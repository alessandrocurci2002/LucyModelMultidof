"""Fit sulle curve di Bode delle stime ibride di ``analisi_dahl_ricampionato.py`` e figure in frequenza per il report.

Richiesta dell'utente (2026-11-01): calcolare il fit anche in Bode. Sui punti misurati in 1-150 Hz, per ogni misura:

* ``complex``: NRMSE su G complessa (lo stesso fit di ``afm.qualita``), intera banda e 1-5, 5-30, 30-150 Hz;
* ``magnitude``: NRMSE sul modulo in dB;
* ``phase``: 100 (1 - ||phi - phi_mis|| / ||phi_mis - media||), con differenze e scarti riportati in (-180, 180] gradi.

Nessuna stima qui: i parametri vengono da ``results_fork_modes_dahl_resampled/summary.json``. Scrive ``bode_fit.csv``,
``sampling.json`` (punti per banda e scarto fra griglia uniforme e griglia misurata) e le figure: Bode con i fit in
legenda, funzione descrittiva dei cuscinetti.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_fit_bode_ricampionato.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))

import matplotlib.pyplot as plt  # noqa: E402

import analisi_dahl_ricampionato as adr  # noqa: E402
import analisi_forcella_lugre as afl  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402

CARTELLA = adr.CARTELLA
VARIANTI = {"report": "ln f width weights (previous report)", "uniform_lnf": "uniform in ln f", "uniform_f": "uniform in f"}
BANDE = {"1-5": (1.0, 5.0), "5-30": (5.0, 30.0), "30-150": (30.0, 150.0)}


def nrmse(y: np.ndarray, y_mod: np.ndarray) -> float:
    return float(100 * (1 - np.linalg.norm(y_mod - y) / np.linalg.norm(y - np.mean(y))))


def avvolgi(gradi: np.ndarray) -> np.ndarray:
    return (gradi + 180.0) % 360.0 - 180.0


def fit_fase(G: np.ndarray, G_mod: np.ndarray) -> float:
    """Fit sulla fase: scarto del modello e dispersione della misura attorno alla media circolare, entrambi avvolti."""
    m = G / np.abs(G)
    scarto = np.rad2deg(np.angle(G_mod / G))
    dispersione = avvolgi(np.rad2deg(np.angle(m)) - np.rad2deg(np.angle(np.mean(m))))
    return float(100 * (1 - np.linalg.norm(scarto) / np.linalg.norm(dispersione)))


def fit_bode(d, G_mod: np.ndarray) -> dict[str, float]:
    out = {"complex": nrmse(d.G, G_mod), "magnitude": nrmse(pt._db(d.G), pt._db(G_mod)), "phase": fit_fase(d.G, G_mod)}
    for nome, (a, b) in BANDE.items():
        k = (d.f >= a) & (d.f <= b)
        out[f"complex_{nome}"] = nrmse(d.G[k], G_mod[k]) if k.sum() >= 3 else float("nan")
        out[f"magnitude_{nome}"] = nrmse(pt._db(d.G[k]), pt._db(G_mod[k])) if k.sum() >= 3 else float("nan")
        out[f"phase_{nome}"] = fit_fase(d.G[k], G_mod[k]) if k.sum() >= 3 else float("nan")
    return out


def figura_cuscinetti(curve, seg):
    fig, assi = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    for j, (asse, (d, A, kk)) in enumerate(curve.items()):
        for i, (et, (K, scala)) in enumerate(kk.items(), 1):
            assi[0, j].plot(d.f, np.abs(K), lw=1.4, label=et, **pt._stile(i))
            assi[1, j].plot(d.f, K.imag / K.real, lw=1.4, **pt._stile(i))
            assi[2, j].axhline(scala, lw=1.0, **pt._stile(i))
        assi[2, j].plot(d.f, 1e6 * A, color=pt.INCHIOSTRO, lw=1.2, label="bearing amplitude from the measured traces")
        assi[0, j].set_title(f"{asse} ({d.etichetta})", fontsize=9)
        for ax in assi[:, j]:
            ax.set_xscale("log"); ax.set_xlim(*afl.BANDA)
        assi[0, j].set_yscale("log"); assi[2, j].set_yscale("log"); assi[2, j].set_xlabel("frequency [Hz]")
        assi[0, j].legend(fontsize=7); assi[2, j].legend(fontsize=7)
    assi[0, 0].set_ylabel("|K_b| [N m/rad]"); assi[1, 0].set_ylabel("loss factor Im K_b / Re K_b")
    assi[2, 0].set_ylabel("A and T_c/sigma0 (horizontal) [urad]")
    pt._titolo(fig, "Bearing describing function at the measured amplitude, for the three hybrid parameter sets", seg)
    return fig


def main() -> None:
    pt.applica_stile()
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    cfg, seg = par.sperimentale, par.segnaposto
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    s = json.loads((CARTELLA / "summary.json").read_text())
    pol = [afm.misura(cfg, fn, "polare", et) for fn, et in afm.MISURE_POLARI]
    az = [afm.misura(cfg, afm.MISURA_AZIMUT[0], "azimut", afm.MISURA_AZIMUT[1])]
    dp, ap = afl.dati_con_ampiezza(cfg, pol)
    da, aa = afl.dati_con_ampiezza(cfg, az)
    righe, campionamento, curve_k = [], {}, {}
    for asse, chiave, cat, mod, dati, amp in (("Elevation", "elevation", cat_pol, fl.polare, dp, ap),
                                              ("Azimuth", "azimuth", cat_az, fl.azimut, da, aa)):
        modelli = {}
        for var, et in VARIANTI.items():
            v = s[chiave]["report_hybrid"]["parameters"] if var == "report" else s[chiave][var]["hybrid"]["parameters"]
            m = mod(cat, v, cfg.Kt, amp)
            q = [fit_bode(d, m.frf_dato(d)) for d in dati]
            righe += [{"axis": chiave, "parameters": var, "measurement": d.etichetta} | x for d, x in zip(dati, q)]
            tot = "/".join(f"{x['complex']:.0f}" for x in q)
            modelli[f"hybrid, {et}: fit {tot} %"] = m
        pt.salva(adr.figura(asse, dati, modelli, seg), CARTELLA / "figures", f"{chiave}_bode_fit")
        plt.close("all")

        d0 = dati[0]
        kk = {}
        for var, et in VARIANTI.items():
            v = s[chiave]["report_hybrid"]["parameters"] if var == "report" else s[chiave][var]["hybrid"]["parameters"]
            kk[et] = (mod(cat, v, cfg.Kt, amp).K_b(d0), 1e6 * v["Tc"] / v["sigma0"])
        curve_k[asse] = (d0, amp[d0.etichetta], kk)

        griglie = {g: [adr.ricampiona(d, amp[d.etichetta], g) for d in dati] for g in adr.GRIGLIE}
        costi = {}
        for var in VARIANTI:
            v = s[chiave]["report_hybrid"]["parameters"] if var == "report" else s[chiave][var]["hybrid"]["parameters"]
            costi[var] = {g: adr.e2(mod(cat, v, cfg.Kt, {u.etichetta: a for u, a in c}), [u for u, _ in c])
                          for g, c in griglie.items()}
        campionamento[chiave] = {"e2": costi}
        for d in dati:
            df_ln = np.geomspace(d.f[0], d.f[-1], d.f.size)
            df_u = np.linspace(d.f[0], d.f[-1], d.f.size)
            campionamento[chiave][d.etichetta] = {
                "n": int(d.f.size), "f_min": float(d.f[0]), "f_max": float(d.f[-1]),
                "old_weight_squared_median": float(np.median(d.peso ** 2)),
                "lnf_max_offset_from_measured_pct": float(100 * np.max(np.abs(df_ln / d.f - 1))),
                "step_lnf_pct": float(100 * (np.exp(np.log(d.f[-1] / d.f[0]) / (d.f.size - 1)) - 1)),
                "step_f_Hz": float(df_u[1] - df_u[0]),
                "points": {g: {n: int(np.sum((x >= a) & (x < b) if b < 150 else (x >= a) & (x <= b)))
                               for n, (a, b) in BANDE.items()} for g, x in (("uniform_lnf", d.f), ("uniform_f", df_u))}}
    pt.salva(figura_cuscinetti(curve_k, seg), CARTELLA / "figures", "bearing_describing_function_resampled")
    plt.close("all")

    with open(CARTELLA / "bode_fit.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(righe[0]))
        wr.writeheader()
        wr.writerows([{k: (f"{x:.2f}" if isinstance(x, float) else x) for k, x in r.items()} for r in righe])
    (CARTELLA / "sampling.json").write_text(json.dumps(campionamento, indent=1))
    for r in righe:
        print(f"{r['axis']:9s} {r['parameters']:11s} {r['measurement'][:34]:34s} complesso {r['complex']:6.1f}  "
              f"modulo {r['magnitude']:6.1f}  fase {r['phase']:6.1f} | 1-5 Hz {r['complex_1-5']:6.1f} {r['magnitude_1-5']:6.1f} "
              f"{r['phase_1-5']:6.1f}")
    print(json.dumps(campionamento, indent=0))


if __name__ == "__main__":
    main()
