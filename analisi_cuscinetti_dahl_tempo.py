"""Funzioni di trasferimento con i cuscinetti di Dahl identificati nel tempo (doppietti V2) al posto di quelli stimati in
frequenza.

Richiesta dell'utente (2026-10-31): confrontare, in azimut e in polare, le funzioni di trasferimento ottenute con i
parametri di Dahl stimati nel dominio del tempo (``closedLoopModel/results_v2_identification``, ``final``: sigma0, T_c,
T_s/T_c, v_s; sigma1 = sigma2 = 0) come caratteristica dei cuscinetti. Due varianti:

* ``sostituiti``: struttura della stima Dahl in frequenza (``results_fork_modes_dahl``), solo i cuscinetti cambiano;
* ``ristimati``: cuscinetti fissi ai valori del tempo, parte strutturale ristimata con gli stessi vincoli (zeta_k <= 0.05,
  CAD +-20 %, K_t +-10 %, minimi quadrati sulle parti reale e immaginaria in 1-150 Hz);
* ``ibrida`` (richiesta successiva): Dahl puro (T_s = T_c) con T_c fisso al valore del tempo, che e' coerente con
  quello stimato in frequenza, e sigma0 stimato sulle FRF insieme alla struttura, con gli stessi vincoli. I doppietti
  misurano lo strisciamento (T_c), non il pre-scorrimento (sigma0, fissato nella selezione V2); le FRF il contrario.

Funzione descrittiva all'ampiezza misurata del cuscinetto, come nelle stime in frequenza.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_cuscinetti_dahl_tempo.py
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))

import matplotlib.pyplot as plt  # noqa: E402

import analisi_forcella_lugre as afl  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import forcella_modale as fm  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402

DAHL_F = RADICE / "results_fork_modes_dahl" / "summary.json"
CARTELLA = RADICE / "results_fork_modes_dahl" / "time_domain_bearings"
NOMI_TEMPO = ("sigma0", "Tc", "Ts_ratio", "vs")


def figura(asse, dati, curve, seg):
    fig, assi = plt.subplots(2, 1, figsize=(11, 7.8), sharex=True)
    for k, d in enumerate(dati):
        col = afm.COLORI_MISURE[k % len(afm.COLORI_MISURE)]
        assi[0].plot(d.f, pt._db(d.G), color=col, marker="o", ms=1.5, lw=0.6, alpha=0.8, label=f"measured: {d.etichetta}")
        assi[1].plot(d.f, np.rad2deg(np.angle(d.G)), color=col, marker="o", ms=1.5, lw=0.6, alpha=0.8)
    for i, (nome, Gs) in enumerate(curve.items(), 1):
        for k, (d, G) in enumerate(zip(dati, Gs)):
            assi[0].plot(d.f, pt._db(G), lw=1.5, label=nome if k == 0 else None, **pt._stile(i))
            assi[1].plot(d.f, pt._fase_avvolta(G), lw=1.2, **pt._stile(i))
    for ax in assi:
        ax.set_xscale("log"); ax.set_xlim(*afl.BANDA)
    assi[0].set_ylabel("magnitude [dB (rad/s)/A]"); assi[1].set_ylabel("phase [deg]")
    assi[1].set_ylim(-180, 180); assi[1].set_yticks(np.arange(-180, 181, 90)); assi[1].set_xlabel("frequency [Hz]")
    assi[0].legend(loc="lower left", fontsize=7)
    pt._titolo(fig, f"{asse}: Dahl bearings identified in the time domain (V2 doublets) vs identified in frequency", seg)
    return fig


def figura_cuscinetti(curve_k, seg):
    fig, assi = plt.subplots(2, 2, figsize=(13, 7), sharex=True)
    for j, (asse, (d, kk)) in enumerate(curve_k.items()):
        for i, (nome, K) in enumerate(kk.items(), 1):
            assi[0, j].plot(d.f, np.abs(K), lw=1.4, label=nome, **pt._stile(i))
            assi[1, j].plot(d.f, K.imag / K.real, lw=1.4, **pt._stile(i))
        assi[0, j].set_title(f"{asse} ({d.etichetta})", fontsize=9)
        for ax in assi[:, j]:
            ax.set_xscale("log"); ax.set_xlim(*afl.BANDA)
        assi[0, j].set_yscale("log"); assi[1, j].set_xlabel("frequency [Hz]")
        assi[0, j].legend(fontsize=7)
    assi[0, 0].set_ylabel("|K_b| [N m/rad]"); assi[1, 0].set_ylabel("loss factor Im K_b / Re K_b")
    pt._titolo(fig, "Bearing describing function at the measured amplitude: time-domain, hybrid and frequency-domain Dahl",
               seg)
    return fig


def main() -> None:
    pt.applica_stile()
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    cfg = par.sperimentale
    Kt, seg = cfg.Kt, par.segnaposto
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    df = json.loads(DAHL_F.read_text())
    lin = json.loads(afl.LINEARE.read_text())
    v2 = json.loads(afl.V2.read_text())["final"]
    pol = [afm.misura(cfg, fn, "polare", et) for fn, et in afm.MISURE_POLARI]
    az = [afm.misura(cfg, afm.MISURA_AZIMUT[0], "azimut", afm.MISURA_AZIMUT[1])]
    dp, ap = afl.dati_con_ampiezza(cfg, pol)
    da, aa = afl.dati_con_ampiezza(cfg, az)
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)

    sommario, curve_k = {"time_domain_dahl": {k: {n: v2[k][n] for n in NOMI_TEMPO} for k in ("azimuth", "elevation")}}, {}
    for asse, chiave, cat, mod, dati, amp in (("Elevation", "elevation", cat_pol, fl.polare, dp, ap),
                                              ("Azimuth", "azimuth", cat_az, fl.azimut, da, aa)):
        tempo = {n: v2[chiave][n] for n in NOMI_TEMPO}
        vf = df[chiave]["parameters"]
        m_f = mod(cat, vf, Kt, amp)
        m_s = mod(cat, vf | tempo, Kt, amp)
        # ristima della sola struttura con i cuscinetti fissi
        parl = [p for p in fl.parametri_polari("dahl", **afl.VINCOLI) if p.nome not in tempo] if chiave == "elevation" \
            else [p for p in fl.parametri_azimut("dahl", **afl.VINCOLI) if p.nome not in tempo]
        e = fm.stima(lambda v, c=cat, a=amp, mo=mod, t=tempo: mo(c, v | t, Kt, a), parl, dati, [vf], max_nfev=800)
        m_r = mod(cat, e.valori | tempo, Kt, amp)
        # stima ibrida: T_c dal tempo, sigma0 (e struttura) dalle FRF
        Tc_t = {"Tc": tempo["Tc"]}
        parh = [p for p in (fl.parametri_polari("dahl", **afl.VINCOLI) if chiave == "elevation"
                            else fl.parametri_azimut("dahl", **afl.VINCOLI)) if p.nome != "Tc"]
        eh = fm.stima(lambda v, c=cat, a=amp, mo=mod, t=Tc_t: mo(c, v | t, Kt, a), parh, dati,
                      [vf, vf | {"sigma0": tempo["sigma0"]}], max_nfev=800)
        m_h = mod(cat, eh.valori | Tc_t, Kt, amp)
        m_lin = (fm.polare(cat, lin[chiave]["parameters"], Kt) if chiave == "elevation"
                 else fm.azimut(cat, lin[chiave]["parameters"], Kt))
        modelli = {"Dahl identified in frequency (sigma0, T_c fitted)": m_f.frf_dato,
                   "hybrid Dahl: T_c from time domain, sigma0 from FRF": m_h.frf_dato,
                   "Dahl from time domain, structure refitted": m_r.frf_dato,
                   "linear hysteretic bearings (saved fit)": lambda d, m=m_lin: m.frf(d.f),
                   "Dahl from time domain, same structure": m_s.frf_dato}
        curve = {n: [g(d) for d in dati] for n, g in modelli.items()}
        in_figura = [n for n in curve if n != "Dahl from time domain, same structure"]
        pt.salva(figura(asse, dati, {n: curve[n] for n in in_figura}, seg), CARTELLA / "figures",
                 f"{chiave}_bode_time_domain_dahl")
        qual = {n: [afm.qualita(n, d, G) for d, G in zip(dati, Gs)] for n, Gs in curve.items()}
        costi = {n: float(0.5 * sum(np.sum((d.peso * np.abs(G - d.G)) ** 2) for d, G in zip(dati, Gs)))
                 for n, Gs in curve.items()}
        sommario[chiave] = {"frequency_domain_dahl": {n: vf[n] for n in ("sigma0", "Tc")},
                            "refit": {"parameters": e.valori, "at_bound": [p.nome for p in e.parametri if e.al_limite(p.nome)],
                                      "cost": e.costo},
                            "hybrid": {"Tc_from_time": tempo["Tc"], "parameters": eh.valori, "std_rel": eh.std_rel,
                                       "at_bound": [p.nome for p in eh.parametri if eh.al_limite(p.nome)],
                                       "cost": eh.costo, "starts": [c for _, c in eh.partenze],
                                       "presliding_scale_urad": 1e6 * tempo["Tc"] / eh.valori["sigma0"]},
                            "cost": costi, "fit_quality": qual,
                            "presliding_scale_urad": {"frequency": 1e6 * vf["Tc"] / vf["sigma0"],
                                                      "time": 1e6 * tempo["Tc"] / tempo["sigma0"]}}
        d0 = dati[0]
        curve_k[asse] = (d0, {"frequency-domain Dahl": m_f.K_b(d0), "hybrid Dahl (T_c time, sigma0 FRF)": m_h.K_b(d0),
                              "time-domain Dahl (V2)": m_s.K_b(d0)})
        print(f"\n{asse}: tempo {tempo} | frequenza sigma0 {vf['sigma0']:.0f}, Tc {vf['Tc']:.3f}", flush=True)
        for n in curve:
            print(f"  {n:52s} costo {costi[n]:.4g}  " + "  ".join(
                f"{q['fit_all']}% ({q['rms_dB']} dB, {q['rms_deg']} deg)" for q in qual[n]), flush=True)
        h = sommario[chiave]["hybrid"]
        print(f"  ibrida: sigma0 {eh.valori['sigma0']:.0f} (std {100 * eh.std_rel['sigma0']:.1f} %), Tc {tempo['Tc']:.3f} fisso, "
              f"g_Kt {eh.valori['g_Kt']:.3f}, pre-scorrimento {h['presliding_scale_urad']:.0f} urad, partenze {h['starts']}, "
              f"al limite {h['at_bound']}")
        print("  ristima, parametri al limite:", sommario[chiave]["refit"]["at_bound"])
        print("  ristima:", {k: round(x, 4) for k, x in e.valori.items() if not k.startswith("zeta")})
    pt.salva(figura_cuscinetti(curve_k, seg), CARTELLA / "figures", "bearing_describing_function_time_vs_frequency")
    plt.close("all")
    (CARTELLA / "summary.json").write_text(json.dumps(sommario, indent=1, default=float))


if __name__ == "__main__":
    main()
