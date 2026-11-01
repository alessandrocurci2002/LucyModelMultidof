"""Stima modale ai minimi quadrati con i cuscinetti ad attrito di LuGre o di Dahl (funzione descrittiva) al posto di
k_g (1 + j eta).

Richiesta dell'utente (2026-10-31): stessi vincoli dell'ultima iterazione salvata
(``saved_fork_modes_modal_20261031_cad20_kt10``): zeta_k <= 0.05, inerzie e rigidezze del CAD entro +-20 %, K_t entro
+-10 % per asse, telescopio rigido sull'anello, minimi quadrati sulle parti reale e immaginaria in 1-150 Hz. I
parametri di LuGre (sigma0, T_c, T_s/T_c, v_s, sigma1, sigma2) sono liberi, come lo erano k_g e eta_g. Modello in
:mod:`gimbal_mdof.forcella_lugre`; ampiezza dei cuscinetti dalla prova (|G_mis| I_exc |S_mis| / omega).

Prima partenza: la stima lineare isteretica salvata (parametri strutturali) con LuGre dall'identificazione nel tempo
dei doppietti V2 (``closedLoopModel/results_v2_identification``); altre partenze variano sigma0 e T_s/T_c (LuGre) o
sigma0 e T_c (Dahl: T_s = T_c, sigma1 = sigma2 = 0, richiesta successiva dello stesso giorno).

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_forcella_lugre.py                      # LuGre -> results_fork_modes_lugre
    .venv/bin/python -I analisi_forcella_lugre.py --attrito dahl       # Dahl  -> results_fork_modes_dahl
"""

from __future__ import annotations

import argparse
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

import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import forcella_modale as fm  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402

NOMI = {"lugre": "LuGre", "dahl": "Dahl"}
LINEARE = RADICE / "saved_fork_modes_modal_20261031_cad20_kt10" / "results" / "summary.json"
V2 = RADICE / "closedLoopModel" / "results_v2_identification" / "summary.json"
BANDA = (1.0, 150.0)
VINCOLI = {"zeta_max": 0.05, "cad": 0.20, "cad_J": 0.20, "kt": 0.10}
ZETA_LUGRE = 0.15
"""Smorzamento della risonanza di pre-scorrimento che fissa sigma1 nelle partenze."""


def dati_con_ampiezza(cfg, misure):
    dati, amp = [], {}
    for m in misure:
        d = mf.dato(m.etichetta, m.f, m.G, BANDA)
        _, tracce, _ = spm.carica_bode(cfg.cartella / m.file)
        dB, deg = tracce["Open Loop"]
        L = 10 ** (np.asarray(dB) / 20) * np.exp(1j * np.deg2rad(deg))
        k = (m.f >= BANDA[0]) & (m.f <= BANDA[1])
        amp[m.etichetta] = fl.ampiezza_misurata(m.G[k], L[k], m.f[k], m.eccitazione_A)
        assert np.allclose(m.f[k], d.f)
        dati.append(d)
    return dati, amp


def partenze(lin: dict, v2: dict, J: float, famiglia: str) -> list[dict[str, float]]:
    base = {k: x for k, x in lin.items() if k not in ("k_g", "eta_g", "zeta_1")}
    if famiglia == "dahl":
        return [base | {"sigma0": v2["sigma0"] * f0, "Tc": v2["Tc"] * ft}
                for f0, ft in ((1.0, 1.0), (3.0, 1.0), (1 / 3, 1.0), (1.0, 0.5))]
    out = []
    for f0, ts in ((1.0, v2["Ts_ratio"]), (3.0, v2["Ts_ratio"]), (1 / 3, v2["Ts_ratio"]), (1.0, 1.5)):
        s0 = v2["sigma0"] * f0
        out.append(base | {"sigma0": s0, "Tc": v2["Tc"], "Ts_ratio": ts, "vs": v2["vs"],
                           "sigma1": 2 * ZETA_LUGRE * np.sqrt(s0 * J), "sigma2": 0.01})
    return out


def figura(misure, dati, mod_l, mod_lin, asse, seg, nome):
    fig, assi = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
    for k, (m, d) in enumerate(zip(misure, dati)):
        col = afm.COLORI_MISURE[k % len(afm.COLORI_MISURE)]
        assi[0].plot(d.f, pt._db(d.G), color=col, marker="o", ms=1.6, lw=0.6, alpha=0.8,
                     label=f"measured: {d.etichetta}")
        assi[1].plot(d.f, np.rad2deg(np.angle(d.G)), color=col, marker="o", ms=1.6, lw=0.6, alpha=0.8)
        G = mod_l.frf_dato(d)
        assi[0].plot(d.f, pt._db(G), lw=1.6, label=f"{nome} bearings" if k == 0 else None, **pt._stile(1))
        assi[1].plot(d.f, pt._fase_avvolta(G), lw=1.3, **pt._stile(1))
    f = np.geomspace(*BANDA, 3000)
    Gl = mod_lin.frf(f)
    assi[0].plot(f, pt._db(Gl), lw=1.2, label="linear hysteretic bearings (saved fit)", **pt._stile(2))
    assi[1].plot(f, pt._fase_avvolta(Gl), lw=1.0, **pt._stile(2))
    for ax in assi:
        ax.set_xscale("log"); ax.set_xlim(*BANDA)
    assi[0].set_ylabel("magnitude [dB (rad/s)/A]"); assi[1].set_ylabel("phase [deg]")
    assi[1].set_ylim(-180, 180); assi[1].set_yticks(np.arange(-180, 181, 90)); assi[1].set_xlabel("frequency [Hz]")
    assi[0].legend(loc="lower left", fontsize=7)
    pt._titolo(fig, f"{asse}: modal least-squares fit with {nome} bearings (describing function at the measured "
                    f"bearing amplitude), 1-150 Hz", seg)
    return fig


def figura_cuscinetti(curve, seg, nome):
    fig, assi = plt.subplots(3, 2, figsize=(13, 8.5), sharex=True)
    for j, (asse, (dati, mod, mod_lin)) in enumerate(curve.items()):
        for k, d in enumerate(dati):
            col = afm.COLORI_MISURE[k % len(afm.COLORI_MISURE)]
            K = mod.K_b(d)
            A = mod.ampiezze[d.etichetta]
            assi[0, j].plot(d.f, np.abs(K), color=col, lw=1.2, label=d.etichetta)
            assi[1, j].plot(d.f, K.imag / K.real, color=col, lw=1.2)
            assi[2, j].plot(d.f, 1e6 * A, color=col, lw=1.2)
        kl, el = mod_lin.k_g, mod_lin.eta_g
        assi[0, j].axhline(kl * np.hypot(1, el), **pt._stile(2), lw=1.0, label="linear fit |k_g (1 + j eta_g)|")
        assi[1, j].axhline(el, **pt._stile(2), lw=1.0)
        assi[2, j].axhline(1e6 * mod.lugre["Tc"] / mod.lugre["sigma0"], color=pt.INCHIOSTRO_SECONDARIO, ls=":",
                           lw=0.9, label="pre-sliding scale T_c / sigma_0")
        assi[0, j].set_title(asse, fontsize=9)
        for ax in assi[:, j]:
            ax.set_xscale("log"); ax.set_xlim(*BANDA)
        assi[0, j].set_yscale("log"); assi[2, j].set_yscale("log")
        assi[0, j].legend(fontsize=7); assi[2, j].legend(fontsize=7)
        assi[2, j].set_xlabel("frequency [Hz]")
    assi[0, 0].set_ylabel("|K_b| [N m/rad]"); assi[1, 0].set_ylabel("loss factor Im K_b / Re K_b")
    assi[2, 0].set_ylabel("bearing amplitude [urad]")
    pt._titolo(fig, f"{nome} bearings: describing function K_b(A, f) at the measured amplitude", seg)
    return fig


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--attrito", choices=tuple(NOMI), default="lugre")
    famiglia = ap.parse_args().attrito
    nome = NOMI[famiglia]
    CARTELLA = RADICE / f"results_fork_modes_{famiglia}"
    pt.applica_stile()
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    cfg = par.sperimentale
    Kt, seg = cfg.Kt, par.segnaposto
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    lin = json.loads(LINEARE.read_text())
    v2 = json.loads(V2.read_text())["final"]
    pol = [afm.misura(cfg, fn, "polare", et) for fn, et in afm.MISURE_POLARI]
    az = [afm.misura(cfg, *afm.MISURA_AZIMUT[:1], "azimut", afm.MISURA_AZIMUT[1])]
    dp, ap = dati_con_ampiezza(cfg, pol)
    da, aa = dati_con_ampiezza(cfg, az)
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)

    esiti, sommario = {}, {"friction": nome, "constraints": VINCOLI, "fit_band_Hz": list(BANDA), "linear_fit": str(LINEARE.relative_to(RADICE)),
                            "amplitude": "|G_mis| I_exc |S_mis| / omega, S_mis = 1 / (1 + L_mis)"}
    curve = {}
    for asse, chiave, cat, mod, parl, dati, amp, misure, J in (
            ("Elevation", "elevation", cat_pol, fl.polare, fl.parametri_polari(famiglia, **VINCOLI), dp, ap, pol,
             float(np.sum(lin["elevation"]["M_diag"][:4]))),
            ("Azimuth", "azimuth", cat_az, fl.azimut, fl.parametri_azimut(famiglia, **VINCOLI), da, aa, az,
             float(np.sum(lin["azimuth"]["M_diag"]))) ):
        p0 = partenze(lin[chiave]["parameters"], v2[chiave], J, famiglia)
        e = fm.stima(lambda v, c=cat, a=amp, mo=mod: mo(c, v, Kt, a), parl, dati, p0, max_nfev=800)
        print(f"{asse}: costo {e.costo:.6g} (lineare {lin[chiave]['cost']:.6g}), partenze "
              f"{[round(c, 7) for _, c in e.partenze]}", flush=True)
        m = mod(cat, e.valori, Kt, amp)
        m_lin = (fm.polare(cat, lin[chiave]["parameters"], Kt) if chiave == "elevation"
                 else fm.azimut(cat, lin[chiave]["parameters"], Kt))
        q = [afm.qualita(f"modal LS, {nome} bearings", d, m.frf_dato(d)) for d in dati]
        ql = {r["measurement"]: r for r in lin["fit_quality"]}
        for r in q:
            print("  ", r, "| lineare:", {k: ql[r["measurement"]][k] for k in ("fit_all", "fit_1-5 Hz", "rms_dB")})
        modi = [{"f_Hz": x.f, "zeta": x.zeta, "residue": x.residuo, "shape": x.forma} for x in m.struttura.modi()]
        sommario[chiave] = {"cost": e.costo, "cost_linear": lin[chiave]["cost"], "nfev": e.nfev,
                            "parameters": e.valori, "std_rel": e.std_rel,
                            "at_bound": [p.nome for p in e.parametri if e.al_limite(p.nome)],
                            "starts": [c for _, c in e.partenze], "fit_quality": q,
                            "fit_quality_linear": [ql[r["measurement"]] for r in q],
                            "structure_modes_without_bearings": modi,
                            "presliding_scale_urad": 1e6 * e.valori["Tc"] / e.valori["sigma0"],
                            "amplitude_range_urad": [float(1e6 * min(a.min() for a in amp.values())),
                                                     float(1e6 * max(a.max() for a in amp.values()))]}
        esiti[chiave] = e
        curve[asse] = (dati, m, m_lin)
        pt.salva(figura(misure, dati, m, m_lin, asse, seg, nome), CARTELLA / "figures", f"{chiave}_bode_{famiglia}")
    pt.salva(figura_cuscinetti(curve, seg, nome), CARTELLA / "figures", "bearing_describing_function")
    plt.close("all")
    (CARTELLA / "summary.json").write_text(json.dumps(sommario, indent=1, default=float))
    for chiave in ("elevation", "azimuth"):
        print(chiave, json.dumps({k: sommario[chiave][k] for k in ("parameters", "at_bound")}, default=float))


if __name__ == "__main__":
    main()
