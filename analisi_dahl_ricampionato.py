"""Stima Dahl sulle FRF ricampionate su griglia uniforme, senza pesi per punto.

Richiesta dell'utente (2026-11-01): al posto del peso w_i = sqrt(dln f_i / 0.01), interpolazione lineare fra i punti
misurati (parti reale e immaginaria), campionamento uniforme con lo stesso numero di punti della misura, poi minimo di

    e^2 = sum_i [Re(G - G_mis)^2 + Im(G - G_mis)^2]

con gli stessi vincoli di ``report_dahl_bearings`` (zeta_k <= 0.05, CAD +-20 %, K_t +-10 %, telescopio rigido
sull'anello, 1-150 Hz). Due griglie: uniforme in f e uniforme in ln f (l'integrale in ln f e' cio' che il peso
approssimava). Ogni griglia ha due stime, come nel report: Dahl libero (sigma0, T_c) e ibrida (T_c dai doppietti V2).

L'ampiezza dei cuscinetti e' interpolata sulla stessa griglia. La qualita' e' calcolata sempre sui punti misurati,
quindi e' confrontabile con la stima del report. Interpolazione quadratica: solo confronto con quella lineare.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_dahl_ricampionato.py
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
from scipy.interpolate import make_interp_spline  # noqa: E402

import analisi_forcella_lugre as afl  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import forcella_modale as fm  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402

DAHL_F = RADICE / "results_fork_modes_dahl" / "summary.json"
IBRIDO = RADICE / "results_fork_modes_dahl" / "time_domain_bearings" / "summary.json"
CARTELLA = RADICE / "results_fork_modes_dahl_resampled"
GRIGLIE = {"uniform_lnf": lambda a, b, n: np.geomspace(a, b, n), "uniform_f": lambda a, b, n: np.linspace(a, b, n)}


def ricampiona(d: mf.DatoFRF, A: np.ndarray, griglia: str) -> tuple[mf.DatoFRF, np.ndarray]:
    """Interpolazione lineare di Re, Im e ampiezza su ``d.f.size`` punti uniformi fra gli estremi misurati; peso 1."""
    f = GRIGLIE[griglia](d.f[0], d.f[-1], d.f.size)
    f[[0, -1]] = d.f[[0, -1]]
    G = np.interp(f, d.f, d.G.real) + 1j * np.interp(f, d.f, d.G.imag)
    return mf.DatoFRF(d.etichetta, f, G, np.ones_like(f)), np.interp(f, d.f, A)


def e2(m, dati) -> float:
    return float(sum(np.sum(np.abs(m.frf_dato(d) - d.G) ** 2) for d in dati))


def interpolazione_quadratica(d: mf.DatoFRF, u: mf.DatoFRF, G_mod: np.ndarray) -> dict[str, float]:
    """Scarto fra interpolante quadratica e lineare, rispetto al segnale e al residuo del fit sulla stessa griglia."""
    Gq = make_interp_spline(d.f, d.G.real, k=2)(u.f) + 1j * make_interp_spline(d.f, d.G.imag, k=2)(u.f)
    dq = np.linalg.norm(Gq - u.G)
    return {"quadratic_vs_linear_over_signal_pct": float(100 * dq / np.linalg.norm(u.G)),
            "quadratic_vs_linear_over_fit_residual_pct": float(100 * dq / np.linalg.norm(G_mod - u.G))}


def figura(asse, dati, curve, seg):
    fig, assi = plt.subplots(2, 1, figsize=(11, 7.8), sharex=True)
    for k, d in enumerate(dati):
        col = afm.COLORI_MISURE[k % len(afm.COLORI_MISURE)]
        assi[0].plot(d.f, pt._db(d.G), color=col, marker="o", ms=1.5, lw=0.6, alpha=0.8, label=f"measured: {d.etichetta}")
        assi[1].plot(d.f, np.rad2deg(np.angle(d.G)), color=col, marker="o", ms=1.5, lw=0.6, alpha=0.8)
    for i, (nome, m) in enumerate(curve.items(), 1):
        for k, d in enumerate(dati):
            G = m.frf_dato(d)
            assi[0].plot(d.f, pt._db(G), lw=1.4, label=nome if k == 0 else None, **pt._stile(i))
            assi[1].plot(d.f, pt._fase_avvolta(G), lw=1.1, **pt._stile(i))
    for ax in assi:
        ax.set_xscale("log"); ax.set_xlim(*afl.BANDA)
    assi[0].set_ylabel("magnitude [dB (rad/s)/A]"); assi[1].set_ylabel("phase [deg]")
    assi[1].set_ylim(-180, 180); assi[1].set_yticks(np.arange(-180, 181, 90)); assi[1].set_xlabel("frequency [Hz]")
    assi[0].legend(loc="lower left", fontsize=7)
    pt._titolo(fig, f"{asse}: hybrid Dahl fit, ln f width weights vs uniform resampling of the measurements", seg)
    return fig


def main() -> None:
    pt.applica_stile()
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    cfg = par.sperimentale
    Kt, seg = cfg.Kt, par.segnaposto
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    lin = json.loads(afl.LINEARE.read_text())
    v2 = json.loads(afl.V2.read_text())["final"]
    df, ib = json.loads(DAHL_F.read_text()), json.loads(IBRIDO.read_text())
    pol = [afm.misura(cfg, fn, "polare", et) for fn, et in afm.MISURE_POLARI]
    az = [afm.misura(cfg, afm.MISURA_AZIMUT[0], "azimut", afm.MISURA_AZIMUT[1])]
    dp, ap = afl.dati_con_ampiezza(cfg, pol)
    da, aa = afl.dati_con_ampiezza(cfg, az)
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)

    sommario = {"constraints": afl.VINCOLI, "fit_band_Hz": list(afl.BANDA),
                "cost": "e^2 = sum(Re(G - G_mis)^2 + Im(G - G_mis)^2), linear interpolation, one per grid point",
                "quality": "always on the measured points"}
    for asse, chiave, cat, mod, parl, dati, amp, J in (
            ("Elevation", "elevation", cat_pol, fl.polare, fl.parametri_polari("dahl", **afl.VINCOLI), dp, ap,
             float(np.sum(lin["elevation"]["M_diag"][:4]))),
            ("Azimuth", "azimuth", cat_az, fl.azimut, fl.parametri_azimut("dahl", **afl.VINCOLI), da, aa,
             float(np.sum(lin["azimuth"]["M_diag"])))):
        Tc = {"Tc": v2[chiave]["Tc"]}
        vf, vh = df[chiave]["parameters"], ib[chiave]["hybrid"]["parameters"] | Tc
        parh = [p for p in parl if p.nome != "Tc"]
        rif = mod(cat, vh, Kt, amp)
        s = {"grid_points": {d.etichetta: int(d.f.size) for d in dati},
             "report_hybrid": {"parameters": vh, "fit_quality": [afm.qualita("report hybrid", d, rif.frf_dato(d)) for d in dati]}}
        curve = {"hybrid, ln f width weights (report)": rif}
        for griglia in GRIGLIE:
            coppie = [ricampiona(d, amp[d.etichetta], griglia) for d in dati]
            du = [c[0] for c in coppie]
            au = {c[0].etichetta: c[1] for c in coppie}
            p0 = afl.partenze(lin[chiave]["parameters"], v2[chiave], J, "dahl") + [vf]
            el = fm.stima(lambda v, c=cat, a=au, mo=mod: mo(c, v, Kt, a), parl, du, p0, max_nfev=800)
            eh = fm.stima(lambda v, c=cat, a=au, mo=mod, t=Tc: mo(c, v | t, Kt, a), parh, du,
                          [el.valori, vh, vh | {"sigma0": v2[chiave]["sigma0"]}], max_nfev=800)
            vu = eh.valori | Tc
            mu_l, mu_h = mod(cat, el.valori, Kt, amp), mod(cat, vu, Kt, amp)
            mg_h = mod(cat, vu, Kt, au)
            s[griglia] = {
                "e2": 2 * eh.costo, "e2_report_parameters_same_grid": e2(mod(cat, vh, Kt, au), du),
                "free_dahl": {"parameters": el.valori, "e2": 2 * el.costo, "starts_e2": [2 * c for _, c in el.partenze],
                              "at_bound": [p.nome for p in el.parametri if el.al_limite(p.nome)],
                              "fit_quality": [afm.qualita(griglia, d, mu_l.frf_dato(d)) for d in dati]},
                "hybrid": {"parameters": vu, "std_rel": eh.std_rel, "starts_e2": [2 * c for _, c in eh.partenze],
                           "at_bound": [p.nome for p in eh.parametri if eh.al_limite(p.nome)],
                           "fit_quality": [afm.qualita(griglia, d, mu_h.frf_dato(d)) for d in dati]},
                "interpolation": {d.etichetta: interpolazione_quadratica(d, u, mg_h.frf_dato(u)) for d, u in zip(dati, du)},
                "uniform_step": {u.etichetta: float(u.f[1] - u.f[0]) if griglia == "uniform_f"
                                 else float(np.log(u.f[1] / u.f[0])) for u in du}}
            curve[f"hybrid, {griglia.replace('_', ' ')}"] = mu_h
            print(f"{asse} {griglia}: e2 {2 * eh.costo:.5g} (parametri del report {s[griglia]['e2_report_parameters_same_grid']:.5g}), "
                  f"sigma0 {vu['sigma0']:.0f} (report {vh['sigma0']:.0f}), g_Kt {vu.get('g_Kt')}, partenze "
                  f"{[round(2 * c, 7) for _, c in eh.partenze]}", flush=True)
            for q in s[griglia]["hybrid"]["fit_quality"]:
                print("  ", q, flush=True)
            print("  ", s[griglia]["interpolation"], flush=True)
        sommario[chiave] = s
        pt.salva(figura(asse, dati, curve, seg), CARTELLA / "figures", f"{chiave}_bode_resampled")
        plt.close("all")
        (CARTELLA / "summary.json").write_text(json.dumps(sommario, indent=1, default=float))


if __name__ == "__main__":
    main()
