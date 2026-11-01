"""Masse modali sullo statore polare e load di azimut flessibile, stimati sulle misure del 28/10/2026 fino a 150 Hz.

Richiesta dell'utente (2026-10-29): le risonanze polari sotto i 150 Hz sono vibrazioni della forcella eccitate dallo
statore del motore e lette dall'encoder, montato sulla stessa struttura. Modello data-driven
(:mod:`gimbal_mdof.modi_forcella`):

* polare: catena H7 con masse modali agganciate allo statore, fittato sulle misure nominali (lente a 12 mm)
  ``..._600_20261028_1perc.csv``, ``..._30_150_20261028.csv`` e
  ``..._78_20261028_2perc_12mm_tightenedlaser_newsupport_4A_80degs_16screws.csv``;
* azimut: catena H7 con il load diviso in base e sommita' della forcella, fittato su
  ``BodeMeasurement_Azimuth_14vit_600_20261028.csv``.

Inerzie fisse al CAD (telescopio, load, rotori); stimati attacco del telescopio, vincolo dei cuscinetti, ritardo
equivalente e i parametri delle masse (polare) o della forcella (azimut). Stima polare in due passi: prima la massa
a ~9 Hz fra 1 e 40 Hz, poi tutte le masse fra 1 e 150 Hz da partenze con tutti i segni dei residui.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_modi_forcella.py
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")  # solve su tante matrici piccole: i thread rallentano
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))

import matplotlib.pyplot as plt  # noqa: E402

from analisi_nodo_azimut import sistema  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402
from gimbal_mdof import vincolo_terra as vt  # noqa: E402

CARTELLA = "results_fork_modes"
OGGI = "Poli_Polar_identification/"
MISURE_POLARI = (
    (OGGI + "BodeMeasurement_Polar_14vit_600_20261028_1perc.csv", "elevation, 1 % grid (1-600 Hz)"),
    (OGGI + "BodeMeasurement_Polar_14vit_30_150_20261028.csv", "elevation, 1 % grid (30-150 Hz)"),
    (OGGI + "BodeMeasurement_Polar_14vit_78_20261028_2perc_12mm_tightenedlaser_newsupport_4A_80degs_16screws.csv",
     "elevation, 16 screws, 4 A (30-300 Hz)"),
)
MISURA_AZIMUT = (OGGI + "BodeMeasurement_Azimuth_14vit_600_20261028.csv", "azimuth, 3 % grid (1-600 Hz)")
F_MODI = (9.0, 56.3, 58.5, 99.2)
"""Frequenze iniziali delle masse sullo statore: risonanze polari misurate sotto i 150 Hz, Hz."""
BANDA_PRIMO_PASSO = (1.0, 40.0)
"""Banda del primo passo della stima polare (solo la massa a ~9 Hz), Hz."""
MODULI_RESIDUO = (0.3, 1.5)
"""Moduli iniziali dei residui delle masse sopra i 20 Hz, 1/(kg m^2): sotto e sopra 1/J dell'asse."""
F_MAX_RADICI = 200.0
COLORI_MISURE = (pt.INCHIOSTRO, pt.SERIE[2], pt.SERIE[3])


def argomenti() -> argparse.Namespace:
    a = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    a.add_argument("--config", type=Path, default=RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    a.add_argument("--f-min", type=float, default=1.0, help="limite inferiore della banda del fit, Hz")
    a.add_argument("--f-max", type=float, default=150.0, help="limite superiore della banda del fit, Hz")
    a.add_argument("--modi", type=float, nargs="+", default=list(F_MODI),
                   help="frequenze iniziali delle masse sullo statore, Hz")
    a.add_argument("--cartella", default=CARTELLA, help="cartella dei risultati (in sim/gimbal_multidof)")
    return a.parse_args()


def misura(cfg: pr.ParametriSperimentali, file: str, asse: str, etichetta: str) -> spm.Misura:
    return spm.carica_misura(cfg, pr.VoceMisura(file, asse, 0.0, True, etichetta, theta0_segnaposto=True))


# ---------------------------------------------------------------------------
# Stime
# ---------------------------------------------------------------------------


def stima_polare(cat: mf.CatenaPolare, par: pr.Parametri, dati_1: list[mf.DatoFRF], dati: list[mf.DatoFRF],
                 Kt: float, f_modi: tuple[float, ...]) -> tuple[mf.Esito | None, mf.Esito]:
    """Se la prima massa e' sotto i 20 Hz, primo passo con solo quella (1-40 Hz); poi tutte le masse: le partenze
    coprono i segni dei residui sopra i 20 Hz e due moduli (:data:`MODULI_RESIDUO`). Senza primo passo le
    partenze coprono anche due valori di k_b."""
    g = par.vincolo_terra.polare if par.vincolo_terra else None
    base = {"k_b": cat.k_b, "c_b": cat.c_b, "k_g": g.k if g else 164.0, "eta_g": g.eta if g else 3.0,
            "tau_d": 2e-4}
    e1 = None
    if f_modi[0] < 20.0:
        p1 = mf.parametri_polari(f_modi[:1], cat)
        partenze_1 = [mf.con_modi({**base, "k_b": kb}, [(f_modi[0], 0.03, a)])
                      for kb in (cat.k_b, 1e6, 1e8) for a in (0.05, 0.3, -0.1)]
        e1 = mf.stima(lambda v: mf.modello_polare(cat, v, 1, Kt), p1, dati_1, partenze_1)
        print(f"polare, primo passo: costo {e1.costo:.3f}", flush=True)
        basi, alti = [e1.valori], f_modi[1:]
    else:
        basi, alti = [{**base, "k_b": kb} for kb in (cat.k_b, 1e6)], f_modi
    n = len(f_modi)
    primo = n - len(alti) + 1
    pn = mf.parametri_polari(f_modi, cat)
    partenze = [mf.con_modi(b, [(f0, 0.01, s * a) for f0, s in zip(alti, segni)], primo=primo)
                for b in basi for segni in itertools.product((1.0, -1.0), repeat=len(alti))
                for a in MODULI_RESIDUO]
    en = mf.stima(lambda v: mf.modello_polare(cat, v, n, Kt), pn, dati, partenze)
    print(f"polare, tutte le masse: costo {en.costo:.3f} ({len(partenze)} partenze)", flush=True)
    return e1, en


def stima_azimut(cat: mf.CatenaAzimut, par: pr.Parametri, dati: list[mf.DatoFRF], Kt: float) -> mf.Esito:
    """Partenze: tre ripartizioni dell'inerzia del load e tre assegnazioni delle risonanze misurate a 49 e 56 Hz
    (base-sommita' e sommita'-telescopio, scambiate, oppure telescopio rigido a 120 Hz)."""
    g = par.vincolo_terra.azimut if par.vincolo_terra else None
    partenze = []
    for alfa in (0.2, 0.4, 0.7):
        J_A = alfa * cat.J_L
        for f_LA, f_P in ((49.0, 56.0), (56.0, 49.0), (56.0, 120.0)):
            partenze.append({"alpha": alfa, "k_LA": J_A * (2 * np.pi * f_LA) ** 2, "c_LA": 1.0,
                             "k_3": cat.J_p * (2 * np.pi * f_P) ** 2, "c_3": cat.c[2], "k_2": cat.k[1],
                             "c_2": cat.c[1], "k_g": g.k if g else 1400.0, "eta_g": g.eta if g else 0.3,
                             "tau_d": 1e-4})
    e = mf.stima(lambda v: mf.modello_azimut(cat, v, Kt), mf.parametri_azimut(cat), dati, partenze)
    print(f"azimut: costo {e.costo:.3f} ({len(partenze)} partenze)", flush=True)
    return e


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def figura_bode(f: np.ndarray, misure: list[spm.Misura], modelli: list[tuple[str, np.ndarray]], titolo: str,
                f_max: float, f_zoom: tuple[float, float], lim_db: tuple[tuple[float, float], tuple[float, float]], *,
                segnaposto: bool):
    fig, assi = plt.subplots(2, 2, figsize=(12.5, 7.2), sharex="col", gridspec_kw={"width_ratios": [1.6, 1.0]})
    for col, f_lim in ((0, (1.0, 600.0)), (1, f_zoom)):
        ax_m, ax_f = assi[0, col], assi[1, col]
        for k, m in enumerate(misure):
            st = {"color": COLORI_MISURE[k % len(COLORI_MISURE)], "marker": "o", "ms": 1.6, "lw": 0.7,
                  "alpha": 0.85}
            ax_m.plot(m.f, m.modulo_db, label=f"measured: {m.etichetta}", **st)
            ax_f.plot(m.f, m.fase_deg, **st)
        for k, (nome, G) in enumerate(modelli):
            st = pt._stile(k)
            ax_m.plot(f, pt._db(G), lw=1.7, label=nome, **st)
            ax_f.plot(f, pt._fase_avvolta(G), lw=1.3, **st)
        for ax in (ax_m, ax_f):
            ax.set_xscale("log")
            ax.set_xlim(*f_lim)
            if f_lim[1] > f_max:
                ax.axvspan(f_max, f_lim[1], color=pt.INCHIOSTRO_TENUE, alpha=0.12, lw=0)
        ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
        ax_m.set_ylim(*lim_db[col])
        ax_f.set_ylabel("phase [deg]")
        ax_f.set_ylim(-180, 180)
        ax_f.set_yticks(np.arange(-180, 181, 90))
        ax_f.set_xlabel("frequency [Hz]")
    assi[0, 0].text(f_max * 1.05, lim_db[0][1] - 3, "outside the\nfit band", fontsize=7, va="top",
                    color=pt.INCHIOSTRO_SECONDARIO)
    assi[0, 0].legend(loc="lower left", fontsize=7)
    assi[0, 1].set_title(f"zoom {f_zoom[0]:g}-{f_zoom[1]:g} Hz", fontsize=9)
    pt._titolo(fig, titolo, segnaposto)
    return fig


def figura_contributi(f: np.ndarray, m: spm.Misura, mod: mf.ModelloPolare, f_max: float, *, segnaposto: bool):
    """Lettura dell'encoder scomposta: disco del rotore e masse sullo statore."""
    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(10.5, 7.0), sharex=True)
    ax_m.plot(m.f, m.modulo_db, color=pt.INCHIOSTRO_TENUE, marker="o", ms=1.4, lw=0.6,
              label=f"measured: {m.etichetta}")
    ax_f.plot(m.f, m.fase_deg, color=pt.INCHIOSTRO_TENUE, marker="o", ms=1.4, lw=0.6)
    G = mod.frf(f)
    ax_m.plot(f, pt._db(G), color=pt.INCHIOSTRO, lw=2.0, label="model: sum (encoder reading)")
    ax_f.plot(f, pt._fase_avvolta(G), color=pt.INCHIOSTRO, lw=1.6)
    # rotore in grigio, masse con la palette (4 colori, 4 tratti)
    for k, (nome, Gk) in enumerate(mod.contributi(f).items()):
        st = {"color": pt.INCHIOSTRO_SECONDARIO, "linestyle": "-"} if k == 0 else pt._stile(k - 1)
        ax_m.plot(f, pt._db(Gk), lw=1.3, label=nome, **st)
        ax_f.plot(f, pt._fase_avvolta(Gk), lw=1.1, **st)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(1.0, 300.0)
        ax.axvspan(f_max, 300.0, color=pt.INCHIOSTRO_TENUE, alpha=0.12, lw=0)
    ax_m.set_ylim(-110, -25)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7, ncol=2)
    pt._titolo(fig, "Elevation model: encoder reading = encoder disc (rotor) minus the head motion driven by the "
                    "stator masses", segnaposto)
    return fig


def figura_errori(curve: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]], f_max: float, *, segnaposto: bool):
    """Errore del modello sulle misure nella banda del fit: modulo (dB) e fase (gradi)."""
    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(10.5, 6.2), sharex=True)
    for k, (nome, f, Gm, G) in enumerate(curve):
        e = np.log(Gm / G)
        st = pt._stile(k)
        ax_m.plot(f, 20 / np.log(10) * e.real, lw=1.0, marker="o", ms=1.4, label=nome, **st)
        ax_f.plot(f, np.rad2deg(e.imag), lw=1.0, marker="o", ms=1.4, **st)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(1.0, f_max)
        ax.axhline(0.0, color=pt.INCHIOSTRO_SECONDARIO, lw=0.8)
    ax_m.set_ylim(-12, 12)
    ax_f.set_ylim(-90, 90)
    ax_m.set_ylabel("model - measured [dB]")
    ax_f.set_ylabel("model - measured [deg]")
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7)
    pt._titolo(fig, f"Error of the fitted models on the measurements, fit band up to {f_max:g} Hz", segnaposto)
    return fig


# ---------------------------------------------------------------------------
# Tabelle
# ---------------------------------------------------------------------------


def righe_parametri(e: mf.Esito, iniziali: dict[str, float]) -> list[list[str]]:
    out = [["parameter", "unit", "estimate", "start_or_H7", "lower", "upper", "std_percent", "at_bound", "note"]]
    for p in e.parametri:
        v = e.valori[p.nome]
        out.append([p.nome, p.unita, f"{v:.5g}", f"{iniziali.get(p.nome, float('nan')):.5g}", f"{p.lo:.4g}",
                    f"{p.hi:.4g}", f"{100 * e.std_rel[p.nome]:.1f}", "yes" if e.al_limite(p.nome) else "",
                    p.nota])
    return out


def qualita(nome: str, d: mf.DatoFRF, G: np.ndarray) -> dict:
    e = np.log(G / d.G)
    return {"model": nome, "measurement": d.etichetta,
            **{f"fit_{k}": round(v, 1) for k, v in spm.fit_bande(d.f, d.G, G).items()},
            "rms_dB": round(float(np.sqrt(np.mean((20 / np.log(10) * e.real) ** 2))), 2),
            "rms_deg": round(float(np.sqrt(np.mean(np.rad2deg(e.imag) ** 2))), 1)}


def caratteristiche(etichetta: str, f: np.ndarray, G: np.ndarray, banda: tuple[float, float]) -> list[list[str]]:
    k = (f >= banda[0]) & (f <= banda[1])
    return [[etichetta, tipo, f"{fx:.2f}", f"{g:.1f}"]
            for tipo, fx, g in spm.estremi(f[k], pt._db(G[k]), prominenza_db=3.0, finestra_dec=0.05)]


def scrivi_csv(percorso: Path, righe: list[list]) -> None:
    with open(percorso, "w", newline="") as fh:
        csv.writer(fh).writerows(righe)


# ---------------------------------------------------------------------------


def main() -> None:
    args = argomenti()
    pt.applica_stile()
    par = pr.parametri_da_config(args.config)
    cfg = par.sperimentale
    assert cfg is not None
    Kt, seg = cfg.Kt, par.segnaposto
    banda = (args.f_min, args.f_max)
    out = RADICE / args.cartella
    f_modi = tuple(args.modi)
    for d in ("figures", "tables"):
        (out / d).mkdir(parents=True, exist_ok=True)

    ml = lz.linearizza(sistema(par).accoppiato, 0.0, 0.0)
    cat_az, cat_pol = mf.catene_h7(ml)
    ml_v = vt.applica(ml, par)
    pol = [misura(cfg, fn, "polare", et) for fn, et in MISURE_POLARI]
    az = misura(cfg, *MISURA_AZIMUT[:1], "azimut", MISURA_AZIMUT[1])
    dati_pol = [mf.dato(m.etichetta, m.f, m.G, banda) for m in pol]
    dati_pol_1 = [d for d in (mf.dato(m.etichetta, m.f, m.G, BANDA_PRIMO_PASSO) for m in pol) if d.f.size >= 3]
    dati_az = [mf.dato(az.etichetta, az.f, az.G, banda)]

    # --- stime ------------------------------------------------------------------------------------------------
    e1, ep = stima_polare(cat_pol, par, dati_pol_1, dati_pol, Kt, f_modi)
    ea = stima_azimut(cat_az, par, dati_az, Kt)
    n = len(f_modi)
    mod_p = mf.modello_polare(cat_pol, ep.valori, n, Kt)
    mod_a = mf.modello_azimut(cat_az, ea.valori, Kt)
    # variante con i cuscinetti a terra (come H7), stessi parametri ristimati: quanto conta il collegamento
    ep_terra = mf.stima(lambda v: mf.modello_polare(cat_pol, v, n, Kt, cuscinetti_su_statore=False),
                        mf.parametri_polari(f_modi, cat_pol), dati_pol, [ep.valori])

    # --- figure -----------------------------------------------------------------------------------------------
    f = np.geomspace(1.0, 600.0, 5000)
    h7 = lambda asse: spm.frf_velocita(ml_v, asse, f, Kt)[spm.canali(asse)[0]]  # noqa: E731
    figure = pt.salva(figura_bode(
        f, pol, [("H7 (configuration)", h7("polare")), (f"fitted: H7 + {n} stator masses", mod_p.frf(f))],
        f"Elevation plant, encoder velocity / current: H7 with modal masses on the motor stator, fit "
        f"{banda[0]:g}-{banda[1]:g} Hz", banda[1], (30.0, 150.0), ((-110, -10), (-75, -35)), segnaposto=seg),
        out / "figures", "bode_polar_fit")
    figure += pt.salva(figura_bode(
        f, [az], [("H7 (configuration, rigid load)", h7("azimut")), ("fitted: flexible load", mod_a.frf(f))],
        f"Azimuth plant, encoder velocity / current: load split into fork base and fork top, fit "
        f"{banda[0]:g}-{banda[1]:g} Hz", banda[1], (30.0, 150.0), ((-95, 0), (-65, -25)), segnaposto=seg),
        out / "figures", "bode_azimuth_fit")
    figure += pt.salva(figura_contributi(f, pol[0], mod_p, banda[1], segnaposto=seg), out / "figures",
                       "polar_stator_masses")
    curve = [(d.etichetta, d.f, mod_p.frf(d.f), d.G) for d in dati_pol] + [
        (d.etichetta, d.f, mod_a.frf(d.f), d.G) for d in dati_az]
    figure += pt.salva(figura_errori(curve, banda[1], segnaposto=seg), out / "figures", "fit_errors")

    # --- tabelle ----------------------------------------------------------------------------------------------
    iniz_p = {"k_b": cat_pol.k_b, "c_b": cat_pol.c_b, **{f"f_{r}": f0 for r, f0 in enumerate(f_modi, 1)}}
    if par.vincolo_terra and par.vincolo_terra.polare:
        iniz_p |= {"k_g": par.vincolo_terra.polare.k, "eta_g": par.vincolo_terra.polare.eta}
    righe_p = righe_parametri(ep, iniz_p)
    righe_p += [[], ["stator mass", "f_Hz", "zeta", "a_1/(kg m^2)", "J_kgm2", "kappa", "k_Nm_rad", "c_Nms_rad"]]
    righe_p += [[r + 1, f"{m.f:.3f}", f"{m.zeta:.4f}", f"{m.a:+.4f}", f"{m.J:.4g}", f"{m.kappa:+.0f}", f"{m.k:.4g}",
                 f"{m.c:.4g}"] for r, m in enumerate(mod_p.modi)]
    scrivi_csv(out / "tables" / "polar_parameters.csv", righe_p)
    iniz_a = {"k_3": cat_az.k[2], "c_3": cat_az.c[2], "k_2": cat_az.k[1], "c_2": cat_az.c[1]}
    if par.vincolo_terra and par.vincolo_terra.azimut:
        iniz_a |= {"k_g": par.vincolo_terra.azimut.k, "eta_g": par.vincolo_terra.azimut.eta}
    modi_a = mod_a.modi(F_MAX_RADICI)
    righe_a = righe_parametri(ea, iniz_a) + [
        [], ["derived", "value"], ["J_L (fork base) kg m^2", f"{mod_a.J_L:.4g}"],
        ["J_A (fork top) kg m^2", f"{mod_a.J_A:.4g}"], ["J_P (telescope CAD) kg m^2", f"{cat_az.J_p:.4g}"],
        [], ["mode", "f_Hz", "zeta"] + [f"shape_{nodo}" for nodo in mf.NODI_AZIMUT]]
    righe_a += [[k + 1, f"{fm:.2f}", f"{z:.4f}"] + [f"{forma[nodo]:+.3f}" for nodo in mf.NODI_AZIMUT]
                for k, (fm, z, forma) in enumerate(modi_a)]
    scrivi_csv(out / "tables" / "azimuth_parameters.csv", righe_a)
    q = [qualita(f"H7 + {n} stator masses", d, mod_p.frf(d.f)) for d in dati_pol]
    q += [qualita("H7 + stator masses, bearings to ground",
                  d, mf.modello_polare(cat_pol, ep_terra.valori, n, Kt, False).frf(d.f)) for d in dati_pol]
    q += [qualita("H7 configuration", d, spm.frf_velocita(ml_v, "polare", d.f, Kt)["theta_e"]) for d in dati_pol]
    q += [qualita("flexible load", d, mod_a.frf(d.f)) for d in dati_az]
    q += [qualita("H7 configuration", d, spm.frf_velocita(ml_v, "azimut", d.f, Kt)["psi_e"]) for d in dati_az]
    scrivi_csv(out / "tables" / "fit_quality.csv", [list(q[0])] + [list(r.values()) for r in q])
    car = [["curve", "kind", "f_Hz", "magnitude_dB"]]
    for m in pol + [az]:
        car += caratteristiche(f"measured: {m.etichetta}", m.f, m.G, banda)
    car += caratteristiche("model: elevation", f, mod_p.frf(f), banda)
    car += caratteristiche("model: azimuth", f, mod_a.frf(f), banda)
    for nome, mod in (("elevation", mod_p), ("azimuth", mod_a)):
        poli, zeri = mod.radici(F_MAX_RADICI)
        car += [[f"model: {nome}", "pole", f"{fp:.2f}", f"zeta {z:.4f}"] for fp, z in poli]
        car += [[f"model: {nome}", "zero", f"{fz:.2f}", f"zeta {z:.4f}"] for fz, z in zeri]
    scrivi_csv(out / "tables" / "features.csv", car)
    scrivi_csv(out / "tables" / "starts.csv", [["axis", "start", "cost"]] + [
        [asse, json.dumps({k: round(v, 6) for k, v in p0.items()}), f"{c:.4f}"]
        for asse, e in (("elevation, first step", e1), ("elevation", ep), ("azimuth", ea)) if e
        for p0, c in e.partenze])

    # --- riepilogo --------------------------------------------------------------------------------------------
    poli_p, zeri_p = mod_p.radici(F_MAX_RADICI)
    sommario = {
        "config": args.config.name, "fit_band_Hz": list(banda), "stator_mass_starts_Hz": list(f_modi), "Kt_Nm_A": Kt,
        "current_loop_Hz": mf.F_CORRENTE, "loss": f"{mf.PERDITA}, f_scale {mf.SCALA_PERDITA}",
        "elevation": {
            "measurements": [fn for fn, _ in MISURE_POLARI], "cost": ep.costo,
            "cost_bearings_to_ground": ep_terra.costo,
            "parameters": ep.valori, "std_rel": ep.std_rel,
            "at_bound": [p.nome for p in ep.parametri if ep.al_limite(p.nome)],
            "stator_masses": [{"f_Hz": m.f, "zeta": m.zeta, "a": m.a, "J_kgm2": m.J, "kappa": m.kappa,
                               "k_Nm_rad": m.k, "c_Nms_rad": m.c} for m in mod_p.modi],
            "axis_inertia_kgm2": cat_pol.J, "poles": poli_p, "zeros": zeri_p},
        "azimuth": {
            "measurement": MISURA_AZIMUT[0], "cost": ea.costo, "parameters": ea.valori, "std_rel": ea.std_rel,
            "at_bound": [p.nome for p in ea.parametri if ea.al_limite(p.nome)],
            "J_L_kgm2": mod_a.J_L, "J_A_kgm2": mod_a.J_A, "J_P_kgm2": cat_az.J_p,
            "modes": [{"f_Hz": fm, "zeta": z, "shape": forma} for fm, z, forma in modi_a],
            "starts": sorted(c for _, c in ea.partenze)},
        "fit_quality": q,
        "h7_chains": {"elevation_roles": list(cat_pol.ruoli), "elevation_M_diag": np.diag(cat_pol.M).tolist(),
                      "k_b_H7": cat_pol.k_b, "c_b_H7": cat_pol.c_b,
                      "azimuth": {"J_e": cat_az.J_e, "J_m": cat_az.J_m, "J_L": cat_az.J_L, "J_p": cat_az.J_p,
                                  "k": list(cat_az.k), "c": list(cat_az.c)}},
    }
    (out / "summary.json").write_text(json.dumps(sommario, indent=2))
    (out / "summary.md").write_text(riepilogo_md(sommario, mod_p, mod_a, ep, ea, q, figure))
    print(json.dumps({k: sommario[k] for k in ("elevation", "azimuth")}, indent=1, default=float)[:6000])


def riepilogo_md(s: dict, mod_p: mf.ModelloPolare, mod_a: mf.ModelloAzimut, ep: mf.Esito, ea: mf.Esito,
                 q: list[dict], figure: list[Path]) -> str:
    vp, va = ep.valori, ea.valori
    md = [f"# Masse modali sullo statore polare e load di azimut flessibile ({s['config']})", "",
          f"Banda del fit {s['fit_band_Hz'][0]:g}-{s['fit_band_Hz'][1]:g} Hz, errore logaritmico complesso, perdita "
          f"{s['loss']}; anello di corrente {s['current_loop_Hz']:g} Hz e ritardo equivalente; inerzie fisse al "
          "CAD. Modello: `src/gimbal_mdof/modi_forcella.py`.", "",
          "## Polare: masse sullo statore", "",
          "Misure: " + ", ".join(f"`{Path(fn).name}`" for fn in s["elevation"]["measurements"]), "",
          "| massa | f [Hz] | zeta | a = kappa/J [1/(kg m^2)] | J [kg m^2] | kappa | k [N m/rad] |",
          "|---|---|---|---|---|---|---|"]
    md += [f"| {r + 1} | {m.f:.2f} | {m.zeta:.4f} | {m.a:+.3f} | {m.J:.3g} | {m.kappa:+.0f} | {m.k:.3g} |"
           for r, m in enumerate(mod_p.modi)]
    md += ["", f"- inerzia rigida dell'asse (CAD): {s['elevation']['axis_inertia_kgm2']:.4f} kg m^2 (1/J = "
           f"{1 / s['elevation']['axis_inertia_kgm2']:.3f})",
           f"- telescopio - mozzo: k_b = {vp['k_b']:.4g} N m/rad (H7 {s['h7_chains']['k_b_H7']:.4g}), "
           f"c_b = {vp['c_b']:.4g} N m s/rad",
           f"- cuscinetti sulla struttura dello statore: k_g = {vp['k_g']:.4g} N m/rad, eta = {vp['eta_g']:.3g}",
           f"- ritardo equivalente: {1e3 * vp['tau_d']:.3f} ms",
           f"- costo {ep.costo:.3f}; con i cuscinetti a terra (come H7) {s['elevation']['cost_bearings_to_ground']:.3f}",
           "- parametri al limite: " + (", ".join(f"`{n}`" for n in s["elevation"]["at_bound"]) or "nessuno"),
           "- poli del modello (Hz, zeta): " + ", ".join(f"{f:.1f} ({z:.3f})" for f, z in s["elevation"]["poles"]),
           "- zeri del modello (Hz, zeta; zeta < 0 = fase non minima): "
           + ", ".join(f"{f:.1f} ({z:+.3f})" for f, z in s["elevation"]["zeros"]), "",
           "## Azimut: load flessibile", "", f"Misura: `{Path(s['azimuth']['measurement']).name}`", "",
           f"- base della forcella J_L = {mod_a.J_L:.4f}, sommita' J_A = {mod_a.J_A:.4f} kg m^2 "
           f"(alpha = {va['alpha']:.3f}); telescopio J_P = {s['azimuth']['J_P_kgm2']:.4f} kg m^2 (CAD)",
           f"- base - sommita': k_LA = {va['k_LA']:.4g} N m/rad, c_LA = {va['c_LA']:.4g} N m s/rad",
           f"- sommita' - telescopio: k_3 = {va['k_3']:.4g} N m/rad (H7 {s['h7_chains']['azimuth']['k'][2]:.4g}), "
           f"c_3 = {va['c_3']:.4g} N m s/rad",
           f"- albero motore - base: k_2 = {va['k_2']:.4g} N m/rad (H7 {s['h7_chains']['azimuth']['k'][1]:.4g}), "
           f"c_2 = {va['c_2']:.4g}",
           f"- cuscinetti: k_g = {va['k_g']:.4g} N m/rad, eta = {va['eta_g']:.3g}; ritardo {1e3 * va['tau_d']:.3f} ms",
           f"- costo {ea.costo:.3f}; costi delle partenze: " + ", ".join(f"{c:.3f}" for c in s["azimuth"]["starts"]),
           "- parametri al limite: " + (", ".join(f"`{n}`" for n in s["azimuth"]["at_bound"]) or "nessuno"),
           "- modi (Hz, zeta, forma e/m/L/A/P): " + "; ".join(
               f"{m['f_Hz']:.1f} ({m['zeta']:.3f}: " + " ".join(f"{v:+.2f}" for v in m["shape"].values()) + ")"
               for m in s["azimuth"]["modes"]), "",
           "## Qualita' del fit (NRMSE %, nella banda del fit; RMS degli errori)", "",
           "| modello | misura | tutta | 1-5 Hz | 5-100 Hz | >100 Hz | RMS dB | RMS deg |", "|---|---|---|---|---|---|---|---|"]
    md += [f"| {r['model']} | {r['measurement']} | {r['fit_all']} | {r['fit_1-5 Hz']} | {r['fit_5-100 Hz']} | "
           f"{r['fit_>100 Hz']} | {r['rms_dB']} | {r['rms_deg']} |" for r in q]
    md += ["", "Figure: " + ", ".join(f"`figures/{p.name}`" for p in figure if p.suffix == ".png"),
           "Tabelle: `tables/polar_parameters.csv`, `tables/azimuth_parameters.csv`, `tables/fit_quality.csv`, "
           "`tables/features.csv`, `tables/starts.csv`"]
    return "\n".join(md) + "\n"


if __name__ == "__main__":
    main()
