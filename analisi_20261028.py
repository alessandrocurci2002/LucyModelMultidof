"""Misure del 28/10/2026 (14 viti anello-telescopio) contro il modello H7 e rigidezza anello-telescopio k_b.

Richiesta dell'utente (2026-10-28):

* confronto fra le risposte in frequenza prese oggi; per il polare quelle con la lente interna a 12 mm (valore
  nominale): ``..._1perc.csv`` (griglia all'1 %, file senza suffisso = 12 mm), ``..._600_20261028.csv`` (griglia
  al 3 %) e ``..._30_150_20261028.csv`` (30-150 Hz all'1 %). Per riferimento anche le misure VTesi del 23/10;
* rigidezza del collegamento anello-telescopio (``polar.hub.telescope_attachment.k_b_Nm_rad``) che porta la
  risonanza polare del modello a quella misurata a ~56 Hz, con l'inerzia del telescopio fissa al valore CAD.

Il polo cercato e' quello del canale coppia polare -> encoder polare che si muove con k_b (gli altri poli che
compaiono nel canale, per esempio il modo del telescopio in yaw, si cancellano con uno zero).

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_20261028.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))

import matplotlib.pyplot as plt  # noqa: E402

from analisi_nodo_azimut import sistema  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402
from gimbal_mdof import tf  # noqa: E402
from gimbal_mdof import vincolo_terra as vt  # noqa: E402

CARTELLA = "results_20261028"
OGGI = "Poli_Polar_identification/"
MISURE_POLARI = (
    (OGGI + "BodeMeasurement_Polar_14vit_600_20261028_1perc.csv", "elevation, lens 12 mm, 1 % grid, 28 Oct"),
    (OGGI + "BodeMeasurement_Polar_14vit_600_20261028.csv", "elevation, lens 12 mm, 3 % grid, 28 Oct"),
    (OGGI + "BodeMeasurement_Polar_14vit_30_150_20261028.csv", "elevation, lens 12 mm, 30-150 Hz, 28 Oct"),
)
"""Misure polari di oggi con la lente a 12 mm; la prima (griglia all'1 %) fissa la risonanza di riferimento."""
MISURE_AZIMUT = ((OGGI + "BodeMeasurement_Azimuth_14vit_600_20261028.csv", "azimuth, 3 % grid, 28 Oct"),)
RIFERIMENTO = {"polare": ("VClosedLoop/VTesi/BodeMeasurement_Polar_14viti.csv", "elevation, 14 screws, 23 Oct"),
               "azimut": ("VClosedLoop/VTesi/BodeMeasurement_Azimuth_14viti.csv", "azimuth, 14 screws, 23 Oct")}
BANDA_RISONANZA = (50.0, 62.0)
"""Banda in cui cercare la risonanza polare a ~56 Hz sulla misura all'1 %, Hz."""
BANDA_POLO = (5.0, 400.0)
"""Banda in cui cercare il polo del telescopio sul giunto nel modello, Hz."""
BANDA_ESTREMI = (20.0, 450.0)
TOLLERANZA_CANCELLAZIONE = 1e-3
CHIAVE_KB = ("polar", "hub", "telescope_attachment", "k_b_Nm_rad")
COLORI_MISURE = (pt.INCHIOSTRO, pt.SERIE[2], pt.SERIE[3], pt.INCHIOSTRO_TENUE)


def argomenti() -> argparse.Namespace:
    a = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    a.add_argument("--config", type=Path, default=RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    a.add_argument("--k-sweep", type=float, nargs=3, default=[1.0e3, 1.0e6, 121], metavar=("K_MIN", "K_MAX", "N"),
                   help="scansione logaritmica di k_b [N m/rad]")
    return a.parse_args()


def misura(cfg: pr.ParametriSperimentali, file: str, asse: str, etichetta: str) -> spm.Misura:
    return spm.carica_misura(cfg, pr.VoceMisura(file, asse, 0.0, True, etichetta, theta0_segnaposto=True))


def picco(m: spm.Misura, banda: tuple[float, float]) -> tuple[float, float]:
    """Massimo del modulo in ``banda``, raffinato con una parabola in (log f, dB) sui tre punti attorno."""
    ok = np.flatnonzero((m.f >= banda[0]) & (m.f <= banda[1]))
    i = ok[int(np.argmax(m.modulo_db[ok]))]
    x, y = np.log(m.f[i - 1: i + 2]), m.modulo_db[i - 1: i + 2]
    a, b, c = np.polyfit(x, y, 2)
    return float(np.exp(-b / (2 * a))), float(c - b * b / (4 * a))


def con_kb(percorso: Path, dati: dict, k_b: float) -> pr.Parametri:
    nodo = dati
    *strada, ultima = CHIAVE_KB
    for s in strada:
        nodo = nodo[s]
    nodo[ultima] = float(k_b)
    return pr.parametri_da_config(percorso, dati=dati)


def radici_polari(ml: lz.ModelloLineare) -> tuple[list[float], list[float]]:
    """Poli e zeri (Hz) poco smorzati del canale coppia polare -> encoder polare in :data:`BANDA_POLO`, senza le
    coppie che si cancellano (modi dell'altro asse)."""
    z = tf.zeri_invarianti(ml.A, ml.Bu[:, 1], ml.Cy[ml.riga("theta_e")])
    sel = lambda r: r[(np.imag(r) > 0) & (np.imag(r) > np.abs(np.real(r)))]  # noqa: E731
    p, z = sel(ml.poli), sel(z)
    vivi_p = [r for r in p if not any(abs(r - s) < TOLLERANZA_CANCELLAZIONE * abs(r) for s in z)]
    vivi_z = [r for r in z if not any(abs(r - s) < TOLLERANZA_CANCELLAZIONE * abs(r) for s in p)]
    hz = lambda v: sorted(f for f in (abs(r) / (2 * np.pi) for r in v) if BANDA_POLO[0] <= f <= BANDA_POLO[1])  # noqa: E731
    return hz(vivi_p), hz(vivi_z)


def coppia_telescopio(par: pr.Parametri) -> tuple[float, float]:
    """``(polo, zero)`` del telescopio sul giunto: il primo polo e il primo zero del canale polare."""
    p, z = radici_polari(lz.linearizza(sistema(par).accoppiato, 0.0, 0.0))
    return (p[0] if p else float("nan")), (z[0] if z else float("nan"))


def frf(par: pr.Parametri, asse: str, f: np.ndarray, Kt: float) -> np.ndarray:
    ml = vt.applica(lz.linearizza(sistema(par).accoppiato, 0.0, 0.0), par)
    return spm.frf_velocita(ml, asse, f, Kt)[spm.canali(asse)[0]]


def figura_bode(f: np.ndarray, modelli: list[tuple[str, np.ndarray]], misure: list[spm.Misura], titolo: str,
                f_zoom: tuple[float, float], *, segnaposto: bool, marcatori: tuple[float, ...] = ()):
    fig, assi = plt.subplots(2, 2, figsize=(12.5, 7.2), sharex="col", gridspec_kw={"width_ratios": [1.6, 1.0]})
    for col, f_lim in ((0, (1.0, 600.0)), (1, f_zoom)):
        ax_m, ax_f = assi[0, col], assi[1, col]
        for k, m in enumerate(misure):
            st = {"color": COLORI_MISURE[k % len(COLORI_MISURE)], "marker": "o", "ms": 1.6, "lw": 0.8,
                  "alpha": 0.5 if k == len(misure) - 1 and len(misure) > 1 else 0.9}
            ax_m.plot(m.f, m.modulo_db, label=f"measured: {m.etichetta}", **st)
            ax_f.plot(m.f, m.fase_deg, **st)
        for k, (nome, G) in enumerate(modelli):
            st = pt._stile(k)
            ax_m.plot(f, pt._db(G), lw=1.6, label=nome, **st)
            ax_f.plot(f, pt._fase_avvolta(G), lw=1.3, **st)
        for ax in (ax_m, ax_f):
            ax.set_xscale("log")
            ax.set_xlim(*f_lim)
            for fm in marcatori:
                ax.axvline(fm, color=pt.INCHIOSTRO_SECONDARIO, lw=0.7, ls=":")
        ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
        ax_f.set_ylabel("phase [deg]")
        ax_f.set_ylim(-180, 180)
        ax_f.set_yticks(np.arange(-180, 181, 90))
        ax_f.set_xlabel("frequency [Hz]")
    assi[0, 0].set_ylim(-100, -20)
    assi[0, 1].set_ylim(-80, -30)
    assi[0, 0].legend(loc="lower left", fontsize=7)
    assi[0, 1].set_title(f"zoom {f_zoom[0]:g}-{f_zoom[1]:g} Hz", fontsize=9)
    pt._titolo(fig, titolo, segnaposto)
    return fig


def figura_scansione(sw: np.ndarray, k_b0: float, k_b: float, f_ris: float, f_antires: float, *, segnaposto: bool):
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(sw[:, 0], sw[:, 1], color=pt.SERIE[0], lw=2.0, label="model pole (resonance)")
    ax.plot(sw[:, 0], sw[:, 2], color=pt.SERIE[1], lw=2.0, ls="--", label="model zero (antiresonance)")
    ax.axhline(f_ris, color=pt.SERIE[0], lw=0.9, alpha=0.6)
    ax.axhline(f_antires, color=pt.SERIE[1], lw=0.9, ls="--", alpha=0.6)
    ax.plot([], [], color=pt.INCHIOSTRO_SECONDARIO, lw=0.9, label="measured, lens 12 mm, 28 Oct (thin lines)")
    for k, testo, dy in ((k_b0, f"configuration\nk_b = {k_b0:.3g}", 25), (k_b, f"pole at {f_ris:.1f} Hz\n"
                                                                              f"k_b = {k_b:.3g}", 35)):
        ax.axvline(k, color=pt.INCHIOSTRO_SECONDARIO, lw=0.9, ls=":")
        ax.annotate(testo, (k, float(np.interp(np.log(k), np.log(sw[:, 0]), sw[:, 2]))), xytext=(8, dy),
                    textcoords="offset points", fontsize=7, color=pt.INCHIOSTRO_SECONDARIO)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("k_b, ring - telescope [N m/rad]")
    ax.set_ylabel("frequency [Hz]")
    ax.legend(loc="upper left", fontsize=7)
    pt._titolo(fig, "Telescope on the ring joint, elevation encoder channel: pole and zero against k_b "
                    "(I_yy fixed at CAD)", segnaposto)
    return fig


def main() -> None:
    args = argomenti()
    pt.applica_stile()
    par0 = pr.parametri_da_config(args.config)
    cfg = par0.sperimentale
    assert cfg is not None and par0.polare.mozzo is not None
    seg = par0.segnaposto
    out = RADICE / CARTELLA
    for d in ("figures", "tables"):
        (out / d).mkdir(parents=True, exist_ok=True)

    # --- misure ----------------------------------------------------------------------------------
    pol = [misura(cfg, fn, "polare", et) for fn, et in MISURE_POLARI]
    az = [misura(cfg, fn, "azimut", et) for fn, et in MISURE_AZIMUT]
    rif = {a: misura(cfg, fn, a, et) for a, (fn, et) in RIFERIMENTO.items()}
    righe = []
    for m in pol + az + list(rif.values()):
        ok = (m.f >= BANDA_ESTREMI[0]) & (m.f <= BANDA_ESTREMI[1])
        for tipo, fm, g in spm.estremi(m.f[ok], m.modulo_db[ok], prominenza_db=3.0, finestra_dec=0.05):
            righe.append([m.file, m.etichetta, tipo, f"{fm:.2f}", f"{g:.1f}"])
    with open(out / "tables" / "measured_features.csv", "w", newline="") as fh:
        csv.writer(fh).writerows([["file", "label", "kind", "f_Hz", "magnitude_dB"]] + righe)
    f_ris, g_ris = picco(pol[0], BANDA_RISONANZA)
    # antirisonanza subito sotto: minimo del modulo fra 40 Hz e la risonanza
    ok = (pol[0].f >= 40.0) & (pol[0].f < f_ris)
    f_antires = float(pol[0].f[ok][np.argmin(pol[0].modulo_db[ok])])
    picchi = {m.etichetta: picco(m, BANDA_RISONANZA)[0] for m in pol}

    # --- k_b per il polo a f_ris, I_yy fissa ------------------------------------------------------------
    dati = pr.carica_config(args.config)
    k_b0 = par0.polare.mozzo.k_b
    k_min, k_max, n_k = args.k_sweep
    ks = np.unique(np.r_[np.geomspace(k_min, k_max, int(n_k)), k_b0])
    sw = np.array([(k, *coppia_telescopio(con_kb(args.config, dati, k))) for k in ks])
    k_b = float(np.exp(brentq(lambda lk: coppia_telescopio(con_kb(args.config, dati, np.exp(lk)))[0] - f_ris,
                              np.log(k_min), np.log(k_max), xtol=1e-6)))
    par1 = con_kb(args.config, dati, k_b)
    s0, s1 = sistema(par0), sistema(par1)
    I_yy = float(s0.payload.I_yy)
    fp0, fz0 = coppia_telescopio(par0)
    fp1, fz1 = coppia_telescopio(par1)
    # inerzia dal lato motore che il modello mette contro il telescopio: (fp/fz)^2 = 1 + I_yy / J_lato
    J_lato = I_yy / ((fp1 / fz1) ** 2 - 1.0)
    with open(out / "tables" / "k_b_sweep.csv", "w", newline="") as fh:
        csv.writer(fh).writerows([["k_b_Nm_rad", "f_pole_Hz", "f_zero_Hz"]]
                                 + [[f"{k:.6g}", f"{a:.3f}", f"{b:.3f}"] for k, a, b in sw])

    # --- Bode ----------------------------------------------------------------------------------------------
    f = np.geomspace(1.0, 600.0, 4000)
    G_pol = [(f"model H7, k_b = {k_b0:.3g} N m/rad (configuration)", frf(par0, "polare", f, cfg.Kt)),
             (f"model H7, k_b = {k_b:.3g} N m/rad (pole at {f_ris:.1f} Hz)", frf(par1, "polare", f, cfg.Kt))]
    G_az = [(f"model H7 (configuration, k_psi,3 = {par0.azimut.attacco_telescopio.k:.3g} N m/rad)",
             frf(par0, "azimut", f, cfg.Kt))]
    fit = {}
    for m in pol + [rif["polare"]]:
        fit[m.etichetta] = {nome: spm.fit_bande(m.f, m.G, frf(p, "polare", m.f, cfg.Kt))
                            for nome, p in ((f"k_b = {k_b0:.3g}", par0), (f"k_b = {k_b:.3g}", par1))}
    for m in az + [rif["azimut"]]:
        fit[m.etichetta] = {"configuration": spm.fit_bande(m.f, m.G, frf(par0, "azimut", m.f, cfg.Kt))}
    figure = pt.salva(figura_bode(f, G_pol, pol + [rif["polare"]], "Elevation plant, encoder velocity / current: "
                                  "28 Oct 2026, lens at 12 mm", (30.0, 150.0), segnaposto=seg,
                                  marcatori=(f_ris,)), out / "figures", "bode_polar_20261028")
    figure += pt.salva(figura_bode(f, G_az, az + [rif["azimut"]], "Azimuth plant, encoder velocity / current: "
                                   "28 Oct 2026", (30.0, 150.0), segnaposto=seg), out / "figures",
                       "bode_azimuth_20261028")
    figure += pt.salva(figura_bode(f, [], pol + az, "Measured plants of 28 Oct 2026 (elevation with the lens at "
                                   "12 mm, azimuth)", (30.0, 150.0), segnaposto=seg), out / "figures",
                       "bode_measured_20261028")
    figure += pt.salva(figura_scansione(sw, k_b0, k_b, f_ris, f_antires, segnaposto=seg), out / "figures",
                       "k_b_sweep")

    # --- riepilogo -------------------------------------------------------------------------------------------
    sommario = {
        "config": args.config.name, "I_yy_kgm2": I_yy,
        "measured": {"resonance_Hz": f_ris, "resonance_dB": g_ris, "antiresonance_below_Hz": f_antires,
                     "resonance_by_file_Hz": picchi},
        "configuration": {"k_b_Nm_rad": k_b0, "pole_Hz": fp0, "zero_Hz": fz0},
        "pole_at_measured_resonance": {"k_b_Nm_rad": k_b, "pole_Hz": fp1, "zero_Hz": fz1,
                                       "J_motor_side_kgm2": J_lato,
                                       "k_b_for_zero_at_antiresonance_Nm_rad": I_yy * (2 * np.pi * f_antires) ** 2},
        "fit_percent": fit,
    }
    (out / "summary.json").write_text(json.dumps(sommario, indent=2))
    md = [f"# Misure del 28/10/2026 e rigidezza anello-telescopio ({args.config.name})", "",
          "Polare con la lente a 12 mm: " + ", ".join(f"`{Path(fn).name}`" for fn, _ in MISURE_POLARI)
          + f"; azimut: `{Path(MISURE_AZIMUT[0][0]).name}`. Riferimento: VTesi del 23/10.", "",
          f"- risonanza polare misurata (griglia all'1 %): **{f_ris:.2f} Hz** ({g_ris:.1f} dB); sulle altre griglie "
          + ", ".join(f"{v:.2f}" for v in list(picchi.values())[1:]) + " Hz",
          f"- antirisonanza subito sotto: {f_antires:.1f} Hz",
          f"- I_yy telescopio (CAD, fissa) = {I_yy:.4f} kg m^2",
          f"- configurazione: k_b = {k_b0:.4g} N m/rad -> zero {fz0:.1f} Hz, polo {fp0:.1f} Hz",
          f"- **k_b per il polo a {f_ris:.2f} Hz: {k_b:.4g} N m/rad** -> zero {fz1:.2f} Hz",
          f"- inerzia dal lato motore contro il telescopio: J = I_yy / ((f_p/f_z)^2 - 1) = {J_lato:.4f} kg m^2",
          f"- k_b che metterebbe invece lo zero a {f_antires:.1f} Hz: I_yy (2 pi f)^2 = "
          f"{I_yy * (2 * np.pi * f_antires) ** 2:.4g} N m/rad", "",
          "Fit NRMSE (%), bande `all`, `5-100 Hz`:", ""]
    md += [f"- {et}: " + "; ".join(f"{n} {v['all']:.1f} / {v['5-100 Hz']:.1f}" for n, v in d.items())
           for et, d in fit.items()]
    md += ["", "Figure: " + ", ".join(f"`figures/{p.name}`" for p in figure if p.suffix == ".png")]
    (out / "summary.md").write_text("\n".join(md) + "\n")
    print(json.dumps({k: sommario[k] for k in ("I_yy_kgm2", "measured", "configuration",
                                              "pole_at_measured_resonance")}, indent=2))
    for et, d in fit.items():
        print("fit", et, {n: {b: round(v, 1) for b, v in x.items()} for n, x in d.items()})


if __name__ == "__main__":
    main()
