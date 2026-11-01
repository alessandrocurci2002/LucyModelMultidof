"""Telescopio appeso al load di azimut (H7, richiesta dell'utente 2026-10-25).

La risonanza misurata a circa 50 Hz sull'azimut dipende dal numero di viti anello-telescopio. In H7 il
load di azimut non e' piu' un corpo unico: (forcella + mozzo + anello) ruota con il nodo L e il telescopio
e' un nodo P appeso a L dalla molla ``azimuth.telescope_attachment.k_Nm_rad``. Lo script confronta H7 con
la configurazione di base (H6) sul canale coppia -> encoder di azimut:

* poli e zeri strutturali del modello concentrato linearizzato (theta0 = 0, Omega0 = 0), senza le coppie
  polo-zero che si cancellano (modi polari visti dall'azimut per i prodotti d'inerzia);
* pianta velocita'/corrente contro le misure di azimut della configurazione, con il vincolo a terra;
* scansione di k: dove cade la coppia polo-zero del telescopio al variare della rigidezza.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_nodo_azimut.py --config config_hub_telescopeonring_alhub_fit_ground_aztel.yaml
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))

import matplotlib.pyplot as plt  # noqa: E402

from gimbal_mdof import concentrato as cc  # noqa: E402
from gimbal_mdof import fem  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import payload as plm  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402
from gimbal_mdof import tf  # noqa: E402
from gimbal_mdof import vincolo_terra as vt  # noqa: E402

TOLLERANZA_CANCELLAZIONE = 1e-3
"""Polo e zero a distanza relativa minore si cancellano (modo dell'altro asse visto per accoppiamento)."""
BANDA_COPPIA = (20.0, 400.0)
"""Banda in cui cercare la coppia polo-zero del telescopio nella scansione di k, Hz."""
MISURE_AZIMUT = ("VClosedLoop/VTesi/BodeMeasurement_Azimuth_14viti.csv", "VClosedLoop/BodeMeasurement_azimuth.csv",
                 "Azimuth_BodeMeasurement_Po_0.csv")
"""Misure di azimut a theta = 0 e 3.5 A, dalla piu' recente; devono stare in ``experimental.measurements``."""
MISURA_POLARE = "VClosedLoop/VTesi/BodeMeasurement_Polar_14viti.csv"
"""Misura polare presa insieme a quella di azimut a 14 viti (23/10/2026)."""
COLORI_MISURE = (pt.INCHIOSTRO, pt.INCHIOSTRO_TENUE, pt.SERIE[2])
TRATTI_MISURE = ("-", "-.", ":")


def argomenti() -> argparse.Namespace:
    a = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    a.add_argument("--config", type=Path, default=RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    a.add_argument("--base", type=Path, default=RADICE / "config_hub_telescopeonring_alhub_fit_ground.yaml",
                   help="configurazione di confronto senza il nodo del telescopio sull'azimut")
    a.add_argument("--k-sweep", type=float, nargs=3, default=[2.0e4, 4.0e5, 41], metavar=("K_MIN", "K_MAX", "N"),
                   help="scansione logaritmica di k [N m/rad]")
    return a.parse_args()


def sistema(par: pr.Parametri):
    """Modello concentrato con il telescopio (come in :func:`gimbal_mdof.sperimentale.sistemi_per_misure`)."""
    assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
            fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    return cc.costruisci_sistema_concentrato(par, plm.costruisci_payload(par, assi[1]), assi)


def radici_azimut(ml: lz.ModelloLineare, f_max: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """``(poli, zeri)`` del canale tau_psi -> psi_e come ``(f [Hz], zeta)``, poco smorzati, fino a ``f_max``,
    senza le coppie che si cancellano."""
    z = tf.zeri_invarianti(ml.A, ml.Bu[:, 0], ml.Cy[ml.riga("psi_e")])
    sel = lambda r: r[(np.imag(r) > 0) & (np.imag(r) > np.abs(np.real(r))) & (np.abs(r) <= 2 * np.pi * f_max)]  # noqa: E731
    p, z = sorted(sel(ml.poli), key=abs), sorted(sel(z), key=abs)
    vivi_p = [r for r in p if not any(abs(r - s) < TOLLERANZA_CANCELLAZIONE * abs(r) for s in z)]
    vivi_z = [r for r in z if not any(abs(r - s) < TOLLERANZA_CANCELLAZIONE * abs(r) for s in p)]
    fz = lambda r: (float(abs(r) / (2 * np.pi)), float(-np.real(r) / abs(r)))  # noqa: E731
    return [fz(r) for r in vivi_p], [fz(r) for r in vivi_z]


def forma_modo(s, ml: lz.ModelloLineare, f: float) -> dict[str, float]:
    """Ampiezze di azimut del modo piu' vicino a ``f`` (Hz), normalizzate sul telescopio (o su L)."""
    w, V = np.linalg.eig(ml.A)
    i = int(np.argmin(np.abs(w - 1j * 2 * np.pi * f)))
    n = ml.n
    st = s.accoppiato.struttura
    v = V[:n, i]
    v = v / v[st.j_psi_p]
    return {nome: float(np.real(v[j])) for nome, j in zip(st.nomi[: st.n_psi], range(st.n_psi))}


def coppia_telescopio(ml: lz.ModelloLineare) -> tuple[float, float]:
    """Polo e zero del canale encoder di azimut nella :data:`BANDA_COPPIA` (``nan`` se assenti)."""
    poli, zeri = radici_azimut(ml, BANDA_COPPIA[1])
    fp = [f for f, _ in poli if f >= BANDA_COPPIA[0]]
    fz = [f for f, _ in zeri if f >= BANDA_COPPIA[0]]
    return (fp[0] if fp else float("nan")), (fz[0] if fz else float("nan"))


def scansione_k(percorso: Path, ks: np.ndarray) -> np.ndarray:
    dati = pr.carica_config(percorso)
    out = []
    for k in ks:
        dati["azimuth"]["telescope_attachment"]["k_Nm_rad"] = float(k)
        par = pr.parametri_da_config(percorso, dati=dati)
        out.append((k, *coppia_telescopio(lz.linearizza(sistema(par).accoppiato, 0.0, 0.0))))
    return np.array(out)


def figura_bode(f: np.ndarray, modelli: list[tuple[str, np.ndarray]], misure: list, *, segnaposto: bool):
    fig, assi = plt.subplots(2, 2, figsize=(12.0, 7.0), sharex="col", gridspec_kw={"width_ratios": [1.6, 1.0]})
    for col, f_lim in ((0, (1.0, 2000.0)), (1, (20.0, 150.0))):
        ax_m, ax_f = assi[0, col], assi[1, col]
        for k, m in enumerate(misure):
            st = {"color": COLORI_MISURE[k], "marker": "o", "ms": 2.0, "lw": 0.8, "alpha": 0.85}
            ax_m.plot(m.f, m.modulo_db, label=f"measured: {m.etichetta}", **st)
            ax_f.plot(m.f, m.fase_deg, **st)
        for k, (nome, G) in enumerate(modelli):
            st = pt._stile(k)
            ax_m.plot(f, pt._db(G), lw=1.6, label=nome, **st)
            ax_f.plot(f, pt._fase_avvolta(G), lw=1.3, **st)
        for ax in (ax_m, ax_f):
            ax.set_xscale("log")
            ax.set_xlim(*f_lim)
        ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
        ax_f.set_ylabel("phase [deg]")
        ax_f.set_ylim(-180, 180)
        ax_f.set_yticks(np.arange(-180, 181, 90))
        ax_f.set_xlabel("frequency [Hz]")
    assi[0, 1].set_ylim(-80, -20)
    assi[0, 0].legend(loc="lower right", fontsize=7)
    assi[0, 1].set_title("zoom on the telescope pair", fontsize=9)
    pt._titolo(fig, "Azimuth plant, encoder velocity / current: telescope attached to the azimuth load", segnaposto)
    return fig


def figura_scansione(sw: np.ndarray, k_cfg: float, misurati: dict[str, list[tuple[str, float]]], *,
                     segnaposto: bool):
    """Polo e zero del modello contro k; le caratteristiche misurate sono linee orizzontali (blu le
    risonanze, arancio le antirisonanze, un tratto per misura)."""
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(sw[:, 0], sw[:, 1], color=pt.SERIE[0], lw=2.0, label="model pole (resonance)")
    ax.plot(sw[:, 0], sw[:, 2], color=pt.SERIE[1], lw=2.0, ls="--", label="model zero (antiresonance)")
    colore = {"risonanza": pt.SERIE[0], "antirisonanza": pt.SERIE[1]}
    for (etichetta, car), ls in zip(misurati.items(), TRATTI_MISURE):
        for tipo, fm in car:
            ax.axhline(fm, color=colore[tipo], lw=0.9, ls=ls, alpha=0.55)
        ax.plot([], [], color=pt.INCHIOSTRO_SECONDARIO, lw=0.9, ls=ls, label=f"measured: {etichetta}")
    ax.axvline(k_cfg, color=pt.INCHIOSTRO_SECONDARIO, lw=1.0, ls=":")
    ax.annotate(f"k = {k_cfg:.3g} N m/rad\n(polar k_b, 12 screws)", (k_cfg, float(np.interp(k_cfg, sw[:, 0], sw[:, 1]))),
                xytext=(-110, 25), textcoords="offset points", fontsize=7, color=pt.INCHIOSTRO_SECONDARIO,
                arrowprops={"arrowstyle": "->", "color": pt.INCHIOSTRO_SECONDARIO, "lw": 0.8})
    ax.set_xscale("log")
    ax.set_xlabel("k, azimuth load - telescope [N m/rad]")
    ax.set_ylabel("frequency [Hz]")
    ax.legend(loc="upper left", fontsize=7, title="measured lines: blue resonance, orange antiresonance",
              title_fontsize=7)
    pt._titolo(fig, "Telescope pole-zero pair on the azimuth encoder channel against k (theta0 = 0)", segnaposto)
    return fig


def main() -> None:
    args = argomenti()
    pt.applica_stile()
    par = pr.parametri_da_config(args.config)
    base = pr.parametri_da_config(args.base)
    at = par.azimut.attacco_telescopio
    if at is None:
        raise SystemExit(f"{args.config.name}: manca azimuth.telescope_attachment (enabled: true)")
    out = RADICE / par.cartella_risultati
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    cfg = par.sperimentale
    assert cfg is not None
    seg = par.segnaposto

    # --- poli e zeri -----------------------------------------------------------------------------
    s7, s6 = sistema(par), sistema(base)
    ml7, ml6 = (lz.linearizza(s.accoppiato, 0.0, 0.0) for s in (s7, s6))
    righe_pz = []
    radici = {}
    for nome, ml in (("H6", ml6), ("H7", ml7)):
        poli, zeri = radici_azimut(ml, cfg.f_max)
        radici[nome] = {"poles": poli, "zeros": zeri}
        righe_pz += [[nome, "pole", f"{fr:.3f}", f"{z:.4f}"] for fr, z in poli]
        righe_pz += [[nome, "zero", f"{fr:.3f}", f"{z:.4f}"] for fr, z in zeri]
    with open(out / "tables" / "azimuth_poles_zeros.csv", "w", newline="") as fh:
        csv.writer(fh).writerows([["model", "kind", "f_Hz", "zeta"]] + righe_pz)
    fp7, fz7 = coppia_telescopio(ml7)
    modo = forma_modo(s7, ml7, fp7)

    # --- misure di azimut ---------------------------------------------------------------------------
    voci = {v.file: v for v in cfg.misure}
    misure = [spm.carica_misura(cfg, voci[nome]) for nome in MISURE_AZIMUT]
    caratteristiche = {}
    for m in misure:
        ok = (m.f >= BANDA_COPPIA[0]) & (m.f <= 150.0)
        car = [(t, fm) for t, fm, _ in spm.estremi(m.f[ok], m.modulo_db[ok], prominenza_db=4.0, finestra_dec=0.15)]
        caratteristiche[m.etichetta] = car

    # --- Bode ------------------------------------------------------------------------------------
    f = np.logspace(0, np.log10(2000.0), 3000)
    G = {nome: spm.frf_velocita(vt.applica(ml, p), "azimut", f, cfg.Kt)["psi_e"]
         for nome, ml, p in (("H6", ml6, base), ("H7", ml7, par))}
    modelli = [(f"H6: telescope rigid on the fork", G["H6"]),
               (f"H7: telescope on k = {at.k:.3g} N m/rad", G["H7"])]
    fit = {}
    for m in misure:
        fit[m.etichetta] = {nome: spm.fit_bande(m.f, m.G, spm.frf_velocita(vt.applica(ml, p), "azimut", m.f, cfg.Kt)
                                                ["psi_e"]) for nome, ml, p in (("H6", ml6, base), ("H7", ml7, par))}
    figure = pt.salva(figura_bode(f, modelli, misure, segnaposto=seg), out / "figures", "bode_azimuth_H6_H7")
    # polare: misura presa insieme a quella di azimut a 14 viti contro il modello della configurazione
    voce_pol = voci[MISURA_POLARE]
    conf_pol = spm.confronta(par, spm.carica_misura(cfg, voce_pol), {True: s7})
    figure += pt.salva(pt.bode_misura(conf_pol, segnaposto=seg, f_lim=(1.0, 2000.0)), out / "figures",
                       "bode_polar_14screws_H7")

    # --- scansione di k ------------------------------------------------------------------------------
    k_min, k_max, n_k = args.k_sweep
    ks = np.unique(np.r_[np.geomspace(k_min, k_max, int(n_k)), at.k])
    sw = scansione_k(args.config, ks)
    with open(out / "tables" / "k_sweep.csv", "w", newline="") as fh:
        csv.writer(fh).writerows([["k_Nm_rad", "f_pole_Hz", "f_zero_Hz"]]
                                 + [[f"{k:.6g}", f"{a:.3f}", f"{b:.3f}"] for k, a, b in sw])
    figure += pt.salva(figura_scansione(sw, at.k, caratteristiche, segnaposto=seg), out / "figures", "k_sweep")

    # --- riepilogo -------------------------------------------------------------------------------------
    pl = s7.payload
    J_L = float(s7.rid_psi.M_r[s7.rid_psi.i_estremo, s7.rid_psi.i_estremo]) + pl.J_ax
    sommario = {
        "config": args.config.name, "base": args.base.name, "k_Nm_rad": at.k,
        "J_L_kgm2": J_L, "J_telescope_psi_kgm2": float(pl.J_psi(0.0)),
        "k_psi2_Nm_rad": float(s7.rid_psi.k_catena[1]),
        "telescope_pair_Hz": {"pole": fp7, "zero": fz7}, "mode_shape_psi": modo,
        "roots": radici, "measured_features_20_150Hz": caratteristiche, "fit_percent": fit,
        "beta_psi": {"H6": s6.beta_psi, "H7": s7.beta_psi},
        "k_sweep": [[float(k), float(a), float(b)] for k, a, b in sw],
    }
    (out / "summary.json").write_text(json.dumps(sommario, indent=2))
    sep = lambda a, b: 100 * (b - a) / a  # noqa: E731
    md = [f"# Telescopio appeso al load di azimut ({args.config.name})", "",
          f"Modello concentrato linearizzato a theta0 = 0, canale coppia -> encoder di azimut. Base: {args.base.name}.",
          "", f"- k (L - telescopio) = {at.k:.4g} N m/rad ({at.fonte})",
          f"- J_L (forcella + alberi polari + mozzo + anello) = {J_L:.4f} kg m^2, J_telescopio(theta0 = 0) = "
          f"{pl.J_psi(0.0):.4f} kg m^2, k_psi,2 = {s7.rid_psi.k_catena[1]:.5g} N m/rad",
          f"- **coppia del telescopio: polo {fp7:.2f} Hz, zero {fz7:.2f} Hz** (zero sopra il polo del "
          f"{sep(fp7, fz7):.1f} %)",
          "- forma del modo (azimut, telescopio = 1): " + ", ".join(f"{k} {v:+.3f}" for k, v in modo.items()),
          f"- beta di azimut (zeta sul primo modo flessibile, ora il modo del telescopio): H6 {s6.beta_psi:.3g}, "
          f"H7 {s7.beta_psi:.3g} s: gli smorzatori di azimut crescono di {s7.beta_psi / s6.beta_psi:.1f} volte",
          "", "| modello | poli [Hz] | zeri [Hz] |", "|---|---|---|"]
    md += [f"| {n} | {', '.join(f'{x:.2f}' for x, _ in r['poles'] if x > 0.5)} | "
           f"{', '.join(f'{x:.2f}' for x, _ in r['zeros'])} |" for n, r in radici.items()]
    md += ["", "Misure (20-150 Hz, prominenza 4 dB):", ""]
    md += [f"- {et}: " + ", ".join(f"{t} {fm:.1f} Hz" for t, fm in car) for et, car in caratteristiche.items()]
    i_cfg = int(np.argmin(np.abs(sw[:, 0] - at.k)))
    md += ["", f"Scansione di k ({ks[0]:.3g}-{ks[-1]:.3g} N m/rad, `tables/k_sweep.csv`): lo scarto fra zero e polo "
           f"resta fra {sep(sw[:, 1], sw[:, 2]).min():.1f} e {sep(sw[:, 1], sw[:, 2]).max():.1f} % "
           f"({sep(sw[i_cfg, 1], sw[i_cfg, 2]):.1f} % a k = {at.k:.3g}).", "",
           "Figure: `figures/bode_azimuth_H6_H7.pdf`, `figures/k_sweep.pdf`."]
    (out / "summary.md").write_text("\n".join(md) + "\n")
    print(json.dumps({k: sommario[k] for k in ("k_Nm_rad", "J_L_kgm2", "J_telescope_psi_kgm2", "k_psi2_Nm_rad",
                                              "telescope_pair_Hz", "mode_shape_psi", "measured_features_20_150Hz",
                                              "beta_psi")}, indent=2))
    for nome in ("H6", "H7"):
        print(nome, "poli", [round(x, 2) for x, _ in radici[nome]["poles"]],
              "zeri", [round(x, 2) for x, _ in radici[nome]["zeros"]])
    for et, d in fit.items():
        print("fit", et, {n: {b: round(v, 1) for b, v in x.items()} for n, x in d.items()})
    print("figure:", *[str(p.relative_to(RADICE)) for p in figure])


if __name__ == "__main__":
    main()
