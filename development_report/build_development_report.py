"""Builds the development report: figures, tables and numbers from the two model versions, then compiles it.

    .venv/bin/python development_report/build_development_report.py [--no-latex]

The two versions compared in the report are kept untouched (traceability):

    7-DOF  config.yaml                                          -> results/
    8-DOF  config_hub_telescopeonring_alhub_stiffring.yaml      -> results_hub_telescopeonring_alhub_stiffring/

Everything that depends on the models is recomputed here from the configuration files; the text of the report
(``development_report.tex``) only reads the generated macros (``generated/numbers.tex``), tables
(``generated/*.tex``) and figures (``figures/``).
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
RADICE = AQUI.parent
sys.path.insert(0, str(RADICE / "src"))

from gimbal_mdof import concentrato as cc  # noqa: E402
from gimbal_mdof import concentrato_analisi as ca  # noqa: E402
from gimbal_mdof import fem  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modal as md  # noqa: E402
from gimbal_mdof import payload as plm  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402
from gimbal_mdof import vincolo_terra as vt  # noqa: E402

FIG, GEN = AQUI / "figures", AQUI / "generated"
CFG_7 = RADICE / "config.yaml"
CFG_8 = RADICE / "config_hub_telescopeonring_alhub_stiffring.yaml"
CFG_8_RING1 = RADICE / "config_hub_telescopeonring_alhub.yaml"   # same as CFG_8 with the ring alone (factor 1)
CFG_FIT = RADICE / "config_hub_telescopeonring_alhub_fit.yaml"  # H5: k_b and kappa_r refitted (19 Oct tests)
CFG_GROUND = RADICE / "config_hub_telescopeonring_alhub_fit_ground.yaml"  # H6: H5 + stiffness to ground, 12 screws
RES_7 = RADICE / "results"
RES_8 = RADICE / "results_hub_telescopeonring_alhub_stiffring"

F_MIN, F_MAX = 0.1, 2000.0
F_GRID = np.logspace(np.log10(F_MIN), np.log10(F_MAX), 2500)


# --------------------------------------------------------------------------- helpers


def sci(x: float, d: int = 3) -> str:
    """LaTeX scientific notation, ``a\\times10^{b}``."""
    if x == 0:
        return "0"
    if not np.isfinite(x):
        return r"\infty" if x > 0 else "--"
    m, e = f"{x:.{d - 1}e}".split("e")
    e = int(e)
    return rf"{m}\times10^{{{e}}}" if e else m


def num(x: float, d: int = 3) -> str:
    """``d`` significant figures, scientific notation for very large or very small values."""
    if x == 0:
        return "0"
    if abs(x) >= 1e4 or abs(x) < 1e-3:
        return sci(x, d)
    s = f"{x:.{d}g}"
    return s


def pct(a: float, b: float) -> str:
    """Relative difference of ``a`` from ``b`` in per cent, without a signed zero."""
    d = 100.0 * (a / b - 1.0)
    return "0.0" if abs(d) < 0.05 else f"{d:+.1f}"


def write_rows(path: Path, rows: list[str]) -> None:
    """Table body: rows joined by ``\\\\``, the last one without it (the report adds it before ``\\bottomrule``)."""
    clean = [r.rstrip().removesuffix("\\\\").rstrip() for r in rows]
    path.write_text(" \\\\\n".join(clean) + "\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as f:
        rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def build_models(cfg: Path):
    par = pr.parametri_da_config(cfg)
    return par, spm.sistemi_per_misure(par)


def pole_zero(s, ch: str, j: int) -> tuple[np.ndarray, np.ndarray]:
    ml = lz.linearizza(s.accoppiato, 0.0, 0.0)
    return spm.poli_zeri_modello(ml, ch, j, 2.0 * F_MAX)


def local_peak(f, g_db, lo, hi) -> tuple[float, float]:
    m = (f >= lo) & (f <= hi)
    k = int(np.argmax(g_db[m]))
    return float(f[m][k]), float(g_db[m][k])


# --------------------------------------------------------------------------- comparison figures


def bode_comparison(conf, model: str, fname: str, title: str) -> None:
    """Measured plant (black) against the model at the encoder (blue, measured channel) and motor (dashed)."""
    mis = conf.misura
    labels = {"psi_e": rf"{model} model, encoder $\dot\psi_e$", "psi_m": rf"{model} model, motor $\dot\psi_m$",
              "theta_e": rf"{model} model, encoder $\dot\theta_e$", "theta_m": rf"{model} model, motor $\dot\theta_m$"}
    fig, (ax_m, ax_f) = pt.plt.subplots(2, 1, figsize=(9.0, 6.2), sharex=True)
    for k, (c, G) in enumerate(conf.G_modello.items()):
        st = pt._stile(k)
        ax_m.plot(conf.f_modello, pt._db(G), label=labels[c], lw=1.5, **st)
        ax_f.plot(conf.f_modello, pt._fase_avvolta(G), lw=1.2, **st)
    ax_m.plot(mis.f, mis.modulo_db, color=pt.INCHIOSTRO, marker="o", ms=2.5, lw=1.0,
              label=f"measurement ({mis.eccitazione_A:g} A)")
    ax_f.plot(mis.f, mis.fase_deg, color=pt.INCHIOSTRO, marker="o", ms=2.5, lw=1.0)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(F_MIN, F_MAX)
        ax.axvspan(F_MIN, mis.f[0], color=pt.GRIGLIA_MINORE, zorder=0)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower right", fontsize=7)
    pt._titolo(fig, title, True)
    fig.tight_layout()
    pt.salva(fig, FIG, fname)
    pt.plt.close(fig)


# --------------------------------------------------------------------------- FEM (all 18 elements)


def full_fem(par):
    """FEM of both axes (all 18 elements, no reduction) coupled through the telescope: with and without it."""
    assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
            fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    pl = plm.costruisci_payload(par, assi[1])
    s_t = md.costruisci_sistema_modale(par, "consistent", payload=pl, assi=assi, valida=False)
    s_n = md.costruisci_sistema_modale(par, "consistent", payload=spm.payload_senza_telescopio(pl),
                                       assi=spm.assi_senza_telescopio(par, assi), valida=False)
    return s_t, s_n


def fem_frf_figure(s_t, s_n, misure: dict, sel: dict, Kt: float, fname: str, title: str) -> dict:
    """FEM (all DOFs) in the units of the measurement, (rad/s)/A, at encoder and motor; measurement in black.
    Rows: azimuth, polar, polar without telescope; columns: magnitude, phase. Returns poles/zeros."""
    casi = (("az", "azimut", s_t, "azimuth (polar at 0 deg)"), ("pol", "polare", s_t, "polar, with telescope"),
            ("nt", "polare", s_n, "polar, without telescope"))
    fig, axs = pt.plt.subplots(3, 2, figsize=(10.0, 9.6), sharex=True)
    out = {}
    for r, (k, asse, s, nome) in enumerate(casi):
        mis = misure[sel[k]]
        ml = s.lineare_completo(mis.theta0, 0.0)
        G = spm.frf_velocita(ml, asse, F_GRID, Kt)
        ce, cm = spm.canali(asse)
        j = 0 if asse == "azimut" else 1
        p, z = spm.poli_zeri_modello(ml, ce, j, 2.0 * F_MAX)
        out[k] = (p[p > 1.0], z[z > 1.0])
        sym = r"\psi" if asse == "azimut" else r"\theta"
        for i, (c, lab) in enumerate(((ce, rf"FEM, encoder $\dot{sym}_e$"), (cm, rf"FEM, motor $\dot{sym}_m$"))):
            st = pt._stile(i)
            axs[r, 0].plot(F_GRID, pt._db(G[c]), lw=1.4, label=lab, **st)
            axs[r, 1].plot(F_GRID, pt._fase_avvolta(G[c]), lw=1.1, **st)
        axs[r, 0].plot(mis.f, mis.modulo_db, color=pt.INCHIOSTRO, marker="o", ms=1.6, lw=0.8,
                       label=f"measurement ({mis.eccitazione_A:g} A)")
        axs[r, 1].plot(mis.f, mis.fase_deg, color=pt.INCHIOSTRO, marker="o", ms=1.6, lw=0.8)
        for fp in p[(p > 1.0) & (p < F_MAX)]:
            for ax in axs[r]:
                ax.axvline(fp, color=pt.GRIGLIA, lw=0.9, ls=":", zorder=0)
        axs[r, 0].set_title(nome, fontsize=9, loc="left")
        axs[r, 0].set_ylabel("magnitude [dB (rad/s)/A]")
        axs[r, 1].set_ylabel("phase [deg]")
        axs[r, 1].set_ylim(-180, 180)
        axs[r, 1].set_yticks(np.arange(-180, 181, 90))
        axs[r, 0].legend(loc="lower right", fontsize=6.5)
        for ax in axs[r]:
            ax.set_xscale("log")
            ax.set_xlim(F_MIN, F_MAX)
            ax.axvspan(F_MIN, mis.f[0], color=pt.GRIGLIA_MINORE, zorder=0)
    for ax in axs[-1]:
        ax.set_xlabel("frequency [Hz]")
    pt._titolo(fig, title, True)
    fig.tight_layout()
    pt.salva(fig, FIG, fname)
    pt.plt.close(fig)
    return out


def mode_shapes_figure(s_t, assi_modi: tuple[tuple[str, str], ...], n_modi: int, fname: str, title: str) -> None:
    """Flexible modes of the FEM of each axis (telescope frozen at theta_ref), shape against the node
    position. A branch node (telescope P hanging from the hub H) is drawn apart, joined to H by a dashed line."""
    fig, axs = pt.plt.subplots(len(assi_modi), n_modi, figsize=(3.4 * n_modi, 2.5 * len(assi_modi)), squeeze=False)
    for r, (asse, nome) in enumerate(assi_modi):
        an, ast = s_t.analisi(asse), s_t.asse(asse)
        mod = ast.modello
        x = np.asarray(mod.x) * 1e3
        ramo = bool(getattr(mod, "ramo", False))
        catena = [i for i in range(len(x)) if not (ramo and i == mod.i_P)]
        for c in range(n_modi):
            ax = axs[r, c]
            k = c + 1
            phi = an.Phi[:, k] / np.max(np.abs(an.Phi[:, k]))
            ax.plot(x[catena], phi[catena], marker="o", ms=3, lw=1.4, color=pt.SERIE[0])
            ax.axhline(0.0, color=pt.GRIGLIA, lw=0.8)
            for i, mk, lab in ((ast.i_e, "s", "encoder"), (ast.i_m, "D", "motor")):
                ax.plot(x[i], phi[i], mk, ms=7, mfc="none", color=pt.INCHIOSTRO, label=lab)
            if ramo:
                iH, iP = mod.i_H, mod.i_P
                ax.plot([x[iH], x[iP]], [phi[iH], phi[iP]], ls="--", color=pt.SERIE[1], lw=1.0)
                ax.plot(x[iP], phi[iP], "*", ms=11, color=pt.SERIE[1], label="telescope (on $k_b$)")
            else:
                ax.plot(x[ast.i_agg], phi[ast.i_agg], "^", ms=7, mfc="none", color=pt.INCHIOSTRO,
                        label="fork L" if asse == "azimut" else "telescope P")
            ax.set_title(f"{nome}, mode {k}: {an.f[k]:.4g} Hz", fontsize=9)
            ax.set_ylim(-1.15, 1.15)
            ax.set_xlabel("position along the shaft [mm]")
            if c == 0:
                ax.set_ylabel("normalised shape [-]")
                ax.legend(loc="best", fontsize=6.5)
    pt._titolo(fig, title, True)
    fig.tight_layout()
    pt.salva(fig, FIG, fname)
    pt.plt.close(fig)


# --------------------------------------------------------------------------- screw test (19 October 2026)


def refined_extremum(f, g, lo, hi, kind, refine: bool = True):
    """Extremum of |G| [dB] in [lo, hi], refined by a parabola in (log f, dB) through the grid point and its
    neighbours: the 19 October records have a 6 % frequency step, as large as the shifts being measured."""
    m = np.where((f >= lo) & (f <= hi))[0]
    i = int(m[np.argmin(g[m])] if kind == "min" else m[np.argmax(g[m])])
    if not refine or i in (m[0], m[-1]):
        return float(f[i])
    c = np.polyfit(np.log10(f[i - 1:i + 2]), g[i - 1:i + 2], 2)
    return float(10 ** (-c[1] / (2 * c[0]))) if c[0] != 0 else float(f[i])


def phase_down_cross(f, ph, lo, hi):
    """Frequency where the phase crosses 0 deg downwards in [lo, hi] (linear interpolation)."""
    m = np.where((f >= lo) & (f <= hi))[0]
    for a, b in zip(m[:-1], m[1:]):
        if ph[a] > 0 >= ph[b]:
            return float(f[a] + (f[b] - f[a]) * ph[a] / (ph[a] - ph[b]))
    return float("nan")


def load(cfg, fname: str, lab: str):
    from types import SimpleNamespace

    return spm.carica_misura(cfg, SimpleNamespace(file=fname, asse="polare", theta0=0.0, telescopio=True,
                                                  etichetta=lab))

# Polar measurements with a different number of screws in the ring--telescope connection. The 4-screw
# configuration is the one of all the earlier polar measurements (16 October included); the 8-screw test is
# still to be done.
SCREWS = (("Elevation_BodeMeasurement_bode_New_16_10_26.csv", "4 screws (16 Oct, 3.5 A)", 4, "k", "-"),
          ("Elevation_BodeMeasurement_4viti.csv", "4 screws (19 Oct, 3 A)", 4, "C0", "-"),
          ("Elevation_BodeMeasurement_2viti.csv", "2 screws (19 Oct, 3 A)", 2, "C1", "--"),
          ("Elevation_BodeMeasurement_0viti.csv", "0 screws (19 Oct, 3 A)", 0, "C3", "-."))


def screw_test(cfg, Jtot: float) -> dict:
    """Overlay of the screw-test measurements and their features; two-mass estimate of k_b (eq. twomass)."""
    window = refined_extremum

    colori = {"k": pt.INCHIOSTRO, "C0": pt.SERIE[0], "C1": pt.SERIE[1], "C3": pt.SERIE[2]}
    fig, (ax_m, ax_f) = pt.plt.subplots(2, 1, figsize=(9.0, 6.6), sharex=True)
    out, righe = {}, []
    for k, (fname, lab, n, c, ls) in enumerate(SCREWS):
        mis = load(cfg, fname, lab)
        f, g = mis.f, mis.modulo_db
        ax_m.plot(f, g, color=colori[c], ls=ls, marker="o", ms=2.0, lw=1.0, label=lab)
        ax_f.plot(f, mis.fase_deg, color=colori[c], ls=ls, marker="o", ms=2.0, lw=1.0)
        if n > 0:
            fz, fp = window(f, g, 12, 25, "min"), window(f, g, 30, 60, "max")
            r = (fp / fz) ** 2 - 1.0
            I = Jtot * r / (1.0 + r)
            kb = (2 * np.pi * fz) ** 2 * I
            feat = dict(fz=fz, fp=fp, ph0=phase_down_cross(f, mis.fase_deg, 30, 70),
                        p90=window(f, g, 75, 110, "max"), z150=window(f, g, 130, 180, "min"),
                        p200=window(f, g, 170, 240, "max"), ratio=fp / fz, I=I, Jsh=Jtot - I, kb=kb)
            righe.append(rf"{lab} & {fz:.2f} & {fp:.1f} & {feat['ph0']:.1f} & {fp / fz:.2f} & {I:.3f} & "
                         rf"{Jtot - I:.3f} & ${sci(kb, 3)}$ & {feat['p90']:.0f} & "
                         rf"{feat['z150']:.0f}\,/\,{feat['p200']:.0f} \\")
        else:
            feat = dict(n1=window(f, g, 25, 32, "min"), n2=window(f, g, 32, 40, "max"), p75=window(f, g, 65, 85, "max"),
                        z170=window(f, g, 150, 190, "min"), p225=window(f, g, 200, 250, "max"))
            righe.append(rf"{lab} & --$^{{a}}$ & -- & -- & -- & -- & -- & -- & "
                         rf"{feat['p75']:.0f}$^{{b}}$ & {feat['z170']:.0f}\,/\,{feat['p225']:.0f} \\")
        out[n if n not in out else 44] = feat
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(1.0, 600.0)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7)
    pt._titolo(fig, "Polar axis with telescope: measured plant against the number of ring-telescope screws", False)
    fig.tight_layout()
    pt.salva(fig, FIG, "screws_polar")
    pt.plt.close(fig)
    write_rows(GEN / "screws.tex", righe)
    return out


# --------------------------------------------------------------------------- 8 and 12 screws (19 October 2026)

# folder LucyIdentification/VAntiRes (README.txt): number of screws in the ring--telescope connection and
# radial gap between ring and tube measured with a feeler gauge (bottom / top of the telescope)
HIGH = (("VAntiRes/BodeMeasurement_4viti.csv", "4 screws", 4, "not recorded"),
        ("VAntiRes/BodeMeasurement_8viti_simmetrico.csv", "8 screws, symmetric", 8, "$<0.1$--0.2 / 0.2--0.4$^{c}$"),
        ("VAntiRes/BodeMeasurement_8viti_1inferiore_2superiore.csv", "8 screws, 0.1/0.2", 8, "0.1 / 0.2"),
        ("VAntiRes/BodeMeasurement_12viti_dueDecimiOmogeneo.csv", "12 screws, 0.2 uniform", 12, "0.2 all round"),
        ("VAntiRes/BodeMeasurement_12viti_1giu_2su.csv", "12 screws, 0.1/0.2", 12, "0.1 / 0.2"),
        ("VAntiRes/BodeMeasurement_12viti_2giu_2su.csv", "12 screws, 0.2/0.2", 12, "0.2 / 0.2"))


def high_screws(cfg, IJJ: float, kg: float, low: dict) -> dict:
    """Overlay of the 4-, 8- and 12-screw records and the table of the whole connection test: antiresonance of
    the telescope pair, next resonance, k_b = (2 pi f_a)^2 I_JJ - k_g with the model I_JJ, features near 400 Hz."""
    stili = [(pt.INCHIOSTRO, "-"), (pt.SERIE[0], "-"), (pt.SERIE[0], "--"), (pt.SERIE[1], "-"), (pt.SERIE[1], "--"),
             (pt.SERIE[2], "-.")]
    fig, (ax_m, ax_f) = pt.plt.subplots(2, 1, figsize=(9.0, 6.6), sharex=True)
    righe, out = [], {}
    kb = lambda fa: (2 * np.pi * fa) ** 2 * IJJ - kg  # noqa: E731
    for fname, lbl, n in (("VAntiRes/BodeMeasurement_2viti.csv", "2 screws", 2),
                          ("VAntiRes/BodeMeasurement_4viti.csv", "4 screws", 4)):
        mis = load(cfg, fname, lbl)
        fz = refined_extremum(mis.f, mis.modulo_db, 12, 25, "min")
        fp = refined_extremum(mis.f, mis.modulo_db, 30, 60, "max")
        if mis.f[-1] > 600:
            quattro = (rf"{refined_extremum(mis.f, mis.modulo_db, 375, 420, 'min', False):.0f}\,/\,"
                       rf"{refined_extremum(mis.f, mis.modulo_db, 395, 470, 'max', False):.0f}")
        else:
            quattro = "--$^{d}$"
        righe.append(rf"{lbl} & {n} & not recorded & {fz:.1f} & {fp:.1f} & ${sci(kb(fz), 3)}$ & {quattro} \\")
    for (fname, lab, n, gap), (c, ls) in zip(HIGH, stili):
        mis = load(cfg, fname, lab)
        f, g = mis.f, mis.modulo_db
        ax_m.plot(f, g, color=c, ls=ls, marker="o", ms=1.8, lw=1.0, label=lab)
        ax_f.plot(f, mis.fase_deg, color=c, ls=ls, marker="o", ms=1.8, lw=1.0)
        if n == 4:
            continue
        fa = refined_extremum(f, g, 40, 90, "min")
        fpk = refined_extremum(f, g, fa * 1.05, 130, "max")
        # grid values: near 400 Hz the pair only marks the presence of the mode, and a parabola on a 6 % step
        # can move the peak by more than the step
        z4, p4 = refined_extremum(f, g, 375, 420, "min", False), refined_extremum(f, g, 395, 470, "max", False)
        out[lab] = dict(fa=fa, fpk=fpk, kb=kb(fa), z4=z4, p4=p4)
        nota = "$^{e}$" if "0.2/0.2" in lab else ""
        righe.append(rf"{lab} & {n} & {gap} & {fa:.1f}{nota} & {fpk:.0f} & ${sci(kb(fa), 3)}$ & "
                     rf"{z4:.0f}\,/\,{p4:.0f} \\")
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(1.0, 2000.0)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7, ncol=2)
    pt._titolo(fig, "Polar axis with telescope: 4, 8 and 12 ring-telescope screws (19 Oct, 3 A)", False)
    fig.tight_layout()
    pt.salva(fig, FIG, "screws_polar_high")
    pt.plt.close(fig)
    write_rows(GEN / "connection.tex", righe)
    return out


# --------------------------------------------------------------------------- refit of k_b and kappa_r (H5)


def variant(cfg: Path, kappa: float | None = None, kb: float | None = None):
    """Hub configuration with the ring factor and/or k_b overridden: (par, reduced models, FEM with and without
    telescope)."""
    import copy

    d = copy.deepcopy(pr.carica_config(cfg))
    h = d["polar"]["hub"]
    if kappa is not None:
        h["ring"]["stiffness_factor"] = float(kappa)
    if kb is not None:
        h["telescope_attachment"]["k_b_Nm_rad"] = float(kb)
    par = pr.parametri_da_config(cfg, dati=d)
    s_t, s_n = full_fem(par)
    return par, spm.sistemi_per_misure(par), s_t, s_n


def polar_pz(m, full: bool):
    """Poles and zeros of tau_theta -> theta_e in band, without the azimuth pole-zero pairs."""
    ml = m.lineare_completo(0.0, 0.0) if full else lz.linearizza(m.accoppiato, 0.0, 0.0)
    p, z = spm.poli_zeri_modello(ml, "theta_e", 1, 2.0 * F_MAX)
    az = (740.0, 790.0)
    keep = lambda a: np.sort(a[(a > 5.0) & ~((a > az[0]) & (a < az[1]))])  # noqa: E731
    return keep(p), keep(z)


def kappa_sweep() -> dict:
    righe, out = [], {}
    for kappa in (1, 10, 30, 100, 300, 1000, 1e6):
        _, sis, f_t, f_n = variant(CFG_FIT, kappa=kappa)
        ln, fn = polar_pz(sis[False], False)[0][0], polar_pz(f_n, True)[0][0]
        lt, ft = polar_pz(sis[True], False)[0][1], polar_pz(f_t, True)[0][1]
        out[kappa] = (ln, fn, lt, ft)
        lab = r"$10^{6}$ (rigid)" if kappa == 1e6 else f"{kappa:g}"
        righe.append(rf"{lab} & {ln:.0f} & {fn:.0f} & {lt:.0f} & {ft:.0f} \\")
    write_rows(GEN / "kappa_sweep.tex", righe)
    return out


def refit_figure(cfg, cases, fname: str) -> list:
    """Measurement against the 8-DOF model and the FEM with hub of H5, k_b set for each configuration."""
    fig, axs = pt.plt.subplots(len(cases), 2, figsize=(10.0, 2.9 * len(cases)), sharex=True)
    righe = []
    for r, (fname_m, lab, kb, tele) in enumerate(cases):
        from types import SimpleNamespace

        mis = spm.carica_misura(cfg, SimpleNamespace(file=fname_m, asse="polare", theta0=0.0, telescopio=tele,
                                                     etichetta=lab))
        _, sis, f_t, f_n = variant(CFG_FIT, kb=kb)
        lum, fe = (sis[True], f_t) if tele else (sis[False], f_n)
        G_l = spm.frf_velocita(lz.linearizza(lum.accoppiato, 0.0, 0.0), "polare", F_GRID, cfg.Kt)["theta_e"]
        G_f = spm.frf_velocita(fe.lineare_completo(0.0, 0.0), "polare", F_GRID, cfg.Kt)["theta_e"]
        for G, st, nm in ((G_l, pt._stile(0), r"8-DOF, encoder $\dot\theta_e$"),
                          (G_f, pt._stile(1), r"FEM with hub, encoder $\dot\theta_e$")):
            axs[r, 0].plot(F_GRID, pt._db(G), lw=1.3, label=nm, **st)
            axs[r, 1].plot(F_GRID, pt._fase_avvolta(G), lw=1.0, **st)
        axs[r, 0].plot(mis.f, mis.modulo_db, color=pt.INCHIOSTRO, marker="o", ms=1.6, lw=0.8, label="measurement")
        axs[r, 1].plot(mis.f, mis.fase_deg, color=pt.INCHIOSTRO, marker="o", ms=1.6, lw=0.8)
        titolo = lab + (f", $k_b$ = {kb:.3g} N m/rad" if tele else "")
        axs[r, 0].set_title(titolo, fontsize=9, loc="left")
        axs[r, 0].set_ylabel("magnitude [dB (rad/s)/A]")
        axs[r, 1].set_ylabel("phase [deg]")
        axs[r, 1].set_ylim(-180, 180)
        axs[r, 1].set_yticks(np.arange(-180, 181, 90))
        axs[r, 0].legend(loc="lower right", fontsize=6.5)
        for ax in axs[r]:
            ax.set_xscale("log")
            ax.set_xlim(1.0, F_MAX)
        if tele:
            pl_, zl_ = polar_pz(lum, False)
            pf_, zf_ = polar_pz(fe, True)
            righe.append((lab, kb, zl_[0], pl_[0], zf_[0], pf_[0]))
    for ax in axs[-1]:
        ax.set_xlabel("frequency [Hz]")
    pt._titolo(fig, "Refitted hub configuration (H5, kappa_r = 1000) against the measurements", True)
    fig.tight_layout()
    pt.salva(fig, FIG, fname)
    pt.plt.close(fig)
    return righe


# --------------------------------------------------------------------------- stiffness to ground (H6)

GROUND_FITS = (("Azimuth_BodeMeasurement_Po_0.csv", "azimuth, 3.5 A", "azimut", True),
               ("Azimuth_BodeMeasurement_Po_0_dep.csv", "azimuth, 2.36 A", "azimut", True),
               ("VAntiRes/BodeMeasurement_12viti_dueDecimiOmogeneo.csv", "polar, 12 screws, 3 A", "polare", True),
               ("Elevation_BodeMeasurement_bode_New_16_10_26.csv", "polar, 4 screws, 3.5 A", "polare", True),
               ("BodeMeasurement_NoTelescope.csv", "polar without telescope, 2.5 A", "polare", False))


def ground_section(cfg) -> dict:
    """Fit of the stiffness to ground below 5 Hz (hysteretic and viscous loss, model rigid inertia fixed) and
    the comparison of H6 (with the constraint) against the measurements."""
    par6 = pr.parametri_da_config(CFG_GROUND)
    sis6 = spm.sistemi_per_misure(par6)
    f6t, f6n = full_fem(par6)
    Kt = cfg.Kt
    J = {(a, t): vt.inerzia_rigida(lz.linearizza(sis6[t].accoppiato, 0.0, 0.0), a, Kt)
         for a, t in (("azimut", True), ("polare", True), ("polare", False))}
    righe, out = [], {}
    for fname, lab, asse, tele in GROUND_FITS:
        mis = load(cfg, fname, lab)
        h = vt.stima(mis.f, mis.modulo_db, mis.fase_deg, J[(asse, tele)], Kt)
        v = vt.stima(mis.f, mis.modulo_db, mis.fase_deg, J[(asse, tele)], Kt, perdita="viscosa")
        out[lab] = (h, v)
        righe.append(rf"{lab} & {h.J:.3g} & {h.k:.0f} & {h.eta:.2f} & {h.f0:.1f} & {h.rms_db:.1f} / {h.rms_deg:.0f} & "
                     rf"{v.k:.0f} & {v.c:.1f} & {v.rms_db:.1f} / {v.rms_deg:.0f} \\")
    write_rows(GEN / "ground_fit.tex", righe)

    casi = (("Azimuth_BodeMeasurement_Po_0.csv", "azimuth (3.5 A)", "azimut", True),
            ("VAntiRes/BodeMeasurement_12viti_dueDecimiOmogeneo.csv", "polar, 12 screws, 0.2 mm uniform (3 A)",
             "polare", True),
            ("BodeMeasurement_NoTelescope.csv", "polar without telescope (2.5 A)", "polare", False))
    fl = np.logspace(np.log10(0.3), np.log10(F_MAX), 2500)
    fig, axs = pt.plt.subplots(len(casi), 2, figsize=(10.0, 3.0 * len(casi)), sharex=True)
    for r, (fname, lab, asse, tele) in enumerate(casi):
        mis = load(cfg, fname, lab)
        if not tele:
            from types import SimpleNamespace

            mis = spm.carica_misura(cfg, SimpleNamespace(file=fname, asse=asse, theta0=0.0, telescopio=False,
                                                         etichetta=lab))
        ce = spm.canali(asse)[0]
        ml = lz.linearizza(sis6[tele].accoppiato, 0.0, 0.0)
        mf = (f6t if tele else f6n).lineare_completo(0.0, 0.0)
        curve = ((spm.frf_velocita(vt.applica(ml, par6, tele), asse, fl, Kt)[ce], pt._stile(0), 1.4,
                  "8-DOF with ground stiffness"),
                 (spm.frf_velocita(vt.applica(mf, par6, tele), asse, fl, Kt)[ce], pt._stile(1), 1.2,
                  "FEM with hub, with ground stiffness"),
                 (spm.frf_velocita(ml, asse, fl, Kt)[ce], {"color": pt.INCHIOSTRO_SECONDARIO, "linestyle": ":"}, 1.2,
                  "8-DOF without it"))
        for G, st, lw, nm in curve:
            axs[r, 0].plot(fl, pt._db(G), lw=lw, label=nm, **st)
            axs[r, 1].plot(fl, pt._fase_avvolta(G), lw=lw * 0.8, **st)
        axs[r, 0].plot(mis.f, mis.modulo_db, color=pt.INCHIOSTRO, marker="o", ms=1.6, lw=0.8, label="measurement")
        axs[r, 1].plot(mis.f, mis.fase_deg, color=pt.INCHIOSTRO, marker="o", ms=1.6, lw=0.8)
        v = par6.vincolo_terra.per(asse, tele)
        G_pts = spm.frf_velocita(vt.applica(ml, par6, tele), asse, mis.f, Kt)[ce]
        fb = spm.fit_bande(mis.f, mis.G, G_pts)
        axs[r, 0].set_title(f"{lab}: k = {v.k:.4g} N m/rad, eta = {v.eta:g}; 8-DOF fit {fb['1-5 Hz']:.0f} % "
                            f"(1-5 Hz), {fb['5-100 Hz']:.0f} % (5-100 Hz)", fontsize=8, loc="left")
        axs[r, 0].set_ylabel("magnitude [dB (rad/s)/A]")
        axs[r, 1].set_ylabel("phase [deg]")
        axs[r, 1].set_ylim(-180, 180)
        axs[r, 1].set_yticks(np.arange(-180, 181, 90))
        axs[r, 0].legend(loc="lower right", fontsize=6.5)
        for ax in axs[r]:
            ax.set_xscale("log")
            ax.set_xlim(fl[0], fl[-1])
            ax.axvspan(fl[0], mis.f[0], color=pt.GRIGLIA_MINORE, zorder=0)
    for ax in axs[-1]:
        ax.set_xlabel("frequency [Hz]")
    pt._titolo(fig, "H6: 8-DOF model and FEM with hub with the stiffness to ground, against the measurements", True)
    fig.tight_layout()
    pt.salva(fig, FIG, "ground_frf")
    pt.plt.close(fig)
    return out


# --------------------------------------------------------------------------- goodness of fit (NRMSE)

FIT_MEAS = (("Azimuth_BodeMeasurement_Po_0.csv", "azimut", True),
            ("VAntiRes/BodeMeasurement_12viti_dueDecimiOmogeneo.csv", "polare", True),
            ("BodeMeasurement_NoTelescope.csv", "polare", False))
FIT_MODELS = (("7-DOF (baseline)", CFG_7, False, False),
              ("8-DOF, H4$^{a}$", CFG_8, False, False),
              ("8-DOF, H6 without ground", CFG_GROUND, False, False),
              ("8-DOF, H6", CFG_GROUND, True, False),
              ("FEM with hub, H6", CFG_GROUND, True, True))


def fit_section(cfg) -> dict:
    """Fit NRMSE (System Identification Toolbox definition) of the model versions on the reference measurements,
    on the complex response and on the magnitude in dB, on all points and by band."""
    from types import SimpleNamespace

    misure = [spm.carica_misura(cfg, SimpleNamespace(file=f, asse=a, theta0=0.0, telescopio=t, etichetta=f))
              for f, a, t in FIT_MEAS]
    bande = [b for b, _, _ in spm.BANDE_FIT]
    righe_c, righe_d, out = [], [], {}
    cache: dict = {}
    for nome, cfgp, terra, fem_ in FIT_MODELS:
        if cfgp not in cache:
            par = pr.parametri_da_config(cfgp)
            cache[cfgp] = (par, spm.sistemi_per_misure(par), full_fem(par) if cfgp == CFG_GROUND else None)
        par, sis, fe = cache[cfgp]
        vc, vd = [], []
        for m, (_, asse, tele) in zip(misure, FIT_MEAS):
            ml = (fe[0] if tele else fe[1]).lineare_completo(0.0, 0.0) if fem_ else lz.linearizza(sis[tele].accoppiato, 0.0, 0.0)
            if terra:
                ml = vt.applica(ml, par, tele)
            G = spm.frf_velocita(ml, asse, m.f, cfg.Kt)[spm.canali(asse)[0]]
            fc = spm.fit_bande(m.f, m.G, G)
            fd = spm.fit_bande(m.f, m.modulo_db, 20 * np.log10(np.abs(G)))
            out[(nome, asse, tele)] = (fc, fd)
            fmt = lambda v: "--" if not np.isfinite(v) else (r"$<-999$" if v < -999 else f"${v:.0f}$")  # noqa: E731
            vc += [fmt(fc[b]) for b in bande]
            vd += [fmt(fd[b]) for b in bande]
        righe_c.append(f"{nome} & " + " & ".join(vc) + r" \\")
        righe_d.append(f"{nome} & " + " & ".join(vd) + r" \\")
    write_rows(GEN / "fit_complex.tex", righe_c)
    write_rows(GEN / "fit_db.tex", righe_d)
    return out


# --------------------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(exist_ok=True)
    GEN.mkdir(exist_ok=True)
    pt.applica_stile()

    par7, sis7 = build_models(CFG_7)
    par8, sis8 = build_models(CFG_8)
    sis8_r1_par, sis8_r1 = build_models(CFG_8_RING1)
    assert cc.variante_di(sis7[True]) == "catena" and cc.variante_di(sis8[True]) == "mozzo"

    # ----- measurements and the two comparisons -------------------------------------------------
    cfg = par7.sperimentale
    misure = {Path(v.file).stem: spm.carica_misura(cfg, v) for v in cfg.misure}
    sel = {"az": "Azimuth_BodeMeasurement_Po_0", "pol": "Elevation_BodeMeasurement_bode_New_16_10_26",
           "nt": "BodeMeasurement_NoTelescope"}
    titles = {"az": "Azimuth axis (polar at 0 deg)", "pol": "Polar axis, with telescope",
              "nt": "Polar axis, without telescope"}
    for tag, par, sis, model in (("7", par7, sis7, "7-DOF"), ("8", par8, sis8, "8-DOF")):
        for k, stem in sel.items():
            conf = spm.confronta(par, misure[stem], sis, F_GRID)
            bode_comparison(conf, model, f"cmp{tag}_{k}", f"{titles[k]}: {model} model against measurement")

    # ----- FEM (all 18 elements): mode shapes and frequency responses in the units of the measurement --------------
    f7t, f7n = full_fem(par7)
    f8t, f8n = full_fem(par8)
    pz_f7 = fem_frf_figure(f7t, f7n, misure, sel, cfg.Kt, "fem_frf_7",
                           "18-element FEM (baseline configuration) against measurement: encoder and motor velocity per ampere")
    pz_f8 = fem_frf_figure(f8t, f8n, misure, sel, cfg.Kt, "fem_frf_8",
                           "FEM with hub (18 elements, telescope on k_b) against measurement: encoder and motor velocity per ampere")
    mode_shapes_figure(f7t, (("azimut", "azimuth"), ("polare", "polar")), 2, "fem_modes_7",
                       "Flexible modes of the 18-element FEM, baseline configuration (telescope frozen at theta = 0)")
    mode_shapes_figure(f8t, (("polare", "polar"),), 3, "fem_modes_8",
                       "Flexible modes of the FEM with hub, polar axis (theta = 0)")

    # ----- reduced models against the FEM they come from (isolated axes, telescope frozen, theta = 0) ------
    def iso_lumped(s, asse):
        f = fem.frequenze_proprie(*s.isolato_ridotto(asse, 0.0)[:2])
        return f[f > 1.0]

    def iso_fem(s, asse):
        f = s.analisi(asse, 0.0).f
        return f[f > 1.0]

    def zero(m, ch, j, full):
        ml = m.lineare_completo(0.0, 0.0) if full else lz.linearizza(m.accoppiato, 0.0, 0.0)
        z = spm.poli_zeri_modello(ml, ch, j, 3.0 * F_MAX)[1]
        return np.sort(z[z > 5.0])

    for tag, sl, sf, sln, sfn in (("7", sis7[True], f7t, sis7[False], f7n), ("8", sis8[True], f8t, sis8[False], f8n)):
        a_l, a_f = iso_lumped(sl, "azimut"), iso_fem(sf, "azimut")
        p_l, p_f = iso_lumped(sl, "polare"), iso_fem(sf, "polare")
        n_l, n_f = iso_lumped(sln, "polare"), iso_fem(sfn, "polare")
        righe = [("azimuth, mode 1", a_f[0], a_l[0]), ("azimuth, mode 2", a_f[1], a_l[1])]
        righe += [(f"polar, mode {k + 1}", p_f[k], p_l[k]) for k in range(3 if tag == "8" else 2)]
        righe += [("polar without telescope, mode 1", n_f[0], n_l[0])]
        za_l, za_f = zero(sl, "psi_e", 0, False), zero(sf, "psi_e", 0, True)
        ze_l, ze_f = zero(sl, "theta_e", 1, False), zero(sf, "theta_e", 1, True)
        # zeros of the polar channel at the azimuth modes are pole-zero pairs of the other axis: not listed
        ze_l, ze_f = ze_l[ze_l < 700.0], ze_f[ze_f < 700.0]
        zm_l, zm_f = zero(sl, "theta_m", 1, False), zero(sf, "theta_m", 1, True)
        righe += [(r"zero $\tau_\psi\to\psi_e$ (load)", za_f[0], za_l[0])]
        righe += [(rf"zero $\tau_\theta\to\theta_e$, {k + 1}", ze_f[k], ze_l[k]) for k in range(min(2, ze_f.size))]
        righe += [(r"zero $\tau_\theta\to\theta_m$, 1", zm_f[0], zm_l[0])]
        write_rows(GEN / f"lumped_vs_fem_{tag}.tex",
                   [rf"{lab} & {num(ff, 4)} & {num(fl_, 4)} & ${pct(fl_, ff)}$ \\" for lab, ff, fl_ in righe])

    # ----- measured features -------------------------------------------------------------------------
    el, nt, az = misure[sel["pol"]], misure[sel["nt"]], misure[sel["az"]]
    f_zero_meas = float(el.f[(el.f > 12) & (el.f < 25)][np.argmin(el.modulo_db[(el.f > 12) & (el.f < 25)])])
    f48, _ = local_peak(el.f, el.modulo_db, 35, 70)
    f90, _ = local_peak(el.f, el.modulo_db, 75, 110)
    f_z152 = float(el.f[(el.f > 120) & (el.f < 175)][np.argmin(el.modulo_db[(el.f > 120) & (el.f < 175)])])
    f193, _ = local_peak(el.f, el.modulo_db, 170, 230)
    f_az_zero = float(az.f[(az.f > 8) & (az.f < 16)][np.argmin(az.modulo_db[(az.f > 8) & (az.f < 16)])])
    f_az_pole, _ = local_peak(az.f, az.modulo_db, 300, 700)
    f_az_lf, _ = local_peak(az.f, az.modulo_db, 1.5, 8)
    nt_zero1 = float(nt.f[(nt.f > 350) & (nt.f < 420)][np.argmin(nt.modulo_db[(nt.f > 350) & (nt.f < 420)])])
    nt_pole, _ = local_peak(nt.f, nt.modulo_db, 415, 470)
    nt_zero2 = float(nt.f[(nt.f > 470) & (nt.f < 560)][np.argmin(nt.modulo_db[(nt.f > 470) & (nt.f < 560)])])

    # ----- poles and zeros of the two models ------------------------------------------------------------
    def first_above(a, x):
        a = np.asarray(a)
        a = a[a > x]
        return float(a.min()) if a.size else float("nan")

    p7, z7e = pole_zero(sis7[True], "theta_e", 1)
    _, z7m = pole_zero(sis7[True], "theta_m", 1)
    p8, z8e = pole_zero(sis8[True], "theta_e", 1)
    _, z8m = pole_zero(sis8[True], "theta_m", 1)
    pa7, za7 = pole_zero(sis7[True], "psi_e", 0)
    pa8, za8 = pole_zero(sis8[True], "psi_e", 0)
    pn7, _ = pole_zero(sis7[False], "theta_e", 1)
    pn8, _ = pole_zero(sis8[False], "theta_e", 1)
    pn8r1, _ = pole_zero(sis8_r1[False], "theta_e", 1)

    s8 = sis8[True]
    Jc = float(s8.payload.I_yy)
    kg = float(s8.payload.k_g(0.0))
    kb = 9.8e3
    kb_fit = (2 * np.pi * f_zero_meas) ** 2 * Jc - kg
    # two-mass estimate from the measured zero and pole (rigid-body line fixes J1 + Jc)
    Jtot = float(read_csv(RES_7 / "tables" / "modal_modes_polare.csv")[0]["massa_modale"])
    ratio = (f48 / f_zero_meas) ** 2 - 1.0
    J1_need = Jtot / (1 + ratio)

    num_ = {
        "fzMeas": f"{f_zero_meas:.2f}", "fpMeas": f"{f48:.0f}", "fNineMeas": f"{f90:.0f}",
        "fzHundredFiftyMeas": f"{f_z152:.0f}", "fTwoHundredMeas": f"{f193:.0f}",
        "fzSevenMotor": f"{first_above(z7m, 5):.2f}", "fzSevenEnc": f"{first_above(z7e, 5):.0f}",
        "fpSevenOne": f"{first_above(p7, 5):.1f}",
        "fzEightEnc": f"{first_above(z8e, 5):.2f}", "fzEightMotorA": f"{first_above(z8m, 5):.2f}",
        "fzEightMotorB": f"{sorted(z8m[z8m > 5])[1]:.1f}",
        "fpEightOne": f"{first_above(p8, 5):.1f}", "fpEightTwo": f"{sorted(p8[p8 > 5])[1]:.0f}",
        "fAzZeroMeas": f"{f_az_zero:.2f}", "fAzZeroSeven": f"{first_above(za7, 5):.2f}",
        "fAzZeroEight": f"{first_above(za8, 5):.2f}", "fAzPoleMeas": f"{f_az_pole:.0f}",
        "fAzPoleSeven": f"{max(pa7[(pa7 > 700) & (pa7 < 800)]):.0f}",
        "fAzPoleEight": f"{max(pa8[(pa8 > 700) & (pa8 < 800)]):.0f}",
        "fAzLFMeas": f"{f_az_lf:.1f}",
        "fNtZeroOne": f"{nt_zero1:.0f}", "fNtPole": f"{nt_pole:.0f}", "fNtZeroTwo": f"{nt_zero2:.0f}",
        "fNtSeven": f"{first_above(pn7, 5):.0f}", "fNtEight": f"{first_above(pn8, 5):.0f}",
        "fNtEightRingOne": f"{first_above(pn8r1, 5):.0f}",
        "Jtot": f"{Jtot:.3f}", "JcEight": f"{Jc:.3f}", "kgEight": f"{kg:.2f}",
        "kgSeven": f"{float(sis7[True].payload.k_g(0.0)):.2f}",
        "kbCur": sci(kb, 3), "kbFit": sci(kb_fit, 3),
        "kbTwoMass": sci((2 * np.pi * f_zero_meas) ** 2 * (Jtot - J1_need), 3),
        "JoneNeed": f"{J1_need:.3f}", "JcNeed": f"{Jtot - J1_need:.3f}", "ratioNeed": f"{ratio:.1f}",
        "fzScaleHalf": f"{f_zero_meas * np.sqrt(0.5):.1f}", "fzScaleDouble": f"{f_zero_meas * np.sqrt(2.0):.1f}",
    }
    def fl(a, lo=5.0, hi=np.inf, k=0):
        a = np.sort(np.asarray(a)[(np.asarray(a) > lo) & (np.asarray(a) < hi)])
        return float(a[k]) if a.size > k else float("nan")

    num_.update({
        "fAzZeroFemSeven": f"{fl(pz_f7['az'][1]):.2f}", "fAzPoleFemSeven": f"{fl(pz_f7['az'][0], 700, 800):.0f}",
        "fzFemSevenEnc": f"{fl(pz_f7['pol'][1]):.0f}", "fpFemSevenOne": f"{fl(pz_f7['pol'][0]):.1f}",
        "fpFemSevenTwo": f"{fl(pz_f7['pol'][0], k=2):.0f}", "fNtFemSeven": f"{fl(pz_f7['nt'][0]):.0f}",
        "fAzZeroFemEight": f"{fl(pz_f8['az'][1]):.2f}", "fAzPoleFemEight": f"{fl(pz_f8['az'][0], 700, 800):.0f}",
        "fzFemEightEnc": f"{fl(pz_f8['pol'][1]):.2f}", "fpFemEightOne": f"{fl(pz_f8['pol'][0]):.1f}",
        "fpFemEightTwo": f"{fl(pz_f8['pol'][0], k=1):.1f}", "fNtFemEight": f"{fl(pz_f8['nt'][0]):.0f}",
        "fzFemSevenRear": f"{fl(pz_f7['pol'][1]):.1f}",
    })
    sc = screw_test(cfg, Jtot)
    s4, s2, s0, s4old = sc[44], sc[2], sc[0], sc[4]
    num_.update({
        "scZfourNew": f"{s4['fz']:.2f}", "scPfourNew": f"{s4['fp']:.1f}", "scZtwo": f"{s2['fz']:.2f}",
        "scPtwo": f"{s2['fp']:.1f}", "scZfourOld": f"{s4old['fz']:.2f}", "scPfourOld": f"{s4old['fp']:.1f}",
        "scPhFour": f"{s4['ph0']:.1f}", "scPhTwo": f"{s2['ph0']:.1f}",
        "scRepZ": f"{100 * (s4['fz'] / s4old['fz'] - 1):+.1f}", "scRepP": f"{100 * (s4['fp'] / s4old['fp'] - 1):+.1f}",
        "scShiftZ": f"{100 * (s2['fz'] / s4['fz'] - 1):+.1f}", "scShiftP": f"{100 * (s2['fp'] / s4['fp'] - 1):+.1f}",
        "scKbRatio": f"{s2['kb'] / s4['kb']:.2f}", "scKbFour": sci(s4['kb'], 3), "scKbTwo": sci(s2['kb'], 3),
        "scJshTwo": f"{s2['Jsh']:.3f}", "scJshFour": f"{s4['Jsh']:.3f}", "scZeroNotchA": f"{s0['n1']:.1f}",
        "scZeroNotchB": f"{s0['n2']:.1f}", "scZeroPk": f"{s0['p75']:.0f}", "scZeroZ": f"{s0['z170']:.0f}",
        "scZeroP": f"{s0['p225']:.0f}",
    })
    # ----- 8 and 12 screws, and refit of the free parameters (H5) ------------------------------------------
    par_fit = pr.parametri_da_config(CFG_FIT)
    s_fit = spm.sistemi_per_misure(par_fit)
    IJJ_fit, kg_fit = float(s_fit[True].payload.I_yy), float(s_fit[True].payload.k_g(0.0))
    hi = high_screws(cfg, IJJ_fit, kg_fit, sc)
    sw = kappa_sweep()
    kb_of = lambda fa: (2 * np.pi * fa) ** 2 * IJJ_fit - kg_fit  # noqa: E731
    kb2, kb12 = round(kb_of(sc[2]["fz"]), -2), round(hi["12 screws, 0.2 uniform"]["fa"] and kb_of(hi["12 screws, 0.2 uniform"]["fa"]), -3)
    casi = (("Elevation_BodeMeasurement_bode_New_16_10_26.csv", "4 screws (16 Oct)", par_fit.mozzo.k_b, True),
            ("Elevation_BodeMeasurement_2viti.csv", "2 screws (19 Oct)", kb2, True),
            ("VAntiRes/BodeMeasurement_12viti_dueDecimiOmogeneo.csv", "12 screws, 0.2 mm uniform (19 Oct)", kb12, True),
            ("BodeMeasurement_NoTelescope.csv", "without telescope", par_fit.mozzo.k_b, False))
    rf = refit_figure(cfg, casi, "refit_frf")
    righe_rf = []
    for lab, kb, zl, pl_, zf, pf in rf:
        righe_rf.append(rf"{lab} & ${sci(kb, 3)}$ & {zl:.1f} & {pl_:.0f} & {zf:.1f} & {pf:.0f} \\")
    write_rows(GEN / "refit.tex", righe_rf)
    e8 = [v["kb"] for k, v in hi.items() if k.startswith("8")]
    e12 = [v["kb"] for k, v in hi.items() if k.startswith("12") and "0.2/0.2" not in k]
    num_.update({
        "kappaNew": "1000", "kbNew": sci(par_fit.mozzo.k_b, 3), "kbTwoFit": sci(kb2, 2), "kbTwelveFit": sci(kb12, 2),
        "kbEightRange": rf"${sci(min(e8), 2)}$--${sci(max(e8), 2)}$",
        "kbTwelveRange": rf"${sci(min(e12), 2)}$--${sci(max(e12), 2)}$",
        "kbFourFitNew": sci(kb_of(sc[44]["fz"]), 3), "kbFourFitOld": sci(kb_of(sc[4]["fz"]), 3),
        "IJJfit": f"{IJJ_fit:.3f}",
        "swLumpHundred": f"{sw[100][0]:.0f}", "swFemHundred": f"{sw[100][1]:.0f}",
        "swLumpThousand": f"{sw[1000][0]:.0f}", "swFemThousand": f"{sw[1000][1]:.0f}",
        "swLumpRigid": f"{sw[1e6][0]:.0f}", "swFemRigid": f"{sw[1e6][1]:.0f}",
        "swLumpOne": f"{sw[1][0]:.0f}", "swLumpThree": f"{sw[300][0]:.0f}",
        "rfZfour": f"{rf[0][2]:.1f}", "rfPfour": f"{rf[0][3]:.0f}", "rfZtwelve": f"{rf[2][2]:.1f}",
        "rfPtwelve": f"{rf[2][3]:.0f}", "rfPtwelveFem": f"{rf[2][5]:.0f}",
        "hiNotchMin": f"{min(v['fa'] for k, v in hi.items() if '0.2/0.2' not in k):.0f}",
        "hiNotchMax": f"{max(v['fa'] for k, v in hi.items() if '0.2/0.2' not in k):.0f}",
        "hiOddNotch": f"{hi['12 screws, 0.2/0.2']['fa']:.0f}",
    })
    fs = fit_section(cfg)
    h6 = {k[1:]: v for k, v in fs.items() if k[0] == "8-DOF, H6"}
    h7 = {k[1:]: v for k, v in fs.items() if k[0] == "7-DOF (baseline)"}
    num_.update({
        "fitAzLow": f"{h6[('azimut', True)][0]['1-5 Hz']:.0f}", "fitPolLow": f"{h6[('polare', True)][0]['1-5 Hz']:.0f}",
        "fitNtLow": f"{h6[('polare', False)][0]['1-5 Hz']:.0f}",
        "fitAzLowSeven": f"{h7[('azimut', True)][0]['1-5 Hz']:.0f}", "fitPolLowSeven": f"{h7[('polare', True)][0]['1-5 Hz']:.0f}",
        "fitPolMid": f"{h6[('polare', True)][0]['5-100 Hz']:.0f}", "fitPolMidDb": f"{h6[('polare', True)][1]['5-100 Hz']:.0f}",
        "fitAzMidDb": f"{h6[('azimut', True)][1]['5-100 Hz']:.0f}", "fitNtMidDb": f"{h6[('polare', False)][1]['5-100 Hz']:.0f}",
    })
    gr = ground_section(cfg)
    ha, va_ = gr["azimuth, 3.5 A"]
    hal, _ = gr["azimuth, 2.36 A"]
    hp, vp = gr["polar, 12 screws, 3 A"]
    hp4, _ = gr["polar, 4 screws, 3.5 A"]
    hn, _ = gr["polar without telescope, 2.5 A"]
    num_.update({
        "gAzK": f"{ha.k:.0f}", "gAzEta": f"{ha.eta:.2f}", "gAzFz": f"{ha.f0:.1f}", "gAzJ": f"{ha.J:.2f}",
        "gAzKlow": f"{hal.k:.0f}", "gAzEtaLow": f"{hal.eta:.2f}",
        "gPolK": f"{hp.k:.0f}", "gPolEta": f"{hp.eta:.2f}", "gPolFz": f"{hp.f0:.1f}", "gPolJ": f"{hp.J:.3f}",
        "gPolKvisc": f"{vp.k:.0f}", "gPolC": f"{vp.c:.1f}", "gPolRmsH": f"{hp.rms_db:.1f}",
        "gPolRmsV": f"{vp.rms_db:.1f}", "gPolDegH": f"{hp.rms_deg:.0f}", "gPolDegV": f"{vp.rms_deg:.0f}",
        "gPolKfour": f"{hp4.k:.0f}", "gNtK": f"{hn.k:.0f}", "gNtEta": f"{hn.eta:.2f}", "gNtFz": f"{hn.f0:.0f}",
        "gNtJ": f"{hn.J:.3f}",
    })
    lump8 = read_csv(RES_8 / "tables" / "lumped_parameters.csv")
    node8 = {(r["asse_o_nodo"], r["elemento_o_ruolo"], r["nodo_a"]): r for r in lump8}
    # hub node H (8-DOF): inertia of the hub node
    num_["JhEight"] = f"{float([r for r in lump8 if r['asse_o_nodo'] == 'nodo' and r['nodo_a'] == 'H'][0]['J_kgm2']):.4f}"
    num_["JpSevenNode"] = f"{float([r for r in read_csv(RES_7 / 'tables' / 'lumped_parameters.csv') if r['asse_o_nodo'] == 'nodo' and r['nodo_a'] == 'P'][0]['J_kgm2']):.2e}"
    num_["JratioEight"] = f"{Jc / float(num_['JhEight']):.1f}"
    from gimbal_mdof import mozzo as mz  # noqa: E402

    rg = mz.rigidezze(par8.mozzo)
    num_.update({"kTube": sci(rg.k_tubo, 3), "kFlange": sci(rg.k_flangia, 3), "kRingZero": f"{rg.k_semianelli:.0f}",
                 "kRingHalf": sci(2 * rg.k_anello, 3), "kHalfHub": sci(rg.k_lato, 3), "kHub": sci(rg.k_totale, 3),
                 "kHubRingAlone": sci(mz.rigidezze(sis8_r1_par.mozzo).k_totale, 3), "ringFactor": f"{rg.fattore:g}"})
    macros = "% generated by build_development_report.py: do not edit\n" + "".join(
        rf"\newcommand{{\{k}}}{{{v}}}" + "\n" for k, v in num_.items())
    (GEN / "numbers.tex").write_text(macros)
    for k, v in num_.items():
        print(f"{k:22s} {v}")

    # ----- tables --------------------------------------------------------------------------------------------
    write_tables(par7, par8)

    # ----- figures that come from the model runs (copied: the report is self-contained) ------------------------
    for src, dst in ((RES_7 / "figures/geometria_azimut.pdf", "geometria_azimut.pdf"),
                     (RES_7 / "figures/geometria_polare.pdf", "geometria_polare.pdf"),
                     (RES_8 / "figures/geometria_mozzo.pdf", "geometria_mozzo.pdf"),
                     (RES_7 / "figures/exp_azimut_sovrapposte.pdf", "exp_azimut_sovrapposte.pdf")):
        shutil.copy(src, FIG / dst)

    # ----- spring--inertia chain of the 8-DOF polar axis, with the symbols of the report -------------------
    fig = pt.catena_ramo_poli_zeri(ca.catena_ramo_modello(sis8[True], 0.0), segnaposto=True, stimata=False)
    simboli = {r"$\theta_P$": r"$\theta_h$", r"$\theta_c$": r"$\theta_p$"}
    for t in fig.findobj(pt.plt.Text):
        if t.get_text() in simboli:
            t.set_text(simboli[t.get_text()])
    # plotting.py marks the node held still to compute each zero with a hatch, which reads as a ground: the
    # hatches are removed (the zeros are explained in the text) together with their note
    from matplotlib.patches import Rectangle
    for ax in fig.axes:
        for r in list(ax.patches):
            if isinstance(r, Rectangle) and r.get_hatch():
                r.remove()
    for t in fig.findobj(pt.plt.Text):
        t.set_text(t.get_text().replace("   (hatched: node clamped to compute the zeros)", ""))
    fig.suptitle("Polar chain of the 8-DOF model, theta0 = 0 deg (J in kg m^2, k in N m/rad)\n"
                 "[RESULT FROM PLACEHOLDER PARAMETERS]", fontsize=10, color=pt.INCHIOSTRO)
    pt.salva(fig, FIG, "chain8_polar")
    pt.plt.close(fig)

    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "development_report.tex"],
                       cwd=AQUI, check=True)


def write_tables(par7, par8) -> None:
    """Element totals, lumped parameters and modal data as LaTeX fragments."""
    # element totals (7-DOF geometry): one row per element, both axes
    for ax, nm in (("azimut", "azimuth"), ("polare", "polar")):
        rows = [r for r in read_csv(RES_7 / "tables" / f"element_breakdown_{ax}.csv") if r["segmento"] == "TOTALE"]
        out = []
        for r in rows:
            k = r["k_Nm_rad"]
            out.append(rf"{r['elemento']} & {float(r['L_mm']):.4g} & {1e3 * float(r['m_kg']):.4g} & "
                       rf"${sci(float(r['J_kgm2']), 3)}$ & " + (rf"${sci(float(k), 3)}$" if k not in ("", "nan", "inf") else "--") + r" \\")
        write_rows(GEN / f"elements_{nm}.tex", out)

    # lumped parameters, with the physical component each element / node stands for
    desc_el = {
        ("7", "azimut", "1"): "shaft with flange, encoder to motor", ("7", "azimut", "2"): "shaft, motor to fork",
        ("7", "polare", "1"): "encoder-side shaft and first bearing seat",
        ("7", "polare", "2"): "second bearing seat and shaft up to the motor",
        ("7", "polare", "3"): "shaft, motor to rear end",
        ("8", "azimut", "1"): "shaft with flange, encoder to motor", ("8", "azimut", "2"): "shaft, motor to fork",
        ("8", "polare", "1"): "encoder-side shaft and first half of the hub",
        ("8", "polare", "2"): "second half of the hub and shaft up to the motor",
        ("8", "polare", "3"): "shaft, motor to rear end",
        ("8", "polare", "4"): "telescope interface (ring--telescope connection)",
    }
    desc_nd = {
        ("7", "azimut", "e"): "encoder end of the azimuth shaft", ("7", "azimut", "m"): "motor node, torque $\\tau_\\psi$",
        ("7", "azimut", "L"): "load: fork and polar shafts (CAD); telescope through its tensor",
        ("7", "polare", "e"): "encoder end of the polar shaft", ("7", "polare", "p"): "telescope seat P, rigid telescope",
        ("7", "polare", "m"): "motor node, torque $\\tau_\\theta$", ("7", "polare", "r"): "rear end of the polar shaft",
        ("8", "azimut", "e"): "encoder end of the azimuth shaft", ("8", "azimut", "m"): "motor node, torque $\\tau_\\psi$",
        ("8", "azimut", "L"): "load: fork and polar shafts (CAD); telescope through its tensor",
        ("8", "polare", "e"): "encoder end of the polar shaft",
        ("8", "polare", "p"): "hub centre H: side hubs, ring, tube section in the ring",
        ("8", "polare", "m"): "motor node, torque $\\tau_\\theta$", ("8", "polare", "r"): "rear end of the polar shaft",
    }
    for tag, res in (("7", RES_7), ("8", RES_8)):
        rows = read_csv(res / "tables" / "lumped_parameters.csv")
        el = [r for r in rows if r["asse_o_nodo"] in ("azimut", "polare")]
        nd = [r for r in rows if r["asse_o_nodo"] == "nodo"]
        out = []
        for r in el:
            ax = r["asse_o_nodo"]
            sym = r"\psi" if ax == "azimut" else r"\theta"
            name = r["elemento_o_ruolo"]
            lab = r"$k_b$" if (tag == "8" and ax == "polare" and name == "4") else rf"$k_{{{sym},{name}}}$"
            sup = rf"{r['nodo_a']}--{r['nodo_b']}"
            out.append(rf"{lab} & {desc_el[(tag, ax, name)]} & {sup} & {r['elementi_FEM'].replace('b', '--')} & "
                       rf"${sci(float(r['J_kgm2']), 3)}$ & ${sci(float(r['k_Nm_rad']), 3)}$ \\")
        write_rows(GEN / f"lumped_elements_{tag}.tex", out)
        out = []
        for i, r in enumerate(nd):
            ax = "azimut" if i < 3 else "polare"
            role = r["elemento_o_ruolo"]
            sym = (r"\psi" if ax == "azimut" else r"\theta")
            if tag == "8" and ax == "polare" and role == "p":
                lab = r"$J_{\theta,h}$"
            else:
                lab = rf"$J_{{{sym},{role}}}$"
            out.append(rf"{lab} & {desc_nd[(tag, ax, role)]} & {r['nodo_a']} & ${sci(float(r['J_kgm2']), 3)}$ \\")
        write_rows(GEN / f"lumped_nodes_{tag}.tex", out)

    # modal table of the 18-element FEM, baseline configuration (kept modes)
    out = []
    for ax, nm in (("azimut", "azimuth"), ("polare", "polar")):
        for r in read_csv(RES_7 / "tables" / f"modal_modes_{ax}.csv")[:5]:
            kept = {"si": "kept", "no": "no"}.get(r["tenuto"], r["tenuto"])
            f = float(r["f_Hz"])
            out.append(rf"{nm} & {r['modo']} & ${num(f, 4)}$ & ${num(float(r['zeta']), 2)}$ & ${num(float(r['hankel_rad_Nm']), 2) if r['hankel_rad_Nm'] != 'inf' else r'\infty'}$ & "
                       rf"{kept} \\")
    write_rows(GEN / "modal_modes.tex", out)


if __name__ == "__main__":
    main()
