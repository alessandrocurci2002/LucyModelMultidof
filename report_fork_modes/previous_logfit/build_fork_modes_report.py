"""Builds the report on the stator modal masses of the elevation axis (two masses) and on the flexible azimuth load.

    ../.venv/bin/python -I build_fork_modes_report.py [--no-latex]

Inputs: the summary.json and tables of analisi_modi_forcella.py in results_fork_modes_2masses (the reported model)
and in results_fork_modes (the four-mass variant, for the comparison). No fit is run here: the models are rebuilt
from the fitted parameters to draw the figures in the fit band (1-150 Hz). The text is in fork_modes_report.tex; this
script writes the figures (figures/), the tables and numbers (generated/) and compiles with latexmk.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt  # noqa: E402

import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
RES2, RES4 = ROOT / "results_fork_modes_2masses", ROOT / "results_fork_modes"
BAND = (1.0, 150.0)
ZOOM = (40.0, 150.0)
TEX_NAME = {"k_b": r"$k_b$", "c_b": r"$c_b$", "k_g": r"$k_g$", "eta_g": r"$\eta_g$", "tau_d": r"$\tau_d$",
            "alpha": r"$\alpha$", "k_LA": r"$k_{LA}$", "c_LA": r"$c_{LA}$", "k_3": r"$k_3$", "c_3": r"$c_3$",
            "k_2": r"$k_2$", "c_2": r"$c_2$"}
UNIT = {"N m/rad": r"N\,m/rad", "N m s/rad": r"N\,m\,s/rad", "s": "s", "Hz": "Hz", "1/(kg m^2)": r"1/(kg\,m$^2$)",
        "": "--"}
KIND = {"antirisonanza": "antiresonance", "risonanza": "resonance"}


def tname(n: str) -> str:
    if n in TEX_NAME:
        return TEX_NAME[n]
    base, _, r = n.rpartition("_")
    return {"f": f"$f_{r}$", "zeta": rf"$\zeta_{r}$", "a": f"$a_{r}$"}[base]


def num(x: float, d: int = 4) -> str:
    """Number in LaTeX, powers of ten as \\times 10^n."""
    if not np.isfinite(x):
        return "--"
    s = f"{x:.{d}g}"
    if "e" in s:
        m, e = s.split("e")
        return rf"\ensuremath{{{m}\times10^{{{int(e)}}}}}"
    return s


def std_txt(v: float) -> str:
    p = 100 * v
    return r"$>10^3$\,\%" if p > 1e3 else f"{p:.1f}\\,\\%"


def neg(t: str) -> str:
    """Minus sign of a number in text mode."""
    return t.replace("-", "$-$", 1) if t.startswith("-") else t


def esc(t: str) -> str:
    return t.replace("%", r"\%")


def rows_tex(rows: list[list[str]]) -> str:
    return " \\\\\n".join(" & ".join(r) for r in rows)


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def starts(path: Path, axis: str) -> dict[str, list[float]]:
    """Distinct start values of each parameter over the starts of ``axis``."""
    out: dict[str, set[float]] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["axis"] == axis:
                for k, v in json.loads(r["start"]).items():
                    out.setdefault(k, set()).add(float(f"{v:.4g}"))
    return {k: sorted(v) for k, v in out.items()}


def costs(path: Path, axis: str) -> list[float]:
    with open(path, newline="") as fh:
        return [float(r["cost"]) for r in csv.DictReader(fh) if r["axis"] == axis]


def param_rows(par: list[mf.Parametro], est: dict, std: dict, st: dict[str, list[float]], h7: dict[str, float],
               at_bound: list[str]) -> list[list[str]]:
    rows = []
    for p in par:
        s0 = st.get(p.nome, [])
        if p.nome.startswith("a_"):
            s0 = sorted({abs(x) for x in s0})
            start = ", ".join(rf"$\pm{x:g}$" for x in s0)
        else:
            start = ", ".join(num(x, 3) for x in s0) or "--"
            if len(s0) > 3:
                start = f"{len(s0)} values, {num(s0[0], 2)}--{num(s0[-1], 2)}"
        rows.append([tname(p.nome), UNIT[p.unita], num(est[p.nome]), num(h7[p.nome]) if p.nome in h7 else "--", start,
                     f"[{num(p.lo, 3)}, {num(p.hi, 3)}]", std_txt(std[p.nome]),
                     "yes" if p.nome in at_bound else ""])
    return rows


# ---------------------------------------------------------------------------
# Figures (fit band only)
# ---------------------------------------------------------------------------


def in_band(m) -> np.ndarray:
    return (m.f >= BAND[0]) & (m.f <= BAND[1])


def fig_bode(f: np.ndarray, misure: list, G: np.ndarray, label: str, ylim: tuple, ylim_zoom: tuple, title: str,
             seg: bool):
    fig, assi = plt.subplots(2, 2, figsize=(12.5, 6.8), sharex="col", gridspec_kw={"width_ratios": [1.6, 1.0]})
    for col, lim in ((0, BAND), (1, ZOOM)):
        ax_m, ax_f = assi[0, col], assi[1, col]
        for k, m in enumerate(misure):
            b = in_band(m)
            st = {"color": afm.COLORI_MISURE[k % len(afm.COLORI_MISURE)], "marker": "o", "ms": 1.6, "lw": 0.7,
                  "alpha": 0.85}
            ax_m.plot(m.f[b], m.modulo_db[b], label=f"measured: {m.etichetta}", **st)
            ax_f.plot(m.f[b], m.fase_deg[b], **st)
        st = pt._stile(1)
        ax_m.plot(f, pt._db(G), lw=1.8, label=label, **st)
        ax_f.plot(f, pt._fase_avvolta(G), lw=1.4, **st)
        for ax in (ax_m, ax_f):
            ax.set_xscale("log")
            ax.set_xlim(*lim)
        ax_m.set_ylim(*(ylim if col == 0 else ylim_zoom))
        ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
        ax_f.set_ylabel("phase [deg]")
        ax_f.set_ylim(-180, 180)
        ax_f.set_yticks(np.arange(-180, 181, 90))
        ax_f.set_xlabel("frequency [Hz]")
    assi[0, 0].legend(loc="lower left", fontsize=7)
    assi[0, 1].set_title(f"zoom {ZOOM[0]:g}-{ZOOM[1]:g} Hz", fontsize=9)
    pt._titolo(fig, title, seg)
    return fig


def fig_contributions(f: np.ndarray, m, mod: mf.ModelloPolare, seg: bool):
    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(10.5, 6.8), sharex=True)
    b = in_band(m)
    ax_m.plot(m.f[b], m.modulo_db[b], color=pt.INCHIOSTRO_TENUE, marker="o", ms=1.4, lw=0.6,
              label=f"measured: {m.etichetta}")
    ax_f.plot(m.f[b], m.fase_deg[b], color=pt.INCHIOSTRO_TENUE, marker="o", ms=1.4, lw=0.6)
    G = mod.frf(f)
    ax_m.plot(f, pt._db(G), color=pt.INCHIOSTRO, lw=2.0, label="model: sum (encoder reading)")
    ax_f.plot(f, pt._fase_avvolta(G), color=pt.INCHIOSTRO, lw=1.6)
    for k, (nome, Gk) in enumerate(mod.contributi(f).items()):
        st = {"color": pt.INCHIOSTRO_SECONDARIO, "linestyle": "-"} if k == 0 else pt._stile(k - 1)
        ax_m.plot(f, pt._db(Gk), lw=1.3, label=nome, **st)
        ax_f.plot(f, pt._fase_avvolta(Gk), lw=1.1, **st)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(*BAND)
    ax_m.set_ylim(-110, -30)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7, ncol=2)
    pt._titolo(fig, "Elevation model: encoder reading split into encoder disc and stator masses", seg)
    return fig


def fig_roots(radici: list[tuple[str, tuple]], seg: bool):
    """Poles and zeros below 150 Hz in the upper half of the s-plane, in Hz."""
    fig, assi = plt.subplots(1, len(radici), figsize=(11.0, 4.6))
    for ax, (nome, (poli, zeri)) in zip(assi, radici):
        for vals, mk, lab in ((poli, "x", "poles"), (zeri, "o", "zeros")):
            v = [(fn, z) for fn, z in vals if fn <= BAND[1]]
            x = [-z * fn for fn, z in v]
            y = [fn * np.sqrt(max(1 - z * z, 0.0)) for fn, z in v]
            ax.plot(x, y, mk, ms=8, mfc="none", mew=1.6, color=pt.INCHIOSTRO if mk == "x" else pt.SERIE[1], label=lab)
            for xi, yi, (fn, z) in zip(x, y, v):
                ax.annotate(f"{fn:.1f} Hz, ζ={z:+.2f}", (xi, yi), xytext=(6, 4), textcoords="offset points",
                            fontsize=7, color=pt.INCHIOSTRO_SECONDARIO)
        ax.axvline(0.0, color=pt.INCHIOSTRO_SECONDARIO, lw=0.9)
        ax.axvspan(0.0, 40.0, color=pt.SERIE[3], alpha=0.07, lw=0)
        ax.set_xlim(-40, 40)
        ax.set_ylim(0, 120)
        ax.set_xlabel("real part, -ζ f_n [Hz]")
        ax.set_ylabel("imaginary part [Hz]")
        ax.set_title(nome, fontsize=9)
        ax.text(20, 112, "right half-plane", ha="center", fontsize=7, color=pt.INCHIOSTRO_SECONDARIO)
    assi[0].legend(loc="lower left", fontsize=7)
    pt._titolo(fig, "Poles and zeros of the fitted models below 150 Hz (upper half of the s-plane)", seg)
    return fig


# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    s2 = json.loads((RES2 / "summary.json").read_text())
    s4 = json.loads((RES4 / "summary.json").read_text())
    assert s2["fit_band_Hz"] == list(BAND), s2["fit_band_Hz"]
    pt.applica_stile()
    par = pr.parametri_da_config(ROOT / s2["config"])
    cfg = par.sperimentale
    Kt, seg = cfg.Kt, par.segnaposto
    ml = lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0)
    cat_az, cat_pol = mf.catene_h7(ml)
    el, az_s = s2["elevation"], s2["azimuth"]
    f_modi = s2["stator_mass_starts_Hz"]
    n = len(f_modi)
    mod_p = mf.modello_polare(cat_pol, el["parameters"], n, Kt)
    mod_a = mf.modello_azimut(cat_az, az_s["parameters"], Kt)
    pol = [afm.misura(cfg, fn, "polare", et) for fn, et in afm.MISURE_POLARI]
    az = afm.misura(cfg, afm.MISURA_AZIMUT[0], "azimut", afm.MISURA_AZIMUT[1])
    M: dict[str, str] = {}

    # ---- figures
    FIG.mkdir(parents=True, exist_ok=True)
    f = np.geomspace(*BAND, 4000)
    pt.salva(fig_bode(f, pol, mod_p.frf(f), f"model: H7 + {n} stator masses", (-75, -30), (-70, -38),
                      "Elevation plant, encoder velocity / current, fit band 1-150 Hz", seg), FIG, "polar_bode")
    pt.salva(fig_bode(f, [az], mod_a.frf(f), "model: flexible load", (-95, -20), (-70, -30),
                      "Azimuth plant, encoder velocity / current, fit band 1-150 Hz", seg), FIG, "azimuth_bode")
    pt.salva(fig_contributions(f, pol[0], mod_p, seg), FIG, "polar_contributions")
    dati_pol = [mf.dato(m.etichetta, m.f, m.G, BAND) for m in pol]
    dati_az = [mf.dato(az.etichetta, az.f, az.G, BAND)]
    curve = [(d.etichetta, d.f, mod_p.frf(d.f), d.G) for d in dati_pol] + [
        (d.etichetta, d.f, mod_a.frf(d.f), d.G) for d in dati_az]
    pt.salva(afm.figura_errori(curve, BAND[1], segnaposto=seg), FIG, "fit_errors")
    pt.salva(fig_roots([("elevation: H7 + 2 stator masses", mod_p.radici(BAND[1])),
                        ("azimuth: flexible load", mod_a.radici(BAND[1]))], seg), FIG, "roots")
    plt.close("all")

    # ---- data
    rows = []
    for (fn, _), m, d in zip(afm.MISURE_POLARI + (afm.MISURA_AZIMUT,), pol + [az], dati_pol + dati_az):
        grid = 100 * np.median(np.diff(np.log(m.f)))
        rows.append(["elevation" if m.asse == "polare" else "azimuth", rf"\path{{{Path(fn).name}}}",
                     f"{m.f[0]:.0f}--{m.f[-1]:.0f}", f"{grid:.0f}", f"{m.eccitazione_A:g}", str(d.f.size)])
    write("tab_data.tex", rows_tex(rows))
    M["nPointsPol"] = str(sum(d.f.size for d in dati_pol))
    M["nPointsAz"] = str(dati_az[0].f.size)

    # ---- elevation parameters
    st_p = starts(RES2 / "tables" / "starts.csv", "elevation")
    h7_p = {"k_b": cat_pol.k_b, "c_b": cat_pol.c_b}
    if par.vincolo_terra and par.vincolo_terra.polare:
        h7_p |= {"k_g": par.vincolo_terra.polare.k, "eta_g": par.vincolo_terra.polare.eta}
    write("tab_polar_params.tex", rows_tex(param_rows(mf.parametri_polari(f_modi, cat_pol), el["parameters"],
                                                     el["std_rel"], st_p, h7_p, el["at_bound"])))
    rows = [[str(r + 1), f"{m['f_Hz']:.2f}", f"{m['zeta']:.4f}", f"{m['a']:+.3f}", f"{m['J_kgm2']:.3f}",
             f"{m['kappa']:+.0f}", num(m["k_Nm_rad"], 3), f"{m['c_Nms_rad']:.1f}"] for r, m in enumerate(el["stator_masses"])]
    write("tab_masses.tex", rows_tex(rows))
    vp = el["parameters"]
    M |= {"kb": num(vp["k_b"], 3), "cb": f"{vp['c_b']:.0f}", "kg": f"{vp['k_g']:.0f}", "etag": f"{vp['eta_g']:.2f}",
          "taudP": f"{1e3 * vp['tau_d']:.2f}", "kbHseven": num(cat_pol.k_b, 3), "cbHseven": f"{cat_pol.c_b:.2f}",
          "Jaxis": f"{el['axis_inertia_kgm2']:.3f}", "invJaxis": f"{1 / el['axis_inertia_kgm2']:.3f}",
          "costP": f"{el['cost']:.2f}", "costPground": f"{el['cost_bearings_to_ground']:.2f}",
          "costPfour": f"{s4['elevation']['cost']:.2f}", "Kt": f"{Kt:g}", "fCL": f"{mf.F_CORRENTE / 1e3:g}",
          "lossScale": f"{mf.SCALA_PERDITA:g}"}
    cp = costs(RES2 / "tables" / "starts.csv", "elevation")
    M["nStartsP"] = str(len(cp))
    M["nStartsPbest"] = str(sum(abs(c - min(cp)) < 1e-3 * min(cp) for c in cp))
    for r, m in enumerate(el["stator_masses"], 1):
        R = "I" * r
        M |= {f"fMass{R}": f"{m['f_Hz']:.1f}", f"zMass{R}": f"{m['zeta']:.3f}", f"aMass{R}": f"{m['a']:+.3f}",
              f"JMass{R}": f"{m['J_kgm2']:.2f}"}
    # residue sign and non-minimum-phase zeros: rigid rotor + two undamped modes
    J = el["axis_inertia_kgm2"]
    (w1, a1), (w2, a2) = [(2 * np.pi * m["f_Hz"], m["a"]) for m in el["stator_masses"]]
    A = 1 / J + a1 + a2
    B = (w1 ** 2 + w2 ** 2) / J + a1 * w2 ** 2 + a2 * w1 ** 2
    C = w1 ** 2 * w2 ** 2 / J
    disc = B ** 2 - 4 * A * C
    u = np.roots([A, B, C]).astype(complex)
    z = np.sqrt(u[0])
    M |= {"nmpA": f"{A:.3f}", "nmpB": num(B, 3), "nmpC": num(C, 3), "nmpDisc": num(disc, 2),
          "nmpF": f"{abs(z) / (2 * np.pi):.1f}", "nmpZeta": f"{abs(z.real) / abs(z):.2f}",
          "aOneJ": f"{a1 * J:+.2f}"}
    poli_p, zeri_p = mod_p.radici(BAND[1])
    M["zerosP"] = ", ".join(f"{fz:.1f}\\,Hz ($\\zeta={zz:+.2f}$)" for fz, zz in zeri_p)
    M["polesP"] = ", ".join(f"{fp:.1f}" for fp, _ in poli_p)

    # ---- azimuth parameters
    st_a = starts(RES2 / "tables" / "starts.csv", "azimuth")
    h7_a = {"k_3": cat_az.k[2], "c_3": cat_az.c[2], "k_2": cat_az.k[1], "c_2": cat_az.c[1]}
    if par.vincolo_terra and par.vincolo_terra.azimut:
        h7_a |= {"k_g": par.vincolo_terra.azimut.k, "eta_g": par.vincolo_terra.azimut.eta}
    write("tab_azimuth_params.tex", rows_tex(param_rows(mf.parametri_azimut(cat_az), az_s["parameters"],
                                                       az_s["std_rel"], st_a, h7_a, az_s["at_bound"])))
    rows = [[str(k + 1), f"{m['f_Hz']:.2f}", f"{m['zeta']:.4f}"] + [f"{m['shape'][q]:+.2f}" for q in mf.NODI_AZIMUT]
            for k, m in enumerate(az_s["modes"])]
    write("tab_azimuth_modes.tex", rows_tex(rows))
    va = az_s["parameters"]
    ca = sorted(costs(RES2 / "tables" / "starts.csv", "azimuth"))
    M |= {"alphaA": f"{va['alpha']:.3f}", "JL": f"{az_s['J_L_kgm2']:.3f}", "JA": f"{az_s['J_A_kgm2']:.3f}",
          "JP": f"{az_s['J_P_kgm2']:.3f}", "JLcad": f"{cat_az.J_L:.3f}", "kLA": num(va["k_LA"], 3),
          "kthree": num(va["k_3"], 3), "cthree": f"{va['c_3']:.0f}", "ktwo": num(va["k_2"], 3),
          "kgA": f"{va['k_g']:.0f}", "etaA": f"{va['eta_g']:.2f}", "taudA": f"{1e3 * va['tau_d']:.2f}",
          "kthreeHseven": num(cat_az.k[2], 3), "ktwoHseven": num(cat_az.k[1], 3), "kone": num(cat_az.k[0], 3),
          "costA": f"{az_s['cost']:.2f}", "nStartsA": str(len(ca)),
          "costsA": ", ".join(f"{c:.2f}" for c in sorted(set(round(c, 2) for c in ca))),
          "fForkAz": f"{az_s['modes'][1]['f_Hz']:.1f}", "zForkAz": f"{az_s['modes'][1]['zeta']:.3f}",
          "fRigidAz": f"{az_s['modes'][0]['f_Hz']:.1f}"}

    # ---- fit quality in the band
    def qrow(q: dict, axis: str) -> list[str]:
        g = lambda k: "--" if q[k] is None or not np.isfinite(q[k]) else neg(f"{q[k]:.1f}")  # noqa: E731
        return [axis, q["model"], esc(q["measurement"]), g("fit_all"), g("fit_1-5 Hz"), g("fit_5-100 Hz"),
                g("fit_>100 Hz"), f"{q['rms_dB']:.2f}", f"{q['rms_deg']:.1f}"]
    qs = [q for q in s2["fit_quality"] if "bearings to ground" not in q["model"]]
    rows = []
    for q in qs:
        axis = "azimuth" if q["measurement"].startswith("azimuth") else "elevation"
        rows.append(qrow(q, axis))
    write("tab_fit.tex", rows_tex(rows).replace("H7 configuration", "H7").replace(
        f"H7 + {n} stator masses", f"H7 + {n} masses"))
    by = {(q["model"], q["measurement"]): q for q in s2["fit_quality"]}
    for k, (fn, et) in enumerate(afm.MISURE_POLARI, 1):
        q = by[(f"H7 + {n} stator masses", et)]
        R = "I" * k
        M |= {f"fitP{R}": f"{q['fit_all']:.0f}", f"dbP{R}": f"{q['rms_dB']:.1f}", f"degP{R}": f"{q['rms_deg']:.1f}",
              f"fitHP{R}": neg(f"{by[('H7 configuration', et)]['fit_all']:.0f}")}
    qa = by[("flexible load", afm.MISURA_AZIMUT[1])]
    M |= {"fitA": f"{qa['fit_all']:.0f}", "dbA": f"{qa['rms_dB']:.1f}", "degA": f"{qa['rms_deg']:.1f}",
          "fitHA": neg(f"{by[('H7 configuration', afm.MISURA_AZIMUT[1])]['fit_all']:.0f}")}
    # two masses against four
    by4 = {(q["model"], q["measurement"]): q for q in s4["fit_quality"]}
    rows = []
    for fn, et in afm.MISURE_POLARI:
        q2, q4 = by[(f"H7 + {n} stator masses", et)], by4[("H7 + 4 stator masses", et)]
        g = lambda q, k: "--" if not np.isfinite(q[k]) else neg(f"{q[k]:.1f}")  # noqa: E731
        rows.append([esc(et)] + [x for k in ("fit_all", "fit_1-5 Hz", "fit_5-100 Hz", "fit_>100 Hz")
                            for x in (g(q4, k), g(q2, k))] + [f"{q4['rms_deg']:.1f}", f"{q2['rms_deg']:.1f}"])
    write("tab_compare.tex", rows_tex(rows))
    m4 = s4["elevation"]["stator_masses"]
    M["massesFour"] = ", ".join(f"{m['f_Hz']:.1f}" for m in m4)

    # ---- resonances and antiresonances
    feat: dict[str, dict[str, list[str]]] = {}
    with open(RES2 / "tables" / "features.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["kind"] in KIND:
                feat.setdefault(r["curve"], {"antiresonance": [], "resonance": []})[KIND[r["kind"]]].append(r["f_Hz"])
    rows = [[esc(c.replace("measured: ", "measured, ").replace("model: ", "model, ")), ", ".join(v["antiresonance"]) or "--",
             ", ".join(v["resonance"]) or "--"] for c, v in feat.items()]
    write("tab_features.tex", rows_tex(rows))

    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("figures, tables and numbers written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "fork_modes_report.tex"],
                       cwd=HERE, check=True)


if __name__ == "__main__":
    main()
