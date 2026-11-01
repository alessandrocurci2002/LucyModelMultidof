"""Builds the report on the study of the October 2026 closed-loop tests (LucyIdentification/VClosedLoop).

    ../../.venv/bin/python -I build_vcl_report.py [--no-latex]

Inputs: the summary.json files of vcl_study.py, current_input_study.py (all iterations) and identifiability.py,
and the figures they wrote (copied into figures/). The text is in vcl_report.tex; this script writes the tables
and numbers (generated/), the summary figures (figures/) and compiles with latexmk.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CL = HERE.parent
sys.path.insert(0, str(CL))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from lucy_cl.controller_tf import TF, describe_roots  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
AXES = ("azimuth", "elevation")
STUDY = CL / "results_vcl_study"
CI = [  # (name, folder, label)
    ("R1", "results_current_input", "first run"),
    ("I6", "results_current_input_I6", "I6"),
    ("I7", "results_current_input_I7", "I7"),
    ("I8", "results_current_input_I8", "I8"),
    ("I9a", "results_current_input_I9a", "I9a"),
    ("I9b", "results_current_input_I9b", "I9b"),
]
COPY = {  # figure in the report <- figure of a run
    "controller_azimuth": STUDY / "figures" / "controller_azimuth.pdf",
    "controller_elevation": STUDY / "figures" / "controller_elevation.pdf",
    "plant_fit": STUDY / "figures" / "plant_fit.pdf",
    "cl_iterations": STUDY / "figures" / "iterations.pdf",
    "cl_validation_I2": STUDY / "figures" / "validation_I2.pdf",
    "cl_validation_I5": STUDY / "figures" / "validation_I5.pdf",
    "imbalance_elevation": STUDY / "imbalance" / "imbalance_elevation.pdf",
    "ci_validation_R1": CL / "results_current_input" / "figures" / "position_validation.pdf",
    "ci_validation_I7": CL / "results_current_input_I7" / "figures" / "position_validation.pdf",
    "ci_validation_I8": CL / "results_current_input_I8" / "figures" / "position_validation.pdf",
    "ci_validation_I9a": CL / "results_current_input_I9a" / "figures" / "position_validation.pdf",
    "corr_I9a_azimuth": CL / "results_identifiability" / "figures" / "corr_results_current_input_I9a_azimuth.pdf",
    "corr_I9a_elevation": CL / "results_identifiability" / "figures" / "corr_results_current_input_I9a_elevation.pdf",
    "corr_I7_elevation": CL / "results_identifiability" / "figures" / "corr_results_current_input_I7_elevation.pdf",
}
TEX_NAME = {"scale": r"$\kappa_i$", "sigma0": r"$\sigma_0$", "Tc": r"$T_c$", "Ts_ratio": r"$T_s/T_c$", "vs": r"$v_s$",
            "d_x_mm": r"$d_x$", "b_{\\theta,h}": r"$b_{\theta,h}$", "H_hyst": r"$H$", "delta_hyst": r"$\delta$",
            "trig_a": r"$a$", "trig_b": r"$b$", "trig_c": r"$c$", "K_cable": r"$K_c$", "theta0_cable_deg": r"$\theta_0$",
            "k_b": r"$k_b$", "ring": r"$\kappa_r$", "kg_az": r"$k_{g,\psi}$", "eta_az": r"$\eta_\psi$",
            "kg_el": r"$k_{g,\theta}$", "eta_el": r"$\eta_\theta$", "Kpp": r"$K_{pp}$", "Kpv": r"$K_{pv}$",
            "Kiv": r"$K_{iv}$"}
UNIT = {"scale": "-", "sigma0": "N m/rad", "Tc": "N m", "Ts_ratio": "-", "vs": "rad/s", "d_x_mm": "mm",
        "b_{\\theta,h}": "N m s/rad", "H_hyst": "N m", "delta_hyst": "rad", "trig_a": "N m", "trig_b": "N m",
        "trig_c": "N m", "K_cable": "N m/rad", "theta0_cable_deg": "deg"}


def tname(n: str) -> str:
    if n in TEX_NAME:
        return TEX_NAME[n]
    if n.startswith("c_{"):
        return "$" + n + "$"
    return n.replace("_", r"\_")


def tt(s: str) -> str:
    return r"\texttt{" + s.replace("_", r"\_") + "}"


def rows_tex(rows: list[list[str]]) -> str:
    return " \\\\\n".join(" & ".join(r) for r in rows)


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def med(x) -> tuple[float, float, float]:
    x = np.array(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.median(x)), float(x.min()), float(x.max())


def mm(x) -> str:
    m, lo, hi = med(x)
    return f"{m:.1f} [{lo:.0f}, {hi:.0f}]"


def fmt(v: float) -> str:
    return f"{v:.4g}".replace("e+0", "e").replace("e-0", "e-")


# --------------------------------------------------------------------------- part A


def part_a(js: dict, macros: dict) -> None:
    # controllers
    rows = []
    for ax in AXES:
        c = js["controllers"][ax]["check"]
        rows.append([ax, f"{c['controller_dB']:.4f}", f"{c['controller_deg']:.3f}", f"{c['feedback_dB']:.4f}",
                     f"{c['feedback_deg']:.2f}"])
    write("tab_ctrl_check.tex", rows_tex(rows))
    for ax in AXES:
        rows = []
        for name in ("C_p", "C_v", "H_fwd", "H_AR34", "H_obs", "C_i", "C_e", "C_fb", "C_WB"):
            d = js["controllers"][ax][name]
            t = TF(np.array(d["num"]), np.array(d["den"]))
            z = "; ".join(describe_roots(t.zeros())) or "--"
            p = "; ".join(describe_roots(t.poles())) or "--"
            k = t.num[np.nonzero(t.num)[0][0]] / t.den[np.nonzero(t.den)[0][0]]
            rows.append(["$" + name.replace("_", "_{") + "}$", f"{t.order[0]}/{t.order[1]}", fmt(k),
                         z.replace("zeta", r"$\zeta$"), p.replace("zeta", r"$\zeta$")])
        write(f"tab_ctrl_{ax}.tex", rows_tex(rows))
    # plant fit
    rows = []
    for ax in ("elevation", "azimuth"):
        for st in js["plant_fit"]["stages"][ax]:
            vals = "; ".join(f"{tname(n)} {fmt(v)}" for n, v in st["values"].items())
            f = st["fits"]
            rows.append([ax, st["stage"].split(":")[0], f"{st['band'][0]:g}--{st['band'][1]:g}", vals,
                         f"{f['1-5 Hz']:.0f}", f"{f['5-100 Hz']:.0f}", f"{f['>100 Hz']:.0f}"])
    write("tab_plant_stages.tex", rows_tex(rows))
    fin = js["plant_fit"]["final"]
    macros.update(kbFit=fmt(fin["k_b"]), ringFit=f"{fin['ring']:.2f}", kgElFit=f"{fin['kg_el']:.0f}",
                  etaElFit=f"{fin['eta_el']:.2f}", kgAzFit=f"{fin['kg_az']:.0f}", etaAzFit=f"{fin['eta_az']:.2f}")
    # closed-loop iterations
    what = {"I1": r"file gains, $K_{pp}=2\pi\,$\texttt{PL.KP} (WorkBench)",
            "I2": r"file gains, $K_{pp}=$\,\texttt{PL.KP} in 1/s",
            "I3": r"$K_{pp},K_{pv},K_{iv}$ identified (grid + least squares), friction of the frequency fit",
            "I4": r"I3 + Coulomb level $T_c$ identified, started from I3",
            "I5": r"grid on the gains with $T_c$ of the moves, then gains and $T_c$ together"}
    rows = []
    for it in js["iterations"]:
        for ax in AXES:
            m = it["medians"]
            v = it["values"][ax]
            rows.append([it["name"] if ax == "azimuth" else "", what[it["name"]] if ax == "azimuth" else "", ax,
                         f"{v['Kpp']:.3g}", f"{v['Kpv']:.3g}", f"{v['Kiv']:.3g}", f"{v['Tc']:.2f}",
                         f"{m[f'{ax}/fit_error'][0]:.0f}", f"{m[f'{ax}/fit_current'][0]:.0f}",
                         "yes" if it["adequate"][ax] else "no"])
    write("tab_cl_iterations.tex", rows_tex(rows))
    last = js["iterations"][-1]
    macros.update(IfiveAzErr=f"{last['medians']['azimuth/fit_error'][0]:.0f}",
                  IfiveAzCur=f"{last['medians']['azimuth/fit_current'][0]:.0f}",
                  IfiveElErr=f"{last['medians']['elevation/fit_error'][0]:.0f}",
                  IfiveElCur=f"{last['medians']['elevation/fit_current'][0]:.0f}")


# --------------------------------------------------------------------------- part B


def part_b(macros: dict) -> dict:
    data = {}
    for key, folder, label in CI:
        js = json.loads((CL / folder / "summary.json").read_text())
        fits = {k: v for k, v in js["fits"].items()}
        if key == "I9b":  # validation with the neutral angles re-estimated by a grid (correction)
            corr = json.loads((CL / folder / "validation_corrected.json").read_text())
            for ax in AXES:
                fits[f"validation/{ax}"] = {"fitted": corr["fits"][f"validation/{ax}"]["fitted"]}
            js["theta0_validation"] = corr["theta0_validation"]
        data[key] = {"js": js, "fits": fits, "label": label}
    # fits table
    rows = []
    for key, d in data.items():
        f = d["fits"]
        rows.append([d["label"]] + [mm(f[f"{s}/{ax}"]["fitted"]) for ax in AXES for s in ("identification", "validation")])
    write("tab_ci_fits.tex", rows_tex(rows))
    # parameters table per axis
    for ax in AXES:
        names = []
        for key, d in data.items():
            for n in d["js"]["final"][ax]:
                if n not in names and n not in ("theta0_by_seq",) and not n.startswith(("c_{", "k_b", "ring")):
                    names.append(n)
        rows = []
        for n in names:
            r = [tname(n), UNIT.get(n, "-")]
            for key, d in data.items():
                fitted = d["js"].get("axes", {}).get(ax, {}).get("names", [])
                v = d["js"]["final"][ax].get(n)
                if n == "theta0_cable_deg" and key == "I9b":
                    th = d["js"]["final"][ax].get("theta0_by_seq", {}) | {k.split("/")[1]: v_ for k, v_ in d["js"]["theta0_validation"].items()
                                                                       if k.startswith(ax)}
                    r.append(r"\shortstack{" + r"\\".join(f"{x:.2f}" for x in th.values()) + "}")
                    continue
                if v is None:
                    r.append("--")
                else:
                    s = fmt(float(v))
                    r.append(r"\textbf{" + s + "}" if n in fitted or (n == "theta0_cable_deg" and "theta0@" in "".join(fitted)) else s)
            rows.append(r)
        write(f"tab_ci_params_{ax}.tex", rows_tex(rows))
    # screening of the first run
    rows = []
    s = data["R1"]["js"]["axes"]
    for ax in AXES:
        for n, rel in s[ax]["screening"].items():
            rows.append([ax, tname(n), f"{100 * rel:.1f}", "yes" if n in s[ax]["kept"] else "no"])
    write("tab_ci_screening.tex", rows_tex(rows))
    # figure: evolution of the fits
    keys = list(data)
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.4), sharey=True)
    for j, ax in enumerate(AXES):
        for o, (s_, c, lab) in enumerate((("identification", "C1", "identification (1 deg moves)"),
                                          ("validation", "C0", "validation (2 deg moves)"))):
            m = np.array([med(data[k]["fits"][f"{s_}/{ax}"]["fitted"]) for k in keys])
            x = np.arange(len(keys)) + (o - 0.5) * 0.18
            axs[j].errorbar(x, m[:, 0], yerr=[m[:, 0] - m[:, 1], m[:, 2] - m[:, 0]], fmt="o-", color=c, capsize=3,
                            ms=4, lw=1.0, label=f"{lab}: median, min-max")
        axs[j].axhline(90, color="k", ls=":", lw=0.8, label="adequacy threshold 90 %")
        axs[j].set_xticks(range(len(keys)), [data[k]["label"] for k in keys])
        axs[j].set_xlabel("run of the current-driven study")
        axs[j].set_ylabel("NRMSE fit of the encoder position [%]")
        axs[j].set_ylim(40, 101)
        axs[j].set_title(ax, fontsize=8.5)
        axs[j].legend(loc="lower left", fontsize=6.5)
        axs[j].grid(alpha=0.3)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fit_evolution.pdf", bbox_inches="tight")
    fig.savefig(FIG / "fit_evolution.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    for key, d in data.items():
        for ax in AXES:
            for s_ in ("identification", "validation"):
                macros[f"ci{key.replace('9', 'Nine').replace('6', 'Six').replace('7', 'Seven').replace('8', 'Eight').replace('1', 'One')}{ax[:2].title()}{s_[:3].title()}"] = \
                    f"{med(d['fits'][f'{s_}/{ax}']['fitted'])[0]:.1f}"
    return data


def identifiability(macros: dict) -> None:
    js = json.loads((CL / "results_identifiability" / "summary.json").read_text())
    rows = []
    for key, d in js.items():
        label, ax = key.split(" / ")
        sd = d["sd"]
        corr = np.array(d["corr"])
        for k, n in enumerate(d["names"]):
            se = f"{sd[k]:.2g} (abs.)" if n in ("theta0_cable_deg", "d_x_mm") else f"{100 * sd[k]:.0f}"
            big = [f"{tname(d['names'][j])} ({corr[k, j]:+.2f})" for j in range(len(d["names"])) if j != k and abs(corr[k, j]) > 0.9]
            lab_tex = label.replace("theta_0", r"$\theta_0$")
            rows.append([lab_tex if k == 0 else "", ax if k == 0 else "",
                         f"{d['condition']:.3g}" if k == 0 else "", f"{d['n_effective']:.0f}" if k == 0 else "",
                         tname(n), fmt(d["values"][k]), se, ", ".join(big) or "--"])
        rows.append(["\\midrule"])
    text = rows_tex(rows[:-1]).replace("\\midrule \\\\", "\\midrule")
    write("tab_identifiability.tex", text)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for dst, src in COPY.items():
        shutil.copyfile(src, FIG / f"{dst}.pdf")
    macros: dict = {}
    js = json.loads((STUDY / "summary.json").read_text())
    part_a(js, macros)
    part_b(macros)
    identifiability(macros)
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in macros.items()))
    print("tables, numbers and figures written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "vcl_report.tex"], cwd=HERE,
                       check=True)


if __name__ == "__main__":
    main()
