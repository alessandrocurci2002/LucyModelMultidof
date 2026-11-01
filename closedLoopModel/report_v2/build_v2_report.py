"""Builds the report on the identification on the doublets of VClosedLoop/V2 and on the resisting torque.

    ../../.venv/bin/python -I build_v2_report.py [--no-latex]

Inputs: the summary.json of v2_identification.py and of cable_torque_analysis.py, and the figures they wrote
(copied into figures/). The text is in v2_report.tex; this script writes the tables and numbers (generated/)
and compiles with latexmk.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CL = HERE.parent
FIG, GEN = HERE / "figures", HERE / "generated"
AXES = ("azimuth", "elevation")
ID = CL / "results_v2_identification"
CAB = ID / "cables"
COPY = {  # figure in the report <- figure of a run
    "model_free": ID / "figures" / "model_free.pdf",
    "cross_validation": ID / "figures" / "cross_validation.pdf",
    "fit_vs_amplitude": ID / "figures" / "fit_vs_amplitude.pdf",
    "position_identification": ID / "figures" / "position_identification.pdf",
    "position_validation": ID / "figures" / "position_validation.pdf",
    "session_timeline": CAB / "figures" / "session_timeline.pdf",
    "cruise_per_file": CAB / "figures" / "cruise_per_file.pdf",
    "holding_vs_amplitude": CAB / "figures" / "holding_vs_amplitude.pdf",
    "relaxation": CAB / "figures" / "relaxation.pdf",
    "relaxation_vs_integral_action": CAB / "figures" / "relaxation_vs_integral_action.pdf",
    "theta0_per_file": CAB / "figures" / "theta0_per_file.pdf",
    "friction_model": CAB / "figures" / "friction_model.pdf",
}
TEX_NAME = {"scale": r"$\kappa_i$", "sigma0": r"$\sigma_0$", "Tc": r"$T_c$", "Ts_ratio": r"$T_s/T_c$", "vs": r"$v_s$",
            "b_{\\theta,h}": r"$b_{\theta,h}$", "K_cable": r"$K_c$", "theta0_cable_deg": r"$\theta_0$"}
UNIT = {"scale": "--", "sigma0": "N\\,m/rad", "Tc": "N\\,m", "Ts_ratio": "--", "vs": "rad/s", "b_{\\theta,h}": "N\\,m\\,s/rad",
        "K_cable": "N\\,m/rad", "theta0_cable_deg": "deg"}
ORDER = ["scale", "Tc", "Ts_ratio", "vs", "sigma0", "b_{\\theta,h}", "K_cable", "theta0_cable_deg"]


def tname(n: str) -> str:
    return TEX_NAME.get(n, n.replace("_", r"\_"))


def short(name: str) -> str:
    return name.replace("V2/", "").replace("MPCLucy_doublets_", "").replace("MPCLucy_", "").replace("_", r"\_")


def rows_tex(rows: list[list[str]]) -> str:
    return " \\\\\n".join(" & ".join(r) for r in rows)


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def mm(x) -> str:
    x = np.array(list(x), dtype=float)
    return f"{np.median(x):.1f} [{x.min():.0f}, {x.max():.0f}]"


def struct_name(k: str) -> str:
    """Readable name of a structure key of v2_identification.py."""
    k = k.replace("S0 no spring", "no cable spring, all free").replace("S1 spring", "cable spring, all free")
    k = k.replace("S1 ", "cable spring, ").replace("S0 ", "no spring, ")
    for n in ("sigma0", "vs", "b_{\\theta,h}"):
        k = k.replace(n, tname(n))
    return k.replace(" x0.333", r" $\times\frac13$").replace(" x1", r" $\times1$").replace(" x3", r" $\times3$")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for dst, src in COPY.items():
        shutil.copyfile(src, FIG / f"{dst}.pdf")
    js = json.loads((ID / "summary.json").read_text())
    cab = json.loads((CAB / "summary.json").read_text())
    final = js["final"]
    M: dict[str, str] = {}

    # ---- data: the V2 session
    amp = js["windows"]["amplitude_deg"]
    rows = []
    for s in cab["session"]:
        n = s["file"]
        ax = s["axis"][0]
        wins = [k for k in amp if k.startswith(n + " ")]
        a = sorted({abs(round(amp[k], 2)) for k in wins})
        excluded = n in ("V2/MPCLucy_doublets_azimut_1deg_15deg_v2",)
        rows.append([short(n), ax, f"{s['start_min']:.1f}", f"{s['duration_s']:.0f}",
                     ", ".join(f"{x:g}" for x in a) or "--", str(len(wins)), "no" if excluded else "yes"])
    write("tab_session.tex", rows_tex(rows))

    # ---- model-free analysis
    rows = []
    for ax in AXES:
        cr, rf = js["model_free"][ax]["cruise"], js["model_free"][ax]["rest_fit"]
        rows.append([ax, f"{cr['friction']:.3f}", f"{cr['gradient']:.2f}", f"{cr['offset']:+.3f}", f"{cr['rms']:.2f}",
                     f"{rf['direction']:.3f}", f"{rf['gradient']:.2f}", f"{rf['rms']:.2f}"])
        M[f"F{ax[:2].title()}"] = f"{cr['friction']:.2f}"
        M[f"Kcr{ax[:2].title()}"] = f"{cr['gradient']:.2f}"
        M[f"Krest{ax[:2].title()}"] = f"{rf['gradient']:.2f}"
        M[f"Hrest{ax[:2].title()}"] = f"{rf['direction']:.2f}"
    write("tab_model_free.tex", rows_tex(rows))
    rows = []
    for ax, d in cab["cruise"].items():
        for n, c in d.items():
            th0 = -np.rad2deg(c["offset"] / c["gradient"]) if abs(c["gradient"]) > 0.3 else float("nan")
            rows.append([ax, short(n), f"{c['friction']:.3f}", f"{c['gradient']:.2f}", f"{c['offset']:+.3f}",
                         f"{th0:+.1f}" if np.isfinite(th0) else "--"])
    for ax, c in cab["cruise_all"].items():
        rows.append([ax, r"\textbf{all files}", f"{c['friction']:.3f}", f"{c['gradient']:.2f}", f"{c['offset']:+.3f}",
                     f"{-np.rad2deg(c['offset'] / c['gradient']):+.1f}"])
        A = ax[:2].title()
        M[f"Fenc{A}"] = f"{c['friction']:.2f}"
        M[f"Kenc{A}"] = f"{c['gradient']:.2f}"
        M[f"cenc{A}"] = f"{c['offset']:+.3f}"
    write("tab_cruise_files.tex", rows_tex(rows))
    rows = []
    for ax, d in cab["sessions"].items():
        for lab, v in d.items():
            rows.append([ax, lab, str(v["rests"]), f"{v['c']:+.3f}", f"{v['h']:.3f}", f"{v['rms']:.3f}"])
        A = ax[:2].title()
        M[f"cSesVtwo{A}"] = f"{d['V2']['c']:+.3f}"
        M[f"cSesOther{A}"] = f"{d['other session']['c']:+.3f}"
        M[f"cSesShift{A}"] = f"{abs(d['V2']['c'] - d['other session']['c']):.2f}"
        M[f"hSesVtwo{A}"] = f"{d['V2']['h']:.2f}"
        M[f"hSesOther{A}"] = f"{d['other session']['h']:.2f}"
    write("tab_sessions.tex", rows_tex(rows))
    for ax in AXES:
        d = cab["cruise"].get(ax, {})
        if d:
            F = [c["friction"] for c in d.values()]
            K = [c["gradient"] for c in d.values()]
            M[f"Fspread{ax[:2].title()}"] = f"{min(F):.2f}--{max(F):.2f}"
            M[f"Kspread{ax[:2].title()}"] = f"{min(K):.2f}--{max(K):.2f}"

    # ---- surrogate
    rows = []
    for ax in AXES:
        s = js["surrogate"][ax]
        rows.append([ax, f"{np.median(s['fit_difference']):.2f}", f"{max(s['fit_difference']):.2f}",
                     f"{100 * max(s['position_difference_rel']):.2f}", f"{s['wall_s']['rigid']:.0f}", f"{s['wall_s']['full']:.0f}"])
    write("tab_surrogate.tex", rows_tex(rows))

    # ---- structures
    rows = []
    for ax in AXES:
        sel = js["selection"][ax]["selected"]
        for k, r in js["structures"][ax].items():
            name = struct_name(k)
            if k == sel:
                name = r"\textbf{" + name + "}"
            rows.append([ax, name, str(len(r["free"])), f"{r['fit_median']:.1f}", f"{r['cv_median']:.1f}",
                         f"{min(r['cv_fits'].values()):.0f}", f"{r['cv_mean_square']:.2f}"])
        rows.append(["\\midrule"])
    write("tab_structures.tex", rows_tex(rows[:-1]).replace("\\midrule \\\\", "\\midrule"))
    for ax in AXES:
        S = js["structures"][ax]
        A = ax[:2].title()
        M[f"cvNoSpring{A}"] = f"{S['S0 no spring']['cv_median']:.1f}"
        M[f"cvSpring{A}"] = f"{S['S1 spring']['cv_median']:.1f}"
        M[f"cvSel{A}"] = f"{S[js['selection'][ax]['selected']]['cv_median']:.1f}"
        M[f"idSel{A}"] = f"{S[js['selection'][ax]['selected']]['fit_median']:.1f}"
        M[f"sel{A}"] = struct_name(js["selection"][ax]["selected"])
        M[f"nStruct{A}"] = str(len(S))
        gap = max(abs(r["fit_median"] - r["cv_median"]) for k, r in S.items() if not k.startswith("S0"))
        M[f"cvGap{A}"] = f"{gap:.1f}"
        M[f"cvGapNoSpring{A}"] = f"{abs(S['S0 no spring']['fit_median'] - S['S0 no spring']['cv_median']):.1f}"
    M["nominalSigmaAz"] = f"{js['nominal']['azimuth']['sigma0']:.0f}"
    M["nominalSigmaEl"] = f"{js['nominal']['elevation']['sigma0']:.0f}"

    # ---- parameters
    rows = []
    for ax in AXES:
        sel = js["structures"][ax][js["selection"][ax]["selected"]]
        for n in ORDER:
            if n not in final[ax] or (n not in sel["free"] and n not in sel["fixed"]):
                continue
            how = "fitted" if n in sel["free"] else "fixed"
            se = js["full"][ax]["se"].get(n)
            jk = "--"
            if n in sel["free"]:
                x = np.array([v[n] for v in sel["cv_values"].values()])
                x = x if n == "theta0_cable_deg" else np.log(x)
                e = float(np.sqrt((x.size - 1) / x.size * np.sum((x - x.mean()) ** 2)))
                jk = f"{e:.2g}$^\\circ$" if n == "theta0_cable_deg" else f"{100 * e:.0f}\\,\\%"
            ses = "--" if se is None else (f"{se:.2g}$^\\circ$" if n == "theta0_cable_deg" else f"{100 * se:.0f}\\,\\%")
            if se is not None and se > 100:
                ses = "$>10^4$\\,\\%"
            rows.append([ax, tname(n), UNIT[n], f"{final[ax][n]:.4g}", how, ses, jk])
        rows.append(["\\midrule"])
    write("tab_params.tex", rows_tex(rows[:-1]).replace("\\midrule \\\\", "\\midrule"))
    for ax in AXES:
        A = ax[:2].title()
        f = final[ax]
        M[f"Tc{A}"] = f"{f['Tc']:.2f}"
        M[f"Ts{A}"] = f"{f['Tc'] * f['Ts_ratio']:.2f}"
        M[f"TsRatio{A}"] = f"{f['Ts_ratio']:.2f}"
        M[f"vs{A}"] = f"{f['vs']:.4f}"
        M[f"vsDeg{A}"] = f"{np.rad2deg(f['vs']):.2f}"
        M[f"sigma{A}"] = f"{f['sigma0']:.0f}"
        M[f"Kc{A}"] = f"{f['K_cable']:.2f}"
        M[f"thz{A}"] = f"{f['theta0_cable_deg']:+.2f}"
        M[f"scale{A}"] = f"{f['scale']:.2f}"
        M[f"pre{A}"] = f"{1e3 * np.rad2deg(f['Tc'] / f['sigma0']):.0f}"
    M["bEl"] = f"{final['elevation'].get('b_{\\theta,h}', 0.0):.1e}"

    # ---- fits
    F = js["fits"]
    rows = []
    for lab, key in (("identification, V2 (49 moves)", "identification"), ("validation, $\\theta_0$ of V2", "validation"),
                     ("validation, $\\theta_0$ of the session from one other file", "validation_theta0_one_file"),
                     ("validation, $\\theta_0$ of each file", "validation_theta0_per_sequence")):
        r = [lab]
        for ax in AXES:
            v = [x for k, x in F[key].items() if f" {ax} " in k]
            r += [str(len(v)), mm(v)]
        rows.append(r)
    write("tab_fits.tex", rows_tex(rows))
    rows = []
    for lo, hi, lab in ((0.0, 0.75, "0.1--0.5"), (0.75, 10.5, "1--10"), (10.5, 99.0, "13--30")):
        r = [lab]
        for ax in AXES:
            v = [x for k, x in F["identification"].items() if f" {ax} " in k and lo <= abs(amp[k]) < hi]
            r += [str(len(v)), mm(v) if v else "--"]
        rows.append(r)
    write("tab_fit_classes.tex", rows_tex(rows))
    for ax in AXES:
        A = ax[:2].title()
        for key, mac in (("identification", "Id"), ("validation", "Val"), ("validation_theta0_one_file", "ValSes"),
                         ("validation_theta0_per_sequence", "ValSeq")):
            v = [x for k, x in F[key].items() if f" {ax} " in k]
            M[f"fit{mac}{A}"] = f"{np.median(v):.1f}"

    # ---- resisting torque: rests at 0 in recording order, relaxation
    rests = cab["rests"]
    rows = []
    for ax in AXES:
        R = [r for r in rests if r["axis"] == ax and abs(r["angle_deg"]) < 0.05]
        for r in R:
            kind = "before the doublet" if r["kind"] == "before the first move" else "after the doublet"
            last = "--" if r["last_move_deg"] is None else f"{r['last_move_deg']:+.2f}"
            rows.append([ax, short(r["file"]), f"{r['t_session_min']:.1f}", kind, last, f"{r['rest_s']:.0f}",
                         f"{r['T_start']:+.3f}", f"{r['T_end']:+.3f}"])
        rows.append(["\\midrule"])
    write("tab_return0.tex", rows_tex(rows[:-1]).replace("\\midrule \\\\", "\\midrule"))
    rows = []
    for ax in AXES:
        R = [r for r in rests if r["axis"] == ax and r["kind"] == "after a move" and r["rest_s"] >= 10.0]
        dT = np.array([r["T_end"] - r["T_start"] for r in R])
        rel = np.array([(r["T_end"] - r["T_start"]) / r["T_start"] for r in R if abs(r["T_start"]) > 0.1])
        tau = np.array([r["tau_s"] for r in R if "tau_s" in r and r["tau_s"] < 400])
        enc = np.array([abs(r["encoder_change_mdeg"]) for r in R])
        ia = np.array([r["integral_action_dT"] for r in R])
        same = np.mean(np.sign(ia) == np.sign(dT)) * 100
        rows.append([ax, str(len(R)), f"{np.median(np.abs(dT)):.3f} [{np.abs(dT).max():.3f}]",
                     f"{100 * np.median(np.abs(rel)):.0f}", f"{np.median(tau):.0f} [{tau.min():.0f}, {tau.max():.0f}]",
                     f"{np.median(enc):.2f} [{enc.max():.2f}]", f"{np.median(np.abs(ia)):.3f}", f"{same:.0f}"])
        A = ax[:2].title()
        M[f"relaxN{A}"] = str(len(R))
        M[f"relaxMax{A}"] = f"{np.abs(dT).max():.2f}"
        M[f"relaxMed{A}"] = f"{np.median(np.abs(dT)):.2f}"
        M[f"relaxRel{A}"] = f"{100 * np.median(np.abs(rel)):.0f}"
        M[f"tau{A}"] = f"{np.median(tau):.0f}"
        M[f"encMax{A}"] = f"{enc.max():.2f}"
        M[f"sameSign{A}"] = f"{same:.0f}"
        M[f"intAct{A}"] = f"{np.median(np.abs(ia)):.2f}"
    write("tab_relax.tex", rows_tex(rows))

    # ---- neutral angle file by file
    rows = []
    for ax, d in cab["theta0_per_file"].items():
        for n, v in d.items():
            rows.append([ax, short(n), f"{v['theta0_deg']:+.2f}", f"{v['fit_median']:.1f}", f"{v['fit_median_session_theta0']:.1f}"])
        for n, v in cab["theta0_other_session"].get(ax, {}).items():
            rows.append([ax, short(n) + " (other session)", f"{v:+.2f}", "--", "--"])
        rows.append(["\\midrule"])
    write("tab_theta0.tex", rows_tex(rows[:-1]).replace("\\midrule \\\\", "\\midrule"))
    for ax, d in cab["theta0_per_file"].items():
        A = ax[:2].title()
        x = np.array([v["theta0_deg"] for v in d.values()])
        M[f"thzFileMin{A}"] = f"{x.min():+.1f}"
        M[f"thzFileMax{A}"] = f"{x.max():+.1f}"
        M[f"thzFileSd{A}"] = f"{x.std(ddof=1):.1f}"
        o = np.array(list(cab["theta0_other_session"].get(ax, {}).values()))
        if o.size:
            M[f"thzOtherMin{A}"] = f"{o.min():+.1f}"
            M[f"thzOtherMax{A}"] = f"{o.max():+.1f}"
            M[f"thzShiftTorque{A}"] = f"{final[ax]['K_cable'] * np.deg2rad(abs(np.mean(o) - final[ax]['theta0_cable_deg'])):.2f}"
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("tables, numbers and figures written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "v2_report.tex"], cwd=HERE,
                       check=True)


if __name__ == "__main__":
    main()
