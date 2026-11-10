"""Builds the report on the Dahl fork model identified with the unweighted cost on the resampled frequency responses.

    ../.venv/bin/python -I build_dahl_resampled_report.py [--no-latex]

Inputs (no fit is run here), all in ``results_fork_modes_dahl_resampled``:

* ``summary.json`` (``analisi_dahl_ricampionato.py``): estimates on the two uniform grids and parameters of the
  previous report;
* ``bode_fit.csv``, ``sampling.json``, ``figures/*_bode_fit.pdf``, ``figures/bearing_describing_function_resampled.pdf``
  (``analisi_fit_bode_ricampionato.py``): fit on the Bode curves, sampling of the grids, figures;
* ``open_loop_fit.csv``, ``closed_loop_fit.csv``, ``time_domain_summary.json`` and the time-domain figures
  (``analisi_tempo_ricampionato.py``).

The doublet measurements in the peak table, the time-domain Dahl parameters and the results with the file gains come
from the results of the previous report (``results_fork_modes_dahl``). The text is in dahl_resampled_report.tex and
model_section.tex; this script copies the figures (figures/), writes tables and numbers (generated/) and compiles with
latexmk.
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
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "report_dahl_bearings"))

import build_dahl_report as old  # noqa: E402
from build_dahl_report import TEST_LABEL, TESTS, UNIT, VINCOLI, klass, num, pct, read_csv, rows_tex, std_txt, tname  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
RES = ROOT / "results_fork_modes_dahl_resampled"
COPY = {"freq_elevation": "elevation_bode_fit", "freq_azimuth": "azimuth_bode_fit",
        "df_resampled": "bearing_describing_function_resampled",
        "ol_azimuth": "azimuth_open_loop_resampled", "ol_elevation": "elevation_open_loop_resampled",
        "cl_azimuth": "azimuth_closed_loop_resampled", "cl_elevation": "elevation_closed_loop_resampled"}
VAR = {"report": "previous weights", "uniform_lnf": r"\textbf{uniform in $\ln f$}", "uniform_f": "uniform in $f$"}
BANDS = ("", "_1-5", "_5-30", "_30-150")
SHORT = {"elevation, 1 % grid (1-600 Hz)": "elevation, 1\\,\\% (1--150\\,Hz)",
         "elevation, 1 % grid (30-150 Hz)": "elevation, 1\\,\\% (30--150\\,Hz)",
         "elevation, 16 screws, 4 A (30-300 Hz)": "elevation, 16 screws, 2\\,\\%",
         "azimuth, 3 % grid (1-600 Hz)": "azimuth, 3\\,\\%"}


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def e3(x: float) -> str:
    return rf"\ensuremath{{{x * 1e4:.3f}\times10^{{-4}}}}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in FIG.glob("*"):
        f.unlink()
    for name, src in COPY.items():
        shutil.copyfile(RES / "figures" / f"{src}.pdf", FIG / f"{name}.pdf")
    s = json.loads((RES / "summary.json").read_text())
    samp = json.loads((RES / "sampling.json").read_text())
    td = json.loads((RES / "time_domain_summary.json").read_text())
    bode = read_csv(RES / "bode_fit.csv")
    ol, cl = read_csv(RES / "open_loop_fit.csv"), read_csv(RES / "closed_loop_fit.csv")
    tb = json.loads((old.TB / "summary.json").read_text())
    fd = json.loads((old.DAHL / "summary.json").read_text())
    gprev = json.loads((old.CL / "stima" / "estimated_gains.json").read_text())
    cl_prev = read_csv(old.CL / "stima" / "closed_loop_fit.csv")
    cl_akd = read_csv(old.CL / "akd" / "closed_loop_fit.csv")
    M: dict[str, str] = {}

    def par(axis, var):
        return s[axis]["report_hybrid"]["parameters"] if var == "report" else s[axis][var]["hybrid"]["parameters"]

    def qual(axis, var):
        return s[axis]["report_hybrid"]["fit_quality"] if var == "report" else s[axis][var]["hybrid"]["fit_quality"]

    # ---- sampling
    rows, total, overlap = [], 0, 0
    for axis in ("elevation", "azimuth"):
        for lab, x in samp[axis].items():
            if lab == "e2":
                continue
            pl, pf = x["points"]["uniform_lnf"], x["points"]["uniform_f"]
            rows.append([SHORT[lab], str(x["n"]), f"{x['old_weight_squared_median']:.2f}", f"{x['step_lnf_pct']:.0f}\\,\\%",
                         f"{x['step_f_Hz']:.2f}\\,Hz", f"{pl['1-5']} / {pl['5-30']} / {pl['30-150']}",
                         f"{pf['1-5']} / {pf['5-30']} / {pf['30-150']}"])
            if axis == "elevation":
                total += x["n"]; overlap += pl["30-150"]
    write("tab_sampling.tex", rows_tex(rows))
    offs = [x["lnf_max_offset_from_measured_pct"] for a in samp.values() for k, x in a.items() if k != "e2"]
    M["overlapShare"] = f"{100 * overlap / total:.0f}"
    M["gridOffset"] = rf"10^{{{int(np.ceil(np.log10(max(offs))))}}}"

    # ---- friction parameters
    rows = []
    for axis in ("azimuth", "elevation"):
        t = tb["time_domain_dahl"][axis]
        rows.append([axis, "doublets (V2, time domain)", num(t["sigma0"], 4), "--", num(t["Tc"], 3), num(t["Ts_ratio"], 3),
                     f"{1e6 * t['Tc'] / t['sigma0']:.0f}"])
        p = s[axis]["report_hybrid"]["parameters"]
        rows.append(["", "hybrid, previous weights", num(p["sigma0"], 4), std_txt(tb[axis]["hybrid"]["std_rel"]["sigma0"]),
                     num(p["Tc"], 3) + " (fixed)", "1", f"{1e6 * p['Tc'] / p['sigma0']:.0f}"])
        f = s[axis]["uniform_lnf"]["free_dahl"]["parameters"]
        rows.append(["", r"free Dahl, uniform in $\ln f$", num(f["sigma0"], 4), "--", num(f["Tc"], 3), "1",
                     f"{1e6 * f['Tc'] / f['sigma0']:.0f}"])
        p = par(axis, "uniform_lnf")
        rows.append(["", r"\textbf{hybrid, uniform in $\ln f$}", r"\textbf{" + num(p["sigma0"], 4) + "}",
                     std_txt(s[axis]["uniform_lnf"]["hybrid"]["std_rel"]["sigma0"]),
                     r"\textbf{" + num(p["Tc"], 3) + "} (fixed)", "1", f"{1e6 * p['Tc'] / p['sigma0']:.0f}"])
        p = par(axis, "uniform_f")
        rows.append(["", "hybrid, uniform in $f$", num(p["sigma0"], 4),
                     std_txt(s[axis]["uniform_f"]["hybrid"]["std_rel"]["sigma0"]), num(p["Tc"], 3) + " (fixed)", "1",
                     f"{1e6 * p['Tc'] / p['sigma0']:.0f}"])
    write("tab_friction.tex", rows_tex(rows))

    # ---- frequency summary and Bode fit
    for axis in ("elevation", "azimuth"):
        rows = []
        for var, lab in VAR.items():
            q = qual(axis, var)
            e = samp[axis]["e2"][var]
            rows.append([lab, e3(e["uniform_lnf"]), e3(e["uniform_f"]), f"{par(axis, var)['g_Kt']:.3f}",
                         " / ".join(f"{r['rms_dB']:.1f}" for r in q), " / ".join(f"{r['rms_deg']:.0f}" for r in q)])
        write(f"tab_cost_{axis}.tex", rows_tex(rows))
    rows = []
    for axis in ("elevation", "azimuth"):
        for lab in [k for k in samp[axis] if k != "e2"]:
            for var, vl in VAR.items():
                r = next(x for x in bode if x["axis"] == axis and x["parameters"] == var and x["measurement"] == lab)
                row = [SHORT[lab] if var == "report" else "", vl]
                for kind in ("complex", "magnitude", "phase"):
                    row += [pct(float(r[f"{kind}{b}"])) for b in BANDS]
                rows.append(row)
            rows[-1][-1] += r" \\ \addlinespace[2pt] %"
    write("tab_bode_fit.tex", rows_tex(rows).removesuffix(r" \\ \addlinespace[2pt] %").replace(r"\\ \addlinespace[2pt] % \\", r"\\ \addlinespace[2pt]"))

    # ---- hybrid parameter tables
    for axis in ("elevation", "azimuth"):
        h = s[axis]["uniform_lnf"]["hybrid"]
        ps = fl.parametri_polari("dahl", **VINCOLI) if axis == "elevation" else fl.parametri_azimut("dahl", **VINCOLI)
        rows = []
        for p in ps:
            if p.nome == "Tc":
                rows.append([r"$T_c$", "bearings, doublets", UNIT["N m"], num(par(axis, "report")["Tc"]),
                             num(par(axis, "uniform_lnf")["Tc"]), num(par(axis, "uniform_f")["Tc"]), "fixed", "--", ""])
                continue
            rows.append([tname(p.nome, axis), klass(p.nome), UNIT.get(p.unita, p.unita), num(par(axis, "report")[p.nome]),
                         r"\textbf{" + num(par(axis, "uniform_lnf")[p.nome]) + "}", num(par(axis, "uniform_f")[p.nome]),
                         f"[{num(p.lo, 3)}, {num(p.hi, 3)}]", std_txt(h["std_rel"].get(p.nome)),
                         "yes" if p.nome in h["at_bound"] else ""])
        write(f"tab_hybrid_{axis}.tex", rows_tex(rows))

    # ---- time domain
    def moves(rows, axis, var, col, t):
        return " / ".join(pct(float(r[col])) for r in rows if r["axis"] == axis and r["parameters"] == var and r["test"] == t)

    write("tab_openloop.tex", rows_tex([[TEST_LABEL[t], moves(ol, "azimut", "uniform_lnf", "fit_motion_pct", t),
                                         moves(ol, "azimut", "uniform_f", "fit_motion_pct", t),
                                         moves(ol, "polare", "uniform_lnf", "fit_motion_pct", t),
                                         moves(ol, "polare", "uniform_f", "fit_motion_pct", t)] for t in TESTS]))
    write("tab_closedloop.tex", rows_tex(
        [[TEST_LABEL[t]] + [moves(cl, ax, "uniform_lnf", c, t) for ax in ("azimut", "polare")
                            for c in ("fit_position_window_pct", "fit_position_motion_pct", "fit_current_motion_pct")]
         for t in TESTS]))
    med = td["medians_without_anomalous_window"]
    rows = []
    for ax, name in (("azimut", "azimuth"), ("polare", "elevation")):
        for var, lab in VAR.items():
            m = med[ax][var]
            rows.append([name if var == "report" else "", lab, pct(m["open_loop_fit_motion"], 1),
                         pct(m["closed_loop_fit_position_motion"], 1), pct(m["closed_loop_fit_current_motion"], 1)])
    write("tab_medians.tex", rows_tex(rows))
    for ax, X in (("azimut", "A"), ("polare", "P")):
        for var, V in (("report", "R"), ("uniform_lnf", "L"), ("uniform_f", "F")):
            m = med[ax][var]
            M |= {f"ol{X}{V}": pct(m["open_loop_fit_motion"], 1), f"clp{X}{V}": pct(m["closed_loop_fit_position_motion"], 1),
                  f"cli{X}{V}": pct(m["closed_loop_fit_current_motion"], 1)}

    def peaks(ax, t):
        meas = [r for r in cl_prev if r["axis"] == ax and r["test"] == t]
        out = [f"{np.median([float(r['max_following_error_meas_mdeg']) for r in meas]):.0f}"]
        sel = {v: [r for r in cl if r["axis"] == ax and r["test"] == t and r["parameters"] == v]
               for v in ("uniform_lnf", "uniform_f")}
        out += [f"{np.median([float(r['max_following_error_model_mdeg']) for r in x]):.0f}" for x in sel.values()]
        out += [f"{np.median([float(r['max_current_meas_A']) for r in meas]):.1f}"]
        out += [f"{np.median([float(r['max_current_model_A']) for r in x]):.1f}" for x in sel.values()]
        return out
    write("tab_closedloop_peaks.tex", rows_tex([[TEST_LABEL[t]] + peaks("azimut", t) + peaks("polare", t) for t in TESTS]))

    rows = []
    for ax, axis in (("azimut", "azimuth"), ("polare", "elevation")):
        g = gprev[ax]
        rows.append([axis, "files (.akd, WorkBench)", f"{g['files']['Kpp']:.2f}", f"{g['files']['Kpv']:.0f}", f"{g['files']['Kiv_hz']:.2f}"])
        rows.append(["", "\\texttt{closedLoopModel} (September)", f"{g['start_closedLoopModel']['Kpp']:.2f}",
                     f"{g['start_closedLoopModel']['Kpv']:.0f}", f"{g['start_closedLoopModel']['Kiv_hz']:.2f}"])
        for var, lab in VAR.items():
            x = td["gains"][axis][var]
            rows.append(["", f"V2, 1$^\\circ$, {lab}", f"{x['Kpp']:.2f}", f"{x['Kpv']:.0f}", f"{x['Kiv_hz']:.2f}"])
    write("tab_gains.tex", rows_tex(rows))

    # ---- numbers
    for axis, X in (("elevation", "P"), ("azimuth", "A")):
        p, r, f = par(axis, "uniform_lnf"), par(axis, "report"), par(axis, "uniform_f")
        fr = s[axis]["uniform_lnf"]["free_dahl"]["parameters"]
        M |= {f"sigmaH{X}": f"{p['sigma0']:.0f}", f"sigmaR{X}": f"{r['sigma0']:.0f}", f"sigmaF{X}": f"{f['sigma0']:.0f}",
              f"sigmaStd{X}": std_txt(s[axis]["uniform_lnf"]["hybrid"]["std_rel"]["sigma0"]),
              f"TcH{X}": f"{p['Tc']:.3f}", f"gKt{X}": f"{p['g_Kt']:.2f}", f"gKtF{X}": f"{f['g_Kt']:.2f}",
              f"presl{X}": f"{1e6 * p['Tc'] / p['sigma0']:.0f}", f"TcFree{X}": f"{fr['Tc']:.3f}",
              f"TcFreeF{X}": f"{s[axis]['uniform_f']['free_dahl']['parameters']['Tc']:.3f}",
              f"TsT{X}": f"{tb['time_domain_dahl'][axis]['Ts_ratio']:.2f}",
              f"fitH{X}": " / ".join(pct(q["fit_all"], 1) for q in qual(axis, "uniform_lnf")),
              f"fitR{X}": " / ".join(pct(q["fit_all"], 1) for q in qual(axis, "report")),
              f"fitF{X}": " / ".join(pct(q["fit_all"], 1) for q in qual(axis, "uniform_f")),
              f"lowH{X}": pct(qual(axis, "uniform_lnf")[0]["fit_1-5 Hz"], 1),
              f"lowR{X}": pct(qual(axis, "report")[0]["fit_1-5 Hz"], 1),
              f"lowF{X}": pct(qual(axis, "uniform_f")[0]["fit_1-5 Hz"], 1),
              f"Arange{X}": f"{fd[axis]['amplitude_range_urad'][0]:.0f}--{fd[axis]['amplitude_range_urad'][1]:.0f}"}
    p = par("elevation", "uniform_lnf")
    a = par("azimuth", "uniform_lnf")
    M |= {"fMassI": f"{p['f_1']:.1f}", "fMassII": f"{p['f_2']:.1f}", "zMassI": f"{p['zeta_2']:.3f}",
          "zMassII": f"{p['zeta_3']:.3f}", "aMassI": num(p["a_1"], 3), "aMassII": num(p["a_2"], 3),
          "zForkA": f"{a['zeta_2']:.3f}", "alphaA": f"{a['alpha']:.3f}", "kLA": num(a["k_LA"], 3), "Kt": "0.258"}
    for ax, axis, X in (("azimut", "azimuth", "A"), ("polare", "elevation", "P")):
        lim = td["limits"][axis]
        M[f"icmdMax{X}"] = f"{max(lim[v]['i_cmd'] for v in lim):.1f}"
        ak = [r for r in cl_akd if r["axis"] == ax]
        M |= {f"akdErr{X}": f"{max(float(r['max_following_error_model_mdeg']) for r in ak):.0f}",
              f"akdErrMeas{X}": f"{max(float(r['max_following_error_meas_mdeg']) for r in ak):.0f}",
              f"akdImax{X}": f"{max(float(r['max_current_model_A']) for r in ak):.0f}",
              f"akdImaxMeas{X}": f"{max(float(r['max_current_meas_A']) for r in ak):.1f}"}
        big = [r for r in cl if r["axis"] == ax and r["parameters"] == "uniform_lnf" and r["test"] != "01deg_minus01deg"
               and not (ax == "azimut" and r["test"] == "8deg_minus8deg" and r["move"] == "#1")]
        M |= {f"clFitMin{X}": f"{min(float(r['fit_position_window_pct']) for r in big):.1f}",
              f"clFitMax{X}": f"{max(float(r['fit_position_window_pct']) for r in big):.1f}"}
        prev = {(r["test"], r["move"]): float(r["max_current_meas_A"]) for r in cl_prev if r["axis"] == ax}
        rat = [float(r["max_current_model_A"]) / prev[(r["test"], r["move"])] for r in big if r["test"] != "05deg_minus05deg"]
        M |= {f"iRatioMed{X}": f"{np.median(rat):.2f}", f"iRatioMin{X}": f"{min(rat):.2f}", f"iRatioMax{X}": f"{max(rat):.2f}"}
    qi = {lab: x for lab, x in s["elevation"]["uniform_f"]["interpolation"].items()} | s["azimuth"]["uniform_f"]["interpolation"]
    M["quadResidual"] = " / ".join(f"{qi[k]['quadratic_vs_linear_over_fit_residual_pct']:.0f}" for k in SHORT)
    M["quadSignal"] = " / ".join(f"{qi[k]['quadratic_vs_linear_over_signal_pct']:.1f}" for k in SHORT)
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("figures, tables and numbers written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "dahl_resampled_report.tex"],
                       cwd=HERE, check=True)


if __name__ == "__main__":
    main()
