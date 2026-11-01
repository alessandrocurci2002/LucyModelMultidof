"""Builds the report on the fork model with Dahl bearing friction (hybrid identification, frequency and time domain).

    ../.venv/bin/python -I build_dahl_report.py [--no-latex]

Inputs (no fit is run here):

* ``results_fork_modes_dahl/time_domain_bearings`` (``analisi_cuscinetti_dahl_tempo.py``): frequency-domain
  comparison of the bearing models and the hybrid identification (T_c from the doublets, sigma0 from the FRF);
* ``results_fork_modes_dahl`` and ``results_fork_modes_lugre`` (``analisi_forcella_lugre.py``): Dahl and LuGre
  identified in frequency;
* ``saved_fork_modes_modal_20261031_cad20_kt10``: linear hysteretic bearings (previous model);
* ``results_fork_modes_dahl/time_domain_check`` (``analisi_tempo_dahl_ibrido.py``): open loop, measured current;
* ``results_fork_modes_dahl/closed_loop_check`` (``analisi_anello_chiuso_ibrido.py``): closed loop, measured position
  command, effective drive gains (``stima``) and file gains (``akd``).

The text is in dahl_bearings_report.tex; this script copies the figures (figures/), writes tables and numbers
(generated/) and compiles with latexmk.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from gimbal_mdof import forcella_lugre as fl  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
DAHL = ROOT / "results_fork_modes_dahl"
TB = DAHL / "time_domain_bearings"
OL = DAHL / "time_domain_check"
CL = DAHL / "closed_loop_check"
LUGRE = ROOT / "results_fork_modes_lugre"
LIN = ROOT / "saved_fork_modes_modal_20261031_cad20_kt10" / "results"
VINCOLI = {"zeta_max": 0.05, "cad": 0.20, "cad_J": 0.20, "kt": 0.10}
COPY = {"freq_elevation": TB / "figures" / "elevation_bode_time_domain_dahl.pdf",
        "freq_azimuth": TB / "figures" / "azimuth_bode_time_domain_dahl.pdf",
        "df_time_vs_freq": TB / "figures" / "bearing_describing_function_time_vs_frequency.pdf",
        "df_amplitude": DAHL / "figures" / "bearing_describing_function.pdf",
        "ol_azimuth": OL / "figures" / "azimuth_time_domain_hybrid.pdf",
        "ol_elevation": OL / "figures" / "elevation_time_domain_hybrid.pdf",
        "cl_azimuth": CL / "stima" / "figures" / "azimuth_closed_loop_hybrid.pdf",
        "cl_elevation": CL / "stima" / "figures" / "elevation_closed_loop_hybrid.pdf"}
MODELS = [("linear hysteretic $k_g(1+j\\eta_g)$ (previous model)", "linear"),
          ("LuGre, all 6 parameters from the FRF", "lugre"),
          ("Dahl, $\\sigma_0$ and $T_c$ from the FRF", "Dahl identified in frequency (sigma0, T_c fitted)"),
          ("\\textbf{Dahl hybrid}: $T_c$ doublets, $\\sigma_0$ FRF", "hybrid Dahl: T_c from time domain, sigma0 from FRF"),
          ("Dahl from the doublets, same structure", "Dahl from time domain, same structure"),
          ("Dahl from the doublets, structure refitted", "Dahl from time domain, structure refitted")]
TESTS = ["01deg_minus01deg", "05deg_minus05deg", "1deg_minus1deg", "3deg_minus3deg", "5deg_minus5deg",
         "8deg_minus8deg"]
TEST_LABEL = {"01deg_minus01deg": "0.1", "05deg_minus05deg": "0.5", "1deg_minus1deg": "1", "3deg_minus3deg": "3",
              "5deg_minus5deg": "5", "8deg_minus8deg": "8"}
TEX_NAME = {"s_J_e": r"$s_{J_e}$", "s_J_hp": r"$s_{J_{hp}}$", "s_J_m": r"$s_{J_m}$", "s_J_r": r"$s_{J_\mathrm{rear}}$",
            "s_J_L": r"$s_{J_L}$", "tau_d": r"$\tau_d$", "alpha": r"$\alpha$", "k_LA": r"$k_{LA}$",
            "g_Kt": r"$g_{K_t}$", "sigma0": r"$\sigma_0$"}
UNIT = {"N m/rad": r"N\,m/rad", "s": "s", "Hz": "Hz", "1/(kg m^2)": r"(kg\,m$^2$)$^{-1}$", "": "--", "N m": r"N\,m"}


def num(x: float, d: int = 4) -> str:
    if x is None or not np.isfinite(x):
        return "--"
    s = f"{x:.{d}g}"
    if "e" in s:
        m, e = s.split("e")
        return rf"\ensuremath{{{m}\times10^{{{int(e)}}}}}"
    return s.replace("-", "$-$", 1) if s.startswith("-") else s


def pct(x: float | None, d: int = 0) -> str:
    if x is None or not np.isfinite(x):
        return "--"
    s = f"{x:.{d}f}"
    return s.replace("-", "$-$", 1) if s.startswith("-") else s


def std_txt(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "--"
    if not np.isfinite(v) or 100 * v > 1e3:
        return "n.d."
    return f"{100 * v:.1f}\\,\\%"


def tname(n: str, axis: str) -> str:
    if n.startswith("s_k_"):
        return rf"$s_{{k_{{\theta{n[-1]}}}}}$" if axis == "elevation" else rf"$s_{{k_{n[-1]}}}$"
    if n in TEX_NAME:
        return TEX_NAME[n]
    base, _, r = n.rpartition("_")
    return {"f": f"$f_{r}$", "zeta": rf"$\zeta_{r}$", "a": f"$a_{r}$"}[base]


def klass(n: str) -> str:
    for k, v in (("s_J", r"CAD, $\pm20$\,\%"), ("s_k", r"CAD, $\pm20$\,\%"), ("g_Kt", r"motor, $\pm10$\,\%"),
                 ("zeta", r"modal, $\leq0.05$"), ("sigma0", "bearings, free"), ("k_LA", "introduced, free"),
                 ("alpha", "introduced, free"), ("f_", "stator mass, free"), ("a_", "stator mass, free"),
                 ("tau_d", "drive")):
        if n.startswith(k):
            return v
    return ""


def rng(v: list[float], d: int) -> str:
    """min--max of a list (one value if all equal at the shown precision)."""
    lo, hi = f"{min(v):.{d}f}", f"{max(v):.{d}f}"
    return lo if lo == hi else f"{lo}--{hi}"


def rows_tex(rows: list[list[str]]) -> str:
    return " \\\\\n".join(" & ".join(r) for r in rows)


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for old in FIG.glob("*"):
        old.unlink()
    for name, src in COPY.items():
        shutil.copyfile(src, FIG / f"{name}.pdf")
    tb = json.loads((TB / "summary.json").read_text())
    fd = json.loads((DAHL / "summary.json").read_text())
    lg = json.loads((LUGRE / "summary.json").read_text())
    lin = json.loads((LIN / "summary.json").read_text())
    gains = json.loads((CL / "stima" / "estimated_gains.json").read_text())
    lim = json.loads((CL / "stima" / "limits.json").read_text())
    ol = read_csv(OL / "position_fit.csv")
    cl = read_csv(CL / "stima" / "closed_loop_fit.csv")
    cl_akd = read_csv(CL / "akd" / "closed_loop_fit.csv")
    M: dict[str, str] = {}
    AX = {"elevation": "P", "azimuth": "A"}

    # ---- frequency-domain comparison of the bearing models
    for axis, X in AX.items():
        rows = []
        for label, key in MODELS:
            if key == "linear":
                cost, q = lin[axis]["cost"], [r for r in lin["fit_quality"]
                                              if r["measurement"].startswith(axis)]
            elif key == "lugre":
                cost, q = lg[axis]["cost"], lg[axis]["fit_quality"]
            else:
                cost, q = tb[axis]["cost"][key], tb[axis]["fit_quality"][key]
            fits = " / ".join(pct(r["fit_all"]) for r in q)
            low = pct(q[0]["fit_1-5 Hz"])
            dbs, degs = rng([r["rms_dB"] for r in q], 1), rng([r["rms_deg"] for r in q], 0)
            rows.append([label, num(cost, 3), fits, low, dbs, degs])
        write(f"tab_freq_{axis}.tex", rows_tex(rows))

    # ---- friction parameters
    rows = []
    for axis in ("azimuth", "elevation"):
        td, h, f = tb["time_domain_dahl"][axis], tb[axis]["hybrid"], fd[axis]
        name = "azimuth" if axis == "azimuth" else "elevation"
        rows.append([name, "doublets (V2, time domain)", num(td["sigma0"], 4), "--", num(td["Tc"], 3), "--",
                     num(td["Ts_ratio"], 3), f"{1e6 * td['Tc'] / td['sigma0']:.0f}"])
        rows.append(["", "FRF (frequency domain)", num(f["parameters"]["sigma0"], 4), std_txt(f["std_rel"]["sigma0"]),
                     num(f["parameters"]["Tc"], 3), std_txt(f["std_rel"]["Tc"]), "1",
                     f"{f['presliding_scale_urad']:.0f}"])
        rows.append(["", r"\textbf{hybrid}", r"\textbf{" + num(h["parameters"]["sigma0"], 4) + "}",
                     std_txt(h["std_rel"]["sigma0"]), r"\textbf{" + num(h["Tc_from_time"], 3) + "} (fixed)", "--", "1",
                     f"{h['presliding_scale_urad']:.0f}"])
    write("tab_friction.tex", rows_tex(rows))

    # ---- hybrid parameter tables
    for axis in ("elevation", "azimuth"):
        h = tb[axis]["hybrid"]
        par = fl.parametri_polari("dahl", **VINCOLI) if axis == "elevation" else fl.parametri_azimut("dahl", **VINCOLI)
        rows = []
        for p in par:
            if p.nome == "Tc":
                rows.append([r"$T_c$", "bearings, doublets", UNIT["N m"], num(h["Tc_from_time"]), "fixed", "--", ""])
                continue
            rows.append([tname(p.nome, axis), klass(p.nome), UNIT.get(p.unita, p.unita), num(h["parameters"][p.nome]),
                         f"[{num(p.lo, 3)}, {num(p.hi, 3)}]", std_txt(h["std_rel"].get(p.nome)),
                         "yes" if p.nome in h["at_bound"] else ""])
        write(f"tab_hybrid_{axis}.tex", rows_tex(rows))

    # ---- open loop and closed loop tables (three moves per test)
    def per_test(rows, axis_key, col, fmt=pct):
        out = []
        for t in TESTS:
            v = [float(r[col]) for r in rows if r["axis"] == axis_key and r["test"] == t]
            out.append(" / ".join(fmt(x) for x in v))
        return out

    rows = []
    for t, a, p in zip(TESTS, per_test(ol, "azimut", "fit_motion_pct"), per_test(ol, "polare", "fit_motion_pct")):
        rows.append([TEST_LABEL[t], a, p])
    write("tab_openloop.tex", rows_tex(rows))
    rows = []
    for t, aw, am, ai, pw, pm, pi in zip(TESTS, per_test(cl, "azimut", "fit_position_window_pct"),
                                         per_test(cl, "azimut", "fit_position_motion_pct"),
                                         per_test(cl, "azimut", "fit_current_motion_pct"),
                                         per_test(cl, "polare", "fit_position_window_pct"),
                                         per_test(cl, "polare", "fit_position_motion_pct"),
                                         per_test(cl, "polare", "fit_current_motion_pct")):
        rows.append([TEST_LABEL[t], aw, am, ai, pw, pm, pi])
    write("tab_closedloop.tex", rows_tex(rows))

    def peaks(rows, axis_key):
        out = []
        for t in TESTS:
            sel = [r for r in rows if r["axis"] == axis_key and r["test"] == t]
            em = np.median([float(r["max_following_error_meas_mdeg"]) for r in sel])
            es = np.median([float(r["max_following_error_model_mdeg"]) for r in sel])
            im = np.median([float(r["max_current_meas_A"]) for r in sel])
            is_ = np.median([float(r["max_current_model_A"]) for r in sel])
            out.append([f"{em:.0f}", f"{es:.0f}", f"{im:.1f}", f"{is_:.1f}"])
        return out
    rows = [[TEST_LABEL[t]] + a + p for t, a, p in zip(TESTS, peaks(cl, "azimut"), peaks(cl, "polare"))]
    write("tab_closedloop_peaks.tex", rows_tex(rows))

    rows = []
    for ax, name in (("azimut", "azimuth"), ("polare", "elevation")):
        g = gains[ax]
        rows.append([name, f"{g['files']['Kpp']:.2f}", f"{g['files']['Kpv']:.0f}", f"{g['files']['Kiv_hz']:.2f}",
                     f"{g['start_closedLoopModel']['Kpp']:.2f}", f"{g['start_closedLoopModel']['Kpv']:.0f}",
                     f"{g['start_closedLoopModel']['Kiv_hz']:.2f}",
                     r"\textbf{" + f"{g['estimated']['Kpp']:.2f}" + "}", r"\textbf{" + f"{g['estimated']['Kpv']:.0f}" + "}",
                     r"\textbf{" + f"{g['estimated']['Kiv_hz']:.2f}" + "}"])
    write("tab_gains.tex", rows_tex(rows))

    # ---- numbers
    for axis, X in AX.items():
        h = tb[axis]["hybrid"]
        hp = h["parameters"]
        q = tb[axis]["fit_quality"]["hybrid Dahl: T_c from time domain, sigma0 from FRF"]
        M |= {f"sigmaH{X}": f"{hp['sigma0']:.0f}", f"sigmaStd{X}": std_txt(h["std_rel"]["sigma0"]),
              f"TcH{X}": f"{h['Tc_from_time']:.3f}", f"gKt{X}": f"{hp['g_Kt']:.2f}",
              f"presl{X}": f"{h['presliding_scale_urad']:.0f}", f"costH{X}": num(h["cost"], 3),
              f"costLin{X}": num(lin[axis]["cost"], 3),
              f"costRed{X}": f"{100 * (1 - h['cost'] / lin[axis]['cost']):.0f}",
              f"sigmaF{X}": f"{fd[axis]['parameters']['sigma0']:.0f}", f"TcF{X}": f"{fd[axis]['parameters']['Tc']:.3f}",
              f"sigmaT{X}": f"{tb['time_domain_dahl'][axis]['sigma0']:.0f}",
              f"TcT{X}": f"{tb['time_domain_dahl'][axis]['Tc']:.3f}",
              f"TsT{X}": f"{tb['time_domain_dahl'][axis]['Ts_ratio']:.2f}",
              f"fitH{X}": " / ".join(pct(r["fit_all"]) for r in q),
              f"dbH{X}": rng([r["rms_dB"] for r in q], 1),
              f"Arange{X}": f"{fd[axis]['amplitude_range_urad'][0]:.0f}--{fd[axis]['amplitude_range_urad'][1]:.0f}"}
        ql = [r for r in lin["fit_quality"] if r["measurement"].startswith(axis)]
        M[f"fitLin{X}"] = " / ".join(pct(r["fit_all"]) for r in ql)
    hpP, hpA = tb["elevation"]["hybrid"]["parameters"], tb["azimuth"]["hybrid"]["parameters"]
    M |= {"fMassI": f"{hpP['f_1']:.1f}", "fMassII": f"{hpP['f_2']:.1f}", "zMassI": f"{hpP['zeta_2']:.3f}",
          "zMassII": f"{hpP['zeta_3']:.3f}", "aMassI": num(hpP["a_1"], 3), "aMassII": num(hpP["a_2"], 3),
          "zForkA": f"{hpA['zeta_2']:.3f}", "alphaA": f"{hpA['alpha']:.3f}", "kLA": num(hpA["k_LA"], 3),
          "Kt": "0.258", "nStartsH": str(len(tb["azimuth"]["hybrid"]["starts"]))}
    for ax, X in (("azimut", "A"), ("polare", "P")):
        g = gains[ax]["estimated"]
        M |= {f"Kpp{X}": f"{g['Kpp']:.1f}", f"Kpv{X}": f"{g['Kpv']:.0f}", f"Kiv{X}": f"{g['Kiv_hz']:.2f}",
              f"icmdMax{X}": f"{lim[ax]['i_cmd']:.1f}", f"vMax{X}": f"{lim[ax]['v_phase']:.1f}"}
        ak = [r for r in cl_akd if r["axis"] == ax]
        M |= {f"akdErr{X}": f"{max(float(r['max_following_error_model_mdeg']) for r in ak):.0f}",
              f"akdErrMeas{X}": f"{max(float(r['max_following_error_meas_mdeg']) for r in ak):.0f}",
              f"akdImax{X}": f"{max(float(r['max_current_model_A']) for r in ak):.0f}",
              f"akdImaxMeas{X}": f"{max(float(r['max_current_meas_A']) for r in ak):.1f}"}
        big = [r for r in cl if r["axis"] == ax and r["test"] not in ("01deg_minus01deg",)
               and not (ax == "azimut" and r["test"] == "8deg_minus8deg" and r["move"] == "#1")]
        M |= {f"clFitMin{X}": f"{min(float(r['fit_position_window_pct']) for r in big):.0f}",
              f"clFitMax{X}": f"{max(float(r['fit_position_window_pct']) for r in big):.0f}"}
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("figures, tables and numbers written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "dahl_bearings_report.tex"],
                       cwd=HERE, check=True)


if __name__ == "__main__":
    main()
