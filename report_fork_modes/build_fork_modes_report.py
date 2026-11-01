"""Builds the report on the modal least-squares identification of the fork modes (elevation and azimuth).

    ../.venv/bin/python -I build_fork_modes_report.py [--no-latex]

Inputs: the summary.json, tables and figures of analisi_forcella_modale.py in results_fork_modes_modal (the reported
model: zeta_k <= 0.05, CAD inertias and stiffnesses +-20 %, K_t +-10 %), results_fork_modes_modal_cad20 (fourth
step, K_t nominal), results_fork_modes_modal_zeta001 (modal fit with
zeta_k <= 0.01), results_fork_modes_modal_zeta005_J15 (zeta_k <= 0.05, CAD +-15 %),
results_fork_modes_modal_zeta005_J20 (intermediate run: inertias +-20 %, stiffnesses +-15 %),
results_fork_modes_modal_diag_cad50 (CAD +-50 %, diagnostic) and of analisi_modi_forcella.py in
results_fork_modes_2masses (first step, logarithmic error). No fit is run here. The text is in fork_modes_report.tex; this script copies the figures
(figures/), writes the tables and numbers (generated/) and compiles with latexmk. The report of the previous step is
kept in previous_logfit/.
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

from gimbal_mdof import forcella_modale as fm  # noqa: E402
from gimbal_mdof.sperimentale import carica_bode  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
RES = ROOT / "results_fork_modes_modal"
DIAG = {"cad": ROOT / "results_fork_modes_modal_diag_cad50"}
STEP2 = ROOT / "results_fork_modes_modal_zeta001"
STEP3 = ROOT / "results_fork_modes_modal_zeta005_J15"
INTER = ROOT / "results_fork_modes_modal_zeta005_J20"
STEP4 = ROOT / "results_fork_modes_modal_cad20"
PREV = ROOT / "results_fork_modes_2masses"
DATA_ROOT = ROOT.parent / "LucyIdentification"
COPY = ["polar_bode", "azimuth_bode", "polar_modal_terms", "azimuth_modal_terms", "fit_errors"]
TEX_NAME = {"s_J_e": r"$s_{J_e}$", "s_J_hp": r"$s_{J_{hp}}$", "s_J_m": r"$s_{J_m}$", "s_J_r": r"$s_{J_\mathrm{rear}}$",
            "s_J_L": r"$s_{J_L}$", "k_g": r"$k_g$", "eta_g": r"$\eta_g$", "tau_d": r"$\tau_d$", "alpha": r"$\alpha$",
            "k_LA": r"$k_{LA}$", "g_Kt": r"$g_{K_t}$"}
UNIT = {"N m/rad": r"N\,m/rad", "s": "s", "Hz": "Hz", "1/(kg m^2)": r"1/(kg\,m$^2$)", "": "--"}
CLASS = {"s_J": "CAD inertia", "s_k": "CAD stiffness", "g_Kt": "motor", "zeta": "modal", "k_g": "ground, free", "eta_g": "ground, free",
         "k_LA": "introduced, free", "alpha": "introduced, free", "f_": "stator mass, free",
         "a_": "stator mass, free", "tau_d": "drive"}


def tname(n: str, axis: str) -> str:
    if n.startswith("s_k_"):
        return rf"$s_{{k_{{\theta{n[-1]}}}}}$" if axis == "elevation" else rf"$s_{{k_{n[-1]}}}$"
    if n in TEX_NAME:
        return TEX_NAME[n]
    base, _, r = n.rpartition("_")
    return {"f": f"$f_{r}$", "zeta": rf"$\zeta_{r}$", "a": f"$a_{r}$"}[base]


def klass(n: str) -> str:
    return next(v for k, v in CLASS.items() if n.startswith(k))


def neg(t: str) -> str:
    return t.replace("-", "$-$", 1) if t.startswith("-") else t


def num(x: float, d: int = 4) -> str:
    if not np.isfinite(x):
        return "--"
    s = f"{x:.{d}g}"
    if "e" in s:
        m, e = s.split("e")
        return rf"\ensuremath{{{m}\times10^{{{int(e)}}}}}"
    return neg(s)


def esc(t: str) -> str:
    return t.replace("%", r"\%")


def std_txt(v: float | None) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "--"
    if not np.isfinite(v) or 100 * v > 1e3:
        return "n.d."
    return f"{100 * v:.1f}\\,\\%"


def rows_tex(rows: list[list[str]]) -> str:
    return " \\\\\n".join(" & ".join(r) for r in rows)


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def starts(path: Path, axis: str) -> dict[str, list[float]]:
    out: dict[str, set[float]] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["axis"] == axis:
                for k, v in json.loads(r["start"]).items():
                    out.setdefault(k, set()).add(float(f"{v:.4g}"))
    return {k: sorted(v) for k, v in out.items()}


def cad_value(n: str, axis: str, cad: dict) -> tuple[float, str] | None:
    """CAD value and unit of the quantity that a factor scales."""
    if axis == "elevation":
        J, roles = cad["elevation_J"], cad["elevation_roles"]
        jj = lambda r: J[roles.index(r)]  # noqa: E731
        k = {frozenset((roles[i], roles[j])): kk for i, j, kk in cad["elevation_springs"]}
        base = {"s_J_e": jj("e"), "s_J_hp": jj("h") + jj("p"), "s_J_m": jj("m"), "s_J_r": jj("r"),
                "s_k_1": k[frozenset(("e", "h"))], "s_k_2": k[frozenset(("h", "m"))],
                "s_k_3": k[frozenset(("m", "r"))]}
    else:
        a = cad["azimuth"]
        base = {"s_J_e": a["J_e"], "s_J_m": a["J_m"], "s_J_L": a["J_L"] + a["J_p"], "s_k_1": a["k"][0],
                "s_k_2": a["k"][1]}
    if n not in base:
        return None
    return base[n], (r"kg\,m$^2$" if "_J_" in n else r"N\,m/rad")


def param_rows(par: list, s: dict, st: dict, axis: str, cad: dict) -> list[list[str]]:
    rows = []
    for p in par:
        x = s["parameters"][p.nome]
        s0 = st.get(p.nome, [])
        if p.nome.startswith("a_"):
            start = ", ".join(rf"$\pm{v:g}$" for v in sorted({abs(y) for y in s0}))
        elif len(s0) > 3:
            start = f"{len(s0)} values"
        else:
            start = ", ".join(num(v, 3) for v in s0) or "--"
        cv = cad_value(p.nome, axis, cad)
        phys = "--" if cv is None else f"{num(x * cv[0], 3)} {cv[1]}"
        cadv = "--" if cv is None else num(cv[0], 3)
        rows.append([tname(p.nome, axis), klass(p.nome), UNIT[p.unita], num(x), cadv, phys, start,
                     f"[{num(p.lo, 3)}, {num(p.hi, 3)}]", std_txt(s["std_rel"][p.nome]),
                     "yes" if p.nome in s["at_bound"] else ""])
    return rows


def mode_rows(s: dict) -> list[list[str]]:
    return [[str(k), f"{m['f_Hz']:.1f}", f"{m['zeta']:.4f}", num(m["residue"], 3), num(m["eta_bearing"], 2)]
            + [neg(f"{m['shape'][n]:+.2f}") for n in s["nodes"]] for k, m in enumerate(s["modes"], 1)]


def nmp(el: dict) -> tuple[float, complex]:
    """B^2 - 4AC and one zero of the rigid axis with the two undamped stator modes."""
    J = sum(el["M_diag"][:4])
    (w1, a1), (w2, a2) = [(2 * np.pi * m["f_Hz"], m["a"]) for m in el["stator_masses"]]
    A, B, C = 1 / J + a1 + a2, (w1 ** 2 + w2 ** 2) / J + a1 * w2 ** 2 + a2 * w1 ** 2, w1 ** 2 * w2 ** 2 / J
    return B ** 2 - 4 * A * C, np.sqrt(np.roots([A, B, C]).astype(complex)[0])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for old in FIG.glob("*"):
        old.unlink()
    for name in COPY:
        shutil.copyfile(RES / "figures" / f"{name}.pdf", FIG / f"{name}.pdf")
    s = json.loads((RES / "summary.json").read_text())
    diag = {k: json.loads((v / "summary.json").read_text()) for k, v in DIAG.items()}
    prev = json.loads((PREV / "summary.json").read_text())
    s2 = json.loads((STEP2 / "summary.json").read_text())
    s3 = json.loads((STEP3 / "summary.json").read_text())
    si = json.loads((INTER / "summary.json").read_text())
    s4 = json.loads((STEP4 / "summary.json").read_text())
    CLASS["g_Kt"] = rf"motor, $\pm{100 * s['kt_variation']:.0f}$\,\%"
    CLASS["zeta"] = rf"modal, $\leq{s['zeta_max']:g}$"
    CLASS["s_J"] = rf"CAD, $\pm{100 * s['cad_inertia_variation']:.0f}$\,\%"
    CLASS["s_k"] = rf"CAD, $\pm{100 * s['cad_variation']:.0f}$\,\%"
    el, az = s["elevation"], s["azimuth"]
    assert s["fit_band_Hz"] == [1.0, 150.0], s["fit_band_Hz"]
    M: dict[str, str] = {}

    # ---- data
    rows = []
    for axis, fn in [("elevation", f) for f in el["measurements"]] + [("azimuth", az["measurement"])]:
        f, _, meta = carica_bode(DATA_ROOT / fn)
        b = (f >= 1.0) & (f <= 150.0)
        grid = 100 * np.median(np.diff(np.log(f)))
        rows.append([axis, rf"\path{{{Path(fn).name}}}", f"{f[0]:.0f}--{f[-1]:.0f}", f"{grid:.0f}",
                     f"{float(meta['MeasurementConfiguration/CurrentExcitation']):g}", str(int(b.sum()))])
    write("tab_data.tex", rows_tex(rows))

    # ---- parameters and modes
    st = {a: starts(RES / "tables" / "starts.csv", a) for a in ("elevation", "azimuth")}
    write("tab_polar_params.tex", rows_tex(param_rows(fm.parametri_polari(s["zeta_max"], s["cad_variation"],
                                                                           s["cad_inertia_variation"],
                                                                           s["kt_variation"]), el,
                                                     st["elevation"], "elevation", s["cad"])))
    write("tab_azimuth_params.tex", rows_tex(param_rows(fm.parametri_azimut(s["zeta_max"], s["cad_variation"],
                                                                             s["cad_inertia_variation"],
                                                                             s["kt_variation"]), az,
                                                       st["azimuth"], "azimuth", s["cad"])))
    write("tab_polar_modes.tex", rows_tex(mode_rows(el)))
    write("tab_azimuth_modes.tex", rows_tex(mode_rows(az)))
    vp, va = el["parameters"], az["parameters"]

    def best(c: list[float]) -> str:
        return str(sum(abs(x - min(c)) < 1e-4 * min(c) for x in c))

    M |= {"costP": num(el["cost"], 3), "costA": num(az["cost"], 3), "nStartsP": str(len(el["starts"])),
          "nStartsA": str(len(az["starts"])), "nStartsPbest": best(el["starts"]), "nStartsAbest": best(az["starts"]),
          "costsA": ", ".join(num(c, 3) for c in sorted({float(f"{c:.3g}") for c in az["starts"]})),
          "nBoundP": str(len(el["at_bound"])), "nParP": str(len(vp)), "nBoundA": str(len(az["at_bound"])),
          "nParA": str(len(va)), "kg": f"{vp['k_g']:.0f}", "etag": f"{vp['eta_g']:.2f}",
          "kgA": f"{va['k_g']:.0f}", "etaA": f"{va['eta_g']:.2f}", "taudA": f"{1e3 * va['tau_d']:.2f}",
          "alphaA": f"{va['alpha']:.3f}", "kLA": num(va["k_LA"], 3), "Kt": f"{s['Kt_Nm_A']:g}",
          "fCL": f"{s['current_loop_Hz'] / 1e3:g}", "JL": f"{az['M_diag'][2]:.3f}", "JAP": f"{az['M_diag'][3]:.3f}",
          "kgPrev": f"{prev['elevation']['parameters']['k_g']:.0f}",
          "etagPrev": f"{prev['elevation']['parameters']['eta_g']:.2f}"}
    for r, m in enumerate(el["stator_masses"], 1):
        R = "I" * r
        M |= {f"fMass{R}": f"{m['f_Hz']:.1f}", f"aMass{R}": neg(f"{m['a']:+.3f}"), f"JMass{R}": f"{m['J_kgm2']:.2f}",
              f"zMass{R}": f"{el['modes'][r]['zeta']:.3f}",
              f"zStd{R}": std_txt(el["std_rel"][f"zeta_{r + 1}"]),
              f"aLog{R}": neg(f"{prev['elevation']['stator_masses'][r - 1]['a']:+.3f}"),
              f"zLog{R}": f"{prev['elevation']['stator_masses'][r - 1]['zeta']:.3f}",
              f"aStep{R}": neg(f"{s2['elevation']['stator_masses'][r - 1]['a']:+.3f}"),
              f"fStep{R}": f"{s2['elevation']['stator_masses'][r - 1]['f_Hz']:.1f}"}
    md, ma = el["modes"], az["modes"]
    M |= {"fBearP": f"{md[0]['f_Hz']:.2f}", "etaBearP": f"{md[0]['eta_bearing']:.2f}",
          "fShaftP": f"{md[3]['f_Hz']:.0f}", "fShaftPb": f"{md[4]['f_Hz']:.0f}", "fBearA": f"{ma[0]['f_Hz']:.2f}",
          "etaBearA": f"{ma[0]['eta_bearing']:.2f}", "fForkA": f"{ma[1]['f_Hz']:.1f}",
          "etaForkA": f"{ma[1]['eta_bearing']:.4f}", "fShaftA": f"{ma[2]['f_Hz']:.0f}",
          "zForkA": f"{ma[1]['zeta']:.3f}", "zForkStd": std_txt(az["std_rel"]["zeta_2"]),
          "zForkStep": f"{s2['azimuth']['modes'][1]['zeta']:.3f}", "fForkStep": f"{s2['azimuth']['modes'][1]['f_Hz']:.1f}",
          "kgStep": f"{s2['elevation']['parameters']['k_g']:.0f}",
          "etagStep": f"{s2['elevation']['parameters']['eta_g']:.2f}",
          "kgAStep": f"{s2['azimuth']['parameters']['k_g']:.0f}", "etaAStep": f"{s2['azimuth']['parameters']['eta_g']:.2f}",
          "alphaStep": f"{s2['azimuth']['parameters']['alpha']:.3f}", "kLAStep": num(s2["azimuth"]["parameters"]["k_LA"], 3),
          "costPStep": num(s2["elevation"]["cost"], 3), "costAStep": num(s2["azimuth"]["cost"], 3),
          "zetaMax": f"{s['zeta_max']:g}", "zetaMaxStep": f"{s2['zeta_max']:g}",
          "cadJ": f"{100 * s['cad_inertia_variation']:.0f}", "cadK": f"{100 * s['cad_variation']:.0f}",
          "sJLA": f"{va['s_J_L']:.2f}", "sJhpP": f"{vp['s_J_hp']:.3f}",
          "sJLAStepJ": f"{s3['azimuth']['parameters']['s_J_L']:.2f}",
          "costPStepJ": num(s3["elevation"]["cost"], 3), "costAStepJ": num(s3["azimuth"]["cost"], 3),
          "costRedA": f"{100 * (1 - s4['azimuth']['cost'] / s3['azimuth']['cost']):.0f}",
          "costRedKt": f"{100 * (1 - az['cost'] / s4['azimuth']['cost']):.0f}",
          "costRedKtP": f"{100 * (1 - el['cost'] / s4['elevation']['cost']):.3f}",
          "gKtA": f"{va['g_Kt']:.2f}", "gKtP": f"{vp['g_Kt']:.3f}", "KtA": f"{s['Kt_Nm_A'] * va['g_Kt']:.3f}",
          "ktVar": f"{100 * s['kt_variation']:.0f}", "ktLo": f"{1 - s['kt_variation']:.2f}",
          "ktHi": f"{1 + s['kt_variation']:.2f}",
          "sktwoAFour": f"{s4['azimuth']['parameters']['s_k_2']:.2f}",
          "sJLAFour": f"{s4['azimuth']['parameters']['s_J_L']:.2f}",
          "JtotAFour": f"{s4['azimuth']['M_diag'][2] + s4['azimuth']['M_diag'][3]:.3f}",
          "costRedAInter": f"{100 * (1 - si['azimuth']['cost'] / s3['azimuth']['cost']):.1f}",
          "fitAInter": f"{si['fit_quality'][3]['fit_all']:.0f}", "sktwoA": f"{va['s_k_2']:.2f}",
          "costRedP": f"{100 * (1 - s4['elevation']['cost'] / s3['elevation']['cost']):.2f}",
          "JtotA": f"{az['M_diag'][2] + az['M_diag'][3]:.3f}",
          "JtotAStepJ": f"{s3['azimuth']['M_diag'][2] + s3['azimuth']['M_diag'][3]:.3f}",
          "cadCostP": f"{100 * abs(1 - diag['cad']['elevation']['cost'] / el['cost']):.1f}",
          "cadCostA": f"{az['cost'] / diag['cad']['azimuth']['cost']:.1f}"}
    # sign of the residues: rigid axis and two undamped stator modes
    for key, src in (("", el), ("Step", s2["elevation"])):
        disc, z = nmp(src)
        M |= {f"nmpDisc{key}": num(disc, 2), f"nmpF{key}": f"{abs(z) / (2 * np.pi):.1f}",
              f"nmpZeta{key}": f"{abs(z.real) / abs(z):.3f}"}
    M["Jaxis"] = f"{sum(el['M_diag'][:4]):.3f}"

    # ---- fit quality
    def g(q: dict, k: str) -> str:
        return "--" if q[k] is None or not np.isfinite(q[k]) else neg(f"{q[k]:.1f}")

    rows = [["elevation" if q["measurement"].startswith("elevation") else "azimuth", esc(q["measurement"]),
             g(q, "fit_all"), g(q, "fit_1-5 Hz"), g(q, "fit_5-100 Hz"), g(q, "fit_>100 Hz"), f"{q['rms_dB']:.2f}",
             f"{q['rms_deg']:.1f}"] for q in s["fit_quality"]]
    write("tab_fit.tex", rows_tex(rows))
    for key, src in (("", s), ("Step", s2), ("StepJ", s3), ("Four", s4), ("Cad", diag["cad"])):
        for k, q in enumerate(src["fit_quality"]):
            R = ["PI", "PII", "PIII", "A"][k]
            M |= {f"fit{R}{key}": f"{q['fit_all']:.0f}", f"db{R}{key}": f"{q['rms_dB']:.1f}",
                  f"deg{R}{key}": f"{q['rms_deg']:.1f}", f"fitLow{R}{key}": g(q, "fit_1-5 Hz")}
    prevq = {q["measurement"]: q for q in prev["fit_quality"] if q["model"] in ("H7 + 2 stator masses", "flexible load")}
    dq = {k: {q["measurement"]: q for q in d["fit_quality"]} for k, d in diag.items()}
    rows = []
    for q in s["fit_quality"]:
        m = q["measurement"]
        row = [esc(m)]
        for src in (prevq[m], {q2["measurement"]: q2 for q2 in s2["fit_quality"]}[m],
                    {q3["measurement"]: q3 for q3 in s3["fit_quality"]}[m],
                    {q4["measurement"]: q4 for q4 in s4["fit_quality"]}[m], q, dq["cad"][m]):
            row += [g(src, "fit_all"), f"{src['rms_dB']:.1f}"]
        rows.append(row)
    write("tab_compare.tex", rows_tex(rows))
    M["boundP"] = ", ".join(tname(n, "elevation") for n in el["at_bound"])
    M["boundA"] = ", ".join(tname(n, "azimuth") for n in az["at_bound"])

    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("figures, tables and numbers written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "fork_modes_report.tex"],
                       cwd=HERE, check=True)


if __name__ == "__main__":
    main()
