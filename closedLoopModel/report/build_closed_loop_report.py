"""Builds the report on the comparison of the closed-loop model of LUCY with the measurements.

    ../../.venv/bin/python -I build_closed_loop_report.py [--no-latex]

Inputs: the results of ``main.py``, ``validate.py`` and ``compare.py`` (summary.json files listed in RESULTS) and
the measurements in sim/LucyIdentification. The time histories of the figures are simulated again here with the
library (a few windows), with the parameters stored in those results. The text is in closed_loop_report.tex;
this script writes the tables and numbers (generated/) and the figures (figures/), then compiles with latexmk.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
CL = HERE.parent
sys.path.insert(0, str(CL))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import main as mn  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import identification as idn  # noqa: E402
from lucy_cl.closed_loop import loop_frf  # noqa: E402
from lucy_cl.data_files import read_akd, read_workbench  # noqa: E402
from lucy_cl.driver import controller_workbench  # noqa: E402
from lucy_cl.validation import read_sequence, rigid_body_regression  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
LI = CL.parent.parent / "LucyIdentification"
SEP_DIR, OCT_DIR = LI / "Validation", LI / "VClosedLoop"
RESULTS = {"sep": CL / "results_validation_H6_zeta10" / "summary.json",
           "oct_time": CL / "results_VClosedLoop_H6_zeta10" / "time_comparison" / "summary.json",
           "oct_main": CL / "results_VClosedLoop_H6_zeta10" / "summary.json"}
CFG_SEP, CFG_OCT = CL / "config_zeta10.yaml", CL / "config_vcl.yaml"
MDEG = 1e3 * 180.0 / np.pi
AMPS = 2.0

plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 6.5,
                     "savefig.bbox": "tight", "lines.linewidth": 1.1})
STYLE = {"measured": dict(color="k", lw=1.6), "akd": dict(color="C3", ls="--"), "estimate": dict(color="C0"),
         "A": dict(color="C3", ls="--"), "B": dict(color="C1", ls="-."), "F": dict(color="C2"),
         "H": dict(color="C0", ls=(0, (1, 1)), lw=1.6)}


def tt(s: str) -> str:
    return r"\texttt{" + s.replace("_", r"\_") + "}"


def save(fig, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=110)
    plt.close(fig)


def rows_tex(rows: list[list[str]]) -> str:
    return " \\\\\n".join(" & ".join(r) for r in rows)


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def f0(x: float, nd: int = 0) -> str:
    if x is None or not np.isfinite(x):
        return "--"
    return f"{x:.{nd}f}"


# --------------------------------------------------------------------------- simulation helpers


def loop_for(cfg_path: Path, seq_dir: Path):
    cfg = yaml.safe_load(cfg_path.read_text())
    idn._init_worker(cfg, str(CL), str(seq_dir), AMPS)
    return cfg


def sim(window_seq: str, axis: str, t_move: float, values: dict, seq_dir: Path, pre=0.3, post=2.5) -> dict:
    seq = read_sequence(seq_dir / f"{window_seq}.csv", AMPS)
    w = [x for x in idn.windows_for_moves(seq, axis, pre=pre, post=post) if abs(x.move - t_move) < 0.05][0]
    return idn.simulate_job((w, values, 200.0))


def window_figure(cases: list[tuple[str, dict[str, dict]]], curves: list[tuple[str, str]], axis_of: list[str],
                  name: str, held: bool = False) -> None:
    """Rows: windows (title, {curve: result}); columns: position error and current of the moving axis
    (and error of the held axis). Every panel has its legend and units."""
    ncol = 3 if held else 2
    fig, axs = plt.subplots(len(cases), ncol, figsize=(6.6 if not held else 9.6, 2.35 * len(cases) + 0.3),
                            squeeze=False)
    for k, ((title, res), ax_m) in enumerate(zip(cases, axis_of)):
        first = res[curves[0][0]]
        t = first["t"]
        cols = [(ax_m, "error", MDEG, "command $-$ encoder [mdeg]"), (ax_m, "current", 1.0, "current [A rms]")]
        if held:
            other = [a for a in AXES if a != ax_m][0]
            cols.append((other, "error", MDEG, "held axis: command $-$ encoder [mdeg]"))
        for j, (axis, key, sc, ylab) in enumerate(cols):
            a = axs[k, j]
            a.plot(t, sc * first[axis][f"{key}_meas"], label="measured", **STYLE["measured"])
            for cname, clabel in curves:
                a.plot(t, sc * res[cname][axis][key], label=clabel, **STYLE[cname])
            a.set_ylabel(ylab)
            a.set_xlabel("time from window start [s]")
            a.set_title(f"{title}: {axis}" + (" (held)" if axis != ax_m else ""), fontsize=7.5)
            a.legend(loc="best")
    fig.tight_layout()
    save(fig, name)


# --------------------------------------------------------------------------- figures and tables


def drive_tables(cfg_sep: dict, cfg_oct: dict) -> dict:
    out = {}
    rows = []
    keys = [("PL.KP", "PL.KP"), ("PL.KI", "PL.KI"), ("VL.KP", "VL.KP [A/(rad/s)]"), ("VL.KI", "VL.KI [Hz]"),
            ("VL.KVFF", "VL.KVFF"), ("OBS.ENABLE", "observer (OBS.ENABLE)"), ("IL.LIMITP", "IL.LIMITP [A rms]")]
    files = {}
    for camp, cfg in (("sep", cfg_sep), ("oct", cfg_oct)):
        for ax in AXES:
            files[(camp, ax)] = read_akd((CL / cfg["drives"][ax]["parameters"]).resolve())
    for k, lab in keys:
        r = [lab]
        for ax in AXES:
            for camp in ("sep", "oct"):
                r.append(f"{float(files[(camp, ax)][f'AXIS1.{k}']):g}")
        rows.append(r)

    def filt(a: dict) -> str:
        out_ = []
        for k in (1, 2, 3, 4):
            if int(float(a[f"AXIS1.VL.ARTYPE{k}"])) == 0:
                continue
            zf, zq, pf, pq = (float(a[f"AXIS1.VL.AR{x}{k}"]) for x in ("ZF", "ZQ", "PF", "PQ"))
            if zf == pf and zq == pq:
                continue
            if zf == pf:
                out_.append(f"AR{k}: {zf:g}\\,Hz, {20 * np.log10(pq / zq):+.0f}\\,dB")
            else:
                out_.append(f"AR{k}: {zf:g}/{pf:g}\\,Hz")
        return "; ".join(out_) or "none"

    r = ["active biquads"]
    for ax in AXES:
        for camp in ("sep", "oct"):
            r.append(filt(files[(camp, ax)]))
    rows.append(r)
    write("tab_drive_files.tex", rows_tex(rows))
    out["files"] = files
    return out


def current_scale(seqs_sep: dict, J: dict) -> None:
    rows = []
    for n, s in seqs_sep.items():
        for ax in AXES:
            rr = rigid_body_regression(s, ax, J[ax])
            rows.append([tt(n), ax, f"{rr['slope']:.2f}", f"{rr['torque_per_unit']:.3f}",
                         f"{rr['torque_per_unit'] * AMPS:.3f}", f"{rr['coulomb_Nm']:.2f}", f"{rr['R2']:.3f}"])
    write("tab_current_scale.tex", rows_tex(rows))
    # figure: current against acceleration (identification sequence)
    from scipy.signal import butter, filtfilt

    s = seqs_sep["MPCLucy_pos_dist_dataset"]
    dt = float(np.median(np.diff(s.t)))
    b, a = butter(2, 3.0 * 2 * dt)
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.2))
    for j, ax in enumerate(AXES):
        th = filtfilt(b, a, s.enc[ax])
        i = filtfilt(b, a, s.current[ax] * AMPS)
        w = np.gradient(th, dt)
        al = np.gradient(w, dt)
        m = np.abs(w) > np.deg2rad(1.0)
        rr = rigid_body_regression(s, ax, J[ax])
        for sg, c in ((1, "C0"), (-1, "C3")):
            mm = m & (np.sign(w) == sg)
            axs[j].plot(al[mm][::20], i[mm][::20], ".", ms=1.5, color=c, alpha=0.5,
                        label=f"measured, {'positive' if sg > 0 else 'negative'} speed")
            xa = np.linspace(al[mm].min(), al[mm].max(), 10)
            axs[j].plot(xa, rr["slope"] * xa + sg * rr["coulomb_units"] + rr["offset_units"], color=c, lw=1.4,
                        label=f"fit: {rr['slope']:.2f} a {'+' if sg > 0 else '-'} {rr['coulomb_units']:.2f}")
        axs[j].set_xlabel(r"encoder acceleration (3 Hz low pass) [rad/s$^2$]")
        axs[j].set_ylabel("AXIS_AMPS channel [units]")
        axs[j].set_title(f"{ax}: J/a = {rr['torque_per_unit']:.3f} N m per unit (J = {J[ax]:.2f} kg m$^2$)",
                         fontsize=7.5)
        axs[j].legend(loc="best", markerscale=5)
    fig.tight_layout()
    save(fig, "fig_current_scale")


def september(sep: dict, seqs_sep: dict) -> dict:
    cfg = loop_for(CFG_SEP, SEP_DIR)
    cl = idn._W["cl"]
    nominal = {ax: idn.current_values(cl, ax) for ax in AXES}
    nominal = {ax: {k: v for k, v in d.items() if k in ("Kpp", "Kpv", "Kiv", "Tc")} for ax, d in nominal.items()}
    est = {ax: sep["estimate"][ax]["values"] for ax in AXES}
    # estimate table
    rows = []
    unit = {"Kpp": "1/s", "Kpv": "A/(rad/s)", "Kiv": "Hz", "Tc": "N m"}
    name = {"Kpp": r"$K_{pp}$", "Kpv": r"$K_{pv}$", "Kiv": r"$K_{iv}$", "Tc": r"$T_c$"}
    for ax in AXES:
        for k in ("Kpp", "Kpv", "Kiv", "Tc"):
            st = sep["estimate"][ax].get("start", {}).get(k, float("nan"))
            rows.append([ax, f"{name[k]} [{unit[k]}]", f"{nominal[ax][k]:.4g}", f"{st:.4g}", f"{est[ax][k]:.4g}",
                         f"{est[ax][k] / nominal[ax][k]:.3g}"])
    write("tab_estimate.tex", rows_tex(rows))
    # validation summary
    rows = []
    fits = {}
    for key, rec in sep["validation"].items():
        n, ax = key.split("/")

        def st(k: str) -> str:
            x = np.array([r[k] for r in rec], dtype=float)
            x = x[np.isfinite(x)]
            return "--" if x.size == 0 else f"{np.median(x):.0f} [{x.min():.0f}, {x.max():.0f}]"

        fits[key] = {k: float(np.nanmedian([r[k] for r in rec])) for k in
                     ("fit_error_akd", "fit_error_estimate", "fit_current_akd", "fit_current_estimate")}
        rows.append([tt(n), ax, str(len(rec)), st("fit_error_akd"), st("fit_error_estimate"), st("fit_current_akd"),
                     st("fit_current_estimate"), st("held_max_meas_mdeg"), st("held_max_estimate_mdeg")])
    write("tab_sep_fit.tex", rows_tex(rows))
    # bar chart
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.0))
    keys = list(fits)
    x = np.arange(len(keys))
    for j, (ch, lab) in enumerate((("error", "position error"), ("current", "current"))):
        for o, (var, c, l) in enumerate((("akd", "C3", "model, .akd parameters (2$\\pi$ PL.KP)"),
                                         ("estimate", "C0", "model, estimated parameters"))):
            v = np.array([fits[k][f"fit_{ch}_{var}"] for k in keys])
            axs[j].bar(x + (o - 0.5) * 0.38, np.clip(v, -100, 100), 0.36, color=c, label=l)
            for xi, vi in zip(x, v):
                if vi < -100:
                    axs[j].text(xi + (o - 0.5) * 0.38, -97, f"{vi:.0f}", ha="center", fontsize=6, rotation=90)
        axs[j].set_xticks(x, [k.replace("MPCLucy_", "").replace("pos_dist_", "").replace("/", "\n") for k in keys],
                          fontsize=7)
        axs[j].set_ylabel(f"median NRMSE fit of the {lab} [%]")
        axs[j].set_ylim(-100, 100)
        axs[j].axhline(0, color="k", lw=0.6)
        axs[j].legend(loc="lower right")
        axs[j].set_title(f"September validation sequences: {lab}", fontsize=7.5)
    fig.tight_layout()
    save(fig, "fig_sep_fit")
    # windows
    cases, axis_of = [], []
    for seqn, ax, tm, title in (("MPCLucy_scotched", "azimuth", 55.68, "scotched, +10 deg at 55.68 s"),
                                ("MPCLucy_pos_dist_test", "elevation", 51.87, "test, -10 deg at 51.87 s"),
                                ("MPCLucy_pos_dist_test", "elevation", 39.42, "test, +1 deg at 39.42 s")):
        res = {"akd": sim(seqn, ax, tm, nominal, SEP_DIR), "estimate": sim(seqn, ax, tm, est, SEP_DIR)}
        cases.append((title, res))
        axis_of.append(ax)
    window_figure(cases, [("akd", "model, .akd parameters (2$\\pi$ PL.KP, no FF)"),
                          ("estimate", "model, estimated parameters")], axis_of, "fig_sep_moves")
    # held axis (September): azimuth move, elevation held
    res = {"akd": sim("MPCLucy_scotched", "azimuth", 51.65, nominal, SEP_DIR),
           "estimate": sim("MPCLucy_scotched", "azimuth", 51.65, est, SEP_DIR)}
    return {"nominal": nominal, "est": est, "fits": fits, "held_sep": res, "cfg": cfg}


def october(oct_t: dict, oct_m: dict) -> dict:
    cfg = loop_for(CFG_OCT, OCT_DIR)
    cl = idn._W["cl"]
    vals = {k: v["values"] for k, v in oct_t["variants"].items()}
    labels = {k: v["label"] for k, v in oct_t["variants"].items()}
    # variants table
    why = {"A": "files as WorkBench reads them (Controller trace)",
           "B": "PL.KP in 1/s, as suggested by the September sequences",
           "C": "A with the velocity feedforward of the files",
           "D": "B with the velocity feedforward of the files",
           "E": "A with the Coulomb level of the September estimate",
           "F": "B with the Coulomb level of the September estimate",
           "G": "position gains of the September estimate, velocity loop of the files",
           "H": "all four parameters of the September estimate (cross-check)"}
    rows = []
    for k, v in vals.items():
        rows.append([k, why[k]] + [f"{v[ax]['Kpp']:.3g}" for ax in AXES]
                    + [f"{v[ax].get('Kpv', cl.drives[ax].p.Kpv):.3g} / {v[ax].get('Kiv', cl.drives[ax].p.Kiv_hz):.3g}"
                       for ax in AXES]
                    + [f"{v[ax]['kvff']:.0f}" for ax in AXES[:1]] + [f"{v['azimuth']['Tc']:.2f} / {v['elevation']['Tc']:.2f}"])
    write("tab_variants.tex", rows_tex(rows))
    # fit table
    moves = oct_t["moves"]
    rows, med = [], {}
    for k in vals:
        r = [k]
        for ax in AXES:
            sel = [m[k] for m in moves if m["axis"] == ax]
            for ch in ("fit_error", "fit_current"):
                x = np.array([s[ch] for s in sel], dtype=float)
                med[(k, ax, ch)] = (float(np.median(x)), float(x.min()), float(x.max()))
                r.append(f"{np.median(x):.0f} [{x.min():.0f}, {x.max():.0f}]")
        rows.append(r)
    write("tab_oct_fit.tex", rows_tex(rows))
    pe = []
    for ax in AXES:
        sel = [m for m in moves if m["axis"] == ax]
        meas = np.array([m["max_error_meas_mdeg"] for m in sel])
        r = [ax, f"{np.median(meas):.0f}"]
        for k in vals:
            x = np.array([m[k]["max_error_mdeg"] for m in sel])
            r.append(f"{np.median(x):.0f}")
        pe.append(r)
    write("tab_oct_peak.tex", rows_tex(pe))
    # bar chart of the fits
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.1))
    ks = list(vals)
    x = np.arange(len(ks))
    for j, ax in enumerate(AXES):
        for o, (ch, c, lab) in enumerate((("fit_error", "C0", "position error"), ("fit_current", "C1", "current"))):
            v = np.array([med[(k, ax, ch)][0] for k in ks])
            lo = np.array([med[(k, ax, ch)][1] for k in ks])
            hi = np.array([med[(k, ax, ch)][2] for k in ks])
            vc = np.clip(v, -100, 100)
            axs[j].bar(x + (o - 0.5) * 0.38, vc, 0.36, color=c, label=f"{lab}: median, min-max")
            axs[j].errorbar(x + (o - 0.5) * 0.38, vc, yerr=[np.clip(vc - np.clip(lo, -100, 100), 0, None),
                                                            np.clip(np.clip(hi, -100, 100) - vc, 0, None)],
                            fmt="none", ecolor="k", lw=0.7, capsize=2)
            for xi, vi in zip(x, v):
                if vi < -100:
                    axs[j].text(xi + (o - 0.5) * 0.38, -95, f"{vi:.0f}", ha="center", fontsize=5.5, rotation=90)
        axs[j].set_xticks(x, ks)
        axs[j].set_xlabel("variant (see the table of the variants)")
        axs[j].set_ylabel("NRMSE fit [%]")
        axs[j].set_ylim(-100, 100)
        axs[j].axhline(0, color="k", lw=0.6)
        axs[j].set_title(f"October step tests, {ax} moves ({sum(m['axis'] == ax for m in moves)})", fontsize=7.5)
        axs[j].legend(loc="lower left")
    fig.tight_layout()
    save(fig, "fig_oct_fit")
    # windows
    curves = [("A", "A: files, 2$\\pi$ PL.KP"), ("B", "B: files, PL.KP in 1/s"),
              ("F", "F: B + Tc of September"), ("H", "H: September estimate")]
    sel = {k: vals[k] for k in ("A", "B", "F", "H")}
    for ax, wins, name in (("azimuth", (("MPCLucy_step_1deg", 5.76, "step 1 deg, +1 deg at 5.76 s"),
                                        ("MPCLucy_2degPol_zeroPol_2degAz_zeroAz", 19.95, "2 deg sequence, +2 deg at 19.95 s")),
                            "fig_oct_az"),
                           ("elevation", (("MPCLucy_step_1deg", 20.16, "step 1 deg, +1 deg at 20.16 s"),
                                          ("MPCLucy_2degAz_2degPol_zeroAz_zeroPol", 23.56, "2 deg sequence, +2 deg at 23.56 s")),
                            "fig_oct_el")):
        cases = []
        for seqn, tm, title in wins:
            cases.append((title, {k: sim(seqn, ax, tm, v, OCT_DIR) for k, v in sel.items()}))
        window_figure(cases, curves, [ax] * len(cases), name)
    held_oct = {k: sim("MPCLucy_2degPol_zeroPol_2degAz_zeroAz", "azimuth", 19.95, vals[k], OCT_DIR) for k in ("B", "H")}
    # frequency domain
    wb = {ax: read_workbench((CL / cfg["drives"][ax]["measurement"]).resolve()) for ax in AXES}
    old = {"azimuth": read_workbench(LI / "Azimuth_BodeMeasurement_Po_0.csv"),
           "elevation": read_workbench(LI / "VAntiRes" / "BodeMeasurement_12viti_dueDecimiOmogeneo.csv")}
    f = np.logspace(0, np.log10(2000), 1500)
    fig, axs = plt.subplots(2, 2, figsize=(9.0, 5.2), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    for j, ax in enumerate(AXES):
        d = cl.drives[ax]
        s = 2j * np.pi * wb[ax].f
        Cm, Cw = controller_workbench(d, 2j * np.pi * f), wb[ax].response("Controller")
        axs[0, j].semilogx(wb[ax].f, 20 * np.log10(np.abs(Cw)), "k.", ms=3, label="WorkBench Controller trace")
        axs[0, j].semilogx(f, 20 * np.log10(np.abs(Cm)), "C0", label="model of the drive")
        axs[1, j].semilogx(wb[ax].f, np.angle(Cw, deg=True), "k.", ms=3, label="WorkBench Controller trace")
        axs[1, j].semilogx(f, np.angle(Cm, deg=True), "C0", label="model of the drive")
        err = controller_workbench(d, s) / Cw
        axs[0, j].set_title(f"{ax}: max error {np.abs(20 * np.log10(np.abs(err))).max():.3f} dB, "
                            f"{np.abs(np.angle(err, deg=True)).max():.2f} deg", fontsize=7.5)
        axs[0, j].set_ylabel("|C| [dB re 1 A/(rad/s)]")
        axs[1, j].set_ylabel("phase [deg]")
        axs[1, j].set_xlabel("frequency [Hz]")
        for a in axs[:, j]:
            a.legend(loc="best")
    fig.tight_layout()
    save(fig, "fig_controller")
    fig, axs = plt.subplots(2, 2, figsize=(9.0, 5.4), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    for j, ax in enumerate(AXES):
        r = loop_frf(cl, ax, f)
        for w, lab, sty in ((old[ax], f"Sept. measurement ({old[ax].path.name.replace('_', ' ')})", dict(color="0.55", marker=".", ls="none", ms=3)),
                            (wb[ax], f"Oct. measurement ({wb[ax].path.name.replace('_', ' ')})", dict(color="k", marker=".", ls="none", ms=3))):
            G = w.response("Plant")
            axs[0, j].semilogx(w.f, 20 * np.log10(np.abs(G)), label=lab, **sty)
            axs[1, j].semilogx(w.f, np.angle(G, deg=True), label=lab, **sty)
        axs[0, j].semilogx(f, 20 * np.log10(np.abs(r.plant)), "C0", label="model (H6, zeta1 = 0.10)")
        axs[1, j].semilogx(f, np.angle(r.plant, deg=True), "C0", label="model (H6, zeta1 = 0.10)")
        axs[0, j].set_title(f"{ax}: plant v_fb / i_cmd", fontsize=7.5)
        axs[0, j].set_ylabel("|P| [dB re 1 (rad/s)/A]")
        axs[1, j].set_ylabel("phase [deg]")
        axs[1, j].set_xlabel("frequency [Hz]")
        for a in axs[:, j]:
            a.legend(loc="best", fontsize=5.8)
    fig.tight_layout()
    save(fig, "fig_bode_plant")
    # Bode fit table and margins
    rows = []
    for ax in AXES:
        m = oct_m["measurement"][ax]
        for tr in ("Plant", "Open Loop", "Closed Loop"):
            rows.append([ax, tr] + [f0(m[tr][b]) for b in ("all", "1-5 Hz", "5-100 Hz", ">100 Hz")])
    write("tab_bode_fit.tex", rows_tex(rows))
    rows = []
    for ax in AXES:
        m = oct_m["measurement"][ax]
        for lab, mm in (("measured", m["margins_measured"]), ("model", m["margins_model"])):
            rows.append([ax, lab, f"{mm['f_c']:.2f}", f"{mm['pm_deg']:.1f} at {mm['f_pm']:.2f}",
                         f"{mm['gm_db']:.1f} at {mm['f_gm']:.0f}"])
    write("tab_margins.tex", rows_tex(rows))
    # damping sensitivity
    rows = []
    for e in oct_m["zeta_sweep"]:
        r = [f"{e['zeta']:g}"]
        for k in ("azimuth_akd", "elevation_akd"):
            mm = e[k]
            r.append(f"{mm['pm_deg']:.0f} / {mm['gm_db']:.1f} at {mm['f_gm']:.0f}")
        r.append(f"{e['eig_max_real']:.3g}")
        rows.append(r)
    write("tab_zeta.tex", rows_tex(rows))
    return {"vals": vals, "labels": labels, "med": med, "held_oct": held_oct, "wb": wb, "old": old}


def held_figure(sep_r: dict, oct_r: dict) -> None:
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.0))
    for a, res, title, curves in ((axs[0], sep_r["held_sep"], "September, scotched: azimuth -2 deg at 51.65 s",
                                   (("akd", "model, .akd parameters"), ("estimate", "model, estimated parameters"))),
                                  (axs[1], oct_r["held_oct"], "October, 2 deg sequence: azimuth +2 deg at 19.95 s",
                                   (("B", "model, B: files, PL.KP in 1/s"), ("H", "model, H: September estimate")))):
        first = res[curves[0][0]]
        e0 = first["elevation"]["error_meas"][0]
        a.plot(first["t"], MDEG * (first["elevation"]["error_meas"] - e0), label="measured", **STYLE["measured"])
        for k, lab in curves:
            r = res[k]
            a.plot(r["t"], MDEG * (r["elevation"]["error"] - r["elevation"]["error"][0]), label=lab, **STYLE[k])
        a.set_xlabel("time from window start [s]")
        a.set_ylabel("elevation (held): change of command $-$ encoder [mdeg]")
        a.set_title(title, fontsize=7.5)
        a.legend(loc="best")
    fig.tight_layout()
    save(fig, "fig_held")


def zeta_figure() -> None:
    cfg = yaml.safe_load(CFG_OCT.read_text())
    cl, _, _ = mn.build(cfg, CL)
    pl = cl.plant
    f = np.logspace(1, np.log10(2000), 2500)
    fig, axs = plt.subplots(2, 1, figsize=(6.4, 4.4), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    z0 = pl.zeta
    for k, z in enumerate((0.01, 0.05, 0.10)):
        pl.set_zeta(z)
        L = loop_frf(cl, "elevation", f).loop
        lab = f"zeta1 = {z:g}" + (" (H6 as identified)" if z == 0.01 else "")
        axs[0].semilogx(f, 20 * np.log10(np.abs(L)), color=f"C{k}", label=lab)
        axs[1].semilogx(f, np.rad2deg(np.unwrap(np.angle(L))), color=f"C{k}", label=lab)
    pl.set_zeta(z0)
    axs[0].axhline(0, color="k", lw=0.6, label="0 dB")
    axs[1].axhline(-540, color="k", lw=0.6, ls=":", label="-540 deg")
    axs[1].axhline(-180, color="k", lw=0.6, ls="--", label="-180 deg")
    axs[0].set_ylabel("|L| [dB]")
    axs[1].set_ylabel("phase (unwrapped) [deg]")
    axs[1].set_xlabel("frequency [Hz]")
    axs[0].set_title("elevation: loop at the current command, October drive parameters", fontsize=7.5)
    for a in axs:
        a.legend(loc="best")
    fig.tight_layout()
    save(fig, "fig_zeta")


def numbers(sep_r: dict, oct_r: dict, oct_m: dict, seqs_sep: dict, seqs_oct: dict) -> None:
    med = oct_r["med"]
    est = sep_r["est"]
    macro = {
        "KppAzEst": f"{est['azimuth']['Kpp']:.2f}", "KppElEst": f"{est['elevation']['Kpp']:.2f}",
        "KpvAzEst": f"{est['azimuth']['Kpv']:.0f}", "KpvElEst": f"{est['elevation']['Kpv']:.0f}",
        "KivAzEst": f"{est['azimuth']['Kiv']:.2f}", "KivElEst": f"{est['elevation']['Kiv']:.2f}",
        "TcAzEst": f"{est['azimuth']['Tc']:.2f}", "TcElEst": f"{est['elevation']['Tc']:.2f}",
        "TcAzHsix": f"{sep_r['nominal']['azimuth']['Tc']:.2f}", "TcElHsix": f"{sep_r['nominal']['elevation']['Tc']:.2f}",
        "SepFitAzEst": f"{np.median([sep_r['fits'][k]['fit_error_estimate'] for k in sep_r['fits'] if 'azimuth' in k]):.0f}",
        "SepFitElEst": f"{np.median([sep_r['fits'][k]['fit_error_estimate'] for k in sep_r['fits'] if 'elevation' in k]):.0f}",
    }
    for k in oct_r["vals"]:
        for ax, a in (("azimuth", "Az"), ("elevation", "El")):
            macro[f"OctFit{a}{k}"] = f"{med[(k, ax, 'fit_error')][0]:.0f}"
            macro[f"OctCur{a}{k}"] = f"{med[(k, ax, 'fit_current')][0]:.0f}"
    dv = oct_m["drive_validation"]
    macro["CtrlErrDb"] = f"{max(dv[ax]['controller_dB'] for ax in AXES):.3f}"
    macro["CtrlErrDeg"] = f"{max(dv[ax]['controller_deg'] for ax in AXES):.2f}"
    macro["LinCheck"] = f"{max(oct_m['linearisation_check'].values()):.0e}"
    macro["NmovesOctAz"] = str(sum(m["axis"] == "azimuth" for m in json.loads(RESULTS["oct_time"].read_text())["moves"]))
    macro["NmovesOctEl"] = str(sum(m["axis"] == "elevation" for m in json.loads(RESULTS["oct_time"].read_text())["moves"]))
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in macro.items()))


def sequences_table(seqs_sep: dict, seqs_oct: dict) -> None:
    role = {"MPCLucy_pos_dist_dataset": "estimate of the effective gains",
            "MPCLucy_pos_dist_test": "validation", "MPCLucy_scotched": "validation"}
    rows = []
    for camp, seqs in (("Sept.", seqs_sep), ("Oct.", seqs_oct)):
        for n, s in seqs.items():
            mv = {ax: s.moves(ax) for ax in AXES}
            amp = sorted({round(abs(np.rad2deg(s.cmd[ax][s.t <= b][-1] - s.cmd[ax][s.t <= a][-1])), 1)
                          for ax in AXES for a, b in mv[ax]})
            rows.append([camp, tt(n), f"{s.t[-1]:.0f}", f"{len(mv['azimuth'])} / {len(mv['elevation'])}",
                         ", ".join(f"{x:g}" for x in amp), role.get(n, "comparison with the parameters of the files")])
    write("tab_sequences.tex", rows_tex(rows))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    res = {k: json.loads(p.read_text()) for k, p in RESULTS.items()}
    cfg_sep, cfg_oct = yaml.safe_load(CFG_SEP.read_text()), yaml.safe_load(CFG_OCT.read_text())
    drive_tables(cfg_sep, cfg_oct)
    seqs_sep = {n: read_sequence(SEP_DIR / f"{n}.csv", AMPS) for n in
                ("MPCLucy_pos_dist_dataset", "MPCLucy_pos_dist_test", "MPCLucy_scotched")}
    seqs_oct = {n: read_sequence(OCT_DIR / f"{n}.csv", AMPS) for n in
                yaml.safe_load((CL / "config_compare_vcl.yaml").read_text())["sequences"]}
    sequences_table(seqs_sep, seqs_oct)
    sep_r = september(res["sep"], seqs_sep)
    cl = idn._W["cl"]
    M, _, _ = cl.plant.linear_matrices()
    n_psi = cl.plant.sim.struttura.n_psi
    J = {"azimuth": float(np.sum(M[:n_psi, :n_psi])), "elevation": float(np.sum(M[n_psi:, n_psi:]))}
    current_scale(seqs_sep, J)
    oct_r = october(res["oct_time"], res["oct_main"])
    held_figure(sep_r, oct_r)
    zeta_figure()
    numbers(sep_r, oct_r, res["oct_main"], seqs_sep, seqs_oct)
    print("tables and figures written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "closed_loop_report.tex"],
                       cwd=HERE, check=True)


if __name__ == "__main__":
    main()
