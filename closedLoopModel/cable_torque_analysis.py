r"""Resisting torque of the axes against the angle and the history (VClosedLoop/V2 doublets).

Questions: does the resisting torque change with the angle, and is the torque at 0 after a doublet different
from the one before it (a symptom of a cable torque that cannot be modelled)? From the measured motor torque
:math:`K_t i` of the V2 files, taken in the order they were recorded (file modification times):

* cruise at 10 deg/s, file by file: half-sum of the two directions at the same angle (position-dependent
  torque) and half-difference (friction);
* holding torque at rest against the amplitude of the move before the rest;
* rests at 0 before and after each doublet, across the whole session (the axis is not moved between files);
* relaxation of the holding torque at constant position, against the change the integral action of the drive
  would give from the measured position error (gains of iteration I5 of vcl_study.py);
* neutral angle :math:`\theta_0` of the position torque fitted file by file with the identified model
  (v2_identification.py), against the one of the other session;

plus the figures of the friction model used (steady-state curve, torque-displacement loops).

Usage (from this folder)::

    ../.venv/bin/python -I cable_torque_analysis.py
"""

from __future__ import annotations

import json
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import current_input_study as cs  # noqa: E402
import main as mn  # noqa: E402
import v2_identification as v2  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl import plant_fit as pf  # noqa: E402
from lucy_cl.validation import read_sequence  # noqa: E402

plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 6.5,
                     "savefig.bbox": "tight"})
CONFIG = "config_v2_identification.yaml"
RES = HERE / "results_v2_identification" / "cables"
GAINS = HERE / "results_vcl_study" / "summary.json"   # drive gains of iteration I5 (last one)
ZERO_DEG = 0.05                                        # a rest within this angle of 0 is "at 0"


def short(name: str) -> str:
    return name.replace("V2/", "").replace("MPCLucy_doublets_", "").replace("MPCLucy_", "")


def amplitude_of(seq, ax: str, a: float, b: float) -> float:
    c = seq.cmd[ax]
    return float(np.rad2deg(c[seq.t <= b][-1] - c[seq.t <= a][-1]))


def relax(t, T_inf, dT, tau):
    return T_inf + dT * np.exp(-t / tau)


def rests_of(name: str, s, ax: str, kt: float, gain: dict, t_file0: float) -> list[dict]:
    """Rests of ``ax`` in a file: before the first move and after each move, with the torque 3 s after the move
    and at the end, the exponential relaxation, the encoder error and the change of torque the integral action
    of the drive gives from it."""
    mv = s.moves(ax)
    rows = []
    other = [a for a, _ in s.moves(v2.OTHER[ax])]
    starts = [a for a, _ in mv] + other
    segs = []
    if mv and mv[0][0] - s.t[0] > 1.5:
        segs.append((None, s.t[0], mv[0][0] - 0.2))
    for k, (a, b) in enumerate(mv):
        segs.append(((a, b), b, min([x for x in starts if x > b] + [s.t[-1]]) - 0.2))
    ki = gain["Kpv"] * 2 * np.pi * gain["Kiv"] * gain["Kpp"]   # current per position error and time, A/(rad s)
    for mvab, t0, t1 in segs:
        t_from = t0 + (3.0 if mvab else 0.0)
        if t1 - t_from < 1.0:
            continue
        m = (s.t >= t_from) & (s.t <= t1)
        T = kt * s.current[ax][m]
        e = s.cmd[ax][m] - s.enc[ax][m]
        tt = s.t[m] - t_from
        n1 = max(10, min(1000, T.size // 10))
        row = {"file": name, "axis": ax, "kind": "after a move" if mvab else "before the first move",
               "t_session_min": (t_file0 + t0) / 60.0, "angle_deg": float(np.rad2deg(np.mean(s.enc[ax][m]))),
               "last_move_deg": amplitude_of(s, ax, *mvab) if mvab else None, "rest_s": float(t1 - t_from),
               "T_start": float(np.mean(T[:n1])), "T_end": float(np.mean(T[-n1:])), "T_mean": float(np.mean(T)),
               "error_start_mdeg": float(out.MDEG * np.mean(e[:n1])), "error_end_mdeg": float(out.MDEG * np.mean(e[-n1:])),
               "encoder_change_mdeg": float(out.MDEG * (np.mean(s.enc[ax][m][-n1:]) - np.mean(s.enc[ax][m][:n1]))),
               "integral_action_dT": float(kt * ki * np.trapezoid(e, tt))}
        if mvab and t1 - t_from >= 10.0:
            try:
                p, _ = curve_fit(relax, tt, T, p0=[T[-1], T[0] - T[-1], 10.0],
                                 bounds=([-3, -3, 0.5], [3, 3, 500.0]), maxfev=5000)
                row.update(T_inf=float(p[0]), dT=float(p[1]), tau_s=float(p[2]))
            except RuntimeError:
                pass
        rows.append(row)
    return rows


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    fig_dir = RES / "figures"
    cc = yaml.safe_load((HERE / CONFIG).read_text())
    cfg = yaml.safe_load((HERE / cc["base_config"]).read_text())
    seq_dir = (HERE / cc["sequences_dir"]).resolve()
    ident = json.loads((HERE / cc["results_dir"] / "summary.json").read_text())
    final = ident["final"]
    gains = json.loads(GAINS.read_text())["iterations"][-1]["values"]
    cl, _, _ = mn.build(cfg, HERE)
    kt = {ax: cl.motors[ax].p.Kt for ax in AXES}
    js: dict = {"gains": gains}

    # ---- files of V2 in the order they were recorded (modification time = end of the recording)
    files = sorted((seq_dir / "V2").glob("MPCLucy_doublets_*.csv"), key=lambda f: f.stat().st_mtime)
    seqs = {f"V2/{f.stem}": read_sequence(f, 2.0) for f in files}
    t_end = {f"V2/{f.stem}": f.stat().st_mtime for f in files}
    t_ref = min(t_end[n] - float(s.t[-1]) for n, s in seqs.items())
    t_start = {n: t_end[n] - float(seqs[n].t[-1]) - t_ref for n in seqs}
    order = list(seqs)
    js["session"] = [{"file": n, "start_min": t_start[n] / 60.0, "duration_s": float(seqs[n].t[-1]),
                      "axis": [ax for ax in AXES if seqs[n].moves(ax)]} for n in order]

    # ---- rests
    rests = [r for n in order for ax in AXES if seqs[n].moves(ax)
             for r in rests_of(n, seqs[n], ax, kt[ax], gains[ax], t_start[n])]
    js["rests"] = rests

    # ---- cruise, file by file
    cruise = {}
    for ax in AXES:
        for n in order:
            if not seqs[n].moves(ax):
                continue
            cr = v2.cruise_torques({n: seqs[n]}, ax, kt[ax], erode_s=0.05, by_encoder=True)
            if not np.isfinite(cr["friction"]) or cr["angle_deg"].size < 500:
                continue
            th, d, T = cr["angle_deg"], cr["direction"], cr["torque"]
            bins = []
            for lo in np.arange(-15.0, 15.0, 1.0):
                m = (th >= lo) & (th < lo + 1.0)
                tp, tm = T[m & (d > 0)], T[m & (d < 0)]
                if tp.size > 30 and tm.size > 30:
                    bins.append([lo + 0.5, float(tp.mean()), float(tm.mean())])
            if len(bins) < 3:  # the two directions do not cross the same angles: friction and position not separable
                continue
            cruise.setdefault(ax, {})[n] = {"friction": cr["friction"], "gradient": cr["gradient"], "offset": cr["offset"],
                                             "rms": cr["rms"], "bins": bins}
    js["cruise"] = cruise
    js["cruise_all"] = {}
    for ax in AXES:
        cr = v2.cruise_torques({n: s for n, s in seqs.items() if n != "V2/MPCLucy_doublets_azimut_1deg_15deg_v2"}, ax,
                               kt[ax], erode_s=0.05, by_encoder=True)
        js["cruise_all"][ax] = {k: cr[k] for k in ("friction", "gradient", "offset", "rms", "bins")}

    # ---- the two sessions on the rests: T = K theta + c + h d with the gradient of the V2 cruise, torque 3 s after
    # the move (the rests of the other session are shorter); c is the position-dependent torque at 0
    old = {n: read_sequence(seq_dir / f"{n}.csv", 2.0) for n in cc["validation_sequences"]}
    js["sessions"] = {}
    for ax in AXES:
        K = js["cruise_all"][ax]["gradient"]
        for lab, S in (("V2", {n: seqs[n] for n in order}), ("other session", old)):
            R = [r for n, s in S.items() if s.moves(ax) for r in rests_of(n, s, ax, kt[ax], gains[ax], 0.0)
                 if r["last_move_deg"] and abs(r["last_move_deg"]) >= 0.9]
            y = np.array([r["T_start"] - K * np.deg2rad(r["angle_deg"]) for r in R])
            A = np.column_stack([np.ones(len(R)), [np.sign(r["last_move_deg"]) for r in R]])
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            js["sessions"].setdefault(ax, {})[lab] = {"rests": len(R), "c": float(c[0]), "h": float(c[1]),
                                                     "rms": float(np.std(y - A @ c)), "K_used": K}

    # ---- neutral angle of the position torque, file by file, with the identified model (full H6)
    rm = pf.RingMap.from_config(cl.plant.config)
    base = {ax: dict(final[ax]) for ax in AXES}
    pre, after = float(cc["window"]["pre_s"]), float(cc["window"]["after_move_s"])
    bounds = {k: tuple(float(x) for x in v) for k, v in cc["bounds"].items()}
    prev = json.loads((RES / "summary.json").read_text()) if (RES / "summary.json").exists() else {}
    th0 = prev.get("theta0_per_file", {})
    with mp.get_context("fork").Pool(int(cc.get("workers", 8)), initializer=cs._init,
                                     initargs=(cfg, str(seq_dir), 2.0, float(cc["input_lowpass_hz"]),
                                               {"a": rm.a, "b": rm.b}, base)) as pool:

        def run(wins, vals, rigid):
            return pool.map(v2.job, [(w, vals, 30.0, rigid) for w in wins])

        for n in order:
            for ax in AXES:
                wins = [w for w in cs.windows(seqs[n], n, pre, after) if w.axis == ax]
                if not wins or final[ax].get("K_cable", 0.0) <= 0.0 or n in th0.get(ax, {}):
                    continue
                grid = np.arange(-10.0, 10.01, 1.0)
                vals = []
                for g0 in grid:
                    v = {a: dict(final[a]) for a in AXES}
                    v[ax]["theta0_cable_deg"] = float(g0)
                    vals.append(v)
                res = pool.map(v2.job, [(w, v, 30.0, False) for v in vals for w in wins])
                nw = len(wins)
                costs = [v2.window_cost(res[k * nw:(k + 1) * nw], 1.0) for k in range(len(grid))]
                v = {a: dict(final[a]) for a in AXES}
                v[ax]["theta0_cable_deg"] = float(grid[int(np.argmin(costs))])
                v_n, *_ = v2.lsq(run, ax, wins, ["theta0_cable_deg"], v, bounds, 10, False,
                                 log=f"theta0 of {short(n)}", loss_scale=1.0)
                fits = [v2.fit_of(r) for r in run(wins, v_n, False)]
                fits_v2 = [v2.fit_of(r) for r in run(wins, final, False)]
                th0.setdefault(ax, {})[n] = {"theta0_deg": v_n[ax]["theta0_cable_deg"], "fit_median": float(np.median(fits)),
                                             "fit_median_session_theta0": float(np.median(fits_v2))}
    js["theta0_per_file"] = th0
    js["theta0_identified"] = {ax: final[ax].get("theta0_cable_deg") for ax in AXES}
    js["theta0_other_session"] = ident.get("theta0_per_sequence", {})
    (RES / "summary.json").write_text(json.dumps(js, indent=1))

    # --------------------------------------------------------------------------- tables
    md = ["# Resisting torque against the angle and the history (V2)", "",
          "Motor torque K_t i of the V2 files in the order they were recorded (file modification times). Rests: torque "
          "3 s after the move (start) and at the end of the rest; relaxation T(t) = T_inf + dT exp(-t/tau) from 3 s after "
          "the move; integral action: change of torque the velocity integrator of the drive gives from the measured "
          "position error, K_t K_pv 2 pi K_iv K_pp int(e dt) (gains of I5).", ""]
    md += ["## Session", "", out.markdown_table(["file", "start [min]", "duration [s]", "axis"],
                                                 [[short(r["file"]), f"{r['start_min']:.1f}", f"{r['duration_s']:.0f}",
                                                   ", ".join(r["axis"])] for r in js["session"]]), ""]
    md += ["## Rests", "", out.markdown_table(
        ["file", "axis", "rest", "angle [deg]", "last move [deg]", "length [s]", "T start [N m]", "T end [N m]",
         "tau [s]", "error start/end [mdeg]", "encoder change [mdeg]", "integral action dT [N m]"],
        [[short(r["file"]), r["axis"], r["kind"], f"{r['angle_deg']:+.2f}",
          "-" if r["last_move_deg"] is None else f"{r['last_move_deg']:+.2f}", f"{r['rest_s']:.0f}",
          f"{r['T_start']:+.3f}", f"{r['T_end']:+.3f}", f"{r['tau_s']:.0f}" if "tau_s" in r else "-",
          f"{r['error_start_mdeg']:+.2f} / {r['error_end_mdeg']:+.2f}", f"{r['encoder_change_mdeg']:+.2f}",
          f"{r['integral_action_dT']:+.3f}"] for r in rests]), ""]
    md += ["## Cruise, file by file", "", "T = F d + K theta + c over the cruise samples of each file (both directions).", "",
           out.markdown_table(["axis", "file", "friction F [N m]", "gradient K [N m/rad]", "offset c [N m]", "rms [N m]"],
                              [[ax, short(n), f"{c['friction']:.3f}", f"{c['gradient']:.3f}", f"{c['offset']:+.3f}",
                                f"{c['rms']:.3f}"] for ax, d in cruise.items() for n, c in d.items()]), ""]
    md += ["Cruise of all files (encoder speed within 3 % of 10 deg/s): " + "; ".join(
        f"{ax} F {c['friction']:.3f} N m, K {c['gradient']:.3f} N m/rad, c {c['offset']:+.3f} N m"
        for ax, c in js["cruise_all"].items()) + ".", ""]
    md += ["## The two sessions on the rests", "", "T(3 s) - K theta = c + h d over the rests after moves of 0.9 deg "
           "or more, K of the V2 cruise: c is the position-dependent torque at 0, h the friction left at the stop.", "",
           out.markdown_table(["axis", "session", "rests", "c [N m]", "h [N m]", "rms [N m]"],
                              [[ax, lab, str(v["rests"]), f"{v['c']:+.3f}", f"{v['h']:.3f}", f"{v['rms']:.3f}"]
                               for ax, d in js["sessions"].items() for lab, v in d.items()]), ""]
    md += ["## Neutral angle file by file", "", "theta0 fitted alone on the moves of each V2 file, all other parameters "
           "as identified on V2 (full H6 model).", "",
           out.markdown_table(["axis", "file", "theta0 [deg]", "fit median [%]", "fit median with theta0 of V2 [%]"],
                              [[ax, short(n), f"{v['theta0_deg']:+.2f}", f"{v['fit_median']:.1f}",
                                f"{v['fit_median_session_theta0']:.1f}"] for ax, d in th0.items() for n, v in d.items()]), "",
           "theta0 identified on all of V2: " + ", ".join(f"{ax} {v:+.2f} deg" for ax, v in js["theta0_identified"].items()
                                                         if v is not None) + "; other session (VClosedLoop steps): "
           + "; ".join(f"{ax} " + ", ".join(f"{x:+.2f}" for x in d.values()) for ax, d in js["theta0_other_session"].items())
           + " deg.", ""]
    (RES / "summary.md").write_text("\n".join(md))

    # --------------------------------------------------------------------------- figures
    # session timeline: torque and angle of the moving axis of each file
    fig, axs = plt.subplots(4, 1, figsize=(11, 9), sharex=False)
    for j, ax in enumerate(AXES):
        a_T, a_th = axs[2 * j], axs[2 * j + 1]
        for n in order:
            s = seqs[n]
            if not s.moves(ax):
                continue
            k = slice(None, None, 50)
            tm = (t_start[n] + s.t[k]) / 60.0
            a_T.plot(tm, kt[ax] * s.current[ax][k], lw=0.5, color="0.4")
            a_th.plot(tm, np.rad2deg(s.enc[ax][k]), lw=0.8, color="C0")
            a_th.text(tm[0], np.rad2deg(np.max(s.enc[ax])) + 1.0, short(n).split("_", 1)[1], fontsize=5.5)
        R = [r for r in rests if r["axis"] == ax]
        for r in R:
            at0 = abs(r["angle_deg"]) < ZERO_DEG
            a_T.plot(r["t_session_min"] + 3.0 / 60 * (r["kind"] == "after a move"), r["T_start"], "o" if at0 else ".",
                     color="C3" if at0 else "C1", ms=5 if at0 else 4)
        a_T.plot([], [], "-", color="0.4", lw=0.5, label="motor torque K_t i (every 50th sample)")
        a_T.plot([], [], "o", color="C3", label="start of a rest at 0 deg")
        a_T.plot([], [], ".", color="C1", label="start of a rest at +/-A")
        a_T.set_ylim(-1.2, 1.2)
        a_T.set_ylabel("torque [N m]")
        a_T.set_title(f"{ax}: V2 session in recording order")
        a_T.legend(loc="upper right", ncol=3)
        a_th.set_ylabel("encoder angle [deg]")
        a_th.set_xlabel("session time [min] (from the file modification times)")
    fig.tight_layout()
    out.save(fig, fig_dir, "session_timeline")

    # cruise file by file
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 7))
    for j, ax in enumerate(AXES):
        for i, (n, c) in enumerate(cruise.get(ax, {}).items()):
            B = np.array(c["bins"])
            if not B.size:
                continue
            col = f"C{i}"
            axs[0, j].plot(B[:, 0], 0.5 * (B[:, 1] + B[:, 2]), "-o", ms=2.5, color=col, label=short(n))
            axs[1, j].plot(B[:, 0], 0.5 * (B[:, 1] - B[:, 2]), "-o", ms=2.5, color=col, label=short(n))
        axs[0, j].set_title(f"{ax}: position-dependent torque (half-sum of the two directions)")
        axs[1, j].set_title(f"{ax}: friction at 10 deg/s (half-difference)")
        for a in axs[:, j]:
            a.set_xlabel("encoder angle [deg]")
            a.set_ylabel("torque [N m]")
            a.legend(loc="best")
        axs[0, j].set_ylim(-0.5, 0.5)
        axs[1, j].set_ylim(0.0, 0.7)
    fig.tight_layout()
    out.save(fig, fig_dir, "cruise_per_file")

    # holding torque against the amplitude of the last move
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4))
    for j, ax in enumerate(AXES):
        a = axs[j]
        cr = v2.cruise_torques({n: s for n, s in seqs.items()}, ax, kt[ax], erode_s=0.05, by_encoder=True)
        for r in [r for r in rests if r["axis"] == ax and r["last_move_deg"]]:
            pos = cr["gradient"] * np.deg2rad(r["angle_deg"]) + cr["offset"] if ax == "elevation" else 0.0
            at0 = abs(r["angle_deg"]) < ZERO_DEG
            a.semilogx(abs(r["last_move_deg"]), np.sign(r["last_move_deg"]) * (r["T_start"] - pos),
                       "o" if at0 else "s", color="C3" if r["last_move_deg"] > 0 else "C0", mfc="none" if at0 else None, ms=5)
        for mk, col, lab in (("s", "C3", "positive move, rest at +/-A"), ("o", "C3", "positive move, rest at 0"),
                             ("s", "C0", "negative move, rest at +/-A"), ("o", "C0", "negative move, rest at 0")):
            a.plot([], [], mk, color=col, mfc="none" if mk == "o" else None, label=lab)
        a.axhline(cr["friction"], color="k", ls="--", lw=0.8, label=f"friction at 10 deg/s ({cr['friction']:.2f} N m)")
        a.set_xlabel("amplitude of the move before the rest [deg]")
        a.set_ylabel("holding torque 3 s after the move x sign of the move [N m]")
        a.set_title(f"{ax}: holding torque against the last move"
                    + (" (position torque of the cruise removed)" if ax == "elevation" else ""))
        a.legend(loc="lower right")
    fig.tight_layout()
    out.save(fig, fig_dir, "holding_vs_amplitude")

    # relaxation during the rests
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 7))
    for j, ax in enumerate(AXES):
        R = [r for r in rests if r["axis"] == ax and r["kind"] == "after a move" and r["rest_s"] >= 10.0]
        amps = np.array([abs(r["last_move_deg"]) for r in R])
        cmap = plt.get_cmap("viridis")
        norm = matplotlib.colors.LogNorm(vmin=0.1, vmax=30)
        for r in R:
            s = seqs[r["file"]]
            t0 = r["t_session_min"] * 60.0 - t_start[r["file"]] + 3.0
            m = (s.t >= t0) & (s.t <= t0 + r["rest_s"])
            k = slice(None, None, 20)
            col = cmap(norm(abs(r["last_move_deg"])))
            axs[0, j].plot(s.t[m][k] - t0, kt[ax] * s.current[ax][m][k] - r["T_start"], color=col, lw=0.7)
            axs[1, j].plot(s.t[m][k] - t0, out.MDEG * (s.enc[ax][m][k] - s.cmd[ax][m][k]), color=col, lw=0.7)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        fig.colorbar(sm, ax=axs[:, j], label="amplitude of the move before the rest [deg]")
        axs[0, j].set_title(f"{ax}: holding torque during the rests ({len(R)} rests of 10 s or more)")
        axs[0, j].set_ylabel("torque change from 3 s after the move [N m]")
        axs[1, j].set_title(f"{ax}: encoder minus command during the same rests")
        axs[1, j].set_ylabel("position error [mdeg]")
        for a in axs[:, j]:
            a.set_xlabel("time from 3 s after the move [s]")
        if amps.size:
            axs[1, j].set_ylim(-2, 2)
    out.save(fig, fig_dir, "relaxation")

    # integral action against the observed relaxation
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4))
    for j, ax in enumerate(AXES):
        R = [r for r in rests if r["axis"] == ax and r["kind"] == "after a move" and r["rest_s"] >= 10.0]
        x = [r["integral_action_dT"] for r in R]
        y = [r["T_end"] - r["T_start"] for r in R]
        axs[j].plot(x, y, "o", color="C0", label="rests of 10 s or more")
        lim = max(0.05, max(np.abs(x + y)) * 1.1) if R else 0.3
        axs[j].plot([-lim, lim], [-lim, lim], "k--", lw=0.8, label="observed = integral action")
        axs[j].set_xlabel("torque change of the drive integral action from the measured error [N m]")
        axs[j].set_ylabel("observed torque change over the rest [N m]")
        axs[j].set_title(f"{ax}: is the relaxation the drive integrator?")
        axs[j].legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "relaxation_vs_integral_action")

    # neutral angle file by file
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4))
    for j, ax in enumerate(AXES):
        a = axs[j]
        d = th0.get(ax, {})
        ks = list(d)
        a.plot(range(len(ks)), [d[k]["theta0_deg"] for k in ks], "o", color="C0", label="V2, fitted file by file")
        if js["theta0_identified"].get(ax) is not None:
            a.axhline(js["theta0_identified"][ax], color="C0", ls="--", lw=0.8, label="V2, fitted on all files")
        o = js["theta0_other_session"].get(ax, {})
        a.plot(np.arange(len(ks), len(ks) + len(o)), list(o.values()), "s", color="C3",
               label="other session (VClosedLoop steps), file by file")
        a.set_xticks(range(len(ks) + len(o)), [short(k) for k in ks] + [short(k) for k in o], rotation=60, fontsize=6)
        a.set_ylabel("neutral angle theta0 [deg]")
        a.set_title(f"{ax}: neutral angle of the position torque (K_c = {final[ax].get('K_cable', 0):.2f} N m/rad)")
        a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "theta0_per_file")

    # friction model: steady-state curve and torque-displacement loops of the identified elements
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    v = np.logspace(-5, 0, 400)
    for j, ax in enumerate(AXES):
        f = final[ax]
        Tc, Ts, vs_ = f["Tc"], f["Tc"] * f["Ts_ratio"], f["vs"]
        b = f.get("b_{\\theta,h}", 0.0)
        g = Tc + (Ts - Tc) * np.exp(-(v / vs_) ** 2) + b * v
        axs[0].semilogx(np.rad2deg(v), g, color=f"C{j}", label=f"{ax}: T_c {Tc:.2f}, T_s {Ts:.2f} N m, v_s {np.rad2deg(vs_):.2f} deg/s")
        for amp_deg, ls in ((0.005, ":"), (0.05, "--"), (1.0, "-")):
            A = np.deg2rad(amp_deg)
            tt = np.linspace(0, 3, 6000)
            x = A * np.sin(2 * np.pi * tt)   # 1 Hz sinusoid, three periods
            xd = A * 2 * np.pi * np.cos(2 * np.pi * tt)
            T = np.zeros_like(tt)
            for k in range(1, tt.size):      # Dahl with the Stribeck level, explicit sub-steps
                dt, Tk = (tt[k] - tt[k - 1]) / 50, T[k - 1]
                for _ in range(50):
                    w = xd[k - 1]
                    lev = Tc + (Ts - Tc) * np.exp(-(w / vs_) ** 2)
                    Tk += dt * f["sigma0"] * (w - abs(w) * Tk / lev)
                T[k] = Tk
            sel = tt >= 1.0
            axs[1 + j].plot(np.rad2deg(x[sel]) / amp_deg, T[sel] + b * xd[sel], ls, color=f"C{j}",
                            label=f"amplitude {amp_deg:g} deg")
        axs[1 + j].set_xlabel("displacement / amplitude [-] (1 Hz sinusoid)")
        axs[1 + j].set_ylabel("friction torque [N m]")
        axs[1 + j].set_title(f"{ax}: Dahl loops, sigma0 {f['sigma0']:.0f} N m/rad")
        axs[1 + j].legend(loc="lower right")
    axs[0].set_xlabel("sliding speed [deg/s]")
    axs[0].set_ylabel("steady-state friction torque [N m]")
    axs[0].set_title("steady-state friction (Stribeck curve + viscous)")
    axs[0].legend(loc="upper right")
    fig.tight_layout()
    out.save(fig, fig_dir, "friction_model")
    print("done", flush=True)


if __name__ == "__main__":
    main()
