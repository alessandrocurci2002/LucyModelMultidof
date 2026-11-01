r"""Plant driven by the measured currents, identified on the doublets of VClosedLoop/V2.

The measured motor currents drive the plant (:class:`lucy_cl.current_driven.CurrentDrivenPlant`, no drive in
the loop); the goodness measure is the NRMSE fit of the encoder position of the moving axis. Steps:

0. model-free analysis of V2: holding torque at rest against the angle and the direction of the last move;
   torque in the cruise at 10 deg/s, whose half-difference between the two directions is the friction at that
   speed and whose half-sum is the position-dependent torque (cables, imbalance); start values of the fits;
1. check of the rigid surrogate (:meth:`~lucy_cl.current_driven.CurrentDrivenPlant.simulate_window_rigid`)
   against the full H6 model on all the identification windows;
2. model structures fitted on all of V2 and scored by leave-one-file-out cross-validation (rigid surrogate):
   cable spring :math:`-K_c(\theta - \theta_0)` yes/no; then each of :math:`\sigma_0`, :math:`v_s`,
   :math:`b_{\theta,h}` fixed at a nominal value (x 1/3, 1, 3) instead of fitted; then the accepted fixings
   together. The structure with the best cross-validation median is selected, a simpler one within
   ``parsimony_points``;
3. selected structure refitted with the full H6 model on all of V2 (warm start), standard errors;
4. validation with the full model on the VClosedLoop steps (another session).

The results are saved after every structure: a stopped run resumes where it stopped.

Usage (from this folder)::

    ../.venv/bin/python -I v2_identification.py --config config_v2_identification.yaml
"""

from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import current_input_study as cs  # noqa: E402
import main as mn  # noqa: E402
from identifiability import effective_samples  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl import plant_fit as pf  # noqa: E402
from lucy_cl.validation import read_sequence  # noqa: E402

MDEG = out.MDEG
OTHER = {"azimuth": "elevation", "elevation": "azimuth"}
SIGNED = {"theta0_cable_deg"}
"""Parameters fitted on their value (the others on their logarithm)."""
plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 6.5,
                     "savefig.bbox": "tight"})


# --------------------------------------------------------------------------- workers


def job(args) -> dict:
    """Position of one window with the parameter values of both axes; ``rigid``: rigid surrogate."""
    win, values, max_cpu, rigid = args
    cs._apply(values)
    cd = cs._W["cd"]
    sim = cd.simulate_window_rigid if rigid else cd.simulate_window
    try:
        r = sim(cs._seq(win.sequence), win.t0, win.t1, max_cpu_s=max_cpu)
        return {"t": r["t"], "ok": r["ok"], "axis": win.axis,
                **{ax: {"enc_meas": r[ax]["enc_meas"], "enc": r[ax]["enc"]} for ax in AXES}}
    finally:
        cs._apply({ax: {} for ax in AXES})


def fit_of(r: dict) -> float:
    ax = r["axis"]
    f = cs.fit_nrmse(r[ax]["enc_meas"], r[ax]["enc"])
    return f if np.isfinite(f) else -999.0


# --------------------------------------------------------------------------- 0. model-free analysis


def rest_torques(seqs: dict, ax: str, kt: float, settle_s: float = 3.0) -> list[dict]:
    """Mean motor torque :math:`K_t i` in the rests of ``ax`` (from ``settle_s`` after a move to the next move),
    with the angle and the direction of the move before the rest."""
    rows = []
    for name, s in seqs.items():
        mv = s.moves(ax)
        if not mv:
            continue
        starts = [a for a, _ in mv] + [a for a, _ in s.moves(OTHER[ax])]
        for a, b in mv:
            t0, t1 = b + settle_s, min([x for x in starts if x > b] + [s.t[-1]]) - 0.2
            if t1 - t0 < 1.0:
                continue
            m = (s.t >= t0) & (s.t <= t1)
            n1 = min(2000, int(m.sum() // 3))
            T = kt * s.current[ax][m]
            rows.append({"sequence": name, "angle_deg": float(np.rad2deg(np.mean(s.enc[ax][m]))),
                         "direction": float(np.sign(s.cmd[ax][s.t <= b][-1] - s.cmd[ax][s.t <= a][-1])),
                         "torque": float(np.mean(T)), "torque_start": float(np.mean(T[:n1])),
                         "torque_end": float(np.mean(T[-n1:])), "rest_s": float(t1 - t0)})
    return rows


def cruise_torques(seqs: dict, ax: str, kt: float, v_deg_s: float = 10.0, erode_s: float = 0.15,
                   by_encoder: bool = False) -> dict:
    """Motor torque in the cruise of the moves (command speed ``v_deg_s``, from ``erode_s`` after its start to
    ``erode_s`` before its end): angle, direction, torque per sample. ``by_encoder``: the cruise is where the
    encoder speed (low-passed at 20 Hz) is within 3 % of ``v_deg_s``, which leaves out the acceleration of the
    axis behind the command (following error of about v/K_pp)."""
    from scipy.ndimage import binary_erosion
    from scipy.signal import butter, filtfilt

    th, d, T = [], [], []
    for s in seqs.values():
        if not s.moves(ax):
            continue
        dt = float(np.median(np.diff(s.t)))
        b, a = butter(2, 20.0 * 2 * dt)
        vc = np.gradient(filtfilt(b, a, s.enc[ax]), s.t) if by_encoder else np.gradient(s.cmd[ax], s.t)
        cr = np.abs(np.abs(vc) - np.deg2rad(v_deg_s)) < np.deg2rad((0.03 if by_encoder else 0.02) * v_deg_s)
        cr = binary_erosion(cr, iterations=int(erode_s / dt))
        if cr.any():
            th.append(np.rad2deg(s.enc[ax][cr]))
            d.append(np.sign(vc[cr]))
            T.append(kt * filtfilt(b, a, s.current[ax])[cr])
    if not th:  # no move long enough to reach the cruise
        e = np.zeros(0)
        return {"angle_deg": e, "direction": e, "torque": e, "friction": float("nan"), "gradient": 0.0, "offset": 0.0,
                "rms": float("nan"), "bins": []}
    th, d, T = (np.concatenate(x) for x in (th, d, T))
    A = np.column_stack([d, np.deg2rad(th), np.ones_like(th)])
    c, *_ = np.linalg.lstsq(A, T, rcond=None)
    edges = np.arange(-16.0, 16.1, 2.0)
    bins = []
    for lo in edges[:-1]:
        m = (th >= lo) & (th < lo + 2.0)
        tp, tm = T[m & (d > 0)], T[m & (d < 0)]
        if tp.size > 50 and tm.size > 50:
            bins.append([lo + 1.0, float(tp.mean()), float(tm.mean())])
    return {"angle_deg": th, "direction": d, "torque": T, "friction": float(c[0]), "gradient": float(c[1]),
            "offset": float(c[2]), "rms": float(np.std(T - A @ c)), "bins": bins}


def window_cost(results: list[dict], loss_scale: float | None) -> float:
    """Sum over the windows of (1 - fit)^2, or of the Cauchy loss of it (see :func:`lsq`)."""
    e2 = np.array([np.sum(cs.residual_of(r) ** 2) for r in results])
    return float(np.sum(e2 if not loss_scale else loss_scale**2 * np.log1p(e2 / loss_scale**2)))


# --------------------------------------------------------------------------- structures


def free_count(st: dict) -> int:
    return len(st["free"])


def make_values(base: dict, ax: str, st: dict, fitted: dict | None = None) -> dict:
    """Values of both axes for structure ``st`` of ``ax``: base, then the fixed values, then ``fitted``."""
    v = {a: dict(base[a]) for a in AXES}
    v[ax].update(st["fixed"])
    v[ax].update(fitted or {})
    return v


def lsq(run, ax: str, wins: list, names: list[str], values: dict, bounds: dict, max_nfev: int, rigid: bool,
        log: str = "", loss_scale: float | None = None):
    """Least squares on the position residuals of ``wins``: ``names`` of ``ax`` free, the rest of ``values``
    fixed. Returns the values, the solution and the standard errors (relative for the logarithmic
    parameters, absolute for the signed ones) from the Jacobian with the effective number of samples."""
    tr = lambda n, x: x if n in SIGNED else np.log(x)  # noqa: E731
    itr = lambda n, x: x if n in SIGNED else np.exp(x)  # noqa: E731
    lo = np.array([tr(n, bounds[n][0]) for n in names])
    hi = np.array([tr(n, bounds[n][1]) for n in names])
    x0 = np.clip(np.array([tr(n, values[ax][n]) for n in names]), lo + 1e-9, hi - 1e-9)

    def to_values(x: np.ndarray) -> dict:
        v = {a: dict(values[a]) for a in AXES}
        v[ax].update({n: float(itr(n, xi)) for n, xi in zip(names, x)})
        return v

    def robust(rs: list[np.ndarray]) -> np.ndarray:
        # per-window Cauchy loss c^2 log(1 + e^2/c^2) of the window error e = 1 - fit: a window the model loses
        # (drift after the move, fit far below zero) does not drive the estimate
        if not loss_scale:
            return np.concatenate(rs)
        out_ = []
        for r in rs:
            e = float(np.linalg.norm(r))
            out_.append(r * (np.sqrt(loss_scale**2 * np.log1p(e * e / loss_scale**2)) / e if e > 0 else 1.0))
        return np.concatenate(out_)

    def resid(x: np.ndarray) -> np.ndarray:
        return robust([cs.residual_of(r) for r in run(wins, to_values(x), rigid)])

    sol = least_squares(resid, x0, bounds=(lo, hi), diff_step=0.03, max_nfev=max_nfev, xtol=1e-3, ftol=1e-4,
                        x_scale=[10.0 if n in SIGNED else 1.0 for n in names])
    rs = [cs.residual_of(r) for r in run(wins, to_values(sol.x), rigid)]
    r0 = robust(rs)
    n_eff = sum(effective_samples(r) for r in rs)
    s2 = float(r0 @ r0) / max(1.0, n_eff - len(names))
    F = sol.jac.T @ sol.jac * (n_eff / r0.size)
    se = np.sqrt(np.clip(np.diag(s2 * np.linalg.pinv(F)), 0, None))
    print(f"    {log} {ax}: " + ", ".join(f"{n} {itr(n, e):.4g}" for n, e in zip(names, sol.x))
          + f" -> cost {float(r0 @ r0):.4f} ({sol.nfev} evaluations)", flush=True)
    return to_values(sol.x), sol, dict(zip(names, se.tolist())), float(r0 @ r0)


# --------------------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config_v2_identification.yaml")
    args = ap.parse_args()
    t_start = time.time()
    cc_path = Path(args.config).resolve()
    cc = yaml.safe_load(cc_path.read_text())
    cfg = yaml.safe_load((HERE / cc["base_config"]).read_text())
    res_dir = HERE / cc["results_dir"]
    fig_dir = res_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    js_path = res_dir / "summary.json"
    js = json.loads(js_path.read_text()) if js_path.exists() else {}
    save = lambda: js_path.write_text(json.dumps(js, indent=1))  # noqa: E731
    seq_dir = (HERE / cc["sequences_dir"]).resolve()
    amps = float(cc.get("amps_per_A", 2.0))
    mc = float(cc.get("max_cpu_per_window_s", 30))
    fl = lambda v: [float(x) for x in v]  # noqa: E731  (PyYAML: 1.0e4 without sign is a string)
    bounds = {k: tuple(fl(v)) for k, v in cc["bounds"].items()}
    pre, after = float(cc["window"]["pre_s"]), float(cc["window"]["after_move_s"])
    ls = cc.get("window_loss_scale")
    ls = None if ls is None else float(ls)

    cl, _, _ = mn.build(cfg, HERE)
    rm = pf.RingMap.from_config(cl.plant.config)
    kt = {ax: cl.motors[ax].p.Kt for ax in AXES}
    dahl = {ax: [e for e in cl.plant.friction if e.axis == ax][0] for ax in AXES}
    seq_id = {n: read_sequence(seq_dir / f"{n}.csv", amps) for n in cc["identification_sequences"]}
    seq_val = {n: read_sequence(seq_dir / f"{n}.csv", amps) for n in cc["validation_sequences"]}

    def amplitude(seq, w) -> float:
        c = seq.cmd[w.axis]
        return float(np.rad2deg(c[seq.t <= w.end][-1] - c[seq.t <= w.move][-1]))

    wid = {ax: [w for n, s in seq_id.items() for w in cs.windows(s, n, pre, after) if w.axis == ax] for ax in AXES}
    wval = [w for n, s in seq_val.items() for w in cs.windows(s, n, pre, after)]
    amp = {w.label: amplitude(seq_id.get(w.sequence) or seq_val[w.sequence], w) for ax in AXES for w in wid[ax]}
    amp.update({w.label: amplitude(seq_val[w.sequence], w) for w in wval})
    short = lambda lab: lab.replace("V2/MPCLucy_doublets_", "").replace("MPCLucy_", "")  # noqa: E731
    print(f"windows: identification azimuth {len(wid['azimuth'])}, elevation {len(wid['elevation'])}; "
          f"validation {len(wval)}", flush=True)

    # ---- 0. model-free analysis of V2
    mf = {}
    for ax in AXES:
        rests = rest_torques(seq_id, ax, kt[ax])
        A = np.array([[np.deg2rad(r["angle_deg"]), 1.0, r["direction"]] for r in rests])
        y = np.array([r["torque"] for r in rests])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        cr = cruise_torques(seq_id, ax, kt[ax])
        K = max(cr["gradient"], 0.05)
        th0 = float(np.clip(-np.rad2deg(cr["offset"] / K), -10.0, 10.0))
        mf[ax] = {"rests": rests, "rest_fit": {"gradient": float(c[0]), "offset": float(c[1]), "direction": float(c[2]),
                                               "rms": float(np.std(y - A @ c))},
                  "cruise": {k: v for k, v in cr.items() if k not in ("angle_deg", "direction", "torque")},
                  "start": {"Tc": cr["friction"] if np.isfinite(cr["friction"]) else float(c[2]), "K_cable": K,
                            "theta0_cable_deg": th0}, "_cruise": cr}
        print(f"  {ax}: cruise friction {cr['friction']:.3f} N m, gradient {cr['gradient']:.3f} N m/rad, offset "
              f"{cr['offset']:+.3f} N m; rests: gradient {c[0]:.3f} N m/rad, direction term {c[2]:.3f} N m", flush=True)
    js["model_free"] = {ax: {k: v for k, v in d.items() if not k.startswith("_")} for ax, d in mf.items()}

    # base values: nominal values of the fixable parameters, start values of the others, no spring
    nominal = {ax: {} for ax in AXES}
    for ax in AXES:
        for n in cc["fixable"][ax]:
            v = cc["nominal"][n]
            nominal[ax][n] = float(dahl[ax].sigma0) if v == "bode" else float(v)
    base = {}
    for ax in AXES:
        b = {"scale": float(cc["start"]["scale"]), "Tc": mf[ax]["start"]["Tc"], "Ts_ratio": float(cc["start"]["Ts_ratio"]),
             "sigma0": nominal[ax]["sigma0"], "vs": nominal[ax]["vs"], "K_cable": 0.0, "theta0_cable_deg": 0.0}
        if "b_{\\theta,h}" in cc["fixable"][ax]:
            b["b_{\\theta,h}"] = nominal[ax]["b_{\\theta,h}"]
        base[ax] = b
    js["nominal"], js["base"] = nominal, base
    save()

    with mp.get_context("fork").Pool(int(cc.get("workers", 8)), initializer=cs._init,
                                     initargs=(cfg, str(seq_dir), amps, float(cc["input_lowpass_hz"]),
                                               {"a": rm.a, "b": rm.b}, base)) as pool:

        def run(wins: list, vals, rigid: bool) -> list[dict]:
            """``vals``: values of both axes, or a function of the window giving them."""
            return pool.map(job, [(w, vals(w) if callable(vals) else vals, mc, rigid) for w in wins])

        # ---- 1. rigid surrogate against the full model
        if "surrogate" not in js:
            js["surrogate"] = {}
            for ax in AXES:
                v = {a: dict(base[a]) for a in AXES}
                if mf[ax]["cruise"]["gradient"] > 0.3:
                    v[ax].update(K_cable=mf[ax]["start"]["K_cable"], theta0_cable_deg=mf[ax]["start"]["theta0_cable_deg"])
                t0 = time.time()
                rr = run(wid[ax], v, True)
                t_r = time.time() - t0
                t0 = time.time()
                rf = run(wid[ax], v, False)
                t_f = time.time() - t0
                df = [abs(fit_of(a) - fit_of(b)) for a, b in zip(rr, rf)]
                dp = [float(np.nanmax(np.abs(a[ax]["enc"] - b[ax]["enc"])) / max(np.ptp(b[ax]["enc_meas"]), 1e-9))
                      for a, b in zip(rr, rf)]
                js["surrogate"][ax] = {"fit_difference": df, "position_difference_rel": dp,
                                       "wall_s": {"rigid": t_r, "full": t_f}}
                print(f"  surrogate {ax}: fit difference median {np.median(df):.2f}, max {max(df):.2f} points; "
                      f"wall time rigid {t_r:.0f} s, full {t_f:.0f} s", flush=True)
            save()

        # ---- 2. structures, fitted on all of V2 and cross-validated by file
        js.setdefault("structures", {ax: {} for ax in AXES})
        par = float(cc["parsimony_points"])

        def evaluate(ax: str, key: str, st: dict, start: dict | None) -> dict:
            """Fit of ``st`` on all of V2 (grid start when ``start`` is None) and leave-one-file-out CV."""
            if key in js["structures"][ax]:
                return js["structures"][ax][key]
            print(f"  [{ax}] {key}: {st['label']}", flush=True)
            names = st["free"]
            if start is None:  # grid on the start
                g = cc["grid"]
                axes_g = [("scale", fl(g["scale"])), ("Tc", [f * mf[ax]["start"]["Tc"] for f in fl(g["Tc_factor"])]),
                          ("Ts_ratio", fl(g["Ts_ratio"]))]
                if "K_cable" in names:
                    axes_g += [("K_cable", [f * mf[ax]["start"]["K_cable"] for f in fl(g["K_cable_factor"])]),
                               ("theta0_cable_deg", [mf[ax]["start"]["theta0_cable_deg"] + o for o in fl(g["theta0_offset_deg"])])]
                pts = list(itertools.product(*[v for _, v in axes_g]))
                vals = [make_values(base, ax, st, dict(zip([k for k, _ in axes_g], p))) for p in pts]
                for v in vals:  # free parameters whose nominal value is zero start from the configured value
                    v[ax].update({n: float(cc["start"][n]) for n in names if n not in SIGNED and v[ax][n] <= 0.0})
                res = pool.map(job, [(w, v, mc, True) for v in vals for w in wid[ax]])
                nw = len(wid[ax])
                costs = [window_cost(res[k * nw:(k + 1) * nw], ls) for k in range(len(pts))]
                kb = int(np.argmin(costs))
                v0 = vals[kb]
                print(f"    grid of {len(pts)} points: best {dict(zip([k for k, _ in axes_g], np.round(pts[kb], 4)))}, "
                      f"cost {costs[kb]:.3f}", flush=True)
            else:
                v0 = make_values(base, ax, st, {k: x for k, x in start.items() if k in names})
            v_all, sol, se, cost = lsq(run, ax, wid[ax], names, v0, bounds, int(cc["max_nfev"]), True, log=key, loss_scale=ls)
            fits_all = [fit_of(r) for r in run(wid[ax], v_all, True)]
            # leave-one-file-out cross-validation, warm start from the fit on all of V2
            cv_fits, cv_vals = {}, {}
            for f_out in sorted({w.sequence for w in wid[ax]}):
                w_tr = [w for w in wid[ax] if w.sequence != f_out]
                w_te = [w for w in wid[ax] if w.sequence == f_out]
                v_tr, *_ = lsq(run, ax, w_tr, names, v_all, bounds, int(cc["max_nfev_cv"]), True,
                                log=f"{key} without {short(f_out)}", loss_scale=ls)
                for w, r in zip(w_te, run(w_te, v_tr, True)):
                    cv_fits[w.label] = fit_of(r)
                cv_vals[f_out] = {n: v_tr[ax][n] for n in names}
            cvf = np.array(list(cv_fits.values()))
            rec = {"label": st["label"], "free": names, "fixed": st["fixed"], "values": v_all[ax], "se": se,
                   "cost": cost, "nfev": int(sol.nfev), "fits": dict(zip([w.label for w in wid[ax]], fits_all)),
                   "cv_fits": cv_fits, "cv_values": cv_vals, "cv_median": float(np.median(cvf)),
                   "cv_mean_square": float(np.mean((1 - np.clip(cvf, -500, 100) / 100) ** 2)),
                   "fit_median": float(np.median(fits_all))}
            js["structures"][ax][key] = rec
            save()
            print(f"    -> identification median {rec['fit_median']:.1f} %, CV median {rec['cv_median']:.1f} %", flush=True)
            return rec

        def choose(ax: str, keys: list[str]) -> str:
            recs = {k: js["structures"][ax][k] for k in keys}
            best = max(recs, key=lambda k: recs[k]["cv_median"])
            near = [k for k in recs if recs[k]["cv_median"] >= recs[best]["cv_median"] - par]
            return min(near, key=lambda k: (len(recs[k]["free"]), -recs[k]["cv_median"]))

        js.setdefault("selection", {})
        for ax in AXES:
            fix = cc["fixable"][ax]
            friction = list(cc["always_free"]) + list(fix)
            # step 1: cable spring or not, all the friction parameters fitted
            s_no = {"label": "Dahl + Stribeck, all fitted, no cable spring", "free": friction, "fixed": {"K_cable": 0.0}}
            s_sp = {"label": "Dahl + Stribeck + cable spring, all fitted", "free": friction + list(cc["spring"]), "fixed": {}}
            evaluate(ax, "S0 no spring", s_no, None)
            evaluate(ax, "S1 spring", s_sp, None)
            k1 = choose(ax, ["S0 no spring", "S1 spring"])
            st1 = s_no if k1.startswith("S0") else s_sp
            print(f"  [{ax}] step 1: {k1}", flush=True)
            # step 2: each fixable parameter fixed at nominal x multiplier
            accepted = {}
            keys2 = [k1]
            for n in fix:
                mults = [1.0] if nominal[ax][n] == 0.0 else fl(cc["multipliers"])
                ks = []
                for mu in mults:
                    k = f"{k1[:2]} {n} fixed x{mu:.3g}" if nominal[ax][n] else f"{k1[:2]} {n} fixed 0"
                    st = {"label": f"{st1['label']}; {n} fixed at {mu:.3g} x nominal ({nominal[ax][n] * mu:.4g})",
                          "free": [x for x in st1["free"] if x != n], "fixed": {**st1["fixed"], n: nominal[ax][n] * mu}}
                    evaluate(ax, k, st, js["structures"][ax][k1]["values"])
                    ks.append(k)
                kb = max(ks, key=lambda k: js["structures"][ax][k]["cv_median"])
                keys2 += ks
                if js["structures"][ax][kb]["cv_median"] >= js["structures"][ax][k1]["cv_median"] - par:
                    accepted[n] = js["structures"][ax][kb]["fixed"][n]
            # step 3: the accepted fixings together
            keys3 = []
            if len(accepted) >= 2:
                k = f"{k1[:2]} " + ", ".join(accepted) + " fixed"
                st = {"label": f"{st1['label']}; " + ", ".join(f"{n} fixed at {v:.4g}" for n, v in accepted.items()),
                      "free": [x for x in st1["free"] if x not in accepted], "fixed": {**st1["fixed"], **accepted}}
                evaluate(ax, k, st, js["structures"][ax][k1]["values"])
                keys3 = [k]
            sel = choose(ax, [k1] + keys2[1:] + keys3)
            js["selection"][ax] = {"step1": k1, "accepted_fixings": accepted, "selected": sel,
                                   "compared": [k1] + keys2[1:] + keys3}
            save()
            print(f"  [{ax}] selected: {sel}", flush=True)

        # ---- 3. selected structure refitted with the full model on all of V2
        final = {a: dict(base[a]) for a in AXES}
        for ax in AXES:
            rec = js["structures"][ax][js["selection"][ax]["selected"]]
            final[ax].update(rec["fixed"])
            final[ax].update(rec["values"])
        if "full" not in js:
            js["full"] = {}
            for ax in AXES:
                rec = js["structures"][ax][js["selection"][ax]["selected"]]
                v_full, sol, se, cost = lsq(run, ax, wid[ax], rec["free"], final, bounds, int(cc["max_nfev_full"]),
                                            False, log="full model", loss_scale=ls)
                js["full"][ax] = {"values": v_full[ax], "se": se, "cost": cost, "nfev": int(sol.nfev)}
                final[ax] = v_full[ax]
                save()
        for ax in AXES:
            final[ax] = dict(js["full"][ax]["values"])
        js["final"] = final

        # ---- 4. fits with the full model: identification (V2) and validation (VClosedLoop steps)
        r_id = {ax: run(wid[ax], final, False) for ax in AXES}
        r_val = run(wval, final, False)
        js["fits"] = {"identification": {w.label: fit_of(r) for ax in AXES for w, r in zip(wid[ax], r_id[ax])},
                      "validation": {w.label: fit_of(r) for w, r in zip(wval, r_val)}}
        # validation with the neutral angle of the cables of each sequence estimated alone (when the spring is in)
        js["fits"]["validation_theta0_per_sequence"] = {}
        js["theta0_per_sequence"] = {}
        for ax in AXES:
            if final[ax].get("K_cable", 0.0) <= 0.0:
                continue
            for sq in cc["validation_sequences"]:
                w_sq = [w for w in wval if w.axis == ax and w.sequence == sq]
                if not w_sq:
                    continue
                grid = np.arange(-10.0, 10.01, 1.0)
                vals = []
                for g0 in grid:
                    v = {a: dict(final[a]) for a in AXES}
                    v[ax]["theta0_cable_deg"] = float(g0)
                    vals.append(v)
                res = pool.map(job, [(w, v, mc, False) for v in vals for w in w_sq])
                nw = len(w_sq)
                costs = [window_cost(res[k * nw:(k + 1) * nw], ls) for k in range(len(grid))]
                kb = int(np.argmin(costs))
                v = {a: dict(final[a]) for a in AXES}
                v[ax]["theta0_cable_deg"] = float(grid[kb])
                v_sq, *_ = lsq(run, ax, w_sq, ["theta0_cable_deg"], v, bounds, 10, False, log=f"theta0 of {short(sq)}",
                               loss_scale=ls)
                js["theta0_per_sequence"].setdefault(ax, {})[sq] = v_sq[ax]["theta0_cable_deg"]
                for w, r in zip(w_sq, run(w_sq, v_sq, False)):
                    js["fits"]["validation_theta0_per_sequence"][w.label] = fit_of(r)
        # neutral angle of the session: estimated on one validation file, used on the other ones
        js["fits"]["validation_theta0_one_file"] = {}
        for ax, by_sq in js["theta0_per_sequence"].items():
            for sq_cal, th_cal in by_sq.items():
                w_te = [w for w in wval if w.axis == ax and w.sequence != sq_cal]
                v = {a: dict(final[a]) for a in AXES}
                v[ax]["theta0_cable_deg"] = th_cal
                for w, r in zip(w_te, run(w_te, v, False)):
                    js["fits"]["validation_theta0_one_file"][f"{short(sq_cal)} -> {w.label}"] = fit_of(r)
        js["windows"] = {"amplitude_deg": amp}
        js["wall_time_s"] = js.get("wall_time_s", 0.0) + time.time() - t_start
        save()

    # --------------------------------------------------------------------------- tables
    thr = float(cc["threshold_position_fit"])
    md = ["# Plant driven by the measured currents, identified on VClosedLoop/V2", "",
          f"Configuration `{cc_path.name}`, model `{cc['base_config']}` (H6 with the plant of vcl_study.py). The measured "
          f"motor currents (low-passed at {cc['input_lowpass_hz']:g} Hz, zero phase) drive the plant without the drives; "
          "the goodness measure is the NRMSE fit of the encoder position of the moving axis (change from the window "
          f"start), windows from {pre:g} s before each move to {after:g} s after the end of its command.", "",
          f"Identification: {len(seq_id)} doublet files of V2 ({len(wid['azimuth'])} azimuth and "
          f"{len(wid['elevation'])} elevation moves). Excluded: "
          + "; ".join(f"`{k}` ({v})" for k, v in (cc.get("excluded_sequences") or {}).items())
          + f". Validation: the {len(seq_val)} step files of VClosedLoop ({len(wval)} moves, another session).", ""]
    md += ["## 0. Model-free analysis of V2", "",
           "Rests: mean motor torque K_t i from 3 s after a move to the next move, regressed on angle and direction of "
           "the last move (T = K theta + c + h d). Cruise at 10 deg/s: T = F d + K theta + c over all the samples; F is "
           "the friction at 10 deg/s, K theta + c the position-dependent torque (torques in nominal K_t units).", ""]
    rows = []
    for ax in AXES:
        rf, cr = js["model_free"][ax]["rest_fit"], js["model_free"][ax]["cruise"]
        rows.append([ax, f"{cr['friction']:.3f}", f"{cr['gradient']:.3f}", f"{cr['offset']:+.3f}", f"{cr['rms']:.3f}",
                     f"{rf['direction']:.3f}", f"{rf['gradient']:.3f}", f"{rf['rms']:.3f}"])
    md += [out.markdown_table(["axis", "cruise: friction F [N m]", "cruise: gradient K [N m/rad]", "cruise: offset c [N m]",
                               "cruise: rms [N m]", "rests: direction term h [N m]", "rests: gradient K [N m/rad]",
                               "rests: rms [N m]"], rows), ""]
    md += ["## 1. Rigid surrogate", "", "Each axis turning as a rigid body (inertia r^T M r, same motor, friction and "
           "extra torques) against the full H6 model, at the start values, on all the identification windows.", ""]
    rows = []
    for ax in AXES:
        s = js["surrogate"][ax]
        rows.append([ax, f"{np.median(s['fit_difference']):.2f}", f"{max(s['fit_difference']):.2f}",
                     f"{100 * max(s['position_difference_rel']):.2f} %", f"{s['wall_s']['rigid']:.0f}", f"{s['wall_s']['full']:.0f}"])
    md += [out.markdown_table(["axis", "fit difference, median [points]", "max [points]", "max position difference / move",
                               "wall time rigid [s]", "wall time full [s]"], rows), ""]
    md += ["## 2. Structures and cross-validation", "",
           "Fit on all of V2 (identification median) and leave-one-file-out cross-validation (median of the fits of the "
           "held-out moves, each fold refitted from the fit on all of V2); rigid surrogate. Selected: best CV median, or "
           f"the structure with fewer free parameters within {par:g} points.", ""]
    rows = []
    for ax in AXES:
        for k, r in js["structures"][ax].items():
            mark = " **(selected)**" if k == js["selection"][ax]["selected"] else ""
            rows.append([ax, k + mark, str(len(r["free"])), f"{r['fit_median']:.1f}", f"{r['cv_median']:.1f}",
                         f"{min(r['cv_fits'].values()):.0f}", f"{r['cv_mean_square']:.3f}"])
    md += [out.markdown_table(["axis", "structure", "free parameters", "identification median [%]", "CV median [%]",
                               "CV worst [%]", "CV mean (1 - fit)^2"], rows), ""]
    md += ["## 3. Selected structure, full H6 model", "",
           "Standard errors from the Jacobian of the position residuals with the effective number of independent samples "
           "(relative for the positive parameters, absolute for theta0).", ""]
    rows = []
    for ax in AXES:
        sel = js["structures"][ax][js["selection"][ax]["selected"]]
        for n in sorted(set(final[ax]) - {"theta0_by_seq"}):
            if n not in sel["free"] and n not in sel["fixed"]:
                continue
            se = js["full"][ax]["se"].get(n)
            src = "fitted" if n in sel["free"] else ("fixed (nominal)" if n not in ("K_cable",) else "fixed (no spring)")
            jk = "-"
            if n in sel["free"]:  # jackknife standard error over the leave-one-file-out fits
                x = np.array([v[n] for v in sel["cv_values"].values()])
                x = x if n in SIGNED else np.log(x)
                e = float(np.sqrt((x.size - 1) / x.size * np.sum((x - x.mean()) ** 2)))
                jk = f"{e:.2g} deg" if n in SIGNED else f"{100 * e:.0f} %"
            rows.append([ax, n, f"{final[ax][n]:.4g}", src,
                         "-" if se is None else (f"{se:.2g} deg" if n in SIGNED else f"{100 * se:.0f} %"), jk,
                         "-" if n not in sel["free"] else f"{sel['values'][n]:.4g}"])
    md += [out.markdown_table(["axis", "parameter", "value", "how", "standard error", "jackknife (CV folds)",
                               "rigid surrogate value"], rows), ""]
    md += ["## 4. Position fit, full model", "", f"NRMSE fit [%], median [min, max]; adequate when the validation median "
           f"reaches {thr:g} %.", ""]
    rows = []
    for lab, d in (("identification (V2)", js["fits"]["identification"]), ("validation (VClosedLoop steps)", js["fits"]["validation"]),
                   ("validation, theta0 of the session from one other file", js["fits"]["validation_theta0_one_file"]),
                   ("validation, theta0 per sequence", js["fits"]["validation_theta0_per_sequence"])):
        for ax in AXES:
            f = np.array([v for k, v in d.items() if f" {ax} " in k])
            if f.size:
                rows.append([lab, ax, str(f.size), f"{np.median(f):.1f} [{f.min():.0f}, {f.max():.0f}]",
                             ("yes" if np.median(f) >= thr else "no") if lab.startswith("validation") else "-"])
    md += [out.markdown_table(["set", "axis", "moves", "fit", "adequate"], rows), "",
           "theta0 of the session: theta0 estimated alone on one validation file (all other parameters from V2) and used "
           "on the moves of the other three files; every file in turn (moves counted once per calibration file). "
           "theta0 per sequence: estimated on the moves it is then scored on (an unknown initial condition).", ""]
    if js.get("theta0_per_sequence"):
        md += [out.markdown_table(["axis", "sequence", "theta0 [deg]"],
                                  [[ax, short(sq), f"{v:+.2f}"] for ax, d in js["theta0_per_sequence"].items()
                                   for sq, v in d.items()]), "",
               "theta0 identified on V2: " + ", ".join(f"{ax} {final[ax]['theta0_cable_deg']:+.2f} deg" for ax in AXES
                                                     if final[ax].get("K_cable", 0.0) > 0) + ".", ""]
    rows = [[short(k), f"{amp[k]:+.2f}", f"{v:.1f}"] for k, v in js["fits"]["validation"].items()]
    md += ["Validation moves:", "", out.markdown_table(["move", "amplitude [deg]", "fit [%]"], rows), ""]
    (res_dir / "summary.md").write_text("\n".join(md))

    # --------------------------------------------------------------------------- figures
    # rests and cruise
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    for j, ax in enumerate(AXES):
        a = axs[0, j]
        R = js["model_free"][ax]["rests"]
        for d, mk, col, lab in ((1, "^", "C3", "after a positive move"), (-1, "v", "C0", "after a negative move")):
            x = [r["angle_deg"] for r in R if r["direction"] == d]
            a.plot(x, [r["torque"] for r in R if r["direction"] == d], mk, color=col, label=f"{lab} (mean over the rest)")
            a.vlines(x, [r["torque_end"] for r in R if r["direction"] == d], [r["torque_start"] for r in R if r["direction"] == d],
                     color=col, lw=0.8, alpha=0.6)
        rf = js["model_free"][ax]["rest_fit"]
        xx = np.linspace(-16, 16, 50)
        for d, col in ((1, "C3"), (-1, "C0")):
            a.plot(xx, rf["gradient"] * np.deg2rad(xx) + rf["offset"] + d * rf["direction"], "--", color=col, lw=0.8,
                   label=f"regression, direction {d:+d}")
        a.set_title(f"{ax}: holding torque at rest (bar: start to end of the rest)")
        a.set_xlabel("encoder angle [deg]")
        a.set_ylabel("motor torque K_t i [N m]")
        a.legend(loc="best")
        a = axs[1, j]
        cr = mf[ax]["_cruise"]
        k = np.arange(0, cr["angle_deg"].size, 20)
        a.plot(cr["angle_deg"][k], cr["torque"][k], ".", color="0.75", ms=1.5, label="cruise samples (every 20th)")
        B = np.array(cr["bins"])
        if B.size:
            a.plot(B[:, 0], B[:, 1], "^-", color="C3", label="mean, positive direction")
            a.plot(B[:, 0], B[:, 2], "v-", color="C0", label="mean, negative direction")
            a.plot(B[:, 0], 0.5 * (B[:, 1] + B[:, 2]), "s-", color="k", label="half-sum: position-dependent torque")
            a.plot(B[:, 0], 0.5 * (B[:, 1] - B[:, 2]), "o-", color="C2", label="half-difference: friction at 10 deg/s")
        a.set_title(f"{ax}: torque in the cruise at 10 deg/s")
        a.set_xlabel("encoder angle [deg]")
        a.set_ylabel("motor torque K_t i [N m]")
        a.set_ylim(-1.2, 1.2)
        a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "model_free")
    # cross-validation
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    for j, ax in enumerate(AXES):
        a = axs[j]
        S = js["structures"][ax]
        ks = list(S)
        y = np.arange(len(ks))
        a.barh(y, [S[k]["cv_median"] for k in ks], color=["C2" if k == js["selection"][ax]["selected"] else "C0" for k in ks],
               alpha=0.8, label="CV median (held-out moves)")
        a.plot([S[k]["fit_median"] for k in ks], y, "kd", label="identification median (all of V2)")
        a.set_yticks(y, [f"{k} ({len(S[k]['free'])} free)" for k in ks], fontsize=6.5)
        a.set_xlabel("position fit NRMSE [%]")
        lo_ = min(min(S[k]["cv_median"] for k in ks), 0.0)
        a.set_xlim(lo_ - 5, 100)
        a.set_title(f"{ax}: structures (green: selected)")
        a.legend(loc="lower left")
    fig.tight_layout()
    out.save(fig, fig_dir, "cross_validation")
    # fit against the amplitude
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for j, ax in enumerate(AXES):
        a = axs[j]
        sel = js["structures"][ax][js["selection"][ax]["selected"]]
        for d, mk, col, lab in ((js["fits"]["identification"], "o", "C0", "identification, V2 (full model)"),
                                (sel["cv_fits"], "x", "C1", "cross-validation, held-out V2 move (rigid)"),
                                (js["fits"]["validation"], "s", "C3", "validation, VClosedLoop steps (full model)")):
            ks = [k for k in d if f" {ax} " in k]
            a.semilogx([abs(amp[k]) for k in ks], np.clip([d[k] for k in ks], -100, 100), mk, color=col, label=lab, mfc="none")
        a.axhline(thr, color="k", lw=0.8, ls=":", label=f"threshold {thr:g} %")
        a.set_xlabel("move amplitude [deg]")
        a.set_ylabel("position fit NRMSE [%] (clipped at -100)")
        a.set_title(f"{ax}: fit of each move, selected structure")
        a.legend(loc="lower right")
    fig.tight_layout()
    out.save(fig, fig_dir, "fit_vs_amplitude")
    # position histories
    for lab, wins, res in (("identification", None, None), ("validation", wval, r_val)):
        sel = []
        for ax in AXES:
            if lab == "identification":
                ws = sorted(zip(wid[ax], r_id[ax]), key=lambda p: abs(amp[p[0].label]))
                pick = [ws[0], ws[len(ws) // 2], ws[-1]]
            else:
                ws = [(w, r) for w, r in zip(wins, res) if w.axis == ax]
                pick = ws[:1] + ws[-1:] if len(ws) > 1 else ws
            sel += pick
        fig, axs = plt.subplots(len(sel), 1, figsize=(7.5, 2.0 * len(sel) + 0.4), squeeze=False)
        for a, (w, r) in zip(axs[:, 0], sel):
            ax = w.axis
            big = abs(amp[w.label]) >= 2.0
            k_, u = (np.rad2deg(1.0), "deg") if big else (MDEG, "mdeg")
            a.plot(r["t"], k_ * r[ax]["enc_meas"], "k", lw=1.6, label="measured")
            a.plot(r["t"], k_ * r[ax]["enc"], "C0", lw=1.1, label=f"model (fit {fit_of(r):.1f} %)")
            a.set_title(f"{lab}: {short(w.label)}, amplitude {amp[w.label]:+.2f} deg", fontsize=7.5)
            a.set_ylabel(f"encoder position change [{u}]")
            a.set_xlabel("time from window start [s]")
            a.legend(loc="best")
        fig.tight_layout()
        out.save(fig, fig_dir, f"position_{lab}")
    print(f"done in {(time.time() - t_start) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
