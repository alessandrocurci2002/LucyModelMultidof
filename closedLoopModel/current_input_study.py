r"""Plant driven by the measured currents, fitted on the encoder position (VClosedLoop sequences).

The measured motor currents drive the model (:class:`lucy_cl.current_driven.CurrentDrivenPlant`); the only
measure of goodness is the NRMSE fit of the encoder position of the moving axis. Starting from iteration I5 of
``vcl_study.py`` (first guess), the parameters of each axis are

* screened: one at a time by +/- ``screening.step``, kept when the cost changes by more than
  ``min_relative_change`` (the others are not seen by the position histories and keep the first guess);
* started from the best point of a small grid (torque scale, break-away ratio);
* fitted by least squares on the identification sequences;

and validated on the other sequences; adequate when the median position fit reaches the threshold.

Usage (from this folder)::

    ../.venv/bin/python -I current_input_study.py --config config_current_input.yaml
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import main as mn  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl import plant_fit as pf  # noqa: E402
from lucy_cl.current_driven import CurrentDrivenPlant  # noqa: E402
from lucy_cl.validation import read_sequence  # noqa: E402

MDEG = out.MDEG
plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 6.5,
                     "savefig.bbox": "tight"})


def fit_nrmse(meas: np.ndarray, model: np.ndarray) -> float:
    from gimbal_mdof.sperimentale import fit_nrmse as f

    ok = np.isfinite(model)
    return f(meas[ok], model[ok]) if ok.sum() >= 10 else float("nan")


# --------------------------------------------------------------------------- windows


class Win:
    def __init__(self, seq: str, axis: str, move: float, end: float, t0: float, t1: float):
        self.sequence, self.axis, self.move, self.end, self.t0, self.t1 = seq, axis, move, end, t0, t1

    @property
    def label(self) -> str:
        return f"{self.sequence} {self.axis} {self.move:.2f} s"


def windows(seq, name: str, pre: float, after: float, min_rest: float = 0.6) -> list[Win]:
    every = [m for ax in AXES for m in seq.moves(ax)]
    out_ = []
    for ax in AXES:
        for a, b in seq.moves(ax):
            if any(s < a - 1e-6 and e > a - min_rest for (s, e) in every) or a - pre < seq.t[0] + 0.1:
                continue
            nxt = min([s for (s, _) in every if s > a + 1e-6] + [seq.t[-1]])
            out_.append(Win(name, ax, float(a), float(b), a - pre, min(b + after, nxt - 0.05)))
    return out_


# --------------------------------------------------------------------------- workers

_W: dict = {}


def _init(cfg: dict, seq_dir: str, amps: float, lowpass: float, ring_ab: dict, base_values: dict) -> None:
    cl, _, _ = mn.build(cfg, HERE)
    if not cfg["plant"].get("balanced_telescope", True):  # gravity: only the offset along the telescope axis
        for n in ("d_y", "d_z"):
            cl.plant.set_value(n, 0.0)
    cd = CurrentDrivenPlant(cl.plant, cl.motors, input_lowpass_hz=lowpass)
    _W.update(cd=cd, seqs={}, dir=Path(seq_dir), amps=amps,
              ring=pf.RingMap(ring_ab["a"], ring_ab["b"]), base=base_values)


def _seq(name: str):
    if name not in _W["seqs"]:
        _W["seqs"][name] = read_sequence(_W["dir"] / f"{name}.csv", _W["amps"])
    return _W["seqs"][name]


SIGNED = {"d_x_mm", "H_hyst", "trig_a", "trig_b", "trig_c", "theta0_cable_deg"}
"""Parameters fitted on their value (not their logarithm)."""


def _apply(values: dict[str, dict[str, float]]) -> None:
    cd = _W["cd"]
    pl = cd.plant
    for ax, v in values.items():
        v = _W["base"][ax] | v
        cd.torque_scale[ax] = float(v["scale"])
        cd.set_friction(ax, sigma0=v["sigma0"], Tc=v["Tc"], Ts=v["Tc"] * v["Ts_ratio"], vs=v["vs"])
        if "H_hyst" in v:  # two-level hysteresis (relay hysteron) on the bearing node
            cd.hysteron[ax] = (float(v["H_hyst"]), float(v["delta_hyst"]))
        if "trig_a" in v:  # torque a sin(theta) + b cos(theta) + c on the bearing node
            cd.position_torque[ax] = (float(v["trig_a"]), float(v["trig_b"]), float(v["trig_c"]))
        if "K_cable" in v:  # elastic torque of the cables -K_c (theta - theta_0)
            cd.cable_spring[ax] = (float(v["K_cable"]), float(np.deg2rad(v["theta0_cable_deg"])))
        for n, x in v.items():
            if isinstance(x, dict):
                continue
            if n == "d_x_mm" and ax == "elevation":  # centre of gravity along the telescope axis (static moment)
                pl.set_value("d_x", x * 1e-3)
            elif n == "k_b":
                pl.set_value(pf.K4, x)
            elif n == "ring":
                for s in (pf.K1, pf.K2):
                    pl.set_value(s, _W["ring"].k(s, x))
            elif n.startswith(("c_{", "k_{", "b_{")):
                pl.set_value(n, x)


def job(args) -> dict:
    win, values, max_cpu = args
    values = {ax: dict(v) for ax, v in values.items()}
    for ax, v in values.items():  # neutral angle of the cables of this sequence, when estimated per sequence
        by_seq = v.get("theta0_by_seq") or (_W["base"][ax].get("theta0_by_seq") or {})
        if win.sequence in by_seq:
            v["theta0_cable_deg"] = by_seq[win.sequence]
    _apply(values)
    r = _W["cd"].simulate_window(_seq(win.sequence), win.t0, win.t1, max_cpu_s=max_cpu)
    try:
        return {"t": r["t"], "ok": r["ok"], "axis": win.axis,
                **{ax: {"enc_meas": r[ax]["enc_meas"], "enc": r[ax]["enc"]} for ax in AXES}}
    finally:
        _apply({ax: {} for ax in AXES})  # back to the base values


def residual_of(r: dict) -> np.ndarray:
    ax = r["axis"]
    y, yh = r[ax]["enc_meas"], r[ax]["enc"]
    d = max(np.linalg.norm(y - np.mean(y)), 1e-12)
    res = (yh - y) / d
    return np.where(np.isfinite(res), res, 1.0 / np.sqrt(max(res.size, 1)))


def cost_of(results: list[dict]) -> float:
    return float(sum(np.sum(residual_of(r) ** 2) for r in results))


# --------------------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config_current_input.yaml")
    args = ap.parse_args()
    t_start = time.time()
    cc_path = Path(args.config).resolve()
    cc = yaml.safe_load(cc_path.read_text())
    cfg = yaml.safe_load((HERE / cc["base_config"]).read_text())
    res_dir = HERE / cc["results_dir"]
    fig_dir = res_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    seq_dir = (HERE / cc["sequences_dir"]).resolve()
    amps = float(cc.get("amps_per_A", 2.0))
    mc = float(cc.get("max_cpu_per_window_s", 60))
    fl = lambda v: [float(x) for x in v]  # noqa: E731  (PyYAML: 1.0e4 without sign is a string)
    bounds = {k: tuple(fl(v)) for k, v in cc["bounds"].items()}

    # first guess: I5 and the fitted plant (or a previous run of this study), gravity when asked
    if cc.get("gravity"):
        cfg["plant"]["balanced_telescope"] = False
    cl, _, _ = mn.build(cfg, HERE)
    i5 = json.loads((HERE / cc["first_guess"]).read_text())["iterations"][-1]
    rm = pf.RingMap.from_config(cl.plant.config)
    base = {}
    for ax in AXES:
        el = [e for e in cl.plant.friction if e.axis == ax][0]
        b = {"scale": float(cc["start"]["scale"]), "sigma0": el.sigma0, "Tc": float(i5["values"][ax]["Tc"]),
             "Ts_ratio": float(cc["start"]["Ts_ratio"]), "vs": float(cc["start"]["vs"])}
        for n in cc["parameters"][ax]:
            if n == "k_b":
                b[n] = cl.plant.value(pf.K4)
            elif n == "ring":
                b[n] = float(rm.b[pf.K1] / (1.0 / cl.plant.value(pf.K1) - rm.a[pf.K1]))
            elif n.startswith(("c_{", "b_{")):
                b[n] = cl.plant.value(n)
        if cc.get("gravity"):
            b["d_x_mm"] = 0.0
        for n, v in (cc.get("new_parameters") or {}).get(ax, {}).items():
            b[n] = float(v)
        base[ax] = b
    if cc.get("first_guess_current_input"):
        prev = json.loads((HERE / cc["first_guess_current_input"]).read_text())["final"]
        for ax in AXES:
            base[ax].update({k: float(v) for k, v in prev[ax].items()})
    for ax, d in (cc.get("start_override") or {}).items():
        base[ax].update({k: float(v) for k, v in d.items()})
    for ax, d in (cc.get("fixed") or {}).items():  # values imposed and not fitted
        base[ax].update({k: float(v) for k, v in d.items()})
    if cc.get("gravity"):
        for ax in AXES:
            cl.plant.set_value("d_x", base[ax].get("d_x_mm", 0.0) * 1e-3)
        cl.plant.set_value("d_z", 0.0)
        cl.plant.set_value("d_y", 0.0)
    first = {ax: dict(b) for ax, b in base.items()}
    for ax in AXES:
        for n in cc["parameters"][ax]:
            if n.startswith("c_{"):
                bounds[n] = (first[ax][n] * float(cc["bounds"]["damper_factor"][0]),
                             first[ax][n] * float(cc["bounds"]["damper_factor"][1]))
    seqs = {n: read_sequence(seq_dir / f"{n}.csv", amps) for n in cc["identification_sequences"] + cc["validation_sequences"]}
    wid = [w for n in cc["identification_sequences"] for w in windows(seqs[n], n, cc["window"]["pre_s"], cc["window"]["after_move_s"])]
    wval = [w for n in cc["validation_sequences"] for w in windows(seqs[n], n, cc["window"]["pre_s"], cc["window"]["after_move_s"])]
    md = ["# Plant driven by the measured currents: fit on the encoder position", "",
          f"Configuration `{cc_path.name}`; model `{cc['base_config']}` with the first guess of iteration I5 "
          f"(`{cc['first_guess']}`). The measured motor currents (low-passed at {cc['input_lowpass_hz']:g} Hz, zero "
          "phase) drive the model without the drives; the only goodness measure is the NRMSE fit of the encoder "
          "position of the moving axis (change from the start of the window). Windows from "
          f"{cc['window']['pre_s']:g} s before each move to {cc['window']['after_move_s']:g} s after the end of its "
          "command. Identification: " + ", ".join(cc["identification_sequences"]) + f" ({len(wid)} moves); validation: "
          + ", ".join(cc["validation_sequences"]) + f" ({len(wval)} moves).", "",
          "Parameters: `scale` torque per measured current over K_t (scale of the current channel, tolerance of K_t, "
          "inertia of the CAD); bearing friction `sigma0` [N m/rad], `Tc` sliding level [N m], `Ts_ratio` break-away "
          "over sliding level, `vs` Stribeck velocity [rad/s] (Dahl element with the LuGre Stribeck curve); `k_b`, "
          "`ring`, dampers `c_*` as in vcl_study.py.", ""]
    js: dict = {"config": cc_path.name, "first_guess": first, "windows": {"identification": [w.label for w in wid],
                                                                        "validation": [w.label for w in wval]}}
    final = {ax: dict(b) for ax, b in first.items()}
    with mp.get_context("fork").Pool(int(cc.get("workers", 4)), initializer=_init,
                                     initargs=(cfg, str(seq_dir), amps, float(cc["input_lowpass_hz"]),
                                               {"a": rm.a, "b": rm.b}, base)) as pool:

        def run(wins: list[Win], vals: dict) -> list[dict]:
            return pool.map(job, [(w, vals, mc) for w in wins])

        def fits(results: list[dict]) -> list[float]:
            return [fit_nrmse(r[r["axis"]]["enc_meas"], r[r["axis"]]["enc"]) for r in results]

        # ---- first guess on both sets
        r_id0, r_val0 = run(wid, first), run(wval, first)
        js["first_guess_fits"] = {"identification": fits(r_id0), "validation": fits(r_val0)}
        print("  first guess: identification", np.round(fits(r_id0), 1), "validation", np.round(fits(r_val0), 1), flush=True)

        js["axes"] = {}
        rows_s, rows_p = [], []
        for ax in cc.get("axes", AXES):
            w_ax = [w for w in wid if w.axis == ax]
            # ---- screening
            c0 = cost_of(run(w_ax, first))
            keep, sens = [], {}
            step = float(cc["screening"]["step"])
            for n in cc["parameters"][ax]:
                dc = []
                for fct in (1 + step, 1 - step):
                    v = {a_: dict(first[a_]) for a_ in AXES}
                    sst = cc["screening"].get("signed_step", 5.0)
                    sst = float(sst.get(n, 0.1) if isinstance(sst, dict) else sst)
                    x = first[ax][n] * fct if n not in SIGNED else first[ax][n] + (fct - 1) / step * sst
                    if n == "Ts_ratio":
                        x = max(1.0, x)
                    v[ax][n] = x
                    dc.append(cost_of(run(w_ax, v)) - c0)
                rel = max(abs(d) for d in dc) / max(c0, 1e-12)
                sens[n] = rel
                if rel > float(cc["screening"]["min_relative_change"]):
                    keep.append(n)
                rows_s.append([ax, n, f"{first[ax][n]:.4g}", f"{100 * rel:.1f} %", "yes" if n in keep else "no"])
            keep += [n for n in cc.get("always_fit", {}).get(ax, []) if n not in keep]
            print(f"  {ax}: cost {c0:.3f}, kept {keep}", flush=True)
            # ---- grid (parameters of the axis listed in the grid)
            import itertools

            gkeys = [k for k in cc["grid"] if k in first[ax]]
            gpts = list(itertools.product(*[fl(cc["grid"][k]) for k in gkeys]))
            gcost = []
            for pt in gpts:
                v = {a_: dict(first[a_]) for a_ in AXES}
                v[ax].update(dict(zip(gkeys, pt)))
                gcost.append(cost_of(run(w_ax, v)))
            kb = int(np.argmin(gcost))
            start = dict(first[ax])
            start.update(dict(zip(gkeys, gpts[kb])))
            print(f"  {ax}: grid of {len(gpts)} points, best {dict(zip(gkeys, gpts[kb]))}, cost {gcost[kb]:.3f}", flush=True)
            names = [n for n in keep] or ["scale"]
            for n in ("scale", "Ts_ratio"):
                if n not in names:
                    names.append(n)
            per_seq = bool(cc.get("theta0_per_sequence")) and "theta0_cable_deg" in names
            if per_seq:  # neutral angle of the cables: one unknown per identification sequence
                names.remove("theta0_cable_deg")
                for sq in cc["identification_sequences"]:
                    n_ = f"theta0@{sq}"
                    names.append(n_)
                    start[n_] = start["theta0_cable_deg"]
                    bounds[n_] = bounds["theta0_cable_deg"]
                    SIGNED.add(n_)
            # ---- least squares
            from scipy.optimize import least_squares

            tr = lambda n, x: x if n in SIGNED else np.log(x)  # noqa: E731
            itr = lambda n, x: x if n in SIGNED else np.exp(x)  # noqa: E731
            lo = np.array([tr(n, bounds[n][0]) for n in names])
            hi = np.array([tr(n, bounds[n][1]) for n in names])
            x0 = np.clip(np.array([tr(n, start[n]) for n in names]), lo + 1e-9, hi - 1e-9)

            def to_values(x: np.ndarray) -> dict:
                v = {a_: dict(first[a_]) for a_ in AXES}
                for n, xi in zip(names, x):
                    if n.startswith("theta0@"):
                        v[ax].setdefault("theta0_by_seq", {})[n.split("@", 1)[1]] = float(itr(n, xi))
                    else:
                        v[ax][n] = float(itr(n, xi))
                return v

            def resid(x: np.ndarray) -> np.ndarray:
                v = to_values(x)
                rr = np.concatenate([residual_of(r) for r in run(w_ax, v)])
                pri = (cc.get("prior") or {}).get(ax, {})
                if pri:  # weak regularisation: weight * (x - x0) / scale
                    rr = np.concatenate([rr, [float(w_) * (v[ax][n] - float(x0_)) / float(sc_)
                                              for n, (x0_, sc_, w_) in pri.items() if n in v[ax]]])
                print(f"    {ax}: " + ", ".join(f"{n} {itr(n, e):.4g}" for n, e in zip(names, x))
                      + f" -> cost {np.sum(rr ** 2):.4f}", flush=True)
                return rr

            sol = least_squares(resid, x0, bounds=(lo, hi), diff_step=0.03, max_nfev=int(cc.get("max_nfev", 40)),
                                xtol=1e-3, ftol=1e-4, x_scale=[10.0 if n in SIGNED else 1.0 for n in names])
            final[ax] = to_values(sol.x)[ax]
            if per_seq:
                # validation sequences: their neutral angle estimated alone, all other parameters fixed
                from scipy.optimize import least_squares as lsq

                for sq in cc["validation_sequences"]:
                    w_sq = [w for w in wval if w.axis == ax and w.sequence == sq]
                    if not w_sq:
                        continue

                    def r_sq(xx: np.ndarray, w_sq=w_sq, sq=sq) -> np.ndarray:
                        v = {a_: dict(final[a_]) for a_ in AXES}
                        v[ax]["theta0_by_seq"] = dict(v[ax].get("theta0_by_seq", {})) | {sq: float(xx[0])}
                        return np.concatenate([residual_of(r) for r in run(w_sq, v)])

                    s1 = lsq(r_sq, [final[ax].get("theta0_cable_deg", 0.0)], bounds=tuple(bounds["theta0_cable_deg"]),
                             diff_step=0.03, x_scale=[1.0], max_nfev=20)
                    final[ax].setdefault("theta0_by_seq", {})[sq] = float(s1.x[0])
                    print(f"  {ax}: theta0 of {sq} (validation) {s1.x[0]:.3f} deg", flush=True)
            js["axes"][ax] = {"screening": sens, "kept": keep, "grid": {"points": gpts, "costs": gcost},
                              "start": start, "names": names, "values": final[ax],
                              "cost": [float(c0), float(np.sum(sol.fun ** 2))]}
            for n in names:
                if n.startswith("theta0@"):
                    sq = n.split("@", 1)[1]
                    rows_p.append([ax, f"theta0_cable_deg ({sq})", f"{first[ax].get('theta0_cable_deg', 0.0):.4g}",
                                   f"{start[n]:.4g}", f"{final[ax]['theta0_by_seq'][sq]:.4g}"])
                else:
                    rows_p.append([ax, n, f"{first[ax][n]:.4g}", f"{start[n]:.4g}", f"{final[ax][n]:.4g}"])
            for sq, v_ in (final[ax].get("theta0_by_seq") or {}).items():
                if sq in cc["validation_sequences"]:
                    rows_p.append([ax, f"theta0_cable_deg ({sq}, validation, estimated alone)", "-", "-", f"{v_:.4g}"])

        # ---- identification and validation with the fitted parameters
        r_id, r_val = run(wid, final), run(wval, final)
    js["final"] = final
    thr = float(cc["threshold_position_fit"])
    rows_f = []
    js["fits"] = {}
    for lab, wins, r0, r1 in (("identification", wid, r_id0, r_id), ("validation", wval, r_val0, r_val)):
        for ax in AXES:
            k = [i for i, w in enumerate(wins) if w.axis == ax]
            f0 = np.array([fits_ for fits_ in [fit_nrmse(r0[i][ax]["enc_meas"], r0[i][ax]["enc"]) for i in k]])
            f1 = np.array([fit_nrmse(r1[i][ax]["enc_meas"], r1[i][ax]["enc"]) for i in k])
            js["fits"][f"{lab}/{ax}"] = {"first_guess": f0.tolist(), "fitted": f1.tolist(),
                                         "moves": [wins[i].label for i in k]}
            rows_f.append([lab, ax, str(len(k)), f"{np.median(f0):.1f} [{f0.min():.0f}, {f0.max():.0f}]",
                           f"{np.median(f1):.1f} [{f1.min():.0f}, {f1.max():.0f}]",
                           ("yes" if np.median(f1) >= thr else "no") if lab == "validation" else "-"])
    # ---- tables
    md += ["## Screening of the parameters", "",
           f"Relative change of the cost (sum of (1 - fit)^2 of the position over the identification moves of the "
           f"axis) for a change of +/-{100 * float(cc['screening']['step']):.0f} % of each parameter from the first guess; "
           f"fitted when above {100 * float(cc['screening']['min_relative_change']):.0f} %.", "",
           out.markdown_table(["axis", "parameter", "first guess", "cost change", "fitted"], rows_s), "",
           "## Fitted parameters", "", out.markdown_table(["axis", "parameter", "first guess (I5)", "start (grid)",
                                                            "fitted"], rows_p), "",
           "## Position fit", "", f"NRMSE fit [%] of the encoder position of the moving axis, median [min, max]; "
           f"adequate when the validation median reaches {thr:g} %.", "",
           out.markdown_table(["set", "axis", "moves", "first guess (I5)", "fitted", "adequate"], rows_f), ""]
    # ---- figures: two identification and two validation moves per axis
    for lab, wins, r0, r1 in (("identification", wid, r_id0, r_id), ("validation", wval, r_val0, r_val)):
        sel = []
        for ax in AXES:
            k = [i for i, w in enumerate(wins) if w.axis == ax][:2]
            sel += [(wins[i], r0[i], r1[i]) for i in k]
        fig, axs = plt.subplots(len(sel), 1, figsize=(7.5, 2.1 * len(sel) + 0.4), squeeze=False)
        for a_, (w, a0, a1) in zip(axs[:, 0], sel):
            ax = w.axis
            a_.plot(a1["t"], MDEG * a1[ax]["enc_meas"], "k", lw=1.6, label="measured")
            a_.plot(a0["t"], MDEG * a0[ax]["enc"], "C3--", lw=1.0,
                    label=f"model, first guess I5 (fit {fit_nrmse(a0[ax]['enc_meas'], a0[ax]['enc']):.0f} %)")
            a_.plot(a1["t"], MDEG * a1[ax]["enc"], "C0", lw=1.1,
                    label=f"model, fitted (fit {fit_nrmse(a1[ax]['enc_meas'], a1[ax]['enc']):.0f} %)")
            a_.set_title(f"{lab}: {w.sequence}, {ax} move at {w.move:.2f} s", fontsize=7.5)
            a_.set_ylabel("encoder position change [mdeg]")
            a_.set_xlabel("time from window start [s]")
            a_.legend(loc="best")
        fig.tight_layout()
        out.save(fig, fig_dir, f"position_{lab}")
    (res_dir / "summary.md").write_text("\n".join(md))
    (res_dir / "summary.json").write_text(json.dumps(js, indent=2, default=lambda o: float(o) if np.isscalar(o) else str(o)))
    print(f"results in {res_dir} ({time.time() - t_start:.0f} s)")


if __name__ == "__main__":
    main()
