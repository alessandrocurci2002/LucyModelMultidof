r"""Closed-loop model of LUCY: mechanical plant + motors + drives.

Builds the three components of :mod:`lucy_cl` from a configuration file and runs, in order:

1. motor: consistency of the datasheet constants (:math:`K_t = \sqrt3 K_b`);
2. drive model against the WorkBench exports (``Controller`` and ``Feedback Biquad`` traces);
3. plant + motor + current loop + drive against the measured ``Plant``, ``Open Loop`` and ``Closed Loop``
   traces, with the tuning active during the test;
4. loops with the drive parameter files (.akd): stability margins, bandwidths, closed-loop eigenvalues;
5. sensitivity of the margins to the damping of the plant (:math:`\zeta_1`, placeholder Q8 of gimbal_mdof);
6. check of the time-domain model: its linearisation against the frequency-domain model;
7. nonlinear time simulations of the scenarios.

Results in ``<results_dir>/``: ``figures/``, ``summary.md``, ``summary.json``.

Usage (from this folder)::

    ../.venv/bin/python main.py --config config.yaml
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lucy_cl import AXES, GIMBAL_NAME  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl.closed_loop import (ClosedLoop, Command, Scenario, bandwidth, loop_frf, margins,  # noqa: E402
                                 state_space_frf)
from lucy_cl.data_files import WorkBenchExport, read_workbench  # noqa: E402
from lucy_cl.driver import Drive, DriveParams, controller_workbench  # noqa: E402
from lucy_cl.motor import Motor, MotorParams  # noqa: E402
from lucy_cl.plant import GimbalPlant  # noqa: E402

TRACES = ("Plant", "Open Loop", "Closed Loop")


# --------------------------------------------------------------------------- construction


def build(cfg: dict, base: Path) -> tuple[ClosedLoop, dict[str, Drive], dict[str, WorkBenchExport]]:
    """Plant, motors and drives from the configuration; returns the closed loop (drives with the .akd
    parameters), the drives with the tuning of the WorkBench tests and the WorkBench exports."""
    pc = cfg["plant"]
    plant = GimbalPlant((base / pc["config"]).resolve(), theta0=float(np.deg2rad(pc.get("theta0_deg", 0.0))),
                        balanced_telescope=bool(pc.get("balanced_telescope", True)),
                        zeta_first_flexible=pc.get("zeta_first_flexible"),
                        yaml_overrides=dict(pc.get("yaml_overrides") or {}),
                        parameter_overrides=dict(pc.get("parameter_overrides") or {}),
                        ground_override=dict(pc.get("ground_override") or {}))
    mc = cfg["motor"]
    mp = MotorParams.from_yaml(base / mc["file"], mc.get("temperature", "cold"))
    motors = {ax: Motor(mp, cogging=bool(mc.get("cogging", False))) for ax in AXES}
    dc = cfg["drives"]
    obs = {"from_parameters": None, "on": True, "off": False}[dc.get("observer", "from_parameters")]
    kw = dict(feedforward=bool(dc.get("feedforward", False)), observer=obs, delay_s=float(dc.get("delay_s", 0.0)),
              pl_kp_unit=str(dc.get("position_gain_unit", "hz")))
    drives, tests, wb = {}, {}, {}
    for ax in AXES:
        akd, meas = base / dc[ax]["parameters"], base / dc[ax]["measurement"]
        drives[ax] = Drive(DriveParams.from_files(ax, akd, **kw))
        tests[ax] = Drive(DriveParams.from_files(ax, akd, meas, **kw))
        wb[ax] = read_workbench(meas)
    fr = pc.get("bearing_friction", {}) or {}
    if fr.get("model", "dahl") == "dahl":
        for ax in AXES:
            v = plant.ground(ax)
            if v is not None:
                # amplitude of the test on which k and eta of the ground constraint were fitted (by default
                # the measurement of the drive section)
                src = (fr.get("amplitude_measurement") or {}).get(ax)
                w_amp = read_workbench(base / src) if src else wb[ax]
                A = w_amp.amplitude_rad(float(fr.get("reference_frequency_Hz", np.sqrt(5.0))))
                plant.add_dahl(ax, v[0], v[1], A, float(fr.get("velocity_smoothing_rad_s", 1e-6)))
    return ClosedLoop(plant, motors, drives), tests, wb


def scenarios(cfg: dict) -> list[Scenario]:
    res = []
    for s in cfg.get("scenarios", []):
        cmds = {}
        for ax, c in s["commands"].items():
            val = np.deg2rad(c["value_deg"] if c["kind"] == "step" else c["value_deg_s"])
            cmds[ax] = Command(c["kind"], float(val), float(c.get("t0_s", 0.0)))
        res.append(Scenario(s["name"], cmds, float(s["t_end_s"]), s.get("description", ""),
                            bool(s.get("compare_linear", False))))
    return res


# --------------------------------------------------------------------------- metrics


def step_metrics(t: np.ndarray, y: np.ndarray, target: float) -> dict:
    yn = y / target
    m = {"final_%": 100 * yn[-1], "overshoot_%": max(0.0, 100 * (yn.max() - 1.0))}
    k10, k90 = np.argmax(yn >= 0.1), np.argmax(yn >= 0.9)
    m["rise_10_90_s"] = float(t[k90] - t[k10]) if yn.max() >= 0.9 else float("nan")
    outside = np.nonzero(np.abs(yn - 1.0) > 0.02)[0]
    m["settling_2%_s"] = float("nan") if outside.size and outside[-1] == t.size - 1 else \
        float(t[outside[-1] + 1]) if outside.size else 0.0
    return m


def linear_response(A: np.ndarray, B: np.ndarray, C: np.ndarray, sc: Scenario, t: np.ndarray) -> dict[str, np.ndarray]:
    u = np.zeros((t.size, 2))
    for j, ax in enumerate(AXES):
        c = sc.command(ax)
        u[:, j] = [c(tk)[0] for tk in t]
    _, y, _ = signal.lsim((A, B, C, np.zeros((2, 2))), u, t, interp=False)
    return {ax: y[:, j] for j, ax in enumerate(AXES)}


# --------------------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--skip-time", action="store_true", help="skip the nonlinear time simulations")
    args = ap.parse_args()
    cfg_path = Path(args.config).resolve()
    base = cfg_path.parent
    cfg = yaml.safe_load(cfg_path.read_text())
    res_dir = base / cfg["results_dir"]
    fig_dir = res_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    md: list[str] = [f"# Closed-loop model of LUCY -- {cfg.get('label', '')}", "",
                     f"Configuration: `{cfg_path.name}`; plant: `{cfg['plant']['config']}`. "
                     f"Generated by `main.py`.", ""]
    js: dict = {"config": cfg_path.name, "label": cfg.get("label", "")}

    cl, tests, wb = build(cfg, base)
    pl = cl.plant
    an = cfg.get("analysis", {})
    f = np.logspace(np.log10(an.get("f_min_Hz", 0.5)), np.log10(an.get("f_max_Hz", 2000.0)), int(an.get("n_points", 4000)))

    # ---- components ----------------------------------------------------------------------------------
    md += ["## 1. Components", "", "### Plant", "",
           f"{pl.n}-DOF lumped model of `gimbal_mdof`, coordinates `{', '.join(pl.names)}`, "
           f"elevation of the operating point {np.rad2deg(pl.theta0):g} deg. Modifications:", ""]
    md += [f"- {n}" for n in pl.notes] + [""]
    M, C, K = pl.linear_matrices()
    import scipy.linalg as sla

    w2 = sla.eigh((K + K.T) / 2, M, eigvals_only=True)
    fn = np.sqrt(np.clip(w2, 0, None)) / (2 * np.pi)
    md += ["Natural frequencies without the ground constraint (rotor inertia included), Hz: "
           + ", ".join(out.fmt(x, 4) for x in fn) + "", ""]
    js["natural_frequencies_Hz"] = fn.tolist()
    mp = cl.motors["azimuth"].p
    md += ["### Motor", "",
           out.markdown_table(["quantity", "value"], [
               ["model", f"{mp.name}, winding temperature {mp.temperature}"],
               ["Kt [N m/A rms]", out.fmt(mp.Kt)], ["Kb [V rms/(rad/s)], line-line", out.fmt(mp.Kb)],
               ["Kt / (sqrt3 Kb) - 1", f"{100 * mp.consistency():+.2f} %"],
               ["R, L line-line", f"{mp.R_ll:g} ohm, {1e3 * mp.L_ll:g} mH (electrical pole {mp.R_ll / mp.L_ll / (2 * np.pi):.0f} Hz)"],
               ["rotor inertia [kg m^2]", out.fmt(mp.J_rotor)],
               ["viscous damping [N m s/rad]", out.fmt(mp.b_viscous)],
               ["Ic / Ip [A rms]", f"{mp.I_cont:g} / {mp.I_peak:g}"]]), ""]
    js["motor"] = {"Kt": mp.Kt, "Kb": mp.Kb, "consistency": mp.consistency()}
    md += ["### Drives", "", "Gains of the parameter files (.akd) and, in brackets when different, of the "
           "WorkBench test used for the validation.", ""]
    rows = []
    for ax in AXES:
        p, q = cl.drives[ax].p, tests[ax].p

        def two(a: float, b: float, nd: int = 4) -> str:
            return out.fmt(a, nd) if np.isclose(a, b) else f"{out.fmt(a, nd)} ({out.fmt(b, nd)})"

        filt = "; ".join(f"AR{k + 1}: {b.describe()}" for k, b in enumerate(p.forward + p.feedback) if b.active) or "none"
        filt_t = "; ".join(f"AR{k + 1}: {b.describe()}" for k, b in enumerate(q.forward + q.feedback) if b.active) or "none"
        rows.append([ax, two(p.Kpp / (2 * np.pi), q.Kpp / (2 * np.pi)), two(p.Kip_hz, q.Kip_hz), two(p.Kpv, q.Kpv),
                     two(p.Kiv_hz, q.Kiv_hz), two(p.Kpi, q.Kpi), f"{p.current_bandwidth_hz:.0f}",
                     filt + ("" if filt == filt_t else f" (test: {filt_t})"),
                     f"{'on, ' + format(p.observer_bw_hz, 'g') + ' Hz' if p.observer else 'off'}",
                     f"{p.i_limit:g} / {p.vel_int_limit:g}", f"{p.v_bus:g}" if p.v_bus else "-"])
    md += [out.markdown_table(["axis", "PL.KP [Hz]", "PL.KI [Hz]", "VL.KP [A/(rad/s)]", "VL.KI [Hz]",
                               "IL.KP [V/A]", "current loop [Hz]", "filters", "observer",
                               "I limit / vel. int. limit [A]", "bus [V]"], rows), "",
           ("Feedforward off (`VL.KVFF` of the files ignored)" if not cfg["drives"].get("feedforward", False)
            else "Velocity feedforward `VL.KVFF` of the files active")
           + "; integrators always active (`VL.KIMODE`); `PL.KI` = 0 on both axes. Position gain: "
           + ("2 pi PL.KP (PL.KP read in Hz, as in the WorkBench Controller trace)"
              if cfg["drives"].get("position_gain_unit", "hz") == "hz" else "PL.KP in 1/s") + ".", ""]

    # ---- friction ------------------------------------------------------------------------------------
    if pl.friction:
        md += ["### Bearing friction (Dahl)", "", out.markdown_table(
            ["axis", "node", "k [N m/rad]", "eta", "test amplitude [mrad]", "sigma0 [N m/rad]", "Tc [N m]"],
            [[el.axis, pl.names[el.index], f"{el.k_eq:g}", f"{el.eta_eq:g}", f"{1e3 * el.amplitude:.3f}",
              f"{el.sigma0:.0f}", f"{el.Tc:.3f}"] for el in pl.friction]), ""]
        out.save(out.dahl_loops(pl.friction), fig_dir, "dahl_friction")
        js["friction"] = [{"axis": el.axis, "sigma0": el.sigma0, "Tc": el.Tc, "A": el.amplitude} for el in pl.friction]

    # ---- 2. drive validation -------------------------------------------------------------------------
    md += ["## 2. Drive model against the WorkBench traces", "",
           "Tuning of the test (WorkBench export), traces on the measurement frequencies.", ""]
    rows = []
    js["drive_validation"] = {}
    for ax in AXES:
        w, d = wb[ax], tests[ax]
        s = 2j * np.pi * w.f
        e_c = controller_workbench(d, s) / w.response("Controller")
        e_f = d.F_ar_fb(s) / w.response("Feedback Biquad")
        v = {"controller_dB": float(np.abs(20 * np.log10(np.abs(e_c))).max()),
             "controller_deg": float(np.abs(np.angle(e_c, deg=True)).max()),
             "feedback_dB": float(np.abs(20 * np.log10(np.abs(e_f))).max()),
             "feedback_deg": float(np.abs(np.angle(e_f, deg=True)).max())}
        js["drive_validation"][ax] = v
        rows.append([ax, w.path.name, f"{w.f[0]:g}-{w.f[-1]:.0f}", f"{v['controller_dB']:.3g} dB / {v['controller_deg']:.2g} deg",
                     f"{v['feedback_dB']:.3g} dB / {v['feedback_deg']:.2g} deg"])
        fg = out.bode([("WorkBench Controller", w.f, w.response("Controller"), dict(color="k", marker=".", ls="none", ms=3)),
                       ("model Controller", f, controller_workbench(d, 2j * np.pi * f), dict(color="C0")),
                       ("WorkBench Feedback Biquad", w.f, w.response("Feedback Biquad"), dict(color="0.5", marker=".", ls="none", ms=3)),
                       ("model AR3 AR4", f, d.F_ar_fb(2j * np.pi * f), dict(color="C1", ls="--"))],
                      f"{ax}: drive model vs WorkBench ({w.path.name})", "C")
        out.save(fg, fig_dir, f"drive_validation_{ax}")
    md += [out.markdown_table(["axis", "export", "band [Hz]", "Controller: max error", "Feedback Biquad: max error"], rows), ""]

    # ---- 3. loops against the measurements -----------------------------------------------------------
    from gimbal_mdof.sperimentale import fit_bande

    md += ["## 3. Loops against the measurements (tuning of the test)", "",
           "Model: plant with the H6 ground constraint k(1 + j eta) (frequency domain), motor electrical "
           "dynamics and current loop. Fit NRMSE (%, complex response, System Identification Toolbox "
           "definition) per band; margins from the WorkBench open loop (Controller x Plant x Feedback Biquad).", ""]
    rows_fit, rows_m = [], []
    js["measurement"] = {}
    for ax in AXES:
        w = wb[ax]
        r_meas = loop_frf(cl, ax, w.f, drive=tests[ax])
        r_test = loop_frf(cl, ax, f, drive=tests[ax])
        r_akd = loop_frf(cl, ax, f)
        jm = {}
        for tr, H in zip(TRACES, (r_meas.plant, r_meas.open_loop, r_meas.closed_loop)):
            fb = fit_bande(w.f, w.response(tr), H)
            jm[tr] = fb
            rows_fit.append([ax, tr] + [out.fmt(fb[k], 3) for k in ("all", "1-5 Hz", "5-100 Hz", ">100 Hz")])
        mm, mo = margins(w.f, w.response("Open Loop")), margins(f, r_test.open_loop)
        jm["margins_measured"], jm["margins_model"] = _clean(mm), _clean(mo)
        js["measurement"][ax] = jm
        for lab, m in (("measured", mm), ("model", mo)):
            rows_m.append([ax, lab, out.fmt(m["f_c"], 3), f"{m['pm_deg']:.1f} @ {out.fmt(m['f_pm'], 3)} Hz",
                           f"{m['gm_db']:.1f} @ {out.fmt(m['f_gm'], 3)} Hz"])
        meas = {tr: w.response(tr) for tr in TRACES}
        same = np.allclose([tests[ax].p.Kpv, tests[ax].p.Kiv_hz], [cl.drives[ax].p.Kpv, cl.drives[ax].p.Kiv_hz]) and \
            tests[ax].n_active_filters() == cl.drives[ax].n_active_filters()
        fg = out.loop_comparison(ax, w.f, meas, r_test, None if same else r_akd, f"measured ({w.path.name})")
        out.save(fg, fig_dir, f"loop_vs_measurement_{ax}")
    md += [out.markdown_table(["axis", "trace", "all", "1-5 Hz", "5-100 Hz", ">100 Hz"], rows_fit), "",
           out.markdown_table(["axis", "source", "crossover [Hz]", "phase margin [deg]", "gain margin [dB]"], rows_m), ""]

    # ---- 4. loops with the .akd parameters ------------------------------------------------------------
    A, B, Cy = cl.linearise()
    ev = np.linalg.eigvals(A)
    md += ["## 4. Loops with the drive parameter files (.akd)", "",
           "Exact loop at the current command (position fed back without the velocity feedback filters, "
           "observer included); plant with the H6 ground constraint.", ""]
    rows = []
    js["akd"] = {}
    for ax in AXES:
        r = loop_frf(cl, ax, f)
        m = margins(f, r.loop)
        bw = bandwidth(f, r.position)
        js["akd"][ax] = {"margins": _clean(m), "position_bandwidth_Hz": bw}
        rows.append([ax, out.fmt(m["f_c"], 3), f"{m['pm_deg']:.1f} @ {out.fmt(m['f_pm'], 3)} Hz",
                     f"{m['gm_db']:.1f} @ {out.fmt(m['f_gm'], 3)} Hz", out.fmt(bw, 3)])
        fg = out.bode([("loop L (exact)", f, r.loop, dict(color="C0")),
                       ("position p_fb/p_cmd", f, r.position, dict(color="C2")),
                       ("plant v_fb/i_cmd x 100", f, 100 * r.plant, dict(color="0.6", lw=0.8))],
                      f"{ax}: loops with the .akd parameters", "H", zero_db=True)
        out.save(fg, fig_dir, f"loop_akd_{ax}")
    unstable = ev[ev.real > 1e-6]
    md += [out.markdown_table(["axis", "crossover [Hz]", "phase margin [deg]", "gain margin [dB]",
                               "position bandwidth -3 dB [Hz]"], rows), "",
           f"Closed-loop eigenvalues of the time-domain model linearised at rest (Dahl friction -> springs "
           f"sigma0): largest real part {ev.real.max():.3g} 1/s"
           + (f"; unstable: " + ", ".join(f"{z.real:.1f} {'+' if z.imag >= 0 else '-'} j{abs(z.imag):.0f} "
                                         f"({abs(z.imag) / (2 * np.pi):.0f} Hz)" for z in unstable if z.imag >= 0)
              if unstable.size else " (stable)") + ".", ""]
    js["eigenvalues_max_real"] = float(ev.real.max())
    js["unstable_eigenvalues"] = [[float(z.real), float(z.imag)] for z in unstable]

    # ---- 5. damping sensitivity ----------------------------------------------------------------------
    zs = [float(z) for z in an.get("zeta_sweep", [])]
    if zs:
        z_ref = pl.zeta
        rows, loops = [], []
        js["zeta_sweep"] = []
        for z in zs:
            pl.set_zeta(z)
            Az, _, _ = cl.linearise()
            row = [f"{z:g}"]
            entry = {"zeta": z, "eig_max_real": float(np.linalg.eigvals(Az).real.max())}
            for ax in AXES:
                for lab, d in (("akd", None), ("test", tests[ax])):
                    m = margins(f, loop_frf(cl, ax, f, drive=d).loop)
                    entry[f"{ax}_{lab}"] = _clean(m)
                    row.append(f"{m['pm_deg']:.0f} deg / {m['gm_db']:.1f} dB @ {out.fmt(m['f_gm'], 3)} Hz")
            loops.append((z, loop_frf(cl, "elevation", f).loop))
            row.append(f"{entry['eig_max_real']:.3g}")
            rows.append(row)
            js["zeta_sweep"].append(entry)
        pl.set_zeta(z_ref)
        md += ["## 5. Sensitivity to the damping of the plant", "",
               f"Stiffness-proportional damping, zeta_1 on the first flexible mode of each axis (gimbal_mdof Q8, "
               f"placeholder; this run: {z_ref:g}). Phase margin / smallest gain margin of the exact loop.", "",
               out.markdown_table(["zeta_1", "azimuth .akd", "azimuth test", "elevation .akd", "elevation test",
                                   "max Re(eig) [1/s]"], rows), ""]
        out.save(out.zeta_sweep(f, loops, "elevation", "elevation, exact loop with the .akd parameters"),
                 fig_dir, "zeta_sweep_elevation")

    # ---- 6. linearisation check ----------------------------------------------------------------------
    rows = []
    js["linearisation_check"] = {}
    fl = f[(f >= 0.5) & (f <= 1000.0)]
    for j, ax in enumerate(AXES):
        H_t = state_space_frf(A, B, Cy, fl, j, j)
        H_f = loop_frf(cl, ax, fl, ground="dahl").position
        err = float(np.max(np.abs(H_t - H_f) / np.maximum(np.abs(H_f), 1e-3)))
        js["linearisation_check"][ax] = err
        rows.append([ax, f"{err:.2e}"])
    md += ["## 6. Check of the time-domain model", "",
           "Position response p_fb/p_cmd of the time-domain model linearised by finite differences against the "
           "frequency-domain model with the same linearised friction (0.5-1000 Hz).", "",
           out.markdown_table(["axis", "max relative difference"], rows), ""]

    # ---- 7. time simulations -------------------------------------------------------------------------
    if not args.skip_time:
        md += ["## 7. Nonlinear time simulations", "",
               "Implicit integration (Radau) of the whole closed loop: plant with all nonlinear terms, Dahl "
               "friction, motor electrical dynamics, drives with current limit, integrator clamp and voltage limit "
               "of the bus. Dotted: linearised model (Dahl friction -> spring sigma0), for the scenarios in the "
               "pre-sliding range and when it is stable.", ""]
        rows = []
        js["scenarios"] = {}
        stable = ev.real.max() < 1e-6
        if not stable:
            md += ["The linearised closed loop is unstable (section 4): the simulations show the growth of the "
                   "unstable mode until the current and voltage limits, and stop when the CPU time per scenario "
                   f"exceeds {an.get('max_cpu_per_scenario_s')} s.", ""]
        for sc in scenarios(cfg):
            r = cl.simulate(sc, max_cpu_s=an.get("max_cpu_per_scenario_s"))
            lin = linear_response(A, B, Cy, sc, r.t) if stable and sc.compare_linear else None
            axes = [ax for ax in AXES if ax in sc.commands]
            out.save(out.scenario(r, axes, lin), fig_dir, f"scenario_{sc.name}")
            entry = {"success": r.success, "message": r.message, "n_rhs": r.n_rhs, "cpu_s": r.cpu_s}
            for ax in AXES:
                s = r.signals[ax]
                c = sc.command(ax)
                e = {"max_abs_i_cmd_A": float(np.abs(s["i_cmd"]).max()), "max_abs_v_V": float(np.abs(s["v"]).max()),
                     "max_abs_error_mdeg": float(out.MDEG * np.abs(s["error"]).max())}
                if c.kind == "step":
                    e |= step_metrics(r.t, s["p_enc"], c.value)
                    if lin is not None:
                        e["nonlinear_minus_linear_%"] = float(100 * np.abs(s["p_enc"] - lin[ax]).max() / abs(c.value))
                elif c.kind == "ramp":
                    h = r.t >= r.t[-1] / 2
                    e["rms_error_2nd_half_mdeg"] = float(out.MDEG * np.sqrt(np.mean(s["error"][h] ** 2)))
                    e["rms_load_error_2nd_half_mdeg"] = float(out.MDEG * np.sqrt(np.mean((s["p_cmd"] - s["p_load"])[h] ** 2)))
                entry[ax] = e
                if c.kind != "hold":
                    rows.append([sc.name, ax, sc.description, "yes" if r.success else f"no ({r.message})",
                                 out.fmt(e.get("rise_10_90_s"), 3), out.fmt(e.get("overshoot_%"), 3),
                                 out.fmt(e.get("settling_2%_s"), 3), out.fmt(e.get("final_%"), 4),
                                 out.fmt(e.get("rms_error_2nd_half_mdeg"), 3), f"{e['max_abs_i_cmd_A']:.3g}",
                                 f"{e['max_abs_v_V']:.3g}", out.fmt(e.get("nonlinear_minus_linear_%"), 3)])
                else:
                    rows.append([sc.name, ax, "(held at the start position)", "yes" if r.success else "no", "-", "-",
                                 "-", "-", f"max error {e['max_abs_error_mdeg']:.3g} mdeg", f"{e['max_abs_i_cmd_A']:.3g}",
                                 f"{e['max_abs_v_V']:.3g}", "-"])
            js["scenarios"][sc.name] = entry
            print(f"  scenario {sc.name}: {r.n_rhs} evaluations, {r.cpu_s:.1f} s CPU, {r.message}")
        md += [out.markdown_table(["scenario", "axis", "description", "completed", "rise 10-90 [s]", "overshoot [%]",
                                   "settling 2% [s]", "final [% of step]", "ramp: rms error, 2nd half [mdeg]",
                                   "max |i_cmd| [A]", "max |v_phase| [V]", "nonlinear - linear [% of step]"], rows), ""]

    md += ["## Assumptions", "",
           "- Drive: series PI loops with the conventions checked on the WorkBench Controller trace; current PI with "
           "its integral zero on the electrical pole R/L (crossover IL.KP/L); continuous-time model, no sampling or "
           "computation delay; observer as a first-order low pass of bandwidth OBS.BW (Velocity Loop page).",
           "- Motor: field-oriented, ideal commutation, per-phase rms model; no cogging; static friction of the "
           "datasheet not modelled (bearing friction is the Dahl element).",
           "- Plant: nonlinear Newton-Euler model of gimbal_mdof; telescope balanced (g = 0); bearing friction in "
           "the time domain as Dahl elements equivalent to the H6 ground constraint at the test amplitude.",
           f"- Measured traces: WorkBench exports; the model of section 3 uses the tuning recorded in them.",
           ""]
    (res_dir / "summary.md").write_text("\n".join(md))
    (res_dir / "summary.json").write_text(json.dumps(js, indent=2, default=float))
    print(f"results in {res_dir} ({time.time() - t_start:.0f} s)")


def _clean(m: dict) -> dict:
    return {k: v for k, v in m.items() if not k.endswith("crossings")} | {
        "gain_crossings": m["gain_crossings"], "phase_crossings": m["phase_crossings"]}


if __name__ == "__main__":
    main()
