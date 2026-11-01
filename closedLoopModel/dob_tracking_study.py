r"""Tracking of a drone on a circle: current control against Hilkert's gyro disturbance observer.

Time-domain closed loop of the LUCY pointer (plant H6, motors, AKD2G drives) with the friction, the cable
spring and the torque scale identified on VClosedLoop/V2, following the reference of Attorrese (2026): drone on
a horizontal circle (radius 10 m, height 5 m, centre 100 m from the pointer, 0.2 rad/s). Two control
structures:

* **current**: position P and velocity PI of the drives closed on the two encoders;
* **proposed**: the same plus a rate gyro on the telescope (ideal) and the single-state disturbance observer of
  Hilkert & Pautler (2011), whose estimate is fed forward to the current command
  (:mod:`lucy_cl.disturbance_observer`).

Metrics on the line of sight of the telescope (fork node for the azimuth, telescope node for the elevation),
after the settling from rest:

* pointing error :math:`E=\lVert\hat n(\psi,\theta)\times\hat r\rVert` (Attorrese, eq. 5.23): rms, mean, max;
* jitter: rms of the error components (cross-elevation, elevation) high-passed above each frequency of
  ``metrics.jitter_bands_hz`` (zero-phase second-order Butterworth).

Also: linear stability of the observer loop (eigenvalues of the linearised model) and disturbance rejection
(bearing torque to line of sight) of the linearised loops.

Usage (from this folder)::

    ../.venv/bin/python -I dob_tracking_study.py --config config_dob_tracking.yaml [--replot]
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import multiprocessing as mp
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from scipy.signal import butter, filtfilt, welch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import main as mn  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl.closed_loop import Scenario, state_space_frf  # noqa: E402
from lucy_cl.disturbance_observer import GyroObserverLoop, ObserverParams, rigid_inertia  # noqa: E402
from lucy_cl.driver import Biquad, Drive  # noqa: E402

MDEG = out.MDEG
plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 6.5,
                     "savefig.bbox": "tight"})


# --------------------------------------------------------------------------- reference


@dataclass(frozen=True)
class Circle:
    """Drone on a horizontal circle; pointer at the origin, azimuth from x, elevation from the horizontal."""

    centre: tuple[float, float]
    radius: float
    height: float
    omega: float
    phase: float = 0.0

    def drone(self, t) -> tuple[np.ndarray, np.ndarray]:
        a = self.omega * np.asarray(t, dtype=float) + self.phase
        p = np.array([self.centre[0] + self.radius * np.cos(a), self.centre[1] + self.radius * np.sin(a),
                      np.full_like(a, self.height)])
        v = np.array([-self.radius * self.omega * np.sin(a), self.radius * self.omega * np.cos(a), np.zeros_like(a)])
        return p, v

    def angles(self, t) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Azimuth, elevation and their rates of the direction to the drone, rad and rad/s."""
        (x, y, z), (vx, vy, _) = self.drone(t)
        rho = np.hypot(x, y)
        drho = (x * vx + y * vy) / rho
        return np.arctan2(y, x), np.arctan2(z, rho), (x * vy - y * vx) / rho**2, -z * drho / (rho**2 + z**2)


@dataclass(frozen=True)
class CircleCommand:
    """Position command of one axis relative to the start, and its rate (the drive feedforward input)."""

    circle: Circle
    axis: str

    def __call__(self, t: float) -> tuple[float, float]:
        psi, th, dpsi, dth = self.circle.angles(t)
        psi0, th0, _, _ = self.circle.angles(0.0)
        return (float(psi - psi0), float(dpsi)) if self.axis == "azimuth" else (float(th - th0), float(dth))


def pointing_error(circle: Circle, t: np.ndarray, psi: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Angle between the line of sight (psi, theta) and the direction to the drone, rad."""
    p, _ = circle.drone(t)
    r = p / np.linalg.norm(p, axis=0)
    n = np.array([np.cos(theta) * np.cos(psi), np.cos(theta) * np.sin(psi), np.sin(theta)])
    return np.arcsin(np.clip(np.linalg.norm(np.cross(n.T, r.T), axis=1), 0.0, 1.0))


# --------------------------------------------------------------------------- model


def build_loop(cc: dict, observer_hz: float | None, feedforward: bool, notch: bool = True,
               friction_structures: dict | None = None, gyro_bandwidth_hz: float | str | None = "config"):
    """Closed loop of the study and its initial state at rest on the reference. ``friction_structures``:
    {axis: structure of v2_identification.py} whose fit replaces the selected one (friction, cable spring,
    torque scale); ``gyro_bandwidth_hz``: first-order gyro, ``"config"`` for the value of the configuration."""
    cfg = yaml.safe_load((HERE / cc["base_config"]).read_text())
    cl, _, _ = mn.build(cfg, HERE)
    pl = cl.plant
    r = cc["reference"]
    circle = Circle(tuple(r["centre_m"]), float(r["radius_m"]), float(r["height_m"]), float(r["omega_rad_s"]),
                    float(r.get("phase_rad", 0.0)))
    psi0, th0, _, _ = circle.angles(0.0)
    pl.theta0 = float(th0)
    gains = json.loads((HERE / cc["drive_gains"]).read_text())["iterations"][-1]["values"]
    drives = {ax: Drive(dataclasses.replace(cl.drives[ax].p, Kpp=gains[ax]["Kpp"], Kpv=gains[ax]["Kpv"],
                                            Kiv_hz=gains[ax]["Kiv"], kvff=1.0 if feedforward else 0.0)) for ax in AXES}
    idj = json.loads((HERE / cc["identified"]).read_text())
    ident = {ax: dict(idj["final"][ax]) for ax in AXES}
    for ax, key in (friction_structures or {}).items():
        ident[ax] = dict(idj["structures"][ax][key]["values"])
    for ax in AXES:
        f = ident[ax]
        pl.friction[:] = [dataclasses.replace(el, sigma0=f["sigma0"], Tc=f["Tc"], Ts=f["Tc"] * f["Ts_ratio"], vs=f["vs"])
                          if el.axis == ax else el for el in pl.friction]
        for n, v in f.items():
            if n.startswith("b_{"):
                pl.set_value(n, float(v))
    oc = cc["observer"]
    obs = {}
    if observer_hz:
        for ax in AXES:
            nc = (oc.get("notch") or {}).get(ax) if notch else None
            filt = (Biquad(float(nc["f_hz"]), float(nc["qz"]), float(nc["f_hz"]), float(nc["qp"])),) if nc else ()
            gb = oc.get("gyro_bandwidth_hz") if gyro_bandwidth_hz == "config" else gyro_bandwidth_hz
            obs[ax] = ObserverParams(J=rigid_inertia(pl, ax), Kt=cl.motors[ax].p.Kt, bandwidth_hz=float(observer_hz),
                                     gyro_bandwidth_hz=gb, ff_filters=filt)
    loop = GyroObserverLoop(pl, cl.motors, drives, obs)
    loop.torque_scale = {ax: float(ident[ax]["scale"]) for ax in AXES}
    loop.cable_spring = {ax: (float(ident[ax]["K_cable"]), float(np.deg2rad(ident[ax]["theta0_cable_deg"]))) for ax in AXES}
    loop.q0[:pl.sim.struttura.n_psi] = psi0
    x0 = loop.initial_state()
    for k, el in enumerate(pl.friction):  # at rest, the friction holds the cable spring
        j = pl.index(el.axis, "bearing")
        K, t0 = loop.cable_spring[el.axis]
        hold = 0.99 * el.level(0.0)
        x0[loop.s_f.start + k] = float(np.clip(-K * (x0[j] - t0), -hold, hold))
    return loop, x0, circle


def run_case(args) -> dict:
    """Simulation of one case (worker): arrays saved to ``runs/<name>.npz``."""
    cc, case, run_dir = args
    loop, x0, circle = build_loop(cc, case.get("observer_hz"), bool(case.get("feedforward")),
                                  friction_structures=case.get("friction_structures"),
                                  gyro_bandwidth_hz=case.get("gyro_bandwidth_hz", "config"))
    ev = np.linalg.eigvals(loop.linearise()[0])
    bad = ev[(ev.real > 1e-6) & (ev.imag >= 0)]
    sim = cc["simulation"]
    t_end = float(sim["t_settle_s"]) + float(sim["orbits"]) * 2 * np.pi / circle.omega
    n_out = int(round(t_end * float(sim["output_rate_hz"]))) + 1
    sc = Scenario(case["name"], {ax: CircleCommand(circle, ax) for ax in AXES}, t_end)
    t0 = time.time()
    r = loop.simulate(sc, n_out=n_out, rtol=float(sim["rtol"]), x0=x0)
    q = r.x[loop.s_q]
    pl = loop.plant
    psi_ref, th_ref, _, _ = circle.angles(r.t)
    data = {"t": r.t, "psi_ref": psi_ref, "th_ref": th_ref, "psi_los": q[loop.jl], "th_los": q[loop.jp],
            "psi_enc": q[loop.ie["azimuth"]], "th_enc": q[loop.ie["elevation"]]}
    for ax in AXES:
        s = r.signals[ax]
        a = ax[:2]
        jb = pl.index(ax, "bearing")
        K, th0 = loop.cable_spring[ax]
        data.update({f"i_{a}": s["i"], f"i_cmd_{a}": s["i_cmd"], f"Td_hat_{a}": s["T_d_hat"], f"w_gyro_{a}": s["w_gyro"],
                     f"T_f_{a}": s["T_f"], f"T_cable_{a}": -K * (q[jb] - th0), f"tau_{a}": s["tau_m"]})
    np.savez_compressed(run_dir / f"{case['name']}.npz", **data)
    return {"name": case["name"], "success": bool(r.success), "message": r.message, "cpu_s": r.cpu_s,
            "wall_s": time.time() - t0, "n_rhs": r.n_rhs, "linear_max_real": float(ev.real.max()),
            "linear_stable": bool(bad.size == 0), "unstable_hz": [float(abs(z.imag) / 2 / np.pi) for z in bad]}


# --------------------------------------------------------------------------- metrics


def metrics(d: dict, circle: Circle, t_settle: float, bands: list[float], kt: dict[str, float]) -> dict:
    t = d["t"]
    m = t >= t_settle
    fs = 1.0 / float(np.median(np.diff(t)))
    E = pointing_error(circle, t, d["psi_los"], d["th_los"])
    dist = np.linalg.norm(circle.drone(t)[0], axis=0)
    e_xel = (d["psi_ref"] - d["psi_los"]) * np.cos(d["th_ref"])
    e_el = d["th_ref"] - d["th_los"]
    res = {"E_rms_mdeg": MDEG * float(np.sqrt(np.mean(E[m] ** 2))), "E_mean_mdeg": MDEG * float(np.mean(E[m])),
           "E_max_mdeg": MDEG * float(np.max(E[m])), "E_p95_mdeg": MDEG * float(np.percentile(E[m], 95)),
           "miss_rms_mm": 1e3 * float(np.sqrt(np.mean((E[m] * dist[m]) ** 2))),
           "xel_rms_mdeg": MDEG * float(np.sqrt(np.mean(e_xel[m] ** 2))), "el_rms_mdeg": MDEG * float(np.sqrt(np.mean(e_el[m] ** 2))),
           "enc_az_rms_mdeg": MDEG * float(np.sqrt(np.mean((d["psi_ref"] - d["psi_enc"])[m] ** 2))),
           "enc_el_rms_mdeg": MDEG * float(np.sqrt(np.mean((d["th_ref"] - d["th_enc"])[m] ** 2))),
           "los_enc_az_rms_mdeg": MDEG * float(np.std((d["psi_los"] - d["psi_enc"])[m])),
           "los_enc_el_rms_mdeg": MDEG * float(np.std((d["th_los"] - d["th_enc"])[m]))}
    for f in bands:
        b, a = butter(2, f / (fs / 2), "high")
        hx, he = filtfilt(b, a, e_xel), filtfilt(b, a, e_el)
        res[f"jitter_{f:g}Hz_urad"] = 1e6 * float(np.sqrt(np.mean(hx[m] ** 2 + he[m] ** 2)))
        res[f"jitter_{f:g}Hz_xel_urad"] = 1e6 * float(np.sqrt(np.mean(hx[m] ** 2)))
        res[f"jitter_{f:g}Hz_el_urad"] = 1e6 * float(np.sqrt(np.mean(he[m] ** 2)))
    for ax in AXES:
        a = ax[:2]
        res[f"i_rms_{a}_A"] = float(np.sqrt(np.mean(d[f"i_{a}"][m] ** 2)))
        res[f"i_max_{a}_A"] = float(np.max(np.abs(d[f"i_{a}"][m])))
        res[f"i_ff_rms_{a}_A"] = float(np.sqrt(np.mean((d[f"Td_hat_{a}"][m] / kt[ax]) ** 2)))
    return res


def stability(cc: dict) -> list[dict]:
    """Largest real part of the eigenvalues of the linearised loop with the observer on one axis, with and
    without the feedforward notch, and the frequency of the unstable modes."""
    rows = []
    for ax in AXES:
        has_notch = bool((cc["observer"].get("notch") or {}).get(ax))
        for notch in ([True, False] if has_notch else [False]):
            for f in cc["observer"]["stability_bandwidths_hz"]:
                loop, _, _ = build_loop(cc, float(f), False, notch=notch)
                loop2 = GyroObserverLoop(loop.plant, loop.motors, loop.drives, {ax: loop.observers[ax]})
                loop2.torque_scale, loop2.cable_spring, loop2.q0 = loop.torque_scale, loop.cable_spring, loop.q0
                e = np.linalg.eigvals(loop2.linearise()[0])
                bad = e[(e.real > 1e-6) & (e.imag >= 0)]
                rows.append({"axis": ax, "notch": notch, "bandwidth_hz": float(f), "max_real": float(e.real.max()),
                             "stable": bool(bad.size == 0),
                             "unstable_hz": [float(abs(z.imag) / 2 / np.pi) for z in bad]})
    return rows


def rejection(cc: dict, f: np.ndarray) -> dict:
    """Line-of-sight angle per N m of torque on the bearing node of each axis, linearised loops with the
    friction sliding (no pre-sliding stiffness): current structure and each observer bandwidth."""
    out_ = {}
    bws = [None] + sorted({c["observer_hz"] for c in cc["cases"] if c.get("observer_hz")})
    for bw in bws:
        loop, _, _ = build_loop(cc, bw, False)
        loop.plant.friction[:] = [dataclasses.replace(el, sigma0=1e-6) for el in loop.plant.friction]
        A, _, _ = loop.linearise()
        M = loop.plant.M(loop.q0)
        for ax in AXES:
            jb = loop.plant.index(ax, "bearing")
            B = np.zeros((loop.n_states, 1))
            e = np.zeros(loop.plant.n)
            e[jb] = 1.0
            B[loop.s_dq, 0] = np.linalg.solve(M, e)
            C = np.zeros((1, loop.n_states))
            C[0, loop.jl if ax == "azimuth" else loop.jp] = 1.0
            out_[(ax, bw)] = state_space_frf(A, B, C, f, 0, 0)
    return out_


# --------------------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config_dob_tracking.yaml")
    ap.add_argument("--replot", action="store_true", help="metrics and figures from the saved runs")
    args = ap.parse_args()
    cc_path = Path(args.config).resolve()
    cc = yaml.safe_load(cc_path.read_text())
    res_dir = HERE / cc["results_dir"]
    run_dir, fig_dir = res_dir / "runs", res_dir / "figures"
    run_dir.mkdir(parents=True, exist_ok=True)
    js_path = res_dir / "summary.json"
    js = json.loads(js_path.read_text()) if js_path.exists() else {}
    cases = cc["cases"]
    t_start = time.time()

    if not args.replot:
        print("stability of the observer loops", flush=True)
        js["stability"] = stability(cc)
        js_path.write_text(json.dumps(js, indent=1))
        todo = [c for c in cases if not (run_dir / f"{c['name']}.npz").exists()]
        print(f"simulating {len(todo)} cases", flush=True)
        with mp.get_context("fork").Pool(min(len(todo), int(cc.get("workers", 8))) or 1) as pool:
            for r in pool.imap_unordered(run_case, [(cc, c, run_dir) for c in todo]):
                js.setdefault("runs", {})[r["name"]] = r
                js_path.write_text(json.dumps(js, indent=1))
                print(f"  {r['name']}: {r['message']} (cpu {r['cpu_s']:.0f} s)", flush=True)

    loop, _, circle = build_loop(cc, None, False)
    kt = {ax: loop.motors[ax].p.Kt for ax in AXES}
    J = {ax: rigid_inertia(loop.plant, ax) for ax in AXES}
    t_settle = float(cc["simulation"]["t_settle_s"])
    bands = [float(b) for b in cc["metrics"]["jitter_bands_hz"]]
    data = {c["name"]: dict(np.load(run_dir / f"{c['name']}.npz")) for c in cases}
    js["metrics"] = {n: metrics(d, circle, t_settle, bands, kt) for n, d in data.items()}
    js["observer"] = {ax: {"J": J[ax], "Kt": kt[ax], "torque_scale": loop.torque_scale[ax]} for ax in AXES}
    f = np.logspace(-1, np.log10(500), 400)
    rej = rejection(cc, f)
    js["rejection_low_frequency"] = {f"{ax} / {bw}": float(np.abs(v[np.argmin(np.abs(f - 1.0))])) for (ax, bw), v in rej.items()}
    js_path.write_text(json.dumps(js, indent=1))

    # --------------------------------------------------------------------------- tables
    lab = {c["name"]: c["label"] for c in cases}
    M = js["metrics"]
    md = ["# Tracking of the circle: current control against the gyro disturbance observer", "",
          f"Configuration `{cc_path.name}`. Reference of Attorrese (2026): circle of radius {cc['reference']['radius_m']:g} m "
          f"at {cc['reference']['height_m']:g} m, centre {cc['reference']['centre_m'][0]:g} m away, "
          f"{cc['reference']['omega_rad_s']:g} rad/s. Metrics on the line of sight of the telescope from "
          f"{t_settle:g} s ({cc['simulation']['orbits']} orbit). Jitter: rms of the line-of-sight error above the given "
          "frequency (1 mdeg = 17.45 urad). Observer: J' = " + ", ".join(f"{ax} {J[ax]:.3f} kg m^2" for ax in AXES) + "; K_t' = "
          + f"{kt['azimuth']:.3f} N m/A (true torque scale " + ", ".join(f"{ax} {loop.torque_scale[ax]:.2f}" for ax in AXES) + ").", ""]
    rows = []
    for n in lab:
        x = M[n]
        rows.append([lab[n], f"{x['E_rms_mdeg']:.2f}", f"{x['E_mean_mdeg']:.2f}", f"{x['E_max_mdeg']:.2f}",
                     f"{x['miss_rms_mm']:.1f}"] + [f"{x[f'jitter_{b:g}Hz_urad']:.3g}" for b in bands]
                    + [f"{x['xel_rms_mdeg']:.2f}", f"{x['el_rms_mdeg']:.2f}", f"{x['i_rms_az_A']:.2f}", f"{x['i_rms_el_A']:.2f}"])
    md += ["## Pointing error and jitter", "", out.markdown_table(
        ["case", "E rms [mdeg]", "E mean [mdeg]", "E max [mdeg]", "miss at the drone, rms [mm]"]
        + [f"jitter > {b:g} Hz [urad]" for b in bands] + ["cross-elevation rms [mdeg]", "elevation rms [mdeg]",
                                                          "current az rms [A]", "current el rms [A]"], rows), ""]
    rows = [[r["axis"], "yes" if r["notch"] else "no", f"{r['bandwidth_hz']:g}", f"{r['max_real']:+.2f}",
             "yes" if r["stable"] else "no", ", ".join(f"{x:.0f}" for x in r["unstable_hz"]) or "-"] for r in js["stability"]]
    md += ["## Stability of the observer loop", "", "Largest real part of the eigenvalues [1/s] of the linearised "
           "loop with the observer on one axis.", "",
           out.markdown_table(["axis", "notch on the feedforward", "observer bandwidth [Hz]", "max real part [1/s]",
                               "stable", "unstable modes [Hz]"], rows), ""]
    (res_dir / "summary.md").write_text("\n".join(md))

    # --------------------------------------------------------------------------- figures
    col = {"encoders": "k", "encoders_ff": "k"}
    # reference
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    tt = np.linspace(0, 2 * np.pi / circle.omega, 600)
    (x, y, z), _ = circle.drone(tt)
    axs[0].plot(x, y, "C0", label="drone on the circle")
    axs[0].plot([0], [0], "k^", label="pointer")
    axs[0].plot([0, x[0]], [0, y[0]], "k:", lw=0.8, label="line of sight at t = 0")
    axs[0].set_aspect("equal", adjustable="datalim")
    axs[0].set_xlabel("x [m]")
    axs[0].set_ylabel("y [m]")
    axs[0].set_title("top view")
    axs[0].legend(loc="upper left")
    psi, th, dpsi, dth = circle.angles(tt)
    axs[1].plot(tt, np.rad2deg(psi), "C0", label="azimuth")
    axs[1].plot(tt, np.rad2deg(th), "C3", label="elevation")
    axs[1].set_xlabel("time [s]")
    axs[1].set_ylabel("angle [deg]")
    axs[1].set_title("reference angles")
    axs[1].legend(loc="best")
    axs[2].plot(tt, np.rad2deg(dpsi), "C0", label="azimuth")
    axs[2].plot(tt, np.rad2deg(dth), "C3", label="elevation")
    axs[2].set_xlabel("time [s]")
    axs[2].set_ylabel("rate [deg/s]")
    axs[2].set_title("reference rates")
    axs[2].legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "reference")
    # pointing error, the compared pairs
    pairs = cc["compare"]
    fig, axs = plt.subplots(len(pairs), 1, figsize=(10, 3.0 * len(pairs)), squeeze=False)
    for a, (n0, n1) in zip(axs[:, 0], pairs):
        for n, c in ((n0, "k"), (n1, "C3")):
            d = data[n]
            E = pointing_error(circle, d["t"], d["psi_los"], d["th_los"])
            a.plot(d["t"], MDEG * E, color=c, lw=0.8, label=f"{lab[n]} (rms {M[n]['E_rms_mdeg']:.2f} mdeg)")
        a.axvline(t_settle, color="0.5", ls=":", lw=0.8, label="end of the settling (start of the metrics)")
        a.set_yscale("log")
        a.set_xlabel("time [s]")
        a.set_ylabel("pointing error [mdeg]")
        a.legend(loc="upper right")
    fig.tight_layout()
    out.save(fig, fig_dir, "pointing_error")
    # error components and jitter of the first pair
    for k, (n0, n1) in enumerate(pairs):
        fig, axs = plt.subplots(3, 2, figsize=(11, 8))
        for n, c in ((n0, "k"), (n1, "C3")):
            d = data[n]
            t = d["t"]
            fs = 1.0 / float(np.median(np.diff(t)))
            m = t >= t_settle
            comps = (("cross-elevation", (d["psi_ref"] - d["psi_los"]) * np.cos(d["th_ref"])),
                     ("elevation", d["th_ref"] - d["th_los"]))
            for j, (nm, e) in enumerate(comps):
                axs[0, j].plot(t, MDEG * e, color=c, lw=0.7, label=lab[n])
                b, a_ = butter(2, bands[0] / (fs / 2), "high")
                axs[1, j].plot(t[m], 1e6 * filtfilt(b, a_, e)[m], color=c, lw=0.6, label=lab[n])
                fw, P = welch(e[m] - np.mean(e[m]), fs=fs, nperseg=min(8192, int(m.sum())))
                axs[2, j].loglog(fw[1:], 1e6 * np.sqrt(P[1:]), color=c, lw=0.8, label=lab[n])
                axs[0, j].set_title(f"{nm}: line-of-sight error")
                axs[1, j].set_title(f"{nm}: error above {bands[0]:g} Hz (jitter)")
                axs[2, j].set_title(f"{nm}: amplitude spectral density of the error")
        for j in range(2):
            axs[0, j].set_ylabel("error [mdeg]")
            axs[1, j].set_ylabel("high-passed error [urad]")
            axs[2, j].set_ylabel("ASD [urad/sqrt(Hz)]")
            axs[2, j].set_xlabel("frequency [Hz]")
            for a in axs[:2, j]:
                a.set_xlabel("time [s]")
            for a in axs[:, j]:
                a.legend(loc="best")
        fig.tight_layout()
        out.save(fig, fig_dir, f"components_{n1}")
    # disturbance estimate against the physical disturbance, main observer case
    n1 = pairs[0][1]
    d = data[n1]
    fig, axs = plt.subplots(2, 1, figsize=(10, 5.5))
    for a, ax in zip(axs, AXES):
        s = ax[:2]
        a.plot(d["t"], -d[f"T_f_{s}"] + d[f"T_cable_{s}"], "k", lw=0.9, label="bearing friction + cable spring (model)")
        a.plot(d["t"], d[f"Td_hat_{s}"], "C3", lw=0.7, label=f"observer estimate T_d' ({lab[n1]})")
        a.set_ylabel("torque [N m]")
        a.set_xlabel("time [s]")
        a.set_title(f"{ax}: disturbance torque on the axis")
        a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "disturbance_estimate")
    # sweep of the observer bandwidth
    sweep = [c for c in cases if not c.get("feedforward")]
    xs = [c.get("observer_hz") or 0.0 for c in sweep]
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.5))
    keys = [("E_rms_mdeg", "pointing error rms [mdeg]")] + [(f"jitter_{b:g}Hz_urad", f"jitter above {b:g} Hz [urad]") for b in bands]
    for a, (k_, yl) in zip(axs, keys):
        a.plot(xs[1:], [M[c["name"]][k_] for c in sweep[1:]], "o-", color="C3", label="with the gyro observer")
        a.axhline(M[sweep[0]["name"]][k_], color="k", ls="--", label="current structure")
        a.set_xscale("log")
        a.set_xticks(xs[1:], [f"{x:g}" for x in xs[1:]])
        a.minorticks_off()
        a.set_xlabel("observer bandwidth [Hz]")
        a.set_ylabel(yl)
        a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "observer_sweep")
    # linear disturbance rejection
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for a, ax in zip(axs, AXES):
        for (ax_, bw), v in rej.items():
            if ax_ != ax:
                continue
            a.loglog(f, 1e6 * np.abs(v), color="k" if bw is None else None, lw=1.2 if bw is None else 0.9,
                     label="current structure" if bw is None else f"with the observer, {bw:g} Hz")
        a.set_xlabel("frequency [Hz]")
        a.set_ylabel("line of sight per bearing torque [urad/(N m)]")
        a.set_title(f"{ax}: disturbance rejection (linearised, sliding friction)")
        a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "disturbance_rejection")
    # motor currents of the first pair
    n0, n1 = pairs[0]
    fig, axs = plt.subplots(2, 1, figsize=(10, 5.5))
    for a, ax in zip(axs, AXES):
        s = ax[:2]
        for n, c in ((n0, "k"), (n1, "C3")):
            a.plot(data[n]["t"], data[n][f"i_{s}"], color=c, lw=0.7, label=f"motor current, {lab[n]}")
        a.plot(data[n1]["t"], -data[n1][f"Td_hat_{s}"] / kt[ax], color="C0", lw=0.7, label="observer feedforward -T_d'/K_t'")
        a.set_ylabel("current [A]")
        a.set_xlabel("time [s]")
        a.set_title(f"{ax}: motor current")
        a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "currents")
    print(f"done in {(time.time() - t_start) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
