r"""Estimate of the effective drive gains (and of the bearing friction level) from the time sequences.

The drive configuration active during the sequences is not available (user, 2026-10-19: "stima dai dati"); the
gains are estimated starting from the .akd files. Free parameters of each axis:

* ``Kpp``  -- position gain, 1/s;
* ``Kpv``  -- velocity proportional gain, A/(rad/s);
* ``Kiv``  -- velocity integral gain, Hz;
* ``Tc``   -- Coulomb level of the Dahl bearing friction, N m (the level fitted on the frequency responses at
  sub-milliradian amplitude, 1.6 N m azimuth and 0.66 N m elevation, is far above the friction during the
  moves of the sequences, about 0.3 N m); :math:`\sigma_0` stays the H6 value.

The filters, the observer, the current loop and the limits stay those of the .akd files. Cost: for each window
and channel (position error ``command - encoder`` and current of the moving axis) the residual normalised as in
the NRMSE fit, :math:`(\hat y - y)/\lVert y - \bar y\rVert`, so that the cost is :math:`\sum (1 - fit)^2`. The
windows are simulated in parallel (one closed loop per worker process).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import AXES
from .closed_loop import ClosedLoop
from .driver import Drive
from .validation import Sequence, read_sequence, simulate_window

PARAMS = ("Kpp", "Kpv", "Kiv", "Tc")
"""Default parameters; lumped plant parameters can be added by symbol (e.g. ``c_{\\psi,2}``)."""
BOUNDS = {"Kpp": (0.3, 200.0), "Kpv": (1.0, 2000.0), "Kiv": (0.02, 60.0), "Tc": (0.02, 3.0)}
DECIMATE = 4
"""Residuals every 4 samples (2 ms)."""


@dataclass(frozen=True)
class Window:
    sequence: str
    t0: float
    t1: float
    axis: str
    """Axis whose move the window contains."""
    move: float = float("nan")
    """Start of the move, s."""

    @property
    def label(self) -> str:
        return f"{self.sequence} {self.axis} move at {self.move:.2f} s (window {self.t0:.2f}-{self.t1:.2f} s)"


def windows_for_moves(seq: Sequence, axis: str, *, pre: float = 0.3, post: float = 2.5,
                      starts: list[float] | None = None, min_rest: float = 0.6) -> list[Window]:
    """One window per move of ``axis`` (or per move starting near ``starts``): from ``pre`` seconds before the
    move to ``post`` seconds after its start, cut 0.05 s before the next move of either axis."""
    moves = seq.moves(axis)
    every = [m for ax in AXES for m in seq.moves(ax)]
    other = [m[0] for m in every]
    out = []
    for a, b in moves:
        if starts is not None and not any(abs(a - s) < 0.05 for s in starts):
            continue
        # at rest before the move: no move of either axis (the other axis included) in [a - min_rest, a)
        busy = any(s < a - 1e-6 and e > a - min_rest for (s, e) in every)
        if busy or a - pre < seq.t[0] + 0.1:
            continue
        nxt = min([s for s in other if s > a + 1e-6] or [seq.t[-1]])
        out.append(Window(seq.name, a - pre, min(a + post, nxt - 0.05), axis, float(a)))
    return out


def apply_params(cl: ClosedLoop, axis: str, values: dict[str, float]) -> ClosedLoop:
    """Closed loop with the parameters of ``axis`` replaced (drive gains, Dahl ``Tc``)."""
    p = cl.drives[axis].p
    rep = {}
    if "Kpp" in values:
        rep["Kpp"] = values["Kpp"]
    if "Kpv" in values:
        rep["Kpv"] = values["Kpv"]
    if "Kiv" in values:
        rep["Kiv_hz"] = values["Kiv"]
    if "kvff" in values:
        rep["kvff"] = values["kvff"]
    drives = dict(cl.drives) | {axis: Drive(dataclasses.replace(p, **rep))}
    pl = cl.plant
    for name, v in values.items():  # lumped plant parameters (dampers, springs) by symbol
        if name.startswith(("c_{", "k_{")):
            pl.set_value(name, float(v))
    if "Tc" in values:
        pl.friction[:] = [dataclasses.replace(el, Tc=values["Tc"]) if el.axis == axis else el for el in pl.friction]
    return ClosedLoop(pl, cl.motors, drives)


def current_values(cl: ClosedLoop, axis: str) -> dict[str, float]:
    p = cl.drives[axis].p
    tc = [el.Tc for el in cl.plant.friction if el.axis == axis]
    return {"Kpp": p.Kpp, "Kpv": p.Kpv, "Kiv": p.Kiv_hz, "Tc": tc[0] if tc else float("nan"), "kvff": p.kvff}


# --------------------------------------------------------------------------- worker processes

_W: dict = {}


def _init_worker(cfg: dict, base: str, seq_dir: str, amps_per_A: float) -> None:
    import main as mn  # the entry point builds the components from the configuration

    cl, _, _ = mn.build(cfg, Path(base))
    _W["cl"] = cl
    _W["seq"] = {}
    _W["dir"] = Path(seq_dir)
    _W["amps"] = amps_per_A


def _sequence(name: str) -> Sequence:
    if name not in _W["seq"]:
        _W["seq"][name] = read_sequence(_W["dir"] / f"{name}.csv", _W["amps"])
    return _W["seq"][name]


def simulate_job(job: tuple[Window, dict[str, dict[str, float]], float]) -> dict:
    """Worker: simulates one window with the given parameters (per axis); returns the measured and simulated
    signals of both axes on the decimated measurement grid."""
    win, values, max_cpu = job
    cl = _W["cl"]
    saved = {ax: current_values(cl, ax) for ax in AXES}
    for ax, v in values.items():  # plant parameters changed by the job are restored afterwards
        saved[ax] |= {n: cl.plant.value(n) for n in v if n.startswith(("c_{", "k_{"))}
    try:
        for ax, v in values.items():
            cl = apply_params(cl, ax, v)
        seq = _sequence(win.sequence)
        r, m = simulate_window(cl, seq, win.t0, win.t1, max_cpu_s=max_cpu)
        tm = seq.t[m][::DECIMATE] - win.t0
        out = {"ok": r.success, "t": tm, "message": r.message, "cpu": r.cpu_s}
        for ax in AXES:
            s = r.signals[ax]
            out[ax] = {
                "error_meas": (seq.cmd[ax][m] - seq.enc[ax][m])[::DECIMATE],
                "current_meas": seq.current[ax][m][::DECIMATE],
                "enc_meas": (seq.enc[ax][m] - seq.enc[ax][m][0])[::DECIMATE],
                "error": np.interp(tm, r.t, s["error"], right=np.nan),
                "current": np.interp(tm, r.t, s["i"], right=np.nan),
                "enc": np.interp(tm, r.t, s["p_enc"] - (seq.enc[ax][m][0] - cl.q0[cl.ie[ax]]), right=np.nan),
            }
        return out
    finally:
        _restore(saved)


def _restore(saved: dict[str, dict[str, float]]) -> None:
    cl = _W["cl"]
    for ax, v in saved.items():
        cl = apply_params(cl, ax, v)
    _W["cl"] = cl


# --------------------------------------------------------------------------- estimate


def normalised_residual(meas: np.ndarray, model: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(meas - np.mean(meas))
    r = (model - meas) / max(d, 1e-12)
    return np.where(np.isfinite(r), r, 1.0 / np.sqrt(max(r.size, 1)))  # unfinished samples count as fit 0


@dataclass
class Estimate:
    axis: str
    start: dict[str, float]
    values: dict[str, float]
    cost_start: float
    cost: float
    n_eval: int
    windows: list[Window]
    history: list[tuple[dict[str, float], float]]


def estimate(pool, axis: str, windows: list[Window], start: dict[str, float], *, weights=(1.0, 1.0),
             max_nfev: int = 40, max_cpu: float = 120.0, log=print,
             fixed: dict[str, dict[str, float]] | None = None) -> Estimate:
    """Least squares on the logarithms of the parameters of ``axis`` (trust region, finite differences).
    ``fixed``: values of the other parameters of both axes during the estimate (default: those of the loop
    built by the workers)."""
    from scipy.optimize import least_squares

    names = [n for n in PARAMS if n in start] + [n for n in start if n.startswith(("c_{", "k_{"))]
    bnd = {n: BOUNDS.get(n, (start[n] * 1e-3, start[n] * 1e3)) for n in names}
    lo = np.log([bnd[n][0] for n in names])
    hi = np.log([bnd[n][1] for n in names])
    x0 = np.clip(np.log([start[n] for n in names]), lo + 1e-9, hi - 1e-9)
    history: list[tuple[dict[str, float], float]] = []

    def residual(x: np.ndarray) -> np.ndarray:
        vals = dict(zip(names, np.exp(x)))
        allv = {ax: dict(v) for ax, v in (fixed or {}).items()}
        allv[axis] = allv.get(axis, {}) | vals
        res = pool.map(simulate_job, [(w, allv, max_cpu) for w in windows])
        r = []
        for w, o in zip(windows, res):
            s = o[axis]
            r.append(weights[0] * normalised_residual(s["error_meas"], s["error"]))
            r.append(weights[1] * normalised_residual(s["current_meas"], s["current"]))
        r = np.concatenate(r)
        cost = float(np.sum(r**2))
        history.append((vals, cost))
        log(f"    {axis}: " + ", ".join(f"{k} {v:.4g}" for k, v in vals.items()) + f" -> cost {cost:.4f}")
        return r

    r0 = residual(x0)
    sol = least_squares(residual, x0, bounds=(lo, hi), diff_step=0.03, x_scale=1.0, max_nfev=max_nfev,
                        xtol=1e-3, ftol=1e-4)
    vals = dict(zip(names, np.exp(sol.x)))
    return Estimate(axis, dict(start), vals, float(np.sum(r0**2)), float(np.sum(sol.fun**2)), sol.nfev + 1,
                    windows, history)


def cruise_gain(seq: Sequence, axis: str, v_min: float = np.deg2rad(10.0)) -> float:
    """Initial guess of the position gain from the following error at constant commanded speed,
    :math:`K_{pp} \\approx \\dot p_{cmd}/(p_{cmd} - p_{enc})`."""
    v = np.gradient(seq.cmd[axis], seq.t)
    e = seq.cmd[axis] - seq.enc[axis]
    a = np.gradient(v, seq.t)
    m = (np.abs(v) > v_min) & (np.abs(a) < 1e-3) & (np.abs(e) > 1e-6)
    return float(np.median(v[m] / e[m])) if m.sum() > 20 else float("nan")
