"""Measured V2 doublets of LUCY: windows, metrics and simulation of the model on the measured position command.

Raw tests: ``MPCLucy_doublets_<axis>_<amplitude>.csv`` exported by the PLC (time in ms, positions in mdeg, current
``AXIS*_AMPS``, which is 2 I_rms; AXIS1 = azimuth, AXIS2 = polar). Each test has three motions (+A, -2A, +A); every
motion is a window. The tuning uses one window per amplitude (the *typical commands*), stored in
``data/doublets/<axis>.npz`` by :func:`save_typical_commands` (``scripts/extract_doublets.py``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.signal import butter, filtfilt

DT = 5e-4
"""Sampling of the windows, s."""
AMPS_PER_A = 2.0
COLUMNS = {"azimuth": ("AzimuthAxis/ActualPosition", "AzimuthAxis/CommandPosition", "VM/VARIABLESTOMONITOR.AXIS1_AMPS"),
           "polar": ("PolarAxis/ActualPosition", "PolarAxis/CommandPosition", "VM/VARIABLESTOMONITOR.AXIS2_AMPS")}
FILES = {"azimuth": "MPCLucy_doublets_azimut_%s.csv", "polar": "MPCLucy_doublets_elevation_%s.csv"}
AMPLITUDES = ("01deg_minus01deg", "05deg_minus05deg", "1deg_minus1deg", "3deg_minus3deg", "5deg_minus5deg", "8deg_minus8deg")
ANOMALOUS = ("azimuth", "8deg_minus8deg", "#1")
"""Window skipped among the typical commands."""
SETTLING_BAND = 5.0
"""Error band of the settling time, mdeg."""
_FILTER = butter(4, 150 * 2 * DT)


@dataclass
class Window:
    name: str
    t: np.ndarray
    I: np.ndarray
    """Current, A rms (filtered at 150 Hz)."""
    pos: np.ndarray
    """Encoder, rad, relative to the start of the window."""
    cmd: np.ndarray
    """Position command, rad, relative to the start of the encoder."""


def load_csv(path: Path, axis: str) -> dict[str, np.ndarray]:
    with open(path) as fh:
        cols = [c.strip() for c in fh.readline().split(",")]
    A = np.genfromtxt(path, delimiter=",", skip_header=1, usecols=range(len(cols)), filling_values=np.nan)
    ip, ic, ii = (cols.index(c) for c in COLUMNS[axis])
    A = A[~np.isnan(A[:, [ip, ic, ii]]).any(axis=1)]
    return {"t": (A[:, 0] - A[0, 0]) * 1e-3, "pos": np.deg2rad(A[:, ip] * 1e-3), "cmd": np.deg2rad(A[:, ic] * 1e-3), "I": A[:, ii] / AMPS_PER_A}


def windows(path: Path, axis: str, before: float = 0.5, after: float = 12.0) -> list[Window]:
    """One window per motion of the command: from ``before`` s before to ``after`` s after (or to the next motion)."""
    d = load_csv(path, axis)
    t = d["t"]
    idx = np.flatnonzero(np.abs(np.diff(d["cmd"])) > 1e-12)
    starts = idx[np.r_[True, np.diff(idx) > 200]]
    out = []
    for k, s in enumerate(starts):
        t0 = t[s] - before
        t1 = min(t[s] + after, t[starts[k + 1]] - before if k + 1 < len(starts) else t[-1])
        m = (t >= t0) & (t <= t1)
        tt = np.arange(t0, t1, DT)
        I = filtfilt(*_FILTER, np.interp(tt, t[m], d["I"][m]))
        p, c = np.interp(tt, t[m], d["pos"][m]), np.interp(tt, t[m], d["cmd"][m])
        out.append(Window(f"{path.stem.removeprefix('MPCLucy_doublets_')} #{k + 1}", tt - tt[0], I, p - p[0], c - p[0]))
    return out


def typical_commands_from_csv(folder: Path, axis: str) -> list[Window]:
    """One doublet window per amplitude, the first one that is not anomalous."""
    out = []
    for amp in AMPLITUDES:
        ok = [w for w in windows(Path(folder) / (FILES[axis] % amp), axis) if (axis, amp, w.name[-2:]) != ANOMALOUS]
        out.append(ok[0])
    return out


def save_typical_commands(ws: list[Window], path: Path) -> None:
    np.savez_compressed(path, names=np.array([w.name for w in ws]), **{f"{k}_{i}": getattr(w, k) for i, w in enumerate(ws) for k in ("t", "I", "pos", "cmd")})


def load_typical_commands(path: Path) -> list[Window]:
    d = np.load(path)
    return [Window(str(n), d[f"t_{i}"], d[f"I_{i}"], d[f"pos_{i}"], d[f"cmd_{i}"]) for i, n in enumerate(d["names"])]


def _motion(w: Window) -> tuple[int, int]:
    moving = np.abs(np.gradient(w.cmd, DT)) > 1e-6
    i0, i1 = np.flatnonzero(moving)[[0, -1]]
    return i0, i1


def doublet_metrics(w: Window, y: np.ndarray, i: np.ndarray) -> dict[str, float]:
    """Tracking error during the motion and 1 s after (rms, max), settling time in the 5 mdeg band after the end of
    the motion, overshoot, current."""
    i0, i1 = _motion(w)
    mv = slice(i0, min(i1 + int(1.0 / DT), len(y)))
    e = np.rad2deg((w.cmd - w.cmd[0]) - y) * 1e3
    out = np.flatnonzero(np.abs(e[i1:]) > SETTLING_BAND)
    final = np.rad2deg(w.cmd[i1] - w.cmd[0]) * 1e3
    pos = np.rad2deg(y[i1:]) * 1e3
    stroke = np.rad2deg(w.cmd[i1] - w.cmd[i0]) * 1e3
    over = float(max(0.0, np.max((pos - final) * np.sign(stroke))))
    return {"stroke_mdeg": abs(stroke), "err_max_mdeg": float(np.abs(e[mv]).max()), "err_rms_mdeg": float(np.sqrt(np.mean(e[mv] ** 2))),
            "t_set_s": float((out[-1] + 1) * DT) if out.size else 0.0, "over_mdeg": over,
            "i_max_A": float(np.abs(i[mv]).max()), "i_rms_A": float(np.sqrt(np.mean(i[mv] ** 2)))}


def simulate_doublet(model, axis: str, c, w: Window, t_max: float | None = None) -> dict:
    """The model on the measured position command of ``w`` (cascade with ``c = None``); outputs on the grid of ``w``,
    starting from the measured holding current."""
    from ..simulation.loops import run_lqg

    a = model.axes[axis]
    k = w.t <= (t_max if t_max else w.t[-1])
    t = np.arange(0.0, w.t[k][-1], model.dt)
    pc, I0 = np.interp(t, w.t, w.cmd - w.cmd[0]), float(w.I[:int(0.4 / DT)].mean())
    r = a.cascade.run(pc, I0) if c is None else run_lqg(a, c, pc, I0, model.dt_master, model.delay)
    return {n: np.interp(w.t[k], t, r[n]) for n in ("p", "i")} | {f"max_{n}": float(np.abs(r[n]).max()) for n in ("i_cmd", "v_phase", "i_int")}
