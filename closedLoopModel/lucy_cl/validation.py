r"""Validation of the closed loop against the time sequences recorded on LUCY (``LucyIdentification/Validation``).

The sequences (2 kHz, about 100 s each) hold, for both axes, the command position and the encoder position
(millidegrees, 360000 per revolution), the current (``AXIS*_AMPS``) and the torque estimated by the drive
(``AXIS*_TORQUE``); ``AXIS1`` is the azimuth, ``AXIS2`` the elevation. The channels start a few samples apart:
the comparison uses the rows where all of them are present.

Current scale. ``AXIS*_AMPS`` saturates at 35.9, twice the current limit of the drives (``IL.LIMITP`` = 18 A
rms), and the rigid-body balance of the moves gives, with the inertias of the model, 0.124-0.158 N m per unit,
about :math:`K_t/2`: the unit is taken as 0.5 A rms (``amps_per_A`` = 2). The torque channel is
``0.1847 x AMPS`` in every file.

The model is driven by the measured command of both axes (linear interpolation) from a state at rest built from
the measurement at the start of the window: positions, holding currents, friction torques and integrators
consistent with them (:func:`state_at_rest`).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import AXES
from .closed_loop import ClosedLoop, Scenario, SimResult

MDEG_TO_RAD = np.pi / 180.0 / 1000.0
COLUMNS = {"azimuth": {"enc": 1, "cmd": 7, "amps": 6, "torque": 3},
           "elevation": {"enc": 2, "cmd": 8, "amps": 5, "torque": 4}}


@dataclass
class Sequence:
    name: str
    t: np.ndarray
    """s, from the first complete row."""
    cmd: dict[str, np.ndarray]
    enc: dict[str, np.ndarray]
    """rad."""
    current: dict[str, np.ndarray]
    """A rms (``AMPS / amps_per_A``)."""
    torque: dict[str, np.ndarray]
    """drive estimate, N m."""
    amps_per_A: float

    def window(self, t0: float, t1: float) -> np.ndarray:
        return (self.t >= t0) & (self.t <= t1)

    def moves(self, axis: str, min_gap: float = 0.3) -> list[tuple[float, float]]:
        """Start and end times of the command moves of ``axis``."""
        v = np.abs(np.diff(self.cmd[axis])) > 1e-9
        idx = np.nonzero(v)[0]
        if idx.size == 0:
            return []
        out, a = [], idx[0]
        for k0, k1 in zip(idx[:-1], idx[1:]):
            if self.t[k1] - self.t[k0] > min_gap:
                out.append((self.t[a], self.t[k0 + 1]))
                a = k1
        out.append((self.t[a], self.t[idx[-1] + 1]))
        return out


def read_sequence(path: str | Path, amps_per_A: float = 2.0) -> Sequence:
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split(",")[:9]
            rows.append([float(x) if x.strip() else np.nan for x in p] + [np.nan] * (9 - len(p)))
    D = np.array(rows)
    ok = ~np.isnan(D).any(axis=1)
    D = D[ok]
    t = (D[:, 0] - D[0, 0]) / 1000.0
    seq = Sequence(Path(path).stem, t, {}, {}, {}, {}, amps_per_A)
    for ax in AXES:
        c = COLUMNS[ax]
        seq.cmd[ax] = D[:, c["cmd"]] * MDEG_TO_RAD
        seq.enc[ax] = D[:, c["enc"]] * MDEG_TO_RAD
        seq.current[ax] = D[:, c["amps"]] / amps_per_A
        seq.torque[ax] = D[:, c["torque"]]
    return seq


class MeasuredCommand:
    """Measured command of one axis as a function of the simulation time (``t0`` = start of the window),
    relative to the start position of the model node ``p_start``: :math:`(p - p_{start}, \\dot p)`."""

    def __init__(self, t: np.ndarray, p: np.ndarray, t0: float, p_start: float):
        self.t, self.p, self.t0, self.p_start = t - t0, p, 0.0, p_start
        self.v = np.gradient(p, t)

    def __call__(self, t: float) -> tuple[float, float]:
        return float(np.interp(t, self.t, self.p) - self.p_start), float(np.interp(t, self.t, self.v))


def state_at_rest(cl: ClosedLoop, p_enc: dict[str, float], i_hold: dict[str, float]) -> np.ndarray:
    r"""Equilibrium at rest at the measured encoder positions: the bearing friction (Dahl state) carries the
    holding torque :math:`K_t i` of the measured current, clipped at :math:`0.99\,T_c` (gravity is off); the
    motor current is the one balancing it (when the measured holding torque exceeds the sliding friction of
    the model the axis would otherwise start moving at the window start); the shafts carry the static twist of
    that torque between the motor node and the bearing node, with the encoder node at the measured position;
    the velocity integrator gives that current command, the current integrator the voltage :math:`R i`."""
    pl = cl.plant
    x = np.zeros(cl.n_states)
    q = np.zeros(pl.n)
    n_psi = pl.sim.struttura.n_psi
    q[:n_psi] = p_enc["azimuth"]
    q[n_psi:] = p_enc["elevation"]
    _, _, K = pl.linear_matrices()
    for j, ax in enumerate(AXES):
        m, d = cl.motors[ax], cl.drives[ax]
        i = float(i_hold[ax])
        tau = m.p.Kt * i
        for kf, el in enumerate(pl.friction):
            if el.axis == ax:
                tau = float(np.clip(tau, -0.99 * el.Tc, 0.99 * el.Tc))
                x[cl.s_f.start + kf] = tau
                i = tau / m.p.Kt
        x[cl.i_cur[ax]] = i
        # static twist of the axis: motor torque on the motor node, reaction of the friction on the bearing
        # node (held there), encoder node at the measured position
        sl = list(range(0, n_psi)) if ax == "azimuth" else list(range(n_psi, pl.n))
        jb = pl.index(ax, "bearing")
        free = [k for k in sl if k != jb]
        dq = np.zeros(pl.n)
        dq[free] = np.linalg.solve(K[np.ix_(free, free)], pl.E[free, j] * tau)
        q[sl] += dq[sl] - dq[pl.index(ax, "encoder")]
        p = d.p
        xd = np.zeros(d.n_states)
        if p.Kiv_hz > 0:
            xd[d.i_zv] = i / (p.Kpv * 2 * np.pi * p.Kiv_hz)
        xd[d.i_zi] = i * p.L_motor / p.Kpi
        x[cl.s_drv[ax]] = xd
    x[cl.s_q] = q
    return x


def rigid_body_regression(seq: Sequence, axis: str, J: float, fc: float = 3.0) -> dict[str, float]:
    """Current against the rigid-body balance during the moves: :math:`i = a\\,\\ddot q + b\\,
    \\mathrm{sgn}\\,\\dot q + c` on the encoder signal low-passed at ``fc`` (moving samples, :math:`|\\dot q|`
    > 1 deg/s). With the inertia ``J`` of the model, :math:`J/a` is the torque per unit of the current
    channel and :math:`bJ/a` the Coulomb friction."""
    from scipy.signal import butter, filtfilt

    dt = float(np.median(np.diff(seq.t)))
    b, a = butter(2, fc * 2 * dt)
    th = filtfilt(b, a, seq.enc[axis])
    i = filtfilt(b, a, seq.current[axis] * seq.amps_per_A)
    w = np.gradient(th, dt)
    al = np.gradient(w, dt)
    m = np.abs(w) > np.deg2rad(1.0)
    X = np.c_[al[m], np.sign(w[m]), np.ones(m.sum())]
    c, *_ = np.linalg.lstsq(X, i[m], rcond=None)
    r = i[m] - X @ c
    return {"slope": float(c[0]), "coulomb_units": float(c[1]), "offset_units": float(c[2]),
            "torque_per_unit": J / c[0], "coulomb_Nm": float(c[1] * J / c[0]), "R2": float(1 - r.var() / i[m].var())}


def simulate_window(cl: ClosedLoop, seq: Sequence, t0: float, t1: float, *, n_avg: int = 200,
                    **kw) -> tuple[SimResult, np.ndarray]:
    """Simulates the window ``[t0, t1]`` of the sequence; returns the result and the measured sample mask.
    The start state is at rest with the encoder positions and currents averaged over the ``n_avg`` samples
    before ``t0``."""
    k0 = int(np.searchsorted(seq.t, t0))
    pre = slice(max(0, k0 - n_avg), max(1, k0))
    p0 = {ax: float(np.mean(seq.enc[ax][pre])) for ax in AXES}
    i0 = {ax: float(np.mean(seq.current[ax][pre])) for ax in AXES}
    cl.plant.theta0 = p0["elevation"]
    cl.q0 = cl.plant.q_equilibrium()
    cl.q0[: cl.plant.sim.struttura.n_psi] = p0["azimuth"]
    x0 = state_at_rest(cl, p0, i0)
    m = seq.window(t0, t1)
    cmds = {ax: MeasuredCommand(seq.t[m], seq.cmd[ax][m], t0, cl.q0[cl.ie[ax]]) for ax in AXES}
    sc = Scenario(f"{seq.name}_{t0:.2f}", cmds, t1 - t0, f"{seq.name}, {t0:.2f}-{t1:.2f} s")
    n_out = int(round((t1 - t0) / 5e-4)) + 1
    r = cl.simulate(sc, n_out=n_out, x0=x0, **kw)
    return r, m
