"""One axis in closed loop: AKD2G drive -> motor -> modal structure, with the Dahl friction of the bearings outside.

State ``[x, xdot, i, drive]``: node rotations and rates of the modal model, motor current (A rms) and drive states.
Inputs ``[u, T_b]``: position command (cascade, the drive in position mode) or velocity command (the drive in velocity
mode, position loop in the master), and the bearing torque. Without saturations the drive is linear: its matrices are
obtained from :meth:`lucy_sim.plant.drive.Drive.step` by finite differences. The linear part is discretised exactly
(ZOH); the Dahl state is updated exactly at constant velocity over the step.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm, solve

from .drive import Drive
from .fork import ModalModel
from .friction import DahlFriction
from .motor import MotorParams


class VelocityModeDrive:
    """The drive in velocity mode: the first input of :class:`ClosedLoopAxis` is the velocity command."""

    def __init__(self, d: Drive):
        self.d, self.p, self.n_states, self.i_zv = d, d.p, d.n_states, d.i_zv

    def step(self, x, v_cmd, _dp, _p_fb, v_fb, i_meas):
        return self.d.step(x, 0.0, v_cmd, 0.0, v_fb, i_meas)


class ClosedLoopAxis:
    """Axis in closed loop: continuous matrices ``A, B``, discrete ``Ad, Bd`` at ``dt`` and output rows."""

    def __init__(self, m: ModalModel, drive: Drive | VelocityModeDrive, motor: MotorParams, friction: DahlFriction, dt: float):
        n, nd = len(m.M), drive.n_states
        self.m, self.drive, self.motor, self.friction = m, drive, motor, friction
        # linear drive: d xd/dt = A xd + B u, v_phase = Cv xd + Dv u, i_cmd = Ci xd + Di u; u = [cmd, p_fb, v_fb, i]
        x0, u0 = np.zeros(nd), np.zeros(4)

        def f(xd, u):
            v, ic, dx = drive.step(xd, u[0], 0.0, u[1], u[2], u[3])
            return np.r_[dx, v, ic]
        y0 = f(x0, u0)
        h = 1e-6
        Jx = np.column_stack([(f(x0 + h * np.eye(nd)[k], u0) - y0) / h for k in range(nd)])
        Ju = np.column_stack([(f(x0, u0 + h * np.eye(4)[k]) - y0) / h for k in range(4)])
        Ad, Bd, Cv, Dv, Ci, Di = Jx[:nd], Ju[:nd], Jx[nd], Ju[nd], Jx[nd + 1], Ju[nd + 1]
        C = m.damping_matrix
        Mi = np.linalg.inv(m.M)
        N = 2 * n + 1 + nd
        ix, iv, ii, idr = slice(0, n), slice(n, 2 * n), 2 * n, slice(2 * n + 1, N)
        A = np.zeros((N, N))
        B = np.zeros((N, 2))
        A[ix, iv] = np.eye(n)
        A[iv, ix] = -Mi @ m.K
        A[iv, iv] = -Mi @ C - motor.b_viscous * np.outer(Mi @ m.e, m.e)
        A[iv, ii] = Mi @ m.e * m.Kt
        B[iv, 1] = -Mi @ m.v
        # drive inputs from the state: p_fb = c x, v_fb = c xdot, i
        G = np.zeros((4, N))
        G[1, ix] = m.c
        G[2, iv] = m.c
        G[3, ii] = 1.0
        g = np.zeros((4, 2))
        g[0, 0] = 1.0
        A[idr, :] += Bd @ G
        A[idr, idr] += Ad
        B[idr, :] += Bd @ g
        # motor: L di/dt = v_phase - R i - Ke omega_rel, omega_rel = e^T xdot (rotor relative to the stator)
        rv = np.zeros(N)
        rv[idr] = Cv
        rv += Dv @ G
        rv[ii] -= motor.R_ph
        rv[iv] -= motor.Ke_ph * m.e
        A[ii, :] = rv / motor.L_ph
        B[ii, :] = (Dv @ g) / motor.L_ph
        self.A, self.B = A, B
        self.out_position = np.r_[m.c, np.zeros(N - n)]
        """Encoder reading."""
        self.out_current = np.eye(N)[ii]
        self.row_i_cmd = np.zeros(N)
        self.row_i_cmd[idr] = Ci
        self.row_i_cmd += Di @ G
        self.in_i_cmd = Di @ g
        self.row_v_phase = rv + motor.R_ph * np.eye(N)[ii] + motor.Ke_ph * np.r_[np.zeros(n), m.e, np.zeros(N - 2 * n)]
        self.in_v_phase = Dv @ g
        self.row_i_int = np.zeros(N)
        wi = 2 * np.pi * drive.p.Kiv_hz
        self.row_i_int[2 * n + 1 + drive.i_zv] = drive.p.Kpv * wi
        E = expm(np.block([[A, B], [np.zeros((2, N + 2))]]) * dt)
        self.Ad, self.Bd, self.N, self.n, self.dt = E[:N, :N], E[:N, N:], N, n, dt

    @property
    def bearing_rate_row(self) -> np.ndarray:
        """Bearing deformation rate :math:`v^T\\dot x` from the state."""
        return np.r_[np.zeros(self.n), self.m.v, np.zeros(self.N - 2 * self.n)]

    def equilibrium(self, Tb: float) -> np.ndarray:
        """State at rest with zero command and bearing torque ``Tb`` (the velocity integrator holds it)."""
        return solve(self.A, -self.B @ np.array([0.0, Tb]))

    def run(self, p_cmd: np.ndarray, I0: float, torque: np.ndarray | None = None, los: np.ndarray | None = None) -> dict[str, np.ndarray]:
        """Cascade (drive in position mode) on the command ``p_cmd`` sampled at ``dt``, relative to the start.

        ``I0``: holding current at the start (the friction holds it, limited to T_c). ``torque``: disturbance torque on
        the bearing port, step by step. ``los``: extra output row on the node positions (line of sight), relative to the
        start. Outputs: encoder ``p``, current ``i``, current command, phase voltage, velocity integrator, ``los``."""
        fr, s0, n = self.friction, self.friction.sigma0, self.n
        Tb0 = float(np.clip(self.m.Kt * I0, -fr.Tc, fr.Tc))
        z = Tb0 / s0
        s = self.equilibrium(Tb0)
        p0 = self.out_position @ s
        P, I, IC, VP, VI, L = (np.empty(len(p_cmd)) for _ in range(6))
        vv = self.bearing_rate_row
        row = None if los is None else np.r_[los[:n], np.zeros(self.N - n)]
        l0 = row @ s if row is not None else 0.0
        for k, pc in enumerate(p_cmd):
            u = np.array([pc, s0 * z + (torque[k] if torque is not None else 0.0)])
            P[k], I[k] = self.out_position @ s - p0, self.out_current @ s
            L[k] = row @ s - l0 if row is not None else P[k]
            IC[k], VP[k], VI[k] = self.row_i_cmd @ s + self.in_i_cmd @ u, self.row_v_phase @ s + self.in_v_phase @ u, self.row_i_int @ s
            s_new = self.Ad @ s + self.Bd @ u
            vd = 0.5 * (vv @ s + vv @ s_new)
            a = s0 * abs(vd) / fr.g(vd)
            if a > 0:
                e = np.exp(-a * self.dt)
                z = z * e + vd / a * (1 - e)
            s = s_new
        return {"p": P, "i": I, "i_cmd": IC, "v_phase": VP, "i_int": VI, "Tb0": Tb0, "los": L}
