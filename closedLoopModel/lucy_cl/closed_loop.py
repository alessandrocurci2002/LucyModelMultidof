r"""Closed loop: mechanical plant + two motors + two drives.

Signal flow of each axis (no algebraic loop: the drive acts on the measured current, the motor on the voltage)::

    p_cmd --> [drive] -- v_phase --> [motor] -- tau_m --> [plant] -- p_fb, v_fb (load-side encoder) --+
                ^  ^                    |  ^                 |                                        |
                |  +------ i (current) -+  +-- omega_m ------+                                        |
                +-----------------------------------------------------------------------------------+

State vector :math:`\mathbf x = [\mathbf q,\ \dot{\mathbf q},\ \mathbf T_f,\ i_\psi,\ i_\theta,\ \mathbf x_{d,\psi},\
\mathbf x_{d,\theta}]`: plant coordinates and velocities, Dahl friction torques, motor currents, drive states.
The rotor inertia of each motor is added to the motor node of the plant when the loop is built.

Frequency domain (linear, about rest at ``theta0``), per axis with :math:`G_{em}`, :math:`G_{mm}` the encoder and
motor rotations per N m of motor torque:

.. math::
    T_i = \frac{i}{i_{cmd}} = \frac{C_i}{L_{ph}s + R_{ph} + C_i + K_{e}K_t\,s\,G_{mm}},\qquad
    P = \frac{v_{fb}}{i_{cmd}} = s\,G_{em}K_t\,T_i,

the WorkBench traces ``Controller`` :math:`C = C_vH_{fwd}(1 + C_p/s)`, ``Open Loop``
:math:`L_{wb} = CPH_{AR3}H_{AR4}` and ``Closed Loop`` :math:`CP/(1 + L_{wb})` (checked on the exports to machine
precision; the observer is not part of them), and the exact loop broken at the current command, with the position fed back without the velocity
feedback path :math:`H_{fb}` (biquads AR3, AR4 and observer): :math:`L_x = C_vH_{fwd}(C_p + sH_{fb})P/s`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from . import AXES
from .driver import Drive, controller_workbench
from .motor import Motor
from .plant import GimbalPlant


# --------------------------------------------------------------------------- commands and scenarios


@dataclass(frozen=True)
class Command:
    """Position command of one axis, relative to the start position: ``step`` (``value`` in rad) or ``ramp``
    (``value`` in rad/s), starting at ``t0``."""

    kind: str = "hold"
    value: float = 0.0
    t0: float = 0.0

    def __call__(self, t: float) -> tuple[float, float]:
        if self.kind == "hold" or t < self.t0:
            return 0.0, 0.0
        if self.kind == "step":
            return self.value, 0.0
        if self.kind == "ramp":
            return self.value * (t - self.t0), self.value
        raise ValueError(f"unknown command {self.kind!r}")


@dataclass(frozen=True)
class Scenario:
    name: str
    commands: Mapping[str, Command]
    t_end: float
    description: str = ""
    compare_linear: bool = False
    """Overlay the response of the linearised model (meaningful in the pre-sliding range of the friction)."""

    def command(self, axis: str) -> Command:
        return self.commands.get(axis, Command())

    @property
    def breakpoints(self) -> list[float]:
        t0s = {getattr(c, "t0", 0.0) for c in self.commands.values()}
        return sorted(t for t in t0s if 0.0 < t < self.t_end)


@dataclass
class SimResult:
    scenario: Scenario
    t: np.ndarray
    x: np.ndarray
    signals: dict[str, dict[str, np.ndarray]]
    """Per axis: ``p_cmd``, ``p_enc``, ``p_motor``, ``p_load`` (azimuth: fork node; elevation: telescope node),
    ``error``, ``i_cmd``, ``i``, ``v``, ``tau_m``, ``T_f`` -- positions relative to the start."""
    n_rhs: int
    cpu_s: float
    success: bool
    message: str


# --------------------------------------------------------------------------- closed loop


class ClosedLoop:
    """The three components connected; time-domain simulation and linearisation."""

    def __init__(self, plant: GimbalPlant, motors: Mapping[str, Motor], drives: Mapping[str, Drive]):
        self.plant, self.motors, self.drives = plant, dict(motors), dict(drives)
        for ax in AXES:
            plant.add_inertia(ax, "motor", self.motors[ax].p.J_rotor, "rotor")
        n, nf = plant.n, plant.n_friction
        self.s_q, self.s_dq, self.s_f = slice(0, n), slice(n, 2 * n), slice(2 * n, 2 * n + nf)
        k = 2 * n + nf
        self.i_cur, self.s_drv = {}, {}
        for ax in AXES:
            self.i_cur[ax] = k
            k += 1
        for ax in AXES:
            self.s_drv[ax] = slice(k, k + self.drives[ax].n_states)
            k += self.drives[ax].n_states
        self.n_states = k
        self.ie = {ax: plant.index(ax, "encoder") for ax in AXES}
        self.im = {ax: plant.index(ax, "motor") for ax in AXES}
        self.q0 = plant.q_equilibrium()
        self.torque_scale = {ax: 1.0 for ax in AXES}
        """Motor torque per :math:`K_t i` of each axis (scale identified by the current-driven model)."""
        self.cable_spring: dict[str, tuple[float, float]] = {}
        """{axis: (K_c [N m/rad], theta_0 [rad])}: torque :math:`-K_c(\\theta - \\theta_0)` on the bearing node."""
        self.signal_names = ("p_cmd", "i_cmd", "v", "i", "tau_m")
        """Internal signals returned by :meth:`evaluate`, in order."""

    # ---- time domain -----------------------------------------------------------------------------------

    def initial_state(self) -> np.ndarray:
        x = np.zeros(self.n_states)
        x[self.s_q] = self.q0
        return x

    def jacobian_steps(self) -> np.ndarray:
        """Perturbations of the finite-difference Jacobian given to the solver: 1e-6 of the state or of a
        typical scale of its group (scipy's adaptive steps grow without bound on the columns that do not act,
        e.g. the velocity integrator while the current is saturated)."""
        typ = np.full(self.n_states, 1e-4)
        typ[self.s_q] = 1e-3
        typ[self.s_dq] = 1e-2
        typ[self.s_f] = 1e-1
        for ax in AXES:
            typ[self.i_cur[ax]] = 1e-1
        return typ

    def atol(self) -> np.ndarray:
        """Absolute tolerances by group: rad, rad/s, N m, A, drive states (units of their inputs)."""
        a = np.full(self.n_states, 1e-10)
        a[self.s_dq] = 1e-8
        a[self.s_f] = 1e-7
        for ax in AXES:
            a[self.i_cur[ax]] = 1e-6
        return a

    def evaluate(self, t: float, x: np.ndarray, cmd: Mapping[str, Callable[[float], tuple[float, float]]],
                 i_exc: Mapping[str, Callable[[float], float]] | None = None) -> tuple[np.ndarray, dict]:
        """State derivative and the internal signals at ``(t, x)``."""
        pl = self.plant
        q, dq, Tf = x[self.s_q], x[self.s_dq], x[self.s_f]
        dx = np.empty(self.n_states)
        tau = np.zeros(2)
        sig = {}
        for j, ax in enumerate(AXES):
            d, m = self.drives[ax], self.motors[ax]
            ie, im = self.ie[ax], self.im[ax]
            p_ref, dp_ref = cmd[ax](t) if ax in cmd else (0.0, 0.0)
            p_cmd = self.q0[ie] + p_ref
            i = x[self.i_cur[ax]]
            v, i_cmd, dxd = d.step(x[self.s_drv[ax]], p_cmd, dp_ref, q[ie], dq[ie], i,
                                   i_exc[ax](t) if i_exc and ax in i_exc else 0.0)
            dx[self.s_drv[ax]] = dxd
            dx[self.i_cur[ax]] = m.current_derivative(i, v, dq[im])
            tau[j] = m.torque(self.torque_scale[ax] * i, q[im], dq[im])
            sig[ax] = (p_ref, i_cmd, v, i, tau[j])
        dx[self.s_q] = dq
        dx[self.s_dq] = pl.accelerations(q, dq, tau, Tf, self.extra_torques(q) if self.cable_spring else None)
        dx[self.s_f] = pl.friction_derivatives(dq, Tf)
        return dx, sig

    def extra_torques(self, q: np.ndarray) -> np.ndarray:
        """Torques on the plant coordinates other than the motors and the friction: the cable springs."""
        Q = np.zeros(self.plant.n)
        for ax, (K, th0) in self.cable_spring.items():
            j = self.plant.index(ax, "bearing")
            Q[j] -= K * (q[j] - th0)
        return Q

    def simulate(self, sc: Scenario, *, n_out: int = 3001, rtol: float = 1e-4, atol_scale: float = 100.0,
                 method: str = "Radau", max_step: float = np.inf, max_cpu_s: float | None = None,
                 x0: np.ndarray | None = None) -> SimResult:
        """Integrates the nonlinear closed loop (implicit solver: the encoder nodes are very stiff), from rest
        at the operating point or from ``x0``. Tolerances: ``rtol`` and :meth:`atol` x ``atol_scale``
        (defaults: 1e-8 rad, 1e-6 rad/s, 1e-5 N m, 1e-4 A); tighter values stall the solver when a bearing
        sticks after sliding (velocity at the scale of the friction smoothing). With ``max_cpu_s`` the
        integration stops when the CPU time exceeds it (e.g. an unstable loop in a limit cycle) and the result
        holds the samples computed so far."""
        import time

        cmd = {ax: sc.command(ax) for ax in AXES}
        count = [0]

        def f(t: float, x: np.ndarray) -> np.ndarray:
            count[0] += 1
            return self.evaluate(t, x, cmd)[0]

        typ = self.jacobian_steps()

        def jac(t: float, x: np.ndarray) -> np.ndarray:
            f0 = f(t, x)
            h = 1e-6 * np.maximum(np.abs(x), typ)
            J = np.empty((x.size, x.size))
            for j in range(x.size):
                xp = x.copy()
                xp[j] += h[j]
                J[:, j] = (f(t, xp) - f0) / h[j]
            return J

        t_out = np.linspace(0.0, sc.t_end, n_out)
        # integration pieces: between the command discontinuities and, with a CPU budget, of at most
        # ``chunk`` seconds so that the budget is checked regularly; each piece ends exactly at its edge
        chunk = sc.t_end if max_cpu_s is None else min(sc.t_end, 0.1)
        edges = np.unique(np.r_[0.0, sc.breakpoints, np.arange(chunk, sc.t_end, chunk), sc.t_end])
        x0 = self.initial_state() if x0 is None else np.asarray(x0, dtype=float).copy()
        ts, xs = [], []
        atol = self.atol() * atol_scale
        ok, msg = True, ""
        t_cpu = time.process_time()
        for a, b in zip(edges[:-1], edges[1:]):
            te = np.unique(np.r_[a, t_out[(t_out > a) & (t_out < b)], b])
            sol = solve_ivp(f, (a, b), x0, method=method, t_eval=te, rtol=rtol, atol=atol, max_step=max_step,
                            **({"jac": jac} if method in ("Radau", "BDF") else {}))
            ok, msg = ok and sol.success, sol.message
            keep = np.isin(sol.t, t_out) & ~np.isin(sol.t, np.concatenate(ts) if ts else [])
            ts.append(sol.t[keep])
            xs.append(sol.y[:, keep])
            if not ok:
                break
            x0 = sol.y[:, -1]
            if max_cpu_s is not None and time.process_time() - t_cpu > max_cpu_s and b < sc.t_end:
                ok, msg = False, f"stopped at t = {b:.3g} s: CPU time above {max_cpu_s:g} s"
                break
        t = np.concatenate(ts)
        x = np.concatenate(xs, axis=1)
        cpu = time.process_time() - t_cpu
        return SimResult(sc, t, x, self._signals(t, x, cmd), count[0], cpu, ok, msg)

    def _signals(self, t: np.ndarray, x: np.ndarray, cmd) -> dict[str, dict[str, np.ndarray]]:
        out = {ax: {k: np.empty(t.size) for k in self.signal_names} for ax in AXES}
        for k in range(t.size):
            _, sig = self.evaluate(t[k], x[:, k], cmd)
            for ax in AXES:
                for name, val in zip(self.signal_names, sig[ax]):
                    out[ax][name][k] = val
        q = x[self.s_q]
        for ax in AXES:
            ie = self.ie[ax]
            o = out[ax]
            o["p_enc"] = q[ie] - self.q0[ie]
            o["error"] = o["p_cmd"] - o["p_enc"]
            o["p_motor"] = q[self.im[ax]] - self.q0[self.im[ax]]
            o["T_f"] = np.zeros(t.size)
            for kf, el in enumerate(self.plant.friction):
                if el.axis == ax:
                    o["T_f"] = x[self.s_f][kf]
        jp = self.plant.index("elevation", "telescope")
        out["elevation"]["p_load"] = q[jp] - self.q0[jp]
        jl = self.plant.index("azimuth", "bearing")
        out["azimuth"]["p_load"] = q[jl] - self.q0[jl]
        return out

    # ---- linearisation of the time-domain model --------------------------------------------------------

    def linearise(self, eps: float = 1e-8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """:math:`\\mathbf A`, :math:`\\mathbf B` (inputs: position commands of azimuth and elevation) and
        :math:`\\mathbf C` (outputs: encoder positions) of the time-domain model about rest, by central
        differences; the Dahl elements linearise to springs :math:`\\sigma_0`."""
        x0 = self.initial_state()
        zero = {ax: (lambda t: (0.0, 0.0)) for ax in AXES}
        A = np.zeros((self.n_states, self.n_states))
        for j in range(self.n_states):
            e = np.zeros(self.n_states)
            e[j] = eps
            A[:, j] = (self.evaluate(0.0, x0 + e, zero)[0] - self.evaluate(0.0, x0 - e, zero)[0]) / (2 * eps)
        B = np.zeros((self.n_states, 2))
        for j, ax in enumerate(AXES):
            up = dict(zero) | {ax: (lambda t, s=+eps: (s, 0.0))}
            dn = dict(zero) | {ax: (lambda t, s=-eps: (s, 0.0))}
            B[:, j] = (self.evaluate(0.0, x0, up)[0] - self.evaluate(0.0, x0, dn)[0]) / (2 * eps)
        C = np.zeros((2, self.n_states))
        for j, ax in enumerate(AXES):
            C[j, self.ie[ax]] = 1.0
        return A, B, C


# --------------------------------------------------------------------------- frequency domain


@dataclass
class LoopFRF:
    """Linear responses of one axis about rest (see the module docstring)."""

    axis: str
    f: np.ndarray
    ground: str
    plant: np.ndarray
    """:math:`v_{fb}/i_{cmd}`, (rad/s)/A."""
    current: np.ndarray
    controller: np.ndarray
    feedback: np.ndarray
    """Feedback biquads AR3, AR4 (WorkBench ``Feedback Biquad``)."""
    open_loop: np.ndarray
    """WorkBench definition, :math:`C P H_{AR3}H_{AR4}` (no observer)."""
    closed_loop: np.ndarray
    """WorkBench definition, :math:`CP/(1 + CPH_{AR3}H_{AR4})`."""
    loop: np.ndarray
    """Exact loop at the current command, observer included."""
    position: np.ndarray
    """:math:`p_{fb}/p_{cmd}`."""
    extra: dict = field(default_factory=dict)


def loop_frf(cl: ClosedLoop, axis: str, f: np.ndarray, *, ground: str = "hysteretic",
             drive: Drive | None = None) -> LoopFRF:
    """Frequency responses of ``axis`` with ``drive`` (default: the drive of the loop). ``ground`` as in
    :meth:`GimbalPlant.frf`."""
    pl, m = cl.plant, cl.motors[axis]
    d = cl.drives[axis] if drive is None else drive
    f = np.asarray(f, dtype=float)
    s = 2j * np.pi * f
    M, C, K = pl.linear_matrices()
    K = K.astype(complex)
    im, ie = cl.im[axis], cl.ie[axis]
    C[im, im] += m.p.b_viscous
    if ground == "hysteretic":
        for ax in AXES:
            g = pl.ground(ax)
            if g is not None:
                j = pl.index(ax, "bearing")
                K[j, j] += g[0] * (1 + 1j * g[1])
                C[j, j] += g[2]
    elif ground == "dahl":
        for el in pl.friction:
            K[el.index, el.index] += el.sigma0
    elif ground != "none":
        raise ValueError(ground)
    col = pl.E[:, 0 if axis == "azimuth" else 1]
    Gem, Gmm = np.empty(f.size, complex), np.empty(f.size, complex)
    for k, sk in enumerate(s):
        x = np.linalg.solve(M * sk**2 + C * sk + K, col)
        Gem[k], Gmm[k] = x[ie], x[im]
    p = m.p
    Ci = d.C_cur(s)
    Ti = Ci / (p.L_ph * s + p.R_ph + Ci + p.Ke_ph * p.Kt * s * Gmm)
    P = s * Gem * p.Kt * Ti * d.delay(s)
    ctrl, fb_wb = controller_workbench(d, s), d.F_ar_fb(s)
    ol = ctrl * P * fb_wb
    Lx = d.C_vel(s) * d.F_fwd(s) * (d.C_pos(s) + s * d.F_fb(s)) * P / s
    pos = d.C_vel(s) * d.F_fwd(s) * d.C_pos(s) * P / s / (1 + Lx)
    return LoopFRF(axis, f, ground, P, Ti, ctrl, fb_wb, ol, ctrl * P / (1 + ol), Lx, pos)


def margins(f: np.ndarray, L: np.ndarray) -> dict:
    """Gain and phase crossings of a loop transfer :math:`L`: phase margins at :math:`|L| = 1`, gain margins
    at :math:`\\angle L = -180^\\circ + k\\,360^\\circ` (log-frequency interpolation)."""
    mag_db = 20 * np.log10(np.abs(L))
    ph = np.rad2deg(np.unwrap(np.angle(L)))
    lf = np.log10(f)
    gc, pc = [], []
    for k in range(f.size - 1):
        if np.sign(mag_db[k]) != np.sign(mag_db[k + 1]):
            r = mag_db[k] / (mag_db[k] - mag_db[k + 1])
            fx = 10 ** (lf[k] + r * (lf[k + 1] - lf[k]))
            phx = ph[k] + r * (ph[k + 1] - ph[k])
            gc.append((float(fx), float((phx + 360.0) % 360.0 - 180.0)))
        n0, n1 = np.floor((ph[k] + 180.0) / 360.0), np.floor((ph[k + 1] + 180.0) / 360.0)
        if n0 != n1:
            target = 360.0 * max(n0, n1) - 180.0
            r = (ph[k] - target) / (ph[k] - ph[k + 1])
            fx = 10 ** (lf[k] + r * (lf[k + 1] - lf[k]))
            gm = -(mag_db[k] + r * (mag_db[k + 1] - mag_db[k]))
            pc.append((float(fx), float(gm)))
    pm = min(gc, key=lambda c: abs(c[1])) if gc else (float("nan"), float("nan"))
    gm = min(pc, key=lambda c: c[1]) if pc else (float("nan"), float("inf"))
    return {"gain_crossings": gc, "phase_crossings": pc, "pm_deg": pm[1], "f_pm": pm[0], "gm_db": gm[1],
            "f_gm": gm[0], "f_c": gc[0][0] if gc else float("nan")}


def bandwidth(f: np.ndarray, T: np.ndarray, level_db: float = -3.0) -> float:
    """First frequency where :math:`|T|` falls below ``level_db`` (relative to its low-frequency value)."""
    g = 20 * np.log10(np.abs(T))
    g = g - g[0]
    idx = np.nonzero(g < level_db)[0]
    if idx.size == 0:
        return float("nan")
    k = idx[0]
    if k == 0:
        return float(f[0])
    r = (g[k - 1] - level_db) / (g[k - 1] - g[k])
    return float(10 ** (np.log10(f[k - 1]) + r * (np.log10(f[k]) - np.log10(f[k - 1]))))


def state_space_frf(A: np.ndarray, B: np.ndarray, C: np.ndarray, f: np.ndarray, i_in: int, i_out: int) -> np.ndarray:
    n = A.shape[0]
    out = np.empty(f.size, complex)
    for k, fk in enumerate(f):
        out[k] = C[i_out] @ np.linalg.solve(2j * np.pi * fk * np.eye(n) - A, B[:, i_in])
    return out
