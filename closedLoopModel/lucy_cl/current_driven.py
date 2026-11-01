r"""Plant driven by the measured motor currents (no drive in the loop).

The motor torque of each axis is :math:`\tau_m = K_t\,i_{meas}(t) - b_m\,\omega_m` (measured current of the
sequence, linear interpolation), and the plant with the Dahl bearing friction gives the encoder positions:

.. math::
    \mathbf M(\theta_p)\ddot{\mathbf q} = \mathbf E\boldsymbol\tau_m + \mathbf Q_f - \mathbf h(\mathbf q,
    \dot{\mathbf q}),\qquad \dot T_f = \sigma_0\bigl(\dot q - |\dot q|_\varepsilon T_f/T_c\bigr).

Optional torques on the bearing node of an axis (``hysteron``, ``position_torque``):

* two-level hysteresis (relay hysteron): :math:`-H\,s`, with :math:`\dot s = (\dot q - |\dot q|_\varepsilon s)/\delta`,
  so that :math:`s` goes to :math:`\pm 1` after a displacement of a few :math:`\delta` in either direction and
  keeps its value at rest (memory of the direction of the last move);
* position-dependent torque :math:`a\sin\theta + b\cos\theta + c` of the absolute angle of the node;
* elastic torque of the cables :math:`-K_c(\theta - \theta_0)`, with :math:`\theta_0` the neutral angle set by
  the layout of the cables (an unknown of the identification).

This tests the mechanics and the friction alone: the controller does not take part. A window starts at rest
from the measurement (encoder positions, friction torques carrying the holding currents, static twist of the
shafts), as in :func:`lucy_cl.validation.state_at_rest`.
"""

from __future__ import annotations

import time
from dataclasses import replace

import numpy as np
from scipy.integrate import solve_ivp

from . import AXES
from .motor import Motor
from .plant import GimbalPlant
from .validation import Sequence


class CurrentDrivenPlant:
    r"""``torque_scale`` per axis multiplies :math:`K_t\,i_{meas}` (scale of the current channel, tolerance of
    :math:`K_t`, inertia of the CAD); ``input_lowpass_hz`` filters the measured current (zero phase) before it
    drives the model."""

    def __init__(self, plant: GimbalPlant, motors: dict[str, Motor], input_lowpass_hz: float | None = 150.0):
        self.plant, self.motors = plant, dict(motors)
        self.torque_scale = {ax: 1.0 for ax in AXES}
        self.input_lowpass_hz = input_lowpass_hz
        self.hysteron: dict[str, tuple[float, float]] = {}
        """{axis: (H [N m], delta [rad])} two-level hysteresis on the bearing node."""
        self.position_torque: dict[str, tuple[float, float, float]] = {}
        """{axis: (a, b, c) [N m]} torque a sin(theta) + b cos(theta) + c on the bearing node."""
        self.cable_spring: dict[str, tuple[float, float]] = {}
        """{axis: (K_c [N m/rad], theta_0 [rad])} elastic torque -K_c (theta - theta_0) on the bearing node."""
        for ax in AXES:
            plant.add_inertia(ax, "motor", self.motors[ax].p.J_rotor, "rotor")
        n = plant.n
        self.s_q, self.s_dq = slice(0, n), slice(n, 2 * n)
        self.s_f = slice(2 * n, 2 * n + plant.n_friction)
        self.i_h = {ax: 2 * n + plant.n_friction + k for k, ax in enumerate(AXES)}
        """State of the hysteron of each axis (unused when the axis has none)."""
        self.n_states = 2 * n + plant.n_friction + len(AXES)
        self.ie = {ax: plant.index(ax, "encoder") for ax in AXES}
        self.im = {ax: plant.index(ax, "motor") for ax in AXES}

    # ---- parameters of the friction ------------------------------------------------------------------

    def set_friction(self, axis: str, **kw) -> None:
        """Changes ``sigma0``, ``Tc``, ``Ts``, ``vs`` of the friction element of ``axis``."""
        pl = self.plant
        kw = {k: (None if v is None else float(v)) for k, v in kw.items()}
        pl.friction[:] = [replace(el, **kw) if el.axis == axis else el for el in pl.friction]

    def friction(self, axis: str):
        return [el for el in self.plant.friction if el.axis == axis][0]

    # ---- simulation ------------------------------------------------------------------------------------

    def rhs(self, t: float, x: np.ndarray, current) -> np.ndarray:
        pl = self.plant
        q, dq, Tf = x[self.s_q], x[self.s_dq], x[self.s_f]
        tau = np.array([self.motors[ax].torque(self.torque_scale[ax] * current[ax](t), q[self.im[ax]], dq[self.im[ax]])
                        for ax in AXES])
        dx = np.zeros(self.n_states)
        dx[self.s_q] = dq
        Q = pl.E @ tau - pl.h(q, dq) + self.extra_torques(q, x)
        for el, tf in zip(pl.friction, Tf):
            Q[el.index] -= tf
        dx[self.s_dq] = np.linalg.solve(pl.M(q), Q)
        dx[self.s_f] = pl.friction_derivatives(dq, Tf)
        for ax, (H, delta) in self.hysteron.items():
            v = dq[pl.index(ax, "bearing")]
            sv = x[self.i_h[ax]]
            dx[self.i_h[ax]] = (v - (np.sqrt(v * v + 1e-12) - 1e-6) * sv) / delta
        return dx

    def extra_torques(self, q: np.ndarray, x: np.ndarray) -> np.ndarray:
        Q = np.zeros(self.plant.n)
        for ax in AXES:
            j = self.plant.index(ax, "bearing")
            Q[j] += self.axis_extra(ax, q[j], x[self.i_h[ax]])[0]
        return Q

    def axis_extra(self, ax: str, th: float, s: float) -> tuple[float, float, float]:
        """Extra torque on the bearing node of ``ax`` at the absolute angle ``th`` (hysteron state ``s``), and
        its derivatives with respect to ``th`` and ``s``."""
        T = dth = ds = 0.0
        if ax in self.hysteron:
            T -= self.hysteron[ax][0] * s
            ds = -self.hysteron[ax][0]
        if ax in self.position_torque:
            a, b, c = self.position_torque[ax]
            T += a * np.sin(th) + b * np.cos(th) + c
            dth += a * np.cos(th) - b * np.sin(th)
        if ax in self.cable_spring:
            K, th0 = self.cable_spring[ax]
            T -= K * (th - th0)
            dth -= K
        return T, dth, ds

    def initial_state(self, p_enc: dict[str, float], i_hold: dict[str, float],
                      s_hyst: dict[str, float] | None = None) -> np.ndarray:
        """At rest at the measured encoder positions. Static balance of each axis: motor torque of the holding
        current on the motor node, gravity (when the model has it) on the telescope, reaction of the bearing
        friction on the bearing node, which carries the net torque (clipped below its break-away level); the
        shafts take the static twist of these loads."""
        pl = self.plant
        x = np.zeros(self.n_states)
        q = np.zeros(pl.n)
        n_psi = pl.sim.struttura.n_psi
        q[:n_psi] = p_enc["azimuth"]
        q[n_psi:] = p_enc["elevation"]
        _, _, K = pl.linear_matrices()
        for ax in AXES:
            x[self.i_h[ax]] = float((s_hyst or {}).get(ax, 0.0)) if ax in self.hysteron else 0.0
        # gravity (and any static term) and the extra torques at the rigid position
        g_load = -pl.h(q, np.zeros(pl.n)) + self.extra_torques(q, x)
        for j, ax in enumerate(AXES):
            sl = list(range(0, n_psi)) if ax == "azimuth" else list(range(n_psi, pl.n))
            tau = self.torque_scale[ax] * self.motors[ax].p.Kt * float(i_hold[ax])
            F = pl.E[:, j] * tau + g_load
            net = float(np.sum(F[sl]))
            for kf, el in enumerate(pl.friction):
                if el.axis == ax:
                    hold = 0.99 * el.level(0.0)
                    x[self.s_f.start + kf] = float(np.clip(net, -hold, hold))
            jb = pl.index(ax, "bearing")
            free = [k for k in sl if k != jb]
            dq = np.zeros(pl.n)
            dq[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
            q[sl] += dq[sl] - dq[self.ie[ax]]
        x[self.s_q] = q
        return x

    def _window_inputs(self, seq: Sequence, t0: float, t1: float, n_avg: int):
        """Rest state before the window (mean encoder positions and currents of the ``n_avg`` samples before
        it, direction of the last move of each axis for the hysteron), samples of the window, filtered currents."""
        k0 = int(np.searchsorted(seq.t, t0))
        pre = slice(max(0, k0 - n_avg), max(1, k0))
        p0 = {ax: float(np.mean(seq.enc[ax][pre])) for ax in AXES}
        i0 = {ax: float(np.mean(seq.current[ax][pre])) for ax in AXES}
        s_h = {}
        for ax in AXES:
            prev = [(a, b) for a, b in seq.moves(ax) if b <= t0]
            if prev:
                a, b = prev[-1]
                s_h[ax] = float(np.sign(seq.cmd[ax][seq.t <= b][-1] - seq.cmd[ax][seq.t <= a][-1]))
        m = seq.window(t0, t1)
        tm = seq.t[m] - t0
        cur = {ax: seq.current[ax][m] for ax in AXES}
        if self.input_lowpass_hz:
            from scipy.signal import butter, filtfilt

            b, a = butter(2, self.input_lowpass_hz * 2 * float(np.median(np.diff(seq.t))))
            cur = {ax: filtfilt(b, a, c) for ax, c in cur.items()}
        return p0, i0, s_h, m, tm, cur

    def simulate_window(self, seq: Sequence, t0: float, t1: float, *, n_avg: int = 200, rtol: float = 1e-4,
                        atol_scale: float = 100.0, max_cpu_s: float | None = 60.0, decimate: int = 4) -> dict:
        """Encoder positions (change from the start, both axes) simulated with the measured currents of the
        window, against the measured ones, on the measurement grid (every ``decimate`` samples)."""
        p0, i0, s_h, m, tm, cur = self._window_inputs(seq, t0, t1, n_avg)
        self.plant.theta0 = p0["elevation"]
        x0 = self.initial_state(p0, i0, s_h)
        current = {ax: (lambda t, ii=cur[ax]: float(np.interp(t, tm, ii))) for ax in AXES}
        atol = np.full(self.n_states, 1e-10)
        atol[self.s_dq] = 1e-8
        atol[self.s_f] = 1e-7
        for ax in AXES:
            atol[self.i_h[ax]] = 1e-6
        atol *= atol_scale
        typ = np.full(self.n_states, 1e-3)
        typ[self.s_dq] = 1e-2
        typ[self.s_f] = 1e-1
        for ax in AXES:
            typ[self.i_h[ax]] = 1e-1
        t_cpu = time.process_time()
        count = [0]

        def f(t, x):
            count[0] += 1
            if max_cpu_s is not None and count[0] % 200 == 0 and time.process_time() - t_cpu > max_cpu_s:
                raise TimeoutError
            return self.rhs(t, x, current)

        def jac(t, x):
            f0 = f(t, x)
            h = 1e-6 * np.maximum(np.abs(x), typ)
            J = np.empty((x.size, x.size))
            for j in range(x.size):
                xp = x.copy()
                xp[j] += h[j]
                J[:, j] = (f(t, xp) - f0) / h[j]
            return J

        te = tm[::decimate]
        ok, msg = True, ""
        try:
            sol = solve_ivp(f, (0.0, te[-1]), x0, method="Radau", t_eval=te, rtol=rtol, atol=atol, jac=jac)
            ok, msg = sol.success, sol.message
            y, ts = sol.y, sol.t
        except TimeoutError:
            ok, msg, y, ts = False, f"CPU time above {max_cpu_s} s", np.zeros((self.n_states, 1)), te[:1]
        out = {"t": te, "ok": ok, "message": msg, "cpu": time.process_time() - t_cpu, "n_rhs": count[0]}
        for ax in AXES:
            enc_meas = (seq.enc[ax][m] - seq.enc[ax][m][0])[::decimate]
            # change from the first measured sample of the window (the model starts at the mean before it)
            sim = np.full(te.size, np.nan)
            n_ok = ts.size if not ok else te.size
            sim[:n_ok] = (y[self.ie[ax]] - seq.enc[ax][m][0])[:n_ok]
            out[ax] = {"enc_meas": enc_meas, "enc": sim, "current_meas": seq.current[ax][m][::decimate],
                       "T_f": y[self.s_f.start + [el.axis for el in self.plant.friction].index(ax)] if ok else None}
        return out

    def simulate_window_rigid(self, seq: Sequence, t0: float, t1: float, *, n_avg: int = 200, rtol: float = 1e-5,
                              max_cpu_s: float | None = 60.0, decimate: int = 4) -> dict:
        r"""As :meth:`simulate_window`, with the plant constrained to the rigid motion of each axis (all its
        nodes turning together, :math:`\mathbf q = \mathbf r\,\theta`): inertia :math:`J = \mathbf r^T\mathbf M
        \mathbf r` at the start of the window, ground damping :math:`\mathbf r^T\mathbf C\mathbf r`, the same
        motor, bearing friction and extra torques; the two axes are uncoupled. A fast surrogate of the position
        on the time scale of the moves (the flexible modes are dropped); balanced telescope only."""
        pl = self.plant
        if pl.value("g") != 0.0:
            raise NotImplementedError("rigid surrogate without gravity only")
        p0, i0, s_h, m, tm, cur = self._window_inputs(seq, t0, t1, n_avg)
        pl.theta0 = p0["elevation"]
        q0 = np.zeros(pl.n)
        n_psi = pl.sim.struttura.n_psi
        q0[:n_psi], q0[n_psi:] = p0["azimuth"], p0["elevation"]
        M, C = pl.M(q0), pl.linear_matrices()[1]
        ax_par = []
        x0 = np.zeros(4 * len(AXES))  # theta, omega, friction torque, hysteron state per axis
        for k, ax in enumerate(AXES):
            r = np.zeros(pl.n)
            r[:n_psi] = ax == "azimuth"
            r[n_psi:] = ax == "elevation"
            el, mot = self.friction(ax), self.motors[ax]
            ax_par.append((ax, float(r @ M @ r), float(r @ C @ r) + mot.p.b_viscous, el, mot.p.Kt * self.torque_scale[ax],
                           self.hysteron.get(ax, (0.0, 1.0))[1]))
            s = float(s_h.get(ax, 0.0)) if ax in self.hysteron else 0.0
            net = self.torque_scale[ax] * mot.p.Kt * i0[ax] + self.axis_extra(ax, p0[ax], s)[0]
            hold = 0.99 * el.level(0.0)
            x0[4 * k: 4 * k + 4] = p0[ax], 0.0, float(np.clip(net, -hold, hold)), s
        i_of = [(lambda t, ii=cur[ax]: float(np.interp(t, tm, ii))) for ax in AXES]
        t_cpu = time.process_time()
        count = [0]

        def f(t, x):
            count[0] += 1
            if max_cpu_s is not None and count[0] % 500 == 0 and time.process_time() - t_cpu > max_cpu_s:
                raise TimeoutError
            dx = np.empty_like(x)
            for k, (ax, J, c, el, kt, delta) in enumerate(ax_par):
                th, w, tf, s = x[4 * k: 4 * k + 4]
                av = np.sqrt(w * w + el.v_eps**2) - el.v_eps
                dx[4 * k] = w
                dx[4 * k + 1] = (kt * i_of[k](t) - c * w - tf + self.axis_extra(ax, th, s)[0]) / J
                dx[4 * k + 2] = el.sigma0 * (w - av * tf / el.level(w))
                dx[4 * k + 3] = (w - (np.sqrt(w * w + 1e-12) - 1e-6) * s) / delta if ax in self.hysteron else 0.0
            return dx

        def jac(t, x):
            Jm = np.zeros((x.size, x.size))
            for k, (ax, J, c, el, kt, delta) in enumerate(ax_par):
                th, w, tf, s = x[4 * k: 4 * k + 4]
                i = 4 * k
                _, dth, ds = self.axis_extra(ax, th, s)
                Jm[i, i + 1] = 1.0
                Jm[i + 1, i: i + 4] = dth / J, -c / J, -1.0 / J, ds / J
                rt = np.sqrt(w * w + el.v_eps**2)
                av, dav = rt - el.v_eps, w / rt
                g = el.level(w)
                dg = 0.0 if el.Ts is None else (el.Ts - el.Tc) * np.exp(-(w / el.vs) ** 2) * (-2.0 * w / el.vs**2)
                Jm[i + 2, i + 1] = el.sigma0 * (1.0 - tf * (dav * g - av * dg) / g**2)
                Jm[i + 2, i + 2] = -el.sigma0 * av / g
                if ax in self.hysteron:
                    rh = np.sqrt(w * w + 1e-12)
                    Jm[i + 3, i + 1] = (1.0 - w / rh * s) / delta
                    Jm[i + 3, i + 3] = -(rh - 1e-6) / delta
            return Jm

        atol = np.tile([1e-8, 1e-6, 1e-5, 1e-4], len(AXES))
        te = tm[::decimate]
        ok, msg = True, ""
        try:
            sol = solve_ivp(f, (0.0, te[-1]), x0, method="Radau", t_eval=te, rtol=rtol, atol=atol, jac=jac)
            ok, msg = sol.success, sol.message
            y, ts = sol.y, sol.t
        except TimeoutError:
            ok, msg, y, ts = False, f"CPU time above {max_cpu_s} s", np.zeros((x0.size, 1)), te[:1]
        out = {"t": te, "ok": ok, "message": msg, "cpu": time.process_time() - t_cpu, "n_rhs": count[0]}
        for k, ax in enumerate(AXES):
            enc_meas = (seq.enc[ax][m] - seq.enc[ax][m][0])[::decimate]
            sim = np.full(te.size, np.nan)
            n_ok = ts.size if not ok else te.size
            sim[:n_ok] = (y[4 * k] - seq.enc[ax][m][0])[:n_ok]
            out[ax] = {"enc_meas": enc_meas, "enc": sim, "current_meas": seq.current[ax][m][::decimate],
                       "T_f": y[4 * k + 2] if ok else None}
        return out
