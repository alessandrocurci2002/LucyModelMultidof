r"""Drives on an EtherCAT network: the master exchanges the process data (PDO) with the drives once per cycle.

Sampled-data closed loop, cycle :math:`T` (1 ms):

* **position targets** (RxPDO, cyclic synchronous position): the master sends each cycle the target for the end
  of the cycle; the drive interpolates linearly between consecutive targets (``interpolation = "linear"``) or holds
  it (``"zoh"``). The velocity feedforward of the drive is the slope of the interpolated target;
* **actual currents and positions** (TxPDO): sampled by the drives at the cycle start and read by the master;
* **rate gyro** on the telescope (:class:`lucy_cl.gyro.MemsGyro`), read by the master over SPI: second-order
  sensor dynamics (states of the continuous model), latency as a first-order Pade delay, noise added at the read;
* **disturbance observer** (Hilkert & Pautler, 2011) in the master, discretised exactly with the inputs held over
  the cycle; its feedforward current, through an optional discrete notch (bilinear, prewarped), reaches the drive
  in the RxPDO of the next cycle (``ff_delay_cycles``) and is held for one cycle.

The drives' own loops (position, velocity, current) stay continuous. Between two cycles the plant, motors, drives
and gyro dynamics are integrated with one Radau solver kept across the cycles (step size and Jacobian are reused).
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Callable, Mapping

import numpy as np
from scipy.integrate import Radau
from scipy.signal import bilinear

from . import AXES
from .closed_loop import ClosedLoop
from .disturbance_observer import ObserverParams
from .driver import Biquad, Drive
from .gyro import MemsGyro
from .motor import Motor
from .plant import GimbalPlant


@dataclass(frozen=True)
class EtherCATParams:
    cycle_s: float = 1e-3
    interpolation: str = "linear"
    """Position target between two cycles: ``linear`` (interpolated by the drive) or ``zoh``."""
    ff_delay_cycles: int = 1
    """Cycles between the master's reading and the drive applying the feedforward current."""


class DiscreteBiquad:
    """Analog biquad by the bilinear transform prewarped at its pole frequency; direct form II transposed."""

    def __init__(self, bq: Biquad, fs: float):
        num, den = bq.coefficients()
        w0 = 2 * np.pi * bq.fp
        fw = w0 / (2 * np.tan(w0 / (2 * fs)))
        self.b, self.a = bilinear(num, den, fs=fw)
        self.b, self.a = self.b / self.a[0], self.a / self.a[0]
        self.z = np.zeros(2)

    def __call__(self, u: float) -> float:
        b, a, z = self.b, self.a, self.z
        y = b[0] * u + z[0]
        z[0] = b[1] * u - a[1] * y + z[1]
        z[1] = b[2] * u - a[2] * y
        return float(y)


class SampledObserverLoop(ClosedLoop):
    """Closed loop with the EtherCAT exchange and, for the axes in ``observers``, the gyro and the observer in
    the master. Extra continuous states per observed axis: the two of the gyro dynamics and the latency."""

    def __init__(self, plant: GimbalPlant, motors: Mapping[str, Motor], drives: Mapping[str, Drive],
                 observers: Mapping[str, ObserverParams], gyro: MemsGyro | None, ecat: EtherCATParams):
        super().__init__(plant, motors, drives)
        self.observers, self.gyro, self.ecat = dict(observers), gyro, ecat
        if self.observers and gyro is None:
            raise ValueError("an observer needs the gyro")
        k = self.n_states
        self.s_gyro = {}
        for ax in self.observers:
            self.s_gyro[ax] = slice(k, k + 3)  # sensor output, its derivative / wn, latency state
            k += 3
        self.n_states = k
        self.jl = plant.index("azimuth", "bearing")
        self.jp = plant.index("elevation", "telescope")
        self._held: dict = {"t0": 0.0, "p0": {ax: 0.0 for ax in AXES}, "slope": {ax: 0.0 for ax in AXES},
                            "v": {ax: 0.0 for ax in AXES}, "i_ff": {ax: 0.0 for ax in AXES}}
        """Inputs held over the current cycle: target (start, slope), velocity feedforward, feedforward current."""

    # ---- continuous part ---------------------------------------------------------------------------

    def true_rate(self, ax: str, q: np.ndarray, dq: np.ndarray) -> float:
        """Rate at the gyro axis: elevation, cross-elevation."""
        return float(dq[self.jp]) if ax == "elevation" else float(dq[self.jl] * np.cos(q[self.jp]))

    def gyro_reading(self, ax: str, x: np.ndarray) -> float:
        """Noise-free reading of the gyro: sensor output through the first-order Pade latency."""
        y, _, d = x[self.s_gyro[ax]]
        return float(-y + 2 * d) if self.gyro.latency_s > 0 else float(y)

    def evaluate(self, t: float, x: np.ndarray, cmd=None, i_exc=None) -> tuple[np.ndarray, dict]:
        h = self._held
        cmd = cmd or {ax: (lambda t, ax=ax: (h["p0"][ax] + h["slope"][ax] * (t - h["t0"]), h["v"][ax])) for ax in AXES}
        exc = i_exc or {ax: (lambda t, v=h["i_ff"][ax]: v) for ax in AXES}
        dx, sig = super().evaluate(t, x, cmd, exc)
        if self.observers:
            q, dq = x[self.s_q], x[self.s_dq]
            g = self.gyro
            wn, z = g.wn, g.damping
            for ax in self.observers:
                y, yd, d = x[self.s_gyro[ax]]
                sl = self.s_gyro[ax]
                dx[sl.start] = wn * yd
                dx[sl.start + 1] = wn * (self.true_rate(ax, q, dq) - y) - 2 * z * wn * yd
                dx[sl.start + 2] = 2 / g.latency_s * (y - d) if g.latency_s > 0 else 0.0
        return dx, sig

    def jacobian_steps(self) -> np.ndarray:
        typ = super().jacobian_steps()
        for sl in self.s_gyro.values():
            typ[sl] = 1e-2
        return typ

    def atol(self) -> np.ndarray:
        a = super().atol()
        for sl in self.s_gyro.values():
            a[sl] = 1e-8
        return a

    # ---- sampled simulation ------------------------------------------------------------------------

    def simulate_sampled(self, refs: Mapping[str, Callable[[float], tuple[float, float]]], t_end: float,
                         x0: np.ndarray, *, noise: Mapping[str, np.ndarray] | None = None, rtol: float = 1e-4,
                         atol_scale: float = 100.0) -> dict:
        """Integrates over ``t_end`` with one master cycle per step of ``ecat.cycle_s``; ``refs`` give the
        target positions relative to the start; ``noise``: gyro noise of each observed axis at the cycle
        times. Returns the states and the master signals at the cycle times."""
        T = self.ecat.cycle_s
        n = int(round(t_end / T))
        tk = np.arange(n + 1) * T
        noise = noise or {}
        typ = self.jacobian_steps()
        count = [0]

        def f(t, x):
            count[0] += 1
            return self.evaluate(t, x)[0]

        def jac(t, x):
            f0 = f(t, x)
            hh = 1e-6 * np.maximum(np.abs(x), typ)
            J = np.empty((x.size, x.size))
            for j in range(x.size):
                xp = x.copy()
                xp[j] += hh[j]
                J[:, j] = (f(t, xp) - f0) / hh[j]
            return J

        obs = self.observers
        notch = {ax: [DiscreteBiquad(b, 1.0 / T) for b in o.ff_filters if b.active] for ax, o in obs.items()}
        w_hat = {ax: 0.0 for ax in obs}
        queue = {ax: [0.0] * max(1, self.ecat.ff_delay_cycles) for ax in AXES}
        log = {k: np.zeros((n + 1, len(AXES))) for k in ("i_cmd", "i", "i_ff", "Td_hat", "w_gyro", "p_cmd")}
        X = np.empty((self.n_states, n + 1))
        x = np.asarray(x0, dtype=float).copy()
        p_prev = {ax: refs[ax](0.0)[0] for ax in AXES}
        solver = None
        t_cpu = time.process_time()
        ok, msg = True, ""
        for k in range(n + 1):
            t = tk[k]
            X[:, k] = x
            q = x[self.s_q]
            # master: observer on the gyro reading and on the actual current of the TxPDO
            for j, ax in enumerate(AXES):
                i_ff_new = 0.0
                if ax in obs:
                    o = obs[ax]
                    w_g = self.gyro_reading(ax, x) + (noise[ax][k] if ax in noise else 0.0)
                    if ax == "azimuth":  # rate about the azimuth axis from the cross-elevation gyro
                        w_g /= np.cos(q[self.ie["elevation"]])
                    i_m = x[self.i_cur[ax]]
                    Td = o.Ko * (w_g - w_hat[ax])
                    a = o.Ko / o.J
                    e = np.exp(-a * T)
                    w_hat[ax] = e * w_hat[ax] + (1 - e) * (w_g + o.Kt * i_m / (o.J * a))
                    u = -Td / o.Kt
                    for nb in notch[ax]:
                        u = nb(u)
                    i_ff_new = u
                    log["Td_hat"][k, j], log["w_gyro"][k, j] = Td, w_g
                queue[ax].append(i_ff_new)
            if k == n:
                break
            # targets of this cycle and the feedforward current applied by the drives
            h = {"t0": t, "p0": {}, "slope": {}, "v": {}, "i_ff": {}}
            for j, ax in enumerate(AXES):
                p_next = refs[ax](tk[k + 1])[0]
                p_now = refs[ax](t)[0]
                if self.ecat.interpolation == "linear":
                    h["p0"][ax], h["slope"][ax], h["v"][ax] = p_now, (p_next - p_now) / T, (p_next - p_now) / T
                else:
                    h["p0"][ax], h["slope"][ax], h["v"][ax] = p_now, 0.0, (p_now - p_prev[ax]) / T
                p_prev[ax] = p_now
                h["i_ff"][ax] = queue[ax].pop(0) if self.ecat.ff_delay_cycles > 0 else queue[ax].pop()
            self._held = h
            _, sig = self.evaluate(t, x)
            for j, ax in enumerate(AXES):
                log["p_cmd"][k, j], log["i_cmd"][k, j], log["i"][k, j] = sig[ax][0], sig[ax][1], sig[ax][3]
                log["i_ff"][k, j] = h["i_ff"][ax]
            # plant, motors, drives and gyro over the cycle
            if solver is None:
                solver = Radau(f, t, x, tk[-1], rtol=rtol, atol=self.atol() * atol_scale, jac=jac, first_step=T / 4)
            else:  # the inputs changed at the cycle start: derivative used by the error estimate of the first step
                solver.f = solver.fun(t, x)
            solver.t_bound = tk[k + 1]
            solver.status = "running"
            while solver.status == "running":
                m = solver.step()
                if solver.status == "failed":
                    ok, msg = False, f"t = {solver.t:.4f} s: {m}"
                    break
            if not ok:
                X, tk = X[:, :k + 1], tk[:k + 1]
                log = {kk: v[:k + 1] for kk, v in log.items()}
                break
            x = solver.y.copy()
        return {"t": tk, "x": X, "log": log, "success": ok, "message": msg or "completed", "n_rhs": count[0],
                "cpu_s": time.process_time() - t_cpu}
