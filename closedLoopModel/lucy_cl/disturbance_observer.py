r"""Rate gyro on the telescope and single-state disturbance observer fed forward to the motor current.

The observer is the one of Hilkert and Pautler (*A reduced-order disturbance observer applied to inertially
stabilized line-of-sight control*, Proc. SPIE 8052, 2011). For each axis a rate gyro on the telescope measures
the inertial rate :math:`\omega_g` of the line of sight; a one-state model of the rigid axis, driven by the motor
torque measured through the current, reconstructs the disturbance torque:

.. math::
    \dot{\hat\omega}=\frac{K_t'\,i_m+\hat T_d}{J'},\qquad \hat T_d=K_o\,(\omega_g-\hat\omega)
    \qquad\Longleftrightarrow\qquad
    \hat T_d=\frac{K_o/J'}{s+K_o/J'}\bigl(J's\,\omega_g-K_t'\,i_m\bigr),

with :math:`J'`, :math:`K_t'` the inertia and torque constant known to the observer and :math:`K_o=2\pi f_oJ'`
(:math:`f_o` the bandwidth of the estimate). The estimate is subtracted at the input of the current loop,
:math:`i_{ff}=-H_n\hat T_d/K_t'`, where :math:`H_n` is an optional chain of biquads (a notch on the structural
mode, as Hilkert's figure 3 passes the estimate through the notch filter of the loop); the position and velocity
loops of the drive, closed on the encoders, are unchanged.

Gyro axes on the telescope: elevation, the rate of the telescope node :math:`\dot\theta_p`; cross-elevation,
:math:`\dot\psi_L\cos\theta_p` (fork node), which the azimuth observer divides by the cosine of the elevation of
its encoder. The gyro is ideal (no noise, no bias) or a first-order lag.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping

import numpy as np

from . import AXES
from .closed_loop import ClosedLoop
from .driver import Biquad, Drive, _Chain
from .motor import Motor
from .plant import GimbalPlant


@dataclass(frozen=True)
class ObserverParams:
    J: float
    """Inertia of the axis known to the observer, kg m^2."""
    Kt: float
    """Torque constant known to the observer, N m/A."""
    bandwidth_hz: float
    """Bandwidth of the disturbance estimate, :math:`K_o/(2\\pi J')`."""
    gyro_bandwidth_hz: float | None = None
    """First-order bandwidth of the gyro; ``None``: ideal."""
    ff_filters: tuple[Biquad, ...] = ()
    """Biquads on the feedforward path, e.g. a notch on the structural mode."""

    @property
    def Ko(self) -> float:
        return 2 * np.pi * self.bandwidth_hz * self.J


def rigid_inertia(plant: GimbalPlant, axis: str, q: np.ndarray | None = None) -> float:
    """Inertia of the rigid motion of ``axis`` (all its nodes turning together), :math:`\\mathbf r^T\\mathbf M\\mathbf r`."""
    q = plant.q_equilibrium() if q is None else q
    n_psi = plant.sim.struttura.n_psi
    r = np.zeros(plant.n)
    if axis == "azimuth":
        r[:n_psi] = 1.0
    else:
        r[n_psi:] = 1.0
    return float(r @ plant.M(q) @ r)


class GyroObserverLoop(ClosedLoop):
    """Closed loop on the encoders plus, for the axes in ``observers``, the gyro and the disturbance observer.
    Extra states: :math:`\\hat\\omega` of each observer, the gyro lag and the feedforward filters when there are."""

    def __init__(self, plant: GimbalPlant, motors: Mapping[str, Motor], drives: Mapping[str, Drive],
                 observers: Mapping[str, ObserverParams]):
        super().__init__(plant, motors, drives)
        self.observers = dict(observers)
        k = self.n_states
        self.i_obs, self.i_gyro, self.s_ff, self._ff = {}, {}, {}, {}
        for ax, o in self.observers.items():
            self.i_obs[ax] = k
            k += 1
            if o.gyro_bandwidth_hz:
                self.i_gyro[ax] = k
                k += 1
            if o.ff_filters:
                self._ff[ax] = _Chain([b.state_space() for b in o.ff_filters if b.active])
                self.s_ff[ax] = slice(k, k + self._ff[ax].n)
                k += self._ff[ax].n
        self.n_states = k
        self.jl = plant.index("azimuth", "bearing")
        self.jp = plant.index("elevation", "telescope")
        self.signal_names = self.signal_names + ("T_d_hat", "w_gyro")

    def gyro_rate(self, ax: str, q: np.ndarray, dq: np.ndarray) -> float:
        """Rate measured by an ideal gyro on the telescope: elevation, cross-elevation."""
        return float(dq[self.jp]) if ax == "elevation" else float(dq[self.jl] * np.cos(q[self.jp]))

    def evaluate(self, t: float, x: np.ndarray, cmd: Mapping[str, Callable[[float], tuple[float, float]]],
                 i_exc: Mapping[str, Callable[[float], float]] | None = None) -> tuple[np.ndarray, dict]:
        q, dq = x[self.s_q], x[self.s_dq]
        est = {}
        exc = dict(i_exc or {})
        for ax, o in self.observers.items():
            w_g = x[self.i_gyro[ax]] if ax in self.i_gyro else self.gyro_rate(ax, q, dq)
            if ax == "azimuth":  # rate about the azimuth axis from the cross-elevation gyro
                w_g /= np.cos(q[self.ie["elevation"]])
            Td = o.Ko * (w_g - x[self.i_obs[ax]])
            i_ff, d_ff = -Td / o.Kt, None
            if ax in self._ff:
                i_ff, d_ff = self._ff[ax].step(x[self.s_ff[ax]], i_ff)
            est[ax] = (Td, w_g, d_ff)
            prev = exc.get(ax)
            exc[ax] = (lambda t, v=i_ff, f=prev: v + (f(t) if f else 0.0))
        dx, sig = super().evaluate(t, x, cmd, exc)
        for ax in AXES:
            if ax in self.observers:
                o = self.observers[ax]
                Td, w_g, d_ff = est[ax]
                dx[self.i_obs[ax]] = (o.Kt * x[self.i_cur[ax]] + Td) / o.J
                if d_ff is not None:
                    dx[self.s_ff[ax]] = d_ff
                if ax in self.i_gyro:
                    dx[self.i_gyro[ax]] = 2 * np.pi * o.gyro_bandwidth_hz * (self.gyro_rate(ax, q, dq) - x[self.i_gyro[ax]])
                sig[ax] = sig[ax] + (Td, w_g)
            else:
                sig[ax] = sig[ax] + (0.0, self.gyro_rate(ax, q, dq))
        return dx, sig

    def jacobian_steps(self) -> np.ndarray:
        typ = super().jacobian_steps()
        for i in list(self.i_obs.values()) + list(self.i_gyro.values()):
            typ[i] = 1e-2
        for sl in self.s_ff.values():
            typ[sl] = 1e-1
        return typ

    def atol(self) -> np.ndarray:
        a = super().atol()
        for i in list(self.i_obs.values()) + list(self.i_gyro.values()):
            a[i] = 1e-8
        for sl in self.s_ff.values():
            a[sl] = 1e-6
        return a
