"""Disturbances: torque on the bearing port of each axis and noise on the position measurement of the master.

A scenario is a :class:`Disturbance`; the generic disturbance torque of the thesis is the zero-mean band-limited
Gaussian torque (``gauss``, ``band``), the model error the scaled Coulomb level of the plant (``friction_factor``).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.signal import butter, sosfiltfilt

LSB = 2 * np.pi / 2 ** 26
"""Resolution of the encoder (26 bit), rad."""


@dataclass(frozen=True)
class Disturbance:
    """All fields optional (zero / absent = off).

    ``gauss``: rms of the zero-mean Gaussian torque [N m], filtered in ``band`` [Hz] (4th-order Butterworth, zero phase);
    ``mean``: static part of the torque (cables, unbalance) [N m], present from the start; ``step``: torque step [N m] at
    ``step_time``; ``sine``: (amplitude [N m], frequency [Hz]) from ``start_time``; ``friction_factor``: Coulomb level of
    the plant multiplied (controllers unchanged); ``noise``: rms of the white noise on the position measurement [rad],
    quantised to the encoder resolution."""

    gauss: float = 0.0
    band: tuple[float, float] = (0.1, 50.0)
    mean: float = 0.0
    step: float = 0.0
    step_time: float = 15.0
    sine: tuple[float, float] | None = None
    start_time: float = 5.0
    friction_factor: float = 1.0
    noise: float = 0.0

    def to_dict(self) -> dict:
        return {k: (list(v) if isinstance(v, tuple) else v) for k, v in asdict(self).items()}


def disturbance_torque(d: Disturbance, t: np.ndarray, seed: int, dt: float) -> np.ndarray:
    """Torque on the bearing port at the times ``t`` (step ``dt``); ``seed`` sets the Gaussian realisation."""
    out = np.zeros_like(t)
    if d.step:
        out += d.step * (t >= d.step_time)
    if d.sine is not None:
        out += d.sine[0] * np.sin(2 * np.pi * d.sine[1] * t) * (t >= d.start_time)
    if d.gauss:
        # margin discarded: the 0.1 Hz filter leaves transients at the edges for tens of seconds
        pad = int(100.0 / dt)
        n = np.random.default_rng(seed).standard_normal(len(t) + 2 * pad)
        n = sosfiltfilt(butter(4, d.band, "bandpass", fs=1 / dt, output="sos"), n)[pad:pad + len(t)]
        out += d.gauss * (n - n.mean()) / n.std()
    if d.mean:
        out += d.mean
    return out


class NoisyMeasurement:
    """White noise and quantisation on the measurement that reaches the estimator of the controller."""

    def __init__(self, c, sigma: float, rng: np.random.Generator):
        self.c, self.sigma, self.rng = c, sigma, rng

    def initial_state(self):
        return self.c.initial_state()

    def output(self, st, r):
        return self.c.output(st, r)

    def __getattr__(self, name):
        # ``measure`` and the other methods of the wrapped controller
        return getattr(self.c, name)

    def update(self, st, y, u):
        self.c.update(st, float(np.round((y + self.sigma * self.rng.standard_normal()) / LSB) * LSB), u)
