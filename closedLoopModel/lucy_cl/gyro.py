r"""MEMS rate gyro: second-order sensor dynamics, Allan-variance noise, digital output read over SPI.

The output of one axis is

.. math::
    \omega_g(t_k)=\Bigl[H(s)\bigl(\omega+n_w+b+r\bigr)\Bigr]_{t=\lfloor (t_k-\tau_l)f_{odr}\rfloor/f_{odr}},\qquad
    H(s)=\frac{\omega_n^2}{s^2+2\zeta\omega_ns+\omega_n^2},

that is the true rate :math:`\omega` plus the noise, through the sensor dynamics (bandwidth :math:`\omega_n/2\pi`),
held at the output data rate :math:`f_{odr}` and read by the controller :math:`\tau_l` after it is sampled (data
ready and SPI burst read). The noise terms come from the Allan deviation:

* angle random walk :math:`N`: white rate noise, :math:`\sigma_A(\tau)=N/\sqrt\tau`;
* bias instability :math:`B` (minimum of the Allan deviation): a first-order Gauss-Markov bias whose Allan
  deviation peaks at :math:`B` at the time of the minimum, :math:`\tau_B\approx1.89\,T_c`, with standard deviation
  :math:`B/0.437`;
* rate random walk :math:`K`: random walk of the bias, :math:`\sigma_A(\tau)=K\sqrt{\tau/3}`;
* a constant turn-on bias.

The default values are those of the Analog Devices ADIS16507-1 (x and y axes): 480 Hz, 2000 SPS, ARW 0.13 deg/sqrt(h),
in-run bias stability 2.3 deg/h; the data sheet gives neither the time of the Allan minimum nor a rate random walk.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import bilinear, lfilter

DEG = np.pi / 180.0


@dataclass(frozen=True)
class MemsGyro:
    name: str = "ADIS16507-1, x and y axes"
    bandwidth_hz: float = 480.0
    """-3 dB bandwidth of the sensor (natural frequency with ``damping`` = 1/sqrt 2)."""
    damping: float = 0.7071
    odr_hz: float = 2000.0
    """Output data rate of the registers read over SPI."""
    latency_s: float = 0.5e-3
    """Age of the sample when the controller uses it (one output period: data ready, SPI burst read)."""
    arw_deg_rthr: float = 0.13
    """Angle random walk, deg/sqrt(h)."""
    bias_instability_deg_hr: float = 2.3
    """In-run bias stability (minimum of the Allan deviation), deg/h."""
    bias_tau_s: float = 50.0
    """Time of the minimum of the Allan deviation, s (assumed: not in the data sheet)."""
    rrw_deg_hr_rthr: float = 0.0
    """Rate random walk, deg/h/sqrt(h)."""
    turn_on_bias_deg_s: float = 0.0
    noise_scale: float = 1.0
    """Multiplies all the noise terms (sensitivity studies); 0: noise-free."""

    # ---- noise coefficients in SI units ----------------------------------------------------------------

    @property
    def N(self) -> float:
        """Angle random walk, rad/sqrt(s) (white rate noise density, rad/s/sqrt(Hz))."""
        return self.noise_scale * self.arw_deg_rthr * DEG / 60.0

    @property
    def B(self) -> float:
        """Bias instability, rad/s."""
        return self.noise_scale * self.bias_instability_deg_hr * DEG / 3600.0

    @property
    def K(self) -> float:
        """Rate random walk, rad/s/sqrt(s)."""
        return self.noise_scale * self.rrw_deg_hr_rthr * DEG / 3600.0 / 60.0

    @property
    def Tc(self) -> float:
        """Correlation time of the Gauss-Markov bias, s."""
        return self.bias_tau_s / 1.89

    @property
    def sigma_gm(self) -> float:
        return self.B / 0.4373

    # ---- sensor dynamics -------------------------------------------------------------------------------

    @property
    def wn(self) -> float:
        return 2 * np.pi * self.bandwidth_hz

    def response(self, f: np.ndarray) -> np.ndarray:
        s = 2j * np.pi * np.asarray(f)
        return self.wn**2 / (s**2 + 2 * self.damping * self.wn * s + self.wn**2)

    def discrete_filter(self, fs: float) -> tuple[np.ndarray, np.ndarray]:
        """H(s) by the bilinear transform prewarped at the natural frequency."""
        fw = self.wn / (2 * np.tan(self.wn / (2 * fs)))
        return bilinear([self.wn**2], [1.0, 2 * self.damping * self.wn, self.wn**2], fs=fw)

    # ---- noise -----------------------------------------------------------------------------------------

    def raw_noise(self, n: int, fs: float, rng: np.random.Generator) -> np.ndarray:
        """Noise at the sensor input on a uniform grid of rate ``fs`` (white + Gauss-Markov + random walk + bias)."""
        x = self.N * np.sqrt(fs) * rng.standard_normal(n)
        if self.B > 0:
            a = np.exp(-1.0 / (fs * self.Tc))
            w = self.sigma_gm * np.sqrt(1 - a * a) * rng.standard_normal(n)
            w[0] = self.sigma_gm * rng.standard_normal()
            x += lfilter([1.0], [1.0, -a], w)
        if self.K > 0:
            x += np.cumsum(self.K / np.sqrt(fs) * rng.standard_normal(n))
        return x + self.noise_scale * self.turn_on_bias_deg_s * DEG

    def output_noise(self, t_read: np.ndarray, rng: np.random.Generator, fs: float = 16000.0) -> np.ndarray:
        """Noise part of the readings at the times ``t_read``: through the sensor dynamics, held at the output
        data rate and read with the latency."""
        t_end = float(t_read[-1]) + 1.0
        n = int(np.ceil(t_end * fs)) + 1
        b, a = self.discrete_filter(fs)
        shaped = lfilter(b, a, self.raw_noise(n, fs, rng))
        t_sample = np.floor((np.asarray(t_read) - self.latency_s) * self.odr_hz) / self.odr_hz
        k = np.clip(np.round(t_sample * fs).astype(int), 0, n - 1)
        return np.where(t_sample >= 0, shaped[k], 0.0)

    def odr_noise(self, duration_s: float, rng: np.random.Generator, fs: float = 8000.0) -> np.ndarray:
        """Noise of the output registers at the output data rate over ``duration_s`` (Allan check)."""
        n = int(np.ceil(duration_s * fs))
        b, a = self.discrete_filter(fs)
        shaped = lfilter(b, a, self.raw_noise(n, fs, rng))
        step = int(round(fs / self.odr_hz))
        return shaped[::step]

    def allan_model(self, tau: np.ndarray) -> dict[str, np.ndarray]:
        """Allan deviation of the noise terms (rad/s) at the averaging times ``tau``."""
        tau = np.asarray(tau, dtype=float)
        out = {"angle random walk": self.N / np.sqrt(tau)}
        if self.B > 0:
            q = 2 * self.sigma_gm**2 / self.Tc
            r = tau / self.Tc
            out["bias instability (Gauss-Markov)"] = np.sqrt(q * self.Tc**2 / tau * (1 - (3 - 4 * np.exp(-r) + np.exp(-2 * r)) / (2 * r)))
        if self.K > 0:
            out["rate random walk"] = self.K * np.sqrt(tau / 3)
        out["total"] = np.sqrt(sum(v**2 for v in out.values()))
        return out


def allan_deviation(x: np.ndarray, fs: float, taus: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Overlapping Allan deviation of the rate samples ``x`` (rate ``fs``) at the averaging times ``taus``."""
    theta = np.concatenate([[0.0], np.cumsum(x)]) / fs
    out_t, out_s = [], []
    for tau in taus:
        m = int(round(tau * fs))
        if m < 1 or 2 * m >= theta.size:
            continue
        d = theta[2 * m:] - 2 * theta[m:-m] + theta[:-2 * m]
        out_t.append(m / fs)
        out_s.append(np.sqrt(np.mean(d**2) / (2 * (m / fs) ** 2)))
    return np.array(out_t), np.array(out_s)
