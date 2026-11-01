r"""Rational transfer functions of the AKD2G controller of one axis.

They are built from the drive parameter file (.akd) and the structure of the loops of the AKD User Manual
(Position Loop, Velocity Loop, Current Loop pages; see ``sim/kollmorgen_control_report``), and checked on the
traces that WorkBench stores in the Bode export of the same test. Polynomials are in :math:`s`, in descending
powers (numpy convention). The blocks are:

.. math::
    C_p = K_{pp},\qquad C_v = K_{pv}\,\frac{s + 2\pi K_{iv}}{s},\qquad H_{fwd} = H_{AR1}H_{AR2},\qquad
    H_{fb} = H_{AR3}H_{AR4}\,\frac{\omega_o}{s + \omega_o},\qquad C_i = K_{pi}\,\frac{s + R/L}{s}.

The control law is

.. math::
    I_{cmd} = C_vH_{fwd}\bigl[C_p\,(P_{cmd} - P_{fb}) - H_{fb}\,s\,P_{fb}\bigr],

so that the transfer function from the position error is :math:`C_e = C_vH_{fwd}C_p`, and the one from the
encoder position (feedback) is :math:`C_{fb} = C_vH_{fwd}(C_p + sH_{fb})`. WorkBench traces the
``Controller`` :math:`C_{WB} = C_vH_{fwd}(1 + C_p^{WB}/s)` (from the velocity feedback, feedback filters
excluded) with :math:`C_p^{WB} = 2\pi\,\texttt{PL.KP}`, and the ``Feedback Biquad`` :math:`H_{AR3}H_{AR4}`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .driver import Biquad, DriveParams


@dataclass(frozen=True)
class TF:
    """Rational transfer function num(s)/den(s)."""

    num: np.ndarray
    den: np.ndarray

    @staticmethod
    def gain(k: float) -> "TF":
        return TF(np.array([float(k)]), np.array([1.0]))

    def __mul__(self, o: "TF") -> "TF":
        return TF(np.polymul(self.num, o.num), np.polymul(self.den, o.den))

    def __add__(self, o: "TF") -> "TF":
        return TF(np.polyadd(np.polymul(self.num, o.den), np.polymul(o.num, self.den)),
                  np.polymul(self.den, o.den))

    def normalised(self) -> "TF":
        """Monic denominator."""
        a = self.den[np.nonzero(np.abs(self.den) > 0)[0][0]]
        return TF(np.trim_zeros(self.num / a, "f"), np.trim_zeros(self.den / a, "f"))

    def __call__(self, s: np.ndarray) -> np.ndarray:
        return np.polyval(self.num, s) / np.polyval(self.den, s)

    def zeros(self) -> np.ndarray:
        return np.roots(np.trim_zeros(self.num, "f"))

    def poles(self) -> np.ndarray:
        return np.roots(np.trim_zeros(self.den, "f"))

    @property
    def order(self) -> tuple[int, int]:
        return len(np.trim_zeros(self.num, "f")) - 1, len(np.trim_zeros(self.den, "f")) - 1


def biquad_tf(b: Biquad) -> TF:
    if not b.active:
        return TF.gain(1.0)
    num, den = b.coefficients()
    return TF(num, den)


def controller_tfs(p: DriveParams) -> dict[str, TF]:
    """Blocks and composite transfer functions of the controller (units: A per rad, A per rad/s, V/A)."""
    wi = 2 * np.pi * p.Kiv_hz
    Cp = TF.gain(p.Kpp) if p.Kip_hz == 0 else TF(np.array([p.Kpp, p.Kpp * 2 * np.pi * p.Kip_hz]), np.array([1.0, 0.0]))
    Cv = TF(np.array([p.Kpv, p.Kpv * wi]), np.array([1.0, 0.0]))
    Hfwd = biquad_tf(p.forward[0]) * biquad_tf(p.forward[1])
    Har = biquad_tf(p.feedback[0]) * biquad_tf(p.feedback[1])
    Hobs = TF.gain(1.0)
    if p.observer:
        wo = 2 * np.pi * p.observer_bw_hz
        Hobs = TF(np.array([wo]), np.array([1.0, wo]))
    Hfb = Har * Hobs
    s_ = TF(np.array([1.0, 0.0]), np.array([1.0]))
    Ci = TF(np.array([p.Kpi, p.Kpi * p.R_motor / p.L_motor]), np.array([1.0, 0.0]))
    kwb = p.Kpp if p.Kpp_workbench is None else p.Kpp_workbench
    Cwb = Cv * Hfwd * TF(np.array([1.0, kwb]), np.array([1.0, 0.0]))
    return {"C_p": Cp, "C_v": Cv, "H_fwd": Hfwd, "H_AR34": Har, "H_obs": Hobs, "H_fb": Hfb, "C_i": Ci,
            "C_e": (Cv * Hfwd * Cp).normalised(), "C_fb": (Cv * Hfwd * (Cp + s_ * Hfb)).normalised(),
            "C_WB": Cwb.normalised()}


def check_workbench(tfs: dict[str, TF], f: np.ndarray, controller: np.ndarray, feedback: np.ndarray) -> dict:
    """Largest magnitude [dB] and phase [deg] differences from the WorkBench traces."""
    s = 2j * np.pi * f
    out = {}
    for name, H, meas in (("controller", tfs["C_WB"](s), controller), ("feedback", tfs["H_AR34"](s), feedback)):
        e = H / meas
        out[f"{name}_dB"] = float(np.max(np.abs(20 * np.log10(np.abs(e)))))
        out[f"{name}_deg"] = float(np.max(np.abs(np.angle(e, deg=True))))
    return out


def describe_roots(r: np.ndarray) -> list[str]:
    """Roots as frequency [Hz] and damping ratio (complex pairs once), real roots as frequency."""
    out, seen = [], []
    for z in sorted(r, key=lambda z: abs(z)):
        if any(abs(z - y) < 1e-6 * max(1.0, abs(z)) for y in seen):
            continue
        if abs(z.imag) > 1e-9 * max(1.0, abs(z)):
            seen += [z, np.conj(z)]
            wn = abs(z)
            out.append(f"{wn / (2 * np.pi):.4g} Hz, zeta {-z.real / wn:.3g}")
        else:
            seen.append(z)
            out.append("0 (integrator)" if abs(z) < 1e-9 else f"{-z.real / (2 * np.pi):.4g} Hz (real)")
    return out
