r"""Kollmorgen AKD2G drive in position mode: position, velocity and current loops.

Structure (Kollmorgen AKD User Manual, "Configuring Axis Settings": Position Loop, Velocity Loop, Current Loop
pages), without feedforward:

.. math::
    v_{cmd} = K_{pp}\Bigl(1 + \frac{2\pi K_{ip}}{s}\Bigr)(p_{cmd} - p_{fb}),\qquad
    i_{cmd} = K_{pv}\Bigl(1 + \frac{2\pi K_{iv}}{s}\Bigr)H_{AR1}H_{AR2}\bigl(v_{cmd} - H_{AR3}H_{AR4}\,v_{fb}\bigr),

.. math::
    u = K_{pi}\Bigl(1 + \frac{R}{L\,s}\Bigr)(i_{cmd} - i),

with :math:`p_{fb}`, :math:`v_{fb}` the position and velocity of the load-side encoder (the feedback of both
outer loops on LUCY), :math:`K_{pp} = 2\pi\,\texttt{PL.KP}` (1/s), :math:`K_{pv} = \texttt{VL.KP}` (A/(rad/s)),
:math:`K_{iv} = \texttt{VL.KI}` (Hz) and :math:`K_{pi} = \texttt{IL.KP}` (V/A, line-line scale).

Conventions not defined in the manual pages and fixed here from the WorkBench exports: the velocity PI is
in series form with the integral gain in Hz, and the position gain in "Hz" multiplies :math:`2\pi`; with them
the ``Controller`` trace that WorkBench stores in every frequency-response export,
:math:`K_{pv}(1 + 2\pi K_{iv}/s)\,H_{AR1}H_{AR2}\,(1 + K_{pp}/s)`, is reproduced within a few tenths of a
dB (:func:`controller_workbench`). The current PI is not described
beyond its proportional gain: its integral zero is placed on the electrical pole :math:`R/L` of the motor, so
that the open current loop is :math:`K_{pi}/(L s)` and its crossover is the idealised one of the manual,
:math:`\texttt{IL.KP}/L` rad/s (Current Loop Gain page), about 1.5 kHz.

The biquads are the continuous second-order sections that reproduce the recipes of the Velocity Loop page
(low pass, lead-lag, resonator):
:math:`H(s) = (s^2/\omega_z^2 + s/(Q_z\omega_z) + 1)/(s^2/\omega_p^2 + s/(Q_p\omega_p) + 1)`.
Limits: user current limit (``IL.LIMITP``/``IL.LIMITN``), output limit of the velocity integrator
(``VL.INTOUTMAX``), phase voltage of the bus (``VBUS.DCNOMINAL``); anti-windup by back-calculation (the pages
do not describe the scheme of the drive). The digital
implementation (sampling, computation delay) is approximated by an optional pure delay in the frequency
domain only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .workbench import WorkBenchExport, read_akd, read_workbench


T_AW = 1e-3
"""Tracking time of the back-calculation anti-windup of the velocity integrator, s."""
T_AW_CURRENT = 2e-4
"""Tracking time of the anti-windup of the current integrator, s."""


# --------------------------------------------------------------------------- biquad filters


@dataclass(frozen=True)
class Biquad:
    """One AR filter: numerator (fz, qz) and denominator (fp, qp), type as in ``VL.ARTYPE``."""

    fz: float
    qz: float
    fp: float
    qp: float
    kind: int = 4

    @property
    def active(self) -> bool:
        return self.kind != 0 and not (self.fz == self.fp and self.qz == self.qp)

    def coefficients(self) -> tuple[np.ndarray, np.ndarray]:
        wz, wp = 2 * np.pi * self.fz, 2 * np.pi * self.fp
        return (np.array([1 / wz**2, 1 / (self.qz * wz), 1.0]), np.array([1 / wp**2, 1 / (self.qp * wp), 1.0]))

    def response(self, s: np.ndarray) -> np.ndarray:
        if not self.active:
            return np.ones_like(s, dtype=complex)
        num, den = self.coefficients()
        return np.polyval(num, s) / np.polyval(den, s)

    def state_space(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Realisation with states in the units of the input: :math:`z_1` is the input through the
        denominator (a low pass of unit static gain) and :math:`z_2 = \\dot z_1/\\omega_p`, so that the
        integrator tolerances apply to all states alike."""
        wz, wp = 2 * np.pi * self.fz, 2 * np.pi * self.fp
        r = wp**2 / wz**2
        A = np.array([[0.0, wp], [-wp, -wp / self.qp]])
        B = np.array([[0.0], [wp]])
        C = np.array([[1.0 - r, wp / (self.qz * wz) - r / self.qp]])
        D = np.array([[r]])
        return A, B, C, D

    def describe(self) -> str:
        if not self.active:
            return "unity"
        if self.fz == self.fp:
            return f"resonator {self.fp:g} Hz, G = {20 * np.log10(self.qp / self.qz):+.1f} dB, Q = {np.sqrt(self.qz * self.qp):.2g}"
        return f"num {self.fz:g} Hz/Q {self.qz:g}, den {self.fp:g} Hz/Q {self.qp:g}"


@dataclass
class _Chain:
    """Biquads in series, as one state-space block."""

    blocks: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = field(default_factory=list)

    @property
    def n(self) -> int:
        return sum(b[0].shape[0] for b in self.blocks)

    def step(self, x: np.ndarray, u: float) -> tuple[float, np.ndarray]:
        """Output and state derivative of the chain for input ``u``."""
        dx = np.empty_like(x)
        k = 0
        y = u
        for A, B, C, D in self.blocks:
            m = A.shape[0]
            xi = x[k:k + m]
            dx[k:k + m] = A @ xi + B[:, 0] * y
            y = float(C[0] @ xi + D[0, 0] * y)
            k += m
        return y, dx


# --------------------------------------------------------------------------- parameters


@dataclass(frozen=True)
class DriveParams:
    """Parameters of one drive, SI units."""

    name: str
    Kpp: float
    """Position gain of the drive, 1/s: 2 pi PL.KP (``pl_kp_unit`` = ``hz``, the convention of the WorkBench
    Controller trace) or PL.KP (``per_s``)."""
    Kip_hz: float
    Kpv: float
    """Velocity proportional gain, A/(rad/s)."""
    Kiv_hz: float
    forward: tuple[Biquad, Biquad]
    feedback: tuple[Biquad, Biquad]
    Kpi: float
    """Current proportional gain, V/A (line-line scale)."""
    R_motor: float
    L_motor: float
    """Motor resistance and inductance known to the drive (MOTOR.R, MOTOR.LQLL), line-line."""
    i_limit: float
    vel_int_limit: float
    v_bus: float | None = None
    """DC bus voltage (VBUS.DCNOMINAL), V; sets the largest phase voltage."""
    kvff: float = 0.0
    """Velocity feedforward gain; LUCY uses no feedforward."""
    observer: bool = False
    observer_bw_hz: float = 0.0
    delay_s: float = 0.0
    source: dict[str, str] = field(default_factory=dict)
    Kpp_workbench: float | None = None
    """Position gain used by WorkBench in its Controller trace, 2 pi PL.KP (``None``: ``Kpp``)."""

    @property
    def current_bandwidth_hz(self) -> float:
        return self.Kpi / self.L_motor / (2 * np.pi)

    @property
    def V_phase_max(self) -> float | None:
        """Largest per-phase rms voltage with sinusoidal modulation of the bus, :math:`V_{dc}/\\sqrt6`."""
        return None if self.v_bus is None else self.v_bus / np.sqrt(6)

    @classmethod
    def from_files(cls, name: str, akd: str | Path, workbench: str | Path | None = None, *,
                   feedforward: bool = False, observer: bool | None = None, delay_s: float = 0.0,
                   pl_kp_unit: str = "hz") -> "DriveParams":
        """Parameters from the ``.akd`` file of the drive; gains and filters are taken from the WorkBench export
        when given (the configuration active during that test), the rest from the ``.akd`` file. ``pl_kp_unit``:
        ``hz`` (position gain 2 pi PL.KP, as WorkBench computes its Controller trace) or ``per_s`` (PL.KP in
        1/s, suggested by the time sequences of September)."""
        a = read_akd(akd)
        g = lambda k: float(a[f"AXIS1.{k}"])  # noqa: E731
        src = {"base": str(akd)}
        Kpp_hz, Kip, Kpv, Kiv = g("PL.KP"), g("PL.KI"), g("VL.KP"), g("VL.KI")
        bq = [Biquad(g(f"VL.ARZF{k}"), g(f"VL.ARZQ{k}"), g(f"VL.ARPF{k}"), g(f"VL.ARPQ{k}"), int(g(f"VL.ARTYPE{k}")))
              for k in (1, 2, 3, 4)]
        Kpi = None
        if workbench is not None:
            w: WorkBenchExport = read_workbench(workbench)
            Kpp_hz, Kip = w.tuning("PositionLoopKP"), w.tuning("PositionLoopKI")
            Kpv, Kiv = w.tuning("VelocityLoopKP"), w.tuning("VelocityLoopKI")
            Kpi = w.tuning("CurrentLoopKP")
            fl = w.filters()
            fwd = [f for f in fl if f[0] == "ArfGroupForward"]
            fb = [f for f in fl if f[0] == "ArfGroupFeedback"]
            bq = [Biquad(*f[1:]) for f in fwd[:2]] + [Biquad(*f[1:]) for f in fb[:2]]
            src["gains and filters"] = str(workbench)
        R, L = g("MOTOR.R"), g("MOTOR.LQLL") * 1e-3
        if Kpi is None:  # IL.KP is read-only and not stored in the .akd file: from IL.BW = IL.KP / L
            Kpi = 2 * np.pi * g("IL.BW") * 1e3 * L
        obs = bool(int(g("OBS.ENABLE"))) if observer is None else observer
        v_bus = float(a["VBUS.DCNOMINAL"]) if "VBUS.DCNOMINAL" in a else None
        if pl_kp_unit not in ("hz", "per_s"):
            raise ValueError(f"pl_kp_unit {pl_kp_unit!r}: hz | per_s")
        Kpp = 2 * np.pi * Kpp_hz if pl_kp_unit == "hz" else Kpp_hz
        return cls(name=name, Kpp=Kpp, Kpp_workbench=2 * np.pi * Kpp_hz, Kip_hz=Kip, Kpv=Kpv, Kiv_hz=Kiv,
                   forward=(bq[0], bq[1]), feedback=(bq[2], bq[3]), Kpi=Kpi, R_motor=R, L_motor=L,
                   i_limit=min(abs(g("IL.LIMITP")), abs(g("IL.LIMITN"))), vel_int_limit=g("VL.INTOUTMAX"), v_bus=v_bus,
                   kvff=g("VL.KVFF") if feedforward else 0.0, observer=obs, observer_bw_hz=g("OBS.BW"),
                   delay_s=delay_s, source=src)


# --------------------------------------------------------------------------- drive


class Drive:
    """AKD2G drive of one axis. States: (position integrator, with PL.KI > 0), forward biquads, velocity
    integrator, feedback biquads, (observer low pass), current integrator."""

    def __init__(self, p: DriveParams, V_phase_max: float | None = None):
        self.p = p
        self.V_phase_max = p.V_phase_max if V_phase_max is None else V_phase_max
        self._fwd = _Chain([b.state_space() for b in p.forward if b.active])
        self._fb = _Chain([b.state_space() for b in p.feedback if b.active])
        # the position integrator exists only with PL.KI > 0 (an unconnected integrator would be a zero
        # column of the Jacobian, on which the implicit solver's finite differences diverge)
        n_zp = 1 if p.Kip_hz != 0 else 0
        n = n_zp + self._fwd.n + 1 + self._fb.n + (1 if p.observer else 0) + 1
        self.n_states = n
        k = 0
        self.i_zp = k if n_zp else None
        k += n_zp
        self.s_fwd = slice(k, k + self._fwd.n); k += self._fwd.n
        self.i_zv = k; k += 1
        self.s_fb = slice(k, k + self._fb.n); k += self._fb.n
        self.i_obs = k if p.observer else None
        k += 1 if p.observer else 0
        self.i_zi = k

    # ---- time domain -----------------------------------------------------------------------------------

    def step(self, x: np.ndarray, p_cmd: float, dp_cmd: float, p_fb: float, v_fb: float, i_meas: float,
             i_exc: float = 0.0) -> tuple[float, float, np.ndarray]:
        """Phase voltage (V rms), current command (A rms) and state derivative."""
        p = self.p
        dx = np.zeros(self.n_states)
        e_p = p_cmd - p_fb
        v_cmd = p.Kpp * e_p + p.kvff * dp_cmd
        if self.i_zp is not None:
            dx[self.i_zp] = e_p
            v_cmd += p.Kpp * 2 * np.pi * p.Kip_hz * x[self.i_zp]
        v_meas = v_fb
        if p.observer:  # the manual describes the observer as a low pass of bandwidth OBS.BW on the feedback
            wo = 2 * np.pi * p.observer_bw_hz
            dx[self.i_obs] = wo * (v_fb - x[self.i_obs])
            v_meas = x[self.i_obs]
        v_f, dx[self.s_fb] = self._fb.step(x[self.s_fb], v_meas)
        e_v = v_cmd - v_f
        e_vf, dx[self.s_fwd] = self._fwd.step(x[self.s_fwd], e_v)
        # velocity PI with the output clamp of the integrator (VL.INTOUTMAX) and the current limit; anti-windup
        # by back-calculation with tracking time T_AW (continuous, unlike conditional integration, so that the
        # implicit solver sees a Lipschitz right-hand side; the AKD pages do not describe their scheme)
        wi = 2 * np.pi * p.Kiv_hz
        i_int_raw = p.Kpv * wi * x[self.i_zv]
        i_int = float(np.clip(i_int_raw, -p.vel_int_limit, p.vel_int_limit))
        i_unsat = p.Kpv * e_vf + i_int
        i_cmd = float(np.clip(i_unsat, -p.i_limit, p.i_limit))
        if wi > 0:
            dx[self.i_zv] = e_vf + ((i_cmd - i_unsat) + (i_int - i_int_raw)) / (p.Kpv * wi * T_AW)
        i_cmd += i_exc
        # current PI (zero on R/L), per-phase voltage limited by the bus, back-calculation as above
        e_i = i_cmd - i_meas
        wr = p.R_motor / p.L_motor
        u = p.Kpi * (e_i + wr * x[self.i_zi])
        v_phase = u / 2.0
        if self.V_phase_max is not None:
            v_phase = float(np.clip(v_phase, -self.V_phase_max, self.V_phase_max))
        dx[self.i_zi] = e_i + (2.0 * v_phase - u) / (p.Kpi * wr * T_AW_CURRENT)
        return v_phase, i_cmd, dx

    # ---- frequency domain ------------------------------------------------------------------------------

    def C_pos(self, s: np.ndarray) -> np.ndarray:
        return self.p.Kpp * (1 + 2 * np.pi * self.p.Kip_hz / s)

    def C_vel(self, s: np.ndarray) -> np.ndarray:
        return self.p.Kpv * (1 + 2 * np.pi * self.p.Kiv_hz / s)

    def F_fwd(self, s: np.ndarray) -> np.ndarray:
        return self.p.forward[0].response(s) * self.p.forward[1].response(s)

    def F_ar_fb(self, s: np.ndarray) -> np.ndarray:
        """Feedback biquads AR3, AR4 (the ``Feedback Biquad`` trace of WorkBench, which excludes the observer)."""
        return self.p.feedback[0].response(s) * self.p.feedback[1].response(s)

    def F_fb(self, s: np.ndarray) -> np.ndarray:
        """Whole velocity feedback path: AR3, AR4 and the observer, which the Velocity Loop page describes as a
        loop that "behaves like a low-pass filter" of bandwidth ``OBS.BW`` (first order here)."""
        h = self.F_ar_fb(s)
        if self.p.observer:
            h = h / (1 + s / (2 * np.pi * self.p.observer_bw_hz))
        return h

    def delay(self, s: np.ndarray) -> np.ndarray:
        return np.exp(-s * self.p.delay_s)

    def C_cur(self, s: np.ndarray) -> np.ndarray:
        """Current PI in per-phase volts per ampere."""
        return self.p.Kpi * (1 + self.p.R_motor / self.p.L_motor / s) / 2.0


    def n_active_filters(self) -> int:
        return sum(b.active for b in self.p.forward + self.p.feedback)


def controller_workbench(d: Drive, s: np.ndarray) -> np.ndarray:
    """``Controller`` trace of WorkBench: from the velocity feedback to the current command, position loop
    included, :math:`C_v H_{AR1}H_{AR2}(1 + C_p/s)`, with the position gain of WorkBench (2 pi PL.KP)."""
    p = d.p
    kpp = p.Kpp if p.Kpp_workbench is None else p.Kpp_workbench
    Cp = kpp * (1 + 2 * np.pi * p.Kip_hz / s)
    return d.C_vel(s) * d.F_fwd(s) * (1 + Cp / s)
