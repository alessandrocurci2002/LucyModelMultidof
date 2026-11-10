r"""Kollmorgen TBM(S)-7646-A frameless brushless torque motor.

Field-oriented model with ideal commutation (the drive keeps the current on the :math:`q` axis, so the
current is the torque-producing component and :math:`i_d = 0`). Per-phase rms equivalent of the
three-phase winding, with the line-line datasheet values split on a star connection:

.. math::
    L_{ph}\,\frac{di}{dt} = v - R_{ph}\,i - E_{ph},\qquad R_{ph} = \frac{R_{LL}}{2},\quad
    L_{ph} = \frac{L_{LL}}{2},\quad E_{ph} = \frac{K_b}{\sqrt3}\,\omega_m,

.. math::
    \tau_m = K_t\,i - b_m\,\omega_m - \tau_{cog}(\theta_m),

with :math:`i` in A rms (the unit of the drive currents), :math:`K_b` the line-line rms back-EMF constant
and :math:`\omega_m` the rotor speed. The power balance :math:`3E_{ph}i = K_t i\,\omega_m` requires
:math:`K_t = \sqrt3\,K_b`, which the datasheet values satisfy (:meth:`MotorParams.consistency`). The rotor
inertia :math:`J_m` belongs to the mechanical model (motor node of :mod:`lucy_sim.plant.fork`).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import yaml

RPM_TO_RADS = 2 * np.pi / 60.0


@dataclass(frozen=True)
class MotorParams:
    """Datasheet parameters, SI units (winding temperature chosen at construction)."""

    name: str
    Kt: float
    """Torque constant, N m/A rms."""
    Kb: float
    """Back-EMF constant, line-line V rms per rad/s."""
    R_ll: float
    L_ll: float
    J_rotor: float
    poles: int
    I_cont: float
    I_peak: float
    V_bus: float
    """Design bus voltage of the winding, V dc."""
    b_viscous: float
    """Viscous damping, N m s/rad."""
    T_static: float
    T_cog_pp: float
    temperature: str

    @classmethod
    def from_yaml(cls, path: str | Path, temperature: Literal["cold", "hot"] = "cold") -> "MotorParams":
        d = yaml.safe_load(Path(path).read_text())
        return cls(name=d["name"],
                   Kt=float(d[f"torque_constant_{temperature}_Nm_per_Arms"]),
                   Kb=float(d[f"back_emf_{temperature}_Vrms_per_krpm"]) / (1000 * RPM_TO_RADS),
                   R_ll=float(d["resistance_ohm"]), L_ll=float(d["inductance_mH"]) * 1e-3,
                   J_rotor=float(d["rotor_inertia_kgm2"]), poles=int(d["poles"]),
                   I_cont=float(d["continuous_current_Arms"]), I_peak=float(d["peak_current_Arms"]),
                   V_bus=float(d["design_bus_voltage_Vdc"]),
                   b_viscous=float(d["viscous_damping_Nm_per_krpm"]) / (1000 * RPM_TO_RADS),
                   T_static=float(d["max_static_friction_Nm"]), T_cog_pp=float(d["cogging_peak_to_peak_Nm"]),
                   temperature=temperature)

    @property
    def R_ph(self) -> float:
        return self.R_ll / 2

    @property
    def L_ph(self) -> float:
        return self.L_ll / 2

    @property
    def Ke_ph(self) -> float:
        """Per-phase rms back-EMF per rad/s."""
        return self.Kb / np.sqrt(3)

    def consistency(self) -> float:
        """Relative difference between :math:`K_t` and :math:`\\sqrt3 K_b` (zero for a consistent datasheet)."""
        return float(self.Kt / (np.sqrt(3) * self.Kb) - 1.0)

    @property
    def V_phase_max(self) -> float:
        """Largest per-phase rms voltage with sinusoidal modulation of the bus, :math:`V_{bus}/\\sqrt6`."""
        return self.V_bus / np.sqrt(6)
