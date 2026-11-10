"""One axis of the plant: identified structure, Dahl bearings, motor and drive, in cascade and in velocity mode, and the
reduced linear model used to design the LQG controller."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from ..control import lqg
from ..structure import AzimuthChain, PolarChain
from .drive import Drive, DriveParams
from .fork import ModalModel, azimuth_model, polar_model
from .friction import DahlFriction
from .loop import ClosedLoopAxis, VelocityModeDrive
from .motor import MotorParams


@dataclass
class Axis:
    """One axis.

    ``cascade``: the drive in position mode (position, velocity and current loops in the AKD2G, the configuration of
    LUCY). ``velocity_loop``: the drive in velocity mode, the plant of the LQG in the master. ``Ap, Bp, Cp``: design model
    of the LQG, the velocity loop without bearings in modal coordinates, truncated, position first. ``los``: selector of
    the node that carries the telescope (line of sight) on the node rotations."""

    name: str
    model: ModalModel
    parameters: dict
    friction: DahlFriction
    drive: Drive
    cascade: ClosedLoopAxis
    velocity_loop: ClosedLoopAxis
    Ap: np.ndarray
    Bp: np.ndarray
    Cp: np.ndarray
    modes: list
    """Modes kept in the design model (the rigid mode first)."""
    all_modes: list
    los: np.ndarray

    @property
    def lightly_damped_modes(self) -> list[lqg.Mode]:
        """Flexible modes of the design model with zeta < 0.1: the modes whose LQG weights ``q_f`` are tuned
        (key ``0`` is the first of them)."""
        return [m for m in self.modes if len(m.indices) == 2 and m.zeta < 0.1]

    def bearing_variants(self) -> dict[str, np.ndarray]:
        """A matrix of the velocity loop with the bearings sliding (no friction stiffness) and stuck (pre-sliding,
        stiffness sigma0)."""
        an = self.velocity_loop
        vv = np.r_[self.model.v, np.zeros(an.N - an.n)]
        return {"sliding": an.A, "stuck": an.A + self.friction.sigma0 * np.outer(an.B[:, 1], vv)}


def build_axis(name: str, chain: AzimuthChain | PolarChain, parameters: Mapping[str, float], friction: DahlFriction, Kt: float,
               akd: str | Path, workbench: str | Path, gains: Mapping[str, float], motor: MotorParams, dt: float,
               los_node: str, f_max: float = 150.0, hinf_min: float = 1e-5) -> Axis:
    """Builds one axis.

    ``parameters``: identified factors of the fork model; ``akd``, ``workbench``: drive parameter file and the
    WorkBench export with the filters active during the tests; ``gains``: effective ``Kpp``, ``Kpv``, ``Kiv_hz`` of the
    cascade, which replace those of the files; ``dt``: simulation step. Modes of the design model: up to ``f_max`` with
    :math:`H_\\infty` norm at least ``hinf_min``."""
    v = dict(parameters)
    m = azimuth_model(chain, v, Kt) if name == "azimuth" else polar_model(chain, v, Kt)  # type: ignore[arg-type]
    p = dataclasses.replace(DriveParams.from_files(name, akd, workbench, pl_kp_unit="hz"), **gains)
    d_c = Drive(p)
    d_v = Drive(dataclasses.replace(p, Kpp=0.0, Kip_hz=0.0, kvff=1.0))
    cascade = ClosedLoopAxis(m, d_c, motor, friction, dt)
    vel = ClosedLoopAxis(m, VelocityModeDrive(d_v), motor, friction, dt)
    Am, Bm, Cm, modes = lqg.modal_form(vel.A, vel.B[:, 0][:, None], vel.out_position[None, :])
    Ar, Br, Cr, kept = lqg.truncate(Am, Bm, Cm, modes, lambda x: x.f_hz <= f_max and x.hinf >= hinf_min)
    Ap, Bp, Cp = lqg.position_first(Ar, Br, Cr)
    los = np.zeros(len(m.M))
    los[m.nodes.index(los_node)] = 1.0
    return Axis(name, m, v, friction, d_c, cascade, vel, Ap, Bp, Cp, kept, modes, los)
