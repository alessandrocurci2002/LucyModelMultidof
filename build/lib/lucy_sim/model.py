"""The whole simulator: structure, plant of the two axes, controllers and coupling, from the files in ``config/``.

    from lucy_sim import GimbalModel, Disturbance
    model = GimbalModel.load()
    t = model.time(36.0)
    ref = {"azimuth": psi(t), "polar": theta(t)}             # absolute line of sight [rad]
    out = model.simulate(ref, t, controller="lqg", disturbance=Disturbance(gauss=0.3, band=(0.5, 5.0)), seed=10)
    out["azimuth"]["los"]                                     # line of sight relative to the start [rad]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

import numpy as np
import yaml

from .control import lqg
from .paths import ROOT
from .plant.axis import Axis, build_axis
from .plant.friction import DahlFriction
from .plant.motor import MotorParams
from .simulation.coupling import LagrangeCoupling
from .simulation.disturbances import Disturbance
from .simulation.loops import run_lqg, simulate_two_axes
from .structure import StructuralModel, build_structure

AXES = ("azimuth", "polar")


def _yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


@dataclass
class GimbalModel:
    """Simulator of LUCY.

    ``axes``: the two :class:`lucy_sim.plant.axis.Axis`; ``controllers[axis]``: ``{"cascade": None, "lqg": <LQG-PI>}``
    (``None`` means the cascade of the drive, the current configuration of LUCY); ``coupling``: torques exchanged by
    the axes through the telescope; ``dt``: simulation step; ``dt_master``, ``delay``: master period and EtherCAT delay
    in master periods."""

    structure: StructuralModel
    axes: dict[str, Axis]
    controllers: dict[str, dict[str, object]]
    coupling: LagrangeCoupling
    dt: float
    dt_master: float
    delay: int
    lqg_settings: dict = field(default_factory=dict)

    @classmethod
    def load(cls, root: str | Path = ROOT, gimbal: str = "config/gimbal.yaml", plant: str = "config/plant.yaml",
             controllers: str = "config/controllers.yaml") -> "GimbalModel":
        """Builds everything (about 15 s: the symbolic derivation of the structure and the LQG synthesis)."""
        root = Path(root)
        S = build_structure(root / gimbal)
        P, C = _yaml(root / plant), _yaml(root / controllers)
        motor = MotorParams.from_yaml(root / P["motor"], P.get("motor_temperature", "cold"))
        dt = float(P["simulation_step_s"])
        dm = C["design_model"]
        axes = {}
        for name in AXES:
            p = P["axes"][name]
            chain = S.azimuth if name == "azimuth" else S.polar
            axes[name] = build_axis(name, chain, p["fork"], DahlFriction(float(p["bearings"]["sigma0"]), float(p["bearings"]["Tc"])),
                                    motor.Kt, root / p["drive_parameters"], root / p["drive_workbench_export"], p["drive_gains"],
                                    motor, dt, p["los_node"], float(dm["f_max_Hz"]), float(dm["hinf_min"]))
        coupling = LagrangeCoupling(S.payload, float(P["axes"]["azimuth"]["fork"]["s_J_L"]))
        model = cls(S, axes, {}, coupling, dt, float(C["master"]["period_s"]), int(C["master"]["delay_periods"]), C["lqg"])
        for name in AXES:
            w = C["lqg"]["axes"][name]
            model.controllers[name] = {"cascade": None, "lqg": model.lqg_controller(name, w["q_p"], w["q_i"], w.get("q_f") or {}, w.get("rho", 1.0))}
        return model

    def lqg_controller(self, axis: str, q_p: float, q_i: float, q_f: Mapping[int, float] | None = None, rho: float = 1.0,
                       measured_pi: bool | None = None, base: float | None = None) -> lqg.LQGController:
        """LQG-PI of one axis for the weights ``q_p``, ``q_i``, ``q_f`` (index of the lightly damped mode -> weight of
        its two states; the other flexible states get ``base``) and ``V = rho Q``."""
        a = self.axes[axis]
        base = float(self.lqg_settings.get("base_flexible_weight", 0.1)) if base is None else base
        q = np.full(len(a.Ap) - 1, base)
        for k, x in (q_f or {}).items():
            q[[i - 1 for i in a.lightly_damped_modes[int(k)].indices]] = float(x)
        if measured_pi is None:
            measured_pi = self.lqg_settings.get("pi_on", "measured") == "measured"
        return lqg.synthesize(a.Ap, a.Bp, a.Cp, self.dt_master, q_i, q_p, q, rho_e=rho, v_continuous=True, measured_pi=measured_pi)

    def time(self, duration: float) -> np.ndarray:
        """Time grid of the simulator, ``[0, duration)`` at the simulation step."""
        return np.arange(0.0, duration, self.dt)

    def simulate(self, reference: Mapping[str, np.ndarray], t: np.ndarray, controller: str | Mapping[str, object] = "lqg",
                 disturbance: Disturbance | None = None, seed: int = 1, coupled: bool = True, **kw) -> dict[str, dict[str, np.ndarray]]:
        """Two coupled axes following ``reference`` (absolute azimuth and elevation of the line of sight [rad] on the
        grid ``t``). ``controller``: ``"lqg"``, ``"cascade"`` or a mapping axis -> controller (``None`` = cascade).
        ``coupled = False``: no coupling torques. Other keywords go to
        :func:`lucy_sim.simulation.loops.simulate_two_axes`."""
        C = {q: self.controllers[q][controller] for q in AXES} if isinstance(controller, str) else dict(controller)
        return simulate_two_axes(self.axes, C, dict(reference), t, disturbance or Disturbance(), seed,
                                 self.coupling if coupled else None, float(reference["polar"][0]), self.dt_master, self.delay, **kw)

    def simulate_axis(self, axis: str, p_cmd: np.ndarray, controller: str | object = "lqg", I0: float = 0.0,
                      torque: np.ndarray | None = None) -> dict[str, np.ndarray]:
        """One axis alone on the position command ``p_cmd`` (relative to the start, at the simulation step); cascade or
        a controller in the master. ``I0``: holding current at the start."""
        a = self.axes[axis]
        c = self.controllers[axis][controller] if isinstance(controller, str) else controller
        if c is None:
            return a.cascade.run(p_cmd, I0, torque)
        return run_lqg(a, c, p_cmd, I0, self.dt_master, self.delay, torque)
