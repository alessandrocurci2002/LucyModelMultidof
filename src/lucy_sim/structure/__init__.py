"""Structural model of the gimbal: geometry -> torsional FEM -> lumped chains -> coupled Lagrange model -> linear chains.

:func:`build_structure` runs the whole chain from ``config/gimbal.yaml``. Its output used by the plant is the pair of
linearised axis chains (:class:`AzimuthChain`, :class:`PolarChain`, at theta = 0, Omega = 0), on which the identified
fork model is built (:mod:`lucy_sim.plant.fork`), and the telescope payload, which gives the coupling torques between
the axes (:mod:`lucy_sim.simulation.coupling`).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import GimbalConfig, load_config
from .fem import AxisFEM, build_axis, stiffness_proportional_beta
from .lagrange import CoupledModel, coupled_model
from .lumped import LumpedChain, lumped_chain
from .payload import Payload, build_payload

BETA_THETA_REF = 0.0
"""Operating point at which the stiffness-proportional damping is tuned (it must not depend on the pointing)."""


@dataclass(frozen=True)
class PolarChain:
    """Polar block (coordinates e, h, m, r, p) of the linearised model at theta = 0."""

    M: np.ndarray
    C: np.ndarray
    K: np.ndarray
    roles: tuple[str, ...]


@dataclass(frozen=True)
class AzimuthChain:
    """Azimuth chain e - m - L - p: inertias, springs (e-m, m-L, L-p) and dampers of the linearised model."""

    J_e: float
    J_m: float
    J_L: float
    J_p: float
    k: tuple[float, float, float]
    c: tuple[float, float, float]


@dataclass(frozen=True)
class StructuralModel:
    config: GimbalConfig
    fem_azimuth: AxisFEM
    fem_polar: AxisFEM
    payload: Payload
    chain_azimuth: LumpedChain
    chain_polar: LumpedChain
    coupled: CoupledModel
    azimuth: AzimuthChain
    polar: PolarChain


def axis_chains(acc: CoupledModel, theta0: float = 0.0) -> tuple[AzimuthChain, PolarChain]:
    """Chains of the two axes from the linearised coupled model (products of inertia between the axes, ~1e-6, neglected)."""
    st = acc.structure
    M, C, K = acc.linearize(theta0, 0.0)
    ip = np.arange(st.n_psi, st.n)
    pol = PolarChain(*(X[np.ix_(ip, ip)].real.copy() for X in (M, C, K)), tuple(st.roles_theta))
    if list(st.roles_psi) != ["e", "m", "L", "p"]:
        raise ValueError(f"unexpected azimuth chain {st.roles_psi}: the telescope must hang from the load")
    az = AzimuthChain(*(float(M[i, i]) for i in range(4)), tuple(float(-K[i, i + 1]) for i in range(3)),  # type: ignore[arg-type]
                      tuple(float(-C[i, i + 1]) for i in range(3)))
    return az, pol


def build_structure(config: GimbalConfig | str | Path) -> StructuralModel:
    """Whole structural chain from the configuration (about 10 s, mostly the SymPy derivation)."""
    cfg = config if isinstance(config, GimbalConfig) else load_config(config)
    fa, fp = build_axis(cfg.azimuth), build_axis(cfg.polar)
    pl = build_payload(cfg, fp)
    ca, cp = lumped_chain(fa), lumped_chain(fp)
    # damping tuned on the first flexible mode of each isolated axis, telescope and fork included
    Ma = ca.M.copy()
    Ma[ca.index("L"), ca.index("L")] += pl.J_ax
    Ma[ca.index("p"), ca.index("p")] += float(pl.J_psi(BETA_THETA_REF))
    Mp = cp.M.copy()
    Mp[cp.index("p"), cp.index("p")] += pl.I_yy
    beta_psi = stiffness_proportional_beta(Ma, ca.K, cfg.zeta)
    beta_theta = stiffness_proportional_beta(Mp, cp.K, cfg.zeta)
    acc = coupled_model(ca, cp, pl, beta_psi=beta_psi, beta_theta=beta_theta)
    az, pol = axis_chains(acc)
    return StructuralModel(cfg, fa, fp, pl, ca, cp, acc, az, pol)


__all__ = ["AzimuthChain", "PolarChain", "StructuralModel", "build_structure", "load_config", "axis_chains"]
