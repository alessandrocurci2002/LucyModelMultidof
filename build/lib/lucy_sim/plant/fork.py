r"""Identified modal model of each axis (structure without the bearings).

The CAD chains of :mod:`lucy_sim.structure` are scaled with the factors identified on the measured frequency responses
(``config/plant.yaml``): inertia and stiffness factors within the CAD tolerance, torque constant within the catalogue
tolerance, free modal damping of each mode. The bearings are not a linear constraint here: the bearing deformation
:math:`v^T x` drives the Dahl friction of :mod:`lucy_sim.plant.friction`, so the first mode is the rigid motion of the
axis.

**Polar** (nodes e, hp, m, r, q1, q2): telescope rigid with the hub (``hp``), motor rotor ``m``, rear end ``r`` and two
modal masses on the motor stator (``q1``, ``q2``). Each is a structural mode normalised on the stator, inertia
:math:`J_r = 1/|a_r|`, spring :math:`J_r\omega_r^2` to ground; the stator reaction :math:`-\tau` acts on them and the
encoder head turns with :math:`\sum\kappa_r q_r`, :math:`\kappa_r = \mathrm{sign}(a_r)`.

**Azimuth** (nodes e, m, L, AP): the load is split into the fork base L (on the bearings) and the fork top A with the
telescope (AP), joined by :math:`k_{LA}`; a single CAD factor scales load and telescope.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np
from scipy.linalg import eigh

from ..structure import AzimuthChain, PolarChain

ZETA_MIN = 1e-4
"""Damping of the rigid-body mode (the bearings are outside the linear model)."""
NODES_POLAR = ("e", "hp", "m", "r", "q1", "q2")
NODES_AZIMUTH = ("e", "m", "L", "AP")


@dataclass(frozen=True)
class ModalModel:
    """Structure of one axis: ``M, K`` with modal viscous damping ``zeta`` (modes in increasing frequency).

    ``e`` maps the motor torque to the nodes, ``c`` gives the encoder reading, ``v`` the bearing deformation (the Dahl
    friction torque acts on the nodes through ``-v``)."""

    M: np.ndarray
    K: np.ndarray
    v: np.ndarray
    e: np.ndarray
    c: np.ndarray
    zeta: np.ndarray
    delay: float
    """Equivalent drive delay of the frequency-response fit, s (not used in the time domain)."""
    Kt: float
    """Torque constant with the identified factor, N m/A."""
    nodes: tuple[str, ...]
    omega: np.ndarray = field(init=False, repr=False)
    Phi: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        w2, Phi = eigh(self.K, self.M)
        object.__setattr__(self, "omega", np.sqrt(np.maximum(w2, 0.0)))
        object.__setattr__(self, "Phi", Phi)

    @property
    def damping_matrix(self) -> np.ndarray:
        """:math:`C = M\\Phi\\,\\mathrm{diag}(2\\zeta\\omega)\\,\\Phi^T M`."""
        return self.M @ self.Phi @ np.diag(2 * self.zeta * self.omega) @ self.Phi.T @ self.M


def springs(K: np.ndarray) -> tuple[list[tuple[int, int, float]], np.ndarray]:
    """Springs ``(i, j, k)`` (i < j) and stiffness to ground per node of a chain stiffness matrix."""
    n = K.shape[0]
    return [(i, j, float(-K[i, j])) for i in range(n) for j in range(i + 1, n) if K[i, j] != 0.0], K.sum(axis=1)


def _spring(K: np.ndarray, i: int, j: int | None, k: float) -> None:
    K[i, i] += k
    if j is not None:
        K[j, j] += k
        K[i, j] -= k
        K[j, i] -= k


def _zeta(v: Mapping[str, float], n: int) -> np.ndarray:
    return np.array([ZETA_MIN if k == 1 else v[f"zeta_{k}"] for k in range(1, n + 1)])


def polar_model(chain: PolarChain, v: Mapping[str, float], Kt: float) -> ModalModel:
    """Polar axis with the telescope rigid with the hub and two modal masses on the motor stator."""
    roles = list(chain.roles)
    ie, ih, im, ir, ip = (roles.index(x) for x in ("e", "h", "m", "r", "p"))
    J = np.diag(chain.M)
    arcs, ground = springs(chain.K)
    k = {frozenset((i, j)): kk for i, j, kk in arcs}
    n = len(NODES_POLAR)
    M, K = np.zeros((n, n)), np.zeros((n, n))
    M[0, 0] = v["s_J_e"] * J[ie]
    M[1, 1] = v["s_J_hp"] * (J[ih] + J[ip])
    M[2, 2] = v["s_J_m"] * J[im]
    M[3, 3] = v["s_J_r"] * J[ir]
    _spring(K, 0, 1, v["s_k_1"] * k[frozenset((ie, ih))])
    _spring(K, 1, 2, v["s_k_2"] * k[frozenset((ih, im))])
    _spring(K, 2, 3, v["s_k_3"] * k[frozenset((im, ir))])
    _spring(K, 1, None, ground[ih] + ground[ip])  # residual to ground of the CAD chain (~0.8 N m/rad), on the merged node
    e, c, vv = np.zeros(n), np.zeros(n), np.zeros(n)
    e[2], c[0], vv[1] = 1.0, 1.0, 1.0
    for r in (1, 2):
        j = 3 + r
        Jr = 1.0 / abs(v[f"a_{r}"])
        M[j, j] = Jr
        _spring(K, j, None, Jr * (2 * np.pi * v[f"f_{r}"]) ** 2)
        e[j], c[j], vv[j] = -1.0, -np.sign(v[f"a_{r}"]), -1.0
    return ModalModel(M, K, vv, e, c, _zeta(v, n), v["tau_d"], Kt * v.get("g_Kt", 1.0), NODES_POLAR)


def azimuth_model(chain: AzimuthChain, v: Mapping[str, float], Kt: float) -> ModalModel:
    """Azimuth axis with the load split into fork base and fork top, the telescope rigid with the top."""
    n = len(NODES_AZIMUTH)
    s = v["s_J_L"]
    M = np.diag([v["s_J_e"] * chain.J_e, v["s_J_m"] * chain.J_m, (1 - v["alpha"]) * s * chain.J_L,
                 v["alpha"] * s * chain.J_L + s * chain.J_p])
    K = np.zeros((n, n))
    _spring(K, 0, 1, v["s_k_1"] * chain.k[0])
    _spring(K, 1, 2, v["s_k_2"] * chain.k[1])
    _spring(K, 2, 3, v["k_LA"])
    vv, e, c = np.zeros(n), np.zeros(n), np.zeros(n)
    vv[2], e[1], c[0] = 1.0, 1.0, 1.0
    return ModalModel(M, K, vv, e, c, _zeta(v, n), v["tau_d"], Kt * v.get("g_Kt", 1.0), NODES_AZIMUTH)
