r"""Lumped-parameter chain of each axis on its master nodes.

Every equivalent element condenses the FEM elements between two consecutive master nodes,

.. math:: J_{eq} = \sum_i J_i, \qquad \frac{1}{k_{eq}} = \sum_i \frac{1}{k_i}

(springs in series); the concentrated inertias of the intermediate nodes go to the equivalent element that contains
them. Each master node carries its own inertia plus half of the inertia of the adjacent equivalent elements. The
branch to the telescope P (massless spring) is an equivalent element of its own, the last one.

Master roles: azimuth ``e`` (encoder), ``m`` (motor), ``L`` (load: fork + hub + ring), ``p`` (telescope); polar
``e``, ``h`` (hub centre), ``m``, ``r`` (rear end), ``p``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fem import AxisFEM, stiffness_matrix

ROLES = {"azimuth": ("e", "m", "L", "p"), "polar": ("e", "h", "m", "r", "p")}


@dataclass(frozen=True)
class EquivalentElement:
    a: str
    b: str
    J: float
    k: float
    ia: int
    ib: int
    """Positions of the two nodes among the masters."""


def guyan_transformation(K: np.ndarray, masters) -> np.ndarray:
    """Static condensation matrix :math:`T` (n x m), rows in global order."""
    n = K.shape[0]
    m = list(masters)
    s = [i for i in range(n) if i not in m]
    T = np.zeros((n, len(m)))
    T[m, np.arange(len(m))] = 1.0
    if s:
        T[s, :] = -np.linalg.solve(K[np.ix_(s, s)], K[np.ix_(s, m)])
    return T


def equivalent_elements(fem: AxisFEM) -> tuple[EquivalentElement, ...]:
    masters = list(fem.masters)
    chain = [m for m in masters if m != fem.i_P]
    out = []
    for a, b in zip(chain[:-1], chain[1:]):
        group = [p for p, (i, j) in zip(fem.elements, fem.links) if a <= i and j <= b]
        if not group:
            raise ValueError(f"axis {fem.name}: no element between nodes {fem.label(a)} and {fem.label(b)}")
        J_int = float(sum(fem.J_nodes[i] for i in range(a + 1, b)))
        out.append(EquivalentElement(fem.label(a), fem.label(b), float(sum(p.J for p in group)) + J_int,
                                     1.0 / float(sum(1.0 / p.k for p in group)), masters.index(a), masters.index(b)))
    group = [p for p, c in zip(fem.elements, fem.links) if c == (fem.i_H, fem.i_P)]
    out.append(EquivalentElement(fem.label(fem.i_H), "P", float(sum(p.J for p in group)) + 0.0,
                                 1.0 / float(sum(1.0 / p.k for p in group)), masters.index(fem.i_H), masters.index(fem.i_P)))
    return tuple(out)


@dataclass(frozen=True)
class LumpedChain:
    """Lumped chain of one axis on the master nodes."""

    axis: str
    roles: tuple[str, ...]
    M: np.ndarray
    K: np.ndarray
    C_b: np.ndarray
    fem: AxisFEM

    @property
    def n(self) -> int:
        return len(self.roles)

    def index(self, role: str) -> int:
        return self.roles.index(role)

    @property
    def arcs(self) -> tuple[tuple[int, int], ...]:
        """Springs between masters: consecutive masters of the chain, then the branch to the telescope (last)."""
        i_p = self.fem.masters.index(self.fem.i_P)
        chain = [c for c in range(self.n) if c != i_p]
        return tuple(zip(chain[:-1], chain[1:])) + ((self.fem.masters.index(self.fem.i_H), i_p),)


def lumped_chain(fem: AxisFEM) -> LumpedChain:
    """Lumped chain on the master nodes; the bearing friction is brought to the masters with the Guyan T."""
    ee = equivalent_elements(fem)
    n = len(fem.masters)
    M = np.diag([float(fem.J_nodes[i]) for i in fem.masters])
    K = np.zeros((n, n))
    for e in ee:
        ids = np.ix_([e.ia, e.ib], [e.ia, e.ib])
        K[ids] += stiffness_matrix(e.k)
        M[ids] += 0.5 * e.J * np.eye(2)
    T = guyan_transformation(fem.K, fem.masters)
    C_b = T.T @ fem.C_b @ T
    expected = (fem.i_encoder, fem.i_motor, fem.i_end, fem.i_P) if fem.name == "azimuth" else (fem.i_encoder, fem.i_H, fem.i_motor, fem.i_end, fem.i_P)
    if tuple(fem.masters) != expected:
        raise ValueError(f"axis {fem.name}: master nodes {fem.master_labels} do not form the expected chain {ROLES[fem.name]}")
    return LumpedChain(fem.name, ROLES[fem.name], M, K, 0.5 * (C_b + C_b.T), fem)
