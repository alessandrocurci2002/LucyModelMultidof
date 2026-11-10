r"""Torsional FEM of each axis, isolated.

Two-node elements with :math:`K_e = k_e[1, -1; -1, 1]` and lumped mass :math:`M_e = \frac{J_e}{2}\mathbb 1`. In the
model with the hub the polar telescope element has no shaft: its two nodes are one node, the hub centre H, and the
telescope is a node P off the chain attached to H by :math:`k_b`; on the azimuth axis the telescope node P hangs from
the load node L. P carries no inertia here: the telescope enters through the payload of the Lagrange model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla

from .config import AxisConfig
from .hub import azimuth_attachment_element, half_hub_element, polar_attachment_element
from .shafts import ElementProperties, Piece, element_properties


@dataclass(frozen=True)
class PlacedPiece:
    """Rigid piece of the axis with the axial position of its centre of mass from node 1."""

    piece: Piece
    x: float


@dataclass(frozen=True)
class AxisFEM:
    """FEM of one axis (without the telescope inertia).

    ``labels`` are the node labels in internal order (``"1"``, ..., ``"H"``, ..., ``"P"``); ``links`` the node
    pairs of ``elements``; ``masters`` the internal indices of the master nodes."""

    name: str
    labels: tuple[str, ...]
    elements: tuple[ElementProperties, ...]
    links: tuple[tuple[int, int], ...]
    K: np.ndarray
    M: np.ndarray
    """Lumped mass matrix (shafts and node inertias)."""
    C_b: np.ndarray
    J_nodes: np.ndarray
    i_encoder: int
    i_motor: int
    i_end: int
    i_P: int
    i_H: int
    """Node the telescope P hangs from (hub centre H, or the azimuth load L)."""
    masters: tuple[int, ...]
    pieces: tuple[PlacedPiece, ...]

    @property
    def J_total(self) -> float:
        return float(sum(p.J for p in self.elements) + self.J_nodes.sum())

    def label(self, i: int) -> str:
        return "P" if i == self.i_P else self.labels[i]

    @property
    def master_labels(self) -> tuple[str, ...]:
        return tuple(self.label(i) for i in self.masters)


def stiffness_matrix(k: float) -> np.ndarray:
    return k * np.array([[1.0, -1.0], [-1.0, 1.0]])


def build_axis(axis: AxisConfig) -> AxisFEM:
    """Assembles the FEM of one axis with the telescope node P off the chain."""
    hub = axis.hub
    originals = tuple(half_hub_element(hub, str(e.id)) if hub is not None and e.id in hub.side_elements
                      else element_properties(e) for e in axis.elements)
    elements: list[ElementProperties] = []
    labels, x_nodes, links = ["1"], [0.0], []
    extra_nodal: dict[str, float] = {}
    alias: dict[str, int] = {}
    i_P: int | None = None
    xu = 0.0
    for el, p in zip(axis.elements, originals):
        i = len(labels) - 1
        if axis.name == "polar" and el.id == axis.payload_element:
            # no shaft: the rigid telescope makes nodes i and i + 1 a single node, the hub centre H
            assert axis.payload_xi is not None
            i_P = i
            alias[str(el.id)] = alias[str(el.id + 1)] = i
            labels[i] = "H"
            x_nodes[i] = xu + axis.payload_xi * el.length
            xu += el.length
            continue
        labels.append(str(el.id + 1))
        xu += p.length
        x_nodes.append(xu)
        links.append((i, i + 1))
        elements.append(p)
    if axis.name == "polar":
        assert hub is not None and i_P is not None
        i_H = i_P
        extra_nodal[str(hub.center_element)] = extra_nodal.get(str(hub.center_element), 0.0) + hub.J_polar
        attachment = polar_attachment_element(hub)
    else:
        assert axis.attachment is not None
        i_H = labels.index(str(axis.end_node))
        attachment = azimuth_attachment_element(axis.attachment)
    labels.append("P")
    x_nodes.append(x_nodes[i_H])
    i_P = len(labels) - 1
    links.append((i_H, i_P))
    elements.append(attachment)

    n = len(labels)
    idx = {e: k for k, e in enumerate(labels)}
    idx.update(alias)
    K, M = np.zeros((n, n)), np.zeros((n, n))
    for p, (a, b) in zip(elements, links):
        ids = np.ix_([a, b], [a, b])
        K[ids] += stiffness_matrix(p.k)
        M[ids] += 0.5 * p.J * np.eye(2)
    J_nodes = np.zeros(n)
    for ni in axis.node_inertias:
        J_nodes[idx[str(ni.node)]] += ni.J
    for node, J in extra_nodal.items():
        J_nodes[idx[node]] += J
    M += np.diag(J_nodes)
    C_b = np.zeros((n, n))
    for el_id, b in axis.bearing_b:
        a_i, b_i = idx[str(el_id)], idx[str(el_id + 1)]
        C_b[a_i, a_i] += b / 2.0
        C_b[b_i, b_i] += b / 2.0
    masters = sorted({(i_P if str(m) == "P" else idx[str(m)]) for m in axis.masters} | {i_H, i_P})

    # rigid pieces placed along the axis (for J^ax)
    x_user = {}
    xc = 0.0
    for k, p in enumerate(originals):
        x_user[str(k + 1)] = xc
        xc += p.length
    x_user[str(axis.n_nodes)] = xc
    pieces = [PlacedPiece(pz, x_user[str(k + 1)] + pz.x_c) for k, p in enumerate(originals) for pz in p.pieces]
    pieces += [PlacedPiece(Piece(f"node {ni.node} '{ni.name}'", ni.m or 0.0, ni.J, ni.I_t or 0.0, 0.0), x_user[str(ni.node)])
               for ni in axis.node_inertias]

    fem = AxisFEM(axis.name, tuple(labels), tuple(elements), tuple(links), K, M, C_b, J_nodes, idx[str(axis.encoder_node)],
                  idx[str(axis.torque_node)], idx[str(axis.end_node)], i_P, i_H, tuple(masters), tuple(pieces))
    check_matrices(fem.M, fem.K, fem.J_total, label=axis.name, massless=(i_P,))
    return fem


def check_matrices(M: np.ndarray, K: np.ndarray, J_total: float, *, label: str = "", massless: tuple[int, ...] = ()) -> None:
    """Symmetry, positive definiteness of M (except the massless nodes), :math:`K 1 = 0`, total inertia."""
    scale_K = max(float(np.max(np.abs(K))), 1e-30)
    assert np.allclose(M, M.T, atol=0.0, rtol=1e-14), f"{label}: M not symmetric"
    assert np.allclose(K, K.T, atol=1e-14 * scale_K), f"{label}: K not symmetric"
    others = [i for i in range(M.shape[0]) if i not in massless]
    assert np.all(np.linalg.eigvalsh(M[np.ix_(others, others)]) > 0.0), f"{label}: M not positive definite"
    one = np.ones(K.shape[0])
    assert np.max(np.abs(K @ one)) <= 1e-10 * scale_K, f"{label}: K 1 != 0"
    assert abs(float(one @ M @ one) - J_total) <= 1e-12 * max(J_total, 1e-30), f"{label}: 1^T M 1 differs from the total inertia"


def natural_frequencies(M: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Natural frequencies in Hz, signed (negative for negative eigenvalues)."""
    w2 = sla.eigh(K, M)[0]
    scale = max(float(np.max(np.abs(w2))), 1.0)
    w2 = np.where(np.abs(w2) < 1e-12 * scale, 0.0, w2)
    return np.sign(w2) * np.sqrt(np.abs(w2)) / (2.0 * np.pi)


def stiffness_proportional_beta(M: np.ndarray, K: np.ndarray, zeta: float) -> float:
    """:math:`\\beta = 2\\zeta/\\omega_1`, with :math:`\\omega_1` the first flexible mode (C = beta K)."""
    w = natural_frequencies(M, K) * 2.0 * np.pi
    flex = w[w > 1e-6 * max(w.max(), 1.0)]
    return 0.0 if flex.size == 0 else 2.0 * zeta / float(flex[0])
