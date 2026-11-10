r"""Properties of the shaft elements from coaxial cylinders.

For an annular layer of length :math:`L`, diameters :math:`D_o`, :math:`D_i`, density :math:`\rho`, shear modulus
:math:`G`:

.. math::
    m = \rho\frac{\pi}{4}(D_o^2 - D_i^2)L,\qquad J = \rho L J_p,\qquad J_p = \frac{\pi}{32}(D_o^4 - D_i^4),
    \qquad I_t = m\left[\frac{D_o^2 + D_i^2}{16} + \frac{L^2}{12}\right].

The coaxial layers of a segment add their inertias and, if they carry torque, their stiffnesses
:math:`k_s = \sum G J_p / L_s`; the segments of an element add their inertias and their compliances in series,
:math:`1/k = \sum 1/k_s + 1/k_{joint}`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import Element


@dataclass(frozen=True)
class Piece:
    """Rigid axisymmetric piece (for the transverse inertia J^ax); ``x_c`` from node ``i`` of its element."""

    name: str
    m: float
    J: float
    I_t: float
    x_c: float


@dataclass(frozen=True)
class ElementProperties:
    """Inertia and torsional stiffness of one element, SI."""

    label: str
    mode: str
    length: float
    J_distributed: float
    J_extra: float
    """Concentrated inertias of the element, split half and half on its nodes."""
    k: float
    m: float
    I_t: float
    pieces: tuple[Piece, ...] = ()

    @property
    def J(self) -> float:
        return self.J_distributed + self.J_extra


def layer_mass(rho: float, d_out: float, d_in: float, L: float) -> float:
    return rho * np.pi / 4.0 * (d_out**2 - d_in**2) * L


def polar_area_moment(d_out: float, d_in: float) -> float:
    return np.pi / 32.0 * (d_out**4 - d_in**4)


def transverse_inertia(m: float, d_out: float, d_in: float, L: float) -> float:
    return m * ((d_out**2 + d_in**2) / 16.0 + L**2 / 12.0)


def _series(compliances) -> float:
    c = float(sum(compliances))
    if c <= 0.0:
        raise ValueError("zero total compliance: the element would be infinitely stiff")
    return 1.0 / c


def element_properties(el: Element, label: str | None = None) -> ElementProperties:
    """Inertia, stiffness and mass of an element (``geometry``, ``manual`` or ``telescope`` mode)."""
    name = label or str(el.id)
    if el.mode == "telescope":
        # no shaft: zero inertia (the telescope is in the payload), rigid connection
        return ElementProperties(name, "telescope", el.length, 0.0, 0.0, float("inf"), 0.0, 0.0)
    pieces: list[Piece] = []
    k_seg: list[float] = []
    if el.mode == "geometry":
        x = 0.0
        for s, seg in enumerate(el.segments, start=1):
            GJ = 0.0
            for j, ly in enumerate(seg.layers, start=1):
                rho, G = ly.material.rho, ly.material.G
                m = layer_mass(rho, ly.d_out, ly.d_in, seg.length)
                Jp = polar_area_moment(ly.d_out, ly.d_in)
                J = rho * seg.length * Jp
                GJ += G * Jp if ly.carries_torque else 0.0
                pieces.append(Piece(f"el.{name} seg.{s} layer {j}", m, J, transverse_inertia(m, ly.d_out, ly.d_in, seg.length),
                                    x + seg.length / 2.0))
            if GJ <= 0.0 and el.k_manual is None:
                raise ValueError(f"element {name}, segment {s}: no layer carries torque")
            if GJ > 0.0:
                k_seg.append(GJ / seg.length)
            x += seg.length
        J_dist = float(sum(p.J for p in pieces))
        compliances = [1.0 / el.k_manual] if el.k_manual is not None else [1.0 / k for k in k_seg]
    else:
        assert el.J_manual is not None and el.k_manual is not None
        J_dist = float(el.J_manual)
        compliances = [1.0 / el.k_manual]
        if None not in (el.m_manual, el.I_t_manual, el.y_c_manual):
            pieces.append(Piece(f"el.{name} (manual)", float(el.m_manual), J_dist, float(el.I_t_manual), float(el.y_c_manual)))  # type: ignore[arg-type]
        else:
            pieces.append(Piece(f"el.{name} (manual, without m/I_t)", 0.0, J_dist, 0.0, el.length / 2.0))
    if el.k_joint is not None:
        compliances.append(1.0 / el.k_joint)
    k = _series(compliances)
    J_extra = float(sum(e.J for e in el.extra))
    for e in el.extra:
        pieces.append(Piece(f"el.{name} extra '{e.name}'", e.m or 0.0, e.J, e.I_t or 0.0, el.length / 2.0))
    m_tot = float(sum(p.m for p in pieces))
    x_G = float(sum(p.m * p.x_c for p in pieces) / m_tot) if m_tot > 0.0 else el.length / 2.0
    I_t = float(sum(p.I_t + p.m * (p.x_c - x_G) ** 2 for p in pieces))
    return ElementProperties(name, el.mode, el.length, J_dist, J_extra, k, m_tot, I_t, tuple(pieces))


def spring(label: str, mode: str, length: float, k: float) -> ElementProperties:
    """Massless spring element (half hub, telescope attachments)."""
    return ElementProperties(label, mode, length, 0.0, 0.0, k, 0.0, 0.0)
