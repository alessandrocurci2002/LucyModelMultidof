r"""Telescope (payload, body 2) and the transverse inertia :math:`J^{ax}` of the polar parts about the azimuth axis.

The Onshape tensor (about the centre of mass, axes parallel to N) is brought to the tensor convention, transported
to :math:`o_2` (Huygens-Steiner) and rotated to the internal pointing frame P (:math:`x_P = \hat K` pointing,
:math:`y_P = \hat\jmath`, :math:`z_P` up at theta = 0). ``Payload.I`` and ``Payload.d`` are in P.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import GimbalConfig
from .fem import AxisFEM
from .hub import telescope_section_inertia
from .inertia import check_physical, normalize_products, rotate, transport


@dataclass(frozen=True)
class Payload:
    """Rigid telescope at node P, and transverse polar contribution J^ax about the azimuth axis."""

    I: np.ndarray
    """Tensor of the telescope about o_2, pointing frame P, kg m^2."""
    m: float
    d: np.ndarray
    """Offset of the centre of mass from o_2, frame P, m."""
    g: float
    J_ax: float

    I_xx = property(lambda self: float(self.I[0, 0]))
    I_yy = property(lambda self: float(self.I[1, 1]))
    I_zz = property(lambda self: float(self.I[2, 2]))
    I_xy = property(lambda self: float(self.I[0, 1]))
    I_xz = property(lambda self: float(self.I[0, 2]))
    I_yz = property(lambda self: float(self.I[1, 2]))
    d_x = property(lambda self: float(self.d[0]))
    d_y = property(lambda self: float(self.d[1]))
    d_z = property(lambda self: float(self.d[2]))

    def J_psi(self, theta: float | np.ndarray) -> float | np.ndarray:
        """:math:`J_{\\psi,p}(\\theta) = I_{xx}\\sin^2\\theta + I_{zz}\\cos^2\\theta + I_{xz}\\sin 2\\theta` (frame P)."""
        return self.I_xx * np.sin(theta) ** 2 + self.I_zz * np.cos(theta) ** 2 + self.I_xz * np.sin(2 * theta)


def build_payload(cfg: GimbalConfig, polar: AxisFEM) -> Payload:
    """Telescope at node P from the ``payload`` section, with the tube section clamped in the ring moved to the hub."""
    pc, fr, ax = cfg.payload, cfg.frames, cfg.polar
    y1 = ax.y_node1_from_o2
    J_ax = float(sum(pp.piece.I_t + pp.piece.m * (0.0 if y1 is None else y1 + ax.nodes_direction * pp.x) ** 2
                     for pp in polar.pieces))
    if any(n.includes_polar_shafts for n in cfg.azimuth.node_inertias):
        # the CAD value of the azimuth load already contains the polar shafts about the azimuth axis
        J_ax = 0.0
    hub = ax.hub
    if hub is not None and not hub.in_azimuth_load:
        J_ax += hub.J_azimuth
    m = float(pc.m)
    I_in = normalize_products(pc.I_G_N, pc.products_sign)
    check_physical(I_in, "payload tensor about the centre of mass")
    r = pc.r_G_N - fr.r_o2
    I = rotate(transport(I_in, m, r), fr.R_NP_0)
    d = fr.R_NP_0.T @ r
    check_physical(I, "payload tensor about o_2")
    section = telescope_section_inertia(hub) if hub is not None else None
    if section is not None:
        # tube clamped in the ring, centred in o_2 along the pointing axis: no Steiner term
        m_t, J_ax_t, J_t = section
        I = I - np.diag([J_ax_t, J_t, J_t])
        d = m * np.asarray(d, dtype=float) / (m - m_t)
        m = m - m_t
    return Payload(I=I, m=m, d=np.asarray(d, dtype=float), g=cfg.g, J_ax=J_ax)
