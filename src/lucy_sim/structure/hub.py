r"""Polar hub and compliant telescope attachment.

Layout::

    (tube + flange) - (telescope support ring) - (flange + tube)

The springs are in series. Tube and flange are shafts in torsion, :math:`k = G J_p / L`. The ring, loaded by
torques about :math:`\hat J` at the two attachment points, works as two half rings in parallel; along each the
moment is constant and splits into out-of-plane bending and wall torsion. With the complementary energy
(Castigliano)

.. math::
    k_{ring} = \frac{4}{\pi R}\left(\frac{1}{E I_r} + \frac{1}{G J_t}\right)^{-1},

with :math:`R` the mean radius, :math:`I_r = t b^3/12`, :math:`J_t = \beta b t^3`, :math:`t` the radial thickness
and :math:`b` the width (ring alone, without the telescope tube). ``stiffness_factor`` multiplies the result.

The hub centre H carries the inertia of the hub about the polar axis (Onshape export). Each half of the hub,
(tube + flange + half ring), is a massless spring between a bearing node and H. The telescope is a node P off
the chain, attached to H by :math:`k_b`.
"""

from __future__ import annotations

import numpy as np

from .config import AzimuthAttachment, Frames, Hub, HubSection
from .inertia import rotate, transport
from .shafts import ElementProperties, spring


def _J_polar(d_out: float, d_in: float) -> float:
    return np.pi * (d_out**4 - d_in**4) / 32.0


def section_stiffness(t: HubSection) -> float:
    """Torsion of a stretch coaxial with the polar axis, :math:`G J_p / L`."""
    return t.material.G * _J_polar(t.d_out, t.d_in) / t.length


def half_hub_stiffness(hub: Hub) -> float:
    """Half hub: tube + flange + half ring (the compliance of the whole ring split in two)."""
    ring = hub.ring
    k_tube, k_flange = section_stiffness(hub.tube), section_stiffness(hub.flange)
    t = (ring.d_out - ring.d_in) / 2.0
    R = (ring.d_out + ring.d_in) / 4.0
    b = ring.width
    E = 2.0 * ring.material.G * (1.0 + ring.poisson)
    EI_r = E * t * b**3 / 12.0
    beta = (1.0 - 0.63 * t / b + 0.052 * (t / b) ** 5) / 3.0
    GJ_t = ring.material.G * beta * b * t**3
    k_half_rings = 4.0 / (np.pi * R * (1.0 / EI_r + 1.0 / GJ_t))
    k_ring = ring.stiffness_factor * k_half_rings
    return 1.0 / (1.0 / k_tube + 1.0 / k_flange + 1.0 / (2.0 * k_ring))


def hub_inertia(hub: Hub, frames: Frames) -> tuple[float, float]:
    """``(J_polar, J_azimuth)``: the Onshape tensor transported to o_2 and rotated to the pointing frame P;
    about the polar axis and about z_hat at theta = 0."""
    I_P = rotate(transport(hub.I_G_N, hub.m, hub.r_G_N - frames.r_o2), frames.R_NP_0)
    return float(I_P[1, 1]), float(I_P[2, 2])


def telescope_section_inertia(hub: Hub) -> tuple[float, float, float] | None:
    """``(m, J_axial, J_transverse)`` of the telescope tube clamped in the ring (``None`` if absent).

    Tube along the pointing axis, centred in o_2: :math:`J_{ax} = m(r_o^2 + r_i^2)/2`, and about a diameter (the
    polar axis, and z_hat at theta = 0) :math:`J_t = m[(r_o^2 + r_i^2)/4 + L^2/12]`."""
    t = hub.telescope_section
    if t is None:
        return None
    ro, ri = t.d_out / 2.0, t.d_in / 2.0
    m = t.material.rho * np.pi * (ro**2 - ri**2) * t.length
    return m, m * (ro**2 + ri**2) / 2.0, m * ((ro**2 + ri**2) / 4.0 + t.length**2 / 12.0)


def half_hub_element(hub: Hub, label: str) -> ElementProperties:
    """Half hub as a massless spring (all the hub inertia is at node H)."""
    length = hub.tube.length + hub.flange.length + hub.ring.d_out / 2.0
    return spring(label, "hub", length, half_hub_stiffness(hub))


def polar_attachment_element(hub: Hub) -> ElementProperties:
    """Spring k_b between the hub centre H and the telescope P."""
    return spring("b", "attachment", 0.0, hub.k_b)


def azimuth_attachment_element(att: AzimuthAttachment) -> ElementProperties:
    """Spring between the azimuth load L (fork + hub + ring) and the telescope P."""
    return spring("b", "attachment", 0.0, att.k)
