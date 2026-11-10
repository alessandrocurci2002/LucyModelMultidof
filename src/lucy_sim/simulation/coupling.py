r"""Velocity-dependent torques exchanged by the two axes through the telescope, from the Lagrange equations.

With :math:`\omega = (\dot\psi\sin\theta, -\dot\theta, \dot\psi\cos\theta)` and the tensor :math:`I` about :math:`o_2`,
:math:`T = \tfrac12 J(\theta)\dot\psi^2 + \tfrac12 I_{yy}\dot\theta^2 + H(\theta)\dot\psi\dot\theta`, with
:math:`J = I_{xx}\sin^2\theta + I_{zz}\cos^2\theta + 2I_{xz}\sin\theta\cos\theta` and :math:`H = -(I_{xy}\sin\theta + I_{yz}\cos\theta)`.
On the right-hand side: azimuth :math:`-(J'\dot\psi\dot\theta + H'\dot\theta^2)`, elevation :math:`+\tfrac12 J'\dot\psi^2`.
The inertial terms :math:`-H\ddot\theta` on the azimuth and :math:`-H\ddot\psi` on the elevation (:math:`H\sim 10^{-6}`) are
added by the co-simulation with :meth:`LagrangeCoupling.H`.

``scale`` is the CAD factor identified on the azimuth load, which the plant already applies to :math:`J(0) = I_{zz}`; it is
applied to the whole tensor so that energy and torques stay consistent.
"""

from __future__ import annotations

import math

import numpy as np


class LagrangeCoupling:
    def __init__(self, payload, scale: float = 1.0):
        I = scale * payload.I
        self.Ixx, self.Izz = float(I[0, 0]), float(I[2, 2])
        self.dIxz, self.Ixz, self.Ixy, self.Iyz = float(I[0, 0] - I[2, 2]), float(I[0, 2]), float(I[0, 1]), float(I[1, 2])

    def H(self, th: float) -> float:
        """Off-diagonal entry of the mass matrix of the telescope."""
        return -(self.Ixy * math.sin(th) + self.Iyz * math.cos(th))

    def delta_J(self, th: float) -> float:
        """:math:`J(\\theta) - J(0)`: change of the telescope inertia about the azimuth axis from the identified one at theta = 0."""
        sn, cs = math.sin(th), math.cos(th)
        return self.Ixx * sn * sn + self.Izz * (cs * cs - 1.0) + 2 * self.Ixz * sn * cs

    def torques(self, th: float, dpsi: float, dth: float) -> tuple[float, float]:
        """``(azimuth, elevation)`` torques for elevation ``th`` and rates ``dpsi``, ``dth``."""
        sn, cs = math.sin(th), math.cos(th)
        dJ = self.dIxz * 2 * sn * cs + 2 * self.Ixz * (cs * cs - sn * sn)
        dH = -(self.Ixy * cs - self.Iyz * sn)
        return -(dJ * dpsi * dth + dH * dth * dth), 0.5 * dJ * dpsi * dpsi


def check_against_lagrange(coupled) -> dict[str, float]:
    """The closed forms against the symbolic Lagrange model (rigid motion of both axes, no dissipation, gravity
    removed); raises if they differ. ``coupled``: :class:`lucy_sim.structure.lagrange.CoupledModel`."""
    acc = coupled.without_dissipation()
    st = acc.structure
    cl = LagrangeCoupling(acc.payload)
    err = 0.0
    for th, dp, dt_ in ((0.05, 0.02, 0.001), (0.3, 0.2, 0.05), (0.8, 0.5, -0.2), (-0.6, -0.3, 0.1)):
        q = np.zeros(st.n)
        q[st.n_psi:] = th
        dq = np.zeros(st.n)
        dq[:st.n_psi], dq[st.n_psi:] = dp, dt_
        f = acc.forces(q, dq) - acc.forces(q, 0 * dq)
        ta, te = cl.torques(th, dp, dt_)
        err = max(err, abs(float(-f[:st.n_psi].sum()) - ta), abs(float(-f[st.n_psi:].sum()) - te))
    if err > 1e-9:
        raise AssertionError(f"coupling torques differ from the Lagrange model: {err:.3e}")
    j, z = st.i_psi_p, np.zeros(st.n)
    M0 = float(acc.M(z)[j, j])
    err_J = err_H = 0.0
    for th in (0.0, 0.05, 0.5, 1.0, -0.8):
        q = z.copy()
        q[st.j_theta_p] = th
        err_J = max(err_J, abs(float(acc.M(q)[j, j]) - M0 - cl.delta_J(th)))
        err_H = max(err_H, abs(float(acc.M(q)[j, st.j_theta_p]) - cl.H(th)))
    if err_J > 1e-9 or err_H > 1e-12:
        raise AssertionError(f"J(theta) or H(theta) differ from the Lagrange model: {err_J:.3e}, {err_H:.3e}")
    return {"torques_N_m": err, "J_kg_m2": err_J, "H_kg_m2": err_H}
