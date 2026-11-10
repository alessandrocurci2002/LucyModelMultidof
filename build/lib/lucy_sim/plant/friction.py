r"""Dahl friction of the bearings, on the bearing deformation :math:`x_b = v^T x`.

.. math::
    \dot z = \dot x_b - \frac{\sigma_0|\dot x_b|}{g(\dot x_b)}\,z,\qquad
    g(v) = T_c\left[1 + (T_s/T_c - 1)\,e^{-(v/v_s)^2}\right],\qquad T_b = \sigma_0 z.

Pure Dahl is :math:`T_s = T_c` (``Ts_ratio = 1``): :math:`\sigma_0` is the pre-sliding stiffness, :math:`T_c` the
Coulomb level. The simulators update :math:`z` exactly at constant velocity over each step:
:math:`z \leftarrow z e^{-\alpha h} + \dot x_b(1 - e^{-\alpha h})/\alpha`, :math:`\alpha = \sigma_0|\dot x_b|/g`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DahlFriction:
    sigma0: float
    """Pre-sliding (bristle) stiffness, N m/rad."""
    Tc: float
    """Coulomb torque, N m."""
    Ts_ratio: float = 1.0
    vs: float = 1.0

    def g(self, v: float) -> float:
        return self.Tc * (1 + (self.Ts_ratio - 1) * np.exp(-(v / self.vs) ** 2))
