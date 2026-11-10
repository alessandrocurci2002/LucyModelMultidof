r"""Pointing geometry and metrics.

Pointing error (Attorrese 2026, eq. 5.20-5.23): :math:`E = \lVert\hat n(\psi,\theta)\times\hat r\rVert`, with
:math:`\hat n` the line of sight and :math:`\hat r` towards the *true* position of the drone; the miss distance is
:math:`\lVert r\rVert E`. The beam is on target when the miss distance is below its radius.
"""

from __future__ import annotations

import numpy as np

BEAM_RADIUS = 0.08
"""Radius of the beam at the drone, m."""
T_TRANSIENT = 5.0
"""Start excluded from the metrics, s."""


class Circle:
    """Test circle of Attorrese: radius 10 m, height 5 m, centre at 100 m, 0.2 rad/s, starting at rest."""

    centre, radius, height, omega = 100.0, 10.0, 5.0, 0.2
    t_start = T_TRANSIENT

    @classmethod
    def duration(cls) -> float:
        """Start plus one orbit."""
        return cls.t_start + 2 * np.pi / cls.omega

    @classmethod
    def reference(cls, t: np.ndarray) -> dict[str, np.ndarray]:
        x, y = cls.centre + cls.radius * np.cos(cls.omega * t), cls.radius * np.sin(cls.omega * t)
        return {"azimuth": np.unwrap(np.arctan2(y, x)), "polar": np.arctan2(cls.height, np.hypot(x, y))}

    @classmethod
    def truth(cls, t: np.ndarray) -> np.ndarray:
        """Position of the drone, frame with LUCY at the origin, m."""
        return np.stack([cls.centre + cls.radius * np.cos(cls.omega * t), cls.radius * np.sin(cls.omega * t), np.full_like(t, cls.height)], axis=1)


def angles(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Azimuth and elevation of the positions ``p`` (n x 3) seen from the origin."""
    return np.arctan2(p[:, 1], p[:, 0]), np.arctan2(p[:, 2], np.hypot(p[:, 0], p[:, 1]))


def line_of_sight(psi: np.ndarray, th: np.ndarray) -> np.ndarray:
    """Unit vector of the line of sight (eq. 5.22)."""
    return np.stack([np.cos(th) * np.cos(psi), np.cos(th) * np.sin(psi), np.sin(th)], axis=1)


def pointing_error(n: np.ndarray, p_true: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(E [rad], miss distance [m])`` of the lines of sight ``n`` towards the positions ``p_true`` (eq. 5.23)."""
    d = np.linalg.norm(p_true, axis=1)
    E = np.linalg.norm(np.cross(n, p_true / d[:, None]), axis=1)
    return E, E * d


def statistics(t: np.ndarray, E: np.ndarray, miss: np.ndarray, t0: float = T_TRANSIENT, beam: float = BEAM_RADIUS) -> dict[str, float]:
    """Pointing error and miss distance after ``t0``; ``on_target_pct`` is the effective pointing."""
    m = t >= t0
    return {"E_rms_mrad": 1e3 * float(np.sqrt(np.mean(E[m] ** 2))), "E_mean_mrad": 1e3 * float(E[m].mean()),
            "E_p95_mrad": 1e3 * float(np.percentile(E[m], 95)), "E_max_mrad": 1e3 * float(E[m].max()),
            "miss_rms_cm": 100 * float(np.sqrt(np.mean(miss[m] ** 2))), "miss_mean_cm": 100 * float(miss[m].mean()),
            "miss_p95_cm": 100 * float(np.percentile(miss[m], 95)), "miss_max_cm": 100 * float(miss[m].max()),
            "on_target_pct": 100 * float(np.mean(miss[m] <= beam))}
