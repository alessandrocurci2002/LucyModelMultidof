"""Inertia tensors and CAD reference frames."""

from __future__ import annotations

import numpy as np

_DIRECTIONS = {"X": (1.0, 0.0, 0.0), "+X": (1.0, 0.0, 0.0), "-X": (-1.0, 0.0, 0.0),
               "Y": (0.0, 1.0, 0.0), "+Y": (0.0, 1.0, 0.0), "-Y": (0.0, -1.0, 0.0),
               "Z": (0.0, 0.0, 1.0), "+Z": (0.0, 0.0, 1.0), "-Z": (0.0, 0.0, -1.0)}
Z_N = np.array([0.0, 0.0, 1.0])
"""Vertical of the inertial frame, up along the azimuth axis."""


def direction(name: str) -> np.ndarray:
    """``"X"``, ``"-Y"``, ... as a unit vector of N."""
    key = str(name).strip().upper()
    if key not in _DIRECTIONS:
        raise ValueError(f"direction '{name}' not recognised; allowed: {sorted(_DIRECTIONS)}")
    return np.array(_DIRECTIONS[key], dtype=float)


def R_y(theta: float) -> np.ndarray:
    """Elementary rotation about y, :math:`[[c, 0, s], [0, 1, 0], [-s, 0, c]]`."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def build_cad_frames(pointing: np.ndarray, polar_axis: np.ndarray, theta_cad: float) -> tuple[np.ndarray, np.ndarray]:
    """Attitudes of body 1 and of the pointing frame P in the CAD pose.

    ``pointing`` is the horizontal pointing direction (i_hat) and ``polar_axis`` the polar axis J_hat in N.
    With j_hat = z_hat x i_hat (it must be -J_hat), :math:`R_{N1} = [\\hat\\imath, \\hat\\jmath, \\hat z]` and
    :math:`R_{NP} = R_{N1} R_y(-\\theta_{cad})` (theta > 0 raises the pointing)."""
    x = np.asarray(pointing, dtype=float)
    x = x / np.linalg.norm(x)
    if abs(float(x @ Z_N)) > 1e-12:
        raise ValueError("the pointing direction of the CAD pose must be horizontal")
    y = np.cross(Z_N, x)
    y = y / np.linalg.norm(y)
    assert np.linalg.norm(np.cross(y, -np.asarray(polar_axis, dtype=float))) < 1e-9, (
        f"j_hat = {y} is not parallel to -J_hat = {-np.asarray(polar_axis)} (needs J_hat = K_hat x I_hat = i_hat x z_hat)")
    R1 = np.column_stack([x, y, Z_N])
    RP = R1 @ R_y(-theta_cad)
    assert np.allclose(R1.T @ R1, np.eye(3), atol=1e-12) and np.isclose(np.linalg.det(R1), 1.0, atol=1e-12)
    assert np.allclose(RP.T @ RP, np.eye(3), atol=1e-12)
    return R1, RP


def normalize_products(I: np.ndarray, convention: str) -> np.ndarray:
    """Tensor convention for the products of inertia: ``tensor`` (off-diagonal terms are -int(x_i x_j dm), unchanged)
    or ``positive`` (+int(x_i x_j dm), negated)."""
    I = np.asarray(I, dtype=float)
    if convention == "tensor":
        return I.copy()
    if convention == "positive":
        J = I.copy()
        off = ~np.eye(3, dtype=bool)
        J[off] = -J[off]
        return J
    raise ValueError(f"products convention '{convention}' not recognised")


def transport(I_G: np.ndarray, m: float, r: np.ndarray) -> np.ndarray:
    """Huygens-Steiner, :math:`I_P = I_G + m(\\|r\\|^2 1 - r r^T)` with :math:`r = r_G - r_P`."""
    I_G, r = np.asarray(I_G, dtype=float), np.asarray(r, dtype=float)
    return I_G + m * (float(r @ r) * np.eye(3) - np.outer(r, r))


def rotate(I: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Change of frame of a tensor, :math:`I|_2 = R^T I|_1 R` with :math:`R = R_{12}`."""
    I, R = np.asarray(I, dtype=float), np.asarray(R, dtype=float)
    return R.T @ I @ R


def check_physical(I: np.ndarray, label: str, tol: float = 1e-9) -> np.ndarray:
    """Symmetry, positive definiteness and triangle inequality of the principal moments; returns them sorted."""
    I = np.asarray(I, dtype=float)
    scale = max(float(np.max(np.abs(I))), 1e-30)
    if float(np.max(np.abs(I - I.T))) / scale > tol:
        raise ValueError(f"{label}: inertia tensor not symmetric")
    p = np.sort(np.linalg.eigvalsh(I))
    if p[0] <= 0.0:
        raise ValueError(f"{label}: inertia tensor not positive definite, principal moments {p}")
    if (p[0] + p[1] - p[2]) / p[2] < -tol:
        raise ValueError(f"{label}: triangle inequality violated, principal moments {p}")
    return p
