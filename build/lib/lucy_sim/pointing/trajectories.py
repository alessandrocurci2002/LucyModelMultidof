"""Recorded drone trajectories (dataset7 and dataset8 of the thesis campaign).

Files in ``data/trajectories/<dataset>/`` (columns ``t x y z qx qy qz qw``): ``truth.txt`` is the TagSLAM ground truth
levelled (rigidly rotated with z up), ``EKF2.txt`` and ``config_v3.txt`` the first run of the two estimators (PX4 EKF2
and OpenVINS with ``config_v3``). The estimates are aligned to the truth with 4 degrees of freedom (position and yaw),
as ``align_posyaw`` of ov_eval, on the common window, which starts and ends with the OpenVINS run. LUCY is placed
100 m from the first true position, 5 m below it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from ..paths import ROOT
from .geometry import angles

DISTANCE = 100.0
HEIGHT = 5.0
"""Distance of LUCY from the first point of the trajectory and height of the drone above LUCY, m (as the circle)."""
ESTIMATORS = ("EKF2", "config_v3")
DATASETS = ("dataset7", "dataset8")


def load_poses(f: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.loadtxt(f, comments="#")
    return d[:, 0], d[:, 1:4], d[:, 4:8]


def align_posyaw(pe: np.ndarray, pg: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Rotation about z and translation that align the estimate ``pe`` to the truth ``pg`` (least squares)."""
    me, mg = pe.mean(0), pg.mean(0)
    rc = ((pg - mg).T @ (pe - me)).T
    th = np.arctan2(rc[0, 1] - rc[1, 0], rc[0, 0] + rc[1, 1])
    c, s = np.cos(th), np.sin(th)
    Rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return Rz, mg - Rz @ me, float(th)


@dataclass
class Estimate:
    t: np.ndarray
    """Times from the start of the window, s."""
    p: np.ndarray
    """Aligned positions, frame with LUCY at the origin, m."""
    ate_m: float
    yaw_deg: float
    rate_Hz: float


@dataclass
class Trajectory:
    dataset: str
    t0: float
    duration: float
    truth: Callable[[np.ndarray], np.ndarray]
    """True position at times from the start of the window (frame with LUCY at the origin)."""
    estimates: dict[str, Estimate]
    origin: np.ndarray
    """Position of LUCY in the frame of the truth, m."""

    def reference(self, estimator: str, t: np.ndarray) -> dict[str, np.ndarray]:
        """Line of sight commanded by the estimator: the latest estimate, held between two estimates."""
        s = self.estimates[estimator]
        psi, th = angles(s.p)
        j = np.clip(np.searchsorted(s.t, t, side="right") - 1, 0, len(psi) - 1)
        return {"azimuth": psi[j], "polar": th[j]}


def load_trajectory(dataset: str, root: str | Path = ROOT) -> Trajectory:
    folder = Path(root) / "data" / "trajectories" / dataset
    tv, pv, _ = load_poses(folder / "truth.txt")
    t_ov = np.loadtxt(folder / "config_v3.txt", comments="#", usecols=0)
    t0, t1 = float(t_ov[0]), float(min(t_ov[-1], tv[-1]))
    origin = pv[0] - np.array([np.sqrt(DISTANCE ** 2 - HEIGHT ** 2), 0.0, HEIGHT])

    def interp(t):
        return np.stack([np.interp(t, tv, pv[:, i]) for i in range(3)], axis=1)
    est = {}
    for name in ESTIMATORS:
        te, pe, _ = load_poses(folder / f"{name}.txt")
        k = (te >= t0) & (te <= t1)
        te, pe = te[k], pe[k]
        pg = interp(te)
        Rz, tz, th = align_posyaw(pe, pg)
        pa = pe @ Rz.T + tz
        err = pa - pg
        est[name] = Estimate(te - t0, pa - origin, float(np.sqrt(np.mean(np.sum(err ** 2, axis=1)))), float(np.degrees(th)),
                             float(1 / np.median(np.diff(te))))
    return Trajectory(dataset, t0, t1 - t0, lambda t: interp(t + t0) - origin, est, origin)
