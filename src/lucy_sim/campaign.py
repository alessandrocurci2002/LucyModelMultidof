"""Simulated campaign of the thesis: simulated pointer (coupled non-linear plant with Dahl bearings) with the cascade and
the tuned LQG-PI, on the test circle and on the recorded drone trajectories.

Trajectories: ``circle`` (Attorrese) and ``(dataset, estimator)`` for dataset7, dataset8 and the estimators EKF2 and
OpenVINS ``config_v3`` (first run): the reference of the loop is the estimate, the error is computed on the TagSLAM
ground truth.

Scenarios (:data:`SCENARIOS`):
    0  nominal: nominal plant, no disturbance;
    1  zero-mean Gaussian disturbance torque on both axes, band 0.5-5 Hz, 0.3 N m rms (generic equivalent of the gusts);
    2  Coulomb friction of the plant x1.5 on both axes (model error: the controllers do not change);
    3  1 and 2 together.
The stochastic scenarios use three realisations, the same for every trajectory and controller (seed + 1 on the polar axis).
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .model import GimbalModel
from .paths import ROOT
from .pointing.geometry import BEAM_RADIUS, T_TRANSIENT, Circle, line_of_sight, pointing_error, statistics
from .pointing.trajectories import DATASETS, DISTANCE, HEIGHT, Trajectory, load_trajectory
from .simulation.disturbances import Disturbance

SCENARIOS = {0: Disturbance(), 1: Disturbance(gauss=0.3, band=(0.5, 5.0)), 2: Disturbance(friction_factor=1.5),
             3: Disturbance(gauss=0.3, band=(0.5, 5.0), friction_factor=1.5)}
SCENARIO_NAMES = {0: "nominal", 1: "torque", 2: "friction", 3: "combined"}
SEEDS = (10, 20, 30)
NOMINAL_SEED = 1
"""Seed of the deterministic scenarios (it does not matter: no random input)."""
TRAJECTORIES = ("circle", ("dataset7", "config_v3"), ("dataset8", "config_v3"), ("dataset7", "EKF2"), ("dataset8", "EKF2"))
CONTROLLERS = ("cascade", "lqg")
DECIMATION = 10
"""Saved series and metrics every DECIMATION simulation steps (1 ms)."""


@dataclass(frozen=True)
class Case:
    trajectory: str | tuple[str, str]
    controller: str
    scenario: int
    seed: int | None = None

    @property
    def tag(self) -> str:
        return self.trajectory if isinstance(self.trajectory, str) else f"{self.trajectory[0]}_{self.trajectory[1]}"

    @property
    def name(self) -> str:
        return f"{self.tag}__{self.controller}__s{self.scenario}" + ("" if self.seed is None else f"_{self.seed}")


def thesis_cases() -> list[Case]:
    stochastic = {k for k, d in SCENARIOS.items() if d.gauss or d.noise}
    return [Case(tr, c, s, seed) for tr in TRAJECTORIES[::-1] for c in CONTROLLERS for s in SCENARIOS
            for seed in (SEEDS if s in stochastic else (None,))]


def run_case(model: GimbalModel, case: Case, trajectories: dict[str, Trajectory] | None = None, root: str | Path = ROOT,
             disturbance: Disturbance | None = None) -> tuple[dict, dict[str, np.ndarray]]:
    """One simulation: ``(metrics, series)``. ``disturbance`` overrides the scenario of the case."""
    if case.trajectory == "circle":
        t = model.time(Circle.duration())
        ref = Circle.reference(t)
        truth = Circle.truth
    else:
        ds, est = case.trajectory
        tr = (trajectories or {}).get(ds) or load_trajectory(ds, root)
        t = model.time(tr.duration)
        ref = tr.reference(est, t)
        truth = tr.truth
    d = SCENARIOS[case.scenario] if disturbance is None else disturbance
    R = model.simulate(ref, t, case.controller, d, NOMINAL_SEED if case.seed is None else case.seed)
    k = np.arange(0, len(t), DECIMATION)
    tk = t[k]
    pv = truth(tk)
    n_los = line_of_sight(ref["azimuth"][0] + R["azimuth"]["los"][k], ref["polar"][0] + R["polar"]["los"][k])
    n_ref = line_of_sight(ref["azimuth"][k], ref["polar"][k])
    E, miss = pointing_error(n_los, pv)
    Es, ms = pointing_error(n_ref, pv)
    Et = np.linalg.norm(np.cross(n_los, n_ref), axis=1)
    m = tk >= T_TRANSIENT
    q = {"traj": case.tag, "controller": case.controller, "scenario": case.scenario, "seed": case.seed} | statistics(tk, E, miss)
    q |= {"estimator_alone": statistics(tk, Es, ms), "loop_E_rms_mrad": 1e3 * float(np.sqrt(np.mean(Et[m] ** 2))),
          "loop_E_max_mrad": 1e3 * float(Et[m].max()), "duration_s": float(tk[-1])}
    for a, b in (("az", "azimuth"), ("el", "polar")):
        q[f"i_max_{a}_A"] = float(np.abs(R[b]["i"][k][m]).max())
        q[f"i_rms_{a}_A"] = float(np.sqrt(np.mean(R[b]["i"][k][m] ** 2)))
        q[f"los_{a}_rms_mrad"] = 1e3 * float(np.sqrt(np.mean((ref[b][k] - ref[b][0] - R[b]["los"][k])[m] ** 2)))
    series = {"t": tk, "E": E, "miss": miss, "miss_est": ms, "los_az": R["azimuth"]["los"][k], "los_el": R["polar"]["los"][k],
              "ref_az": ref["azimuth"][k] - ref["azimuth"][0], "ref_el": ref["polar"][k] - ref["polar"][0],
              "i_az": R["azimuth"]["i"][k], "i_el": R["polar"]["i"][k]}
    return q, series


_MODEL: GimbalModel | None = None
_TRAJ: dict[str, Trajectory] = {}


def _init(root: str) -> None:
    global _MODEL, _TRAJ
    _MODEL = GimbalModel.load(root)
    _TRAJ = {ds: load_trajectory(ds, root) for ds in DATASETS}


def _worker(args) -> dict:
    case, out = args
    q, s = run_case(_MODEL, case, _TRAJ)  # type: ignore[arg-type]
    np.savez_compressed(Path(out) / "series" / f"{case.name}.npz", **{k: v.astype(np.float32) for k, v in s.items()})
    print(case.name, "| E rms %.3f mrad | miss rms %.2f cm | on target %.1f %% | i max %.1f / %.1f A"
          % (q["E_rms_mrad"], q["miss_rms_cm"], q["on_target_pct"], q["i_max_az_A"], q["i_max_el_A"]), flush=True)
    return q


def run_campaign(out: str | Path, cases: list[Case] | None = None, workers: int = 3, root: str | Path = ROOT) -> dict:
    """Runs the cases (default: the whole thesis campaign, 80 simulations) in ``workers`` processes; writes
    ``summary.json`` and the series (float32, one per case) to ``out/series``."""
    out = Path(out)
    (out / "series").mkdir(parents=True, exist_ok=True)
    cases = thesis_cases() if cases is None else cases
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with ProcessPoolExecutor(workers, initializer=_init, initargs=(str(root),)) as ex:
        rows = list(ex.map(_worker, [(c, str(out)) for c in cases]))
    geo = {}
    for ds in DATASETS:
        tr = load_trajectory(ds, root)
        geo[ds] = {"origin_m": tr.origin.tolist(), "window_s": tr.duration,
                   "estimators": {n: {"ate_m": e.ate_m, "yaw_deg": e.yaw_deg, "rate_Hz": e.rate_Hz} for n, e in tr.estimates.items()}}
    summary = {"cases": rows, "geometry": geo, "scenarios": {str(k): v.to_dict() for k, v in SCENARIOS.items()}, "seeds": list(SEEDS),
               "beam_radius_m": BEAM_RADIUS, "distance_m": DISTANCE, "height_m": HEIGHT}
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary
