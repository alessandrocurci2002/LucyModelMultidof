"""Short co-simulations on the test circle against the original thesis code, and one full case of the thesis campaign."""

import json

import numpy as np
import pytest

from conftest import REF, rel
from lucy_sim import Disturbance
from lucy_sim.pointing.geometry import Circle

CASES = {"cascade_nominal": ("cascade", Disturbance(), 1), "lqg_torque": ("lqg", Disturbance(gauss=0.3, band=(0.5, 5.0)), 10),
         "lqg_friction": ("lqg", Disturbance(friction_factor=1.5), 1)}


@pytest.mark.parametrize("tag", CASES)
def test_circle_4s(model, tag):
    ref = np.load(REF / "circle_4s.npz")
    controller, d, seed = CASES[tag]
    t = model.time(4.0)
    R = model.simulate(Circle.reference(t), t, controller, d, seed)
    for ax in ("azimuth", "polar"):
        for k in ("los", "i", "p"):
            assert rel(R[ax][k][::10], ref[f"{tag}__{ax}__{k}"]) < 1e-6, (ax, k)


@pytest.mark.slow
def test_thesis_campaign_cases(model):
    """Four cases of the thesis campaign (circle, dataset7 with EKF2, dataset8 with config_v3; both controllers;
    nominal, torque, friction, combined): the metrics of the thesis."""
    from lucy_sim.campaign import Case, run_case

    for r in json.loads((REF / "thesis_hil_cases.json").read_text()):
        traj = "circle" if r["traj"] == "circle" else tuple(r["traj"].split("_", 1))
        q, _ = run_case(model, Case(traj, r["controller"], r["scenario"], r["seed"]))
        for k in ("E_rms_mrad", "miss_rms_cm", "on_target_pct", "i_max_az_A", "i_max_el_A", "loop_E_rms_mrad"):
            assert q[k] == pytest.approx(r[k], rel=1e-6, abs=1e-9), (r["traj"], k)
