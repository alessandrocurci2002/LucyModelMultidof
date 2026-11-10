"""The model built from config/ against the thesis model (reference values computed with the original thesis code)."""

import numpy as np
import pytest

from conftest import rel
from lucy_sim.control import lqg
from lucy_sim.simulation.coupling import check_against_lagrange

AXES = ("azimuth", "polar")


@pytest.mark.parametrize("axis", AXES)
def test_structure_of_the_axes(model, ref_model, axis):
    m, r = model.axes[axis].model, ref_model[axis]
    assert rel(m.M, r["M"]) < 1e-12
    assert rel(m.K, r["K"]) < 1e-12
    assert m.Kt == pytest.approx(r["Kt"], rel=1e-14)
    assert rel(m.omega, r["omega"]) < 1e-9


@pytest.mark.parametrize("axis", AXES)
def test_closed_loops(model, ref_model, axis):
    a, r = model.axes[axis], ref_model[axis]
    assert np.trace(a.cascade.A) == pytest.approx(r["cascade_A_trace"], rel=1e-10)
    assert np.linalg.norm(a.velocity_loop.A) == pytest.approx(r["velocity_A_fro"], rel=1e-10)
    assert len(a.Ap) == r["n_design"]


@pytest.mark.parametrize("axis", AXES)
def test_lqg_gains(model, ref_model, axis):
    c, r = model.controllers[axis]["lqg"], ref_model[axis]
    assert isinstance(c, lqg.LQGPIMeasured)
    assert c.k_p == pytest.approx(r["k_p"], rel=1e-8)
    assert c.k_i == pytest.approx(r["k_i"], rel=1e-8)
    assert rel(c.K_e[:, 0], r["K_e"]) < 1e-8
    if r["K_f"]:
        assert rel(c.K_f, r["K_f"]) < 1e-8


def test_payload(model, ref_model):
    p, r = model.structure.payload, ref_model["payload"]
    assert rel(p.I, r["I"]) < 1e-12 and p.m == pytest.approx(r["m"]) and p.J_ax == pytest.approx(r["J_ax"])


def test_coupling_torques_match_lagrange(model):
    err = check_against_lagrange(model.structure.coupled)
    assert err["torques_N_m"] < 1e-9


@pytest.mark.parametrize("axis", AXES)
def test_tuned_loop_is_stable_with_bearings_sliding_and_stuck(model, axis):
    a, c = model.axes[axis], model.controllers[axis]["lqg"]
    an = a.velocity_loop
    for A in a.bearing_variants().values():
        assert np.abs(lqg.closed_loop_poles(c, A, an.B[:, 0], an.out_position, model.delay)).max() < 1.0
