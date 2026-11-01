"""Checks of the closed-loop library: components against their references and interconnection consistency.

Run from ``closedLoopModel``: ``../.venv/bin/python -m pytest tests``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from lucy_cl import AXES, GIMBAL_NAME  # noqa: E402
from lucy_cl.closed_loop import Command, Scenario, loop_frf, margins, state_space_frf  # noqa: E402
from lucy_cl.driver import Biquad, controller_workbench  # noqa: E402
from lucy_cl.motor import MotorParams  # noqa: E402
from lucy_cl.plant import GimbalPlant, dahl_describing, dahl_from_hysteretic  # noqa: E402

import main as mn  # noqa: E402


@pytest.fixture(scope="module")
def loop():
    cfg = yaml.safe_load((ROOT / "config_zeta10.yaml").read_text())
    return mn.build(cfg, ROOT)


# --------------------------------------------------------------------------- components


def test_motor_datasheet_consistency():
    for temp in ("cold", "hot"):
        p = MotorParams.from_yaml(ROOT / "data" / "motor_TBM7646A.yaml", temp)
        assert abs(p.consistency()) < 0.01  # Kt = sqrt(3) Kb within 1 %
        assert p.Ke_ph == pytest.approx(p.Kt / 3, rel=0.01)  # power balance 3 E i = Kt i omega


def test_biquad_state_space_equals_transfer_function():
    for b in (Biquad(65, 5.62, 65, 17.8), Biquad(500, 3, 14, 3), Biquad(5000, 0.707, 800, 0.707)):
        A, B, C, D = b.state_space()
        for fk in (1.0, 14.0, 65.0, 800.0, 3000.0):
            s = 2j * np.pi * fk
            H = (C @ np.linalg.solve(s * np.eye(2) - A, B) + D)[0, 0]
            assert H == pytest.approx(b.response(np.array([s]))[0], rel=1e-9)


def test_drive_model_reproduces_workbench_traces(loop):
    _, tests, wb = loop
    for ax in AXES:
        w, d = wb[ax], tests[ax]
        s = 2j * np.pi * w.f
        e = controller_workbench(d, s) / w.response("Controller")
        assert np.abs(20 * np.log10(np.abs(e))).max() < 0.01
        assert np.abs(np.angle(e, deg=True)).max() < 0.5
        e = d.F_ar_fb(s) / w.response("Feedback Biquad")
        assert np.abs(20 * np.log10(np.abs(e))).max() < 0.1


def test_workbench_open_loop_definition(loop):
    _, _, wb = loop
    for w in wb.values():
        ol = w.response("Controller") * w.response("Plant") * w.response("Feedback Biquad")
        assert np.allclose(ol, w.response("Open Loop"), rtol=1e-6)
        cl = w.response("Controller") * w.response("Plant") / (1 + ol)
        assert np.allclose(cl, w.response("Closed Loop"), rtol=1e-6)


def test_dahl_equivalent_stiffness():
    from scipy.integrate import solve_ivp

    k, eta, A = 1418.0, 0.32, 7e-4
    s0, Tc = dahl_from_hysteretic(k, eta, A)
    assert dahl_describing(s0, Tc, A) == pytest.approx(k * (1 + 1j * eta), rel=1e-6)
    w = 2 * np.pi * 2.0
    sol = solve_ivp(lambda t, T: s0 * (A * w * np.cos(w * t) - abs(A * w * np.cos(w * t)) * T / Tc), (0, 3.0), [0.0],
                    max_step=1e-3, rtol=1e-9, atol=1e-12, dense_output=True)
    t = np.linspace(2.0, 2.5, 20001)[:-1]
    T = sol.sol(t)[0]
    a, b = 2 * np.mean(T * np.sin(w * t)), 2 * np.mean(T * np.cos(w * t))
    assert a / A == pytest.approx(k, rel=1e-3)
    assert b / a == pytest.approx(eta, rel=1e-3)


def test_plant_equals_gimbal_mdof_model():
    """Without modifications the plant is the lumped model of gimbal_mdof, ground constraint included."""
    from gimbal_mdof import linearize as lz
    from gimbal_mdof import sperimentale as spm
    from gimbal_mdof import vincolo_terra as vt

    cfg = ROOT.parent / "config_hub_telescopeonring_alhub_fit_ground.yaml"
    pl = GimbalPlant(cfg, balanced_telescope=False)
    ml = vt.applica(lz.linearizza(pl.system.accoppiato, 0.0, 0.0), pl.par, True)
    f = np.array([1.0, 3.0, 18.0, 58.0, 150.0, 400.0])
    for ax in AXES:
        ref = spm.frf_velocita(ml, GIMBAL_NAME[ax], f, 1.0)[spm.canali(GIMBAL_NAME[ax])[0]] / (2j * np.pi * f)
        assert np.allclose(pl.frf(ax, f)["encoder"], ref, rtol=1e-5)


# --------------------------------------------------------------------------- interconnection


def test_time_domain_linearisation_equals_frequency_model(loop):
    cl, _, _ = loop
    A, B, C = cl.linearise()
    f = np.logspace(0, 3, 60)
    for j, ax in enumerate(AXES):
        Ht = state_space_frf(A, B, C, f, j, j)
        Hf = loop_frf(cl, ax, f, ground="dahl").position
        assert np.max(np.abs(Ht - Hf) / np.abs(Hf)) < 1e-5


def test_axes_nearly_decoupled_at_rest(loop):
    """At rest the axes couple only through the CAD product of inertia I_yz (~1e-6 kg m^2)."""
    cl, _, _ = loop
    A, B, C = cl.linearise()
    f = np.array([1.0, 10.0, 100.0])
    for j_in in (0, 1):
        direct = np.abs(state_space_frf(A, B, C, f, j_in, j_in))
        cross = np.abs(state_space_frf(A, B, C, f, j_in, 1 - j_in))
        assert np.all(cross < 1e-3 * direct)


def test_small_step_nonlinear_close_to_linear(loop):
    from scipy import signal

    cl, _, _ = loop
    val = np.deg2rad(1e-3)
    sc = Scenario("t", {"azimuth": Command("step", val)}, 0.3)
    r = cl.simulate(sc, n_out=301)
    assert r.success
    A, B, C = cl.linearise()
    u = np.zeros((r.t.size, 2))
    u[:, 0] = val
    _, y, _ = signal.lsim((A, B, C, np.zeros((2, 2))), u, r.t, interp=False)
    err = np.abs(r.signals["azimuth"]["p_enc"] - y[:, 0]).max() / val
    assert err < 0.05  # the Dahl element is close to its linearisation for small motions


def test_current_limit_respected(loop):
    cl, _, _ = loop
    sc = Scenario("t", {"azimuth": Command("step", np.deg2rad(0.5))}, 0.15)
    r = cl.simulate(sc, n_out=151)
    s = r.signals["azimuth"]
    lim = cl.drives["azimuth"].p.i_limit
    assert np.abs(s["i_cmd"]).max() <= lim + 1e-9
    assert np.abs(s["i_cmd"]).max() == pytest.approx(lim)


def test_h6_elevation_stability_depends_on_damping(loop):
    """Finding of 2026-10-19: with zeta_1 = 0.01 (placeholder Q8 of H6) the elevation loop is unstable at the
    hub-telescope mode (~300 Hz); with zeta_1 = 0.10 it is stable."""
    cl, _, _ = loop
    pl = cl.plant
    z0 = pl.zeta
    f = np.logspace(0, np.log10(2000), 3000)
    try:
        pl.set_zeta(0.01)
        assert np.linalg.eigvals(cl.linearise()[0]).real.max() > 1.0
        assert margins(f, loop_frf(cl, "elevation", f).loop)["gm_db"] < 0
        pl.set_zeta(0.10)
        assert np.linalg.eigvals(cl.linearise()[0]).real.max() < 1e-6
        assert margins(f, loop_frf(cl, "elevation", f).loop)["gm_db"] > 3
    finally:
        pl.set_zeta(z0)
