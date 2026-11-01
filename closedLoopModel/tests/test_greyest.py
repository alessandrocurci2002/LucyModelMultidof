"""Checks of the grey-box estimator (lucy_cl.greyest) on systems with known parameters.

Run from ``closedLoopModel``: ``../.venv/bin/python -m pytest tests/test_greyest.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from lucy_cl.greyest import (FrequencyData, GreyestOptions, GreyModel, GreyParameter, TimeData,  # noqa: E402
                             compare, greyest, state_space_response, zoh)

TRUE = {"m": 1.0, "k": 4.0e4, "c": 8.0}   # mass-spring-damper: 31.8 Hz, zeta 0.02


def msd(p, Ts):
    """Force -> velocity of a mass on a spring and a damper."""
    A = [[0.0, 1.0], [-p["k"] / p["m"], -p["c"] / p["m"]]]
    return A, [[0.0], [1.0 / p["m"]]], [[0.0, 1.0]], [[0.0]]


def msd_model(k0=1.5e4, c0=20.0, **kw) -> GreyModel:
    return GreyModel.create(msd, {"m": GreyParameter("m", 1.0, free=False),
                                  "k": dict(value=k0, minimum=1e2, maximum=1e7),
                                  "c": dict(value=c0, minimum=1e-3, maximum=1e3)}, **kw)


def msd_data(noise: float = 0.01, seed: int = 0, f=None) -> FrequencyData:
    rng = np.random.default_rng(seed)
    f = np.geomspace(5, 200, 300) if f is None else f
    G = state_space_response(*(np.asarray(X, float) for X in msd(TRUE, 0)), 2 * np.pi * f)[:, 0, 0]
    G = G * (1 + noise * (rng.standard_normal(f.size) + 1j * rng.standard_normal(f.size)))
    return FrequencyData(f, G, name="msd")


def test_frequency_complex_error_recovers_parameters():
    est = greyest(msd_data(), msd_model())
    p, rep = est.parameters, est.report
    assert p["m"] == 1.0                                   # fixed parameter untouched
    assert p["k"] == pytest.approx(TRUE["k"], rel=0.005)
    assert p["c"] == pytest.approx(TRUE["c"], rel=0.03)
    assert rep.fit_percent[0, 0] > 95
    assert rep.status == "converged"
    assert 0 < rep.std["k"] < 0.01 * p["k"]
    assert rep.covariance.shape == (2, 2)


def test_frequency_log_error_and_band():
    # starting far away (resonance at 3x the true frequency): the log error with a band still converges
    opt = GreyestOptions(frequency_error="log", weighting_filter=[10, 100])
    est = greyest(msd_data(), msd_model(k0=4e5, c0=1.0), opt)
    assert est.parameters["k"] == pytest.approx(TRUE["k"], rel=0.01)
    assert est.parameters["c"] == pytest.approx(TRUE["c"], rel=0.05)
    _, fit = compare(msd_data(), est.model, options=opt)
    assert fit[0, 0] > 95


def test_multi_experiment():
    d1 = msd_data(seed=1, f=np.geomspace(5, 40, 150))
    d2 = msd_data(seed=2, f=np.geomspace(25, 200, 150))
    est = greyest([d1, d2], msd_model())
    assert est.parameters["k"] == pytest.approx(TRUE["k"], rel=0.005)
    assert est.report.fit_percent.shape == (2, 1)


def test_bounds_respected_and_flagged():
    model = GreyModel.create(msd, {"m": GreyParameter("m", 1.0, free=False),
                                   "k": dict(value=1.5e4, minimum=1e2, maximum=2.0e4),
                                   "c": dict(value=20.0, minimum=1e-3, maximum=1e3)})
    est = greyest(msd_data(), model)
    k = next(p for p in est.report.parameters if p.name == "k")
    assert k.value <= 2.0e4 * (1 + 1e-9)
    assert k.at_bound


def test_frequency_only_model_hysteretic():
    # receptance of a mass on a hysteretic spring k (1 + j eta): not rational in s, frfun only
    def frf(p, w):
        return 1.0 / (p["k"] * (1 + 1j * p["eta"]) - p["m"] * w**2)

    f = np.geomspace(2, 50, 200)
    G = frf({"m": 2.0, "k": 1.2e4, "eta": 0.3}, 2 * np.pi * f)
    model = GreyModel.create(None, {"m": GreyParameter("m", 2.0, free=False),
                                    "k": dict(value=5e3, minimum=1.0, maximum=1e6),
                                    "eta": dict(value=1.0, minimum=1e-3, maximum=10.0)}, frfun=frf)
    est = greyest(FrequencyData(f, G), model, GreyestOptions(frequency_error="log"))
    assert est.parameters["k"] == pytest.approx(1.2e4, rel=1e-4)
    assert est.parameters["eta"] == pytest.approx(0.3, rel=1e-4)


# --------------------------------------------------------------------------- time domain


def dc_motor(p, Ts):
    """DC motor of the MATLAB greyest example: states angle and rate, outputs both."""
    tau, G = p["tau"], p["G"]
    return [[0.0, 1.0], [0.0, -1.0 / tau]], [[0.0], [G / tau]], np.eye(2), np.zeros((2, 1))


def motor_data(seed=0, x0=(0.0, 0.0), noise=0.01, Ts=0.01, N=1500) -> TimeData:
    rng = np.random.default_rng(seed)
    t = np.arange(N) * Ts
    u = np.repeat(rng.choice([-1.0, 1.0], N // 25 + 1), 25)[:N]   # PRBS-like
    true = GreyModel.create(dc_motor, {"tau": 0.5, "G": 2.0})
    y = true.simulate(t, u, np.array(x0))
    y += noise * y.std(axis=0) * rng.standard_normal(y.shape)
    return TimeData(t, u, y, name="dc motor")


def motor_model() -> GreyModel:
    return GreyModel.create(dc_motor, {"tau": dict(value=1.5, minimum=1e-3, maximum=100.0),
                                       "G": dict(value=0.8, minimum=1e-3, maximum=100.0)})


def test_time_domain_dc_motor():
    est = greyest(motor_data(), motor_model(), GreyestOptions(output_weight="normalize"))
    assert est.parameters["tau"] == pytest.approx(0.5, rel=0.01)
    assert est.parameters["G"] == pytest.approx(2.0, rel=0.01)
    assert (est.report.fit_percent > 95).all()


def test_time_domain_initial_state_estimate():
    data = motor_data(x0=(0.3, -1.0), noise=0.002)
    est = greyest(data, motor_model(), GreyestOptions(output_weight="normalize", initial_state="estimate"))
    assert est.parameters["tau"] == pytest.approx(0.5, rel=0.01)
    assert est.x0[0] == pytest.approx([0.3, -1.0], abs=0.02)


def test_discrete_model_time_data():
    Ts = 0.01

    def dc_discrete(p, Ts):
        A, B, C, D = (np.asarray(X, float) for X in dc_motor(p, 0.0))
        Ad, Bd = zoh(A, B, Ts)
        return Ad, Bd, C, D

    model = GreyModel.create(dc_discrete, {"tau": dict(value=1.5, minimum=1e-3, maximum=100.0),
                                           "G": dict(value=0.8, minimum=1e-3, maximum=100.0)}, Ts=Ts)
    est = greyest(motor_data(Ts=Ts), model, GreyestOptions(output_weight="normalize"))
    assert est.parameters["tau"] == pytest.approx(0.5, rel=0.01)


# --------------------------------------------------------------------------- statistics


@pytest.mark.parametrize("error, covariance", [("complex", "robust"), ("log", "classic"), ("log", "robust")])
def test_covariance_matches_monte_carlo(error, covariance):
    """The reported standard deviation agrees with the spread of the estimates over noise realisations. The noise
    is proportional to the response: with the complex error its variance changes along the data, and only the
    robust covariance holds (the classic one underestimates the spread about 7 times)."""
    opt = GreyestOptions(frequency_error=error, covariance=covariance)
    ks, stds = [], []
    for seed in range(30):
        est = greyest(msd_data(noise=0.03, seed=100 + seed), msd_model(k0=3.5e4, c0=10.0), opt)
        ks.append(est.parameters["k"])
        stds.append(est.report.std["k"])
    ratio = np.std(ks) / np.mean(stds)
    assert 0.6 < ratio < 1.6


def test_regularization_pulls_unidentified_parameter_to_nominal():
    # data far below the resonance: the damper is not identified, the regularization keeps it at its nominal
    data = msd_data(f=np.geomspace(0.5, 3, 50))
    opt = GreyestOptions(regularization_lambda=1.0, regularization_R=np.diag([0.0, 1.0]))
    est = greyest(data, msd_model(k0=3e4, c0=20.0), opt)
    assert est.parameters["k"] == pytest.approx(TRUE["k"], rel=0.02)
    assert est.parameters["c"] == pytest.approx(20.0, rel=0.05)
    assert "c" in est.report.summary()


def test_soft_constraint_limits_damping():
    """penalty_fun: with zeta <= 0.01 (true 0.02) the damper stops at the limit c = 2 zeta sqrt(k m)."""
    def penalty(p):
        zeta = p["c"] / (2 * np.sqrt(p["k"] * p["m"]))
        return 1e3 * np.array([max(zeta - 0.01, 0.0)])

    base = msd_model()
    model = GreyModel(base.odefun, base.parameters, penalty_fun=penalty)
    est = greyest(msd_data(), model, GreyestOptions(frequency_error="log"))
    zeta = est.parameters["c"] / (2 * np.sqrt(est.parameters["k"]))
    assert zeta == pytest.approx(0.01, abs=5e-4)
    assert est.parameters["k"] == pytest.approx(TRUE["k"], rel=0.01)


# --------------------------------------------------------------------------- adapter of the LUCY plant


def test_gimbal_adapter_reproduces_loop_frf():
    """At the nominal values the adapter returns the plant of loop_frf; the hub springs split into a geometric
    part and the ring reproduce the lumped values."""
    import yaml

    import main as mn
    from lucy_cl.closed_loop import loop_frf
    from lucy_cl.gimbal_grey import Bounds, PlantGreyModel, axis_parameters
    from lucy_cl.plant_fit import K1, RingMap

    cfg = yaml.safe_load((ROOT / "config_vtesi_greybox.yaml").read_text())
    cl, tests, _ = mn.build(cfg, ROOT)
    rm = RingMap.from_config((ROOT / cfg["plant"]["config"]).resolve())
    f = np.geomspace(1, 400, 50)
    k1 = cl.plant.value(K1)
    for ax in ("azimuth", "elevation"):
        ref = loop_frf(cl, ax, f, drive=tests[ax]).plant
        m = PlantGreyModel(cl, ax, tests[ax], rm, axis_parameters(cl, ax, rm))
        G = m.grey_model(Bounds()).frequency_response(2 * np.pi * f)[:, 0, 0]
        assert np.max(np.abs(G / ref - 1)) < 1e-12
    assert cl.plant.value(K1) == pytest.approx(k1, rel=1e-12)
