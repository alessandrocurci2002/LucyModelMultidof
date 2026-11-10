"""One point of the LQG tuning: the final weights of the thesis give the final point of the thesis tuning."""

import pytest

from lucy_sim.paths import ROOT
from lucy_sim.tuning.doublets import load_typical_commands
from lucy_sim.tuning.gawronski import Tuner, load_settings

FINAL = {"azimuth": {"k_p": 25.213396150325185, "bandwidth_Hz": 3.0736451954962756, "modulus": 0.7814428144195608, "i_max_A": 17.39081233891924},
         "polar": {"k_p": 644.416035809674, "bandwidth_Hz": 19.334810474075084, "modulus": 0.6806092177675414, "i_max_A": 9.075482744975883}}


@pytest.mark.slow
@pytest.mark.parametrize("axis", ("azimuth", "polar"))
def test_final_point(model, axis):
    T = Tuner(model, axis, load_settings(), load_typical_commands(ROOT / "data" / "doublets" / f"{axis}.npz"))
    w = model.lqg_settings["axes"][axis]
    x = T.point(w["q_p"], w["q_i"], w.get("q_f") or {}, w["rho"])
    assert not T.violations(x, True)
    r = FINAL[axis]
    for k in ("k_p", "bandwidth_Hz", "modulus"):
        assert x[k] == pytest.approx(r[k], rel=1e-8)
    assert x["doublets"]["i_max_A"] == pytest.approx(r["i_max_A"], rel=1e-6)
