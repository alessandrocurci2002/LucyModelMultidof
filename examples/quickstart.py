#!/usr/bin/env python3
"""Cascade and LQG-PI on the test circle, nominal and with the generic disturbance torque of the thesis.

    python examples/quickstart.py
"""

import numpy as np

from lucy_sim import Disturbance, GimbalModel
from lucy_sim.pointing.geometry import Circle, line_of_sight, pointing_error, statistics

model = GimbalModel.load()                                    # config/ and data/ of the repository
t = model.time(Circle.duration())                             # 5 s start + one orbit, 0.1 ms step
ref = Circle.reference(t)                                     # absolute azimuth and elevation of the drone [rad]
k = np.arange(0, len(t), 10)                                  # metrics at 1 ms

for name, d in (("nominal", Disturbance()), ("torque 0.3 N m rms, 0.5-5 Hz", Disturbance(gauss=0.3, band=(0.5, 5.0)))):
    for controller in ("cascade", "lqg"):
        out = model.simulate(ref, t, controller, d, seed=10)
        n = line_of_sight(ref["azimuth"][0] + out["azimuth"]["los"][k], ref["polar"][0] + out["polar"]["los"][k])
        E, miss = pointing_error(n, Circle.truth(t[k]))
        s = statistics(t[k], E, miss)
        print(f"{name:32s} {controller:8s} E rms {s['E_rms_mrad']:.2f} mrad, on target {s['on_target_pct']:.1f} %, "
              f"max current {np.abs(out['azimuth']['i']).max():.1f} / {np.abs(out['polar']['i']).max():.1f} A")
