"""Closed-loop model of the LUCY mount: mechanical plant, motors and drives.

The library has three independent components and one interconnection:

* :mod:`lucy_cl.plant`  -- nonlinear mechanical model (by default the 8-DOF H6 lumped model of
  ``gimbal_mdof``, any configuration of that package can be used), with the bearing friction;
* :mod:`lucy_cl.motor`  -- Kollmorgen TBM(S)-7646-A frameless motor (datasheet parameters);
* :mod:`lucy_cl.driver` -- Kollmorgen AKD2G drive: position, velocity and current loops with the
  biquad filters, parameters from the drive files and the WorkBench exports;
* :mod:`lucy_cl.closed_loop` -- connects the three: time-domain simulation of the nonlinear closed
  loop and frequency responses of the linearised loops.

``main.py`` (one level up) builds the three components from ``config.yaml`` and runs the analyses.
"""

from __future__ import annotations

import sys
from pathlib import Path

# the plant uses the gimbal_mdof package of the parent project (sim/gimbal_multidof/src)
_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

AXES = ("azimuth", "elevation")
"""Axis names used by the closed-loop library (``gimbal_mdof`` calls them ``azimut`` and ``polare``)."""
GIMBAL_NAME = {"azimuth": "azimut", "elevation": "polare"}
