# LUCY gimbal simulator

Closed-loop simulator of the LUCY alt-azimuth pointer: the version used in the thesis of A. Curci (2026), chapters on
the LQG controller and on the hardware-in-the-loop tests.

- **Plant.** Two axes, azimuth and polar (elevation). Each is a modal structure identified on the measured frequency
  responses, with Dahl friction on the bearings. The axes are coupled through the telescope by the Coriolis and
  gyroscopic torques of the Lagrange equations, and the azimuth inertia depends on the elevation.
- **Drive.** Kollmorgen AKD2G model with position, velocity and current loops, the biquad filters and the TBM(S)-7646-A
  motor.
- **Controllers.**
  - The *cascade* is the current configuration of LUCY: the position loop runs in the drive, with the effective gains
    identified on the doublets.
  - The *LQG-PI* is Gawronski's controller with the PI on the measured error (Fig. 8.4b). It runs in the master at
    1 kHz with one EtherCAT period of delay and replaces the position loop.
- **Tuning.** The thesis procedure for the LQG weights (Gawronski 2008, par. 8.2) on the measured doublets of LUCY.
- **Campaign.** The simulated campaign of the thesis: the test circle and the recorded drone trajectories (dataset7 and
  dataset8, estimators EKF2 and OpenVINS `config_v3`), with disturbance-torque and friction scenarios.

The plant with the Dahl bearings is the reference ("proxy of reality") on which the controllers are compared. The
measurements were used to identify it, not to compare controllers.

The regression tests check that the code reproduces the thesis results exactly: plant matrices, LQG gains, tuning
points and the metrics of the HIL campaign.

## Install

Python 3.10 or newer.

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest              # fast tests, about 10 s
.venv/bin/python -m pytest -m slow      # full thesis cases and tuning points, a few minutes
```

## Quick start

```python
from lucy_sim import Disturbance, GimbalModel
from lucy_sim.pointing.geometry import Circle

model = GimbalModel.load()                     # builds everything from config/ and data/ (a few seconds)
t = model.time(Circle.duration())              # 0.1 ms grid
ref = Circle.reference(t)                      # {"azimuth": psi(t), "polar": theta(t)}, absolute, rad
out = model.simulate(ref, t, controller="lqg", # or "cascade"
                     disturbance=Disturbance(gauss=0.3, band=(0.5, 5.0)), seed=10)
out["azimuth"]["los"], out["polar"]["i"]       # line of sight relative to the start [rad], motor current [A rms]
```

The same example with the pointing metrics is in `examples/quickstart.py`. Other entry points:

| | |
|---|---|
| `model.simulate_axis("polar", p_cmd, "lqg")` | one axis alone on a position command (relative, at the simulation step) |
| `model.lqg_controller(axis, q_p, q_i, q_f, rho)` | an LQG-PI with other weights, to pass as `controller={"azimuth": c_az, "polar": c_el}` |
| `python scripts/run_hil.py [--quick]` | the thesis campaign (80 simulations) → `results/hil/summary.json` and the series |
| `python scripts/tune_lqg.py` | re-tunes the LQG-PI (about 30 min per axis) → `results/tuning/` |
| `python scripts/extract_doublets.py <dir>` | rebuilds `data/doublets/` from the raw PLC exports of the V2 doublets |

`results/` is ignored by git.

## How the model is built

1. **Structure** (`config/gimbal.yaml`). The shafts, the hub and the CAD inertias give a torsional FEM. It is reduced
   to lumped chains on the master nodes and coupled through the telescope with the Lagrange equations (SymPy). The
   result is linearised at theta = 0 into one chain per axis.
2. **Plant** (`config/plant.yaml`, `data/drive`, `data/motor`). The identified factors scale those chains into the
   modal fork model of each axis. Adding the Dahl bearings, the motor and the AKD2G drive gives the cascade and the
   velocity loop of each axis, plus the reduced design model of the LQG.
3. **Control** (`config/controllers.yaml`). LQG weights and master timing give the LQG-PI of each axis.
4. **Simulation.** The two axes run together with the coupling torques and the disturbances, following a reference: the
   circle or a drone trajectory.

| Package | Content |
|---|---|
| `lucy_sim.structure` | `config` (reads `gimbal.yaml`), `shafts` and `hub` (inertia and stiffness of the parts), `fem`, `payload` (telescope tensor), `lumped` (chains on the master nodes), `lagrange` (coupled equations of motion, linearisation) |
| `lucy_sim.plant` | `fork` (identified modal model), `friction` (Dahl), `motor`, `drive` (AKD2G) and `workbench` (its files), `loop` (axis in closed loop), `axis` (assembled axis + LQG design model) |
| `lucy_sim.control.lqg` | modal form and truncation, synthesis (CARE/FARE), LQG-PI with the PI on the measured (thesis) or estimated error, closed-loop poles and transfer functions |
| `lucy_sim.simulation` | `loops` (single axis, two coupled axes), `coupling` (Lagrange torques through the telescope), `disturbances` |
| `lucy_sim.pointing` | `geometry` (pointing error, miss distance, effective pointing, test circle), `trajectories` (recorded drone flights) |
| `lucy_sim.tuning` | `doublets` (measured V2 doublets, metrics), `gawronski` (tuning procedure) |
| `lucy_sim.campaign` | scenarios and cases of the thesis campaign |
| `lucy_sim.model` | `GimbalModel`: loads and assembles everything |

## Changing the model

- **Mechanics** (`config/gimbal.yaml`). Shaft geometry, materials, CAD inertias, hub, telescope tensor. The file
  explains every field. The topology is fixed: on both axes the telescope hangs from the shaft chain through a spring
  (hub model on the polar axis, attachment to the load on the azimuth).
- **Identified plant** (`config/plant.yaml`):
  - `fork` holds the factors on the CAD chains;
  - `bearings` holds the Dahl friction;
  - `drive_gains` holds the effective gains of the cascade.

  The drive filters and the rest of the drive parameters come from the `.akd` file and from the WorkBench export in
  `data/drive/`.
- **Controllers** (`config/controllers.yaml`). LQG weights per axis, the master period and the delay. The PI can act
  on the measured error (thesis) or on the estimated one. `q_f: {0: w}` sets the weight of the first lightly damped
  mode of the design model.
- **Scenarios** (`lucy_sim.simulation.disturbances.Disturbance`). The available terms are:
  - a Gaussian torque with its band;
  - a static mean torque;
  - a step;
  - a sine;
  - a factor on the Coulomb friction (model error);
  - noise on the position measurement.

  The thesis scenarios are in `lucy_sim.campaign.SCENARIOS`.
- **A new controller in the master.** It must provide `initial_state()`, `output(state, r)` and
  `update(state, y, u_applied)`, plus `measure(state, y)` if it needs the measured position before computing the
  command. Pass it through `controller={...}`; `None` means the cascade.
- **New trajectories.** Add `data/trajectories/<name>/` with `truth.txt` and one file per estimator (columns
  `t x y z qx qy qz qw`), or pass any reference to `model.simulate`.

## Conventions

- **Units and angles.** SI units and radians inside the code. Currents are A rms, as in the drive.
- **Axis keys.** `azimuth` and `polar`. Polar = elevation; theta > 0 raises the pointing.
- **References.** Absolute angles of the line of sight. The loops follow them relative to their first value, and the
  outputs `los` and `p` are relative to the start.
- **Simulation step.** 0.1 ms: exact ZOH of the linear part, exact update of the Dahl state. The master runs at 1 ms.
- **No saturations.** The simulation has no current limit. Check the current in the outputs: the tuning keeps it below
  18 A on the doublets.
- **Coupling.** The coupling torques use the elevation and the rates at the start of each step. The modal damping of
  the azimuth stays the one at theta = 0.

## Data in this repository

| | |
|---|---|
| `data/drive/` | AKD2G parameter files and WorkBench exports of the V2 session (gains and filters active during the tests) |
| `data/motor/TBM7646A.yaml` | datasheet of the Kollmorgen TBM(S)-7646-A |
| `data/doublets/` | typical V2 doublets (one window per amplitude) used by the tuning |
| `data/trajectories/` | dataset7 and dataset8: TagSLAM ground truth (levelled), EKF2 and OpenVINS `config_v3`, first run |
| `tests/reference/` | values computed with the original thesis code, for the regression tests |

The identification itself (fit of the fork model on the frequency responses, fit of the effective gains on the
doublets) is not part of this repository: its results are the numbers in `config/plant.yaml`.
