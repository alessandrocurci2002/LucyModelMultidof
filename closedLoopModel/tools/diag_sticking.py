"""Diagnosis of the slow integration in the sticking phase of el_step_500mdeg (open issue, see STATO_LAVORI.md).

From closedLoopModel:
    ../.venv/bin/python tools/diag_sticking.py save 3.0        # integrates to t = 3.0 s and saves the state
    ../.venv/bin/python tools/diag_sticking.py Radau:1e-6:1 3.0 4.0   # method:rtol:atol_scale, from t0 to t1
"""

import signal as sg
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import main as mn  # noqa: E402

STATE = ROOT / "tools" / "state_el_step_500mdeg.npy"
cfg = yaml.safe_load((ROOT / "config_zeta10.yaml").read_text())
cl, tests, wb = mn.build(cfg, ROOT)
sc = [s for s in mn.scenarios(cfg) if s.name == "el_step_500mdeg"][0]
cmd = {ax: sc.command(ax) for ax in ("azimuth", "elevation")}
cnt = [0]


def f(t, x):
    cnt[0] += 1
    return cl.evaluate(t, x, cmd)[0]


mode = sys.argv[1]
if mode == "save":
    t_save = float(sys.argv[2])
    sol = solve_ivp(f, (0, t_save), cl.initial_state(), method="Radau", rtol=1e-6, atol=cl.atol())
    np.save(STATE, sol.y[:, -1])
    print("saved state at", t_save, "s,", cnt[0], "evaluations")
    sys.exit()

t0, t1 = float(sys.argv[2]), float(sys.argv[3])
x0 = np.load(STATE)
method, rtol, ascale = mode.split(":")


def handler(*_):
    raise TimeoutError


sg.signal(sg.SIGALRM, handler)
sg.alarm(120)
tc = time.process_time()
n, jb = cl.plant.n, cl.plant.index("elevation", "bearing")
try:
    sol = solve_ivp(f, (t0, t1), x0, method=method, rtol=float(rtol), atol=cl.atol() * float(ascale))
    v = sol.y[n + jb]
    print(mode, "success", sol.success, "steps", sol.t.size, "nfev", cnt[0], "cpu", round(time.process_time() - tc, 1),
          "v_bearing", v.min(), v.max(), "min step", np.diff(sol.t).min())
except TimeoutError:
    print(mode, "TIMEOUT after 120 s, nfev", cnt[0])
