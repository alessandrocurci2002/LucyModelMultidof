#!/usr/bin/env python3
"""Simulated campaign of the thesis: cascade and LQG-PI on the test circle and on the recorded drone trajectories,
scenarios 0-3 (nominal, disturbance torque, higher friction, both), three realisations of the torque.

    python scripts/run_hil.py [--out results/hil] [--workers 3] [--quick]

Writes ``<out>/summary.json`` (metrics of every case) and ``<out>/series/<case>.npz`` (time series at 1 ms, float32).
``--quick``: one case per scenario on the circle with the LQG-PI. The whole campaign is 80 simulations of 30-90 s
each; with 4 or more workers watch the memory (each worker holds the model, about 1 GB).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from lucy_sim.campaign import SCENARIOS, SEEDS, Case, run_campaign
from lucy_sim.paths import ROOT


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "hil")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    cases = [Case("circle", "lqg", s, SEEDS[0] if SCENARIOS[s].gauss else None) for s in SCENARIOS] if args.quick else None
    run_campaign(args.out, cases, args.workers)


if __name__ == "__main__":
    main()
