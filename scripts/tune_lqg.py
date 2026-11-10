#!/usr/bin/env python3
"""Re-tunes the LQG-PI of LUCY with the procedure of the thesis (Gawronski 2008, par. 8.2).

    python scripts/tune_lqg.py [--out results/tuning] [--axes azimuth polar]

Writes ``<out>/summary.json`` (every evaluated point, as in the thesis) and ``<out>/controllers_weights.yaml``, the
``lqg.axes`` block to copy into ``config/controllers.yaml``. About 30 minutes per axis; the axes run in parallel.
"""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

import yaml  # noqa: E402

from lucy_sim import GimbalModel  # noqa: E402
from lucy_sim.paths import ROOT  # noqa: E402
from lucy_sim.tuning.gawronski import tune, weights  # noqa: E402


def _one(axis: str) -> dict:
    return tune(GimbalModel.load(ROOT), (axis,))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "tuning")
    ap.add_argument("--axes", nargs="+", default=["azimuth", "polar"], choices=["azimuth", "polar"])
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(len(args.axes)) as ex:
        summary = {k: v for part in ex.map(_one, args.axes) for k, v in part.items()}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1))
    (args.out / "controllers_weights.yaml").write_text(yaml.safe_dump({"axes": weights(summary)}, sort_keys=False))
    for ax, w in weights(summary).items():
        print(ax, w)


if __name__ == "__main__":
    main()
