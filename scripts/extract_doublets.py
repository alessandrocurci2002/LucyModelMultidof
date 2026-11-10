#!/usr/bin/env python3
"""Extracts the typical doublet commands used by the LQG tuning from the raw PLC exports of the V2 doublets.

    python scripts/extract_doublets.py <folder with MPCLucy_doublets_*.csv>

Writes ``data/doublets/azimuth.npz`` and ``data/doublets/polar.npz`` (one window per amplitude: time, filtered current,
encoder and command). The raw files (about 18 MB each) are not part of the repository.
"""

from __future__ import annotations

import sys
from pathlib import Path

from lucy_sim.paths import ROOT
from lucy_sim.tuning.doublets import save_typical_commands, typical_commands_from_csv


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = ROOT / "data" / "doublets"
    out.mkdir(parents=True, exist_ok=True)
    for axis in ("azimuth", "polar"):
        ws = typical_commands_from_csv(Path(sys.argv[1]), axis)
        save_typical_commands(ws, out / f"{axis}.npz")
        print(axis, [w.name for w in ws])


if __name__ == "__main__":
    main()
