"""Confronto fra l'LQG-PI della tesi (errore stimato) e quello della Fig. 8.4b (errore misurato): stesse simulazioni, stessi scenari.

    .venv/bin/python -I lqg_pi_errore_misurato/confronto.py

Scrive results/confronto.md.
"""

from __future__ import annotations

import json
from pathlib import Path

QUI = Path(__file__).resolve().parent
TESI = QUI.parent / "results_lqg_tuning" / "no_bearings" / "thesis_hil" / "numbers.json"
NUOVO = QUI / "results" / "hil" / "numbers.json"
TRAIETTORIE = {"circle": "circle", "dataset7_config_v3": "dataset7, config_v3", "dataset8_config_v3": "dataset8, config_v3",
               "dataset7_EKF2": "dataset7, EKF2", "dataset8_EKF2": "dataset8, EKF2"}
SCENARI = {0: "nominal", 1: "torque", 2: "friction", 3: "combined"}


def main() -> None:
    a, b = json.loads(TESI.read_text()), json.loads(NUOVO.read_text())
    righe = ["| Trajectory | Scenario | eta cascade | eta LQG est. | eta LQG meas. | E rms LQG est. / meas. [mrad] | effect az est. / meas. [urad] | effect pol est. / meas. [urad] | i max pol est. / meas. [A] |",
             "|---|---|---|---|---|---|---|---|---|"]
    for t, nome in TRAIETTORIE.items():
        for s, sn in SCENARI.items():
            k = lambda c, q: f"{t}|{c}|{s}|{q}"  # noqa: E731
            if k("lqg", "E_rms_mrad") not in b:
                continue
            eta = "on_target_pct" if t == "circle" else "eta"
            nan = (float("nan"), float("nan"))
            ea, eb = a.get(k("lqg", "effect"), nan), b.get(k("lqg", "effect"), nan)
            righe.append(f"| {nome} | {sn} | {b[k('cascade', eta)]:.1f} | {a[k('lqg', eta)]:.1f} | {b[k('lqg', eta)]:.1f} | "
                         f"{a[k('lqg', 'E_rms_mrad')]:.3f} / {b[k('lqg', 'E_rms_mrad')]:.3f} | {ea[0]:.0f} / {eb[0]:.0f} | {ea[1]:.0f} / {eb[1]:.0f} | "
                         f"{a[k('lqg', 'i_max_el_A')]:.1f} / {b[k('lqg', 'i_max_el_A')]:.1f} |")
    (QUI / "results" / "confronto.md").write_text("\n".join(righe) + "\n")
    print("\n".join(righe))


if __name__ == "__main__":
    main()
