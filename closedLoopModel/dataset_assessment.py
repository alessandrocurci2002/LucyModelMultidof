r"""Information content of a set of sequences for the current-driven model (dataset assessment).

At the parameters of a fitted run, computes the Jacobian of the normalised position residuals on the moves of a
set of sequences (as ``identifiability.py``) and reports condition number, effective independent samples,
standard errors and correlations. Used to compare the 1 deg sequences of VClosedLoop with the doublets of
VClosedLoop/V2 as identification sets, before identifying anything.

Usage (from this folder)::

    ../.venv/bin/python -I dataset_assessment.py [azimuth] [elevation]

With axes given, only those are recomputed; the other entries of information.json are kept.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import current_input_study as cs  # noqa: E402
import main as mn  # noqa: E402
from identifiability import REL, STEP, SIGNED, effective_samples  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl import plant_fit as pf  # noqa: E402
from lucy_cl.validation import read_sequence  # noqa: E402

LI = HERE.parent.parent / "LucyIdentification" / "VClosedLoop"
SETS = {
    "VClosedLoop 1 deg (old identification set)": [LI / "MPCLucy_step_1deg.csv", LI / "MPCLucy_1degAz_zeroAz_1degPol__zeroPol.csv"],
    "VClosedLoop/V2 doublets": sorted((LI / "V2").glob("MPCLucy_doublets_*.csv")),
}
MODEL = ("I9a (cable spring)", "config_current_input_I9a.yaml", "results_current_input_I9a",
         {"azimuth": ["scale", "sigma0", "Tc", "Ts_ratio", "K_cable", "theta0_cable_deg"],
          "elevation": ["scale", "sigma0", "Tc", "Ts_ratio", "vs", "b_{\\theta,h}", "K_cable", "theta0_cable_deg"]})


def main() -> None:
    res_dir = HERE / "results_dataset_v2"
    res_dir.mkdir(exist_ok=True)
    label, cfg_name, res_name, axes = MODEL
    cc = yaml.safe_load((HERE / cfg_name).read_text())
    cfg = yaml.safe_load((HERE / cc["base_config"]).read_text())
    final = json.loads((HERE / res_name / "summary.json").read_text())["final"]
    cl, _, _ = mn.build(cfg, HERE)
    rm = pf.RingMap.from_config(cl.plant.config)
    only = set(sys.argv[1:]) or set(axes)  # axes to (re)compute; the others are kept from information.json
    js_path = res_dir / "information.json"
    js = json.loads(js_path.read_text()) if js_path.exists() and only != set(axes) else {}
    s2_ref: dict[str, float] = {}  # noise level of the first set, to compare the sets at the same noise
    for set_name, files in SETS.items():
        seq_dir = files[0].parent
        seqs = {f.stem: read_sequence(f) for f in files}
        wins = [w for n, s in seqs.items() for w in cs.windows(s, n, 0.3, 0.5)]
        with mp.get_context("fork").Pool(8, initializer=cs._init,
                                         initargs=(cfg, str(seq_dir), 2.0, float(cc["input_lowpass_hz"]),
                                                   {"a": rm.a, "b": rm.b}, final)) as pool:
            for ax, names in axes.items():
                if ax not in only:
                    continue
                w_ax = [w for w in wins if w.axis == ax]
                base = {a_: dict(final[a_]) for a_ in AXES}

                def resid(vals: dict) -> list[np.ndarray]:
                    return [cs.residual_of(r) for r in pool.map(cs.job, [(w, vals, 120.0) for w in w_ax])]

                r0s = resid(base)
                r0 = np.concatenate(r0s)
                J = np.empty((r0.size, len(names)))
                for k, n in enumerate(names):
                    cols = []
                    for sg in (+1, -1):
                        v = {a_: dict(base[a_]) for a_ in AXES}
                        v[ax][n] = base[ax][n] + sg * STEP[n] if n in SIGNED else base[ax][n] * np.exp(sg * REL)
                        cols.append(np.concatenate(resid(v)))
                    J[:, k] = (cols[0] - cols[1]) / (2 * (STEP[n] if n in SIGNED else REL))
                n_eff = sum(effective_samples(r) for r in r0s)
                s2 = float(r0 @ r0) / max(1.0, n_eff - len(names))
                F = J.T @ J * (n_eff / r0.size)
                Js = J / np.linalg.norm(J, axis=0)
                sv = np.linalg.svd(Js, compute_uv=False)
                s2_ref.setdefault(ax, s2)
                cov = s2 * np.linalg.pinv(F)
                sd = np.sqrt(np.clip(np.diag(cov), 0, None))
                corr = cov / np.outer(sd, sd)
                # standard errors at the noise level of the first set: the information of the moves alone,
                # independent of how well the fitted parameters extrapolate to this set
                sd_ref = np.sqrt(np.clip(np.diag(s2_ref[ax] * np.linalg.pinv(F)), 0, None))
                key = f"{set_name} / {ax}"
                js[key] = {"set": set_name, "axis": ax, "names": names, "moves": len(w_ax), "n_effective": n_eff,
                           "condition": float(sv[0] / sv[-1]), "s2": s2, "sd": sd.tolist(),
                           "sd_common_noise": sd_ref.tolist(), "corr": corr.tolist(),
                           "fit_median": float(np.median([100 * (1 - np.linalg.norm(r)) for r in r0s]))}
                print(f"{key}: {len(w_ax)} moves, n_eff {n_eff:.0f}, condition {sv[0] / sv[-1]:.3g}, "
                      + ", ".join(f"{n} {('%.2g' % sd_ref[k]) if n in SIGNED else ('%.0f%%' % (100 * sd_ref[k]))}" for k, n in enumerate(names)),
                      flush=True)
                js_path.write_text(json.dumps(js, indent=2))  # after every case: a stopped run keeps what it computed
                (res_dir / "information.md").write_text(report(js, label))


def report(js: dict, label: str) -> str:
    md = ["# Information content of the identification sets", "",
          f"Model and parameters of {label}; Jacobian of the normalised position residuals on the moves of "
          "each set (windows from 0.3 s before a move to 0.5 s after the end of its command).", ""]
    for key, e in js.items():
        names, sd, sd_ref, corr = e["names"], e["sd"], e["sd_common_noise"], np.array(e["corr"])
        rows = []
        for k, n in enumerate(names):
            se, se_r = ((f"{x[k]:.2g} deg" if n in SIGNED else f"{100 * x[k]:.0f} %") for x in (sd, sd_ref))
            big = [f"{names[j]} ({corr[k, j]:+.2f})" for j in range(len(names)) if j != k and abs(corr[k, j]) > 0.9]
            rows.append([n, se, se_r, ", ".join(big) or "-"])
        md += [f"## {key.replace(' / ', ': ')}", "",
               f"{e['moves']} moves, effective independent samples {e['n_effective']:.0f}, condition number "
               f"{e['condition']:.3g}, median fit at the {label} parameters {e['fit_median']:.1f} %.", "",
               out.markdown_table(["parameter", "standard error", "standard error at the noise of the first set",
                                   "correlations |r| > 0.9"], rows), ""]
    return "\n".join(md)

if __name__ == "__main__":
    main()
