r"""Identifiability of the parameters of the current-driven model at a fitted point.

Jacobian of the normalised position residuals on the identification moves (central differences in the
coordinates of the fit: logarithm of positive parameters, value of signed ones), then

* singular values of the Jacobian scaled per column, condition number;
* correlation matrix of the estimates, :math:`\mathrm{corr} = D^{-1/2}(J^TJ)^{-1}D^{-1/2}`;
* standard error of each parameter, :math:`\sigma^2\,\mathrm{diag}(J^TJ)^{-1}` with :math:`\sigma^2` from the
  residuals and an effective number of independent samples per move (consecutive samples of a smooth position
  history are strongly correlated: the effective count is estimated from the autocorrelation of the residuals).

Usage (from this folder)::

    ../.venv/bin/python -I identifiability.py
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

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import current_input_study as cs  # noqa: E402
import main as mn  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl import plant_fit as pf  # noqa: E402
from lucy_cl.validation import read_sequence  # noqa: E402

CASES = [  # (label, config, results, parameters per axis)
    ("first run (no cables)", "config_current_input.yaml", "results_current_input",
     {"azimuth": ["scale", "sigma0", "Tc", "Ts_ratio"]}),
    ("I7 (no cables)", "config_current_input_I7.yaml", "results_current_input_I7",
     {"elevation": ["scale", "sigma0", "Tc", "Ts_ratio", "vs", "d_x_mm", "b_{\\theta,h}"]}),
    ("I9a (cable spring, one theta_0)", "config_current_input_I9a.yaml", "results_current_input_I9a",
     {"azimuth": ["scale", "sigma0", "Tc", "Ts_ratio", "K_cable", "theta0_cable_deg"],
      "elevation": ["scale", "sigma0", "Tc", "Ts_ratio", "vs", "b_{\\theta,h}", "K_cable", "theta0_cable_deg"]}),
]
SIGNED = {"d_x_mm", "theta0_cable_deg"}
STEP = {"theta0_cable_deg": 0.05, "d_x_mm": 0.2}   # absolute steps of the signed parameters
REL = 0.03                                          # relative step of the others (on the logarithm)


def effective_samples(r: np.ndarray) -> float:
    """Number of independent samples of a residual history: n / (1 + 2 sum of its autocorrelation up to the
    first zero crossing)."""
    r = r - r.mean()
    if np.allclose(r, 0):
        return 1.0
    ac = np.correlate(r, r, "full")[r.size - 1:] / np.dot(r, r)
    k0 = np.argmax(ac < 0) if np.any(ac < 0) else ac.size
    return float(max(1.0, r.size / (1 + 2 * np.sum(ac[1:k0]))))


def main() -> None:
    res_dir = HERE / "results_identifiability"
    res_dir.mkdir(exist_ok=True)
    md = ["# Identifiability of the parameters of the current-driven model", "",
          "Jacobian of the normalised position residuals on the identification moves (1 deg sequences) at the "
          "fitted point; fit coordinates: ln p for positive parameters, p for signed ones. Condition number of the "
          "Jacobian with columns scaled to unit norm; correlations and standard errors of the estimates from "
          "(J^T J)^-1, with the residual variance and the effective number of independent samples (autocorrelation "
          "of the residuals). Relative standard error: of p for ln-parameters, absolute for signed ones.", ""]
    js = {}
    for label, cfg_name, res_name, axes in CASES:
        cc = yaml.safe_load((HERE / cfg_name).read_text())
        cfg = yaml.safe_load((HERE / cc["base_config"]).read_text())
        if cc.get("gravity"):
            cfg["plant"]["balanced_telescope"] = False
        final = json.loads((HERE / res_name / "summary.json").read_text())["final"]
        seq_dir = (HERE / cc["sequences_dir"]).resolve()
        cl, _, _ = mn.build(cfg, HERE)
        rm = pf.RingMap.from_config(cl.plant.config)
        seqs = {n: read_sequence(seq_dir / f"{n}.csv") for n in cc["identification_sequences"]}
        wid = [w for n in cc["identification_sequences"] for w in cs.windows(seqs[n], n, cc["window"]["pre_s"],
                                                                                cc["window"]["after_move_s"])]
        with mp.get_context("fork").Pool(8, initializer=cs._init,
                                         initargs=(cfg, str(seq_dir), 2.0, float(cc["input_lowpass_hz"]),
                                                   {"a": rm.a, "b": rm.b}, final)) as pool:
            for ax, names in axes.items():
                w_ax = [w for w in wid if w.axis == ax]
                base = {a_: dict(final[a_]) for a_ in AXES}

                def resid(vals: dict) -> np.ndarray:
                    return np.concatenate([cs.residual_of(r) for r in pool.map(cs.job, [(w, vals, 60.0) for w in w_ax])])

                r0 = resid(base)
                J = np.empty((r0.size, len(names)))
                for k, n in enumerate(names):
                    cols = []
                    for sg in (+1, -1):
                        v = {a_: dict(base[a_]) for a_ in AXES}
                        if n in SIGNED:
                            v[ax][n] = base[ax][n] + sg * STEP[n]
                        else:
                            v[ax][n] = base[ax][n] * np.exp(sg * REL)
                        cols.append(resid(v))
                    h = 2 * (STEP[n] if n in SIGNED else REL)
                    J[:, k] = (cols[0] - cols[1]) / h
                # effective samples: residuals of each move
                n_win = len(w_ax)
                per = r0.size // n_win
                n_eff = sum(effective_samples(r0[i * per:(i + 1) * per]) for i in range(n_win))
                s2 = float(r0 @ r0) / max(1.0, n_eff - len(names))
                # information scaled by the effective count (independent samples carry the information)
                F = J.T @ J * (n_eff / r0.size)
                Js = J / np.linalg.norm(J, axis=0)
                sv = np.linalg.svd(Js, compute_uv=False)
                cov = s2 * np.linalg.pinv(F)
                sd = np.sqrt(np.clip(np.diag(cov), 0, None))
                corr = cov / np.outer(sd, sd)
                key = f"{label} / {ax}"
                js[key] = {"names": names, "singular_values": sv.tolist(), "condition": float(sv[0] / sv[-1]),
                           "n_samples": int(r0.size), "n_effective": n_eff, "corr": corr.tolist(), "sd": sd.tolist(),
                           "values": [base[ax][n] for n in names]}
                rows = []
                for k, n in enumerate(names):
                    v = base[ax][n]
                    se = f"{sd[k]:.2g} (abs.)" if n in SIGNED else f"{100 * sd[k]:.0f} %"
                    big = [names[j] for j in range(len(names)) if j != k and abs(corr[k, j]) > 0.9]
                    rows.append([n, f"{v:.4g}", se, ", ".join(f"{b} ({corr[k, names.index(b)]:+.2f})" for b in big) or "-"])
                md += [f"## {label}: {ax}", "",
                       f"{n_win} moves, {r0.size} samples, effective independent samples {n_eff:.0f}; condition number "
                       f"{sv[0] / sv[-1]:.3g}; scaled singular values " + ", ".join(f"{x:.3g}" for x in sv) + ".", "",
                       out.markdown_table(["parameter", "fitted value", "standard error", "correlations |r| > 0.9"], rows), ""]
                fig, a = plt.subplots(figsize=(4.8, 4.0))
                im = a.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
                a.set_xticks(range(len(names)), [n.replace("_{\\theta,h}", "_th") for n in names], rotation=60, fontsize=7)
                a.set_yticks(range(len(names)), [n.replace("_{\\theta,h}", "_th") for n in names], fontsize=7)
                for i in range(len(names)):
                    for j in range(len(names)):
                        a.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize=5.5)
                fig.colorbar(im, ax=a, label="correlation of the estimates [-]")
                a.set_title(f"{label}: {ax}", fontsize=8)
                fig.tight_layout()
                out.save(fig, res_dir / "figures", f"corr_{res_name}_{ax}")
                print(f"{key}: condition {sv[0] / sv[-1]:.3g}, n_eff {n_eff:.0f}", flush=True)
    (res_dir / "summary.md").write_text("\n".join(md))
    (res_dir / "summary.json").write_text(json.dumps(js, indent=2))
    print("done")


if __name__ == "__main__":
    main()
