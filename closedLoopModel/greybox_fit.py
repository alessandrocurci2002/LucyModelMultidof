r"""Grey-box estimate of the lumped plant on measured plant responses, with lucy_cl.greyest.

For each axis the parameters of :mod:`lucy_cl.gimbal_grey` are estimated on the WorkBench "Plant" trace (encoder
velocity over current command) of the measurement given in ``drives.<axis>.measurement``, within the bounds of the
``fit`` section (inertias and geometric springs within a relative band, the other springs and the damping free).
With ``fit.max_damping_ratio`` the damping ratio of every mode up to ``fit.damping_check_fmax_Hz`` is limited (soft
constraint, :meth:`PlantGreyModel.damping_penalty`). The search starts from every combination of the ``starts``;
the estimate with the lowest loss is kept.

Outputs in ``results_dir``: ``summary.md``/``summary.json`` (estimates, standard deviations, fit by band, modes,
starts), ``tables/parameters_<axis>.csv``, ``figures/bode_<axis>.pdf``, and ``<config>_fit.yaml``, a closed-loop
configuration with the estimated plant.

Run (from this folder)::

    ../.venv/bin/python -I greybox_fit.py --config config_vtesi_greybox.yaml
    ../.venv/bin/python -I greybox_fit.py --config config_vcl_greybox.yaml
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import main as mn  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl.gimbal_grey import GROUND, Bounds, PlantGreyModel, axis_parameters  # noqa: E402
from lucy_cl.greyest import FrequencyData, GreyestOptions, compare, greyest  # noqa: E402
from lucy_cl.plant_fit import RingMap  # noqa: E402

LABEL = {"azimuth": "azimuth", "elevation": "elevation"}


def measurement(path: Path) -> FrequencyData:
    f, tr, _ = spm.carica_bode(path)
    g, p = tr["Plant"]
    return FrequencyData(f, 10 ** (g / 20) * np.exp(1j * np.deg2rad(p)), name=path.name)


def check_fmax(fit: dict) -> float:
    """Upper frequency of the modes whose damping is limited (``null`` or absent: all the modes)."""
    v = fit.get("damping_check_fmax_Hz")
    return float("inf") if v is None else float(v)


def estimate_axis(cl, tests, axis: str, data: FrequencyData, ring_map, fit: dict, log) -> dict:
    model = PlantGreyModel(cl, axis, tests[axis], ring_map, axis_parameters(cl, axis, ring_map))
    # PyYAML reads 1.0e9 (no sign in the exponent) as a string: every number of the section is converted
    bounds = Bounds(**{k: (tuple(float(x) for x in v) if isinstance(v, list) else float(v))
                       for k, v in fit["bounds"].items()})
    opt = GreyestOptions(weighting_filter=list(fit["band_Hz"]), frequency_error=fit["frequency_error"],
                         loss=fit["loss"], f_scale=float(fit["f_scale"]), max_iterations=int(fit["max_iterations"]),
                         covariance=fit.get("covariance", "robust"))
    penalty = None
    if fit.get("max_damping_ratio") is not None:
        penalty = model.damping_penalty(float(fit["max_damping_ratio"]), check_fmax(fit),
                                        float(fit.get("damping_penalty_weight", 1.0e3)))
    nominal = model.grey_model(bounds, penalty=penalty)
    starts = {k: [float(x) for x in v] for k, v in fit.get("starts", {}).get(axis, {}).items()}
    combos = [dict(zip(starts, c)) for c in itertools.product(*starts.values())] or [{}]
    runs = []
    for c in combos:
        t0 = time.time()
        est = greyest(data, model.grey_model(bounds, start=c, penalty=penalty), opt)
        runs.append((est, c))
        log(f"  {axis}: start {c} -> loss {est.report.loss:.5g}, fit {est.report.fit_percent[0, 0]:.1f} %, "
            f"{est.report.n_evaluations} evaluations, {time.time() - t0:.0f} s")
    best, start = min(runs, key=lambda r: r[0].report.loss)
    model.apply(best.model.values)  # leave the plant at the estimate
    return {"model": model, "nominal": nominal, "best": best, "start": start, "runs": runs, "options": opt}


def bode_figure(axis: str, data: FrequencyData, curves: dict[str, np.ndarray], fits: dict[str, dict], label: str,
                band: tuple[float, float]):
    fig, (ax_m, ax_p) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax_m.plot(data.f, 20 * np.log10(np.abs(data.G[:, 0, 0])), color=pt.INCHIOSTRO, marker="o", ms=1.5, lw=0.8,
              label=f"measured: {data.name}")
    ax_p.plot(data.f, np.angle(data.G[:, 0, 0], deg=True), color=pt.INCHIOSTRO, marker="o", ms=1.5, lw=0.8)
    for k, (name, G) in enumerate(curves.items()):
        st = pt._stile(k)
        fb = fits[name]
        ax_m.plot(data.f, 20 * np.log10(np.abs(G)), lw=1.5, label=f"{name} (fit in band {fb['band']:.0f} %; 1-5 Hz "
                  f"{fb['1-5 Hz']:.0f}, 5-100 Hz {fb['5-100 Hz']:.0f} %)", **st)
        ax_p.plot(data.f, np.angle(G, deg=True), lw=1.2, **st)
    for ax in (ax_m, ax_p):
        ax.set_xscale("log")
        ax.set_xlim(data.f[0], data.f[-1])
        ax.axvspan(band[1], data.f[-1], color=pt.GRIGLIA, alpha=0.5, zorder=0)  # outside the estimation band
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_p.set_ylabel("phase [deg]")
    ax_p.set_ylim(-180, 180)
    ax_p.set_yticks(np.arange(-180, 181, 90))
    ax_p.set_xlabel("frequency [Hz]")
    ax_m.legend(fontsize=7, loc="lower left")
    fig.suptitle(f"{LABEL[axis].capitalize()} plant, encoder velocity / current: {label}\n(shaded: outside the "
                 f"estimation band {band[0]:g}-{band[1]:g} Hz)", fontsize=10)
    return fig


def main() -> None:
    a = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    a.add_argument("--config", type=Path, default=HERE / "config_vtesi_greybox.yaml")
    a.add_argument("--axes", nargs="+", default=list(AXES), choices=list(AXES))
    args = a.parse_args()
    cfg = yaml.safe_load(args.config.read_text())
    base = args.config.resolve().parent
    out = base / cfg["results_dir"]
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        log_lines.append(msg)

    pt.applica_stile()
    cl, tests, _ = mn.build(cfg, base)
    ring_map = RingMap.from_config((base / cfg["plant"]["config"]).resolve())
    fit = cfg["fit"]
    summary: dict = {"config": args.config.name, "plant": cfg["plant"]["config"], "fit": fit, "axes": {}}
    md = [f"# Grey-box estimate: {cfg['label']} ({args.config.name})", "",
          f"Model `{cfg['plant']['config']}`, estimator `lucy_cl.greyest`: error {fit['frequency_error']}, loss "
          f"{fit['loss']} (f_scale {fit['f_scale']}), band {fit['band_Hz'][0]:g}-{fit['band_Hz'][1]:g} Hz, "
          f"covariance {fit.get('covariance', 'robust')}. Bounds: inertias and geometric springs "
          f"+-{100 * float(fit['bounds']['relative']):.0f} %, the other springs free; damping "
          + (f"limited to zeta <= {fit['max_damping_ratio']} on every mode up to "
             f"{check_fmax(fit):g} Hz (ground loss factor excluded)."
             if fit.get("max_damping_ratio") is not None else "free."), ""]
    overrides: dict[str, float] = {}
    ground: dict[str, dict] = {}
    for axis in args.axes:
        data = measurement((base / cfg["drives"][axis]["measurement"]).resolve())
        log(f"{axis}: {data.name}, {len(data.f)} points")
        r = estimate_axis(cl, tests, axis, data, ring_map, fit, log)
        best, rep, model = r["best"], r["best"].report, r["model"]
        # responses: nominal (configuration) and estimate, on the measurement points
        G_nom = compare(data, r["nominal"])[0][0][:, 0, 0]
        model.apply(best.model.values)
        G_est = best.model.frequency_response(data.w)[:, 0, 0]
        b0, b1 = (float(x) for x in fit["band_Hz"])
        inb = (data.f >= b0) & (data.f <= b1)
        fits = {n: {**spm.fit_bande(data.f, data.G[:, 0, 0], G), "band": spm.fit_nrmse(data.G[inb, 0, 0], G[inb])}
                for n, G in (("nominal", G_nom), ("grey-box estimate", G_est))}
        pt.salva(bode_figure(axis, data, {"nominal": G_nom, "grey-box estimate": G_est}, fits, cfg["label"], (b0, b1)),
                 out / "figures", f"bode_{axis}")
        f_chk = check_fmax(fit)
        modes = model.modes(best.model.values, f_chk)
        bnd = {p.name: (p.minimum, p.maximum) for p in r["nominal"].parameters}
        kinds = {p.name: (p.kind, p.nominal, p.note) for p in model.params}
        rows = []
        for p in rep.parameters:
            kind, nom, note = kinds[p.name]
            rel = p.std / abs(p.value) if np.isfinite(p.std) and p.value else float("nan")
            rows.append({"name": p.name, "kind": kind, "nominal": nom, "min": bnd[p.name][0], "max": bnd[p.name][1],
                         "start": p.initial, "estimate": p.value, "ratio_to_nominal": p.value / nom, "std": p.std,
                         "std_rel": rel, "at_bound": p.at_bound, "note": note})
        with open(out / "tables" / f"parameters_{axis}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        undetermined = [x["name"] for x in rows if not (x["std_rel"] < 0.5)]
        summary["axes"][axis] = {
            "measurement": data.name, "start": r["start"], "loss": rep.loss, "loss_start": rep.loss_start,
            "fit_percent": fits, "status": rep.status, "evaluations": rep.n_evaluations,
            "jacobian_condition": rep.jacobian_condition, "parameters": rows,
            "modes_Hz_zeta": [[round(f, 3), round(z, 4)] for f, z in modes],
            "correlation": {"names": list(rep.free_names), "matrix": np.round(rep.correlation, 3).tolist()},
            "runs": [{"start": c, "loss": e.report.loss, "fit": float(e.report.fit_percent[0, 0])} for e, c in r["runs"]],
            "undetermined": undetermined}
        md += [f"## {LABEL[axis].capitalize()}", "",
               f"Measurement `{data.name}`; best start {r['start']}; {rep.status}, {rep.n_evaluations} evaluations; "
               f"loss {rep.loss_start:.4g} -> {rep.loss:.4g}; Jacobian condition {rep.jacobian_condition:.3g}.", "",
               f"| fit [%] | band {b0:g}-{b1:g} Hz | 1-5 Hz | 5-100 Hz | >100 Hz |", "|---|---|---|---|---|"]
        md += [f"| {n} | " + " | ".join(f"{v[b]:.1f}" for b in ("band", "1-5 Hz", "5-100 Hz", ">100 Hz")) + " |"
               for n, v in fits.items()]
        md += ["", f"Modes of the axis in the estimated model up to {f_chk:g} Hz (ground stiffness included; zeta 1 = "
               "overdamped): "
               + ", ".join(f"{f:.1f} Hz (zeta {z:.3f})" for f, z in modes) + "."]
        md += ["", "| parameter | kind | nominal | bounds | estimate | estimate/nominal | std/estimate | note |",
               "|---|---|---|---|---|---|---|---|"]
        for x in rows:
            md.append(f"| `{x['name']}` | {x['kind']} | {x['nominal']:.4g} | {x['min']:.3g} - {x['max']:.3g} | "
                      f"{x['estimate']:.4g}{' **(at bound)**' if x['at_bound'] else ''} | {x['ratio_to_nominal']:.3g} | "
                      f"{x['std_rel']:.1%} | {x['note']} |")
        md += ["", "Starts: " + "; ".join(f"{c} -> loss {e.report.loss:.4g}" for e, c in r["runs"]), "",
               "Not determined by the data (std above 50 % of the estimate): "
               + (", ".join(f"`{n}`" for n in undetermined) or "none") + ".", ""]
        log(rep.summary())
        # closed-loop configuration with the estimate: lumped values written as overrides
        pl = cl.plant
        tag = "psi" if axis == "azimuth" else "theta"
        for s in pl._sym:
            if s.startswith((rf"J_{{\{tag}", rf"k_{{\{tag}", rf"c_{{\{tag}")) or s in (
                    ("J^{ax}", "I_zz") if axis == "azimuth" else ("I_yy",)):
                overrides[s] = float(pl.value(s))
        g = pl.ground(axis)
        ground[axis] = {"k": float(g[0]), "eta": float(g[1])}

    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    (out / "summary.md").write_text("\n".join(md) + "\n")
    (out / "run.log").write_text("\n".join(log_lines) + "\n")
    fit_cfg = {k: v for k, v in cfg.items() if k != "fit"}
    fit_cfg["results_dir"] = cfg["results_dir"] + "_closed_loop"
    fit_cfg["label"] = cfg["label"] + ", estimated values (greybox_fit.py)"
    fit_cfg["plant"] = dict(cfg["plant"], parameter_overrides=overrides, ground_override=ground)
    (base / f"{args.config.stem}_fit.yaml").write_text(
        f"# Closed-loop configuration with the plant estimated by greybox_fit.py from {args.config.name}.\n"
        + yaml.safe_dump(fit_cfg, sort_keys=False, allow_unicode=True))
    log(f"results in {out}")


if __name__ == "__main__":
    main()
