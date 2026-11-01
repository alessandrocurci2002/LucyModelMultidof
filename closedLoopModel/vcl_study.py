r"""Study on the October closed-loop tests (LucyIdentification/VClosedLoop as the only source of truth).

1-2. Transfer functions of the elevation and azimuth controllers from ``Parameters_*.akd`` and the structure of
     the AKD2G loops, checked on the WorkBench traces of ``BodeMeasurement_*.csv``.
3.   Elevation plant: fit of the hub - telescope spring k_b, of the ring stiffness factor and of the ground spring.
4.   Azimuth plant: fit of the ground spring.
5.   Both axes: fit of the dampers.
6.   Validation of the closed loop on ``MPCLucy_*.csv`` with the NRMSE fit of the time histories; when it is not
     adequate the controller is reconsidered (back to 1) through the iterations of the configuration.

Results in ``<results_dir>/``: ``figures/``, ``summary.md``, ``summary.json``, ``plant_fitted.yaml`` and the
closed-loop configuration ``config_vcl_fit.yaml`` (fitted plant) in this folder.

Usage (from this folder)::

    ../.venv/bin/python -I vcl_study.py --config config_vcl_study.yaml
"""

from __future__ import annotations

import argparse
import copy
import json
import multiprocessing as mp
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
from lucy_cl import AXES  # noqa: E402
from lucy_cl import controller_tf as ctf  # noqa: E402
from lucy_cl import identification as idn  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl import plant_fit as pf  # noqa: E402
from lucy_cl.closed_loop import loop_frf  # noqa: E402
from lucy_cl.data_files import read_workbench  # noqa: E402
from lucy_cl.driver import DriveParams  # noqa: E402
from lucy_cl.validation import read_sequence  # noqa: E402

MDEG = out.MDEG
plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 6.5,
                     "savefig.bbox": "tight"})


def fit_nrmse(meas: np.ndarray, model: np.ndarray) -> float:
    from gimbal_mdof.sperimentale import fit_nrmse as f

    ok = np.isfinite(model)
    return f(meas[ok], model[ok]) if ok.sum() >= 10 else float("nan")


def fit_bands(f: np.ndarray, G: np.ndarray, Gm: np.ndarray) -> dict[str, float]:
    from gimbal_mdof.sperimentale import fit_bande

    return {k: float(v) for k, v in fit_bande(f, G, Gm).items()}


def poly_str(c: np.ndarray) -> str:
    return "[" + ", ".join(f"{x:.6g}" for x in c) + "]"


# --------------------------------------------------------------------------- steps 1-2


def controllers(cl, cfg: dict, wb: dict, fig_dir: Path, md: list, js: dict) -> None:
    md += ["## 1-2. Controller transfer functions", "",
           "Built from `Parameters_*.akd` with the structure of the AKD2G loops (series PI velocity loop, biquads "
           "AR1-AR2 forward and AR3-AR4 on the feedback, observer as a first-order low pass of bandwidth OBS.BW, "
           "proportional position loop with PL.KI = 0, current PI with its zero on R/L). Control law "
           "I_cmd = C_v H_fwd [C_p (P_cmd - P_fb) - H_fb s P_fb]; C_e = C_v H_fwd C_p from the position error, "
           "C_fb = C_v H_fwd (C_p + s H_fb) from the encoder position, C_WB = C_v H_fwd (1 + 2 pi PL.KP / s) the "
           "WorkBench Controller trace. Polynomials in s, descending powers; A per rad, A per rad/s.", ""]
    js["controllers"] = {}
    rows_chk = []
    for ax in AXES:
        akd_path = (HERE / cfg["drives"][ax]["parameters"]).resolve()
        w = wb[ax]
        # tuning recorded in the Bode export against the parameter file
        p_file = DriveParams.from_files(ax, akd_path, pl_kp_unit="hz")
        p_test = DriveParams.from_files(ax, akd_path, w.path, pl_kp_unit="hz")
        same = {k: (getattr(p_file, k), getattr(p_test, k)) for k in ("Kpp", "Kpv", "Kiv_hz")}
        tfs = ctf.controller_tfs(p_file)
        chk = ctf.check_workbench(tfs, w.f, w.response("Controller"), w.response("Feedback Biquad"))
        rows_chk.append([ax, f"{p_file.Kpp / (2 * np.pi):g}", f"{p_file.Kpv:g}", f"{p_file.Kiv_hz:g}",
                         "yes" if all(np.isclose(a, b) for a, b in same.values()) else f"no: {same}",
                         f"{chk['controller_dB']:.4f} dB / {chk['controller_deg']:.3f} deg",
                         f"{chk['feedback_dB']:.4f} dB / {chk['feedback_deg']:.2f} deg"])
        blocks = []
        for name in ("C_p", "C_v", "H_fwd", "H_AR34", "H_obs", "C_i", "C_e", "C_fb", "C_WB"):
            t = tfs[name]
            blocks.append([name, f"{t.order[0]}/{t.order[1]}", poly_str(t.num), poly_str(t.den),
                           "; ".join(ctf.describe_roots(t.zeros())) or "-", "; ".join(ctf.describe_roots(t.poles())) or "-"])
        md += [f"### {ax}", "", out.markdown_table(["block", "order num/den", "numerator", "denominator", "zeros",
                                                    "poles"], blocks), ""]
        js["controllers"][ax] = {name: {"num": t.num.tolist(), "den": t.den.tolist()} for name, t in tfs.items()}
        js["controllers"][ax]["check"] = chk
        # figure
        f = np.logspace(-1, np.log10(2000), 1500)
        s = 2j * np.pi * f
        p_s = DriveParams.from_files(ax, akd_path, pl_kp_unit="per_s")
        tfs_s = ctf.controller_tfs(p_s)
        fig, axs = plt.subplots(2, 2, figsize=(10, 5.6), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
        for j, (title, curves, meas) in enumerate((
                ("WorkBench traces", [("C_WB, model", tfs["C_WB"], "C0"), ("H_AR3 H_AR4, model", tfs["H_AR34"], "C1")],
                 [("Controller trace", w.response("Controller")), ("Feedback Biquad trace", w.response("Feedback Biquad"))]),
                ("controller of the loop", [("C_e, K_pp = 2 pi PL.KP", tfs["C_e"], "C0"), ("C_e, K_pp = PL.KP [1/s]", tfs_s["C_e"], "C3"),
                                             ("C_fb / s (from encoder velocity)", tfs["C_fb"] * ctf.TF(np.array([1.0]), np.array([1.0, 0.0])), "C2"),
                                             ("H_fb (AR3 AR4 observer)", tfs["H_fb"], "C1")], []))):
            for lab, t, c in curves:
                H = t(s)
                axs[0, j].semilogx(f, 20 * np.log10(np.abs(H)), color=c, label=lab)
                axs[1, j].semilogx(f, np.angle(H, deg=True), color=c, label=lab)
            for k, (lab, H) in enumerate(meas):
                axs[0, j].semilogx(w.f, 20 * np.log10(np.abs(H)), ".", color=f"C{k}", ms=4, mfc="none", label=lab)
                axs[1, j].semilogx(w.f, np.angle(H, deg=True), ".", color=f"C{k}", ms=4, mfc="none", label=lab)
            axs[0, j].set_title(f"{ax}: {title}", fontsize=8)
            axs[0, j].set_ylabel("magnitude [dB re 1 A/rad, A/(rad/s) or 1]")
            axs[1, j].set_ylabel("phase [deg]")
            axs[1, j].set_xlabel("frequency [Hz]")
            for a in axs[:, j]:
                a.legend(loc="best")
        fig.tight_layout()
        out.save(fig, fig_dir, f"controller_{ax}")
    md += ["### Check against the WorkBench traces", "",
           out.markdown_table(["axis", "PL.KP", "VL.KP [A/(rad/s)]", "VL.KI [Hz]", "tuning of the export = file",
                               "Controller: max error", "Feedback Biquad: max error"], rows_chk), "",
           "The Bode export records the same gains and filters as the parameter file. The units of the files differ "
           "(azimuth positions in degrees, elevation in millidegrees); the gains of the loops do not depend on them "
           "(VL.KP in A/(rad/s), PL.KP in velocity unit per position unit).", ""]


# --------------------------------------------------------------------------- steps 3-5


def plant_fit(cl, sc: dict, wb: dict, fig_dir: Path, md: list, js: dict) -> dict:
    pfc = sc["plant_fit"]
    rm = pf.RingMap.from_config(cl.plant.config)
    dm = pf.DampingMap.from_plant(cl.plant)
    pl = cl.plant
    # float(): PyYAML reads exponents without sign (1.0e4) as strings
    bounds = {k: tuple(float(x) for x in v) for k, v in pfc["bounds"].items() if k != "damper_factor"}
    dfac = [float(x) for x in pfc["bounds"]["damper_factor"]]
    c0 = {n: pl.value(n) for n in pl._sym if n.startswith("c_")}
    for n, v in c0.items():
        bounds[n] = (v * dfac[0], v * dfac[1])
    prior = {n: (v, float(pfc.get("damper_prior_weight", 0.3))) for n, v in c0.items()}
    f_plot = np.logspace(0, np.log10(2000), 1500)
    G0 = {ax: loop_frf(cl, ax, f_plot).plant for ax in AXES}
    fits0 = {ax: fit_bands(wb[ax].f, wb[ax].response("Plant"), loop_frf(cl, ax, wb[ax].f).plant) for ax in AXES}
    js["plant_fit"] = {"start": {}, "stages": {}}
    rows = []
    for ax in ("elevation", "azimuth"):
        a = pfc[ax]
        w = wb[ax]
        G = w.response("Plant")
        tag = "theta" if ax == "elevation" else "psi"
        dampers = [n for n in c0 if f"\\{tag}" in n]
        js["plant_fit"]["start"][ax] = pf.current_values(cl, a["stiffness"] + dampers, rm)
        stages = (("A: stiffness and ground (steps 3-4)", a["stiffness"], tuple(a["band_stiffness"]), None),
                  ("B: dampers (step 5)", dampers, tuple(a["band_damping"]), prior),
                  ("C: joint refinement", a["stiffness"] + dampers, tuple(a["band_damping"]), prior))
        js["plant_fit"]["stages"][ax] = []
        for label, names, band, pr in stages:
            r = pf.fit(cl, ax, w.f, G, names, bounds, band, rm, loss=pfc.get("loss", "soft_l1"),
                       f_scale=float(pfc.get("f_scale", 0.5)), prior=pr)
            fb = fit_bands(w.f, G, loop_frf(cl, ax, w.f).plant)
            js["plant_fit"]["stages"][ax].append({"stage": label, "band": band, "values": r.values,
                                                  "cost": [r.cost_start, r.cost], "n_eval": r.n_eval, "fits": fb})
            rows.append([ax, label, f"{band[0]:g}-{band[1]:g}",
                         "; ".join(f"{n} {r.start[n]:.4g} -> {v:.4g}" for n, v in r.values.items()),
                         f"{r.cost_start:.1f} -> {r.cost:.1f}"] + [f"{fb[k]:.0f}" for k in ("all", "1-5 Hz", "5-100 Hz", ">100 Hz")])
    md += ["## 3-5. Plant parameters on the October frequency responses", "",
           "Measured plant: encoder velocity over current command of `BodeMeasurement_*.csv`; model: H6 with the "
           "motor, the current loop and the ground constraint k(1 + j eta). Residuals ln(G_model/G_measured) "
           "(magnitude and phase) on the points of the band, robust loss; parameters on their logarithms. Ring: "
           "stiffness factor kappa_r of gimbal_mdof (1 = ring alone, H6 = 1000), acting on k_theta,1 (encoder - "
           "hub) and k_theta,2 (hub - motor) as a compliance in series. Dampers: element dampers of the lumped "
           "model, with a weak prior on their start values (H6, zeta_1 = 0.10) for those the band does not "
           "determine.", "",
           "Fit NRMSE of the plant before the fit (H6): " + "; ".join(
               f"{ax} " + ", ".join(f"{k} {v:.0f} %" for k, v in fits0[ax].items()) for ax in AXES), "",
           out.markdown_table(["axis", "stage", "band [Hz]", "parameters (start -> fit)", "cost",
                               "fit all [%]", "1-5 Hz", "5-100 Hz", ">100 Hz"], rows), ""]
    # figure: measured, H6, fitted
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 5.8), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    for j, ax in enumerate(AXES):
        w = wb[ax]
        G = w.response("Plant")
        Gf = loop_frf(cl, ax, f_plot).plant
        curves = [(f"measured (BodeMeasurement_{ax}.csv)", w.f, G, dict(color="k", marker=".", ls="none", ms=3)),
                  ("model H6 (zeta1 = 0.10)", f_plot, G0[ax], dict(color="C3", ls="--")),
                  ("model, fitted parameters", f_plot, Gf, dict(color="C0"))]
        for lab, f_, H, sty in curves:
            axs[0, j].semilogx(f_, 20 * np.log10(np.abs(H)), label=lab, **sty)
            axs[1, j].semilogx(f_, np.angle(H, deg=True), label=lab, **sty)
        a = sc["plant_fit"][ax]
        for x in a["band_damping"]:
            axs[0, j].axvline(x, color="0.6", lw=0.6, ls=":")
        axs[0, j].set_title(f"{ax}: plant v_fb / i_cmd (dotted: fit band)", fontsize=8)
        axs[0, j].set_ylabel("|P| [dB re 1 (rad/s)/A]")
        axs[1, j].set_ylabel("phase [deg]")
        axs[1, j].set_xlabel("frequency [Hz]")
        for a_ in axs[:, j]:
            a_.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "plant_fit")
    final = {}
    for ax in AXES:
        for st in js["plant_fit"]["stages"][ax]:
            final.update(st["values"])
    js["plant_fit"]["final"] = final
    return {"values": final, "ring_map": rm, "c0": c0}


def fitted_config(cfg: dict, cl, final: dict, rm) -> dict:
    """Closed-loop configuration with the fitted plant (lumped overrides and ground constraint)."""
    c = copy.deepcopy(cfg)
    pl = cl.plant
    over = {pf.K4: float(final["k_b"]), pf.K1: float(rm.k(pf.K1, final["ring"])), pf.K2: float(rm.k(pf.K2, final["ring"]))}
    for n, v in final.items():
        if n.startswith("c_"):
            over[n] = float(v)
    c["plant"]["parameter_overrides"] = over
    c["plant"]["ground_override"] = {"azimuth": {"k": float(final["kg_az"]), "eta": float(final["eta_az"])},
                                     "elevation": {"k": float(final["kg_el"]), "eta": float(final["eta_el"])}}
    # the Dahl elements take the amplitude of the October tests, on which the ground spring is fitted
    c["plant"]["bearing_friction"].pop("amplitude_measurement", None)
    c["results_dir"] = "results_vcl_fit"
    c["label"] = "H6 with the plant parameters fitted on the October frequency responses (vcl_study.py)"
    _ = pl
    return c


# --------------------------------------------------------------------------- step 6


def windows(seqs: dict, names: list[str], pre: float, post: float) -> list:
    out_ = []
    for n in names:
        for ax in AXES:
            out_ += idn.windows_for_moves(seqs[n], ax, pre=pre, post=post)
    return out_


def evaluate(pool, wins: list, values: dict, mc: float) -> list[dict]:
    res = pool.map(idn.simulate_job, [(w, values, mc) for w in wins])
    rec = []
    for w, o in zip(wins, res):
        ax = w.axis
        rec.append({"sequence": w.sequence, "axis": ax, "move": w.move,
                    "fit_error": fit_nrmse(o[ax]["error_meas"], o[ax]["error"]),
                    "fit_current": fit_nrmse(o[ax]["current_meas"], o[ax]["current"]),
                    "fit_position": fit_nrmse(o[ax]["enc_meas"], o[ax]["enc"]), "ok": o["ok"], "_res": o})
    return rec


def medians(rec: list[dict]) -> dict:
    out_ = {}
    for ax in AXES:
        sel = [r for r in rec if r["axis"] == ax]
        for k in ("fit_error", "fit_current", "fit_position"):
            x = np.array([r[k] for r in sel], dtype=float)
            out_[(ax, k)] = (float(np.nanmedian(x)), float(np.nanmin(x)), float(np.nanmax(x)), len(sel))
    return out_


def window_figure(sel: list[tuple[str, dict]], name: str, fig_dir: Path, title: str) -> None:
    fig, axs = plt.subplots(len(sel), 2, figsize=(9.0, 2.3 * len(sel) + 0.4), squeeze=False)
    for k, (lab, r) in enumerate(sel):
        o, ax = r["_res"], r["axis"]
        t = o["t"]
        for j, (key, sc, ylab) in enumerate((("error", MDEG, "command - encoder [mdeg]"), ("current", 1.0, "current [A rms]"))):
            a = axs[k, j]
            a.plot(t, sc * o[ax][f"{key}_meas"], "k", lw=1.5, label="measured")
            a.plot(t, sc * o[ax][key], "C0", lw=1.1, label=f"model (fit {r['fit_' + key]:.0f} %)")
            a.set_title(f"{lab}: {ax}, move at {r['move']:.2f} s", fontsize=7.5)
            a.set_xlabel("time from window start [s]")
            a.set_ylabel(ylab)
            a.legend(loc="best")
    fig.suptitle(title, fontsize=8.5)
    fig.tight_layout()
    out.save(fig, fig_dir, name)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config_vcl_study.yaml")
    args = ap.parse_args()
    t_start = time.time()
    sc_path = Path(args.config).resolve()
    sc = yaml.safe_load(sc_path.read_text())
    cfg = yaml.safe_load((HERE / sc["base_config"]).read_text())
    res_dir = HERE / sc["results_dir"]
    fig_dir = res_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    seq_dir = (HERE / sc["sequences_dir"]).resolve()
    amps = float(sc.get("amps_per_A", 2.0))
    md = ["# Study on the October closed-loop tests (VClosedLoop)", "",
          f"Configuration `{sc_path.name}`; base model `{sc['base_config']}`. Generated by `vcl_study.py`. Only the "
          "files of `LucyIdentification/VClosedLoop` are used.", ""]
    js: dict = {"config": sc_path.name}

    cl, _, _ = mn.build(cfg, HERE)
    wb = {ax: read_workbench((HERE / cfg["drives"][ax]["measurement"]).resolve()) for ax in AXES}
    controllers(cl, cfg, wb, fig_dir, md, js)
    print("  steps 1-2 done", flush=True)
    pfr = plant_fit(cl, sc, wb, fig_dir, md, js)
    cfg_fit = fitted_config(cfg, cl, pfr["values"], pfr["ring_map"])
    (res_dir / "plant_fitted.yaml").write_text(yaml.safe_dump({"values": {k: float(v) for k, v in pfr["values"].items()},
                                                               "plant": {k: cfg_fit["plant"][k] for k in ("parameter_overrides", "ground_override")}},
                                                              sort_keys=False, allow_unicode=True))
    header = ("# Closed-loop configuration with the plant fitted on the October frequency responses\n"
              f"# (written by vcl_study.py from {sc_path.name}; base {sc['base_config']}).\n")
    (HERE / "config_vcl_fit.yaml").write_text(header + yaml.safe_dump(cfg_fit, sort_keys=False, allow_unicode=True))
    print("  steps 3-5 done", flush=True)

    # ---- 6. validation in the time domain, iterations on the controller -------------------------------
    vc = sc["validation"]
    seqs = {n: read_sequence(seq_dir / f"{n}.csv", amps) for n in vc["sequences"]}
    thr = vc["threshold"]
    mc = float(vc.get("max_cpu_per_window_s", 120))
    cl_fit, _, _ = mn.build(cfg_fit, HERE)
    friction = {el.axis: (el.sigma0, el.Tc, el.amplitude) for el in cl_fit.plant.friction}
    js["friction"] = friction
    eig = {}
    for unit in ("hz", "per_s"):
        cfg_u = copy.deepcopy(cfg_fit)
        cfg_u["drives"]["position_gain_unit"] = unit
        c_u, _, _ = mn.build(cfg_u, HERE)
        A, _, _ = c_u.linearise()
        eig[unit] = float(np.linalg.eigvals(A).real.max())
    js["eig_max_real"] = eig
    md += ["## 6. Validation of the closed loop on the MPCLucy sequences", "",
           f"Fitted plant (`config_vcl_fit.yaml`); bearing friction as Dahl elements equivalent to the fitted ground "
           f"springs at the amplitude of the October Bode tests: " + "; ".join(
               f"{ax} sigma0 {v[0]:.0f} N m/rad, Tc {v[1]:.2f} N m (amplitude {1e3 * v[2]:.2f} mrad)" for ax, v in friction.items())
           + f". Largest real part of the closed-loop eigenvalues: {eig['hz']:.2g} 1/s (2 pi PL.KP), {eig['per_s']:.2g} 1/s "
           "(PL.KP in 1/s).", "",
           f"Every move simulated from the measured state at rest with the measured commands of both axes; NRMSE fit "
           f"of the position error command - encoder and of the current of the moving axis (and of the encoder "
           f"position, for reference). Adequate: median fit of the position error >= {thr['error']:g} % and of the "
           f"current >= {thr['current']:g} % on each axis.", ""]
    js["iterations"] = []
    rows = []
    best = None
    with mp.get_context("fork").Pool(int(sc.get("workers", 4)), initializer=idn._init_worker,
                                     initargs=(cfg_fit, str(HERE), str(seq_dir), amps)) as pool:
        base_vals = {ax: idn.current_values(cl_fit, ax) for ax in AXES}
        prev_vals, prev_adequate = None, {ax: True for ax in AXES}
        for it in vc["iterations"]:
            vals = {}
            for ax in AXES:
                p = cl_fit.drives[ax].p
                kpp_file = p.Kpp_workbench / (2 * np.pi)
                vals[ax] = {"Kpp": 2 * np.pi * kpp_file if it["position_gain_unit"] == "hz" else kpp_file,
                            "Kpv": p.Kpv, "Kiv": p.Kiv_hz, "kvff": 0.0, "Tc": base_vals[ax]["Tc"]}
            test_names = list(vc["sequences"])
            ident = {}
            if "identify" in it:
                assert prev_vals is not None
                bad = it.get("axes") or [ax for ax in AXES if not prev_adequate[ax]]
                ids = it["identification_sequences"]
                test_names = [n for n in vc["sequences"] if n not in ids]
                mc_id = float(it.get("max_cpu_per_window_s", mc))
                for ax in bad:
                    wins = [w for w in windows(seqs, ids, vc["pre_s"], vc["post_s"]) if w.axis == ax]
                    combos, costs = [], []
                    tc0 = None
                    if "Tc" in it["identify"] and it.get("start_Tc") == "regression":
                        # Coulomb level of the moves (rigid-body balance, model inertia) as the start
                        from lucy_cl.validation import rigid_body_regression
                        Mx, _, _ = cl_fit.plant.linear_matrices()
                        n_psi = cl_fit.plant.sim.struttura.n_psi
                        sl = slice(0, n_psi) if ax == "azimuth" else slice(n_psi, cl_fit.plant.n)
                        Jax = float(np.sum(Mx[sl, sl]))
                        tcs = [rigid_body_regression(seqs[n], ax, Jax)["coulomb_Nm"] for n in ids]
                        tc0 = float(np.clip(np.median(tcs), 0.05, 2.0))
                        print(f"  {it['name']}: {ax} start Tc {tc0:.3f} N m from the rigid-body balance {tcs}",
                              flush=True)
                    if "grid" in it:
                        # start of the least squares: best point of a coarse grid of the AKD2G gains (the gains
                        # of the file can give an unstable loop, where the local search does not move), with the
                        # friction level of the start when it is identified too
                        grid = it["grid"][ax]
                        combos = [dict(zip(grid, c)) for c in np.array(np.meshgrid(*grid.values())).T.reshape(-1, len(grid))]
                        if tc0 is not None:
                            combos = [c | {"Tc": tc0} for c in combos]
                        jobs = [(w, {a_: dict(vals[a_]) | (c if a_ == ax else {}) for a_ in AXES}, mc_id)
                                for c in combos for w in wins]
                        resg = pool.map(idn.simulate_job, jobs)
                        for k, c in enumerate(combos):
                            r_ = resg[k * len(wins):(k + 1) * len(wins)]
                            costs.append(sum(float(np.sum(idn.normalised_residual(o[ax]["error_meas"], o[ax]["error"]) ** 2)
                                                   + np.sum(idn.normalised_residual(o[ax]["current_meas"], o[ax]["current"]) ** 2))
                                             for o in r_))
                        k_best = int(np.argmin(costs))
                        start = {k: float(v) for k, v in combos[k_best].items()}
                        print(f"  {it['name']}: {ax} grid of {len(combos)} points, best {start} cost {costs[k_best]:.3f}",
                              flush=True)
                    else:
                        # start from the previous iteration
                        start = {k: float(prev_vals[ax][k]) for k in it["identify"] if k in prev_vals[ax]}
                        if tc0 is not None:
                            start["Tc"] = tc0
                    e = idn.estimate(pool, ax, wins, start, max_nfev=int(it.get("max_nfev", 30)), max_cpu=mc_id,
                                     log=lambda s_: print(s_, flush=True), fixed=vals)
                    ident[ax] = {"grid": {"points": [{k: float(v) for k, v in c.items()} for c in combos],
                                          "costs": costs}, "start": e.start, "values": e.values,
                                 "cost": [e.cost_start, e.cost], "windows": [w.label for w in wins]}
                    vals[ax].update(e.values)
            wins = windows(seqs, test_names, vc["pre_s"], vc["post_s"])
            rec = evaluate(pool, wins, vals, mc)
            med = medians(rec)
            adequate = {ax: med[(ax, "fit_error")][0] >= thr["error"] and med[(ax, "fit_current")][0] >= thr["current"]
                        for ax in AXES}
            prev_adequate, prev_vals = adequate, vals
            js["iterations"].append({"name": it["name"], "label": it["label"], "values": vals, "identified": ident,
                                     "test_sequences": test_names,
                                     "medians": {f"{a}/{k}": v for (a, k), v in med.items()},
                                     "adequate": adequate,
                                     "moves": [{k: v for k, v in r.items() if k != "_res"} for r in rec]})
            for ax in AXES:
                e_, c_, p_ = med[(ax, "fit_error")], med[(ax, "fit_current")], med[(ax, "fit_position")]
                rows.append([it["name"], ax, f"{vals[ax]['Kpp']:.3g}", f"{vals[ax]['Kpv']:.3g}", f"{vals[ax]['Kiv']:.3g}",
                             str(e_[3]), f"{e_[0]:.0f} [{e_[1]:.0f}, {e_[2]:.0f}]", f"{c_[0]:.0f} [{c_[1]:.0f}, {c_[2]:.0f}]",
                             f"{p_[0]:.1f}", "yes" if adequate[ax] else "no"])
            # figure: a 1 deg and a 2 deg move of each axis
            sel = []
            for ax in AXES:
                for amp_seq in [n for n in test_names if "step_1deg" in n or "1deg" in n][:1] + [n for n in test_names if "2deg" in n][:1]:
                    cand = [r for r in rec if r["axis"] == ax and r["sequence"] == amp_seq]
                    if cand:
                        sel.append((amp_seq.replace("MPCLucy_", ""), cand[0]))
            window_figure(sel, f"validation_{it['name']}", fig_dir, f"{it['name']}: {it['label']}")
            print(f"  {it['name']}: " + "; ".join(f"{ax} error {med[(ax, 'fit_error')][0]:.0f} %, current "
                                                   f"{med[(ax, 'fit_current')][0]:.0f} %" for ax in AXES), flush=True)
            best = it["name"]
            if all(adequate.values()):
                break
    md += [out.markdown_table(["iteration", "axis", "K_pp [1/s]", "K_pv [A/(rad/s)]", "K_iv [Hz]", "moves",
                               "fit error [%]", "fit current [%]", "fit position [%]", "adequate"], rows), ""]
    for itr in js["iterations"]:
        md += [f"- **{itr['name']}**: {itr['label']}; sequences: {', '.join(itr['test_sequences'])}."]
        for ax, d in itr["identified"].items():
            md += [f"  - {ax} identified on {', '.join(sorted({w.split(' ')[0] for w in d['windows']}))}: "
                   + ", ".join(f"{k} {d['start'][k]:.4g} -> {v:.4g}" for k, v in d["values"].items())
                   + f" (cost {d['cost'][0]:.2f} -> {d['cost'][1]:.2f})"]
    md += [""]
    # bar chart of the iterations
    its = js["iterations"]
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.0))
    for j, ax in enumerate(AXES):
        x = np.arange(len(its))
        for o, (k, c, lab) in enumerate((("fit_error", "C0", "position error"), ("fit_current", "C1", "current"))):
            v = np.array([i_["medians"][f"{ax}/{k}"][0] for i_ in its])
            axs[j].bar(x + (o - 0.5) * 0.38, np.clip(v, -100, 100), 0.36, color=c, label=f"{lab}, median")
            for xi, vi in zip(x, v):
                axs[j].text(xi + (o - 0.5) * 0.38, min(max(vi, -95), 95), f"{vi:.0f}", ha="center", fontsize=6)
        axs[j].axhline(thr["error"], color="C0", ls=":", lw=0.8, label=f"threshold error {thr['error']:g} %")
        axs[j].axhline(thr["current"], color="C1", ls=":", lw=0.8, label=f"threshold current {thr['current']:g} %")
        axs[j].set_xticks(x, [i_["name"] for i_ in its])
        axs[j].set_ylim(-100, 100)
        axs[j].set_ylabel("NRMSE fit [%]")
        axs[j].set_xlabel("iteration")
        axs[j].set_title(f"{ax}: MPCLucy moves", fontsize=8)
        axs[j].legend(loc="lower left")
    fig.tight_layout()
    out.save(fig, fig_dir, "iterations")
    js["last_iteration"] = best
    (res_dir / "summary.md").write_text("\n".join(md))
    (res_dir / "summary.json").write_text(json.dumps(js, indent=2, default=lambda o: float(o) if np.isscalar(o) else str(o)))
    print(f"results in {res_dir} ({time.time() - t_start:.0f} s)")


if __name__ == "__main__":
    main()
