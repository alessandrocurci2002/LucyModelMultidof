r"""Tracking of the circle with the signal constraints: EtherCAT at 1 ms and a MEMS gyro with Allan-variance noise.

The study of ``dob_tracking_study.py`` (current control on the encoders against Hilkert's gyro disturbance
observer fed forward to the current) repeated with:

* the drives on an EtherCAT network, process data exchanged every 1 ms (:mod:`lucy_cl.ethercat`): position
  targets interpolated by the drive, actual currents read by the master, observer in the master, feedforward
  current applied at the next cycle;
* the gyro of an Analog Devices ADIS16507-1 read over SPI (:mod:`lucy_cl.gyro`): second-order dynamics,
  2000 SPS, latency, noise from the Allan deviation (angle random walk, bias instability).

Same metrics on the line of sight of the telescope (pointing error, jitter above each band), compared with the
continuous loop and ideal gyro of ``dob_tracking_study.py``; Allan-deviation check of the gyro noise model.

Usage (from this folder)::

    ../.venv/bin/python -I dob_ethercat_study.py --config config_dob_ethercat.yaml [--replot]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy.signal import butter, filtfilt, welch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import dob_tracking_study as dt  # noqa: E402
from lucy_cl import AXES  # noqa: E402
from lucy_cl import output as out  # noqa: E402
from lucy_cl.ethercat import EtherCATParams, SampledObserverLoop  # noqa: E402
from lucy_cl.gyro import DEG, MemsGyro, allan_deviation  # noqa: E402

MDEG = out.MDEG
plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 6.5,
                     "savefig.bbox": "tight"})


def gyro_of(cc: dict, noise_scale: float = 1.0) -> MemsGyro:
    return MemsGyro(**{k: (v if k == "name" else float(v)) for k, v in cc["gyro"].items()}, noise_scale=float(noise_scale))


def build(cc: dict, base: dict, case: dict):
    """Sampled loop of one case, its initial state and the reference."""
    loop0, x0b, circle = dt.build_loop(base, None, bool(case.get("feedforward")))
    obs = dt.build_loop(base, case["observer_hz"], False)[0].observers if case.get("observer_hz") else {}
    e = cc["ethercat"]
    ecat = EtherCATParams(float(e["cycle_s"]), str(case.get("interpolation", e["interpolation"])), int(e["ff_delay_cycles"]))
    gyro = gyro_of(cc, case.get("noise_scale", 1.0)) if obs else None
    loop = SampledObserverLoop(loop0.plant, loop0.motors, loop0.drives, obs, gyro, ecat)
    loop.torque_scale, loop.cable_spring = loop0.torque_scale, loop0.cable_spring
    loop.q0[:] = loop0.q0
    x0 = loop.initial_state()
    x0[:loop0.n_states] = x0b
    return loop, x0, circle


def run_case(args) -> dict:
    """Simulation of one case (worker): arrays saved to ``runs/<name>.npz``."""
    cc, base, case, run_dir = args
    loop, x0, circle = build(cc, base, case)
    sim = base["simulation"]
    t_end = float(sim["t_settle_s"]) + float(sim["orbits"]) * 2 * np.pi / circle.omega
    T = loop.ecat.cycle_s
    t_end = round(t_end / T) * T
    tk = np.arange(int(round(t_end / T)) + 1) * T
    rng = np.random.default_rng(int(cc["seed"]))
    noise = {ax: loop.gyro.output_noise(tk, rng) for ax in loop.observers} if loop.observers else {}
    r = loop.simulate_sampled({ax: dt.CircleCommand(circle, ax) for ax in AXES}, t_end, x0, noise=noise,
                              rtol=float(cc["rtol"]))
    t, X, lg = r["t"], r["x"], r["log"]
    q = X[loop.s_q]
    psi_ref, th_ref, _, _ = circle.angles(t)
    data = {"t": t, "psi_ref": psi_ref, "th_ref": th_ref, "psi_los": q[loop.jl], "th_los": q[loop.jp],
            "psi_enc": q[loop.ie["azimuth"]], "th_enc": q[loop.ie["elevation"]]}
    pl = loop.plant
    for j, ax in enumerate(AXES):
        a = ax[:2]
        jb = pl.index(ax, "bearing")
        K, th0 = loop.cable_spring[ax]
        kf = [el.axis for el in pl.friction].index(ax)
        data.update({f"i_{a}": X[loop.i_cur[ax]], f"i_cmd_{a}": lg["i_cmd"][:, j], f"i_ff_{a}": lg["i_ff"][:, j],
                     f"Td_hat_{a}": lg["Td_hat"][:, j], f"w_gyro_{a}": lg["w_gyro"][:, j],
                     f"T_f_{a}": X[loop.s_f][kf], f"T_cable_{a}": -K * (q[jb] - th0),
                     f"noise_{a}": noise.get(ax, np.zeros(t.size))[:t.size]})
    np.savez_compressed(run_dir / f"{case['name']}.npz", **data)
    return {"name": case["name"], "success": bool(r["success"]), "message": r["message"], "cpu_s": r["cpu_s"],
            "n_rhs": r["n_rhs"]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config_dob_ethercat.yaml")
    ap.add_argument("--replot", action="store_true", help="metrics and figures from the saved runs")
    args = ap.parse_args()
    cc_path = Path(args.config).resolve()
    cc = yaml.safe_load(cc_path.read_text())
    base = yaml.safe_load((HERE / cc["base"]).read_text())
    res_dir = HERE / cc["results_dir"]
    run_dir, fig_dir = res_dir / "runs", res_dir / "figures"
    run_dir.mkdir(parents=True, exist_ok=True)
    js_path = res_dir / "summary.json"
    js = json.loads(js_path.read_text()) if js_path.exists() else {}
    cases = cc["cases"]
    t_start = time.time()

    if not args.replot:
        todo = [c for c in cases if not (run_dir / f"{c['name']}.npz").exists()]
        # longest cases first (noisy observers)
        todo.sort(key=lambda c: -(c.get("observer_hz") or 0) * (1.0 if c.get("noise_scale", 1.0) else 0.1))
        print(f"simulating {len(todo)} cases", flush=True)
        if todo:
            with mp.get_context("fork").Pool(min(len(todo), int(cc.get("workers", 8)))) as pool:
                for r in pool.imap_unordered(run_case, [(cc, base, c, run_dir) for c in todo]):
                    js.setdefault("runs", {})[r["name"]] = r
                    js_path.write_text(json.dumps(js, indent=1))
                    print(f"  {r['name']}: {r['message']} (cpu {r['cpu_s']:.0f} s)", flush=True)

    # ---- gyro noise model: Allan deviation and spectrum
    g = gyro_of(cc)
    rng = np.random.default_rng(int(cc["seed"]) + 1)
    x = g.odr_noise(float(cc["allan_check_s"]), rng, fs=4000.0)
    taus = np.logspace(np.log10(2 / g.odr_hz), np.log10(cc["allan_check_s"] / 10), 40)
    ta, ad = allan_deviation(x, g.odr_hz, taus)
    mdl = g.allan_model(ta)
    js["gyro"] = {"name": g.name, "output_noise_deg_s": float(np.std(x) / DEG), "N_rad_rts": g.N, "B_rad_s": g.B,
                  "allan_tau_s": ta.tolist(), "allan_sim_deg_h": (ad / DEG * 3600).tolist(),
                  "allan_model_deg_h": (mdl["total"] / DEG * 3600).tolist()}

    # ---- metrics
    loop, _, circle = build(cc, base, {"observer_hz": None})
    kt = {ax: loop.motors[ax].p.Kt for ax in AXES}
    t_settle = float(base["simulation"]["t_settle_s"])
    bands = [float(b) for b in cc["jitter_bands_hz"]]
    data = {c["name"]: dict(np.load(run_dir / f"{c['name']}.npz")) for c in cases}
    M = {}
    for n, d in data.items():
        d = {k: v[:-1] for k, v in d.items()}  # the last sample has no master signals
        data[n] = d
        M[n] = dt.metrics(d, circle, t_settle, bands, kt)
        m = d["t"] >= t_settle
        for ax in AXES:
            M[n][f"i_ff_applied_rms_{ax[:2]}_A"] = float(np.sqrt(np.mean(d[f"i_ff_{ax[:2]}"][m] ** 2)))
    js["metrics"] = M
    cont = {}
    cdir = HERE / cc["continuous_results"]
    if (cdir / "summary.json").exists():
        cm = json.loads((cdir / "summary.json").read_text()).get("metrics", {})
        for n in ("encoders", "dob20", "encoders_ff", "dob20_ff"):
            if n in cm and (cdir / "runs" / f"{n}.npz").exists():
                dc = dict(np.load(cdir / "runs" / f"{n}.npz"))
                cont[n] = dt.metrics(dc, circle, t_settle, bands, kt)
    js["continuous_metrics"] = cont
    js_path.write_text(json.dumps(js, indent=1))

    # ---- tables
    lab = {c["name"]: c["label"] for c in cases}
    md = ["# Tracking of the circle with EtherCAT at 1 ms and the ADIS16507 gyro", "",
          f"Configuration `{cc_path.name}` on `{cc['base']}`. EtherCAT cycle {1e3 * float(cc['ethercat']['cycle_s']):g} ms, "
          f"targets {cc['ethercat']['interpolation']}, feedforward current applied {cc['ethercat']['ff_delay_cycles']} cycle "
          f"later. Gyro {g.name}: {g.bandwidth_hz:g} Hz second order, {g.odr_hz:g} SPS, latency {1e3 * g.latency_s:g} ms, "
          f"ARW {cc['gyro']['arw_deg_rthr']} deg/sqrt(h), bias instability {cc['gyro']['bias_instability_deg_hr']} deg/h "
          f"(output noise of the model {js['gyro']['output_noise_deg_s']:.3f} deg/s rms; data sheet 0.068). Metrics on the "
          f"line of sight from {t_settle:g} s, one orbit; jitter: rms above each frequency.", ""]
    head = ["case", "E rms [mdeg]", "E max [mdeg]"] + [f"jitter > {b:g} Hz [urad]" for b in bands] + \
        ["FF current rms az / el [A]", "motor current rms az / el [A]"]
    rows = []
    for n in lab:
        x_ = M[n]
        rows.append([lab[n], f"{x_['E_rms_mdeg']:.2f}", f"{x_['E_max_mdeg']:.2f}"]
                    + [f"{x_[f'jitter_{b:g}Hz_urad']:.3g}" for b in bands]
                    + [f"{x_['i_ff_applied_rms_az_A']:.2f} / {x_['i_ff_applied_rms_el_A']:.2f}",
                       f"{x_['i_rms_az_A']:.2f} / {x_['i_rms_el_A']:.2f}"])
    md += ["## EtherCAT and gyro", "", out.markdown_table(head, rows), ""]
    if cont:
        rows = []
        for n, x_ in cont.items():
            y_ = M[n]
            rows.append([n, f"{x_['E_rms_mdeg']:.2f} -> {y_['E_rms_mdeg']:.2f}", f"{x_['E_max_mdeg']:.2f} -> {y_['E_max_mdeg']:.2f}"]
                        + [f"{x_[f'jitter_{b:g}Hz_urad']:.3g} -> {y_[f'jitter_{b:g}Hz_urad']:.3g}" for b in bands])
        md += ["## Continuous loop with an ideal gyro -> EtherCAT and ADIS16507", "",
               out.markdown_table(["case", "E rms [mdeg]", "E max [mdeg]"] + [f"jitter > {b:g} Hz [urad]" for b in bands], rows), ""]
    rows = [[n, r["message"], f"{r['cpu_s']:.0f}"] for n, r in js.get("runs", {}).items()]
    md += ["## Runs", "", out.markdown_table(["case", "result", "cpu [s]"], rows), ""]
    (res_dir / "summary.md").write_text("\n".join(md))

    # --------------------------------------------------------------------------- figures
    # gyro model
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    a = axs[0]
    a.loglog(ta, ad / DEG * 3600, "k.", ms=4, label=f"simulated noise ({cc['allan_check_s']:g} s at {g.odr_hz:g} SPS)")
    for k_, v in mdl.items():
        a.loglog(ta, v / DEG * 3600, "-" if k_ == "total" else "--", lw=1.2 if k_ == "total" else 0.8,
                 label=f"model: {k_}")
    a.loglog([1.0], [cc["gyro"]["arw_deg_rthr"] * 60], "rs", label="data sheet: ARW at 1 s")
    a.loglog([g.bias_tau_s], [cc["gyro"]["bias_instability_deg_hr"]], "r^", label="data sheet: bias instability (time assumed)")
    a.set_xlabel("averaging time [s]")
    a.set_ylabel("Allan deviation [deg/h]")
    a.set_title(f"{g.name}: Allan deviation of the rate noise")
    a.legend(loc="lower left")
    a = axs[1]
    f = np.logspace(0, np.log10(g.odr_hz / 2), 300)
    a.semilogx(f, 20 * np.log10(np.abs(g.response(f))), "C0", label="magnitude [dB]")
    a2 = a.twinx()
    a2.semilogx(f, np.angle(g.response(f), deg=True) - 360 * f * g.latency_s, "C3", label="phase with the latency [deg]")
    a.set_xlabel("frequency [Hz]")
    a.set_ylabel("magnitude [dB]")
    a2.set_ylabel("phase [deg]")
    a.set_title("gyro dynamics: second order + latency")
    h1, l1 = a.get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, loc="lower left")
    fig.tight_layout()
    out.save(fig, fig_dir, "gyro_model")
    # pointing error of the pairs, with the continuous ideal case
    pairs = cc["compare"]
    fig, axs = plt.subplots(len(pairs), 1, figsize=(10, 3.0 * len(pairs)), squeeze=False)
    for a, (n0, n1) in zip(axs[:, 0], pairs):
        for n, c in ((n0, "k"), (n1, "C3")):
            d = data[n]
            E = dt.pointing_error(circle, d["t"], d["psi_los"], d["th_los"])
            a.plot(d["t"], MDEG * E, color=c, lw=0.6, label=f"{lab[n]}, EtherCAT + ADIS16507 (rms {M[n]['E_rms_mdeg']:.2f} mdeg)")
        if n1 in cont:
            dc = dict(np.load(cdir / "runs" / f"{n1}.npz"))
            E = dt.pointing_error(circle, dc["t"], dc["psi_los"], dc["th_los"])
            a.plot(dc["t"], MDEG * E, color="C0", lw=0.8, ls="--", label=f"{lab[n1]}, continuous, ideal gyro (rms {cont[n1]['E_rms_mdeg']:.2f} mdeg)")
        a.axvline(t_settle, color="0.5", ls=":", lw=0.8, label="start of the metrics")
        a.set_yscale("log")
        a.set_ylim(1e-3, None)
        a.set_xlabel("time [s]")
        a.set_ylabel("pointing error [mdeg]")
        a.legend(loc="lower right")
    fig.tight_layout()
    out.save(fig, fig_dir, "pointing_error")
    # components, jitter and spectra of the main pair with velocity feedforward
    n0, n1 = pairs[-1]
    n2 = f"{n1}_nonoise" if f"{n1}_nonoise" in data else None
    fig, axs = plt.subplots(3, 2, figsize=(11, 8))
    for n, c in [(n0, "k"), (n1, "C3")] + ([(n2, "C0")] if n2 else []):
        d = data[n]
        t = d["t"]
        fs = 1.0 / float(np.median(np.diff(t)))
        m = t >= t_settle
        comps = (("cross-elevation", (d["psi_ref"] - d["psi_los"]) * np.cos(d["th_ref"])), ("elevation", d["th_ref"] - d["th_los"]))
        for j, (nm, e) in enumerate(comps):
            axs[0, j].plot(t, MDEG * e, color=c, lw=0.6, label=lab[n])
            b, a_ = butter(2, bands[0] / (fs / 2), "high")
            axs[1, j].plot(t[m], 1e6 * filtfilt(b, a_, e)[m], color=c, lw=0.5, label=lab[n])
            fw, P = welch(e[m] - np.mean(e[m]), fs=fs, nperseg=min(8192, int(m.sum())))
            axs[2, j].loglog(fw[1:], 1e6 * np.sqrt(P[1:]), color=c, lw=0.8, label=lab[n])
            axs[0, j].set_title(f"{nm}: line-of-sight error")
            axs[1, j].set_title(f"{nm}: error above {bands[0]:g} Hz (jitter)")
            axs[2, j].set_title(f"{nm}: amplitude spectral density of the error")
    for j in range(2):
        axs[0, j].set_ylabel("error [mdeg]")
        axs[1, j].set_ylabel("high-passed error [urad]")
        axs[2, j].set_ylabel("ASD [urad/sqrt(Hz)]")
        axs[2, j].set_xlabel("frequency [Hz]")
        for a in axs[:2, j]:
            a.set_xlabel("time [s]")
        for a in axs[:, j]:
            a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, f"components_{n1}")
    # sweep of the observer bandwidth: with and without gyro noise
    sw = sorted({c["observer_hz"] for c in cases if c.get("observer_hz") and c.get("feedforward") and "zoh" not in c["name"]})
    keys = [("E_rms_mdeg", "pointing error rms [mdeg]")] + [(f"jitter_{b:g}Hz_urad", f"jitter above {b:g} Hz [urad]") for b in bands]
    fig, axs = plt.subplots(1, len(keys), figsize=(3.4 * len(keys), 3.4))
    for a, (k_, yl) in zip(axs, keys):
        a.plot(sw, [M[f"dob{f:g}_ff"][k_] for f in sw], "o-", color="C3", label="gyro with noise (ADIS16507)")
        a.plot(sw, [M[f"dob{f:g}_ff_nonoise"][k_] for f in sw], "s--", color="C0", label="gyro without noise")
        a.axhline(M["encoders_ff"][k_], color="k", ls="--", lw=0.9, label="current structure (+ velocity FF)")
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xticks(sw, [f"{f:g}" for f in sw])
        a.minorticks_off()
        a.set_xlabel("observer bandwidth [Hz]")
        a.set_ylabel(yl)
        a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "observer_sweep")
    # feedforward current: noise-equivalent current
    fig, axs = plt.subplots(2, 2, figsize=(11, 6))
    for j, ax in enumerate(AXES):
        s = ax[:2]
        for n, c in [(n1, "C3")] + ([(n2, "C0")] if n2 else []):
            d = data[n]
            t = d["t"]
            m = t >= t_settle
            axs[0, j].plot(t, d[f"i_ff_{s}"], color=c, lw=0.4, label=f"feedforward current, {lab[n]}")
            fw, P = welch(d[f"i_ff_{s}"][m] - np.mean(d[f"i_ff_{s}"][m]), fs=1.0 / float(np.median(np.diff(t))), nperseg=4096)
            axs[1, j].loglog(fw[1:], np.sqrt(P[1:]), color=c, lw=0.8, label=lab[n])
        axs[0, j].set_title(f"{ax}: observer feedforward current")
        axs[0, j].set_xlabel("time [s]")
        axs[0, j].set_ylabel("current [A]")
        axs[1, j].set_title(f"{ax}: its amplitude spectral density")
        axs[1, j].set_xlabel("frequency [Hz]")
        axs[1, j].set_ylabel("ASD [A/sqrt(Hz)]")
        for a in axs[:, j]:
            a.legend(loc="best")
    fig.tight_layout()
    out.save(fig, fig_dir, "feedforward_current")
    print(f"done in {(time.time() - t_start) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
