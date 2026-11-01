"""Builds the report on the gyro disturbance observer (Hilkert) for the tracking of the circle of Attorrese.

    ../../.venv/bin/python -I build_dob_report.py [--no-latex]

Inputs: the summary.json and the saved runs of dob_tracking_study.py (continuous loop, ideal gyro) and of
dob_ethercat_study.py (EtherCAT at 1 ms, ADIS16507 gyro), and their figures (copied into figures/). The text is
in dob_report.tex; this script writes the tables and numbers (generated/), the figures of the jitter analysis
(figures/) and compiles with latexmk.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.signal import welch

HERE = Path(__file__).resolve().parent
CL = HERE.parent
sys.path.insert(0, str(CL))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import dob_tracking_study as dt  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
AXES = ("azimuth", "elevation")
CONT, ECAT = CL / "results_dob_tracking", CL / "results_dob_ethercat"
COPY = {
    "reference": CONT / "figures" / "reference.pdf",
    "cont_pointing_error": CONT / "figures" / "pointing_error.pdf",
    "cont_components_ff": CONT / "figures" / "components_dob20_ff.pdf",
    "cont_disturbance_estimate": CONT / "figures" / "disturbance_estimate.pdf",
    "cont_disturbance_rejection": CONT / "figures" / "disturbance_rejection.pdf",
    "cont_observer_sweep": CONT / "figures" / "observer_sweep.pdf",
    "ecat_gyro_model": ECAT / "figures" / "gyro_model.pdf",
    "ecat_pointing_error": ECAT / "figures" / "pointing_error.pdf",
    "ecat_components_ff": ECAT / "figures" / "components_dob20_ff.pdf",
    "ecat_observer_sweep": ECAT / "figures" / "observer_sweep.pdf",
    "ecat_feedforward_current": ECAT / "figures" / "feedforward_current.pdf",
}
URAD_DEG = 180.0 / np.pi * 1e-6   # degrees per microradian
plt.rcParams.update({"font.size": 8.5, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 7,
                     "savefig.bbox": "tight"})


def rows_tex(rows: list[list[str]]) -> str:
    return " \\\\\n".join(" & ".join(r) for r in rows)


def write(name: str, text: str) -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def sci(x: float, digits: int = 2) -> str:
    """Number in LaTeX scientific notation."""
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    if -2 <= e <= 2:
        return f"{x:.{max(0, digits - 1 - e)}f}"
    return f"${x / 10**e:.{digits - 1}f}\\cdot10^{{{e}}}$"


def band(m: dict, lo: float, hi: float | None) -> float:
    """Rms of the line-of-sight error between two jitter bands (from the nested values above each frequency)."""
    a = m[f"jitter_{lo:g}Hz_urad"]
    b = m[f"jitter_{hi:g}Hz_urad"] if hi else 0.0
    return float(np.sqrt(max(a * a - b * b, 0.0)))


MACRO_AX = {"azimuth": "Az", "elevation": "El"}
T_CONT, T_PEAK = 2.18, 6.74   # continuous and peak stall torque of the TBM(S)-7646-A, N m (motor data sheet)


def torque_profile(cc_c: dict, fin: dict, M: dict) -> None:
    """Motor torque kappa_i K_t i of the current structure and of the observer (ideal and ADIS16507 gyro),
    EtherCAT runs with the velocity feedforward: figure, table and numbers."""
    from scipy.signal import butter, filtfilt

    cases = [("encoders_ff", "current (encoders)", "k"), ("dob20_ff_nonoise", "observer 20 Hz, ideal gyro", "C0"),
             ("dob20_ff", "observer 20 Hz, ADIS16507", "C3")]
    t_settle = float(cc_c["simulation"]["t_settle_s"])
    kt = 0.258
    r = cc_c["reference"]
    circle = dt.Circle(tuple(r["centre_m"]), float(r["radius_m"]), float(r["height_m"]), float(r["omega_rad_s"]))
    data = {}
    for n, _, _ in cases:
        d = dict(np.load(ECAT / "runs" / f"{n}.npz"))
        data[n] = {k: v[:-1] for k, v in d.items()}
    t = data["encoders_ff"]["t"]
    m = t >= t_settle
    fs = 1.0 / float(np.median(np.diff(t)))
    b_lp, a_lp = butter(2, 10.0 / (fs / 2))
    # inertial torque of the reference, J times the angular acceleration of the pointing angles
    tt = np.linspace(t_settle, t[-1], 20001)
    psi, th, dpsi, dth = circle.angles(tt)
    acc = {"azimuth": np.gradient(dpsi, tt), "elevation": np.gradient(dth, tt)}
    J = {"azimuth": 2.162, "elevation": 0.900}
    rows = []
    rev = tt[np.nonzero(np.diff(np.sign(dpsi)))[0][0]]  # first velocity reversal of the azimuth
    rev_el = tt[np.nonzero(np.diff(np.sign(dth)))[0][0]]
    fig, axs = plt.subplots(3, 2, figsize=(11, 9))
    for j, ax in enumerate(AXES):
        a = ax[:2]
        kap = fin[ax]["scale"]
        for n, lab, col in cases:
            tau = kap * kt * data[n][f"i_{a}"]
            lo = filtfilt(b_lp, a_lp, tau)
            rms, pk = float(np.sqrt(np.mean(tau[m] ** 2))), float(np.max(np.abs(tau[m])))
            rms_lo, rms_hi = float(np.sqrt(np.mean(lo[m] ** 2))), float(np.std((tau - lo)[m]))
            rows.append([ax if n == "encoders_ff" else "", lab, f"{rms:.3f}", f"{pk:.2f}", f"{rms_lo:.3f}", f"{rms_hi:.3f}",
                         f"{100 * rms / T_CONT:.0f}", f"{100 * pk / T_PEAK:.0f}"])
            key = {"encoders_ff": "Cur", "dob20_ff_nonoise": "Ideal", "dob20_ff": "Dob"}[n]
            A = MACRO_AX[ax]
            M[f"tauRms{key}{A}"], M[f"tauPk{key}{A}"] = f"{rms:.2f}", f"{pk:.2f}"
            M[f"tauLo{key}{A}"], M[f"tauHi{key}{A}"] = f"{rms_lo:.2f}", f"{rms_hi:.2f}"
            lw = 0.4 if n == "dob20_ff" else 0.9
            axs[0, j].plot(t[m], tau[m], color=col, lw=lw, label=lab)
            tr = rev if ax == "azimuth" else rev_el
            w = (t >= tr - 1.5) & (t <= tr + 2.5)
            axs[1, j].plot(t[w], tau[w], color=col, lw=lw, label=lab)
            fw, P = welch(tau[m] - tau[m].mean(), fs=fs, nperseg=4096)
            axs[2, j].loglog(fw[1:], np.sqrt(P[1:]), color=col, lw=0.8, label=lab)
        for n, key in (("encoders_ff", "Cur"), ("dob20_ff_nonoise", "Ideal"), ("dob20_ff", "Dob")):
            M[f"fric{key}{MACRO_AX[ax]}"] = f"{np.mean(np.abs(data[n][f'T_f_{a}'][m])):.2f}"
        M[f"cable{MACRO_AX[ax]}"] = f"{np.mean(data['encoders_ff'][f'T_cable_{a}'][m]):+.2f}"
        ia = np.max(np.abs(J[ax] * acc[ax]))
        M[f"tauInertia{MACRO_AX[ax]}"] = f"{1e3 * ia:.1f}"
        axs[0, j].set_title(f"{ax}: motor torque over one orbit")
        axs[1, j].set_title(f"{ax}: around the first velocity reversal ({(rev if ax == 'azimuth' else rev_el):.1f} s)")
        axs[2, j].set_title(f"{ax}: amplitude spectral density of the torque")
        for a_ in axs[:2, j]:
            a_.set_xlabel("time [s]")
            a_.set_ylabel("motor torque [N m]")
        axs[2, j].set_xlabel("frequency [Hz]")
        axs[2, j].set_ylabel("ASD [N m/sqrt(Hz)]")
        for a_ in axs[:, j]:
            a_.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIG / "torque_profile.pdf")
    plt.close(fig)
    write("tab_torque.tex", rows_tex(rows))
    M["tCont"], M["tPeak"] = f"{T_CONT:.2f}", f"{T_PEAK:.2f}"
    M["revAz"], M["revEl"] = f"{rev:.1f}", f"{rev_el:.1f}"


def rejection_analysis(cc_c: dict, cj: dict, fin: dict, M: dict) -> None:
    """Disturbance rejection, bearing torque to line of sight: numerical transfer functions of the linearised
    loops, rigid-body model with the transfer functions of the drive, asymptotes and compliance floor."""
    f = np.logspace(-1, np.log10(500), 500)
    num = dt.rejection(cc_c, f)
    loop, _, _ = dt.build_loop(cc_c, None, False)
    pl = loop.plant
    _, _, K = pl.linear_matrices()
    n_psi = pl.sim.struttura.n_psi
    s = 2j * np.pi * f
    bws = sorted({c["observer_hz"] for c in cc_c["cases"] if c.get("observer_hz")})
    rows = []
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
    for j, ax in enumerate(AXES):
        A = MACRO_AX[ax]
        d = loop.drives[ax]
        J = cj["observer"][ax]["J"]
        kt = loop.motors[ax].p.Kt * fin[ax]["scale"]
        # static compliance between the encoder and the line of sight for a bearing torque balanced by the motor
        sl = list(range(0, n_psi)) if ax == "azimuth" else list(range(n_psi, pl.n))
        ie, im, ib = loop.ie[ax], loop.im[ax], pl.index(ax, "bearing")
        il = loop.jl if ax == "azimuth" else loop.jp
        free = [k for k in sl if k != ie]
        F = np.zeros(pl.n)
        F[ib] += 1.0
        F[im] -= 1.0
        q = np.zeros(pl.n)
        q[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        floor = float(q[il] - q[ie])
        # rigid-body loop with the transfer functions of the drive (current loop ideal)
        Cl = d.F_fwd(s) * d.C_vel(s) * (d.C_pos(s) + s * d.F_fb(s))
        G_cur = 1.0 / (J * s**2 + kt * Cl)
        p = d.p
        kI = p.Kpv * p.Kpp * 2 * np.pi * p.Kiv_hz
        M[f"floor{A}"] = f"{1e6 * abs(floor):.1f}"
        M[f"kI{A}"] = f"{kI:.0f}"
        M[f"ktk{A}"] = f"{kt:.3f}"
        M[f"lowCur{A}"] = f"{1e6 * 2 * np.pi * 0.1 / (kt * kI):.0f}"
        M[f"highHund{A}"] = f"{1e6 / (J * (2 * np.pi * 100) ** 2):.2f}"
        a = axs[j]
        a.loglog(f, 1e6 * np.abs(num[(ax, None)]), color="k", lw=1.6, label="current structure (linearised H6 loop)")
        a.loglog(f, 1e6 * np.abs(G_cur + floor), color="k", ls="--", lw=0.9, label="current, rigid model with the drive filters")
        for k, bw in enumerate(bws):
            col = f"C{k}"
            a.loglog(f, 1e6 * np.abs(num[(ax, bw)]), color=col, lw=1.1, label=f"observer {bw:g} Hz (linearised H6 loop)")
            if bw == 20.0:
                wo = 2 * np.pi * bw
                a.loglog(f, 1e6 * np.abs(G_cur * s / (s + wo) + floor), color=col, ls="--", lw=0.9,
                         label="observer 20 Hz, rigid model x s/(s+w_o)")
        w = 2 * np.pi * f
        a.loglog(f, 1e6 * w / (kt * kI), ":", color="0.4", lw=0.9, label="low f, current: w/(kappa K_t k_I)")
        a.loglog(f, 1e6 * w**2 / (kt * kI * 2 * np.pi * 20), ":", color="C2", lw=0.9,
                 label="low f, observer 20 Hz: w^2/(kappa K_t k_I w_o)")
        a.loglog(f, 1e6 / (J * w**2), ":", color="C4", lw=0.9, label="high f: 1/(J w^2)")
        if abs(floor) * 1e6 > 0.5:
            a.axhline(1e6 * abs(floor), color="C5", ls="-.", lw=0.9,
                      label=f"compliance encoder-LOS: {1e6 * abs(floor):.0f} urad/(N m)")
        a.set_ylim(1e-2, 2e4)
        a.set_xlabel("frequency [Hz]")
        a.set_ylabel("LOS angle per bearing torque [urad/(N m)]")
        a.set_title(f"{ax}: disturbance rejection, bearing torque to line of sight")
        a.legend(loc="lower left", fontsize=5.8)
        # table: magnitude at selected frequencies, ratio to the current structure and Hilkert's |s/(s+w_o)|
        for fq in (0.1, 1.0, 3.0, 10.0, 30.0):
            k_ = int(np.argmin(np.abs(f - fq)))
            g0 = abs(num[(ax, None)][k_])
            r_ = [ax if fq == 0.1 else "", f"{fq:g}", f"{1e6 * g0:.3g}"]
            for bw in bws:
                g = abs(num[(ax, bw)][k_])
                r_.append(f"{1e6 * g:.3g} ({g / g0:.2g}; {fq / np.hypot(fq, bw):.2g})")
            rows.append(r_)
        k3 = int(np.argmax(np.abs(num[(ax, None)])))
        M[f"peakF{A}"], M[f"peakG{A}"] = f"{f[k3]:.1f}", f"{1e6 * abs(num[(ax, None)][k3]):.0f}"
        k1 = int(np.argmin(np.abs(f - 1.0)))
        M[f"oneCur{A}"] = f"{1e6 * abs(num[(ax, None)][k1]):.0f}"
        M[f"oneDob{A}"] = f"{1e6 * abs(num[(ax, 20.0)][k1]):.0f}"
        rows.append(["\\midrule"])
    fig.tight_layout()
    fig.savefig(FIG / "rejection_annotated.pdf")
    plt.close(fig)
    write("tab_rejection.tex", rows_tex(rows[:-1]).replace("\\midrule \\\\", "\\midrule"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for dst, src in COPY.items():
        shutil.copyfile(src, FIG / f"{dst}.pdf")
    cj = json.loads((CONT / "summary.json").read_text())
    ej = json.loads((ECAT / "summary.json").read_text())
    cc_c = yaml.safe_load((CL / "config_dob_tracking.yaml").read_text())
    cc_e = yaml.safe_load((CL / "config_dob_ethercat.yaml").read_text())
    Mc, Me = cj["metrics"], ej["metrics"]
    M: dict[str, str] = {}

    # ---- model parameters
    gains = json.loads((CL / "results_vcl_study" / "summary.json").read_text())["iterations"][-1]["values"]
    fin = json.loads((CL / "results_v2_identification" / "summary.json").read_text())["final"]
    rows = []
    for ax in AXES:
        g, f = gains[ax], fin[ax]
        rows.append([ax, f"{g['Kpp']:.2f}", f"{g['Kpv']:.0f}", f"{g['Kiv']:.2f}", f"{f['scale']:.2f}", f"{f['Tc']:.2f}",
                     f"{f['Tc'] * f['Ts_ratio']:.2f}", f"{np.rad2deg(f['vs']):.2f}", f"{f['sigma0']:.0f}",
                     f"{f['K_cable']:.2f}", f"{f['theta0_cable_deg']:+.2f}", f"{cj['observer'][ax]['J']:.3f}"])
    write("tab_model.tex", rows_tex(rows))

    # ---- reference
    r = cc_c["reference"]
    circle = dt.Circle(tuple(r["centre_m"]), float(r["radius_m"]), float(r["height_m"]), float(r["omega_rad_s"]))
    tt = np.linspace(0, 2 * np.pi / circle.omega, 4000)
    psi, th, dpsi, dth = circle.angles(tt)
    M.update(psiMax=f"{np.rad2deg(np.abs(psi).max()):.2f}", thMin=f"{np.rad2deg(th.min()):.2f}", thMax=f"{np.rad2deg(th.max()):.2f}",
             dpsiMax=f"{np.rad2deg(np.abs(dpsi).max()):.2f}", dthMax=f"{np.rad2deg(np.abs(dth).max()):.3f}",
             period=f"{2 * np.pi / circle.omega:.1f}", followAz=f"{1e3 * np.rad2deg(np.abs(dpsi).max()) / gains['azimuth']['Kpp']:.0f}")

    # ---- stability
    rows = []
    st = cj["stability"]
    for ax, notch in (("azimuth", False), ("elevation", True), ("elevation", False)):
        rr = [x for x in st if x["axis"] == ax and x["notch"] == notch]
        if not rr:
            continue
        rows.append([ax, "yes" if notch else "no"] + ["stable" if x["stable"] else
                                                      f"\\textbf{{{x['unstable_hz'][0]:.0f}\\,Hz}}" for x in rr])
    bws = [x["bandwidth_hz"] for x in st if x["axis"] == "azimuth"]
    write("tab_stability.tex", rows_tex(rows))
    write("tab_stability_head.tex", " & ".join(["axis", "notch"] + [f"{b:g}\\,Hz" for b in bws]))
    write("tab_stability_cols.tex", "ll" + "c" * len(bws))

    # ---- results tables (continuous and EtherCAT)
    lab_c = {c["name"]: c["label"] for c in cc_c["cases"]}
    rows = []
    for n, lab in lab_c.items():
        x = Mc[n]
        rows.append([lab.replace("%", "\\%").replace("_", "\\_"), f"{x['E_rms_mdeg']:.2f}", f"{x['E_max_mdeg']:.1f}", f"{x['jitter_1Hz_urad']:.3g}",
                     f"{x['jitter_10Hz_urad']:.3g}", f"{x['xel_rms_mdeg']:.2f}", f"{x['el_rms_mdeg']:.2f}"])
    write("tab_cont.tex", rows_tex(rows))
    lab_e = {c["name"]: c["label"] for c in cc_e["cases"]}
    rows = []
    for n, lab in lab_e.items():
        x = Me[n]
        rows.append([lab.replace("_", "\\_"), f"{x['E_rms_mdeg']:.2f}", f"{x['E_max_mdeg']:.1f}", f"{x['jitter_1Hz_urad']:.3g}",
                     f"{x['jitter_10Hz_urad']:.3g}", f"{x['jitter_100Hz_urad']:.3g}",
                     f"{x['i_ff_applied_rms_az_A']:.2f} / {x['i_ff_applied_rms_el_A']:.2f}"])
    write("tab_ecat.tex", rows_tex(rows))

    # ---- focus: current against the observer with the noisy gyro (EtherCAT)
    cases = [("encoders_ff", "current (encoders)"), ("dob20_ff_nonoise", "observer 20\\,Hz, ideal gyro"),
             ("dob20_ff", "observer 20\\,Hz, ADIS16507")]
    rows_m, rows_d = [], []
    for n, lab in cases:
        x = Me[n]
        rows_m.append([lab, f"{x['E_rms_mdeg']:.2f}", f"{x['E_p95_mdeg']:.2f}", f"{x['E_max_mdeg']:.1f}",
                       f"{x['miss_rms_mm']:.1f}", sci(x["jitter_1Hz_urad"], 3), sci(band(x, 1, 10), 3),
                       sci(band(x, 10, 100), 3), sci(x["jitter_100Hz_urad"], 3)])
        rows_d.append([lab, sci(x["E_rms_mdeg"] / 1e3), sci(x["E_max_mdeg"] / 1e3), sci(x["jitter_1Hz_urad"] * URAD_DEG),
                       sci(band(x, 1, 10) * URAD_DEG), sci(band(x, 10, 100) * URAD_DEG), sci(x["jitter_100Hz_urad"] * URAD_DEG)])
    write("tab_focus.tex", rows_tex(rows_m))
    write("tab_focus_deg.tex", rows_tex(rows_d))
    b0, bi, bn = Me["encoders_ff"], Me["dob20_ff_nonoise"], Me["dob20_ff"]
    M.update(
        ErmsCur=f"{b0['E_rms_mdeg']:.1f}", ErmsDob=f"{bn['E_rms_mdeg']:.2f}", ErmsIdeal=f"{bi['E_rms_mdeg']:.2f}",
        EmaxCur=f"{b0['E_max_mdeg']:.0f}", EmaxDob=f"{bn['E_max_mdeg']:.1f}",
        EratioRms=f"{b0['E_rms_mdeg'] / bn['E_rms_mdeg']:.0f}", EratioMax=f"{b0['E_max_mdeg'] / bn['E_max_mdeg']:.0f}",
        JoneCur=f"{b0['jitter_1Hz_urad']:.0f}", JoneDob=f"{bn['jitter_1Hz_urad']:.1f}", JoneIdeal=f"{bi['jitter_1Hz_urad']:.1f}",
        JoneRatio=f"{b0['jitter_1Hz_urad'] / bn['jitter_1Hz_urad']:.0f}",
        JtenCur=f"{b0['jitter_10Hz_urad']:.2f}", JtenDob=f"{bn['jitter_10Hz_urad']:.2f}", JtenIdeal=f"{bi['jitter_10Hz_urad']:.3f}",
        JhunCur=sci(b0["jitter_100Hz_urad"]), JhunDob=f"{bn['jitter_100Hz_urad']:.2f}",
        BandOneTenCur=f"{band(b0, 1, 10):.0f}", BandOneTenDob=f"{band(bn, 1, 10):.1f}",
        BandTenHunCur=f"{band(b0, 10, 100):.2f}", BandTenHunDob=f"{band(bn, 10, 100):.2f}",
        NoiseJone=f"{np.sqrt(max(bn['jitter_1Hz_urad'] ** 2 - bi['jitter_1Hz_urad'] ** 2, 0)):.1f}",
        IffNoise=f"{np.sqrt(max(bn['i_ff_applied_rms_az_A'] ** 2 - bi['i_ff_applied_rms_az_A'] ** 2, 0)):.1f}",
        IffAz=f"{bn['i_ff_applied_rms_az_A']:.2f}", IffAzIdeal=f"{bi['i_ff_applied_rms_az_A']:.2f}",
        ErmsNoFFCur=f"{Me['encoders']['E_rms_mdeg']:.0f}", ErmsNoFFDob=f"{Me['dob20']['E_rms_mdeg']:.0f}",
        ErmsZoh=f"{Me['dob20_ff_zoh']['E_rms_mdeg']:.2f}", ErmsCurZoh=f"{Me['encoders_ff_zoh']['E_rms_mdeg']:.1f}",
        gyroNoise=f"{ej['gyro']['output_noise_deg_s']:.3f}",
        twist=f"{Mc['dob20_ff']['los_enc_az_rms_mdeg']:.2f}", encAzDob=f"{Mc['dob20_ff']['enc_az_rms_mdeg']:.2f}")
    for f in (5, 10, 20, 40):
        M[f"Jsweep{['Five', 'Ten', 'Twenty', 'Forty'][[5, 10, 20, 40].index(f)]}"] = f"{Me[f'dob{f}_ff']['jitter_1Hz_urad']:.1f}"
    # noise floor of the line of sight following the gyro between the position-loop and the observer bandwidths
    N = cc_e["gyro"]["arw_deg_rthr"] * np.pi / 180 / 60
    M["noiseFloor"] = f"{1e6 * N / (2 * np.pi) * np.sqrt(1 / 1.0 - 1 / 20.0):.1f}"
    M["noiseFloorTen"] = f"{1e6 * N / (2 * np.pi) * np.sqrt(1 / 10.0 - 1 / 20.0):.1f}"
    M["arwSI"] = f"{1e5 * N:.2f}"

    # ---- figure: jitter by band
    fig, ax = plt.subplots(figsize=(8, 3.6))
    bands_ = [("1--10 Hz", lambda x: band(x, 1, 10)), ("10--100 Hz", lambda x: band(x, 10, 100)),
              ("above 100 Hz", lambda x: x["jitter_100Hz_urad"]), ("total above 1 Hz", lambda x: x["jitter_1Hz_urad"])]
    w = 0.26
    for k, ((n, lab), col) in enumerate(zip(cases, ("k", "C0", "C3"))):
        v = [fn(Me[n]) for _, fn in bands_]
        ax.bar(np.arange(len(bands_)) + (k - 1) * w, v, w, color=col, alpha=0.8, label=lab.replace("\\,", " "))
    ax.set_yscale("log")
    ax.set_xticks(range(len(bands_)), [b for b, _ in bands_])
    ax.set_ylabel("rms line-of-sight error in the band [urad]")
    ax.set_title("jitter by frequency band, EtherCAT 1 ms, velocity feedforward, one orbit")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "jitter_bands.pdf")
    plt.close(fig)

    # ---- figure: cumulative jitter (rms of the error above each frequency) from the spectra
    t_settle = float(cc_c["simulation"]["t_settle_s"])
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for (n, lab), col in zip(cases, ("k", "C0", "C3")):
        d = dict(np.load(ECAT / "runs" / f"{n}.npz"))
        d = {k: v[:-1] for k, v in d.items()}
        t = d["t"]
        m = t >= t_settle
        fs = 1.0 / float(np.median(np.diff(t)))
        for j, (nm, e) in enumerate((("cross-elevation", (d["psi_ref"] - d["psi_los"]) * np.cos(d["th_ref"])),
                                     ("elevation", d["th_ref"] - d["th_los"]))):
            fw, P = welch(e[m] - e[m].mean(), fs=fs, nperseg=8192)
            df = fw[1] - fw[0]
            cum = np.sqrt(np.cumsum((P * df)[::-1])[::-1])
            axs[j].loglog(fw[1:], 1e6 * cum[1:], color=col, lw=1.0, label=lab.replace("\\,", " "))
            axs[j].set_title(f"{nm}: rms of the error above each frequency")
    for a in axs:
        a.set_xlabel("frequency [Hz]")
        a.set_ylabel("cumulative rms [urad]")
        a.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "cumulative_jitter.pdf")
    plt.close(fig)

    torque_profile(cc_c, fin, M)
    rejection_analysis(cc_c, cj, fin, M)
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("tables, numbers and figures written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "dob_report.tex"], cwd=HERE,
                       check=True)


if __name__ == "__main__":
    main()
