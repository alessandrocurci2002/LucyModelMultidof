"""Builds the report on the Gawronski LQG controller (PI + flexible mode controller) on the LUCY doublets.

    ../.venv/bin/python -I build_lqg_report.py [--no-latex]

Inputs (no fit and no simulation here): ``results_lqg_gawronski`` (``analisi_lqg_gawronski.py``) and
``results_fork_modes_dahl_resampled`` (structure and effective gains). The text is in lqg_report.tex; this script copies
the figures (figures/), writes tables and numbers (generated/) and compiles with latexmk.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "closedLoopModel"), str(ROOT)]

import analisi_anello_chiuso_ibrido as aac  # noqa: E402
from lucy_cl.driver import DriveParams  # noqa: E402

FIG, GEN = HERE / "figures", HERE / "generated"
RES = ROOT / "results_lqg_gawronski"
RIC = ROOT / "results_fork_modes_dahl_resampled"
COPY = {"freq_elevation": RIC / "figures" / "elevation_bode_fit.pdf", "freq_azimuth": RIC / "figures" / "azimuth_bode_fit.pdf",
        "vl_azimuth": RES / "figures" / "azimuth_velocity_loop.pdf", "vl_elevation": RES / "figures" / "elevation_velocity_loop.pdf",
        "scan": RES / "figures" / "tuning_scan.pdf", "tf": RES / "figures" / "closed_loop_tf.pdf",
        "loop": RES / "figures" / "position_loop.pdf", "db_azimuth": RES / "figures" / "azimuth_doublets.pdf",
        "db_elevation": RES / "figures" / "elevation_doublets.pdf"}
TESTS = {"01deg_minus01deg": "0.1", "05deg_minus05deg": "0.5", "1deg_minus1deg": "1", "3deg_minus3deg": "3", "5deg_minus5deg": "5", "8deg_minus8deg": "8"}
SCHEMES = ("measured", "cascade", "lqg")
AX = (("azimuth", "azimut", "A"), ("elevation", "polare", "P"))


def rows_tex(rows):
    return " \\\\\n".join(" & ".join(str(x) for x in r) for r in rows)


def write(name, text):
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def sci(x, d=2):
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    if -2 <= e <= 3:
        s = f"{x:.{max(d - e, 0)}f}"
        return s.replace("-", "$-$", 1)
    m = x / 10 ** e
    return rf"\ensuremath{{{m:.{d}f}\times10^{{{e}}}}}"


def mod_label(f, om):
    return "structural" if np.min(np.abs(om - f)) < 3.0 else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for f in FIG.glob("*"):
        f.unlink()
    for n, src in COPY.items():
        shutil.copyfile(src, FIG / f"{n}.pdf")
    Sm = json.loads((RES / "summary.json").read_text())
    Tu = json.loads((RES / "tuning.json").read_text())
    Me = json.loads((RES / "medians.json").read_text())
    Hu = json.loads((RES / "hunting.json").read_text())
    Ri = json.loads((RIC / "summary.json").read_text())
    Td = json.loads((RIC / "time_domain_summary.json").read_text())
    with open(RES / "doublet_metrics.csv", newline="") as fh:
        R = list(csv.DictReader(fh))
    M = {}

    # ---- structure and drive
    rows = []
    for ch, asse, X in AX:
        v = Ri[ch]["uniform_lnf"]["hybrid"]["parameters"]
        rows.append([ch, f"{v['sigma0']:.0f}", f"{v['Tc']:.3f}", f"{v['g_Kt']:.3f}", f"{v['tau_d'] * 1e3:.2f}",
                     " / ".join(f"{m['f_hz']:.1f}" for m in Sm[ch]["modes"] if len(m["indici"]) == 2 and m["zeta"] < 0.1),
                     " / ".join(f"{m['zeta']:.3f}" for m in Sm[ch]["modes"] if len(m["indici"]) == 2 and m["zeta"] < 0.1)])
    write("tab_structure.tex", rows_tex(rows))
    rows = []
    for ch, asse, X in AX:
        akd, wb = aac.DRIVE[asse]
        p = DriveParams.from_files(ch, aac.VCL / akd, aac.VCL / wb, pl_kp_unit="hz")
        g = Td["gains"][ch]["uniform_lnf"]
        fl = [b.describe() for b in p.forward if b.active]
        fb = [b.describe() for b in p.feedback if b.active] + ([f"low pass {p.observer_bw_hz:.0f} Hz (observer)"] if p.observer else [])
        rows.append([ch, f"{g['Kpp']:.2f}", f"{g['Kpv']:.1f}", f"{g['Kiv_hz']:.2f}", f"{p.Kpi / p.L_motor / (2 * np.pi):.0f}", f"{p.i_limit:.0f}",
                     "; ".join(fl).replace("%", r"\%") or "none", "; ".join(fb).replace("%", r"\%") or "none"])
        M[f"Kpp{X}"], M[f"Kpv{X}"], M[f"Kiv{X}"] = f"{g['Kpp']:.2f}", f"{g['Kpv']:.0f}", f"{g['Kiv_hz']:.2f}"
    write("tab_drive.tex", rows_tex(rows))

    # ---- modes of the velocity loop
    for ch, asse, X in AX:
        rows, tutti = [], Tu[ch]["modes_all"]
        scelti = {round(m["f_hz"], 4) for m in Sm[ch]["modes"]}
        for m in tutti:
            rig = not np.isfinite(m["hinf"])
            rows.append(["rigid" if rig else f"{m['f_hz']:.2f}", "--" if rig or len(m["indici"]) == 1 else f"{m['zeta']:.3f}",
                         "$\\infty$" if rig else sci(m["hinf"]), str(len(m["indici"])), "yes" if round(m["f_hz"], 4) in scelti else ""])
        write(f"tab_modes_{ch}.tex", rows_tex(rows))

    # ---- tuning
    rows = []
    for ch, asse, X in AX:
        s, t = Sm[ch], Tu[ch]["best"]
        rows.append([ch, sci(t["q_i"]), sci(t["q_p"]), sci(t["q_f"]) if t["q_f"] else "0", sci(t["rho"]), f"{s['k_p']:.2f}", f"{s['k_i']:.2f}",
                     f"{0.25 * s['k_p'] ** 2:.0f}", f"{s['Kpp_cascade']:.2f}", f"{np.linalg.norm(s['K_f']):.0f}"])
        M |= {f"kp{X}": f"{s['k_p']:.1f}", f"ki{X}": f"{s['k_i']:.1f}", f"qi{X}": sci(t["q_i"]), f"qp{X}": sci(t["q_p"]),
              f"nOk{X}": str(sum(r["ok"] for r in Tu[ch]["scan"])), f"nScan{X}": str(len(Tu[ch]["scan"])), f"nRed{X}": str(s["states_reduced"]), f"nFull{X}": str(s["states_full"]),
              f"tauI{X}": f"{s['k_p'] / s['k_i']:.1f}"}
    write("tab_weights.tex", rows_tex(rows))
    for ch, asse, X in AX:
        rows = []
        for r in Tu[ch]["scan_flex"]:
            sel = r["q_f"] == Tu[ch]["best"]["q_f"] and r["rho"] == Tu[ch]["best"]["rho"]
            b = (lambda x: rf"\textbf{{{x}}}") if sel else str
            rows.append([b(sci(r["q_f"]) if r["q_f"] else "0"), b(sci(r["rho"])), b(f"{r['K_f_flex']:.0f}"), b(" / ".join(f"{z:.3f}" for z in r["zeta"].values())), b(f"{r['mod_min']:.2f}"),
                         b(f"{r['indice']:.2f}"), b("yes" if r["ok"] else "no")])
        write(f"tab_flex_{ch}.tex", rows_tex(rows))
    rows = []
    for ch, asse, X in AX:
        for r in Hu[ch]:
            rows.append([ch if r is Hu[ch][0] else "", r"\textbf{final}" if r["final"] else "", f"{r['k_p']:.1f}", f"{r['k_i']:.1f}", f"{r['e_max_last_3s_mdeg']:.1f}", f"{100 * r['moving_fraction']:.0f}", str(r["sign_changes"]), f"{r['modulus_stuck']:.2f}"])
    write("tab_hunting.tex", rows_tex(rows))

    # ---- linear
    rows = []
    for ch, asse, X in AX:
        s = Sm[ch]
        for sch, lab in (("cascade", "cascade"), ("lqg", r"\textbf{LQG}")):
            m = s["margins"][sch]
            z = s["zeta_flex"][sch]
            rows.append([ch if sch == "cascade" else "", lab, f"{s['bandwidth_Hz'][sch]:.2f}", f"{s['peak_T_dB'][sch]:.2f}", f"{s['T_wy_at_1Hz_dB'][sch]:.1f}",
                         f"{m['sliding']['modulus']:.2f}", f"{m['stuck']['modulus']:.2f}", f"{m['measured']['modulus']:.2f}",
                         f"{min(x['gm_dB'] for x in m.values()):.1f}", f"{min(x['pm_deg'] for x in m.values()):.0f}", " / ".join(f"{x:.3f}" for x in z.values())])
        rob = s["robustness"]
        M[f"rhoMax{X}"] = f"{max(v['rho_max'] for v in rob.values()):.5f}"
        M[f"robMod{X}"] = f"{min(v['modulus'] for v in rob.values()):.2f}"
        M[f"iCmd{X}"] = f"{s['limits']['lqg']['i_cmd']:.1f}"
        M[f"iCmdC{X}"] = f"{s['limits']['cascade']['i_cmd']:.1f}"
        M[f"bwC{X}"], M[f"bwL{X}"] = f"{s['bandwidth_Hz']['cascade']:.2f}", f"{s['bandwidth_Hz']['lqg']:.2f}"
        M[f"zOne{X}"] = f"{list(s['zeta_flex']['lqg'].values())[0]:.3f}"
        M[f"zOneC{X}"] = f"{list(s['zeta_flex']['cascade'].values())[0]:.3f}"
    write("tab_linear.tex", rows_tex(rows))
    rows = []
    for ch, asse, X in AX:
        reg = sorted(Sm[ch]["regulator_poles_hz_zeta"])
        est = sorted(Sm[ch]["estimator_poles_hz_zeta"])
        rows.append([ch, "; ".join(f"{a:.2f}" + (f" ({b:.2f})" if b < 0.999 else "") for a, b in reg), "; ".join(f"{a:.2f}" + (f" ({b:.2f})" if b < 0.999 else "") for a, b in est)])
    write("tab_poles.tex", rows_tex(rows))

    # ---- doublets
    KEYS = ("err_max_mdeg", "err_rms_mdeg", "over_mdeg", "t_set_s", "i_max_A")
    FMT = {"err_max_mdeg": ".0f", "err_rms_mdeg": ".0f", "over_mdeg": ".0f", "t_set_s": ".1f", "i_max_A": ".1f"}
    for ch, asse, X in AX:
        rows = []
        for t, lab in TESTS.items():
            row = [lab]
            for k in KEYS:
                for sch in SCHEMES:
                    x = [float(r[k]) for r in R if r["axis"] == asse and r["test"] == t and r["scheme"] == sch and not (asse == "azimut" and t == "8deg_minus8deg" and r["move"] == "#1")]
                    row.append(format(np.median(x), FMT[k]))
            rows.append(row)
        write(f"tab_doublets_{ch}.tex", rows_tex(rows))
    rows = []
    for ch, asse, X in AX:
        for sch, lab in (("measured", "measured (cascade, real axis)"), ("cascade", "model, cascade"), ("lqg", r"\textbf{model, LQG}")):
            m = Me[asse][sch]
            rows.append([ch if sch == "measured" else "", lab] + [format(m[k], FMT[k]) for k in KEYS] + [f"{m['i_rms_A']:.2f}"])
        c, l = Me[asse]["cascade"], Me[asse]["lqg"]
        M |= {f"errC{X}": f"{c['err_max_mdeg']:.0f}", f"errL{X}": f"{l['err_max_mdeg']:.0f}", f"rmsC{X}": f"{c['err_rms_mdeg']:.0f}", f"rmsL{X}": f"{l['err_rms_mdeg']:.0f}",
              f"rmsD{X}": f"{100 * (1 - l['err_rms_mdeg'] / c['err_rms_mdeg']):.0f}", f"errD{X}": f"{100 * (1 - l['err_max_mdeg'] / c['err_max_mdeg']):.0f}",
              f"tC{X}": f"{c['t_set_s']:.1f}", f"tL{X}": f"{l['t_set_s']:.1f}", f"oL{X}": f"{l['over_mdeg']:.0f}", f"iC{X}": f"{c['i_max_A']:.1f}", f"iL{X}": f"{l['i_max_A']:.1f}",
              f"iD{X}": f"{100 * (l['i_max_A'] / c['i_max_A'] - 1):.0f}"}
    write("tab_medians.tex", rows_tex(rows))
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("figures, tables and numbers written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "lqg_report.tex"], cwd=HERE, check=True)


if __name__ == "__main__":
    main()
