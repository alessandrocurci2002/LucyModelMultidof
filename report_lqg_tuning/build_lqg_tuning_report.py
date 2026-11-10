"""Builds the living report on the tuning of the LQG-PI position controller (one section per step).

    ../.venv/bin/python -I build_lqg_tuning_report.py [--no-latex]

Inputs (nothing is computed here): ``results_lqg_tuning`` (``analisi_lqg_punto_lavoro.py``, ``analisi_maxwell_tempo.py``, ``analisi_lqg_passo1_flessibili.py``) and
``results_fork_modes_dahl_resampled`` (structure and effective gains). The text is in lqg_tuning_report.tex.
"""

from __future__ import annotations

import argparse
import json
import numpy as np
import math
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT / "results_lqg_tuning"
RIC = ROOT / "results_fork_modes_dahl_resampled"
FIG, GEN = HERE / "figures", HERE / "generated"
COPY = {"plant_operating_points": RES / "figures" / "plant_operating_points.pdf", "circle_bearing_state": RES / "figures" / "circle_bearing_state.pdf"}
AX = (("azimuth", "A"), ("elevation", "P"))
S24 = RES / "steps_2_4"
CM = RES / "comparison_measured"
NB = RES / "no_bearings"
COPY |= {f"nb1_{k}": NB / "step1_flexible" / "figures" / f"{k}.pdf" for k in ("first_guess", "process_uniform")}
COPY |= {f"nb2_{k}": NB / "steps_2_4" / "figures" / f"{k}.pdf" for k in ("step_2", "step_3", "final_tf_step", "effect_of_V", "final_doublet")}
COPY |= {f"cc_{k}": NB / "circle" / "figures" / f"{k}.pdf" for k in ("circle_nominal", "circle_scenarios", "circle_torque_step", "circle_all_disturbances", "disturbance_rejection", "twy_modes", "loop_decomposition")}
COPY |= {f"tp_{k}": NB / "thesis_pointing" / "figures" / f"{k}.pdf" for k in ("trajectories", "disturbance", "circle_miss", "dataset_7_miss", "dataset_8_miss", "effective_pointing", "miss_cdf", "speed")}
COPY |= {f"cm_{k}": CM / "figures" / f"{k}.pdf" for k in ("measured_doublets", "azimuth_chain", "elevation_checks")}
COPY |= {f"s24_{k}": S24 / "figures" / f"{k}.pdf" for k in ("step_2", "step_3", "final_tf_step", "effect_of_V", "final_doublet")}
TD = RES / "time_domain_maxwell"
COPY |= {f"td_{c}_{k}": TD / "figures" / f"{c}_{k}.pdf" for c in ("azimuth", "elevation") for k in ("open_loop", "closed_loop_estimated", "closed_loop_akd")}
S1 = RES / "step1_flexible"
COPY |= {f"s1_{k}": S1 / "figures" / f"{k}.pdf" for k in ("first_guess", "process_uniform", "process_mode_ladder", "process_trace", "dahl_check")}
LABEL = {"none (sliding)": "no bearings (Taylor, steady sliding)", "stuck (sigma0)": r"stuck, $\sigma_0$ (Taylor, at rest)",
         "describing function, measured amplitude": "describing function, test amplitude (not state space)",
         "viscous k_eq + c_eq": r"viscous, $k_{eq}+c_{eq}$",
         "Maxwell (spring-damper)": r"\textbf{Maxwell, $k_1s/(s+a)$}",
         "Maxwell with free k_0": r"Maxwell with free $k_0$ (not adopted)",
         "two Maxwell branches": "two Maxwell branches (not adopted)",
         "reduced synthesis model": r"\textbf{reduced synthesis model}"}
SINGLE = {"none (sliding)": "no bearings", "stuck (sigma0)": r"stuck, $\sigma_0$", "describing function, measured amplitude": "describing function, test amplitude",
          "viscous k_eq + c_eq": r"viscous, $k_{eq}+c_{eq}$", "Maxwell (spring-damper)": r"\textbf{Maxwell}", "measured": "measured"}


def damping(modes, lo=50.0, hi=105.0):
    return [m for m in modes if lo < m["f_hz"] < hi and m["zeta"] < 0.1]


def write(name, text):
    GEN.mkdir(parents=True, exist_ok=True)
    (GEN / name).write_text(text + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-latex", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    for n, src in COPY.items():
        shutil.copyfile(src, FIG / f"{n}.pdf")
    O = json.loads((RES / "operating_point.json").read_text())
    Ri = json.loads((RIC / "summary.json").read_text())
    Td = json.loads((RIC / "time_domain_summary.json").read_text())
    M = {}
    rows = []
    for key, lab in list(LABEL.items()) + [(f"describing function, A = {x} urad", rf"describing function, $A={x}\,\mu$rad") for x in (10, 100, 1000, 10000)]:
        row = [lab]
        for ch, X in AX:
            for e in O[ch]["errors"][key][:2]:
                row += [f"{e['rms_dB']:.1f}", f"{e['rms_deg']:.0f}"]
        rows.append(" & ".join(row))
    write("tab_candidates.tex", " \\\\\n".join(rows))
    rows = []
    for ch, X in AX:
        o, c, mx = O[ch], O[ch]["circle"], O[ch]["maxwell"]
        a = o["frf_amplitude_urad_1_10Hz"]
        rows.append(" & ".join([ch, f"{o['presliding_scale_urad']:.0f}", f"{a[0]:.0f}--{a[2]:.0f} ({a[1]:.0f})", f"{mx['k_1']:.0f}", f"{mx['f_a_Hz']:.1f}",
                                f"{mx['c']:.1f}", f"{o['c_eq']:.1f}", f"{o['origin_residue']:.2f}", str(o["synthesis"]["order"])]))
        M |= {f"ko{X}": f"{mx['k_1']:.0f}", f"fa{X}": f"{mx['f_a_Hz']:.1f}", f"pre{X}": f"{100 * c['presliding_fraction']:.0f}", f"rev{X}": str(c["reversals"]),
              f"vmax{X}": f"{c['velocity_max_rad_s']:.4f}", f"res{X}": f"{o['origin_residue']:.2f}", f"cm{X}": f"{mx['c']:.0f}", f"gapM{X}": f"{o['gain_at_1Hz_dB']['maxwell_over_measured']:.1f}",
              f"gapS{X}": f"{o['gain_at_1Hz_dB']['no_bearings_over_maxwell']:.0f}", f"tq{X}": f"{mx['c'] * math.radians(10):.0f}",
              f"tc{X}": f"{Ri[ch]['uniform_lnf']['hybrid']['parameters']['Tc']:.2f}", f"kzfree{X}": f"{o['free_k_0']['k_0']:.0f}",
              f"pkM{X}": f"{o['peaks_dB']['Maxwell (spring-damper)'][0] - o['measured_peak_Hz_dB'][1]:.1f}",
              f"pkV{X}": f"{o['measured_peak_Hz_dB'][1] - o['peaks_dB']['viscous k_eq + c_eq'][0]:.1f}",
              f"zM{X}": f"{damping(o['modes_with_bearings'])[0]['zeta']:.3f}", f"zV{X}": f"{damping(o['modes_viscous'])[0]['zeta']:.3f}",
              f"zN{X}": f"{damping(o['modes_no_bearings'])[0]['zeta']:.3f}",
              f"lfM{X}": f"{o['errors']['Maxwell (spring-damper)'][0]['rms_dB']:.1f}", f"lfN{X}": f"{o['errors']['none (sliding)'][0]['rms_dB']:.1f}"}
    write("tab_point.tex", " \\\\\n".join(rows))
    rows = []
    for ch, X in AX:
        for q in O[ch]["weight_scan"]:
            e, sel = q["errors"][0], q["weight"] == 3.0
            r = [f"{q['weight']:.0f}", f"{q['k_1']:.0f}", f"{q['f_a_Hz']:.2f}", f"{q['c']:.0f}", f"{e['rms_dB']:.1f}", f"{e['rms_deg']:.0f}", f"{q['peak_Hz']:.1f}", f"$-${abs(q['peak_dB']):.1f}", f"{q['zeta_mode']:.4f}"]
            rows.append(" & ".join([ch if q["weight"] == 1.0 else ""] + [rf"\textbf{{{x}}}" if sel else x for x in r]))
        rows.append(" & ".join([r"\emph{measured}", "", "", "", "", "", f"{O[ch]['measured_peak_Hz_dB'][0]:.1f}", f"$-${abs(O[ch]['measured_peak_Hz_dB'][1]):.1f}", ""]))
    write("tab_weights.tex", " \\\\\n".join(rows))
    rows = []
    for key, lab in (("none (sliding)", "no bearings"), ("viscous k_eq + c_eq", r"viscous, $k_{eq}+c_{eq}$"), ("describing function, measured amplitude", "describing function (not state space)"), ("Maxwell (spring-damper)", r"\textbf{Maxwell}"), ("reduced synthesis model", r"\textbf{reduced synthesis model}")):
        row = [lab]
        for ch, X in AX:
            q = O[ch]["resonance_band_errors"][key]
            row += [f"{q['rms_dB']:.1f}", f"{q['rms_deg']:.0f}"]
        rows.append(" & ".join(row))
    write("tab_resbands.tex", " \\\\\n".join(rows))
    rows = []
    for ch, X in AX:
        o = O[ch]
        rows.append(ch + " & " + "; ".join(f"{m['f_hz']:.2f}" + (f" ({m['zeta']:.3f})" if m["zeta"] < 0.999 else "") for m in o["synthesis"]["modes"] if m["f_hz"] > 0))
    write("tab_poles.tex", " \\\\\n".join(rows))
    rows = []
    for key, lab in SINGLE.items():
        row = [lab]
        for ch, X in AX:
            row += [f"{x:.1f}" for x in O[ch]["peaks_dB"][key]]
        rows.append(" & ".join(row))
    write("tab_peaks.tex", " \\\\\n".join(rows))
    rows = []
    for ch, X in AX:
        for z0, zv, zm in zip(damping(O[ch]["modes_no_bearings"]), damping(O[ch]["modes_viscous"]), damping(O[ch]["modes_with_bearings"])):
            rows.append(" & ".join([ch, f"{z0['f_hz']:.1f}", f"{z0['zeta']:.3f}", f"{zv['zeta']:.3f}", rf"\textbf{{{zm['zeta']:.3f}}}"]))
    write("tab_damping.tex", " \\\\\n".join(rows))
    rows = []
    for ch, X in AX:
        v, g = Ri[ch]["uniform_lnf"]["hybrid"]["parameters"], Td["gains"][ch]["uniform_lnf"]
        rows.append(" & ".join([ch, f"{v['sigma0']:.0f}", f"{v['Tc']:.3f}", f"{v['g_Kt']:.3f}", f"{g['Kpp']:.2f}", f"{g['Kpv']:.1f}", f"{g['Kiv_hz']:.2f}"]))
    write("tab_model.tex", " \\\\\n".join(rows))
    W = json.loads((TD / "summary.json").read_text())
    Med = W["medians_without_anomalous_window"]
    rows = []
    for ch, ax, nome in (("azimuth", "azimut", "azimuth"), ("elevation", "polare", "elevation")):
        m = Med[ax]
        for n, lab in (("dahl", "Dahl"), ("maxwell", r"\textbf{Maxwell}")):
            row = [nome if n == "dahl" else "", lab, f"{m['open_loop'][n]['fit_motion']:.0f}"]
            for g in ("estimated", "akd"):
                q = m["closed_loop"][g][n]
                row += [f"{q['fit_position_motion']:.0f}", f"{q['fit_current_motion']:.0f}", f"{q['i_max_model']:.1f}"]
            rows.append(" & ".join(x.replace("-", "$-$", 1) if x.startswith("-") else x for x in row))
    write("tab_td_medians.tex", " \\\\\n".join(rows))
    rows = []
    for ch in ("azimuth", "elevation"):
        g = W[ch]["gains"]
        for lab, x in (("files (.akd, WorkBench)", g["files"]), ("estimated, Dahl", g["dahl_estimated"]), (r"\textbf{estimated, Maxwell}", g["maxwell_estimated"])):
            rows.append(" & ".join([ch if x is g["files"] else "", lab, f"{x['Kpp']:.1f}", f"{x['Kpv']:.1f}", f"{x['Kiv_hz']:.2f}"]))
    write("tab_td_gains.tex", " \\\\\n".join(rows))
    M |= {f"meas{X}": f"{Med[a]['measured']['i_max']:.1f}" for a, X in (("azimut", "A"), ("polare", "P"))}
    M |= {f"peakM{X}": f"{W[c]['peak_i_cmd_A']['maxwell_estimated']:.0f}" for c, X in AX}
    M |= {f"peakDG{X}": f"{W[c]['peak_i_cmd_A']['maxwell_dahl_gains']:.0f}" for c, X in AX}
    neg = lambda x: f"$-${abs(x):.0f}" if x < 0 else f"{x:.0f}"
    for a, c, X in (("azimut", "azimuth", "A"), ("polare", "elevation", "P")):
        q, g = Med[a]["closed_loop"], W[c]["gains"]
        M |= {f"olD{X}": neg(Med[a]["open_loop"]["dahl"]["fit_motion"]), f"olM{X}": neg(Med[a]["open_loop"]["maxwell"]["fit_motion"]),
              f"clE{X}": neg(q["estimated"]["maxwell"]["fit_position_motion"]), f"clC{X}": neg(q["estimated"]["maxwell"]["fit_current_motion"]),
              f"clG{X}": neg(q["dahl_gains"]["maxwell"]["fit_position_motion"]), f"clGC{X}": neg(q["dahl_gains"]["maxwell"]["fit_current_motion"]),
              f"clF{X}": neg(q["akd"]["maxwell"]["fit_position_motion"]), f"clFD{X}": neg(q["akd"]["dahl"]["fit_position_motion"]),
              f"kivM{X}": f"{g['maxwell_estimated']['Kiv_hz']:.1f}", f"kivD{X}": f"{g['dahl_estimated']['Kiv_hz']:.2f}",
              f"kpvM{X}": f"{g['maxwell_estimated']['Kpv']:.1f}", f"kpvD{X}": f"{g['dahl_estimated']['Kpv']:.1f}"}
    F = json.loads((S1 / "summary.json").read_text())
    neg1 = lambda x: f"$-${abs(x):.1f}" if round(x, 1) < 0 else f"{abs(x):.1f}"
    polo = lambda q: "stable" if q["abs_z"] < 1 + 1e-9 else rf"\textbf{{unstable}}, {q['f_Hz']:.1f}\,Hz"
    exp10 = lambda x: f"{x:g}" if x < 1e3 else rf"$10^{{{math.log10(x):.0f}}}$" if abs(math.log10(x) % 1) < 1e-9 else rf"$3\cdot10^{{{math.floor(math.log10(x)):.0f}}}$"
    r_f, r_m, r_u, r_l, r_a, r_d, r_v, r_q = [], [], [], [], [], [], [], []
    for ch, X in AX:
        o = F[ch]
        fg, st = o["first_guess"], o["first_guess_step"]
        im = o["order_most_evident_first"][0]
        r_f.append(" & ".join([ch, f"{fg['k_p']:.2f}", f"{fg['k_i']:.2f}", f"{fg['bandwidth_Hz']:.2f}", f"{st['overshoot_pct']:.0f}", f"{st['t_set_2pct_s']:.1f}", f"{st['vib_peak_urad']:.1f}", f"{st['vib_after_100ms_urad']:.2f}"]))
        for m0, m in zip(o["modes"], fg["modes"]):
            r_m.append(" & ".join([ch, f"{m0['f_Hz']:.1f}", f"{m0['zeta']:.3f}", neg1(m["band_max_dB"]), f"{m['peak']['height_dB']:.1f}"]))
        for u in o["uniform_ladder"]:
            z = [f"{u['q_f']:g}", f"{u['k_p']:.2f}", f"{u['bandwidth_Hz']:.3f}", f"{u['step']['overshoot_pct']:.0f}", f"{u['step']['t_set_2pct_s']:.1f}" if u["step"]["t_set_2pct_s"] < 19 else r"$>$19", neg1(u["modes"][im]["band_max_dB"]), f"{u['modes'][im]['zeta_regulator']:.4f}", f"{u['modulus']:.2f}"]
            if u["q_f"] == 0.1:
                z = [rf"\textbf{{{q}}}" for q in z]
            r_u.append(" & ".join([ch if u["q_f"] == 0.0 else ""] + z))
        for x, v in zip(o["mode_ladder_levels"], o["mode_ladder"][str(im)]):
            m = v["modes"][im]
            r_l.append(" & ".join([ch if x == 0.1 else "", exp10(x), f"{m['zeta_regulator']:.4f}", neg1(m["band_max_dB"]), f"{m['peak']['height_dB']:.1f}" if m["peak"] else "none", f"{v['modulus']:.2f}", f"{v['bandwidth_Hz']:.3f}", f"{v['k_p']:.3f}"]))
        for lab, k in (("$q_f=0.1$", "q_f = 0.1 on all flexible states"), ("$q_f=0$", "q_f = 0")):
            q = o["dahl_check"][k]
            y, lp = q["y_100_120s_mdeg"], q["largest_pole"]
            r_d.append(" & ".join([ch if k.startswith("q_f = 0.1") else "", lab, polo(lp["full_maxwell"]), polo(lp["sliding"]), polo(lp["stuck"]), f"{q['step_20s']['overshoot_pct']:.0f}", f"{y[0]:.1f}--{y[1]:.1f}", f"{q['y_minus_y_hat_end_mdeg']:.2f}", f"{q['i_max_A']:.1f}"]))
        for w in o["V_discretisation"]:
            lp = w["largest_pole"]
            r_v.append(" & ".join([ch if w["q_f"] == 0 and w["estimator"].endswith("sample") else "", exp10(w["q_f"]) if w["q_f"] else "0", w["estimator"], f"{w['K_e_position_per_s']:.2f}", polo(lp["sliding"]), polo(lp["stuck"])]))
        for p in o["preview_higher_gains"]:
            r_a.append(" & ".join([ch if p["q_p"] == 25.0 else "", f"{p['q_p']:.0f}", f"{p['q_i']:.0f}", f"{p['k_p']:.1f}", f"{p['k_i']:.1f}", neg1(p["modes"][im]["band_max_dB"]), f"{p['modulus']:.2f}"]))
        fr = [v["state"].split(",")[0] for v in o["q_vector"][2:]]
        vis = list(dict.fromkeys(fr))
        r_q.append(" & ".join([ch, "0.1", "0.3"]) + " & " + r"; ".join(rf"{x}: {'$[0.1,\ 0.1]$' if fr.count(x) == 2 else '0.1'}" for x in vis))
        lv = o["mode_ladder_levels"].index(1e4)
        v4 = o["mode_ladder"][str(im)][lv]
        q0 = o["uniform_ladder"][0]
        M |= {f"fgKp{X}": f"{q0['k_p']:.2f}", f"fgPk{X}": neg1(fg["modes"][im]["band_max_dB"]), f"fgZ{X}": f"{fg['modes'][im]['zeta_regulator']:.4f}", f"fgBw{X}": f"{fg['bandwidth_Hz']:.2f}",
              f"vib{X}": f"{st['vib_peak_urad']:.1f}", f"kpDrop{X}": f"{100 * (1 - v4['k_p'] / fg['k_p']):.1f}", f"hiPk{X}": neg1(o["preview_higher_gains"][-1]["modes"][im]["band_max_dB"]), f"hiKp{X}": f"{o['preview_higher_gains'][-1]['k_p']:.0f}",
              f"bias{X}": f"{o['dahl_check']['q_f = 0.1 on all flexible states']['y_minus_y_hat_end_mdeg']:.2f}", f"zFour{X}": f"{v4['modes'][im]['zeta_regulator']:.3f}", f"pkFour{X}": neg1(v4["modes"][im]["band_max_dB"])}
    righe = []
    for r in json.loads((S1 / "weights_to_gains.json").read_text()):
        z = [f"{r['q_p']:g}", f"{r['q_i']:g}"]
        for x in (r["formula"], r["rigid_dlqr"], r["azimuth_maxwell"], r["elevation_maxwell"]):
            z += [f"{x['k_p']:.2f}", f"{x['k_i']:.2f}"]
        righe.append(" & ".join(z + [f"{r['azimuth_no_bearings']['k_p']:.2f}", f"{r['elevation_no_bearings']['k_p']:.2f}"]))
    write("tab_q_gains.tex", " \\\\\n".join(righe))
    for n, x in (("first", r_f), ("modes", r_m), ("uniform", r_u), ("mode", r_l), ("preview", r_a), ("dahl", r_d), ("v", r_v), ("vector", r_q)):
        write(f"tab_s1_{n}.tex", " \\\\\n".join(x))
    F2 = json.loads((S24 / "summary.json").read_text())
    def en(x):
        if x == 0:
            return "0"
        e = math.floor(math.log10(abs(x)))
        return f"{x:.3g}" if -2 <= e < 3 else rf"${x / 10 ** e:.2g}\cdot10^{{{e}}}$"
    lim_t = {"current": "current", "unstable on a full plant": "unstable on a full plant", "bandwidth at the fundamental": "fundamental", "min |1 + L|": r"$\min|1+L|$", "critical gain": "(7.9)", "oscillations": "oscillations", "resonance peak": "peak"}
    lim = lambda x: ", ".join(lim_t[q] for q in x.get("limits", [])) or "--"
    def riga(ch0, x, chiave, bold=False):
        d = x["doublets"]
        z = [en(x[chiave]), f"{x['k_p']:.4g}", f"{x['k_i']:.3g}", f"{x['bandwidth_Hz']:.3g}", f"{x['step']['t_set_2pct_s']:.3g}", f"{x['step']['overshoot_pct']:.1f}", f"{d['i_max_A']:.1f}", f"{x['modulus']:.2f}", neg1(x["modes"][0]["band_max_dB"]), lim(x)]
        if bold:
            z = [rf"\textbf{{{q}}}" for q in z]
        return " & ".join([ch0] + z)
    r2, r3, r4, rf_, rr, rv, rx, rq2 = [], [], [], [], [], [], [], []
    for ch, X in AX:
        o = F2[ch]
        for k, (rows_, key, sel) in {"p": (r2, "q_p", o["q_p_selected"]), "i": (r3, "q_i", o["q_i_selected"])}.items():
            pts = sorted(o[key], key=lambda x: x[key])
            giu = [x for x in pts if round(math.log10(x[key]) * 2) == math.log10(x[key]) * 2 or x[key] == sel or abs(math.log10(x[key]) - round(math.log10(x[key]) * 2) / 2) < 0.03 or x["limits"]]
            vista = []
            for x in giu:
                if x["limits"] and any(q["limits"] for q in vista):
                    continue
                vista.append(x)
            if len(vista) > 11:
                vista = vista[:3] + vista[-8:]
            for n, x in enumerate(vista):
                rows_.append(riga(ch if n == 0 else "", x, key, x[key] == sel))
        for x in o["step_4"]:
            rq = x["q_f"].get("0", "0.1")
            r4.append(" & ".join([ch, en(float(rq)) if rq != "0.1" else "0.1", neg1(x["modes"][0]["band_max_dB"]), neg1(x["modes"][1]["band_max_dB"]), f"{x['modes'][0]['zeta_regulator']:.3f}", f"{x['modulus']:.2f}", f"{x['bandwidth_Hz']:.3g}", f"{x['k_p']:.4g}"]))
        pt_g = sorted([x for x in o["fine_grid"] if not x["limits"]], key=lambda x: x["index"])[:4]
        for n, x in enumerate(pt_g):
            d = x["doublets"]
            rf_.append(" & ".join([ch if n == 0 else "", f"{x['index']:.3f}".replace("-", "$-$"), en(x["q_p"]), en(x["q_i"]), f"{x['k_p']:.3g}", f"{x['k_i']:.3g}", f"{x['bandwidth_Hz']:.3g}", f"{x['step']['t_set_2pct_s']:.3g}", f"{x['step']['overshoot_pct']:.1f}", f"{d['i_max_A']:.1f}"]))
        for n, x in enumerate(o["q_p_raised"]):
            rr.append(riga(ch if n == 0 else "", x, "q_p", x["q_p"] == o["q_p_raised_selected"]))
        dec = {e["rho"]: e for e in o["V_decision"]}
        for n, q in enumerate(o["V"]):
            x, m = q["at_same_q_p"], q["at_q_p_max"]
            d = x["doublets"]
            e = dec.get(q["rho"])
            rv.append(" & ".join([ch if n == 0 else "", f"{q['rho']:g}", f"{x['bandwidth_Hz']:.3g}", f"{x['step']['t_set_2pct_s']:.3g}", f"{x['modulus']:.2f}", f"{d['i_max_A']:.1f}", f"{d['t_set_median_s']:.2f}", f"{d['err_rms_median_mdeg']:.0f}", f"{d['over_max_mdeg']:.0f}", en(m["q_p"]) if m else "--", f"{m['bandwidth_Hz']:.3g}" if m else "--", f"{e['index']:.3f}".replace("-", "$-$") if e else "--"]))
        fin, c = o["final"], o["cascade"]
        dfin, dc, cd = fin["doublets"], o["dahl_check"], c["doublets"]
        for lab, v1, v2 in (("$k_p$ [1/s]; $K_{pp}$", f"{fin['k_p']:.4g}", f"{c['K_pp']:.2f}"), ("$k_i$ [1/s$^2$]", f"{fin['k_i']:.4g}", "--"), ("band, Maxwell [Hz]", f"{fin['bandwidth_Hz']:.3g}", f"{c['bandwidth_Hz']:.3g}"),
                            ("$t_s$, 1$^\\circ$, 2\\,\\%, Maxwell [s]", f"{fin['step']['t_set_2pct_s']:.3g}", f"{c['step']['t_set_2pct_s']:.3g}"), ("overshoot, Maxwell [\\%]", f"{fin['step']['overshoot_pct']:.1f}", f"{c['step']['overshoot_pct']:.1f}"),
                            ("largest current, doublets, Dahl [A]", f"{dfin['i_max_A']:.1f}", f"{cd['i_max_A']:.1f}"), ("median rms error, doublets, Dahl [mdeg]", f"{dfin['err_rms_median_mdeg']:.0f}", f"{cd['err_rms_median_mdeg']:.0f}"), ("median $t_\\mathrm{set}$, doublets, Dahl [s]", f"{dfin['t_set_median_s']:.2f}", f"{cd['t_set_median_s']:.2f}"),
                            ("largest overshoot, doublets, Dahl [mdeg]", f"{dfin['over_max_mdeg']:.0f}", f"{cd['over_max_mdeg']:.0f}")):
            rx.append((ch, lab, v1, v2))
        vis = [v["state"].split(",")[0] for v in F[ch]["q_vector"][2:]]
        qf56 = o["q_f_final"].get("0")
        f56 = f"{o['fundamental_Hz']:.2f} Hz"
        parti = []
        for x in dict.fromkeys(vis):
            val = en(qf56) if (qf56 and x == f56) else "0.1"
            parti.append(rf"{x}: $[{val.strip('$')},\ {val.strip('$')}]$" if vis.count(x) == 2 else f"{x}: {val}")
        rq2.append(" & ".join([ch, en(fin["q_i"]), en(fin["q_p"]), "; ".join(parti)]))
        m_ = o["final_measured_plant"]
        M |= {f"finKp{X}": f"{fin['k_p']:.4g}", f"finKi{X}": f"{fin['k_i']:.3g}", f"finBw{X}": f"{fin['bandwidth_Hz']:.3g}", f"finTs{X}": f"{fin['step']['t_set_2pct_s']:.3g}", f"finOv{X}": f"{fin['step']['overshoot_pct']:.0f}", f"finMod{X}": f"{fin['modulus']:.2f}", f"finI{X}": f"{dfin['i_max_A']:.1f}",
              f"finErr{X}": f"{dfin['err_rms_median_mdeg']:.0f}", f"finOvD{X}": f"{dfin['over_max_mdeg']:.0f}", f"casBw{X}": f"{c['bandwidth_Hz']:.3g}", f"casTs{X}": f"{c['step']['t_set_2pct_s']:.3g}", f"casErr{X}": f"{cd['err_rms_median_mdeg']:.0f}", f"casOvD{X}": f"{cd['over_max_mdeg']:.0f}", f"casI{X}": f"{cd['i_max_A']:.1f}",
              f"mMod{X}": f"{m_['modulus']:.2f}", f"mHz{X}": f"{m_['modulus_Hz']:.0f}", f"mL{X}": neg1(m_["L_max_dB"]), f"stepVib{X}": f"{dc['step_20s']['vib_after_100ms_urad']:.0f}", f"stepOv{X}": f"{dc['step_20s']['overshoot_pct']:.0f}", f"stepI{X}": f"{dc['i_max_A']:.0f}",
              f"sTwo{X}": en(o["q_p_selected"]), f"sThree{X}": en(o["q_i_selected"]), f"sRaised{X}": en(o["q_p_raised_selected"]), f"biasD{X}": f"{abs(dc['y_minus_y_hat_end_mdeg']):.2f}"}
    write("tab_s2.tex", " \\\\\n".join(r2)); write("tab_s3.tex", " \\\\\n".join(r3)); write("tab_s4.tex", " \\\\\n".join(r4)); write("tab_fine.tex", " \\\\\n".join(rf_)); write("tab_raise.tex", " \\\\\n".join(rr)); write("tab_v2.tex", " \\\\\n".join(rv))
    write("tab_fin_vector.tex", " \\\\\n".join(rq2))
    lab_ = list(dict.fromkeys(x[1] for x in rx))
    write("tab_final.tex", " \\\\\n".join(" & ".join([l] + [next(x[k] for x in rx if x[0] == ch and x[1] == l) for ch, _ in AX for k in (2, 3)]) for l in lab_))
    C = json.loads((CM / "summary.json").read_text())
    Nz = json.loads((CM / "noise_check.json").read_text())
    rm, rc, ra, re = [], [], [], []
    amp = {"01deg": "0.1", "05deg": "0.5", "1deg": "1", "3deg": "3", "5deg": "5", "8deg": "8"}
    for ch, X in AX:
        o = C[ch]
        fin_key = "Final LQG (Maxwell)"
        for n, (pm, pc, pf) in enumerate(zip(o["measured"]["per_amplitude"], o["chain"]["Cascade"]["per_amplitude"], o["chain"][fin_key]["per_amplitude"])):
            a_ = next(v for k, v in amp.items() if k in pm["name"])
            rm.append(" & ".join([ch if n == 0 else "", a_] + [f"{q:.0f}" for q in (pm["err_rms_mdeg"], pc["err_rms_mdeg"], pf["err_rms_mdeg"], pm["err_max_mdeg"], pc["err_max_mdeg"], pf["err_max_mdeg"])] + [f"{q:.1f}" for q in (pm["i_max_A"], pc["i_max_A"], pf["i_max_A"])]))
        righe_c = [("Measured (real cascade)", o["measured"])] + list(o["chain"].items())
        for n, (nome, x) in enumerate(righe_c):
            lenti = x.get("slow_poles_sliding", [])
            pol = "; ".join(f"{q['f_Hz']:.1f} ({q['zeta']:.2f})" for q in lenti if q["f_Hz"] < 8.0 and q["zeta"] < 0.99)[:40]
            rc.append(" & ".join([ch if n == 0 else "", nome, f"{x['i_max_A']:.1f}", f"{x['err_rms_median_mdeg']:.0f}", f"{x['over_max_mdeg']:.0f}", f"{x['t_set_median_s']:.2f}", f"{x['K_e_position_per_s']:.3g}" if "K_e_position_per_s" in x else "--", pol or "--"]))
    for x in C["azimuth"]["configurations"]:
        d = x["doublets"]
        ra.append(" & ".join([x["internal_model"], en(x["rho"]) if x["rho"] != 1.0 else "1", f"{x['q_i']:g}", en(x["q_p"]), f"{x['k_p']:.1f}", f"{x['k_i']:.2f}", f"{d['i_max_A']:.1f}", f"{d['err_rms_median_mdeg']:.0f}", f"{d['over_max_mdeg']:.0f}", f"{d['t_set_median_s']:.2f}", f"{x['bandwidth_sliding_Hz']:.2f}", f"{x['modulus_measured']:.2f}"]))
    for v, n_ in zip(C["elevation"]["checks"], Nz["elevation"]["ladder"]):
        re.append(" & ".join([en(v["q_p"]), f"{v['k_p']:.0f}", f"{v['nominal']['err_rms_median_mdeg']:.0f}", f"{v['nominal']['i_max_A']:.1f}", f"{v['nominal']['over_max_mdeg']:.0f}", f"{n_['noise']['10']['err_rms_median_mdeg']:.0f}", f"{n_['noise']['10']['i_max_A']:.1f}", f"{v['two_more_cycles']['err_rms_median_mdeg']:.0f}", "yes" if v["g_Kt_0.9_stable"] and v["g_Kt_1.1_stable"] else "no", f"{v['modulus_measured']:.2f}", neg1(v["T_ry_measured_max_dB"])]))
    Eq = json.loads((CM / "cascade_equivalence.json").read_text())
    re2 = []
    def req(ch0, nome, x, k=None):
        return " & ".join([ch0, nome, f"{k:.2f}" if k else "--", f"{x['i_max_A']:.1f}", f"{x['err_rms_median_mdeg']:.0f}", f"{x['over_max_mdeg']:.0f}", f"{x['t_set_median_s']:.2f}", f"{x['t_set_max_s']:.2f}"])
    for ch, X in AX:
        o = Eq[ch]
        re2.append(req(ch, "Cascade", o["cascade"], o["K_pp"]))
        re2.append(req("", "P on the encoder", o["P_encoder"], o["K_pp"]))
        x = o["same_k_p"]["0.1"]
        re2.append(req("", r"LQG, $q_p=K_{pp}^2$, $q_i=0.1$", x, x["k_p"]))
        for y in o["same_error"]:
            lab = rf"\textbf{{LQG, same error}}, $q_i={en(y['q_i']).strip('$')}$, $\rho={en(y['rho']).strip('$')}$"
            re2.append(req("", lab, y["lqg"], y["k_p"]))
        y = o["same_error"][0]
        re2.append(req("", r"same, $K_f=0$", y["without_K_f"], y["k_p"]))
    write("tab_cm_equiv.tex", " \\\\\n".join(re2))
    N1 = json.loads((NB / "step1_flexible" / "summary.json").read_text())
    N2 = json.loads((NB / "steps_2_4" / "summary.json").read_text())
    Ci = json.loads((NB / "circle" / "summary.json").read_text())
    Ex = json.loads((NB / "circle" / "extra.json").read_text())
    r_n1, r_n2, r_n3, r_n4, r_nr, r_nf, r_nq, r_c, r_cq, r_x = [], [], [], [], [], [], [], [], [], []
    for ch, X in AX:
        o, q = N1[ch], N2[ch]
        fg, st = o["first_guess"], o["first_guess_step"]
        r_n1.append(" & ".join([ch, f"{fg['k_p']:.2f}", f"{fg['k_i']:.2f}", f"{fg['bandwidth_Hz']:.2f}", f"{st['overshoot_pct']:.0f}", f"{st['t_set_2pct_s']:.1f}", neg1(fg["modes"][0]["band_max_dB"]), f"{fg['modes'][0]['peak']['height_dB']:.1f}", neg1(fg["modes"][1]["band_max_dB"]), f"{st['vib_peak_urad']:.1f}"]))
        for key, rows_, sel in (("q_p", r_n2, q["q_p_selected"]), ("q_i", r_n3, q["q_i_selected"])):
            pts = sorted(q[key], key=lambda x: x[key])
            vista = [x for x in pts if x[key] == sel or abs(math.log10(x[key]) * 2 - round(math.log10(x[key]) * 2)) < 0.05]
            pf = next((x for x in pts if x["limits"]), None)
            if pf is not None and pf not in vista:
                vista.append(pf)
            vista = sorted({id(x): x for x in vista if not (x["limits"] and x is not pf)}.values(), key=lambda x: x[key])
            if len(vista) > 11:
                vista = vista[:3] + vista[-8:]
            for n, x in enumerate(vista):
                rows_.append(riga(ch if n == 0 else "", x, key, x[key] == sel))
        for x in q["step_4"]:
            rq = x["q_f"].get("0", "0.1")
            r_n4.append(" & ".join([ch, en(float(rq)) if rq != "0.1" else "0.1", neg1(x["modes"][0]["band_max_dB"]), neg1(x["modes"][1]["band_max_dB"]), f"{x['modes'][0]['zeta_regulator']:.3f}", f"{x['modulus']:.2f}", f"{x['bandwidth_Hz']:.3g}", f"{x['k_p']:.4g}"]))
        for n, x in enumerate(q["q_p_raised"]):
            r_nr.append(riga(ch if n == 0 else "", x, "q_p", x["q_p"] == q["q_p_raised_selected"]))
        fin, c = q["final"], q["cascade"]
        df, cd = fin["doublets"], c["doublets"]
        r_nf.append((ch, ("$k_p$ [1/s]; $K_{pp}$", f"{fin['k_p']:.4g}", f"{c['K_pp']:.2f}"), ("$k_i$ [1/s$^2$]", f"{fin['k_i']:.3g}", "--"), ("band [Hz]", f"{fin['bandwidth_Hz']:.3g}", f"{c['bandwidth_Hz']:.3g}"), (r"$t_s$, 1$^\circ$, 2\,\% [s]", f"{fin['step']['t_set_2pct_s']:.3g}", f"{c['step']['t_set_2pct_s']:.3g}"), (r"overshoot [\%]", f"{fin['step']['overshoot_pct']:.1f}", f"{c['step']['overshoot_pct']:.1f}"),
                     ("largest current, doublets, Dahl [A]", f"{df['i_max_A']:.1f}", f"{cd['i_max_A']:.1f}"), ("median rms error, doublets, Dahl [mdeg]", f"{df['err_rms_median_mdeg']:.0f}", f"{cd['err_rms_median_mdeg']:.0f}"), (r"median $t_\mathrm{set}$, doublets, Dahl [s]", f"{df['t_set_median_s']:.2f}", f"{cd['t_set_median_s']:.2f}"), ("largest overshoot, doublets, Dahl [mdeg]", f"{df['over_max_mdeg']:.0f}", f"{cd['over_max_mdeg']:.0f}"), (r"$\min|1+L|$, internal model", f"{fin['modulus']:.2f}", "--")))
        vis = [v["state"].split(",")[0] for v in o["q_vector"][2:]]
        qf56 = q["q_f_final"].get("0")
        f56 = f"{q['fundamental_Hz']:.2f} Hz"
        parti = []
        for x in dict.fromkeys(vis):
            val = en(qf56) if (qf56 and x == f56) else "0.1"
            parti.append(rf"{x}: $[{val.strip('$')},\ {val.strip('$')}]$" if vis.count(x) == 2 else f"{x}: {val}")
        r_nq.append(" & ".join([ch, en(fin["q_i"]), en(fin["q_p"]), "; ".join(parti)]))
        M |= {f"nbKp{X}": f"{fin['k_p']:.4g}", f"nbKi{X}": f"{fin['k_i']:.3g}", f"nbQp{X}": en(q["q_p_selected"]), f"nbQi{X}": en(q["q_i_selected"]), f"nbQr{X}": en(q["q_p_raised_selected"]), f"nbMm{X}": f"{q['final_measured_plant']['modulus']:.2f}"}
        ee = Ex[ch]["cases"]
        for et, lab in (("cascade", "cascade"), ("LQG-PI tuned", "LQG-PI tuned"), ("LQG-PI, k_i = 0.5 k_p", r"$k_i\approx0.5k_p$"), ("LQG-PI, V = 100 Q", "$V=100Q$"), ("LQG-PI, q_p = 1e5", "$q_p=10^5$")):
            if et in ee:
                v = ee[et]
                d = v["doublets"]
                r_x.append(" & ".join([ch if et == "cascade" else "", lab, neg1(v["L_at_1Hz_dB"]), neg1(v["T_wy_at_1Hz_dB"]), neg1(v["T_wy_peak_dB"]), f"{v['bandwidth_T_ry_sliding_Hz']:.3g}", f"{d['err_rms_median_mdeg']:.0f}", f"{d['i_max_A']:.1f}", f"{d['over_max_mdeg']:.0f}", f"{d['t_set_median_s']:.2f}"]))
    lf = [x[0] for x in r_nf[0][1:]]
    write("tab_nb_final.tex", " \\\\\n".join(" & ".join([l, r_nf[0][1 + i][1], r_nf[0][1 + i][2], r_nf[1][1 + i][1], r_nf[1][1 + i][2]]) for i, l in enumerate(lf)))
    sc_l = {"nominal": "nominal", "noise 1 urad": r"noise 1\,$\mu$rad", "noise 10 urad": r"noise 10\,$\mu$rad", "torque step 0.5 N m": r"step 0.5\,N\,m", "Gaussian 0.1 N m rms": r"Gaussian 0.1\,N\,m", "Gaussian 0.1 N m rms, mean 0.2 N m": r"same, mean 0.2\,N\,m", "sine 0.2 N m, 1 Hz": r"sine 0.2\,N\,m, 1\,Hz", "sine 0.2 N m, 10 Hz": r"sine 0.2\,N\,m, 10\,Hz", "sine 0.2 N m, 56 Hz": r"sine 0.2\,N\,m, 56\,Hz", "all: 1 urad, step, Gaussian with mean": "noise, step, Gaussian with mean"}
    for sn, lab in sc_l.items():
        x = Ci["scenarios"][sn]
        a_, b_, c_, d_ = x["cascade"], x["LQG-PI tuned"], x["LQG-PI, V = 100 Q"], x["LQG-PI, q_p = 1e5"]
        fm = lambda y: f"{y:.3g}"
        r_c.append(" & ".join([lab] + [fm(q["E_rms_mdeg"]) for q in (a_, b_, d_)] + [fm(q["dist_total_rms_urad"]) for q in (a_, b_, c_, d_)] + [fm(q["dist_total_peak_urad"]) for q in (a_, b_, d_)] + [fm(q["jitter_10Hz_urad"]) for q in (a_, b_, d_)]))
    r_s, r_l2, r_qf, r_cp = [], [], [], []
    for ch, X in AX:
        sn = Ex[ch]["sines"]
        for n, f_ in enumerate(q["f_Hz"] for q in sn["cascade"]):
            c0, l0 = sn["cascade"][n], sn["LQG-PI tuned"][n]
            r_s.append(" & ".join([ch if n == 0 else "", f"{f_:g}"] + [f"{y:.3g}" for y in (c0["encoder_urad"], c0["line_of_sight_urad"], c0["linear_sliding_urad"], c0["linear_stuck_urad"], l0["encoder_urad"], l0["line_of_sight_urad"], l0["linear_sliding_urad"], l0["linear_stuck_urad"])]))
        d = Ex[ch]["loop_decomposition"]
        for lab, k in (("cascade", "cascade"), ("LQG-PI", "LQG"), ("$K_f=0$", "K_f = 0"), ("fast estimator", "fast estimator"), ("fast estimator, $K_f=0$: $k_p$ alone", "fast estimator, K_f = 0")):
            r_l2.append(" & ".join([ch if k == "cascade" else "", lab] + [neg1(y) for y in d[k]]))
    sa = Ex["azimuth"]["sines"]
    for et, lab in (("LQG-PI tuned", "0.1"), ("LQG-PI, q_f = 10000", "$10^4$"), ("LQG-PI, q_f = 30000", r"$3\cdot10^4$"), ("LQG-PI, q_f = 100000", "$10^5$")):
        q1 = next(y for y in sa[et] if y["f_Hz"] == 1.0)
        q56 = next(y for y in sa[et] if y["f_Hz"] == 56.0)
        st_ = Ex["azimuth"]["sines_largest_pole"][et]
        r_qf.append(" & ".join([lab, f"{q56['encoder_urad']:.1f}", f"{q56['line_of_sight_urad']:.1f}", f"{q1['encoder_urad']:.0f}", "yes" if max(st_.values()) < 1 + 1e-9 else "no"]))
    for y in Ex["coupling"]["rows"]:
        r_cp.append(" & ".join([f"{y['f_Hz']:g}"] + [neg1(y[k]) for k in ("az_torque_to_az", "el_torque_to_el", "el_torque_to_az", "az_torque_to_el")]))
    Dd = json.loads((NB / "circle" / "dither.json").read_text())
    r_d2 = []
    for ch, X in AX:
        for et, lab in (("cascade", "cascade"), ("LQG-PI tuned", "LQG-PI tuned")):
            y = Dd[ch][et]
            r_d2.append(" & ".join([ch if et == "cascade" else "", lab] + [f"{y[k]:.3g}" for k in ("low_urad", "low_far_urad", "high_urad", "error_rms_nominal_mdeg", "error_rms_with_sine_mdeg")]))
    r_ce = []
    for nome, lab in (("test circle", "test circle"), ("fast, high circle (centre 15 m, height 15 m, 0.5 rad/s)", "fast, high circle")):
        q = Ci["coupling_effect"][nome]
        for et, l2 in (("cascade", "cascade"), ("LQG-PI tuned", "LQG-PI tuned")):
            y = q[et]
            r_ce.append(" & ".join([lab if et == "cascade" else "", f"{q['psi_rate_max']:.3g}" if et == "cascade" else "", f"{q['theta_mean_deg']:.0f}" if et == "cascade" else "", l2, en(y["torque_az_max_N_m"]), en(y["torque_el_max_N_m"]), en(y["torques_only_los_az_rms_urad"]), en(y["torques_only_los_el_rms_urad"]), en(y["J_only_los_az_rms_urad"]), en(y["H_only_los_az_rms_urad"]), en(y["H_only_los_el_rms_urad"]), en(y["effect_los_az_peak_urad"]), f"{y['i_az_max_fixed_J_A']:.1f}", f"{y['i_az_max_A']:.1f}"]))
    r_st = []
    for x in Ci["stability_vs_elevation"]:
        z = [f"{x['theta_deg']:.0f}", f"{x['J_ratio']:.3f}"]
        for et in ("cascade", "LQG-PI tuned", "LQG-PI, q_p = 1e5"):
            for n in ("sliding", "stuck"):
                v = x[f"{et}, {n}"]
                z.append("stable" if v < 1 + 1e-9 else rf"\textbf{{{v:.4f}}}")
        r_st.append(" & ".join(z))
    write("tab_cc_stab_theta.tex", " \\\\\n".join(r_st))
    Tp = json.loads((NB / "thesis_pointing" / "summary.json").read_text())
    Td_ = json.loads((NB / "thesis_pointing" / "derived.json").read_text())
    r_p, r_q, r_m = [], [], []
    for tag, tr_, st_ in (("circle", "circle", "--"), ("dataset7_EKF2", r"\texttt{dataset\_7}", "EKF2"), ("dataset7_config_v3", r"\texttt{dataset\_7}", r"\texttt{config\_v3}"), ("dataset8_EKF2", r"\texttt{dataset\_8}", "EKF2"), ("dataset8_config_v3", r"\texttt{dataset\_8}", r"\texttt{config\_v3}")):
        for et, lab, cod in (("cascade", "cascade", "cascade"), ("LQG-PI tuned", "LQG-PI", "lqg")):
            sel = [c for c in Tp["cases"] if c["traj"] == tag and c["controller"] == et]
            n0 = next(c for c in sel if c["seed"] is None)
            ds_ = [c for c in sel if c["seed"] is not None and not c.get("mean_N_m")]
            dm_ = [c for c in sel if c.get("mean_N_m")]
            if dm_:
                em_ = [c["on_target_pct"] for c in dm_]
                e0_ = np.mean([c["on_target_pct"] for c in ds_])
                mf = Td_["mean_effect"][f"{tag}__{cod}"]
                mv = lambda k: float(np.mean([y[k] for y in mf]))
                r_m.append(" & ".join([*((tr_, st_) if cod == "cascade" else ("", "")), lab, f"{e0_:.1f}", rf"\textbf{{{np.mean(em_):.1f}}} ({min(em_):.1f}--{max(em_):.1f})", "0.0" if abs(np.mean(em_) - e0_) < 0.05 else f"{np.mean(em_) - e0_:+.1f}",
                                       f"{mv('los_az_mean_urad'):.1f}", f"{abs(mv('los_el_mean_urad')):.1f}", f"{max(y['total_peak_urad'] for y in mf):.1f}", f"{mv('i_az_mean_A'):.2f}", f"{mv('i_el_mean_A'):.2f}"]))
            med = lambda k: float(np.mean([c[k] for c in ds_]))
            et_ = [c["on_target_pct"] for c in ds_]
            prim = (tr_, st_) if cod == "cascade" else ("", "")
            ef = Td_["effect_urad"][f"{tag}__{cod}"]
            r_p.append(" & ".join([*prim, lab, f"{n0['E_rms_mrad']:.2f}", f"{n0['E_max_mrad']:.2f}", f"{n0['miss_rms_cm']:.1f}", f"{n0['miss_p95_cm']:.1f}", rf"\textbf{{{n0['on_target_pct']:.1f}}}", f"{med('E_rms_mrad'):.2f}", f"{med('miss_rms_cm'):.1f}", rf"\textbf{{{np.mean(et_):.1f}}} ({min(et_):.1f}--{max(et_):.1f})", f"{np.mean(ef):.0f}"]))
            ea = n0["estimator_alone"]
            ate = "--" if tag == "circle" else f"{100 * Tp['geometry'][tag[:8]]['estimators'][st_.replace(chr(92) + 'texttt{', '').replace(chr(92) + '_', '_').replace('}', '')]['ate_m']:.1f}"
            r_q.append(" & ".join([*prim, lab, ate if cod == "cascade" else "", f"{ea['miss_rms_cm']:.1f}" if cod == "cascade" else "", f"{ea['on_target_pct']:.1f}" if cod == "cascade" else "", f"{n0['loop_E_rms_mrad']:.2f}", f"{n0['los_az_rms_mrad']:.2f}", f"{n0['los_el_rms_mrad']:.2f}", f"{n0['i_max_az_A']:.1f}", f"{n0['i_max_el_A']:.1f}", f"{med('i_max_az_A'):.1f}", f"{med('i_max_el_A'):.1f}"]))
    write("tab_tp_main.tex", " \\\\\n".join(r_p))
    write("tab_tp_loop.tex", " \\\\\n".join(r_q))
    write("tab_tp_mean.tex", " \\\\\n".join(r_m))
    Vs = json.loads((NB / "thesis_pointing" / "speed.json").read_text())
    r_v = []
    for ds, lab in (("dataset7", r"\texttt{dataset\_7}"), ("dataset8", r"\texttt{dataset\_8}")):
        righe_v = [(1, {c["controller"]: c for c in Tp["cases"] if c["traj"] == f"{ds}_config_v3" and c["seed"] is None})]
        righe_v += [(fv, {c["controller"]: c for c in Vs["cases"] if c["traj"] == f"{ds}_config_v3_x{fv}"}) for fv in Vs["factors"]]
        il = lambda c: rf"\textbf{{{c['i_max_el_A']:.1f}}}" if max(c["i_max_az_A"], c["i_max_el_A"]) > 18.0 else f"{c['i_max_el_A']:.1f}"
        for n_, (fv, d) in enumerate(righe_v):
            a_, b_ = d["cascade"], d["LQG-PI tuned"]
            r_v.append(" & ".join([lab if n_ == 0 else "", f"{fv}", f"{fv * Vs['mean_speed_m_s'][ds]:.2f}", f"{a_['estimator_alone']['on_target_pct']:.1f}", f"{a_['on_target_pct']:.1f}", f"{b_['on_target_pct']:.1f}",
                                   f"{b_['on_target_pct'] - a_['on_target_pct']:+.1f}", f"{a_['loop_E_rms_mrad']:.2f}", f"{b_['loop_E_rms_mrad']:.2f}", il(a_), il(b_)]))
    write("tab_tp_speed.tex", " \\\\\n".join(r_v))
    Bn = json.loads((NB / "circle" / "benefits.json").read_text())
    Bm = json.loads((NB / "circle" / "benefits_mid_band.json").read_text())
    def nome_caso(et):
        return et.replace("P, k = ", "P, $k=$").replace("LQG-PI, q_p = 1e5", "LQG-PI, $q_p=10^5$").replace("LQG-PI, V_f = 1e6", "LQG-PI, $V_f=10^6$")
    r_b1, r_b2 = [], []
    for asse, tit in (("azimut", "azimuth"), ("polare", "elevation")):
        for n_, (et, y) in enumerate(Bn[asse].items()):
            l_, d_ = y["linear"]["sliding"], y["doublets"]
            mz = Bm[asse][et]
            pm = mz["margins"]["sliding"]["pm_deg"]
            r_b1.append(" & ".join([tit if n_ == 0 else "", nome_caso(et), f"{pm:.0f}", f"{l_['band_Hz']:.1f}", f"{l_['T_ry_peak_dB']:.1f}", f"{d_['err_rms_median_mdeg']:.0f}", f"{d_['i_max_A']:.1f}", f"{d_['over_max_mdeg']:.0f}", f"{d_['t_set_median_s']:.2f}", f"{y['error_rms_mdeg']:.1f}"]))
            sn = y["scenarios"]
            se = lambda am, fs: next(z for z in mz["sines"] if z["amplitude_N_m"] == am and z["f_Hz"] == fs)["line_of_sight_urad"]
            r_b2.append(" & ".join([tit if n_ == 0 else "", nome_caso(et)] + [f"{sn[k]['effect_rms_urad']:.1f}" for k in ("noise 1 urad", "torque step 0.5 N m", "process noise 0.1 N m rms", "sine 0.2 N m, 1 Hz")] + [f"{se(0.2, 6.6):.0f}", f"{se(0.8, 3.0):.0f}", f"{se(0.8, 6.6):.0f}"]))
    Bw_ = json.loads((NB / "circle" / "benefits_wind.json").read_text())
    r_w = []
    for asse, tit in (("azimut", "azimuth"), ("polare", "elevation")):
        cw = Bw_[asse]["cascade"]
        media = lambda y, liv, k: (y[f"{liv}_N_m_seed_1"][k] + y[f"{liv}_N_m_seed_2"][k]) / 2
        for n_, (et, y) in enumerate(Bw_[asse].items()):
            pm = Bm[asse][et]["margins"]["sliding"]["pm_deg"] if et in Bm[asse] else None
            r_w.append(" & ".join([tit if n_ == 0 else "", nome_caso(et).replace("LQG-PI, k_i = 0.5 k_p", "LQG-PI, $k_i=0.5k_p$"), "--" if pm is None else f"{pm:.0f}", f"{media(y, '0.1', 'los_urad'):.1f}", f"{media(y, '0.1', 'encoder_urad'):.1f}", f"{media(cw, '0.1', 'los_urad') / media(y, '0.1', 'los_urad'):.1f}", f"{media(y, '0.5', 'los_urad'):.1f}", f"{media(cw, '0.5', 'los_urad') / media(y, '0.5', 'los_urad'):.1f}", f"{media(y, '0.5', 'i_max_A'):.1f}"]))
    write("tab_cc_benefits_wind.tex", " \\\\\n".join(r_w))
    write("tab_cc_benefits_cmd.tex", " \\\\\n".join(r_b1))
    write("tab_cc_benefits_dist.tex", " \\\\\n".join(r_b2))
    Tm = json.loads((NB / "circle" / "twy_modes.json").read_text())
    Tv = json.loads((NB / "circle" / "twy_modes_vf.json").read_text())
    Tr = json.loads((NB / "circle" / "twy_modes_rest.json").read_text())
    To = json.loads((NB / "circle" / "twy_modes_orbit.json").read_text())
    def polo_56(lista):
        p = [x for x in lista if min(abs(x[0] - y) for y in (49.4, 96.4, 98.3)) > 0.5]
        return f"{p[0][0]:.1f}, {p[0][1]:.3f}" if p else r"$>400$"
    r_tm = []
    for et, q_, v_ in (("Q_f = V_f = 0.1", "0.1", "0.1"), ("tuned, Q_f = V_f = 1e4", "$10^4$", "$10^4$"), ("Q_f = 0.1, V_f = 1e6", "0.1", "$10^6$"), ("Q_f = 1e4, V_f = 1e6", "$10^4$", "$10^6$"), ("Q_f = 1e4, V_f = 1e8", "$10^4$", "$10^8$"), ("Q_f = 1e4, V_f = 1e12", "$10^4$", "$10^{12}$"), ("Q_f = V_f = 1e6", "$10^6$", "$10^6$")):
        x = Tm["cases"][et]
        zs = x["stuck"]["largest_z"]
        r_tm.append(" & ".join([q_, v_, polo_56(x["regulator"]), polo_56(x["estimator"]), f"{x['T_ry_peak_40_70_dB']:.1f}", f"{x['sliding']['T_wy_56_dB'][0]:.1f}", f"{x['sliding']['T_wy_56_dB'][1]:.1f}", f"{x['stuck']['T_wy_56_dB'][0]:.1f}", f"{x['stuck']['T_wy_56_dB'][1]:.1f}", "stable" if zs < 1 + 1e-9 else rf"\textbf{{{zs:.4f}}}"]))
    write("tab_cc_twy_modes.tex", " \\\\\n".join(r_tm))
    r_td = []
    for et, lab in (("cascade", "cascade"), ("tuned, V_f = 1e4", r"tuned, $V_f=10^4$"), ("V_f = 1e6", r"$V_f=10^6$")):
        mo, ri = Tv[et], Tr[et]
        m_ = lambda fs, am: next(r for r in mo["sine"] if r["f_Hz"] == fs and r["amplitude_N_m"] == am)["dB"]
        r_ = lambda fs, am: next(r for r in ri if r["f_Hz"] == fs and r["amplitude_N_m"] == am)
        d_ = mo["doublets"]
        r_td.append(" & ".join([lab, f"{m_(57.0, 0.01):.1f}", f"{m_(57.0, 0.2):.1f}", f"{r_(57.8, 0.002)['dB']:.1f}", f"{r_(57.8, 0.2)['dB']:.1f}", f"{r_(58.5, 0.002)['dB']:.1f}", f"{r_(57.8, 0.002)['line_of_sight_urad']:.3f}", f"{d_['err_rms_median_mdeg']:.1f}", f"{d_['i_max_A']:.1f}", f"{d_['t_set_median_s']:.2f}", f"{To[et]['error_rms_mdeg']:.2f}", f"{To[et]['noise_1_urad_los_effect_urad']:.1f}", f"{To[et]['noise_10_urad_los_effect_urad']:.1f}"]))
    write("tab_cc_twy_dahl.tex", " \\\\\n".join(r_td))
    write("tab_cc_coupled.tex", " \\\\\n".join(r_ce))
    write("tab_cc_sines.tex", " \\\\\n".join(r_s)); write("tab_cc_loop.tex", " \\\\\n".join(r_l2)); write("tab_cc_qf.tex", " \\\\\n".join(r_qf)); write("tab_cc_coupling.tex", " \\\\\n".join(r_cp)); write("tab_cc_dither.tex", " \\\\\n".join(r_d2))
    write("tab_nb_step1.tex", " \\\\\n".join(r_n1)); write("tab_nb_s2.tex", " \\\\\n".join(r_n2)); write("tab_nb_s3.tex", " \\\\\n".join(r_n3)); write("tab_nb_s4.tex", " \\\\\n".join(r_n4)); write("tab_nb_raise.tex", " \\\\\n".join(r_nr)); write("tab_nb_vector.tex", " \\\\\n".join(r_nq)); write("tab_cc.tex", " \\\\\n".join(r_c)); write("tab_cc_x.tex", " \\\\\n".join(r_x))
    write("tab_cm_doublets.tex", " \\\\\n".join(rm)); write("tab_cm_chain.tex", " \\\\\n".join(rc)); write("tab_cm_conf.tex", " \\\\\n".join(ra)); write("tab_cm_el.tex", " \\\\\n".join(re))
    write("numbers.tex", "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in M.items()))
    print("written")
    if not args.no_latex:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "lqg_tuning_report.tex"], cwd=HERE, check=True)


if __name__ == "__main__":
    main()
