"""Doppietti V2 nel tempo con i parametri stimati sulle FRF ricampionate, a confronto con quelli del report.

Richiesta dell'utente (2026-11-01): confronto nel tempo, in anello aperto (corrente misurata in ingresso) e in anello
chiuso (riferimento di posizione in ingresso), delle stime ibride di ``analisi_dahl_ricampionato.py`` (griglia uniforme
in ln f e in f) con quella del report ``report_dahl_bearings``.

Simulatori e metriche sono quelli di ``analisi_tempo_dahl_ibrido.py`` e ``analisi_anello_chiuso_ibrido.py``. In anello
chiuso Kpp, Kpv, Kiv sono stimati per ogni variante con la stessa procedura del report (doppietto da 1 grado, plant
fisso); per la variante del report si riusano i guadagni salvati, e il confronto con i CSV del report controlla che il
calcolo sia lo stesso. Le stime sulle FRF non vengono ripetute.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_tempo_ricampionato.py
"""

from __future__ import annotations

import csv
import dataclasses
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))
sys.path.insert(0, str(RADICE / "closedLoopModel"))

import matplotlib.pyplot as plt  # noqa: E402
from scipy.signal import savgol_filter  # noqa: E402

import analisi_anello_chiuso_ibrido as aac  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import doppietti_tempo as dt  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl.driver import Drive, DriveParams  # noqa: E402
from lucy_cl.motor import MotorParams  # noqa: E402

CARTELLA = RADICE / "results_fork_modes_dahl_resampled"
SOMMARIO = CARTELLA / "summary.json"
REPORT_OL = RADICE / "results_fork_modes_dahl" / "time_domain_check" / "position_fit.csv"
REPORT_CL = RADICE / "results_fork_modes_dahl" / "closed_loop_check" / "stima"
VARIANTI = {"report": "report: ln f width weights", "uniform_lnf": "uniform in ln f", "uniform_f": "uniform in f"}
ANOMALA = ("azimut", "8deg_minus8deg", "#1")
"""Finestra esclusa dalle mediane: fit di posizione negativo (-136 %) gia' nei risultati del report."""


def titolo_prova(key) -> str:
    return f"{key[0].split('_')[0].replace('01deg', '0.1deg').replace('05deg', '0.5deg')} doublet, move {key[1]}"


def figura_aperto(nome, curve, seg):
    fig, ax = plt.subplots(3, 3, figsize=(15, 9), sharex="col", gridspec_kw={"height_ratios": [2, 1.2, 1]})
    for j, key in enumerate(aac.FIGURA):
        w, i0 = curve["report"][key][0], curve["report"][key][2]
        k = w.t <= w.t[i0] + 2.5
        t = w.t[k]
        md = lambda q: np.rad2deg(q[k]) * 1e3  # noqa: E731
        ax[0, j].plot(t, md(w.cmd), color=pt.INCHIOSTRO, lw=0.8, ls=":", label="command")
        ax[0, j].plot(t, md(w.pos), color=pt.INCHIOSTRO, lw=1.6, label="measured encoder")
        ax[1, j].plot(t, savgol_filter(w.pos, 41, 3, deriv=1, delta=dt.DT)[k], color=pt.INCHIOSTRO, lw=1.2)
        for i, (var, et) in enumerate(VARIANTI.items(), 1):
            y = curve[var][key][1]
            ax[0, j].plot(t, md(y), lw=1.3, label=f"model, {et}", **pt._stile(i))
            ax[1, j].plot(t, np.gradient(y, dt.DT)[k], lw=1.0, **pt._stile(i))
        ax[2, j].plot(t, w.I[k], color=pt.INCHIOSTRO, lw=1.0)
        ax[0, j].set_title(titolo_prova(key), fontsize=9)
        ax[2, j].set_xlabel("time from window start [s]")
    ax[0, 0].set_ylabel("position [mdeg]"); ax[1, 0].set_ylabel("velocity [rad/s]")
    ax[2, 0].set_ylabel("input current [A rms]")
    ax[0, 1].legend(loc="best", fontsize=7)
    pt._titolo(fig, f"{nome}: open loop, measured current input, hybrid Dahl parameters from the three FRF costs", seg)
    return fig


def figura_chiuso(nome, curve, seg):
    fig, ax = plt.subplots(3, 3, figsize=(15, 9.5), sharex="col", gridspec_kw={"height_ratios": [2, 1.3, 1.6]})
    for j, key in enumerate(aac.FIGURA):
        w, i0 = curve["report"][key][0], curve["report"][key][3]
        k = (w.t >= w.t[i0] - 0.2) & (w.t <= w.t[i0] + 2.0)
        t = w.t[k]
        md = lambda q: np.rad2deg(q[k]) * 1e3  # noqa: E731
        ax[0, j].plot(t, md(w.cmd), color=pt.INCHIOSTRO, lw=0.8, ls=":", label="position command (model input)")
        ax[0, j].plot(t, md(w.pos), color=pt.INCHIOSTRO, lw=1.6, label="measured encoder")
        ax[1, j].plot(t, md(w.pos - w.cmd), color=pt.INCHIOSTRO, lw=1.3)
        ax[2, j].plot(t, w.I[k], color=pt.INCHIOSTRO, lw=1.1, label="measured (AMPS / 2)")
        for i, (var, et) in enumerate(VARIANTI.items(), 1):
            _, y, im, _ = curve[var][key]
            ax[0, j].plot(t, md(y), lw=1.3, label=f"model, {et}", **pt._stile(i))
            ax[1, j].plot(t, md(y - (w.cmd - w.cmd[0])), lw=1.1, **pt._stile(i))
            ax[2, j].plot(t, im[k], lw=1.0, **pt._stile(i))
        ax[0, j].set_title(titolo_prova(key), fontsize=9)
        ax[2, j].set_xlabel("time from window start [s]")
    ax[0, 0].set_ylabel("position [mdeg]"); ax[1, 0].set_ylabel("following error [mdeg]")
    ax[2, 0].set_ylabel("motor current [A rms]")
    ax[0, 1].legend(loc="lower right", fontsize=7); ax[2, 0].legend(fontsize=7)
    pt._titolo(fig, f"{nome}: closed loop, measured position command input, gains estimated for each parameter set", seg)
    return fig


def controllo_report(righe, percorso, colonna, indice) -> float:
    """Massimo scarto, in punti di fit, fra la variante del report e il CSV del report."""
    with open(percorso, newline="") as fh:
        rif = {(r["axis"], r["test"], r["move"]): float(r[colonna]) for r in csv.DictReader(fh)}
    return max(abs(r[indice] - rif[(r[0], r[1], r[2])]) for r in righe if r[3] == "report")


def mediana(righe, asse, var, indice) -> float:
    return float(np.median([r[indice] for r in righe if r[0] == asse and r[3] == var and tuple(r[:3]) != ANOMALA]))


def main() -> None:
    pt.applica_stile()
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    Kt, seg = par.sperimentale.Kt, par.segnaposto
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    s = json.loads(SOMMARIO.read_text())
    salvati = json.loads((REPORT_CL / "estimated_gains.json").read_text())
    mot = MotorParams.from_yaml(RADICE / "closedLoopModel" / "data" / "motor_TBM7646A.yaml", "cold")
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    r_ol, r_cl, guadagni, limiti = [], [], {}, {}
    for asse, chiave, nome, cat, mod in (("azimut", "azimuth", "Azimuth", cat_az, fl.azimut),
                                         ("polare", "elevation", "Elevation", cat_pol, fl.polare)):
        akd, wb = aac.DRIVE[asse]
        p_file = DriveParams.from_files(chiave, aac.VCL / akd, aac.VCL / wb, pl_kp_unit="hz")
        finestre = {a: dt.finestre(aac.VCL / "V2" / (aac.PROVE[asse] % a), asse) for a in aac.AMPIEZZE}
        c_ol, c_cl = {}, {}
        for var in VARIANTI:
            v = s[chiave]["report_hybrid"]["parameters"] if var == "report" else s[chiave][var]["hybrid"]["parameters"]
            m = mod(cat, v, Kt, {}).struttura
            att = dt.Attrito(v["sigma0"], v["Tc"])
            if var == "report":
                g = salvati[asse]["estimated"]
            else:
                g, _ = aac.stima_guadagni(m, dataclasses.replace(p_file, **aac.SETTEMBRE[asse]), mot, att,
                                          finestre[aac.STIMA])
            guadagni.setdefault(chiave, {})[var] = g
            print(f"{asse} {var}: guadagni {g}", flush=True)
            drive = Drive(dataclasses.replace(p_file, **g))
            an = aac.AnelloChiuso(m, drive, mot, att)
            sim = dt.Simulatore(m)
            c_ol[var], c_cl[var] = {}, {}
            picchi = {"i_cmd": 0.0, "v_phase": 0.0, "i_int": 0.0}
            for a, ww in finestre.items():
                for w in ww:
                    vc = np.abs(np.gradient(w.cmd, dt.DT)) > 1e-6
                    i0, i1 = np.flatnonzero(vc)[[0, -1]]
                    # anello aperto: fase di movimento fino a 0.5 s dopo la fine del comando
                    y = sim.corri(w.I, att)
                    mv = slice(i0, min(i1 + int(0.5 / dt.DT), len(y)))
                    r_ol.append([asse, a, w.nome[-2:], var, dt.fit(y, w.pos), dt.fit(y[mv], w.pos[mv])])
                    c_ol[var][(a, w.nome[-2:])] = (w, y, i0)
                    # anello chiuso: fino a 1 s dopo la fine del comando
                    t = np.arange(0.0, w.t[-1], aac.PASSO)
                    r = an.corri(np.interp(t, w.t, w.cmd - w.cmd[0]), float(w.I[:int(0.4 / dt.DT)].mean()))
                    for k in picchi:
                        picchi[k] = max(picchi[k], float(np.abs(r[k]).max()))
                    y = np.interp(w.t, t, r["p"])
                    im = np.interp(w.t, t, r["i"])
                    mv = slice(i0, min(i1 + int(1.0 / dt.DT), len(y)))
                    r_cl.append([asse, a, w.nome[-2:], var, dt.fit(y, w.pos), dt.fit(y[mv], w.pos[mv]),
                                 dt.fit(im[mv], w.I[mv]), np.rad2deg(np.max(np.abs(y - (w.cmd - w.cmd[0]))[mv])) * 1e3,
                                 float(np.abs(im[mv]).max())])
                    c_cl[var][(a, w.nome[-2:])] = (w, y, im, i0)
            limiti.setdefault(chiave, {})[var] = picchi | {"i_limit": drive.p.i_limit,
                                                           "vel_int_limit": drive.p.vel_int_limit,
                                                           "V_phase_max": drive.V_phase_max}
        pt.salva(figura_aperto(nome, c_ol, seg), CARTELLA / "figures", f"{chiave}_open_loop_resampled")
        pt.salva(figura_chiuso(nome, c_cl, seg), CARTELLA / "figures", f"{chiave}_closed_loop_resampled")
        plt.close("all")

    with open(CARTELLA / "open_loop_fit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["axis", "test", "move", "parameters", "fit_window_pct", "fit_motion_pct"])
        wr.writerows([r[:4] + [f"{x:.1f}" for x in r[4:]] for r in r_ol])
    with open(CARTELLA / "closed_loop_fit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["axis", "test", "move", "parameters", "fit_position_window_pct", "fit_position_motion_pct",
                     "fit_current_motion_pct", "max_following_error_model_mdeg", "max_current_model_A"])
        wr.writerows([r[:4] + [f"{x:.2f}" for x in r[4:]] for r in r_cl])
    mediane = {asse: {var: {"open_loop_fit_motion": mediana(r_ol, asse, var, 5),
                            "closed_loop_fit_position_motion": mediana(r_cl, asse, var, 5),
                            "closed_loop_fit_current_motion": mediana(r_cl, asse, var, 6)} for var in VARIANTI}
               for asse in ("azimut", "polare")}
    controllo = {"open_loop_max_fit_difference": controllo_report(r_ol, REPORT_OL, "fit_motion_pct", 5),
                 "closed_loop_max_fit_difference": controllo_report(r_cl, REPORT_CL / "closed_loop_fit.csv",
                                                                    "fit_position_motion_pct", 5)}
    print(json.dumps({"medians_without_anomalous_window": mediane, "report_reproduction": controllo}, indent=1))
    print("limiti:", json.dumps(limiti, indent=1))
    (CARTELLA / "time_domain_summary.json").write_text(json.dumps(
        {"gains": guadagni, "limits": limiti, "medians_without_anomalous_window": mediane,
         "report_reproduction": controllo}, indent=1))


if __name__ == "__main__":
    main()
