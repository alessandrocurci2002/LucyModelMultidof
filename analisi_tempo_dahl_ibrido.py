"""Verifica nel tempo del modello con i cuscinetti di Dahl ibridi (T_c dai doppietti, sigma0 dalle FRF).

Richiesta dell'utente (2026-10-31): confronto nel dominio del tempo con gli ultimi parametri, cioe' la stima ibrida di
``analisi_cuscinetti_dahl_tempo.py`` (``results_fork_modes_dahl/time_domain_bearings/summary.json``, sezione
``hybrid``): struttura e sigma0 dalle FRF, T_c dai doppietti V2, Dahl puro (T_s = T_c), g_Kt del modello. Ingresso: la
corrente misurata dei doppietti V2 (AMPS = 2 I_rms); uscita: l'encoder, confrontato con la misura. Nessun parametro e'
stimato qui. Simulatore in :mod:`gimbal_mdof.doppietti_tempo`.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_tempo_dahl_ibrido.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))

import matplotlib.pyplot as plt  # noqa: E402
from scipy.signal import savgol_filter  # noqa: E402

import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import doppietti_tempo as dt  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402

IBRIDO = RADICE / "results_fork_modes_dahl" / "time_domain_bearings" / "summary.json"
V2 = RADICE.parent / "LucyIdentification" / "VClosedLoop" / "V2"
CARTELLA = RADICE / "results_fork_modes_dahl" / "time_domain_check"
PROVE = {"azimut": "MPCLucy_doublets_azimut_%s.csv", "polare": "MPCLucy_doublets_elevation_%s.csv"}
AMPIEZZE = ("01deg_minus01deg", "05deg_minus05deg", "1deg_minus1deg", "3deg_minus3deg", "5deg_minus5deg",
            "8deg_minus8deg")
FIGURA = (("01deg_minus01deg", "#1"), ("1deg_minus1deg", "#1"), ("5deg_minus5deg", "#2"))


def main() -> None:
    pt.applica_stile()
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    Kt, seg = par.sperimentale.Kt, par.segnaposto
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    ib = json.loads(IBRIDO.read_text())
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    righe, sommario = [], {}
    for asse, chiave, cat, mod in (("azimut", "azimuth", cat_az, fl.azimut), ("polare", "elevation", cat_pol, fl.polare)):
        h = ib[chiave]["hybrid"]
        v = h["parameters"] | {"Tc": h["Tc_from_time"]}
        m = mod(cat, v, Kt, {}).struttura
        att = dt.Attrito(v["sigma0"], v["Tc"])
        sim = dt.Simulatore(m)
        curve = {}
        for a in AMPIEZZE:
            for w in dt.finestre(V2 / (PROVE[asse] % a), asse):
                y = sim.corri(w.I, att)
                # fase di movimento: dall'inizio del comando a 0.5 s dopo la sua fine
                vc = np.abs(np.gradient(w.cmd, dt.DT)) > 1e-6
                i0, i1 = np.flatnonzero(vc)[[0, -1]]
                mv = slice(i0, min(i1 + int(0.5 / dt.DT), len(y)))
                k3 = min(int((w.t[i0] + 3.0) / dt.DT), len(y) - 1)
                righe.append([asse, a, w.nome[-2:], dt.fit(y, w.pos), dt.fit(y[mv], w.pos[mv]),
                              np.rad2deg(w.pos[i1]) * 1e3, np.rad2deg(y[i1]) * 1e3,
                              np.rad2deg(w.pos[k3]) * 1e3, np.rad2deg(y[k3]) * 1e3])
                curve[(a, w.nome[-2:])] = (w, y, i0)
        sommario[chiave] = {"parameters": v, "g_Kt": v["g_Kt"]}

        fig, ax = plt.subplots(3, 3, figsize=(15, 9), sharex="col", gridspec_kw={"height_ratios": [2, 1.2, 1]})
        for j, key in enumerate(FIGURA):
            w, y, i0 = curve[key]
            k = w.t <= w.t[i0] + 2.5
            t = w.t[k]
            md = lambda q: np.rad2deg(q[k]) * 1e3  # noqa: E731
            ax[0, j].plot(t, md(w.cmd), color=pt.INCHIOSTRO, lw=0.8, ls=":", label="command")
            ax[0, j].plot(t, md(w.pos), color=pt.INCHIOSTRO, lw=1.6, label="measured encoder")
            ax[0, j].plot(t, md(y), lw=1.4, label="model: hybrid Dahl bearings, measured current input", **pt._stile(1))
            ax[1, j].plot(t, savgol_filter(w.pos, 41, 3, deriv=1, delta=dt.DT)[k], color=pt.INCHIOSTRO, lw=1.2)
            ax[1, j].plot(t, np.gradient(y, dt.DT)[k], lw=1.1, **pt._stile(1))
            ax[2, j].plot(t, w.I[k], color=pt.INCHIOSTRO, lw=1.0)
            ax[0, j].set_title(f"{key[0].split('_')[0].replace('01deg', '0.1deg').replace('05deg', '0.5deg')} doublet, "
                               f"move {key[1]}", fontsize=9)
            ax[2, j].set_xlabel("time from window start [s]")
        ax[0, 0].set_ylabel("position [mdeg]"); ax[1, 0].set_ylabel("velocity [rad/s]")
        ax[2, 0].set_ylabel("input current [A rms]")
        ax[0, 0].legend(loc="lower right", fontsize=7)
        nome = "Azimuth" if asse == "azimut" else "Elevation"
        pt._titolo(fig, f"{nome}: hybrid Dahl model (sigma0 {v['sigma0']:.0f} N m/rad from FRF, T_c {v['Tc']:.3f} N m "
                        f"from doublets, g_Kt {v['g_Kt']:.2f}) driven by the measured current", seg)
        pt.salva(fig, CARTELLA / "figures", f"{chiave}_time_domain_hybrid")
        plt.close(fig)

    with open(CARTELLA / "position_fit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["axis", "test", "move", "fit_window_pct", "fit_motion_pct", "end_of_command_meas_mdeg",
                     "end_of_command_model_mdeg", "at_3s_meas_mdeg", "at_3s_model_mdeg"])
        wr.writerows([r[:3] + [f"{x:.1f}" for x in r[3:]] for r in righe])
    for r in righe:
        print(f"{r[0]:7s} {r[1]:18s} {r[2]}: fit finestra {r[3]:8.1f} %, fit movimento {r[4]:7.1f} % | fine comando "
              f"mis {r[5]:8.1f} mod {r[6]:8.1f} | a 3 s mis {r[7]:8.1f} mod {r[8]:8.1f} mdeg")
    (CARTELLA / "summary.json").write_text(json.dumps(sommario, indent=1))


if __name__ == "__main__":
    main()
