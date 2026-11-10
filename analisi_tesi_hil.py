#!/usr/bin/env python3
r"""Scenari definitivi per la tesi (capitolo "Hardware in the Loop testing: Simulated Pointer - Real Drone").

Puntatore simulato (plant non lineare accoppiato con cuscinetti di Dahl, come in ``analisi_tesi_puntamento``) con la cascata tarata
sui doppietti e l'LQG-PI tarato. Traiettorie: circonferenza di Attorrese e traiettorie registrate del drone (dataset7 e dataset8,
riferimento dagli stimatori EKF2 e OpenVINS ``config_v3``, prima run; errore calcolato sulla ground truth TagSLAM).

Scenari:
    0  nominale: plant nominale, nessun disturbo;
    1  coppia di disturbo gaussiana a media nulla su entrambi gli assi, banda 0.5-5 Hz (equivalente generico delle raffiche);
    2  attrito di Coulomb del plant x1.5 su entrambi gli assi (errore di modello: i controllori non cambiano);
    3  1 e 2 insieme.

Figure nello stile del capitolo della campagna sperimentale (``style.py`` della tesi: PNG a 200 dpi, 15 cm, Computer Modern)
in ``Thesis-Writing/plots/hil/``; tabelle LaTeX in ``results_lqg_tuning/no_bearings/thesis_hil/tables/``.
Richiesta dell'utente (2026-11-03): la circonferenza va nel capitolo LQG (sezione dei risultati simulati, scenari 0-2:
figura ``plots/lqg/lqg_circle_miss.png`` e tabella ``lqg_circle.tex``), le traiettorie registrate restano nel capitolo HIL.
Errore di inseguimento in elevazione contro quello in azimut sulla circonferenza (come Gawronski, Fig. 8.14):
``plots/lqg/lqg_circle_servo_error.png``; in ``numbers.json`` il puntamento efficace con la sola componente in azimut o
polare dell'effetto della coppia (chiavi ``circle|<cod>|1|eta_az_only`` ed ``eta_el_only``), per spiegare perche' la coppia a
media nulla alza eta dell'LQG-PI e abbassa quello della cascata.

    .venv/bin/python -I analisi_tesi_hil.py             simulazioni (LAVORI processi), poi figure e tabelle
    .venv/bin/python -I analisi_tesi_hil.py --prova     un caso per scenario sulla circonferenza
    .venv/bin/python -I analisi_tesi_hil.py --figure    solo figure e tabelle
"""

from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
SOLO_FIGURE = "--figure" in sys.argv
PROVA = "--prova" in sys.argv
sys.argv = [sys.argv[0]] + (["--figure"] if SOLO_FIGURE else [])
sys.path[:0] = [str(RADICE)]

import analisi_tesi_puntamento as tp  # noqa: E402  (senza --figure costruisce assi e controllori in tp.G)

cc, g, apl = tp.cc, tp.g, tp.apl

CARTELLA = RADICE / "results_lqg_tuning" / "no_bearings" / "thesis_hil"
FIGURE = Path("/home/alessandro/Desktop/Thesis/Thesis-Writing/plots/hil")
FIGURE_LQG = FIGURE.parent / "lqg"
SCENARI_CERCHIO = (0, 1, 2)
"""Scenari della circonferenza nel capitolo LQG: nominale, coppia, attrito."""
STILE = Path("/home/alessandro/Desktop/Thesis/OpenVins/OpenVincenzo/colcon_ws/src/open_vins/ov_data/sc_cantina01_thesis/output/plots/scripts")
"""``style.py`` del capitolo della campagna sperimentale: stessi font, dimensioni, formato e colori degli stimatori."""

SIGMA_TW, BANDA = 0.3, (0.5, 5.0)
"""Scenario 1: rms per asse [N m] e banda [Hz]. I picchi superano l'attrito di Coulomb dei due assi (0.23 e 0.41 N m), dove il guadagno d'anello conta."""
F_ATTRITO = 1.5
"""Scenario 2: livello di Coulomb del plant moltiplicato; in azimut e' il valore piu' alto delle strutture provate nell'identificazione (0.35 N m)."""
SCENARI = {0: {}, 1: {"gauss": SIGMA_TW, "banda": BANDA}, 2: {"attrito": F_ATTRITO}, 3: {"gauss": SIGMA_TW, "banda": BANDA, "attrito": F_ATTRITO}}
NOMI_SCENARI = {0: "Nominal", 1: "Torque", 2: "Friction", 3: "Combined"}
STOCASTICI = (1, 3)
SEMI = (10, 20, 30)
"""Tre realizzazioni della coppia, le stesse per tutte le traiettorie e i controllori (seme + 1 sull'asse polare)."""
TRAIETTORIE = ["circle", ("dataset7", "config_v3"), ("dataset8", "config_v3"), ("dataset7", "EKF2"), ("dataset8", "EKF2")]
CONTROLLORI = {"cascade": "cascade", "LQG-PI tuned": "lqg"}


def nome(tag: str, cod: str, scen: int, seme) -> str:
    return f"{tag}__{cod}__s{scen}" + ("" if seme is None else f"_{seme}")


def un_caso(arg) -> dict:
    """Una simulazione: (traiettoria, controllore, scenario, seme). Metriche come ``tp.un_caso``."""
    traj, et, scen, seme = arg
    G = tp.G
    C = {q: G["Cc"][q][et] for q in G["Ax"]}
    if traj == "circle":
        t = np.arange(0.0, cc.T_FINE, g.PASSO_SIM)
        ref = apl.circonferenza(t)
        vero = lambda x: np.stack([apl.CENTRO + apl.RAGGIO * np.cos(apl.OMEGA * x), apl.RAGGIO * np.sin(apl.OMEGA * x), np.full_like(x, apl.QUOTA)], axis=1)
        tag = "circle"
    else:
        ds, st = traj
        tr = G["TR"][ds]
        s = tr["stime"][st]
        t = np.arange(0.0, tr["durata"], g.PASSO_SIM)
        psi, th = tp.angoli(s["p"])
        j = np.clip(np.searchsorted(s["t"], t, side="right") - 1, 0, len(psi) - 1)
        ref = {"azimut": psi[j], "polare": th[j]}
        vero = tr["vero"]
        tag = f"{ds}_{st}"
    R = cc.simula_coppia(G["Ax"], C, ref, t, SCENARI[scen], G["righe"], cc.SEME if seme is None else seme, G["cl"], float(ref["polare"][0]))
    k = np.arange(0, len(t), tp.DECIMA)
    tk = t[k]
    pv = vero(tk)
    psi_l, th_l = ref["azimut"][0] + R["azimut"]["los"][k], ref["polare"][0] + R["polare"]["los"][k]
    E, miss = tp.errore(tp.versore(psi_l, th_l), pv)
    Es, ms = tp.errore(tp.versore(ref["azimut"][k], ref["polare"][k]), pv)
    Et = np.linalg.norm(np.cross(tp.versore(psi_l, th_l), tp.versore(ref["azimut"][k], ref["polare"][k])), axis=1)
    m = tk >= tp.T_TRANS
    cod = CONTROLLORI[et]
    q = {"traj": tag, "controller": cod, "scenario": scen, "seed": seme} | tp.statistiche(tk, E, miss)
    q |= {"estimator_alone": tp.statistiche(tk, Es, ms), "loop_E_rms_mrad": 1e3 * float(np.sqrt(np.mean(Et[m] ** 2))), "loop_E_max_mrad": 1e3 * float(Et[m].max()), "duration_s": float(tk[-1])}
    q |= {f"i_max_{a}_A": float(np.abs(R[b]["i"][k][m]).max()) for a, b in (("az", "azimut"), ("el", "polare"))}
    q |= {f"i_rms_{a}_A": float(np.sqrt(np.mean(R[b]["i"][k][m] ** 2))) for a, b in (("az", "azimut"), ("el", "polare"))}
    q |= {f"los_{a}_rms_mrad": 1e3 * float(np.sqrt(np.mean((ref[b][k] - ref[b][0] - R[b]["los"][k])[m] ** 2))) for a, b in (("az", "azimut"), ("el", "polare"))}
    np.savez_compressed(CARTELLA / "series" / f"{nome(tag, cod, scen, seme)}.npz", t=tk.astype(np.float32), E=E.astype(np.float32), miss=miss.astype(np.float32), miss_est=ms.astype(np.float32),
                        los_az=R["azimut"]["los"][k].astype(np.float32), los_el=R["polare"]["los"][k].astype(np.float32), ref_az=(ref["azimut"][k] - ref["azimut"][0]).astype(np.float32),
                        ref_el=(ref["polare"][k] - ref["polare"][0]).astype(np.float32), i_az=R["azimut"]["i"][k].astype(np.float32), i_el=R["polare"]["i"][k].astype(np.float32))
    print(nome(tag, cod, scen, seme), "| E rms %.3f mrad | miss rms %.2f cm | on target %.1f %% | i max %.1f / %.1f A" % (q["E_rms_mrad"], q["miss_rms_cm"], q["on_target_pct"], q["i_max_az_A"], q["i_max_el_A"]), flush=True)
    return q


# ----------------------------------------------------------------------------------------------- figure (stile della tesi)
def stile():
    sys.path.insert(0, str(STILE))
    import style  # noqa: E402
    style.setup()
    return style


CTRL = {"cascade": dict(label="Cascade", color="#BF7168", ls="--", marker="s", lw=1.3),   # mattone pastello deciso (prima #CC79A7)
        "lqg": dict(label="LQG-PI", color="#5E8A71", ls="-", marker="o", lw=1.3)}       # verde bosco pastello (prima #009E73)
"""Colori dei controllori (Okabe-Ito), distinti dal blu di OpenVINS e dall'arancio di EKF2 usati nel resto della tesi; il tratto differisce sempre."""
IDEALE = dict(label="Perfect controller (upper bound)", color="#777777", ls=":", lw=1.1)
FASCIO = dict(color="#000000", ls=(0, (4, 2)), lw=0.8)
TAG = {"circle": "circle", "dataset7_config_v3": "dataset7, OpenVINS config v3", "dataset8_config_v3": "dataset8, OpenVINS config v3", "dataset7_EKF2": "dataset7, EKF2", "dataset8_EKF2": "dataset8, EKF2"}


def serie(tag, cod, scen, seme=None) -> dict:
    return dict(np.load(CARTELLA / "series" / f"{nome(tag, cod, scen, seme)}.npz"))


def tag_pannello(st, a_, testo, x=0.01, y=0.97):
    a_.text(x, y, testo, transform=a_.transAxes, va="top", fontsize=st.FONT_SIZE - 1, bbox=dict(facecolor="white", edgecolor="none", pad=1))


def legenda(fig, h, l, ncol):
    fig.legend(h, l, loc="lower center", ncol=ncol, bbox_to_anchor=(0.5, 1.0))


def fig_traiettorie(st) -> None:
    """Vista dall'alto delle traiettorie registrate, con la direzione di Lucy."""
    import matplotlib.pyplot as plt
    fig, ax = st.figure(aspect=0.5, ncols=2)
    vero = {**st.line_kw("reference"), "lw": 1.0, "label": "True position (TagSLAM)"}
    ax[0].set_ylabel("Distance along the line of sight [m]")
    for a_, ds in zip(ax, ("dataset7", "dataset8")):
        tr = tp.traiettoria(ds)
        v = tr["vero"](np.arange(0.0, tr["durata"], 0.02))
        p0 = v[0]
        a_.plot(v[:, 1] - p0[1], v[:, 0] - p0[0], **vero)
        for alg, kw in (("ekf2", "EKF2"), ("openvins_dv2", "config_v3")):
            p = tr["stime"][kw]["p"]
            a_.plot(p[:, 1] - p0[1], p[:, 0] - p0[0], **{**st.line_kw(alg), "lw": 0.9})
        a_.plot(0, 0, "o", color="#000000", ms=3.5)
        tag_pannello(st, a_, ds)
    for a_ in ax:
        a_.set_xlim(-4.8, 1.8)
        a_.set_ylim(-4.6, 2.4)
    for a_ in ax:
        a_.set_aspect("equal", adjustable="box")
        a_.set_xlabel("Lateral position [m]")
    # solo la direzione di Lucy: la freccia esce dal bordo inferiore, Lucy e' 100 m piu' in basso, fuori scala
    for a_ in ax:
        a_.annotate("", xy=(0.12, -0.005), xycoords="axes fraction", xytext=(0.12, 0.17), textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8))
        a_.text(0.16, 0.085, "towards Lucy\n(100 m, off scale)", transform=a_.transAxes, ha="left", va="center", fontsize=st.FONT_SIZE - 2, color="#555555")
    h, l = ax[0].get_legend_handles_labels()
    legenda(fig, h, l, 3)
    fig.tight_layout()
    st.save(fig, FIGURE / "hil_trajectories.png")
    plt.close(fig)


def fig_disturbo(st) -> None:
    """Scenario 1: una realizzazione della coppia nel tempo e densita' spettrale di ampiezza."""
    from scipy.signal import welch
    t = np.arange(0.0, 120.0, g.PASSO_SIM)
    sc = SCENARI[1]
    d = {n: cc.disturbo(sc, t, SEMI[0] + i) for i, n in enumerate(("azimuth", "polar"))}
    kw = {"azimuth": dict(color="#333333", ls="-", lw=0.9, label="Azimuth"), "polar": dict(color="#8C8C8C", ls="--", lw=0.9, label="Polar")}
    fig, ax = st.figure(aspect=0.36, ncols=2, gridspec_kw={"width_ratios": [1.4, 1]})
    k = (t >= 30.0) & (t <= 34.0)
    for n, x in d.items():
        ax[0].plot(t[k] - 30.0, x[k], **kw[n])
        f, P = welch(x, fs=1 / g.PASSO_SIM, nperseg=int(20.0 / g.PASSO_SIM))
        ax[1].loglog(f[1:], np.sqrt(P[1:]), **{**kw[n], "label": None})
    for tc, ls, lab in ((0.2316, ":", "Coulomb level, azimuth"), (0.4073, "-.", "Coulomb level, polar")):
        for s in (1, -1):
            ax[0].axhline(s * tc, color="#B0B0B0", lw=0.7, ls=ls, label=lab if s == 1 else None, zorder=1)
    ax[0].set_xlabel("Time [s]")
    ax[0].set_ylabel(r"Disturbance torque $\tau_w$ [N m]")
    # sotto 0.1 Hz la stima di Welch e' dominata dalla dispersione dei segmenti, non dal segnale (il filtro attenua oltre 100 dB)
    ax[1].set_xlim(0.1, 50)
    ax[1].set_ylim(1e-4, 1)
    for fb in BANDA:
        ax[1].axvline(fb, color="#999999", lw=0.6, ls=":")
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylabel(r"ASD of $\tau_w$ [N m/$\sqrt{\mathrm{Hz}}$]")
    h, l = ax[0].get_legend_handles_labels()
    legenda(fig, h, l, 4)
    fig.tight_layout()
    st.save(fig, FIGURE / "hil_disturbance.png")


def raccogli(S: list, tag: str, cod: str, scen: int, chiave: str) -> list:
    return [y[chiave] for y in S if y["traj"] == tag and y["controller"] == cod and y["scenario"] == scen]


def fig_efficace(st, S: list) -> None:
    """Puntamento efficace per scenario: media sulle realizzazioni, barre d'errore sull'intervallo."""
    pannelli = ("dataset7_config_v3", "dataset8_config_v3")
    fig, ax = st.figure(aspect=0.45, ncols=2, sharey=True)
    x = np.arange(len(SCENARI))
    w = 0.36
    for a_, tag in zip(ax, pannelli):
        for i, cod in enumerate(("cascade", "lqg")):
            v = [raccogli(S, tag, cod, sc, "on_target_pct") for sc in SCENARI]
            m = np.array([np.mean(y) for y in v])
            err = np.array([[np.mean(y) - min(y) for y in v], [max(y) - np.mean(y) for y in v]])
            s = CTRL[cod]
            a_.bar(x + (i - 0.5) * w, m, width=w * 0.95, color=s["color"], edgecolor="#777777", linewidth=0.2, label=s["label"])
            a_.errorbar(x + (i - 0.5) * w, m, yerr=err, fmt="none", ecolor="#555555", elinewidth=0.4, capsize=1.2, capthick=0.4)
        ide = next(y for y in S if y["traj"] == tag and y["scenario"] == 0)["estimator_alone"]["on_target_pct"]
        a_.axhline(ide, **{k: IDEALE[k] for k in ("color", "ls", "lw")}, label=IDEALE["label"])
        a_.set_xticks(x)
        a_.set_xticklabels([NOMI_SCENARI[s] for s in SCENARI], rotation=30, ha="right")
        a_.grid(axis="x", visible=False)
        a_.set_ylim(0, 124)
        a_.set_yticks(range(0, 101, 20))
        tag_pannello(st, a_, TAG[tag].replace(", OpenVINS ", "\n"))
    ax[0].set_ylabel("Effective pointing [%]")
    h, l = ax[0].get_legend_handles_labels()
    legenda(fig, h, l, 3)
    fig.tight_layout()
    st.save(fig, FIGURE / "hil_effective_pointing.png")


def fig_mancato(st) -> None:
    """Distanza di mancato puntamento nel tempo: dataset8 con config_v3, scenari 0 e 3 (prima realizzazione)."""
    righe = (("dataset8_config_v3", "dataset8, OpenVINS config v3", 38),)
    fig, ax = st.figure(aspect=0.36, nrows=1, ncols=2, sharey="row", squeeze=False)
    for i, (tag, lab, ymax) in enumerate(righe):
        for j, (scen, seme) in enumerate(((0, None), (3, SEMI[0]))):
            a_ = ax[i, j]
            ogni = 20
            x = serie(tag, "cascade", scen, seme)
            a_.plot(x["t"][::ogni], 100 * x["miss_est"][::ogni], **IDEALE, zorder=1)
            for cod in ("cascade", "lqg"):
                x = serie(tag, cod, scen, seme)
                s = CTRL[cod]
                a_.plot(x["t"][::ogni], 100 * x["miss"][::ogni], color=s["color"], ls=s["ls"], lw=0.8 if scen == 0 else 0.6, label=s["label"], zorder=2 + (cod == "lqg"))
            a_.axhline(100 * tp.R_FASCIO, **FASCIO, label="Beam radius, 8 cm", zorder=4)
            a_.axvspan(0, tp.T_TRANS, color="#EEEEEE", zorder=0, lw=0)
            a_.set_xlim(0, float(x["t"][-1]))
            a_.set_ylim(0, ymax)
            tag_pannello(st, a_, f"{lab}, {NOMI_SCENARI[scen].lower()}")
        ax[i, 0].set_ylabel("Miss distance [cm]")
    for a_ in ax[-1]:
        a_.set_xlabel("Time [s]")
    h, l = ax[-1, 0].get_legend_handles_labels()
    legenda(fig, h, l, 4)
    fig.tight_layout()
    st.save(fig, FIGURE / "hil_miss_distance.png")


def errori_servo(tag: str, cod: str, scen: int, seme=None) -> tuple:
    """Errore di inseguimento in azimut e in elevazione [rad], distanza dal drone [m], dopo il transitorio."""
    x = serie(tag, cod, scen, seme)
    m = x["t"] >= tp.T_TRANS
    r = x["miss"][m] / np.maximum(x["E"][m], 1e-12)
    return (x["ref_az"] - x["los_az"])[m], (x["ref_el"] - x["los_el"])[m], r


def eta_componenti(cod: str) -> dict:
    """Puntamento efficace sulla circonferenza con la sola componente in azimut o polare dell'effetto della coppia
    (scenario 1 meno scenario 0, punto per punto), media sulle realizzazioni."""
    az0, el0, r = errori_servo("circle", cod, 0)
    c = np.cos(np.arctan2(tp.QUOTA, r))
    out = {"eta_az_only": [], "eta_el_only": []}
    for s in SEMI:
        az1, el1, _ = errori_servo("circle", cod, 1, s)
        out["eta_az_only"].append(100 * np.mean(r * np.hypot(c * az1, el0) <= tp.R_FASCIO))
        out["eta_el_only"].append(100 * np.mean(r * np.hypot(c * az0, el1) <= tp.R_FASCIO))
    return {k: float(np.mean(v)) for k, v in out.items()}


def fig_servo_cerchio(st) -> None:
    """Capitolo LQG: errore in elevazione contro errore in azimut sulla circonferenza (Gawronski, Fig. 8.14),
    nominale e con la prima realizzazione della coppia, con il fascio a 100 m."""
    fig, ax = st.figure(aspect=0.5, ncols=2, sharex=True, sharey=True)
    ph = np.linspace(0, 2 * np.pi, 200)
    rb = 1e3 * tp.R_FASCIO / tp.D_LUCY
    for a_, cod, tag in zip(ax, ("cascade", "lqg"), ("(a) Cascade", "(b) LQG-PI")):
        az1, el1, _ = errori_servo("circle", cod, 1, SEMI[0])
        a_.plot(1e3 * az1, 1e3 * el1, ".", ms=1.0, color=CTRL[cod]["color"], alpha=0.5, label="Disturbance torque")
        az0, el0, _ = errori_servo("circle", cod, 0)
        a_.plot(1e3 * az0, 1e3 * el0, "-", lw=1.0, color="#000000", label="Nominal")
        a_.plot(rb * np.cos(ph), rb * np.sin(ph), **FASCIO, label="Beam radius, 8 cm at 100 m")
        a_.set_aspect("equal")
        a_.set_xlim(-3.2, 3.2)
        a_.set_ylim(-2.6, 2.6)
        a_.set_title(tag, loc="left", fontsize=st.FONT_SIZE)
        a_.set_xlabel("Azimuth servo error [mrad]")
    ax[0].set_ylabel("Elevation servo error [mrad]")
    h, l = [], []
    for a_ in ax:
        for hh, ll in zip(*a_.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh), l.append(ll)
    # nella legenda i punti del disturbo come marcatori leggibili
    h[0] = h[0].__class__([], [], ls="none", marker="o", ms=3, color=CTRL["cascade"]["color"])
    h.insert(1, h[0].__class__([], [], ls="none", marker="o", ms=3, color=CTRL["lqg"]["color"]))
    l = ["Torque, cascade", "Torque, LQG-PI"] + l[1:]
    legenda(fig, h, l, 4)
    fig.tight_layout()
    FIGURE_LQG.mkdir(parents=True, exist_ok=True)
    st.save(fig, FIGURE_LQG / "lqg_circle_servo_error.png")


def fig_cerchio(st) -> None:
    """Capitolo LQG: distanza di mancato puntamento sulla circonferenza, scenari 0-2 (coppia: prima realizzazione)."""
    fig, ax = st.figure(aspect=0.95, nrows=3, ncols=1, sharex=True, sharey=True)
    for a_, scen in zip(ax, SCENARI_CERCHIO):
        seme = SEMI[0] if scen in STOCASTICI else None
        for cod in ("cascade", "lqg"):
            x = serie("circle", cod, scen, seme)
            s = CTRL[cod]
            a_.plot(x["t"][::10], 100 * x["miss"][::10], color=s["color"], ls=s["ls"], lw=0.8, label=s["label"], zorder=2 + (cod == "lqg"))
        a_.axhline(100 * tp.R_FASCIO, **FASCIO, label="Beam radius, 8 cm", zorder=4)
        a_.axvspan(0, tp.T_TRANS, color="#EEEEEE", zorder=0, lw=0)
        a_.set_xlim(0, float(x["t"][-1]))
        a_.set_ylim(0, 44)
        a_.set_ylabel("Miss distance [cm]")
        tag_pannello(st, a_, "(%s) %s" % ("abc"[scen], ("Nominal", "Disturbance torque", "Higher friction")[scen]))
    ax[-1].set_xlabel("Time [s]")
    h, l = ax[0].get_legend_handles_labels()
    legenda(fig, h, l, 3)
    fig.tight_layout()
    FIGURE_LQG.mkdir(parents=True, exist_ok=True)
    st.save(fig, FIGURE_LQG / "lqg_circle_miss.png")


# ----------------------------------------------------------------------------------------------- tabelle
def effetto(tag: str, cod: str, scen: int) -> tuple[float, float]:
    """Spostamento della linea di vista rispetto allo scenario 0 dello stesso caso, rms per asse [urad], medio sulle realizzazioni."""
    b = serie(tag, cod, 0)
    m = b["t"] >= tp.T_TRANS
    semi = SEMI if scen in STOCASTICI else (None,)
    az, el = [], []
    for s in semi:
        x = serie(tag, cod, scen, s)
        az.append(1e6 * np.sqrt(np.mean((x["los_az"] - b["los_az"])[m] ** 2)))
        el.append(1e6 * np.sqrt(np.mean((x["los_el"] - b["los_el"])[m] ** 2)))
    return float(np.mean(az)), float(np.mean(el))


def tabelle(S: list) -> dict:
    (CARTELLA / "tables").mkdir(parents=True, exist_ok=True)
    etichette = {"dataset7_config_v3": ("dataset7", r"\texttt{config\_v3}"), "dataset8_config_v3": ("dataset8", r"\texttt{config\_v3}"),
                 "dataset7_EKF2": ("dataset7", "EKF2"), "dataset8_EKF2": ("dataset8", "EKF2")}
    nc = {"cascade": "cascade", "lqg": "LQG-PI"}
    eta, err, ax_, num = [], [], [], {}
    for n_, (tag, (tr_, est)) in enumerate(etichette.items()):
        if n_:
            for lst in (eta, err, ax_):
                lst.append(r"\addlinespace")
        ide = next(y for y in S if y["traj"] == tag and y["scenario"] == 0)["estimator_alone"]
        for i, cod in enumerate(("cascade", "lqg")):
            # dataset, stimatore e limite del controllore perfetto solo sulla riga della cascata
            prim = [tr_, est, f"{ide['on_target_pct']:.1f}"] if i == 0 else ["", "", ""]
            c_e = []
            for sc in SCENARI:
                v = raccogli(S, tag, cod, sc, "on_target_pct")
                c_e.append(f"{v[0]:.1f}" if len(v) == 1 else f"{np.mean(v):.1f} ({min(v):.1f}--{max(v):.1f})")
                num[f"{tag}|{cod}|{sc}|eta"] = float(np.mean(v))
            eta.append(" & ".join(prim + [nc[cod]] + c_e) + r" \\")
            c_r = [f"{np.mean(raccogli(S, tag, cod, sc, 'E_rms_mrad')):.2f}" for sc in SCENARI]
            c_m = [f"{max(raccogli(S, tag, cod, sc, 'E_max_mrad')):.2f}" for sc in (0, 3)]
            i_az = max(max(raccogli(S, tag, cod, sc, "i_max_az_A")) for sc in SCENARI)
            i_el = max(max(raccogli(S, tag, cod, sc, "i_max_el_A")) for sc in SCENARI)
            err.append(" & ".join(prim[:2] + [nc[cod]] + c_r + c_m + [f"{i_az:.1f}", f"{i_el:.1f}"]) + r" \\")
            c_a = []
            for sc in (1, 2, 3):
                a, e = effetto(tag, cod, sc)
                c_a += [f"{a:.0f}", f"{e:.0f}"]
                num[f"{tag}|{cod}|{sc}|effect"] = (a, e)
            ax_.append(" & ".join(prim[:2] + [nc[cod]] + c_a) + r" \\")
            for sc in SCENARI:
                for k in ("E_rms_mrad", "miss_rms_cm", "loop_E_rms_mrad", "los_az_rms_mrad", "los_el_rms_mrad", "i_max_az_A", "i_max_el_A"):
                    num[f"{tag}|{cod}|{sc}|{k}"] = float(np.mean(raccogli(S, tag, cod, sc, k)))
    cer = []
    for i, cod in enumerate(("cascade", "lqg")):
        for sc in SCENARI_CERCHIO:
            r = lambda k: raccogli(S, "circle", cod, sc, k)  # noqa: E731
            eta_c = r("on_target_pct")
            c_eta = f"{eta_c[0]:.1f}" if len(eta_c) == 1 else f"{np.mean(eta_c):.1f} ({min(eta_c):.1f}--{max(eta_c):.1f})"
            a, e = effetto("circle", cod, sc) if sc else (None, None)
            eff = ["--", "--"] if sc == 0 else [f"{a:.0f}", f"{e:.0f}"]
            cer.append(" & ".join([nc[cod] if sc == 0 else "", NOMI_SCENARI[sc], f"{np.mean(r('E_rms_mrad')):.2f}", f"{max(r('E_max_mrad')):.2f}", c_eta] + eff
                                  + [f"{max(r('i_max_az_A')):.1f}", f"{max(r('i_max_el_A')):.1f}"]) + r" \\")
            for k in ("E_rms_mrad", "E_max_mrad", "on_target_pct", "i_max_az_A", "i_max_el_A", "los_az_rms_mrad", "los_el_rms_mrad"):
                num[f"circle|{cod}|{sc}|{k}"] = float(np.mean(r(k)))
            num[f"circle|{cod}|{sc}|E_max_max_mrad"] = float(max(r("E_max_mrad")))
            if sc:
                num[f"circle|{cod}|{sc}|effect"] = (a, e)
            if sc == 1:
                for k, v in eta_componenti(cod).items():
                    num[f"circle|{cod}|1|{k}"] = v
        if i == 0:
            cer.append(r"\addlinespace")
    for n_, lst in (("hil_eta", eta), ("hil_error", err), ("hil_axis", ax_), ("lqg_circle", cer)):
        (CARTELLA / "tables" / f"{n_}.tex").write_text("\n".join(lst) + "\n")
    (CARTELLA / "numbers.json").write_text(json.dumps(num, indent=1))
    return num


def main() -> None:
    (CARTELLA / "series").mkdir(parents=True, exist_ok=True)
    FIGURE.mkdir(parents=True, exist_ok=True)
    if PROVA:
        for sc in SCENARI:
            un_caso(("circle", "LQG-PI tuned", sc, SEMI[0] if sc in STOCASTICI else None))
        return
    if not SOLO_FIGURE:
        compiti = [(tr, et, sc, s) for tr in TRAIETTORIE[::-1] for et in CONTROLLORI for sc in SCENARI for s in (SEMI if sc in STOCASTICI else (None,))]
        with ProcessPoolExecutor(int(os.environ.get("LAVORI", "3"))) as ex:
            righe = list(ex.map(un_caso, compiti))
        geo = {ds: {"origin_m": tp.G["TR"][ds]["origine"].tolist(), "window_s": tp.G["TR"][ds]["durata"],
                    "estimators": {n: {k: s[k] for k in ("ate_m", "yaw_deg", "rate_Hz")} for n, s in tp.G["TR"][ds]["stime"].items()}} for ds in tp.DATASET}
        (CARTELLA / "summary.json").write_text(json.dumps({"cases": righe, "geometry": geo, "scenarios": {str(k): {kk: (list(vv) if isinstance(vv, tuple) else vv) for kk, vv in v.items()} for k, v in SCENARI.items()},
                                                           "seeds": SEMI, "beam_radius_m": tp.R_FASCIO, "distance_m": tp.D_LUCY, "height_m": tp.QUOTA}, indent=1))
    S = json.loads((CARTELLA / "summary.json").read_text())["cases"]
    st = stile()
    fig_traiettorie(st)
    fig_disturbo(st)
    fig_efficace(st, S)
    fig_mancato(st)
    fig_cerchio(st)
    fig_servo_cerchio(st)
    tabelle(S)


if __name__ == "__main__":
    main()
