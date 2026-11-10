#!/usr/bin/env python3
r"""Errore di puntamento di LUCY per la tesi: cascata e LQG-PI tarato sul plant non lineare accoppiato con cuscinetti di Dahl.

Definizione dell'errore (Attorrese 2026, eq. 5.20-5.23): :math:`E=\lVert\hat n(\psi,\theta)\times\hat r\rVert`, con :math:`\hat r` verso la posizione *vera*
del drone. Il riferimento dell'anello viene dagli stimatori (EKF2 e OpenVINS con ``config_drone_v2``, ``config_v3`` nella tesi, prima run) sui dataset 7 e 8,
tenuto costante fra due stime. Scenari: circonferenza e traiettorie, senza e con coppia di disturbo gaussiana su entrambi gli assi, a media nulla e,
con ``--media``, a media ``MEDIA_TW`` (aggiunge questi casi a ``summary.json``). Con ``--figure`` rifa' solo le figure. Risultati in ``results_lqg_tuning/no_bearings/thesis_pointing/``.
"""

from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]
SOLO_FIGURE = "--figure" in sys.argv
PROVA = "--prova" in sys.argv
MEDIA = "--media" in sys.argv
VELOCITA = "--velocita" in sys.argv
sys.argv = [sys.argv[0], "--senza-cuscinetti"]

import analisi_lqg_circonferenza as cc  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_punto_lavoro as apl  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning" / "no_bearings" / "thesis_pointing"
DATI = Path("/home/alessandro/Desktop/Thesis/OpenVins/OpenVincenzo/colcon_ws/src/open_vins/ov_data/sc_cantina01_thesis/OpenVINS_eval")
DATASET = {"dataset7": "dataset_7", "dataset8": "dataset_8"}
STIMATORI = {"EKF2": "ekf2", "config_v3": "openvins_dv2"}
CONTROLLORI = ("cascade", "LQG-PI tuned")
D_LUCY = 100.0
QUOTA = 5.0
"""Distanza di LUCY dal primo punto della traiettoria e quota del drone rispetto a LUCY, la stessa della circonferenza."""
R_FASCIO = 0.08
SIGMA_TW = 0.1
SEMI = (10, 20, 30)
"""Coppia di disturbo: rms per asse, banda 0.1-50 Hz (``cc.disturbo``), tre realizzazioni indipendenti, le stesse per tutti i controllori e le traiettorie."""
FATTORI_V = (2, 3, 5, 10)
"""Velocita' delle traiettorie reali moltiplicata per questi fattori sullo stesso percorso (``scala``), solo con ``config_v3``; 1 e' il caso base."""
I_LIMITE = 18.0
"""Limite di corrente del drive usato sui doppietti; il modello non satura, quindi i casi oltre il limite si segnalano."""
MEDIA_TW = 0.2
"""Parte statica della coppia (cavi, sbilanciamento), uguale sui due assi e sommata alle stesse realizzazioni: gia' presente all'avvio, come in ``cc.SCENARI``."""
T_TRANS = 5.0
DECIMA = 10


def carica(f: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.loadtxt(f, comments="#")
    return d[:, 0], d[:, 1:4], d[:, 4:8]


def allinea_posyaw(pe: np.ndarray, pg: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Allineamento a 4 gradi di liberta' (posizione e imbardata), come ``align_posyaw`` di ``evaluate_thesis.py`` e ov_eval."""
    me, mg = pe.mean(0), pg.mean(0)
    rc = (pg - mg).T @ (pe - me)
    rc = rc.T
    th = np.arctan2(rc[0, 1] - rc[1, 0], rc[0, 0] + rc[1, 1])
    c, s = np.cos(th), np.sin(th)
    Rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return Rz, mg - Rz @ me, float(th)


def traiettoria(ds: str) -> dict:
    """Vero (``truths_leveled``: ``truths`` ruotato in modo rigido con z in alto) e stime allineate, in un riferimento con LUCY nell'origine, sulla finestra comune."""
    tv, pv, _ = carica(DATI / "truths_leveled" / f"{ds}.txt")
    tov = np.loadtxt(DATI / "algorithms" / "openvins_dv2" / ds / "run_1.txt", comments="#", usecols=0)
    t0, t1 = float(tov[0]), float(min(tov[-1], tv[-1]))
    origine = pv[0] - np.array([np.sqrt(D_LUCY ** 2 - QUOTA ** 2), 0.0, QUOTA])
    interp = lambda t: np.stack([np.interp(t, tv, pv[:, i]) for i in range(3)], axis=1)
    stime = {}
    for n, alg in STIMATORI.items():
        te, pe, _ = carica(DATI / "algorithms" / alg / ds / "run_1.txt")
        k = (te >= t0) & (te <= t1)
        te, pe = te[k], pe[k]
        pg = interp(te)
        Rz, tz, th = allinea_posyaw(pe, pg)
        pa = pe @ Rz.T + tz
        err = pa - pg
        stime[n] = {"t": te - t0, "p": pa - origine, "ate_m": float(np.sqrt(np.mean(np.sum(err ** 2, axis=1)))), "yaw_deg": float(np.degrees(th)), "rate_Hz": float(1 / np.median(np.diff(te)))}
    return {"t0": t0, "durata": t1 - t0, "vero": lambda t: interp(t + t0) - origine, "stime": stime, "origine": origine, "ds": ds}


def angoli(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.arctan2(p[:, 1], p[:, 0]), np.arctan2(p[:, 2], np.hypot(p[:, 0], p[:, 1]))


def versore(psi: np.ndarray, th: np.ndarray) -> np.ndarray:
    """Eq. 5.22."""
    return np.stack([np.cos(th) * np.cos(psi), np.cos(th) * np.sin(psi), np.sin(th)], axis=1)


def errore(n: np.ndarray, pv: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Eq. 5.23 (a 100 m la regolarizzazione della 5.21 non conta) e distanza del drone dalla linea di vista, che e' :math:`\\lVert r\\rVert E`."""
    d = np.linalg.norm(pv, axis=1)
    E = np.linalg.norm(np.cross(n, pv / d[:, None]), axis=1)
    return E, E * d


def scala(tr: dict, st: str, fv: float) -> tuple:
    """Stesso percorso ``fv`` volte piu' veloce. Vero e stima si comprimono nel tempo; lo stimatore resta alla sua frequenza, quindi a ogni suo periodo
    il riferimento e' la stima piu' recente. L'errore della stima resta quello della corsa originale: e' un'ipotesi."""
    s = tr["stime"][st]
    ts = s["t"] / fv
    tick = np.arange(0.0, tr["durata"] / fv, 1.0 / s["rate_Hz"])
    j = np.clip(np.searchsorted(ts, tick, side="right") - 1, 0, len(ts) - 1)
    return tr["durata"] / fv, tick, s["p"][j], (lambda x: tr["vero"](fv * x))


def statistiche(t: np.ndarray, E: np.ndarray, miss: np.ndarray, t0: float = T_TRANS) -> dict:
    m = t >= t0
    return {"E_rms_mrad": 1e3 * float(np.sqrt(np.mean(E[m] ** 2))), "E_mean_mrad": 1e3 * float(E[m].mean()), "E_p95_mrad": 1e3 * float(np.percentile(E[m], 95)), "E_max_mrad": 1e3 * float(E[m].max()),
            "miss_rms_cm": 100 * float(np.sqrt(np.mean(miss[m] ** 2))), "miss_mean_cm": 100 * float(miss[m].mean()), "miss_p95_cm": 100 * float(np.percentile(miss[m], 95)), "miss_max_cm": 100 * float(miss[m].max()),
            "on_target_pct": 100 * float(np.mean(miss[m] <= R_FASCIO))}


def imposta() -> dict:
    par, Kt, S, T, mot, assi = g.carica()
    Fi = json.loads(cc.TUNING.read_text())
    Ax, Cc, righe = {}, {}, {}
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        m0 = next(m for m in a.modi if len(m.indici) == 2 and m.zeta < 0.1)
        Ax[asse], Cc[asse] = a, cc.casi(a, Fi[ch]["final"], m0.f_hz)
        riga = np.zeros(a.cascata.N)
        riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
        righe[asse] = riga
    scala_J = float(S["azimuth"]["uniform_lnf"]["hybrid"]["parameters"]["s_J_L"])
    return {"Ax": Ax, "Cc": Cc, "righe": righe, "cl": cc.CoppieLagrange(afm.sistema(par).payload, scala_J), "TR": {ds: traiettoria(ds) for ds in DATASET}}


if not SOLO_FIGURE:
    G = imposta()


def un_caso(arg) -> tuple:
    """Una simulazione. ``traj``: ``circle`` o (dataset, stimatore). ``seme``: None senza disturbo. ``media``: parte statica della coppia [N m]."""
    traj, et, seme, *resto = arg
    media = resto[0] if resto else 0.0
    sc = {} if seme is None else {"gauss": SIGMA_TW} | ({"media": media} if media else {})
    C = {q: G["Cc"][q][et] for q in G["Ax"]}
    fv = 1
    if traj == "circle":
        t = np.arange(0.0, cc.T_FINE, g.PASSO_SIM)
        ref = apl.circonferenza(t)
        vero = lambda x: np.stack([apl.CENTRO + apl.RAGGIO * np.cos(apl.OMEGA * x), apl.RAGGIO * np.sin(apl.OMEGA * x), np.full_like(x, apl.QUOTA)], axis=1)
        tag = "circle"
    else:
        ds, st, *kv = traj
        fv = kv[0] if kv else 1
        tr = G["TR"][ds]
        if fv == 1:
            durata, ts, ps, vero = tr["durata"], tr["stime"][st]["t"], tr["stime"][st]["p"], tr["vero"]
        else:
            durata, ts, ps, vero = scala(tr, st, fv)
        t = np.arange(0.0, durata, g.PASSO_SIM)
        psi, th = angoli(ps)
        j = np.clip(np.searchsorted(ts, t, side="right") - 1, 0, len(psi) - 1)
        ref = {"azimut": psi[j], "polare": th[j]}
        tag = f"{ds}_{st}" + ("" if fv == 1 else f"_x{fv:g}")
    R = cc.simula_coppia(G["Ax"], C, ref, t, sc, G["righe"], cc.SEME if seme is None else seme, G["cl"], float(ref["polare"][0]))
    k = np.arange(0, len(t), DECIMA)
    tk = t[k]
    pv = vero(tk)
    psi_l, th_l = ref["azimut"][0] + R["azimut"]["los"][k], ref["polare"][0] + R["polare"]["los"][k]
    E, miss = errore(versore(psi_l, th_l), pv)
    Es, ms = errore(versore(ref["azimut"][k], ref["polare"][k]), pv)
    Et = np.linalg.norm(np.cross(versore(psi_l, th_l), versore(ref["azimut"][k], ref["polare"][k])), axis=1)
    # stesso tratto di percorso escluso a ogni velocita'
    t_tr = T_TRANS / fv
    m = tk >= t_tr
    q = {"traj": tag, "controller": et, "seed": seme, "mean_N_m": media, "speed_factor": fv} | statistiche(tk, E, miss, t_tr)
    q |= {"estimator_alone": statistiche(tk, Es, ms, t_tr), "loop_E_rms_mrad": 1e3 * float(np.sqrt(np.mean(Et[m] ** 2))), "loop_E_max_mrad": 1e3 * float(Et[m].max()), "duration_s": float(tk[-1]), "evaluated_s": float(tk[-1] - t_tr)}
    q |= {f"i_max_{a}_A": float(np.abs(R[b]["i"][k][m]).max()) for a, b in (("az", "azimut"), ("el", "polare"))}
    q |= {f"i_rms_{a}_A": float(np.sqrt(np.mean(R[b]["i"][k][m] ** 2))) for a, b in (("az", "azimut"), ("el", "polare"))}
    q |= {f"los_{a}_rms_mrad": 1e3 * float(np.sqrt(np.mean((ref[b][k] - ref[b][0] - R[b]["los"][k])[m] ** 2))) for a, b in (("az", "azimut"), ("el", "polare"))}
    nome = f"{tag}__{'cascade' if et == 'cascade' else 'lqg'}__{'none' if seme is None else seme}{'_mean' if media else ''}"
    np.savez_compressed(CARTELLA / "series" / f"{nome}.npz", t=tk.astype(np.float32), E=E.astype(np.float32), miss=miss.astype(np.float32), miss_est=ms.astype(np.float32), los_az=R["azimut"]["los"][k].astype(np.float32), los_el=R["polare"]["los"][k].astype(np.float32),
                        ref_az=(ref["azimut"][k] - ref["azimut"][0]).astype(np.float32), ref_el=(ref["polare"][k] - ref["polare"][0]).astype(np.float32), i_az=R["azimut"]["i"][k].astype(np.float32), i_el=R["polare"]["i"][k].astype(np.float32))
    print(nome, "| E rms %.3f mrad | miss rms %.2f cm | on target %.1f %%" % (q["E_rms_mrad"], q["miss_rms_cm"], q["on_target_pct"]), "| est %.1f %%" % q["estimator_alone"]["on_target_pct"], flush=True)
    return q


COL = {"cascade": "#eb6834", "lqg": "#2a78d6"}
NOME = {"cascade": "Cascade", "lqg": "LQG-PI"}
GRIGIO = "#6b6b6b"
SCENARI_FIG = [("circle", "Circle"), ("dataset7_EKF2", "dataset_7\nEKF2"), ("dataset7_config_v3", "dataset_7\nconfig_v3"), ("dataset8_EKF2", "dataset_8\nEKF2"), ("dataset8_config_v3", "dataset_8\nconfig_v3")]


def serie(tag: str, ctrl: str, seme, media: bool = False) -> dict:
    return dict(np.load(CARTELLA / "series" / f"{tag}__{ctrl}__{'none' if seme is None else seme}{'_mean' if media else ''}.npz"))


def stile_tesi() -> None:
    import matplotlib as mpl
    mpl.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 10, "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 8.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e3e2de", "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": True, "legend.framealpha": 0.92, "legend.edgecolor": "#d0cfca", "pdf.fonttype": 42, "lines.linewidth": 1.1, "figure.dpi": 150})


def salva_fig(fig, nome: str) -> None:
    import matplotlib.pyplot as plt
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(CARTELLA / "figures" / f"{nome}.pdf", bbox_inches="tight")
    fig.savefig(CARTELLA / "figures" / f"{nome}.png", bbox_inches="tight", dpi=170)
    plt.close(fig)


def fascio(ax, etichetta: bool = True) -> None:
    ax.axhline(100 * R_FASCIO, color="#0b0b0b", lw=0.9, ls=(0, (4, 2)), label="Beam radius, 8 cm" if etichetta else None, zorder=4)


def fig_geometria(TR: dict) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 4.6), sharex=True)
    for q, (ds, lab) in zip(ax, DATASET.items()):
        tr = TR[ds]
        t = np.arange(0.0, tr["durata"], 0.02)
        v = tr["vero"](t)
        p0 = v[0]
        q.plot(v[:, 1] - p0[1], v[:, 0] - p0[0], color="#0b0b0b", lw=1.5, label="Ground truth", zorder=3)
        for n, c, ls in (("EKF2", "#9a6fb0", "--"), ("config_v3", "#1baf7a", "-")):
            p = tr["stime"][n]["p"]
            q.plot(p[:, 1] - p0[1], p[:, 0] - p0[0], color=c, ls=ls, lw=1.0, label=n, zorder=2)
        q.plot(0, 0, "o", color="#0b0b0b", ms=5, zorder=5, label="Start of the evaluation")
        q.set_title(lab)
        q.set_ylabel("Along the line of sight, $x$ [m]")
        q.set_aspect("equal", adjustable="box")
        q.set_xlim(-4.8, 2.6)
        q.annotate("", xy=(2.05, -1.25), xytext=(2.05, -0.2), arrowprops={"arrowstyle": "->", "color": GRIGIO, "lw": 1.0})
        q.text(2.05, -0.1, "to LUCY,\n100 m", ha="center", va="bottom", fontsize=8.5, color=GRIGIO)
    ax[1].set_xlabel("Lateral position, $y$ [m]")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    salva_fig(fig, "trajectories")


def fig_disturbo() -> None:
    import matplotlib.pyplot as plt
    t = np.arange(0.0, 60.0, g.PASSO_SIM)
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.8), gridspec_kw={"width_ratios": [1.5, 1]})
    d = {n: cc.disturbo({"gauss": SIGMA_TW}, t, SEMI[0] + i) for i, n in enumerate(("Azimuth", "Elevation"))}
    k = (t >= 20.0) & (t <= 22.0)
    for (n, x), c in zip(d.items(), ("#2a78d6", "#eb6834")):
        ax[0].plot(t[k] - 20.0, x[k], color=c, lw=0.9, label=n)
    ax[0].set_xlabel("Time [s]")
    ax[0].set_ylabel("Disturbance torque $\\tau_w$ [N m]")
    ax[0].legend(loc="upper right", ncol=2)
    z = np.linspace(-0.4, 0.4, 400)
    ax[1].hist(d["Azimuth"], bins=60, density=True, color="#2a78d6", alpha=0.45, label="Azimuth samples")
    ax[1].hist(d["Elevation"], bins=60, density=True, color="#eb6834", alpha=0.35, label="Elevation samples")
    ax[1].plot(z, np.exp(-0.5 * (z / SIGMA_TW) ** 2) / (SIGMA_TW * np.sqrt(2 * np.pi)), color="#0b0b0b", lw=1.1, label="$\\mathcal{N}(0,\\,0.1^2)$")
    ax[1].set_xlabel("$\\tau_w$ [N m]")
    ax[1].set_ylabel("Probability density [1/(N m)]")
    ax[1].set_ylim(0, 5.6)
    ax[1].legend(loc="upper right", fontsize=7.5)
    fig.tight_layout()
    salva_fig(fig, "disturbance")


def fig_circonferenza() -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0), sharey=True)
    for q, (se, tit) in zip(ax, ((None, "No disturbance"), (SEMI[0], "Disturbance, $\\sigma_{\\tau}$ = 0.1 N m"))):
        for c in ("cascade", "lqg"):
            x = serie("circle", c, se)
            q.plot(x["t"][::20], 100 * x["miss"][::20], color=COL[c], label=NOME[c])
        fascio(q)
        q.axvspan(0, T_TRANS, color="#e9e8e4", zorder=0, label="Excluded transient")
        q.set_title(tit)
        q.set_xlabel("Time [s]")
        q.set_xlim(0, 36.4)
    ax[0].set_ylabel("Miss distance at the drone [cm]")
    ax[0].set_ylim(0, None)
    ax[0].legend(loc="upper right", ncol=2)
    fig.tight_layout()
    salva_fig(fig, "circle_miss")


def fig_dataset(ds: str) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.0), sharex=True, sharey="row")
    for i, st in enumerate(STIMATORI):
        tag = f"{ds}_{st}"
        for j, se in enumerate((None, SEMI[0])):
            q = ax[i, j]
            x0 = serie(tag, "cascade", se)
            q.plot(x0["t"][::20], 100 * x0["miss_est"][::20], color=GRIGIO, lw=0.8, label="Estimator alone (ideal pointer)", zorder=1)
            for c in ("cascade", "lqg"):
                x = serie(tag, c, se)
                q.plot(x["t"][::20], 100 * x["miss"][::20], color=COL[c], lw=0.8, label=NOME[c], zorder=2)
            fascio(q)
            q.axvspan(0, T_TRANS, color="#e9e8e4", zorder=0, label="Excluded transient")
            q.set_xlim(0, float(x["t"][-1]))
            q.set_title(f"{st}, " + ("no disturbance" if se is None else "disturbance, $\\sigma_{\\tau}$ = 0.1 N m"))
        ax[i, 0].set_ylabel("Miss distance [cm]")
    for q in ax[1]:
        q.set_xlabel("Time from the start of the evaluation [s]")
    ax[0, 0].set_ylim(0, None)
    ax[1, 0].set_ylim(0, None)
    h, l = ax[1, 1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.02), fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    salva_fig(fig, f"{DATASET[ds]}_miss")


def fig_efficace(C: list) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 5.2), sharex=True)
    w = 0.34
    x = np.arange(len(SCENARI_FIG))
    for q, dist, tit in ((ax[0], False, "No disturbance"), (ax[1], True, "Disturbance, $\\sigma_{\\tau}$ = 0.1 N m: mean of three realisations, whiskers show the range")):
        ide = [next(y for y in C if y["traj"] == tag)["estimator_alone"]["on_target_pct"] for tag, _ in SCENARI_FIG]
        for off, et, c in ((-w / 2 - 0.01, "cascade", "cascade"), (w / 2 + 0.01, "LQG-PI tuned", "lqg")):
            v, lo, hi = [], [], []
            for tag, _ in SCENARI_FIG:
                p = [y["on_target_pct"] for y in C if y["traj"] == tag and y["controller"] == et and (y["seed"] is not None) == dist and not y.get("mean_N_m")]
                v.append(np.mean(p)); lo.append(np.mean(p) - min(p)); hi.append(max(p) - np.mean(p))
            q.bar(x + off, v, w, color=COL[c], label=NOME[c], zorder=2, yerr=[lo, hi] if dist else None, error_kw={"elinewidth": 0.9, "ecolor": "#0b0b0b", "capsize": 2})
            for xi, vi, id_ in zip(x + off, v, ide):
                q.text(xi, max(vi, id_ if id_ < vi + 7 else vi) + 2.5, f"{vi:.1f}", ha="center", fontsize=8, color="#0b0b0b")
        q.hlines(ide, x - 0.37, x + 0.37, color="#0b0b0b", lw=1.5, label="Estimator alone (ideal pointer)", zorder=3)
        q.set_ylim(0, 112)
        q.set_yticks(range(0, 101, 20))
        q.set_ylabel("Effective pointing [%]")
        q.set_title(tit)
        q.grid(axis="x", visible=False)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels([n for _, n in SCENARI_FIG])
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    salva_fig(fig, "effective_pointing")


def fig_cdf() -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 3.3), sharey=True)
    def cdf(q, x, **kw):
        m = x["t"] >= T_TRANS
        v = np.sort(100 * x["miss"][m])
        q.plot(v, np.linspace(0, 100, len(v)), **kw)
    for c in ("cascade", "lqg"):
        cdf(ax[0], serie("circle", c, None), color=COL[c], label=NOME[c])
    ax[0].set_title("Circle")
    for q, ds in zip(ax[1:], DATASET):
        for c in ("cascade", "lqg"):
            cdf(q, serie(f"{ds}_config_v3", c, None), color=COL[c], label=f"{NOME[c]}, config_v3")
            cdf(q, serie(f"{ds}_EKF2", c, None), color=COL[c], ls="--", label=f"{NOME[c]}, EKF2")
        q.set_title(DATASET[ds])
    for q in ax:
        q.axvline(100 * R_FASCIO, color="#0b0b0b", lw=0.9, ls=(0, (4, 2)), label="Beam radius, 8 cm")
        q.set_xscale("log")
        q.set_xlim(0.3, 300)
        q.set_ylim(0, 100)
        q.set_xlabel("Miss distance [cm]")
    ax[0].set_ylabel("Fraction of time below [%]")
    h, l = ax[2].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.1), fontsize=8.5)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    salva_fig(fig, "miss_cdf")


def fig_velocita(S: dict, V: dict) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.4), sharex=True, sharey="row")
    for j, (ds, lab) in enumerate(DATASET.items()):
        v0 = V["mean_speed_m_s"][ds]
        righe = [c for c in S["cases"] if c["traj"] == f"{ds}_config_v3" and c["seed"] is None] + [c for c in V["cases"] if c["traj"].startswith(f"{ds}_config_v3")]
        for et, cod in (("cascade", "cascade"), ("LQG-PI tuned", "lqg")):
            r = sorted((c for c in righe if c["controller"] == et), key=lambda c: c.get("speed_factor", 1))
            x = np.array([v0 * c.get("speed_factor", 1) for c in r])
            fuori = np.array([max(c["i_max_az_A"], c["i_max_el_A"]) > I_LIMITE for c in r])
            for q, y in ((ax[0, j], np.array([c["on_target_pct"] for c in r])), (ax[1, j], np.array([c["loop_E_rms_mrad"] for c in r]))):
                q.plot(x, y, color=COL[cod], lw=1.2, label=NOME[cod], zorder=3)
                q.plot(x[~fuori], y[~fuori], "o", color=COL[cod], ms=4.5, zorder=4)
                if fuori.any():
                    q.plot(x[fuori], y[fuori], "o", mfc="white", mec=COL[cod], mew=1.2, ms=5, zorder=5, label=f"{NOME[cod]}, current above 18 A")
            if et == "cascade":
                ax[0, j].plot(x, [c["estimator_alone"]["on_target_pct"] for c in r], color=GRIGIO, ls="--", lw=1.1, marker="s", ms=3.5, label="Estimator alone (ideal pointer)", zorder=2)
        ax[1, j].axhline(1e3 * R_FASCIO / D_LUCY, color="#0b0b0b", lw=0.9, ls=(0, (4, 2)), label="Beam half-angle at 100 m, 0.8 mrad", zorder=1)
        ax[0, j].set_title(lab)
        ax[1, j].set_xlabel("Mean speed of the drone [m/s]")
    for q in ax.flat:
        q.set_xscale("log")
        q.set_xticks([0.2, 0.5, 1.0, 2.0])
        q.set_xticklabels(["0.2", "0.5", "1", "2"])
        q.minorticks_off()
    ax[0, 0].set_ylim(0, 105)
    ax[0, 0].set_ylabel("Effective pointing [%]")
    ax[1, 0].set_ylim(0, None)
    ax[1, 0].set_ylabel("Loop error, rms [mrad]")
    h, l = [], []
    for q in (ax[0, 1], ax[1, 1]):
        for hh, ll in zip(*q.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh); l.append(ll)
    fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.04), fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    salva_fig(fig, "speed")


def sensibilita(ds: str, st: str) -> dict:
    """Stimatore da solo (puntatore ideale) con altre scelte: allineamento dai soli primi 5 s (posizione e imbardata dall'assetto), e LUCY lungo y invece che lungo x."""
    from scipy.spatial.transform import Rotation as Rt
    tv, pv, qv = carica(DATI / "truths_leveled" / f"{ds}.txt")
    tr = traiettoria(ds)
    te, pe, qe = carica(DATI / "algorithms" / STIMATORI[st] / ds / "run_1.txt")
    k = (te >= tr["t0"]) & (te <= tr["t0"] + tr["durata"])
    te, pe, qe = te[k], pe[k], qe[k]
    pg = np.stack([np.interp(te, tv, pv[:, i]) for i in range(3)], axis=1)
    k0 = te <= te[0] + 5.0
    dR = Rt.from_quat(qv[np.clip(np.searchsorted(tv, te[k0]), 0, len(tv) - 1)]) * Rt.from_quat(qe[k0]).inv()
    z = dR.as_euler("zyx")[:, 0]
    yaw = float(np.arctan2(np.sin(z).mean(), np.cos(z).mean()))
    Rz0 = np.array([[np.cos(yaw), -np.sin(yaw), 0.0], [np.sin(yaw), np.cos(yaw), 0.0], [0.0, 0.0, 1.0]])
    tz0 = (pg[k0] - pe[k0] @ Rz0.T).mean(0)
    Rz, tz, _ = allinea_posyaw(pe, pg)
    pa, pa0 = pe @ Rz.T + tz, pe @ Rz0.T + tz0
    m = te >= te[0] + T_TRANS
    o = tr["origine"]
    eta = lambda p: 100 * float(np.mean(errore(versore(*angoli(p - o)), pg - o)[1][m] <= R_FASCIO))
    # LUCY lungo y: stessa distanza e quota
    oy = pv[0] - np.array([0.0, np.sqrt(D_LUCY ** 2 - QUOTA ** 2), QUOTA])
    etay = 100 * float(np.mean(errore(versore(*angoli(pa - oy)), pg - oy)[1][m] <= R_FASCIO))
    e0 = np.linalg.norm(pa0 - pg, axis=1)
    return {"eta_full_alignment_pct": eta(pa), "eta_start_alignment_pct": eta(pa0), "ate_start_alignment_m": float(np.sqrt(np.mean(e0 ** 2))), "eta_lucy_along_y_pct": etay}


def derivati(C: list) -> dict:
    """Effetto del disturbo: differenza della linea di vista dallo stesso caso senza disturbo, rms su due assi, per realizzazione; e rms della coppia usata."""
    out = {}
    for tag, _ in SCENARI_FIG:
        for c in ("cascade", "lqg"):
            b = serie(tag, c, None)
            m = b["t"] >= T_TRANS
            out[f"{tag}__{c}"] = [1e6 * float(np.sqrt(np.mean((x["los_az"] - b["los_az"])[m] ** 2 + (x["los_el"] - b["los_el"])[m] ** 2))) for x in (serie(tag, c, s) for s in SEMI)]
    t = np.arange(0.0, 60.0, g.PASSO_SIM)
    d = {n: cc.disturbo({"gauss": SIGMA_TW}, t, SEMI[0] + i) for i, n in enumerate(("azimuth", "elevation"))}
    stat = {n: {"mean_N_m": float(x.mean()), "rms_N_m": float(x.std()), "peak_N_m": float(np.abs(x).max()), "excess_kurtosis": float(((x - x.mean()) ** 4).mean() / x.var() ** 2 - 3)} for n, x in d.items()}
    sens = {f"{ds}_{st}": sensibilita(ds, st) for ds in DATASET for st in STIMATORI}
    # Effetto della parte statica: stessa realizzazione con e senza media
    mm = {}
    if any(y.get("mean_N_m") for y in C):
        for tag, _ in SCENARI_FIG:
            for c in ("cascade", "lqg"):
                v = []
                for s in SEMI:
                    x0, x1 = serie(tag, c, s), serie(tag, c, s, True)
                    m = x0["t"] >= T_TRANS
                    da, de = (x1["los_az"] - x0["los_az"])[m], (x1["los_el"] - x0["los_el"])[m]
                    v.append({"los_az_mean_urad": 1e6 * float(da.mean()), "los_el_mean_urad": 1e6 * float(de.mean()), "total_rms_urad": 1e6 * float(np.sqrt(np.mean(da ** 2 + de ** 2))),
                              "total_peak_urad": 1e6 * float(np.sqrt(da ** 2 + de ** 2).max()), "i_az_mean_A": float((x1["i_az"] - x0["i_az"])[m].mean()), "i_el_mean_A": float((x1["i_el"] - x0["i_el"])[m].mean()),
                              "miss_rms_change_cm": 100 * float(np.sqrt(np.mean(x1["miss"][m] ** 2)) - np.sqrt(np.mean(x0["miss"][m] ** 2)))})
                mm[f"{tag}__{c}"] = v
    return {"effect_urad": out, "torque": stat, "sensitivity": sens, "mean_effect": mm}


def figure() -> None:
    stile_tesi()
    S = json.loads((CARTELLA / "summary.json").read_text())
    TR = {ds: traiettoria(ds) for ds in DATASET}
    fig_geometria(TR)
    fig_disturbo()
    fig_circonferenza()
    for ds in DATASET:
        fig_dataset(ds)
    fig_efficace(S["cases"])
    fig_cdf()
    if (CARTELLA / "speed.json").exists():
        fig_velocita(S, json.loads((CARTELLA / "speed.json").read_text()))
    (CARTELLA / "derived.json").write_text(json.dumps(derivati(S["cases"]), indent=1))


def main() -> None:
    if SOLO_FIGURE:
        figure()
        return
    (CARTELLA / "series").mkdir(parents=True, exist_ok=True)
    if PROVA:
        import time
        t0 = time.time()
        q = un_caso(("circle", "LQG-PI tuned", None))
        print(q, time.time() - t0)
        return
    traj = ["circle"] + [(ds, st) for ds in DATASET for st in STIMATORI]
    if VELOCITA:
        compiti = [((ds, "config_v3", fv), et, None) for ds in DATASET for fv in FATTORI_V for et in CONTROLLORI]
        with ProcessPoolExecutor(int(os.environ.get("LAVORI", "5"))) as ex:
            righe = list(ex.map(un_caso, compiti))
        v = {}
        for ds in DATASET:
            tr = G["TR"][ds]
            p = tr["vero"](np.arange(0.0, tr["durata"], 0.02))
            v[ds] = float(np.mean(np.linalg.norm(np.diff(p, axis=0), axis=1)) / 0.02)
        (CARTELLA / "speed.json").write_text(json.dumps({"factors": FATTORI_V, "mean_speed_m_s": v, "cases": righe}, indent=1))
        return
    if MEDIA:
        compiti = [(x, et, s, MEDIA_TW) for x in traj for et in CONTROLLORI for s in SEMI]
        with ProcessPoolExecutor(int(os.environ.get("LAVORI", "5"))) as ex:
            righe = list(ex.map(un_caso, compiti))
        S = json.loads((CARTELLA / "summary.json").read_text())
        S["cases"] = [y for y in S["cases"] if not y.get("mean_N_m")] + righe
        S["mean_N_m"] = MEDIA_TW
        (CARTELLA / "summary.json").write_text(json.dumps(S, indent=1))
        return
    compiti = [(x, et, s) for x in traj for et in CONTROLLORI for s in (None, *SEMI)]
    with ProcessPoolExecutor(int(os.environ.get("LAVORI", "5"))) as ex:
        righe = list(ex.map(un_caso, compiti))
    geo = {ds: {"origin_m": G["TR"][ds]["origine"].tolist(), "window_s": G["TR"][ds]["durata"], "estimators": {n: {k: s[k] for k in ("ate_m", "yaw_deg", "rate_Hz")} for n, s in G["TR"][ds]["stime"].items()}} for ds in DATASET}
    (CARTELLA / "summary.json").write_text(json.dumps({"cases": righe, "geometry": geo, "sigma_N_m": SIGMA_TW, "beam_radius_m": R_FASCIO, "seeds": SEMI, "distance_m": D_LUCY, "height_m": QUOTA}, indent=1))


if __name__ == "__main__":
    main()
