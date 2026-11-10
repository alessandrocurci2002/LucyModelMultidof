"""Passi 2-4 della taratura dell'LQG-PI (Gawronski 2008, par. 8.2), iterazione e prova su V, dopo il passo 1 di ``analisi_lqg_passo1_flessibili.py``.

Funzioni di trasferimento e gradino sul modello ridotto con Maxwell, Q = V. Obiettivo dell'utente: assestamento basso, banda alta.
Passo 2: si alza q_p finche' la banda raggiunge la fondamentale. Passo 3: si alza q_i finche' compaiono oscillazioni, entro (7.9).
Passo 4: si alzano i pesi dei modi i cui picchi sono troppo alti (D14). Poi taratura fine sull'indice (8.27) ridotto a t_s e banda, si
rialza q_p, e V = rho Q: su Maxwell T_ry non dipende da V, quindi V si giudica anche sul plant Dahl. Le soglie non date dal libro sono in
``limiti`` e nel registro delle decisioni del report.

    .venv/bin/python -I analisi_lqg_passi_2_4.py [--senza-cuscinetti]

Risultati in results_lqg_tuning/steps_2_4/; con ``--senza-cuscinetti``, modello interno senza cuscinetti, in results_lqg_tuning/no_bearings/steps_2_4/.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]

import matplotlib.pyplot as plt  # noqa: E402
from scipy.linalg import expm  # noqa: E402
from scipy.signal import butter, sosfiltfilt  # noqa: E402

import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_passo1_flessibili as p1  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning" / ("no_bearings" if p1.SENZA else "") / "steps_2_4"
SCALA_QP = (0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8)
SCALA_QF = (1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7)
RHO = (1e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1e3)
I_LIMITE = 18.0
MOD_MIN = 0.5
PICCO_MAX = -6.0
"""Picco lasciato dal libro in Fig. 8.5 (D14)."""
MIGLIORAMENTO = 0.05
"""Variazione dell'indice sotto la quale un'altra iterazione non vale."""
T_STEP = 20.0


def comandi_tipici(asse: str) -> list:
    """Un doppietto per ampiezza; per l'8 gradi in azimut si salta la finestra anomala."""
    fin = g.finestre_asse(asse)
    out = []
    for t, f in fin.items():
        ok = [w for w in f if (asse, t, w.nome[-2:]) != g.ANOMALA]
        out.append(ok[0])
    return out


def doppietti(b, c, W) -> dict:
    ms, picco = [], 0.0
    for w in W:
        r = g.simula(b.a, c, w)
        ms.append(g.metriche(w, r["p"], r["i"]))
        picco = max(picco, r["max_i_cmd"])
    return {"i_max_A": float(picco), "t_set_median_s": float(np.median([m["t_set_s"] for m in ms])), "t_set_max_s": float(max(m["t_set_s"] for m in ms)),
            "err_rms_median_mdeg": float(np.median([m["err_rms_mdeg"] for m in ms])), "over_max_mdeg": float(max(m["over_mdeg"] for m in ms))}


def oscillazioni(r: dict) -> int:
    """Estremi dell'errore fuori dalla banda del 2 % dopo il primo passaggio sul comando, sotto 20 Hz: le oscillazioni del guadagno
    integrale hanno frequenza minore della fondamentale (par. 8.1.7), la vibrazione dei modi e' del passo 4."""
    if not np.all(np.isfinite(r["y"])) or np.abs(r["y"]).max() > 1e3:
        return 999
    y = sosfiltfilt(butter(4, 20.0, "lowpass", fs=1 / g.PASSO_MASTER, output="sos"), r["y"])
    t, e = r["t"], (y - p1.GRADINO) / p1.GRADINO
    e = e[t >= p1.T_AVVIO]
    i0 = np.flatnonzero(e >= 0)
    if not i0.size:
        return 0
    e = e[i0[0]:]
    d = np.diff(e)
    est = np.flatnonzero(d[:-1] * d[1:] < 0) + 1
    return int(np.sum(np.abs(e[est]) > 0.02))


class Asse2:
    """Un asse: ``p1.Banco``, comandi tipici e plant completi."""

    def __init__(self, b, mx):
        self.b, self.W = b, comandi_tipici(b.asse)
        a, an = b.a, b.a.anello_v
        self.f_fond = b.modi[0]["f_Hz"]
        self.piante = {n: (A, an.B[:, 0], an.uscite["p"]) for n, A in g.varianti_cuscinetto(a).items()}
        if mx:
            self.piante["full_maxwell"] = lqg.con_cuscinetti(an.A, an.B[:, 0], an.B[:, 1], an.uscite["p"], np.r_[a.m.v, np.zeros(an.N - an.n)], np.r_[np.zeros(an.n), a.m.v, np.zeros(an.N - 2 * an.n)], 0.0, mx["k_1"], mx["a"])

    def controllore(self, q_p, q_i, q_f=None, rho=1.0):
        b = self.b
        q = np.full(b.n - 1, p1.Q_F)
        for i, x in (q_f or {}).items():
            q[[j - 1 for j in b.modi[int(i)]["stati"]]] = x
        return lqg.sintesi(b.a.Ap, b.a.Bp, b.a.Cp, g.PASSO_MASTER, q_i, q_p, q, rho_e=rho, v_continua=True)

    def punto(self, q_p, q_i, q_f=None, rho=1.0, con_doppietti=True) -> dict:
        c = self.controllore(q_p, q_i, q_f, rho)
        v = self.b.valuta(c)
        r = self.b.gradino(c, T_STEP)
        out = {"q_p": float(q_p), "q_i": float(q_i), "q_f": {str(k): x for k, x in (q_f or {}).items()}, "rho": rho} | v | {"step": r["metriche"], "extrema_outside_2pct": oscillazioni(r), "k_i_over_critical": float(v["k_i"] / (0.25 * v["k_p"] ** 2)),
                                                                                                                              "full_plants": {n: p1.peggiore(p1.poli_chiusi(c, *x)) for n, x in self.piante.items()}}
        if con_doppietti:
            out["doublets"] = doppietti(self.b, c, self.W)
        return out


def limiti(x: dict, f_fond: float, integrale: bool = False) -> list[str]:
    """Le condizioni violate: banda alla fondamentale e corrente (libro), min |1 + L| e stabilita' sui plant completi (scelte mie), piu' (7.9) e oscillazioni per il passo 3."""
    v = []
    if x["bandwidth_Hz"] >= f_fond:
        v.append("bandwidth at the fundamental")
    if "doublets" in x and x["doublets"]["i_max_A"] > I_LIMITE:
        v.append("current")
    if x["modulus"] < MOD_MIN:
        v.append("min |1 + L|")
    if any(q["abs_z"] > 1 + 1e-9 for q in x["full_plants"].values()):
        v.append("unstable on a full plant")
    if integrale and x["k_i_over_critical"] > 1:
        v.append("critical gain")
    if integrale and x["extrema_outside_2pct"] >= 2:
        v.append("oscillations")
    return v


def alza(asse: Asse2, scala, fun, integrale: bool, passi: int = 5) -> tuple[list, float, list[str]]:
    """Sale lungo ``scala``, poi biseca fra l'ultimo punto buono e il primo che viola un limite. ``fun(q)`` da' il punto."""
    punti, buono, cattivo, motivo = [], None, None, []
    for q in scala:
        x = fun(q)
        x["limits"] = limiti(x, asse.f_fond, integrale)
        punti.append(x)
        print(asse.b.chiave, "%.4g" % q, round(x["k_p"], 2), round(x["k_i"], 2), round(x["bandwidth_Hz"], 3), round(x["step"]["t_set_2pct_s"], 3), round(x["doublets"]["i_max_A"], 1), round(x["modulus"], 3), x["limits"], flush=True)
        if x["limits"]:
            cattivo, motivo = q, x["limits"]
            break
        buono = q
    if buono is None:
        return punti, float(scala[0]), motivo
    if cattivo is None:
        return punti, float(buono), []
    lo, hi = buono, cattivo
    for _ in range(passi):
        m = float(np.sqrt(lo * hi))
        x = fun(m)
        x["limits"] = limiti(x, asse.f_fond, integrale)
        punti.append(x)
        print(asse.b.chiave, " bis %.4g" % m, round(x["k_p"], 2), round(x["doublets"]["i_max_A"], 1), round(x["modulus"], 3), x["limits"], flush=True)
        if x["limits"]:
            hi, motivo = m, x["limits"]
        else:
            lo = m
    return punti, float(lo), motivo


def passo_2(asse):
    return alza(asse, SCALA_QP, lambda q: asse.punto(q, p1.Q_I), False)


def passo_3(asse, q_p, q_f=None):
    kp = asse.punto(q_p, p1.Q_I, q_f, con_doppietti=False)["k_p"]
    q_c = (0.25 * kp ** 2) ** 2
    n = int(np.ceil(np.log10(max(q_c, 1.0) / p1.Q_I) * 2)) + 2
    scala = [p1.Q_I * 10 ** (0.5 * k) for k in range(n)]
    return alza(asse, scala, lambda q: asse.punto(q_p, q, q_f), True)


def passo_4(asse, q_p, q_i, q_f=None) -> tuple[list, dict]:
    """Alza il peso dei soli modi con il picco oltre PICCO_MAX, uno alla volta, finche' il picco scende o un limite lo impedisce."""
    q_f = dict(q_f or {})
    storia = []
    x = asse.punto(q_p, q_i, q_f)
    x["limits"] = limiti(x, asse.f_fond)
    storia.append(x)
    for im in range(len(asse.b.modi)):
        if x["modes"][im]["band_max_dB"] < PICCO_MAX:
            continue
        for q in SCALA_QF:
            prova = q_f | {im: q}
            y = asse.punto(q_p, q_i, prova)
            y["limits"] = limiti(y, asse.f_fond)
            y["raised_mode"] = im
            storia.append(y)
            print(asse.b.chiave, "passo 4 modo", im, "%.3g" % q, round(y["modes"][im]["band_max_dB"], 1), round(y["modulus"], 3), y["limits"], flush=True)
            if y["limits"]:
                break
            q_f, x = prova, y
            if y["modes"][im]["band_max_dB"] < PICCO_MAX:
                break
    return storia, q_f


def indice(x: dict, rif: dict) -> float:
    """(8.27) ridotta a t_s e banda, ciascuno rispetto al punto di riferimento."""
    return float(x["step"]["t_set_2pct_s"] / rif["step"]["t_set_2pct_s"] - x["bandwidth_Hz"] / rif["bandwidth_Hz"])


def fine(asse, rif, q_f) -> tuple[list, dict]:
    griglia = []
    q_p0, q_i0 = rif["q_p"], rif["q_i"]
    kp = rif["k_p"]
    q_c = (0.25 * kp ** 2) ** 2
    assi_i = sorted(set([float(x) for x in np.geomspace(q_i0, max(q_c, 10 * q_i0), 7)]))
    for fq in (0.15, 0.3, 0.5, 0.7, 1.0, 1.4):
        for qi in assi_i:
            x = asse.punto(q_p0 * fq, qi, q_f, con_doppietti=False)
            lim = limiti(x, asse.f_fond, True)
            if not lim:
                x["doublets"] = doppietti(asse.b, asse.controllore(x["q_p"], qi, q_f), asse.W)
                lim = limiti(x, asse.f_fond, True)
            x["limits"] = lim
            x["index"] = indice(x, rif)
            griglia.append(x)
    ok = [x for x in griglia if not x["limits"]]
    return griglia, min(ok, key=lambda x: x["index"]) if ok else rif


def rialza(asse, q_p, q_i, q_f, rho=1.0) -> tuple[list, float, list[str]]:
    """Da q_p sale finche' un limite lo ferma; conta anche il picco dei modi, perche' qui il passo 4 e' gia' stato fatto."""
    punti, qs, motivo = [], None, []
    for fat in (1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0):
        x = asse.punto(q_p * fat, q_i, q_f, rho)
        x["limits"] = limiti(x, asse.f_fond) + (["resonance peak"] if any(m["band_max_dB"] >= PICCO_MAX for m in x["modes"]) else [])
        punti.append(x)
        if x["limits"]:
            motivo = x["limits"]
            break
        qs = q_p * fat
    return punti, qs, motivo


def prova_v(asse, q_p, q_i, q_f) -> list[dict]:
    """V = rho Q, a q_p fisso e rialzandolo per ogni rho."""
    out = []
    for rho in RHO:
        punti, qs, motivo = rialza(asse, q_p, q_i, q_f, rho)
        out.append({"rho": rho, "at_same_q_p": punti[0], "q_p_max": qs, "stop": motivo, "at_q_p_max": next((x for x in punti if x["q_p"] == qs), None)})
        print(asse.b.chiave, "rho", rho, "base", round(punti[0]["modulus"], 3), round(punti[0]["doublets"]["i_max_A"], 1), round(punti[0]["doublets"]["t_set_median_s"], 2), round(punti[0]["doublets"]["err_rms_median_mdeg"], 1), punti[0]["limits"], "q_p max", qs, motivo, flush=True)
    return out


def scegli_v(V: list[dict], rif: dict) -> tuple[float, float, list[dict]]:
    """Si cambia V solo se l'indice sul modello migliora e sul plant Dahl ne' t_set ne' l'errore rms peggiorano."""
    esiti = []
    for q in V:
        x = q["at_q_p_max"]
        if q["rho"] == 1.0 or x is None:
            continue
        d, d0 = x["doublets"], rif["doublets"]
        esiti.append({"rho": q["rho"], "index": indice(x, rif), "t_set_ratio": d["t_set_median_s"] / d0["t_set_median_s"], "err_ratio": d["err_rms_median_mdeg"] / d0["err_rms_median_mdeg"],
                      "accepted": bool(indice(x, rif) < -MIGLIORAMENTO and d["t_set_median_s"] <= d0["t_set_median_s"] and d["err_rms_median_mdeg"] <= d0["err_rms_median_mdeg"])})
    ok = [e for e in esiti if e["accepted"]]
    if not ok:
        return 1.0, rif["q_p"], esiti
    b = min(ok, key=lambda e: e["index"])
    return b["rho"], next(q["q_p_max"] for q in V if q["rho"] == b["rho"]), esiti


def misurato(a, x, cfg) -> dict:
    """Fuori dalla procedura: anello con il plant dalla FRF misurata, 1-150 Hz."""
    mis = g.misura_asse(cfg, a.b.asse)
    k = (mis.f >= 1.0) & (mis.f <= 150.0)
    L = g.anello_lqg(a.controllore(x["q_p"], x["q_i"], x["q_f"], x["rho"]), g.plant_misurato(a.b.a, mis.f[k], mis.G[k]), mis.f[k])[0]
    i = int(np.argmin(np.abs(1 + L)))
    return {"modulus": float(np.abs(1 + L[i])), "modulus_Hz": float(mis.f[k][i]), "L_max_dB": float(pt._db(L).max())}


def cascata_su_maxwell(a) -> dict:
    """La cascata attuale sul modello con Maxwell completo, per confrontare banda e gradino a parita' di plant."""
    an, cas, mx = a.b.a, a.b.a.cascata, a.b.mx
    C = cas.uscite["p"]
    A2, B2, C2 = lqg.con_cuscinetti(cas.A, cas.B[:, 0], cas.B[:, 1], C, np.r_[an.m.v, np.zeros(cas.N - cas.n)], np.r_[np.zeros(cas.n), an.m.v, np.zeros(cas.N - 2 * cas.n)], 0.0, mx["k_1"], mx["a"]) if mx else (cas.A, cas.B[:, 0], C)
    T = g.frf(A2, B2, C2, p1.F)
    N = len(A2)
    E = expm(np.block([[A2, B2[:, None]], [np.zeros((1, N + 1))]]) * g.PASSO_MASTER)
    t = np.arange(0.0, T_STEP, g.PASSO_MASTER)
    x, y = np.zeros(N), np.empty(len(t))
    for k, tk in enumerate(t):
        y[k] = C2 @ x
        x = E[:N, :N] @ x + E[:N, N] * (p1.GRADINO if tk >= p1.T_AVVIO else 0.0)
    return {"bandwidth_Hz": g.banda(p1.F, T), "step": p1.metriche_gradino(t, y, 1 / g.PASSO_MASTER)["metriche"], "doublets": doppietti(a.b, None, a.W), "K_pp": an.drive.p.Kpp}


def fig_scala(B, R, chiave, xlab, seg, nome):
    fig, ax = plt.subplots(4, 2, figsize=(13, 10), sharex="col")
    for j, (a, r) in enumerate(zip(B, R)):
        pts = sorted(r[chiave], key=lambda x: x[xlab])
        x = [q[xlab] for q in pts]
        ax[0, j].plot(x, [q["bandwidth_Hz"] for q in pts], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[0, j].axhline(a.f_fond, color=pt.INCHIOSTRO_SECONDARIO, lw=0.8, ls=":"); ax[0, j].set_yscale("log")
        ax[1, j].plot(x, [q["step"]["t_set_2pct_s"] for q in pts], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[1, j].set_yscale("log")
        ax[2, j].plot(x, [q["doublets"]["i_max_A"] for q in pts], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[2, j].axhline(I_LIMITE, color=pt.INCHIOSTRO, lw=0.7, ls="--"); ax[2, j].set_ylim(0, 45)
        ax[3, j].plot(x, [q["modulus"] for q in pts], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[3, j].axhline(MOD_MIN, color=pt.INCHIOSTRO, lw=0.7, ls="--")
        ax[3, j].set_xscale("log"); ax[3, j].set_xlabel(nome)
        for q in ax[:, j]:
            q.axvline(r[chiave + "_selected"], color=pt.INCHIOSTRO, lw=0.8, ls=":")
        ax[0, j].set_title("Azimuth" if a.b.asse == "azimut" else "Elevation", fontsize=9)
    ax[0, 0].set_ylabel("bandwidth [Hz]; dotted: fundamental"); ax[1, 0].set_ylabel("t_s, 2 %, 1 deg step [s]"); ax[2, 0].set_ylabel("largest current, typical commands [A]"); ax[3, 0].set_ylabel("min |1 + L|")
    pt._titolo(fig, f"{p1.MODELLO.capitalize()}; current on the Dahl plant. Vertical dotted line: selected {nome}", seg)
    return fig


def fig_finale(B, R, seg):
    fig, ax = plt.subplots(3, 2, figsize=(13, 9.5))
    for j, (a, r) in enumerate(zip(B, R)):
        for n, (et, x) in enumerate(r["curve"].items(), 1):
            c = a.controllore(x["q_p"], x["q_i"], x["q_f"], x["rho"])
            T = a.b.tf(c)[0]
            st = a.b.gradino(c, 6.0)
            ax[0, j].plot(p1.F, pt._db(T), lw=1.2, label=et, **pt._stile(n))
            k = p1.F >= 30.0
            ax[1, j].plot(p1.F[k], pt._db(T[k]), lw=1.1, **pt._stile(n))
            ax[2, j].plot(st["t"], np.rad2deg(st["y"]), lw=1.2, **pt._stile(n))
        ax[0, j].axhline(PICCO_MAX, color=pt.INCHIOSTRO_SECONDARIO, lw=0.7, ls="--"); ax[1, j].axhline(PICCO_MAX, color=pt.INCHIOSTRO_SECONDARIO, lw=0.7, ls="--")
        ax[0, j].set_xscale("log"); ax[1, j].set_xscale("log"); ax[0, j].set_xlabel("frequency [Hz]"); ax[1, j].set_xlabel("frequency [Hz]"); ax[2, j].set_xlabel("time [s]")
        ax[0, j].set_title("Azimuth" if a.b.asse == "azimut" else "Elevation", fontsize=9); ax[0, j].legend(fontsize=7)
    ax[0, 0].set_ylabel("T_ry [dB]; dashed: -6 dB"); ax[1, 0].set_ylabel("T_ry, zoom on the modes [dB]"); ax[2, 0].set_ylabel("1 deg step [deg]")
    pt._titolo(fig, f"Tuning on the {p1.MODELLO}, Q = V", seg)
    return fig


def fig_v(B, R, seg):
    fig, ax = plt.subplots(4, 2, figsize=(13, 10), sharex="col")
    for j, (a, r) in enumerate(zip(B, R)):
        x = [q["rho"] for q in r["V"]]
        base = [q["at_same_q_p"] for q in r["V"]]
        ax[0, j].plot(x, [q["modulus"] for q in base], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[0, j].axhline(MOD_MIN, color=pt.INCHIOSTRO, lw=0.7, ls="--")
        ax[1, j].plot(x, [q["doublets"]["i_max_A"] for q in base], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[1, j].axhline(I_LIMITE, color=pt.INCHIOSTRO, lw=0.7, ls="--")
        ax[2, j].plot(x, [q["doublets"]["err_rms_median_mdeg"] for q in base], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO)
        mx = [q["at_q_p_max"]["bandwidth_Hz"] if q["at_q_p_max"] else np.nan for q in r["V"]]
        ax[3, j].plot(x, mx, marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[3, j].set_yscale("log")
        ax[3, j].set_xscale("log"); ax[3, j].set_xlabel("rho, V = rho Q")
        ax[0, j].set_title("Azimuth" if a.b.asse == "azimut" else "Elevation", fontsize=9)
    ax[0, 0].set_ylabel("min |1 + L|, internal model"); ax[1, 0].set_ylabel("largest current, Dahl [A]"); ax[2, 0].set_ylabel("median rms error, Dahl [mdeg]"); ax[3, 0].set_ylabel("bandwidth at the largest q_p [Hz]")
    pt._titolo(fig, "Effect of V: at the weights selected with Q = V (rows 1-3), and raising q_p again for each V (row 4)", seg)
    return fig


def fig_doppietti(B, R, seg):
    fig, ax = plt.subplots(2, 2, figsize=(13, 7), sharex="col")
    for j, (a, r) in enumerate(zip(B, R)):
        w = a.W[3]
        c = a.controllore(r["final"]["q_p"], r["final"]["q_i"], r["final"]["q_f"], r["final"]["rho"])
        for n, (et, cc) in enumerate((("LQG-PI", c), ("present cascade", None)), 1):
            x = g.simula(a.b.a, cc, w)
            e = np.rad2deg((w.cmd - w.cmd[0]) - x["p"]) * 1e3
            ax[0, j].plot(w.t, e, lw=1.1, label=et, **pt._stile(n)); ax[1, j].plot(w.t, x["i"], lw=0.9, **pt._stile(n))
        ax[0, j].set_title(("Azimuth" if a.b.asse == "azimut" else "Elevation") + ", " + w.nome, fontsize=9); ax[0, j].legend(fontsize=7); ax[1, j].set_xlabel("time [s]")
    ax[0, 0].set_ylabel("following error [mdeg]"); ax[1, 0].set_ylabel("current [A]")
    pt._titolo(fig, "Final weights on the Dahl plant, measured command of a doublet, no current limit", seg)
    return fig


def senza_array(o):
    if isinstance(o, dict):
        return {str(k): senza_array(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [senza_array(v) for v in o]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def main() -> None:
    par, Kt, S, T, mot, assi = g.carica()
    seg = par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    O = json.loads((RADICE / "results_lqg_tuning" / "operating_point.json").read_text())
    B, R, out = [], [], {}
    for asse, ch, cat, mod in assi:
        mx = None if p1.SENZA else O[ch]["maxwell"]
        a = Asse2(p1.Banco(asse, ch, cat, mod, Kt, S, T, mot, mx), mx)
        p2, qp, m2 = passo_2(a)
        p3, qi, m3 = passo_3(a, qp)
        st4, qf = passo_4(a, qp, qi)
        fine1 = a.punto(qp, qi, qf)
        griglia, migliore = fine(a, fine1, qf)
        secondo = bool(migliore is not fine1 and migliore["index"] < -MIGLIORAMENTO)
        qp2, qi2 = (migliore["q_p"], migliore["q_i"]) if secondo else (qp, qi)
        st4b, qf2 = passo_4(a, qp2, qi2, qf) if secondo else ([], qf)
        punti_r, qp_r, stop_r = rialza(a, qp2, qi2, qf2)
        qp_r = qp_r if qp_r is not None else qp2
        finale_q = a.punto(qp_r, qi2, qf2)
        V = prova_v(a, qp_r, qi2, qf2)
        rho, qp3, esiti_v = scegli_v(V, finale_q)
        finale = a.punto(qp3, qi2, qf2, rho)
        finale["limits"] = limiti(finale, a.f_fond, True)
        casc = cascata_su_maxwell(a)
        check = {k: p1.verifica_dahl(a.b, a.controllore(finale["q_p"], finale["q_i"], finale["q_f"], finale["rho"])) for k in ("final",)}
        r = {"fundamental_Hz": a.f_fond, "typical_commands": [w.nome for w in a.W], "q_p": p2, "q_p_selected": qp, "q_p_stop": m2, "q_i": p3, "q_i_selected": qi, "q_i_stop": m3, "step_4": st4, "q_f_after_step_4": {str(k): v for k, v in qf.items()},
             "after_steps": fine1, "fine_grid": griglia, "fine_best": migliore if migliore is not fine1 else None, "second_iteration": secondo, "step_4_second": st4b, "q_f_final": {str(k): v for k, v in qf2.items()}, "at_Q_equal_V": finale_q,
             "q_p_raised": punti_r, "q_p_raised_selected": qp_r, "q_p_raised_stop": stop_r, "V": V, "V_decision": esiti_v, "rho_selected": rho, "final": finale, "final_measured_plant": misurato(a, finale, par.sperimentale), "cascade": casc, "dahl_check": {q: check["final"][q] for q in ("largest_pole", "y_100_120s_mdeg", "y_minus_y_hat_end_mdeg", "step_20s", "i_max_A")},
             "curve": {"after Step 1": a.punto(p1.Q_P, p1.Q_I, con_doppietti=False), "after Steps 2-4": fine1, "final": finale}}
        out[ch] = r
        B.append(a); R.append(r)
        (CARTELLA / "summary.json").write_text(json.dumps(senza_array(out), indent=1))
    salva = lambda f, n: (pt.salva(f, CARTELLA / "figures", n), plt.close(f))
    salva(fig_scala(B, R, "q_p", "q_p", seg, "q_p"), "step_2")
    salva(fig_scala(B, R, "q_i", "q_i", seg, "q_i"), "step_3")
    salva(fig_finale(B, R, seg), "final_tf_step")
    salva(fig_v(B, R, seg), "effect_of_V")
    salva(fig_doppietti(B, R, seg), "final_doublet")


if __name__ == "__main__":
    main()
