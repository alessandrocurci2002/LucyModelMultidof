#!/usr/bin/env python3
r"""In quali condizioni l'LQG-PI porta un vantaggio: confronto con la cascata e con il P sull'encoder a guadagno alzato (solo asse, plant con Dahl).

Il P a guadagno alzato separa cio' che viene dall'LQG da cio' che viene dal solo :math:`k_p` piu' alto. Per ogni caso: banda e picco di
:math:`T_{ry}`, :math:`T_{wy}`, stabilita', doppietti, e sull'orbita errore della linea di vista ed effetto di disturbi e rumore, come
differenza dallo stesso controllore senza disturbi. Risultati in ``results_lqg_tuning/no_bearings/circle/benefits.json``.
"""

from __future__ import annotations

import dataclasses
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]
SOLO_MEDIE = "--medie" in sys.argv
SOLO_KF = "--kf" in sys.argv
SOLO_VENTO = "--vento" in sys.argv
sys.argv = [sys.argv[0], "--senza-cuscinetti"]

import analisi_lqg_circonferenza as cc  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_punto_lavoro as apl  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402

cm = cc.cm
# Raffiche di vento (par. 5.2): spettro di Davenport, velocita' media e livelli di coppia ipotizzati, come rms della coppia all'asse
V_MEDIA = 5.0
T_VENTO = 120.0
RMS_VENTO = (0.1, 0.5)
SEMI_VENTO = (1, 2)
F_MEDIE = (3.0, 5.0, 6.6, 8.0)
AMP_MEDIE = (0.2, 0.8)
PK = {"azimut": (12.0, 25.0, 60.0), "polare": (12.0, 60.0, 100.0)}
F = np.geomspace(0.01, 300.0, 3000)


COSTRUITI_KI: dict = {}


def costruisci() -> dict:
    par, Kt, S, T, mot, assi = g.carica()
    Fi = json.loads(cc.TUNING.read_text())
    out = {}
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        fi = Fi[ch]["final"]
        m0 = next(m for m in a.modi if len(m.indici) == 2 and m.zeta < 0.1)
        ca = cc.casi(a, fi, m0.f_hz)
        COSTRUITI_KI[asse] = {"ki": ca["LQG-PI, k_i = 0.5 k_p"]}
        C = {"cascade": None} | {f"P, k = {k:g}": cm.PEncoder(k) for k in PK[asse]} | {"LQG-PI tuned": ca["LQG-PI tuned"], "LQG-PI, q_p = 1e5": ca["LQG-PI, q_p = 1e5"]}
        if asse == "polare":
            q, v = np.full(len(a.Ap) - 1, 0.1), np.full(len(a.Ap) - 1, 0.1)
            i = [x - 1 for x in m0.indici]
            q[i], v[i] = 1e4, 1e6
            C["LQG-PI, V_f = 1e6"] = lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, fi["q_i"], fi["q_p"], q, v_f=v, v_continua=True)
        out[asse] = (a, C)
    return out


TUTTI = costruisci()
CI = {asse: x["ki"] for asse, x in COSTRUITI_KI.items()}


def lineare(a, c) -> dict:
    an, var = a.anello_v, g.varianti_cuscinetto(a)
    C, Bu, Bw = an.uscite["p"], an.B[:, 0], an.B[:, 1]
    ritardo = np.exp(-2j * np.pi * F * g.PASSO_MASTER * g.RITARDO) * g.zoh(F)
    out = {}
    for n, A in var.items():
        P, Pw = g.frf(A, Bu, C, F), g.frf(A, Bw, C, F)
        if c is None:
            k = a.drive.p.Kpp
            L, N = k * ritardo * P, k * ritardo * P
        elif isinstance(c, cm.PEncoder):
            L, N = c.k * ritardo * P, c.k * ritardo * P
        else:
            L, N = g.anello_lqg(c, P, F)
        Tw, Tr = pt._db(Pw / (1 + L)), N / (1 + L)
        out[n] = {"band_Hz": g.banda(F, Tr), "T_ry_peak_dB": float(pt._db(Tr).max()), "T_wy_peak_dB": float(Tw.max()), "T_wy_peak_Hz": float(F[Tw.argmax()]), "T_wy_0.1_dB": float(np.interp(0.1, F, Tw)), "T_wy_1_dB": float(np.interp(1.0, F, Tw)), "T_wy_3_dB": float(np.interp(3.0, F, Tw)), "T_wy_10_dB": float(np.interp(10.0, F, Tw)), "T_wy_56_dB": float(np.interp(56.12, F, Tw))}
    return out


def un_caso(arg) -> tuple:
    asse, et = arg
    a, C = TUTTI[asse]
    c = C[et]
    t = np.arange(0.0, cc.T_FINE, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    r = ref[asse] - ref[asse][0]
    riga = np.zeros(a.cascata.N)
    riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
    m = t >= cc.T_AVVIO
    rms = lambda x: float(np.sqrt(np.mean(x[m] ** 2)))
    base = cc.corri(a, c, r, t, {}, riga, cc.SEME)
    q = {"linear": lineare(a, c), "error_rms_mdeg": cc.MDEG * rms(r - base["los"]), "i_max_circle_A": float(np.abs(base["i"][m]).max()), "scenarios": {}}
    for nome, sc in cc.SCENARI:
        if nome == "nominal":
            continue
        x = cc.corri(a, c, r, t, sc, riga, cc.SEME)
        ok = np.isfinite(x["los"][m]).all() and np.abs(x["los"][m]).max() < 10.0
        q["scenarios"][nome] = {"effect_rms_urad": 1e6 * rms(x["los"] - base["los"]) if ok else None, "error_rms_mdeg": cc.MDEG * rms(r - x["los"]) if ok else None}
    d = cm.prova(a, c, cm.p24.comandi_tipici(asse))
    q["doublets"] = {n: d[n] for n in ("err_rms_median_mdeg", "i_max_A", "over_max_mdeg", "t_set_median_s")}
    print(asse, et, round(q["error_rms_mdeg"], 2), q["doublets"], flush=True)
    return asse, et, q


def margini(a, c) -> dict:
    an, var = a.anello_v, g.varianti_cuscinetto(a)
    C, Bu = an.uscite["p"], an.B[:, 0]
    ritardo = np.exp(-2j * np.pi * F * g.PASSO_MASTER * g.RITARDO) * g.zoh(F)
    out = {}
    for n, A in var.items():
        P = g.frf(A, Bu, C, F)
        k = a.drive.p.Kpp if c is None else getattr(c, "k", None)
        L = k * ritardo * P if k is not None else g.anello_lqg(c, P, F)[0]
        out[n] = g.margini(L, F)
    return out


def caso_medio(arg) -> tuple:
    """Seni in banda media: ampiezza alla frequenza del seno, encoder e linea di vista, sulla differenza dal caso senza disturbi."""
    asse, et = arg
    a, C = TUTTI[asse]
    c = C[et]
    t = np.arange(0.0, 12.0, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    r = ref[asse] - ref[asse][0]
    riga = np.zeros(a.cascata.N)
    riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
    k = t >= 8.0
    base = cc.corri(a, c, r, t, {}, riga, cc.SEME)
    righe = []
    for am in AMP_MEDIE:
        for fs in F_MEDIE:
            x = cc.corri(a, c, r, t, {"sin": (am, fs)}, riga, cc.SEME)
            X = np.c_[np.sin(2 * np.pi * fs * t[k]), np.cos(2 * np.pi * fs * t[k])]
            ampiezza = lambda n: 1e6 * float(np.linalg.norm(np.linalg.lstsq(X, (x[n] - base[n])[k], rcond=None)[0]))
            righe.append({"f_Hz": fs, "amplitude_N_m": am, "encoder_urad": ampiezza("p"), "line_of_sight_urad": ampiezza("los")})
    q = {"sines": righe, "margins": margini(a, c)}
    print(asse, et, [(z["amplitude_N_m"], z["f_Hz"], round(z["line_of_sight_urad"], 1)) for z in righe], flush=True)
    return asse, et, q


def medie() -> None:
    compiti = [(asse, et) for asse, (_, C) in TUTTI.items() for et in C]
    out = {asse: {} for asse in TUTTI}
    with ProcessPoolExecutor(int(os.environ.get("LAVORI", "6"))) as ex:
        for asse, et, q in ex.map(caso_medio, compiti):
            out[asse][et] = q
    (cc.CARTELLA / "benefits_mid_band.json").write_text(json.dumps(cm.p24.senza_array(out), indent=1))


def raffica(t: np.ndarray, rms: float, seme: int) -> np.ndarray:
    """Coppia di raffica: rumore bianco sagomato in frequenza con la radice dello spettro di Davenport (5.14), a media nulla e con il rms dato."""
    bianco = np.fft.rfft(np.random.default_rng(seme).standard_normal(len(t)))
    om = 2 * np.pi * np.fft.rfftfreq(len(t), g.PASSO_SIM)
    beta = 600.0 / (np.pi * V_MEDIA)
    Sv = beta * om / (1 + (beta * om) ** 2) ** (4.0 / 3.0)
    x = np.fft.irfft(bianco * np.sqrt(Sv), len(t))
    return rms * (x - x.mean()) / x.std()


def caso_vento(arg) -> tuple:
    asse, et = arg
    a, C = TUTTI[asse]
    c = C[et] if et != "LQG-PI, k_i = 0.5 k_p" else CI[asse]
    t = np.arange(0.0, T_VENTO, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    r = ref[asse] - ref[asse][0]
    riga = np.zeros(a.cascata.N)
    riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
    m = t >= cc.T_AVVIO
    rms = lambda x: float(np.sqrt(np.mean(x[m] ** 2)))
    corsa = (lambda d: a.cascata.corri(r, 0.0, d, riga)) if c is None else (lambda d: g.corri_lqg(a, c, r, 0.0, g.RITARDO, d, riga))
    base = corsa(np.zeros_like(t))
    q = {"error_nominal_mdeg": cc.MDEG * rms(r - base["los"])}
    for liv in RMS_VENTO:
        for seme in SEMI_VENTO:
            x = corsa(raffica(t, liv, seme))
            q[f"{liv:g}_N_m_seed_{seme}"] = {"los_urad": 1e6 * rms(x["los"] - base["los"]), "encoder_urad": 1e6 * rms(x["p"] - base["p"]), "error_mdeg": cc.MDEG * rms(r - x["los"]), "i_max_A": float(np.abs(x["i"][m]).max())}
    print(asse, et, {k: round(v["los_urad"], 1) for k, v in q.items() if isinstance(v, dict)}, flush=True)
    return asse, et, q


def vento() -> None:
    compiti = [(asse, et) for asse, (_, C) in TUTTI.items() for et in C] + [(asse, "LQG-PI, k_i = 0.5 k_p") for asse in TUTTI]
    out = {asse: {} for asse in TUTTI}
    with ProcessPoolExecutor(int(os.environ.get("LAVORI", "7"))) as ex:
        for asse, et, q in ex.map(caso_vento, compiti):
            out[asse][et] = q
    t = np.arange(0.0, T_VENTO, g.PASSO_SIM)
    X = np.abs(np.fft.rfft(raffica(t, 1.0, SEMI_VENTO[0]))) ** 2
    fr = np.fft.rfftfreq(len(t), g.PASSO_SIM)
    out["gust_energy_below_Hz"] = {f"{x:g}": float(X[(fr > 0) & (fr <= x)].sum() / X[fr > 0].sum()) for x in (0.1, 1.0, 3.0, 10.0, 50.0)}
    (cc.CARTELLA / "benefits_wind.json").write_text(json.dumps(cm.p24.senza_array(out), indent=1))


def parti_kf() -> None:
    """Margine e picco di :math:`T_{ry}` con :math:`K_f` completo, solo sui poli del velocity loop (sotto 10 Hz), solo sui modi flessibili, nullo."""
    out = {}
    for asse, (a, C) in TUTTI.items():
        an = a.anello_v
        P = g.frf(g.varianti_cuscinetto(a)["sliding"], an.B[:, 0], an.uscite["p"], F)
        lenti = [i - 1 for m in a.modi if 0 < m.f_hz < 10 for i in m.indici]
        out[asse] = {}
        for et in ("LQG-PI tuned", "LQG-PI, q_p = 1e5"):
            c = C[et]
            vel = np.zeros_like(c.K_f)
            vel[lenti] = c.K_f[lenti]
            out[asse][et] = {"k_p": c.k_p}
            for lab, K in (("full", c.K_f), ("velocity loop only", vel), ("flexible only", c.K_f - vel), ("none", 0 * c.K_f)):
                L, N = g.anello_lqg(dataclasses.replace(c, K_f=K), P, F)
                out[asse][et][lab] = {"pm_deg": g.margini(L, F)["pm_deg"], "T_ry_peak_dB": float(pt._db(N / (1 + L)).max()), "band_Hz": g.banda(F, N / (1 + L))}
    (cc.CARTELLA / "benefits_kf.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out)[:600])


def main() -> None:
    if SOLO_VENTO:
        vento()
        return
    if SOLO_KF:
        parti_kf()
        return
    if SOLO_MEDIE:
        medie()
        return
    compiti = [(asse, et) for asse, (_, C) in TUTTI.items() for et in C]
    out = {asse: {} for asse in TUTTI}
    with ProcessPoolExecutor(int(os.environ.get("LAVORI", "6"))) as ex:
        for asse, et, q in ex.map(un_caso, compiti):
            out[asse][et] = q
    (cc.CARTELLA / "benefits.json").write_text(json.dumps(cm.p24.senza_array(out), indent=1))


if __name__ == "__main__":
    main()
