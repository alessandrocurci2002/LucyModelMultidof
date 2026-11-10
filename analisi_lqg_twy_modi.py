#!/usr/bin/env python3
r"""Perche' la :math:`T_{wy}` in elevazione non scende ai modi a 56 e 98 Hz, con :math:`q_f = 10^4` a 56 Hz.

Pesi di regolatore e stimatore separati sui modi (:math:`Q_f`, :math:`V_f`, p. 119 del libro): poli, zeri di :math:`F_y`, picchi di
:math:`T_{wy}` con cuscinetti in strisciamento e bloccati, e seni sul plant con Dahl, sulla differenza dal caso senza seno.
Solo diagnosi: il controllore tarato non cambia. Risultati in ``results_lqg_tuning/no_bearings/circle/``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]
SOLO_VF = "--vf" in sys.argv
SOLO_RIPOSO = "--riposo" in sys.argv
SOLO_ORBITA = "--orbita" in sys.argv
sys.argv = [sys.argv[0], "--senza-cuscinetti"]

import matplotlib.pyplot as plt  # noqa: E402
from scipy.linalg import eigvals  # noqa: E402

import analisi_lqg_circonferenza as cc  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_passo1_flessibili as p1  # noqa: E402
import analisi_lqg_punto_lavoro as apl  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402

F_56 = np.linspace(50.0, 64.0, 561)
F_98 = np.linspace(90.0, 110.0, 201)
F_SENI = (50.0, 54.0, 56.0, 57.0, 57.8, 58.5, 60.0, 98.0)
# (Q_f a 56 Hz, V_f a 56 Hz, Q_f a 98 Hz, V_f a 98 Hz)
VARIANTI = {"tuned, Q_f = V_f = 1e4": (1e4, 1e4, 0.1, 0.1), "Q_f = V_f = 0.1": (0.1, 0.1, 0.1, 0.1), "Q_f = 0.1, V_f = 1e6": (0.1, 1e6, 0.1, 0.1), "Q_f = 1e4, V_f = 1e6": (1e4, 1e6, 0.1, 0.1),
            "Q_f = 1e4, V_f = 1e8": (1e4, 1e8, 0.1, 0.1), "Q_f = 1e4, V_f = 1e12": (1e4, 1e12, 0.1, 0.1), "Q_f = V_f = 1e6": (1e6, 1e6, 0.1, 0.1), "98 Hz: Q_f = V_f = 1e6": (1e4, 1e4, 1e6, 1e6), "98 Hz: Q_f = 1e6, V_f = 1e10": (1e4, 1e4, 1e6, 1e10)}
IN_DAHL = ("cascade", "tuned, Q_f = V_f = 1e4", "Q_f = 1e4, V_f = 1e6")


def poli_modi(z: np.ndarray, lo: float, hi: float) -> list[list[float]]:
    s = np.log(np.asarray(z, complex)) / g.PASSO_MASTER
    s = s[s.imag > 0]
    return [[float(abs(x) / (2 * np.pi)), float(-x.real / abs(x))] for x in sorted(s, key=abs) if lo < abs(x) / (2 * np.pi) < hi]


def zeri_poli_fy(c) -> tuple[list, list]:
    """Zeri e poli di :math:`F_y`, senza ritardo e senza integrale, che sono lontani dai modi."""
    n = len(c.Ad)
    Kx = np.r_[c.k_p, c.K_f]
    Ac = c.Ad - np.outer(c.K_e[:, 0], np.eye(n)[0]) - np.outer(c.Bd[:, 0], Kx)
    M = np.block([[Ac, c.K_e[:, [0]]], [Kx[None, :], np.zeros((1, 1))]])
    E = np.zeros_like(M)
    E[:n, :n] = np.eye(n)
    z = eigvals(M, E)
    return poli_modi(z[np.isfinite(z)], 40.0, 110.0), poli_modi(np.linalg.eigvals(Ac), 40.0, 110.0)


def picco(f, T):
    i = int(np.argmax(np.abs(T)))
    return float(pt._db(T[i])), float(f[i])


def fig_twy(a, C, seg):
    an, var = a.anello_v, g.varianti_cuscinetto(a)
    Cy, Bu, Bw = an.uscite["p"], an.B[:, 0], an.B[:, 1]
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex="col")
    voci = {k: C[k] for k in ("Q_f = V_f = 0.1", "tuned, Q_f = V_f = 1e4", "Q_f = 1e4, V_f = 1e6")}
    for j, n in enumerate(("sliding", "stuck")):
        Pw, Pu = g.frf(var[n], Bw, Cy, F_56), g.frf(var[n], Bu, Cy, F_56)
        ax[0, j].plot(F_56, pt._db(Pw), lw=1.2, label="no position loop (= cascade)", **pt._stile(1))
        ax[1, j].plot(F_56, pt._db(0.0 * Pu + 1.0), lw=0.8, color=pt.INCHIOSTRO_SECONDARIO, ls=":")
        for k, (et, c) in enumerate(voci.items(), 2):
            L = g.anello_lqg(c, Pu, F_56)[0]
            ax[0, j].plot(F_56, pt._db(Pw / (1 + L)), lw=1.2, label=et, **pt._stile(k))
            ax[1, j].plot(F_56, pt._db(L), lw=1.2, label=et, **pt._stile(k))
        ax[0, j].set_title(f"Bearings {n}", fontsize=9)
        ax[1, j].set_xlabel("frequency [Hz]")
    # Punti sul plant con Dahl: seno piccolo con l'asse in moto a sinistra, con l'asse fermo a destra
    Mo, Ri = (json.loads((cc.CARTELLA / n).read_text()) for n in ("twy_modes_vf.json", "twy_modes_rest.json"))
    for et, m in (("cascade", 1), ("tuned, V_f = 1e4", 3), ("V_f = 1e6", 4)):
        for j, (righe, am, tit) in enumerate(((Mo[et]["sine"], 0.01, "moving axis, 0.01 N m"), (Ri[et], 0.002, "axis at rest, 0.002 N m"))):
            v = [r for r in righe if r["amplitude_N_m"] == am]
            ax[0, j].plot([r["f_Hz"] for r in v], [r["dB"] for r in v], ls="none", marker="o", ms=6, mfc="none", color=pt._stile(m)["color"], label=f"Dahl, {tit}: {et}")
    for q in ax[0]:
        q.set_ylim(-125.0, -55.0)
        q.legend(fontsize=7, loc="lower center", ncol=2)
    ax[0, 0].set_ylabel("T_wy, bearing torque to encoder [dB rad/(N m)]")
    ax[1, 0].set_ylabel("|L| [dB]")
    pt._titolo(fig, "Elevation, 56 Hz mode: disturbance rejection and loop gain", seg)
    return fig


def prova_vf(a, sint, asse) -> dict:
    """Ampiezza del seno sul plant con Dahl, rumore di misura, errore sulla circonferenza e doppietti, al variare del solo :math:`V_f` a 56 Hz."""
    t = np.arange(0.0, 12.0, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    r = ref[asse] - ref[asse][0]
    riga = np.zeros(a.cascata.N)
    riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
    k, m = t >= 8.0, t >= cc.T_AVVIO
    W = cc.cm.p24.comandi_tipici(asse)
    rms = lambda x: float(np.sqrt(np.mean(x[m] ** 2)))
    out = {}
    for et, c in {"cascade": None, "tuned, V_f = 1e4": sint(1e4, 1e4, 0.1, 0.1), "V_f = 1e5": sint(1e4, 1e5, 0.1, 0.1), "V_f = 1e6": sint(1e4, 1e6, 0.1, 0.1)}.items():
        base = cc.corri(a, c, r, t, {}, riga, cc.SEME)
        q = {"error_rms_mdeg": cc.MDEG * rms(r - base["los"]), "sine": []}
        for fs, am in ((57.0, 0.2), (57.8, 0.2), (57.0, 0.05), (57.8, 0.05), (57.0, 0.01), (57.8, 0.01), (57.0, 0.002), (57.8, 0.002)):
            x = cc.corri(a, c, r, t, {"sin": (am, fs)}, riga, cc.SEME)
            X = np.c_[np.sin(2 * np.pi * fs * t[k]), np.cos(2 * np.pi * fs * t[k])]
            A = float(np.linalg.norm(np.linalg.lstsq(X, (x["p"] - base["p"])[k], rcond=None)[0]))
            q["sine"].append({"f_Hz": fs, "amplitude_N_m": am, "dB": 20 * np.log10(A / am)})
        for sg in (1e-6, 1e-5):
            x = cc.corri(a, c, r, t, {"sigma": sg}, riga, cc.SEME)
            q[f"noise_{sg * 1e6:g}_urad_los_effect_urad"] = 1e6 * rms(x["los"] - base["los"])
        d = cc.cm.prova(a, c, W)
        q["doublets"] = {n: d[n] for n in ("err_rms_median_mdeg", "i_max_A", "over_max_mdeg", "t_set_median_s")}
        out[et] = q
        print(et, json.dumps(q)[:900], flush=True)
    (cc.CARTELLA / "twy_modes_vf.json").write_text(json.dumps(out, indent=1))
    return out


def a_riposo(a, sint, asse) -> dict:
    """Asse fermo, riferimento costante: con seni piccoli il cuscinetto resta bloccato. Il polo meno smorzato del tarato ha un tempo di
    assestamento di circa 3 s, per cui si misura dopo 9 s dall'inizio del seno."""
    t = np.arange(0.0, 20.0, g.PASSO_SIM)
    r = np.zeros_like(t)
    riga = np.zeros(a.cascata.N)
    riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
    k = t >= cc.T_AVVIO + 9.0
    out = {}
    for et, c in {"cascade": None, "tuned, V_f = 1e4": sint(1e4, 1e4, 0.1, 0.1), "V_f = 1e6": sint(1e4, 1e6, 0.1, 0.1)}.items():
        righe = []
        for fs in (56.0, 57.8, 58.5):
            for am in (0.002, 0.01, 0.05, 0.2):
                x = cc.corri(a, c, r, t, {"sin": (am, fs)}, riga, cc.SEME)
                X = np.c_[np.sin(2 * np.pi * fs * t[k]), np.cos(2 * np.pi * fs * t[k])]
                A = float(np.linalg.norm(np.linalg.lstsq(X, x["p"][k], rcond=None)[0]))
                L = float(np.linalg.norm(np.linalg.lstsq(X, x["los"][k], rcond=None)[0]))
                righe.append({"f_Hz": fs, "amplitude_N_m": am, "dB": 20 * np.log10(A / am), "encoder_urad": 1e6 * A, "line_of_sight_urad": 1e6 * L})
        out[et] = righe
        print(et, [(q["f_Hz"], q["amplitude_N_m"], round(q["dB"], 1), round(q["encoder_urad"], 2)) for q in righe], flush=True)
    (cc.CARTELLA / "twy_modes_rest.json").write_text(json.dumps(out, indent=1))
    return out


def orbita(a, sint, asse) -> dict:
    """Orbita intera, comprese le inversioni di velocita', per il solo asse: errore e effetto del rumore sulla linea di vista."""
    t = np.arange(0.0, cc.T_FINE, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    r = ref[asse] - ref[asse][0]
    riga = np.zeros(a.cascata.N)
    riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
    m = t >= cc.T_AVVIO
    rms = lambda x: float(np.sqrt(np.mean(x[m] ** 2)))
    out = {}
    for et, c in {"cascade": None, "tuned, V_f = 1e4": sint(1e4, 1e4, 0.1, 0.1), "V_f = 1e6": sint(1e4, 1e6, 0.1, 0.1)}.items():
        base = cc.corri(a, c, r, t, {}, riga, cc.SEME)
        q = {"error_rms_mdeg": cc.MDEG * rms(r - base["los"])}
        for sg in (1e-6, 1e-5):
            x = cc.corri(a, c, r, t, {"sigma": sg}, riga, cc.SEME)
            q[f"noise_{sg * 1e6:g}_urad_los_effect_urad"] = 1e6 * rms(x["los"] - base["los"])
        out[et] = q
        print(et, q, flush=True)
    (cc.CARTELLA / "twy_modes_orbit.json").write_text(json.dumps(out, indent=1))
    return out


def main() -> None:
    par, Kt, S, T, mot, assi = g.carica()
    Fi = json.loads(cc.TUNING.read_text())
    asse, ch, cat, mod = next(x for x in assi if x[0] == "polare")
    a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
    fi = Fi[ch]["final"]
    m56 = next(m for m in a.modi if len(m.indici) == 2 and m.zeta < 0.1)
    m98 = next(m for m in a.modi if abs(m.f_hz - 98.25) < 1.0)
    i56, i98 = [i - 1 for i in m56.indici], [i - 1 for i in m98.indici]
    an, var = a.anello_v, g.varianti_cuscinetto(a)
    Cy, Bu, Bw = an.uscite["p"], an.B[:, 0], an.B[:, 1]

    def sint(qf, vf, q98, v98):
        q, v = np.full(len(a.Ap) - 1, 0.1), np.full(len(a.Ap) - 1, 0.1)
        q[i56], v[i56], q[i98], v[i98] = qf, vf, q98, v98
        return lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, fi["q_i"], fi["q_p"], q, v_f=v, v_continua=True)

    if SOLO_VF:
        prova_vf(a, sint, asse)
        return
    if SOLO_RIPOSO:
        a_riposo(a, sint, asse)
        return
    if SOLO_ORBITA:
        orbita(a, sint, asse)
        return
    C = {k: sint(*x) for k, x in VARIANTI.items()}
    f = np.r_[F_56, F_98]
    nf = len(F_56)
    out = {"model_modes": [{"f_Hz": m.f_hz, "zeta": m.zeta, "hinf": m.hinf} for m in a.modi if np.isfinite(m.hinf)], "plant": {}, "cases": {}}
    PU, PW = {n: g.frf(A, Bu, Cy, f) for n, A in var.items()}, {n: g.frf(A, Bw, Cy, f) for n, A in var.items()}
    for n, A in var.items():
        ev = np.linalg.eigvals(A)
        ev = ev[ev.imag > 0]
        out["plant"][n] = {"poles": [[float(abs(x) / (2 * np.pi)), float(-x.real / abs(x))] for x in sorted(ev, key=abs) if abs(x) / (2 * np.pi) < 110], "peak_56_dB": picco(F_56, PW[n][:nf]), "peak_98_dB": picco(F_98, PW[n][nf:])}
    fr = np.geomspace(0.1, 150.0, 600)
    PUr = g.frf(var["sliding"], Bu, Cy, fr)
    for et, c in C.items():
        r = {"regulator": poli_modi(c.poli_regolatore, 40.0, 110.0), "estimator": poli_modi(c.poli_stimatore, 40.0, 400.0)}
        r["F_y_zeros"], r["F_y_poles"] = zeri_poli_fy(c)
        r["F_y_56_Hz"] = float(abs(c.funzione_trasferimento(np.exp(2j * np.pi * np.array([56.12]) * g.PASSO_MASTER), g.RITARDO)[1][0]))
        Lr, Tr = g.anello_lqg(c, PUr, fr)
        r["T_ry_peak_40_70_dB"] = float(pt._db(Tr / (1 + Lr))[(fr > 40) & (fr < 70)].max())
        for n, A in var.items():
            L = g.anello_lqg(c, PU[n], f)[0]
            Tw = PW[n] / (1 + L)
            z = p1.poli_chiusi(c, A, Bu, Cy)
            r[n] = {"largest_z": float(np.abs(z).max()), "closed_loop": poli_modi(z, 40.0, 110.0), "T_wy_56_dB": picco(F_56, Tw[:nf]), "T_wy_98_dB": picco(F_98, Tw[nf:]), "L_peak_56": float(np.abs(L[:nf]).max()), "L_peak_98": float(np.abs(L[nf:]).max())}
        out["cases"][et] = r
        print(et, r["regulator"][1], r["estimator"][1:2], r["sliding"]["T_wy_56_dB"], r["stuck"]["T_wy_56_dB"], r["stuck"]["largest_z"], flush=True)

    # Seni sul plant con Dahl, sulla circonferenza
    t = np.arange(0.0, 12.0, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    r0_ = ref[asse] - ref[asse][0]
    riga = np.zeros(a.cascata.N)
    riga[a.m.nodi.index(cc.NODO_LOS[asse])] = 1.0
    k = t >= 8.0
    seni = {}
    for et in IN_DAHL:
        c = None if et == "cascade" else C[et]
        base = cc.corri(a, c, r0_, t, {}, riga, cc.SEME)
        righe = []
        for fs in F_SENI:
            x = cc.corri(a, c, r0_, t, {"sin": (cc.AMP_SENI, fs)}, riga, cc.SEME)
            X = np.c_[np.sin(2 * np.pi * fs * t[k]), np.cos(2 * np.pi * fs * t[k])]
            amp = lambda n: 1e6 * float(np.linalg.norm(np.linalg.lstsq(X, (x[n] - base[n])[k], rcond=None)[0]))
            righe.append({"f_Hz": fs, "encoder_urad": amp("p"), "line_of_sight_urad": amp("los"), "dB": 20 * np.log10(1e-6 * amp("p") / cc.AMP_SENI)})
        seni[et] = righe
        print(et, [(q["f_Hz"], round(q["encoder_urad"], 2), round(q["dB"], 1)) for q in righe], flush=True)
    out["sines_on_Dahl"] = seni
    (cc.CARTELLA / "twy_modes.json").write_text(json.dumps(out, indent=1))
    pt.applica_stile()
    fig_twy(a, C, par.segnaposto).savefig(cc.CARTELLA / "figures" / "twy_modes.pdf")


if __name__ == "__main__":
    main()
