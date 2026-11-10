"""Inseguimento della circonferenza di Attorrese con l'LQG-PI tarato, confrontato con la cascata sul plant Dahl, proxy della realta'.

Riferimento e metriche dello studio di Hilkert (``closedLoopModel/dob_tracking_study.py``): raggio 10 m, quota 5 m, centro a 100 m, 0.2 rad/s; 5 s di
avvio esclusi, un'orbita; errore di puntamento della linea di vista, jitter sopra 1 e 10 Hz. Modello interno senza cuscinetti, simulatore Dahl.
Poi disturbi di coppia sulla porta dei cuscinetti e rumore sulla misura di posizione. I livelli sono assunti: la misura a riposo non da' il rumore.
Il rumore non entra nella velocita' del drive, uguale per tutti i controllori. I due assi si scambiano a ogni passo le coppie di Coriolis e giroscopiche delle
equazioni di Lagrange (``CoppieLagrange``, sui nodi del telescopio), e l'inerzia del telescopio attorno all'azimut dipende dall'elevazione. Nessun limite di corrente nella simulazione. L'errore totale e' dominato dall'inseguimento in azimut, quindi l'effetto di un disturbo si misura come differenza
dal caso senza disturbi, a parita' di controllore.

    .venv/bin/python -I analisi_lqg_circonferenza.py [--extra | --dither | --senza-extra]

``--extra`` ricalcola solo la reiezione lineare e i doppietti delle varianti; ``--dither`` solo la separazione in banda dell'effetto del seno a 56 Hz.

Risultati in results_lqg_tuning/no_bearings/circle/.
"""

from __future__ import annotations

import dataclasses
import json
import math
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]

import matplotlib.pyplot as plt  # noqa: E402
from scipy.linalg import expm  # noqa: E402
from scipy.signal import butter, filtfilt, sosfiltfilt  # noqa: E402

import analisi_lqg_confronto_misure as cm  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402
import analisi_lqg_punto_lavoro as apl  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning" / "no_bearings" / "circle"
TUNING = RADICE / "results_lqg_tuning" / "no_bearings" / "steps_2_4_measured_error" / "summary.json"
"""Tuning dell'LQG-PI della tesi. Dal 2026-11-05 e' quello con il PI sull'errore misurato (Gawronski Fig. 8.4b,
``lqg_pi_errore_misurato/``); il precedente, con il PI sull'errore stimato, e' in ``steps_2_4/summary.json``."""
PI_MISURATO = True
"""``casi``: l'LQG-PI tarato con il PI sull'errore misurato (Fig. 8.4b) invece che sull'errore stimato (8.15)."""
T_AVVIO, MDEG, SEME = apl.T_AVVIO, 1e3 * 180 / np.pi, 1
T_FINE = T_AVVIO + 2 * np.pi / apl.OMEGA
NODO_LOS = {"azimut": "AP", "polare": "hp"}
"""Il nodo che porta il telescopio, rigido sull'anello."""
T_GRADINO = 15.0
PASSO_THETA = np.deg2rad(0.1)
THETA_POLI = (0.0, 3.0, 30.0, 45.0, 60.0, 80.0)
FREQ_SENI = (0.3, 1.0, 3.0, 10.0, 30.0, 56.0, 100.0)
AMP_SENI = 0.2
SCALA_QF_AZ = (1e4, 3e4, 1e5)
F_ANELLO = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 56.0])
TUTTI_DIST = "all: 1 urad, step, Gaussian with mean"
SCENARI = (("nominal", {}), ("noise 1 urad", {"sigma": 1e-6}), ("noise 10 urad", {"sigma": 1e-5}), ("torque step 0.5 N m", {"step": 0.5}), ("Gaussian 0.1 N m rms", {"gauss": 0.1}),
           ("Gaussian 0.1 N m rms, mean 0.2 N m", {"gauss": 0.1, "media": 0.2}), ("sine 0.2 N m, 1 Hz", {"sin": (0.2, 1.0)}), ("sine 0.2 N m, 10 Hz", {"sin": (0.2, 10.0)}), ("sine 0.2 N m, 56 Hz", {"sin": (0.2, 56.0)}),
           (TUTTI_DIST, {"sigma": 1e-6, "step": 0.5, "gauss": 0.1, "media": 0.2}))
"""La media (cavi, sbilanciamento) c'e' gia' prima dell'inseguimento: la co-simulazione parte nell'equilibrio in cui l'integratore di velocita' la regge."""


class CoppieLagrange:
    """Coppie dipendenti dalle velocita' scambiate fra gli assi dal telescopio, dalle equazioni di Lagrange di ``gimbal_mdof.coupled``.

    Con omega = (psi' sin th, -th', psi' cos th) e il tensore I rispetto a O, T = 1/2 J(th) psi'^2 + 1/2 I_yy th'^2 + H(th) psi' th', dove
    J = I_xx sin^2 th + I_zz cos^2 th + 2 I_xz sin th cos th e H = -(I_xy sin th + I_yz cos th). Da qui, a secondo membro:
    azimut -(J' psi' th' + H' th'^2), elevazione +1/2 J' psi'^2. I termini inerziali -H th'' sull'azimut e -H psi'' sull'elevazione (H ~ 1e-6) sono in ``H``.
    ``scala``: fattore CAD stimato sul carico di azimut, che il banco applica gia' a J(0) = I_zz; va applicato a tutto il tensore perche' energia e coppie restino coerenti."""

    def __init__(self, pl, scala: float = 1.0):
        I = scala * pl.I
        self.Ixx, self.Izz = float(I[0, 0]), float(I[2, 2])
        self.dIxz, self.Ixz, self.Ixy, self.Iyz = float(I[0, 0] - I[2, 2]), float(I[0, 2]), float(I[0, 1]), float(I[1, 2])

    def H(self, th: float) -> float:
        """Elemento fuori diagonale della matrice di massa del telescopio."""
        return -(self.Ixy * math.sin(th) + self.Iyz * math.cos(th))

    def delta_J(self, th: float) -> float:
        """J(th) - J(0): quanto cambia l'inerzia del telescopio attorno all'azimut rispetto al banco, identificato a th = 0."""
        sn, cs = math.sin(th), math.cos(th)
        return self.Ixx * sn * sn + self.Izz * (cs * cs - 1.0) + 2 * self.Ixz * sn * cs

    def coppie(self, th: float, dpsi: float, dth: float) -> tuple[float, float]:
        sn, cs = math.sin(th), math.cos(th)
        dJ = self.dIxz * 2 * sn * cs + 2 * self.Ixz * (cs * cs - sn * sn)
        dH = -(self.Ixy * cs - self.Iyz * sn)
        return -(dJ * dpsi * dth + dH * dth * dth), 0.5 * dJ * dpsi * dpsi


def verifica_coppie(par) -> dict:
    """Le forme chiuse contro le forze del modello simbolico di Lagrange, con moto rigido dei due assi, senza dissipazione e gravita' sottratta."""
    sy = afm.sistema(par)
    acc = sy.accoppiato.senza_dissipazione()
    st = acc.struttura
    cl = CoppieLagrange(sy.payload)
    errore, righe = 0.0, []
    for th, dp, dt_ in ((0.05, 0.02, 0.001), (0.3, 0.2, 0.05), (0.8, 0.5, -0.2), (-0.6, -0.3, 0.1)):
        q = np.zeros(st.n)
        q[st.n_psi:] = th
        dq = np.zeros(st.n)
        dq[:st.n_psi], dq[st.n_psi:] = dp, dt_
        f = np.asarray(acc.funzioni.forze(q, dq, acc.parametri), dtype=float).reshape(-1) - np.asarray(acc.funzioni.forze(q, 0 * dq, acc.parametri), dtype=float).reshape(-1)
        ta, te = cl.coppie(th, dp, dt_)
        fa, fe = float(-f[:st.n_psi].sum()), float(-f[st.n_psi:].sum())
        errore = max(errore, abs(fa - ta), abs(fe - te))
        righe.append({"theta_rad": th, "psi_rate": dp, "theta_rate": dt_, "azimuth_closed_form": ta, "azimuth_lagrange": fa, "elevation_closed_form": te, "elevation_lagrange": fe})
    if errore > 1e-9:
        raise AssertionError(f"coppie diverse dal modello di Lagrange: {errore:.3e}")
    # J(th) contro la matrice di massa simbolica, nodo del telescopio in azimut
    j = st.j_psi_p
    z = np.zeros(st.n)
    M0 = float(acc.M(z)[j, j])
    err_J = 0.0
    for th in (0.05, 0.5, 1.0, -0.8):
        q = z.copy()
        q[st.j_theta_p] = th
        err_J = max(err_J, abs(float(acc.M(q)[j, j]) - M0 - cl.delta_J(th)))
    if err_J > 1e-9:
        raise AssertionError(f"J(theta) diversa dal modello di Lagrange: {err_J:.3e}")
    # H(th) contro lo stesso elemento della matrice di massa
    err_H = 0.0
    for th in (0.0, 0.05, 0.5, 1.0, -0.8):
        q = z.copy()
        q[st.j_theta_p] = th
        err_H = max(err_H, abs(float(acc.M(q)[j, st.j_theta_p]) - cl.H(th)))
    if err_H > 1e-12:
        raise AssertionError(f"H(theta) diversa dal modello di Lagrange: {err_H:.3e}")
    return {"largest_difference_N_m": errore, "largest_difference_J_kg_m2": err_J, "largest_difference_H_kg_m2": err_H, "rows": righe}


def simula_coppia(Ax, C, ref, t, sc, righe, seme, cl, th0, inerzia: bool = True, con_H: bool = True) -> dict:
    """Co-simulazione dei due assi, ciascuno come ``AnelloChiuso.corri`` o ``g.corri_lqg``, con le coppie ``cl`` (``None``: assi indipendenti) calcolate a inizio passo
    dalle velocita' dei nodi del telescopio; esse entrano in quei nodi, come in Lagrange. Con ``inerzia`` l'inerzia del nodo del telescopio in azimut segue J(theta):
    cambia solo la riga della sua accelerazione, e il modello discreto si ricalcola a passi di ``PASSO_THETA``. M e' diagonale in azimut; lo smorzamento modale resta quello di theta = 0.
    Con ``con_H`` si aggiungono le coppie -H th'' e -H psii'', con l'accelerazione dell'altro nodo sul passo precedente: H^2 / (J I_yy) ~ 1e-12, quindi il ritardo di un passo non conta."""
    assi_ = ("azimut", "polare")
    n_t = len(t)
    d = {}
    for i, q in enumerate(assi_):
        a, c = Ax[q], C[q]
        an = a.cascata if c is None else a.anello_v
        Mi = np.linalg.inv(a.m.M)
        bc = np.zeros(an.N)
        bc[an.n:2 * an.n] = Mi @ righe[q][:an.n]
        B3 = np.c_[an.B, bc]
        E = expm(np.block([[an.A, B3], [np.zeros((3, an.N + 3))]]) * g.PASSO_SIM)
        sigma = sc.get("sigma", 0.0)
        r = ref[q] - ref[q][0]
        if c is None and sigma:
            n = np.random.default_rng(seme + i).standard_normal(int(n_t * g.PASSO_SIM / g.PASSO_MASTER) + 2) * sigma
            r = r - n[np.floor(t / g.PASSO_MASTER).astype(int)]
        x = None if c is None else (cm.ConRumore(c, sigma, np.random.default_rng(seme + i)) if sigma else c)
        # errore di modello: livello di Coulomb del plant scalato, sigma0 e controllori invariati
        att = dataclasses.replace(a.att, Tc=a.att.Tc * sc["attrito"]) if "attrito" in sc else a.att
        d[q] = {"an": an, "Ad": E[:an.N, :an.N], "Bd": E[:an.N, an.N:], "s": a.cascata.equilibrio(sc.get("media", 0.0)), "z": 0.0, "s0": a.att.sigma0, "att": att, "r": r, "dist": disturbo(sc, t, seme + i), "c": x, "st": x.stato_iniziale() if x else None,
                "coda": [0.0] * g.RITARDO, "u": 0.0, "ys": 0.0, "vv": np.r_[np.zeros(an.n), a.m.v, np.zeros(an.N - 2 * an.n)], "rp": an.uscite["p"], "ri": an.uscite["i"], "rl": np.r_[righe[q][:an.n], np.zeros(an.N - an.n)],
                "rv": np.r_[np.zeros(an.n), righe[q][:an.n], np.zeros(an.N - 2 * an.n)], "P": np.empty(n_t), "I": np.empty(n_t), "L": np.empty(n_t), "T": np.empty(n_t)}
    sotto = int(round(g.PASSO_MASTER / g.PASSO_SIM))
    Da, De = d["azimut"], d["polare"]
    a_az = Ax["azimut"]
    an_az = Da["an"]
    i_v = an_az.n + int(np.argmax(righe["azimut"][:an_az.n]))
    J0 = float(a_az.m.M[i_v - an_az.n, i_v - an_az.n])
    bc0 = np.zeros(an_az.N)
    bc0[an_az.n:2 * an_az.n] = np.linalg.inv(a_az.m.M) @ righe["azimut"][:an_az.n]
    cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    def modello(th: float) -> tuple[np.ndarray, np.ndarray]:
        key = int(round(th / PASSO_THETA))
        if key not in cache:
            fatt = J0 / (J0 + cl.delta_J(key * PASSO_THETA))
            A_ = an_az.A.copy()
            B_ = np.c_[an_az.B, bc0]
            A_[i_v, :] *= fatt
            B_[i_v, :] *= fatt
            E_ = expm(np.block([[A_, B_], [np.zeros((3, an_az.N + 3))]]) * g.PASSO_SIM)
            cache[key] = (E_[:an_az.N, :an_az.N], E_[:an_az.N, an_az.N:])
        return cache[key]

    v_a0 = v_e0 = 0.0
    for k in range(n_t):
        if cl is not None:
            th_k = th0 + float(De["rl"] @ De["s"])
            v_a, v_e = float(Da["rv"] @ Da["s"]), float(De["rv"] @ De["s"])
            ta, te = cl.coppie(th_k, v_a, v_e)
            if con_H:
                h = cl.H(th_k)
                ta -= h * (v_e - v_e0) / g.PASSO_SIM
                te -= h * (v_a - v_a0) / g.PASSO_SIM
                v_a0, v_e0 = v_a, v_e
            if inerzia:
                Da["Ad"], Da["Bd"] = modello(th_k)
        else:
            ta, te = 0.0, 0.0
        Da["T"][k], De["T"][k] = ta, te
        for D, tq in ((Da, ta), (De, te)):
            s = D["s"]
            y = float(D["rp"] @ s)
            if D["c"] is None:
                u0 = D["r"][k]
            else:
                if k % sotto == 0:
                    D["ys"] = y
                    if hasattr(D["c"], "misura"):  # aggancio: un controllore con ``misura`` (PI sull'errore misurato, Gawronski Fig. 8.4b) riceve y(k) prima di ``uscita``
                        D["c"].misura(D["st"], y)
                    D["coda"].append(D["c"].uscita(D["st"], float(D["r"][k])))
                    D["u"] = D["coda"].pop(0)
                u0 = D["u"]
            ing = np.array([u0, D["s0"] * D["z"] + D["dist"][k], tq])
            D["P"][k], D["I"][k], D["L"][k] = y, float(D["ri"] @ s), float(D["rl"] @ s)
            s_new = D["Ad"] @ s + D["Bd"] @ ing
            vd = 0.5 * float(D["vv"] @ s + D["vv"] @ s_new)
            al = D["s0"] * abs(vd) / D["att"].g(vd)
            if al > 0:
                e = math.exp(-al * g.PASSO_SIM)
                D["z"] = D["z"] * e + vd / al * (1 - e)
            D["s"] = s_new
            if D["c"] is not None and k % sotto == sotto - 1:
                D["c"].stima(D["st"], D["ys"], D["u"])
    return {q: {"p": d[q]["P"], "i": d[q]["I"], "los": d[q]["L"], "coupling_torque": d[q]["T"]} for q in assi_}


def stabilita_theta(Ax, Cc, cl) -> list[dict]:
    """Polo piu' grande dell'anello in azimut con J(theta), controllore fisso, identificato a theta = 0; cuscinetti in strisciamento e bloccati."""
    a = Ax["azimut"]
    an = a.anello_v
    casc = a.cascata
    i = int(np.argmax(np.array([1.0 if n == NODO_LOS["azimut"] else 0.0 for n in a.m.nodi])))
    J0 = float(a.m.M[i, i])
    righe_ = []
    for gr in THETA_POLI:
        f = J0 / (J0 + cl.delta_J(np.deg2rad(gr)))
        r = {"theta_deg": gr, "J_ratio": 1 / f}
        for et, c in Cc["azimut"].items():
            if et not in ("cascade", "LQG-PI tuned", "LQG-PI, q_p = 1e5"):
                continue
            for n, A in g.varianti_cuscinetto(a).items():
                A_ = A.copy()
                A_[an.n + i, :] *= f
                B_ = an.B[:, 0].copy()
                B_[an.n + i] *= f
                if c is None:
                    Ac = casc.A.copy()
                    if n == "stuck":
                        Ac = Ac + a.att.sigma0 * np.outer(casc.B[:, 1], np.r_[a.m.v, np.zeros(casc.N - casc.n)])
                    Ac[casc.n + i, :] *= f
                    r[f"{et}, {n}"] = float(np.exp(np.linalg.eigvals(Ac).real.max() * g.PASSO_MASTER))
                else:
                    r[f"{et}, {n}"] = float(np.abs(cm.p1.poli_chiusi(c, A_, B_, an.uscite["p"])).max())
        righe_.append(r)
        print("theta", gr, {k: round(v, 5) for k, v in r.items()}, flush=True)
    return righe_


def prova_coppie(Ax, Cc, righe, cl) -> dict:
    """Effetto del cross coupling: la circonferenza del test e una piu' veloce e alta, dove i termini sono piu' grandi. Solo coppie, e coppie con J(theta)."""
    out = {}
    for nome, (ce, ra, qo, om, tf) in {"test circle": (apl.CENTRO, apl.RAGGIO, apl.QUOTA, apl.OMEGA, T_FINE), "fast, high circle (centre 15 m, height 15 m, 0.5 rad/s)": (15.0, 10.0, 15.0, 0.5, 18.0)}.items():
        t = np.arange(0.0, tf, g.PASSO_SIM)
        x, y = ce + ra * np.cos(om * t), ra * np.sin(om * t)
        ref = {"azimut": np.unwrap(np.arctan2(y, x)), "polare": np.arctan2(qo, np.hypot(x, y))}
        m = t >= T_AVVIO
        out[nome] = {"psi_rate_max": float(np.abs(np.gradient(ref["azimut"], g.PASSO_SIM)).max()), "theta_rate_max": float(np.abs(np.gradient(ref["polare"], g.PASSO_SIM)).max()), "theta_mean_deg": float(np.rad2deg(ref["polare"].mean()))}
        for et in ("cascade", "LQG-PI tuned"):
            C = {q: Cc[q][et] for q in Ax}
            r0 = simula_coppia(Ax, C, ref, t, {}, righe, SEME, None, float(ref["polare"][0]))
            r2 = simula_coppia(Ax, C, ref, t, {}, righe, SEME, cl, float(ref["polare"][0]), False, False)
            r3 = simula_coppia(Ax, C, ref, t, {}, righe, SEME, cl, float(ref["polare"][0]), True, False)
            r1 = simula_coppia(Ax, C, ref, t, {}, righe, SEME, cl, float(ref["polare"][0]))
            rms = lambda v: float(np.sqrt(np.mean(v[m] ** 2)))
            out[nome][et] = {"torque_az_max_N_m": float(np.abs(r1["azimut"]["coupling_torque"][m]).max()), "torque_el_max_N_m": float(np.abs(r1["polare"]["coupling_torque"][m]).max()),
                             "torques_only_los_az_rms_urad": 1e6 * rms(r2["azimut"]["los"] - r0["azimut"]["los"]), "torques_only_los_el_rms_urad": 1e6 * rms(r2["polare"]["los"] - r0["polare"]["los"]),
                             "J_only_los_az_rms_urad": 1e6 * rms(r3["azimut"]["los"] - r2["azimut"]["los"]), "J_only_los_el_rms_urad": 1e6 * rms(r3["polare"]["los"] - r2["polare"]["los"]),
                             "H_only_los_az_rms_urad": 1e6 * rms(r1["azimut"]["los"] - r3["azimut"]["los"]), "H_only_los_el_rms_urad": 1e6 * rms(r1["polare"]["los"] - r3["polare"]["los"]),
                             "i_az_max_A": float(np.abs(r1["azimut"]["i"][m]).max()), "i_az_max_fixed_J_A": float(np.abs(r2["azimut"]["i"][m]).max()),
                             "effect_los_az_rms_urad": 1e6 * rms(r1["azimut"]["los"] - r0["azimut"]["los"]), "effect_los_el_rms_urad": 1e6 * rms(r1["polare"]["los"] - r0["polare"]["los"]),
                             "effect_los_az_peak_urad": 1e6 * float(np.abs((r1["azimut"]["los"] - r0["azimut"]["los"])[m]).max()), "effect_los_el_peak_urad": 1e6 * float(np.abs((r1["polare"]["los"] - r0["polare"]["los"])[m]).max()),
                             "error_az_rms_mdeg": MDEG * rms(ref["azimut"] - ref["azimut"][0] - r1["azimut"]["los"]), "error_el_rms_mdeg": MDEG * rms(ref["polare"] - ref["polare"][0] - r1["polare"]["los"])}
            print(nome, et, {k: round(v, 5) for k, v in out[nome][et].items()}, flush=True)
    return out


def disturbo(sc: dict, t: np.ndarray, seme: int) -> np.ndarray:
    d = np.zeros_like(t)
    if "step" in sc:
        d += sc["step"] * (t >= T_GRADINO)
    if "sin" in sc:
        d += sc["sin"][0] * np.sin(2 * np.pi * sc["sin"][1] * t) * (t >= T_AVVIO)
    if "proc" in sc:
        n = np.random.default_rng(seme).standard_normal(len(t))
        n = sosfiltfilt(butter(4, (0.1, 50.0), "bandpass", fs=1 / g.PASSO_SIM, output="sos"), n)
        d += sc["proc"] * n / n.std()
    if "gauss" in sc:
        # Margine scartato: il filtro a 0.1 Hz lascia transitori ai bordi per decine di secondi, e il segnale non sarebbe gaussiano
        pad = int(100.0 / g.PASSO_SIM)
        n = np.random.default_rng(seme).standard_normal(len(t) + 2 * pad)
        n = sosfiltfilt(butter(4, sc.get("banda", (0.1, 50.0)), "bandpass", fs=1 / g.PASSO_SIM, output="sos"), n)[pad:pad + len(t)]
        d += sc["gauss"] * (n - n.mean()) / n.std()
    if "media" in sc:
        # Parte statica (cavi, sbilanciamento): costante dall'inizio della simulazione
        d += sc["media"]
    return d


def corri(a, c, r, t, sc, riga, seme):
    """Riferimento relativo all'inizio; il rumore entra nella misura di posizione del solo anello esterno, a 1 kHz."""
    d = disturbo(sc, t, seme)
    sigma = sc.get("sigma", 0.0)
    if c is None:
        n = np.random.default_rng(seme).standard_normal(int(len(t) * g.PASSO_SIM / g.PASSO_MASTER) + 2) * sigma
        pc = r - (n[np.floor(t / g.PASSO_MASTER).astype(int)] if sigma else 0.0)
        return a.cascata.corri(pc, 0.0, d, riga)
    x = cm.ConRumore(c, sigma, np.random.default_rng(seme)) if sigma else c
    return g.corri_lqg(a, x, r, 0.0, g.RITARDO, d, riga)


def errore_puntamento(t, psi, th) -> np.ndarray:
    """Angolo fra la linea di vista e la direzione del drone, come ``pointing_error`` di dob_tracking_study."""
    x, y = apl.CENTRO + apl.RAGGIO * np.cos(apl.OMEGA * t), apl.RAGGIO * np.sin(apl.OMEGA * t)
    p = np.array([x, y, np.full_like(x, apl.QUOTA)])
    p = p / np.linalg.norm(p, axis=0)
    n = np.array([np.cos(th) * np.cos(psi), np.cos(th) * np.sin(psi), np.sin(th)])
    return np.arcsin(np.clip(np.linalg.norm(np.cross(n.T, p.T), axis=1), 0.0, 1.0))


def metriche(t, ref, R) -> dict:
    """R: uscite dei due assi. Errore di puntamento della linea di vista, jitter, encoder e corrente."""
    m = t >= T_AVVIO
    psi0, th0 = ref["azimut"][0], ref["polare"][0]
    psi, th = psi0 + R["azimut"]["los"], th0 + R["polare"]["los"]
    E = errore_puntamento(t, psi, th)
    ex, ee = (ref["azimut"] - psi) * np.cos(ref["polare"]), ref["polare"] - th
    fs = 1 / (t[1] - t[0])
    out = {"E_rms_mdeg": MDEG * float(np.sqrt(np.mean(E[m] ** 2))), "E_max_mdeg": MDEG * float(E[m].max()), "E_p95_mdeg": MDEG * float(np.percentile(E[m], 95)), "E_mean_mdeg": MDEG * float(E[m].mean()),
           "miss_rms_mm": 1e3 * float(np.sqrt(np.mean((E[m] * 100.0) ** 2))), "xel_rms_mdeg": MDEG * float(np.sqrt(np.mean(ex[m] ** 2))), "el_rms_mdeg": MDEG * float(np.sqrt(np.mean(ee[m] ** 2)))}
    for f in (1.0, 10.0):
        b, a = butter(2, f / (fs / 2), "high")
        hx, he = filtfilt(b, a, ex), filtfilt(b, a, ee)
        out[f"jitter_{f:g}Hz_urad"] = 1e6 * float(np.sqrt(np.mean(hx[m] ** 2 + he[m] ** 2)))
    for ax, q in (("az", "azimut"), ("el", "polare")):
        out[f"i_rms_{ax}_A"] = float(np.sqrt(np.mean(R[q]["i"][m] ** 2)))
        out[f"i_max_{ax}_A"] = float(np.abs(R[q]["i"][m]).max())
        out[f"enc_{ax}_rms_mdeg"] = MDEG * float(np.sqrt(np.mean((ref[q] - ref[q][0] - R[q]["p"])[m] ** 2)))
        out[f"los_{ax}_rms_mdeg"] = MDEG * float(np.sqrt(np.mean((ref[q] - ref[q][0] - R[q]["los"])[m] ** 2)))
    return out


def casi(a, fi, f_modo) -> dict:
    """Tarato, integrale al valore tipico del libro (k_i = 0.5 k_p), V piu' grande, e in elevazione un guadagno minore."""
    def sint(qp, qi, rho=1.0):
        q = np.full(len(a.Ap) - 1, 0.1)
        if fi["q_f"]:
            m = next(x for x in a.modi if abs(x.f_hz - f_modo) < 3.0)
            q[[i - 1 for i in m.indici]] = float(fi["q_f"]["0"])
        return lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, qi, qp, q, rho_e=rho, v_continua=True)
    tarato = sint(fi["q_p"], fi["q_i"])
    return {"cascade": None, "LQG-PI tuned": lqg.pi_misurato(tarato) if PI_MISURATO else tarato, "LQG-PI, k_i = 0.5 k_p": sint(fi["q_p"], (0.5 * fi["k_p"]) ** 2), "LQG-PI, V = 100 Q": sint(fi["q_p"], fi["q_i"], 100.0), "LQG-PI, q_p = 1e5": sint(1e5, fi["q_i"])}


def t_wy(a, c, f) -> dict:
    an = a.anello_v
    C, Bu, Bw = an.uscite["p"], an.B[:, 0], an.B[:, 1]
    out = {}
    for n, A in g.varianti_cuscinetto(a).items():
        if c is None:
            if n == "sliding":
                T = g.frf(a.cascata.A, a.cascata.B[:, 1], a.cascata.uscite["p"], f)
            else:
                vx = np.r_[a.m.v, np.zeros(a.cascata.N - a.cascata.n)]
                T = g.frf(a.cascata.A + a.att.sigma0 * np.outer(a.cascata.B[:, 1], vx), a.cascata.B[:, 1], a.cascata.uscite["p"], f)
        else:
            T = g.frf(A, Bw, C, f) / (1 + g.anello_lqg(c, g.frf(A, Bu, C, f), f)[0])
        out[n] = T
    return out


def fig_nominale(t, ref, D, nomi, seg):
    fig, ax = plt.subplots(3, 2, figsize=(13, 8.5), sharex=True)
    d = 50
    for j, (q, tit) in enumerate((("azimut", "Azimuth"), ("polare", "Elevation"))):
        for n, et in enumerate(nomi, 1):
            r = D[et][q]
            e = (ref[q] - ref[q][0] - r["los"]) * MDEG
            ax[0, j].plot(t[::d], e[::d], lw=1.1, label=et, **pt._stile(n))
            ax[1, j].plot(t[::d], np.abs(e[::d]), lw=1.1, **pt._stile(n))
            ax[2, j].plot(t[::d], r["i"][::d], lw=0.9, **pt._stile(n))
        ax[0, j].set_title(tit, fontsize=9); ax[1, j].set_yscale("log"); ax[1, j].set_ylim(1e-3, None); ax[2, j].set_xlabel("time [s]")
        for q2 in ax[:, j]:
            q2.axvline(T_AVVIO, color=pt.INCHIOSTRO_SECONDARIO, lw=0.6, ls=":")
    ax[0, 0].set_ylabel("line-of-sight error [mdeg]"); ax[1, 0].set_ylabel("|error| [mdeg]"); ax[2, 0].set_ylabel("current [A]"); ax[0, 0].legend(fontsize=7)
    pt._titolo(fig, "Circle on the Dahl plant, no disturbances. Dotted: end of the transient left out of the metrics", seg)
    return fig


def fig_scenari(M, nomi, seg):
    fig, ax = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
    sc = [s for s, _ in SCENARI if s != "nominal"]
    w = 0.8 / len(nomi)
    x = np.arange(len(sc))
    for n, et in enumerate(nomi):
        st = pt._stile(n + 1)
        for q, k in zip(ax, ("dist_total_rms_urad", "E_rms_mdeg", "i_max_az_A")):
            q.bar(x + (n - (len(nomi) - 1) / 2) * w, [max(M[s][et][k], 1e-3) for s in sc], width=w * 0.9, color=st["color"], label=et)
    ax[0].set_yscale("log"); ax[1].set_yscale("log")
    ax[0].set_ylabel("effect of the disturbance on the line of sight, rms [urad]"); ax[1].set_ylabel("total pointing error, rms [mdeg]"); ax[2].set_ylabel("largest azimuth current [A]")
    ax[2].axhline(18.0, color=pt.INCHIOSTRO, lw=0.7, ls="--")
    corto = {"torque step 0.5 N m": "step 0.5 N m", "Gaussian 0.1 N m rms": "Gaussian 0.1 N m", "Gaussian 0.1 N m rms, mean 0.2 N m": "Gaussian, mean 0.2 N m", TUTTI_DIST: "all, with mean"}
    for q in ax[:2]:
        q.tick_params(labelbottom=False)
    ax[2].set_xticks(x); ax[2].set_xticklabels([corto.get(s, s) for s in sc], rotation=15, ha="right", fontsize=8); ax[0].legend(fontsize=7)
    pt._titolo(fig, "Circle on the Dahl plant. Top: difference from the same controller with no disturbance. Dashed: 18 A", seg)
    return fig


def fig_anello(E, seg):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for j, (ch, tit) in enumerate((("azimuth", "Azimuth"), ("elevation", "Elevation"))):
        d = E[ch]["loop_decomposition"]
        for n, k in enumerate(("cascade", "LQG", "K_f = 0", "fast estimator"), 1):
            ax[j].plot(d["f_Hz"], d[k], marker="o", ms=4, lw=1.2, label=k, **pt._stile(n))
        ax[j].plot(d["f_Hz"], d["fast estimator, K_f = 0"], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO, label="fast estimator, K_f = 0 (k_p alone)")
        ax[j].axhline(0.0, color=pt.INCHIOSTRO_SECONDARIO, lw=0.7, ls=":"); ax[j].set_xscale("log"); ax[j].set_xlabel("frequency [Hz]"); ax[j].set_title(tit, fontsize=9); ax[j].legend(fontsize=7)
    ax[0].set_ylabel("|L| on the plant with bearings sliding [dB]")
    pt._titolo(fig, "Where the feedback gain goes: loop gain of the tuned LQG-PI, decomposed", seg)
    return fig


def fig_gradino(t, ref, D, nomi, seg):
    fig, ax = plt.subplots(2, 2, figsize=(13, 7), sharex=True)
    k = (t >= T_GRADINO - 0.5) & (t <= T_GRADINO + 8.0)
    for j, (q, tit) in enumerate((("azimut", "Azimuth"), ("polare", "Elevation"))):
        for n, et in enumerate(nomi, 1):
            r = D[et][q]
            ax[0, j].plot(t[k], ((ref[q] - ref[q][0] - r["los"]) * MDEG)[k], lw=1.1, label=et, **pt._stile(n))
            ax[1, j].plot(t[k], r["i"][k], lw=0.9, **pt._stile(n))
        ax[0, j].set_title(tit + ": torque step of 0.5 N m at 15 s", fontsize=9); ax[1, j].set_xlabel("time [s]")
    ax[0, 0].set_ylabel("line-of-sight error [mdeg]"); ax[1, 0].set_ylabel("current [A]"); ax[0, 0].legend(fontsize=7)
    pt._titolo(fig, "Circle on the Dahl plant, torque step", seg)
    return fig


def fig_twy(A, C, seg):
    f = np.geomspace(0.01, 150.0, 800)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for j, (a, cc, tit) in enumerate(zip(A, C, ("Azimuth", "Elevation"))):
        for n, (et, c) in enumerate(cc.items(), 1):
            T = t_wy(a, c, f)
            ax[j].plot(f, pt._db(T["sliding"]), lw=1.2, label=et + ", sliding", **pt._stile(n))
        ax[j].set_xscale("log"); ax[j].set_xlabel("frequency [Hz]"); ax[j].set_title(tit, fontsize=9); ax[j].legend(fontsize=7)
    ax[0].set_ylabel("T_wy, bearing torque to encoder [dB rad/(N m)]")
    pt._titolo(fig, "Disturbance rejection on the plant with bearings sliding", seg)
    return fig


def effetto_disturbo(t, ref, R, R0) -> dict:
    """Differenza dal caso senza disturbi con lo stesso controllore: linea di vista ed encoder, in urad."""
    m = t >= T_AVVIO
    dx = (R["azimut"]["los"] - R0["azimut"]["los"]) * np.cos(ref["polare"])
    de = R["polare"]["los"] - R0["polare"]["los"]
    rms = lambda x: 1e6 * float(np.sqrt(np.mean(x[m] ** 2)))
    return {"dist_total_rms_urad": 1e6 * float(np.sqrt(np.mean(dx[m] ** 2 + de[m] ** 2))), "dist_total_peak_urad": 1e6 * float(np.sqrt(dx[m] ** 2 + de[m] ** 2).max()), "dist_los_az_rms_urad": rms(dx), "dist_los_el_rms_urad": rms(de),
            "dist_enc_az_rms_urad": rms(R["azimut"]["p"] - R0["azimut"]["p"]), "dist_enc_el_rms_urad": rms(R["polare"]["p"] - R0["polare"]["p"])}


def identifica_seni(a, casi_, riga, r, t, asse) -> dict:
    """Ampiezza alla frequenza del disturbo, sulla differenza dal caso senza disturbi, contro la T_wy lineare con cuscinetti in strisciamento e bloccati."""
    A0 = g.varianti_cuscinetto(a)
    an = a.anello_v
    out = {}
    k = t >= 8.0
    for et, c in casi_.items():
        r0 = corri(a, c, r, t, {}, riga, SEME)
        righe = []
        for f in FREQ_SENI:
            x = corri(a, c, r, t, {"sin": (AMP_SENI, float(f))}, riga, SEME)
            X = np.c_[np.sin(2 * np.pi * f * t[k]), np.cos(2 * np.pi * f * t[k])]
            amp = lambda n: 1e6 * float(np.linalg.norm(np.linalg.lstsq(X, (x[n] - r0[n])[k], rcond=None)[0]))
            ln = {}
            ff = np.array([float(f)])
            for n, A in A0.items():
                if c is None:
                    Ac = a.cascata.A if n == "sliding" else a.cascata.A + a.att.sigma0 * np.outer(a.cascata.B[:, 1], np.r_[a.m.v, np.zeros(a.cascata.N - a.cascata.n)])
                    ln[n] = 1e6 * AMP_SENI * float(abs(g.frf(Ac, a.cascata.B[:, 1], a.cascata.uscite["p"], ff))[0])
                else:
                    Tw = g.frf(A, an.B[:, 1], an.uscite["p"], ff) / (1 + g.anello_lqg(c, g.frf(A, an.B[:, 0], an.uscite["p"], ff), ff)[0])
                    ln[n] = 1e6 * AMP_SENI * float(abs(Tw)[0])
            righe.append({"f_Hz": float(f), "encoder_urad": amp("p"), "line_of_sight_urad": amp("los"), "linear_sliding_urad": ln["sliding"], "linear_stuck_urad": ln["stuck"]})
        out[et] = righe
        print(asse, et, [(q["f_Hz"], round(q["encoder_urad"], 1), round(q["line_of_sight_urad"], 1)) for q in righe], flush=True)
    return out


def scomponi_anello(a, c, f) -> dict:
    """Guadagno d'anello: LQG; K_f = 0; stimatore veloce (rho = 1e6, recupero dell'anello del regolatore); entrambi, cioe' il solo k_p."""
    an = a.anello_v
    P = g.frf(g.varianti_cuscinetto(a)["sliding"], an.B[:, 0], an.uscite["p"], f)
    cf = lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, c["q_i"], c["q_p"], c["q"], rho_e=1e6, v_continua=True)
    c0 = lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, c["q_i"], c["q_p"], c["q"], v_continua=True)
    db = lambda x: pt._db(g.anello_lqg(x, P, f)[0])
    return {"f_Hz": f.tolist(), "cascade": pt._db(a.drive.p.Kpp * P).tolist(), "LQG": db(c0).tolist(), "K_f = 0": db(dataclasses.replace(c0, K_f=0 * c0.K_f)).tolist(), "fast estimator": db(cf).tolist(), "fast estimator, K_f = 0": db(dataclasses.replace(cf, K_f=0 * cf.K_f)).tolist(),
            "T_wy_no_position_loop_dB": pt._db(g.frf(g.varianti_cuscinetto(a)["sliding"], an.B[:, 1], an.uscite["p"], f)).tolist(), "k_p": c0.k_p, "K_e_position_per_s": float(c0.K_e[0, 0] / g.PASSO_MASTER)}


def accoppiamento(par) -> dict:
    """Modello accoppiato del progetto (senza cuscinetti e masse di statore) nel punto della circonferenza: dalle coppie agli encoder, |G| in dB."""
    acc = afm.sistema(par).accoppiato
    th, Om = float(np.arctan2(apl.QUOTA, apl.CENTRO)), float(apl.RAGGIO * apl.OMEGA / apl.CENTRO)
    ml = lz.linearizza(acc, th, Om)
    righe = []
    for f in (0.3, 1.0, 10.0, 56.0, 100.0):
        s_ = 2j * np.pi * f
        Y = ml.Cq @ np.linalg.solve(s_ * s_ * ml.M0 + s_ * ml.C0 + ml.K0, ml.E)
        db = 20 * np.log10(np.abs(Y[:2, :2]))
        righe.append({"f_Hz": f, "az_torque_to_az": float(db[0, 0]), "el_torque_to_el": float(db[1, 1]), "el_torque_to_az": float(db[0, 1]), "az_torque_to_el": float(db[1, 0])})
    return {"theta_rad": th, "Omega_rad_s": Om, "rows": righe}


def extra(assi, Kt, S, T, mot, Fi) -> dict:
    """Banda di comando e reiezione sul plant in strisciamento, guadagno d'anello e doppietti, per ogni caso."""
    f = np.geomspace(0.01, 150.0, 1500)
    out = {}
    t = np.arange(0.0, 15.0, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        m0 = next(m for m in a.modi if len(m.indici) == 2 and m.zeta < 0.1)
        riga = np.zeros(a.cascata.N)
        riga[a.m.nodi.index(NODO_LOS[asse])] = 1.0
        fi = Fi[ch]["final"]
        q = np.full(len(a.Ap) - 1, 0.1)
        idx = [i - 1 for i in m0.indici]
        if fi["q_f"]:
            q[idx] = float(fi["q_f"]["0"])
        id_c = {"cascade": None, "LQG-PI tuned": casi(a, fi, m0.f_hz)["LQG-PI tuned"], "LQG-PI, V = 100 Q": casi(a, fi, m0.f_hz)["LQG-PI, V = 100 Q"]}
        for x in SCALA_QF_AZ if asse == "azimut" else ():
            q2 = q.copy()
            q2[idx] = x
            id_c[f"LQG-PI, q_f = {x:g}"] = lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, fi["q_i"], fi["q_p"], q2, v_continua=True)
        seni = identifica_seni(a, id_c, riga, ref[asse] - ref[asse][0], t, asse)
        stabili = {et: {n: float(np.abs(cm.p1.poli_chiusi(c, A, a.anello_v.B[:, 0], a.anello_v.uscite["p"])).max()) for n, A in g.varianti_cuscinetto(a).items()} for et, c in id_c.items() if c is not None}
        scom = scomponi_anello(a, {"q_p": fi["q_p"], "q_i": fi["q_i"], "q": q}, F_ANELLO)
        W = cm.p24.comandi_tipici(asse)
        r = {}
        an = a.anello_v
        P = g.frf(g.varianti_cuscinetto(a)["sliding"], an.B[:, 0], an.uscite["p"], f)
        for et, c in casi(a, Fi[ch]["final"], m0.f_hz).items():
            Tw = t_wy(a, c, f)["sliding"]
            db = pt._db(Tw)
            q = {"T_wy_peak_dB": float(db.max()), "T_wy_peak_Hz": float(f[db.argmax()]), "T_wy_at_1Hz_dB": float(np.interp(1.0, f, db)), "T_wy_at_10Hz_dB": float(np.interp(10.0, f, db)), "T_wy_at_0.1Hz_dB": float(np.interp(0.1, f, db)),
                 "doublets": cm.prova(a, c, W)}
            if c is not None:
                L, N = g.anello_lqg(c, P, f)
                q |= {"L_at_1Hz_dB": float(np.interp(1.0, f, pt._db(L))), "bandwidth_T_ry_sliding_Hz": g.banda(f, N / (1 + L))}
            else:
                q |= {"L_at_1Hz_dB": float(np.interp(1.0, f, pt._db(a.drive.p.Kpp * P))), "bandwidth_T_ry_sliding_Hz": g.banda(f, g.frf(a.cascata.A, a.cascata.B[:, 0], a.cascata.uscite["p"], f))}
            r[et] = q
            print(ch, et, {k: round(v, 2) for k, v in q.items() if k != "doublets"}, round(q["doublets"]["err_rms_median_mdeg"]), round(q["doublets"]["i_max_A"], 1), round(q["doublets"]["over_max_mdeg"]), round(q["doublets"]["t_set_median_s"], 2), flush=True)
        out[ch] = {"cases": r, "sines": seni, "sines_largest_pole": stabili, "loop_decomposition": scom}
    out["coupling"] = accoppiamento(cm.pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml"))
    (CARTELLA / "extra.json").write_text(json.dumps(cm.p24.senza_array(out), indent=1))
    return out


def dither(assi, Kt, S, T, mot, Fi) -> dict:
    """Effetto del seno a 56 Hz sulla linea di vista: sotto 10 Hz, sopra 30 Hz, e lontano dalle inversioni del riferimento."""
    t = np.arange(0.0, T_FINE, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    m = t >= T_AVVIO
    lp, hp = (butter(4, f, b, fs=1 / g.PASSO_SIM, output="sos") for f, b in ((10.0, "low"), (30.0, "high")))
    rms = lambda x, k: float(np.sqrt(np.mean(x[k] ** 2)))
    out = {}
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        m0 = next(x for x in a.modi if len(x.indici) == 2 and x.zeta < 0.1)
        C = casi(a, Fi[ch]["final"], m0.f_hz)
        riga = np.zeros(a.cascata.N)
        riga[a.m.nodi.index(NODO_LOS[asse])] = 1.0
        r = ref[asse] - ref[asse][0]
        v = np.gradient(r, g.PASSO_SIM)
        inversioni = t[1:][np.diff(np.sign(v)) != 0]
        lontano = m.copy()
        for q in inversioni:
            lontano &= np.abs(t - q) > 3.0
        out[ch] = {"reversals_s": inversioni.tolist()}
        for et in ("cascade", "LQG-PI tuned"):
            r0, x = corri(a, C[et], r, t, {}, riga, SEME), corri(a, C[et], r, t, {"sin": (AMP_SENI, 56.0)}, riga, SEME)
            d = x["los"] - r0["los"]
            dl, dh = sosfiltfilt(lp, d), sosfiltfilt(hp, d)
            out[ch][et] = {"low_urad": 1e6 * rms(dl, m), "high_urad": 1e6 * rms(dh, m), "low_far_urad": 1e6 * rms(dl, lontano), "error_rms_nominal_mdeg": MDEG * rms(r - r0["los"], m), "error_rms_with_sine_mdeg": MDEG * rms(r - x["los"], m)}
    (CARTELLA / "dither.json").write_text(json.dumps(cm.p24.senza_array(out), indent=1))
    return out


def main() -> None:
    par, Kt, S, T, mot, assi = g.carica()
    if "--dither" in sys.argv:
        dither(assi, Kt, S, T, mot, json.loads(TUNING.read_text()))
        return
    if "--extra" in sys.argv:
        extra(assi, Kt, S, T, mot, json.loads(TUNING.read_text()))
        return
    seg = par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    Fi = json.loads(TUNING.read_text())
    t = np.arange(0.0, T_FINE, g.PASSO_SIM)
    ref = apl.circonferenza(t)
    scala_J = float(S["azimuth"]["uniform_lnf"]["hybrid"]["parameters"]["s_J_L"])
    cl = CoppieLagrange(afm.sistema(par).payload, scala_J)
    verifica = verifica_coppie(par) | {"azimuth_load_factor": scala_J}
    print("coppie contro Lagrange:", verifica["largest_difference_N_m"], flush=True)
    Ax, Cc, righe = {}, {}, {}
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        m0 = next(m for m in a.modi if len(m.indici) == 2 and m.zeta < 0.1)
        Ax[asse], Cc[asse] = a, casi(a, Fi[ch]["final"], m0.f_hz)
        riga = np.zeros(a.cascata.N)
        riga[a.m.nodi.index(NODO_LOS[asse])] = 1.0
        righe[asse] = riga
    tutti = list(Cc["polare"])
    M, serie, nominali = {}, {}, {}
    for sc_n, sc in SCENARI:
        M[sc_n] = {}
        for et in tutti:
            R = simula_coppia(Ax, {q: Cc[q][et] for q in Ax}, ref, t, sc, righe, SEME, cl, float(ref["polare"][0]))
            nominali.setdefault(et, R)
            M[sc_n][et] = metriche(t, ref, R) | effetto_disturbo(t, ref, R, nominali[et])
            if sc_n in ("nominal", "torque step 0.5 N m", TUTTI_DIST):
                serie.setdefault(sc_n, {})[et] = R
            x = M[sc_n][et]
            print(sc_n, "|", et, "| E %.2f xel %.2f el %.2f dist rms %.1f peak %.1f az %.1f el %.1f jit10 %.2f iaz %.1f iel %.1f" % (x["E_rms_mdeg"], x["xel_rms_mdeg"], x["el_rms_mdeg"], x["dist_total_rms_urad"], x["dist_total_peak_urad"], x["dist_los_az_rms_urad"], x["dist_los_el_rms_urad"], x["jitter_10Hz_urad"], x["i_max_az_A"], x["i_max_el_A"]), flush=True)
    st = {q: {et: {n: p1 for n, p1 in ((n, float(np.abs(cm.p1.poli_chiusi(c, A, Ax[q].anello_v.B[:, 0], Ax[q].anello_v.uscite["p"])).max())) for n, A in g.varianti_cuscinetto(Ax[q]).items())} for et, c in Cc[q].items() if c is not None} for q in Ax}
    # a coppie spente la co-simulazione deve dare gli stessi risultati dei due assi separati
    controllo = {}
    for et in ("cascade", "LQG-PI tuned"):
        R0 = simula_coppia(Ax, {q: Cc[q][et] for q in Ax}, ref, t, {}, righe, SEME, None, float(ref["polare"][0]))
        controllo[et] = {q: 1e6 * float(np.abs(R0[q]["los"] - corri(Ax[q], Cc[q][et], ref[q] - ref[q][0], t, {}, righe[q], SEME)["los"]).max()) for q in Ax}
    print("differenza da assi separati [urad]:", controllo, flush=True)
    gu = {q: {et: {"k_p": c.k_p, "k_i": c.k_i} for et, c in Cc[q].items() if c is not None} for q in Ax}
    (CARTELLA / "summary.json").write_text(json.dumps(cm.p24.senza_array({"scenarios": M, "largest_pole": st, "gains": gu, "K_pp": {q: Ax[q].drive.p.Kpp for q in Ax}, "coupling_check_against_lagrange": verifica, "cosimulation_vs_separate_urad": controllo, "coupling_effect": prova_coppie(Ax, Cc, righe, cl), "stability_vs_elevation": stabilita_theta(Ax, Cc, cl)}), indent=1))
    salva = lambda f, n: (pt.salva(f, CARTELLA / "figures", n), plt.close(f))
    principali = ["cascade", "LQG-PI tuned", "LQG-PI, k_i = 0.5 k_p"]
    salva(fig_nominale(t, ref, serie["nominal"], principali, seg), "circle_nominal")
    salva(fig_scenari(M, ["cascade", "LQG-PI tuned", "LQG-PI, V = 100 Q", "LQG-PI, q_p = 1e5"], seg), "circle_scenarios")
    salva(fig_gradino(t, ref, serie["torque step 0.5 N m"], principali, seg), "circle_torque_step")
    salva(fig_nominale(t, ref, serie[TUTTI_DIST], principali, seg), "circle_all_disturbances")
    salva(fig_twy([Ax["azimut"], Ax["polare"]], [{k: Cc[q][k] for k in principali} for q in ("azimut", "polare")], seg), "disturbance_rejection")
    if "--senza-extra" not in sys.argv:
        salva(fig_anello(extra(assi, Kt, S, T, mot, Fi), seg), "loop_decomposition")
        dither(assi, Kt, S, T, mot, Fi)


if __name__ == "__main__":
    main()
