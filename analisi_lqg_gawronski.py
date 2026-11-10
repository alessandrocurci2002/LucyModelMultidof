"""LQG di Gawronski (PI + flexible mode controller) sui doppietti V2, a confronto con la cascata posizione-velocita'-corrente.

Richiesta dell'utente (2026-11-01). Il plant della sintesi e' la maglia di velocita' chiusa del drive (come nel libro,
Fig. 8.1): corrente e velocita' restano nell'AKD, l'LQG sostituisce l'anello di posizione nel master a 1 kHz. La struttura
e' la stima ibrida su griglia uniforme in ln f, senza attrito nella sintesi; K_pv e K_iv sono quelli efficaci della
cascata, cosi' il confronto riguarda il solo anello di posizione.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_lqg_gawronski.py
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
sys.path[:0] = [str(RADICE / "src"), str(RADICE), str(RADICE / "closedLoopModel")]

import matplotlib.pyplot as plt  # noqa: E402
from scipy.linalg import expm  # noqa: E402

import analisi_anello_chiuso_ibrido as aac  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import doppietti_tempo as dt  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402
from lucy_cl.driver import Drive, DriveParams  # noqa: E402
from lucy_cl.motor import MotorParams  # noqa: E402

CARTELLA = RADICE / "results_lqg_gawronski"
RIC = RADICE / "results_fork_modes_dahl_resampled"
PASSO_MASTER = 1e-3
PASSO_SIM = 1e-4
RITARDO = 1
"""Periodi di ritardo fra il calcolo nel master e il drive (EtherCAT)."""
F_MAX_MODELLO = 150.0
HINF_MIN = 1e-5
ANOMALA = ("azimut", "8deg_minus8deg", "#1")
FIGURA = (("01deg_minus01deg", "#1"), ("1deg_minus1deg", "#1"), ("5deg_minus5deg", "#2"))
BANDA_MODULO = 5.0
"""Banda d'errore dell'assestamento, mdeg."""


class DriveVelocita:
    """Il drive in modo velocita': il primo ingresso di ``AnelloChiuso`` e' il comando di velocita'."""

    def __init__(self, d: Drive):
        self.d, self.p, self.n_states, self.i_zv = d, d.p, d.n_states, d.i_zv

    def step(self, x, v_cmd, _dp, _p_fb, v_fb, i_meas):
        return self.d.step(x, 0.0, v_cmd, 0.0, v_fb, i_meas)


@dataclasses.dataclass
class Asse:
    nome: str
    chiave: str
    m: object
    v: dict
    att: dt.Attrito
    cascata: aac.AnelloChiuso
    anello_v: aac.AnelloChiuso
    drive: Drive
    Ap: np.ndarray
    Bp: np.ndarray
    Cp: np.ndarray
    modi: list
    modi_tutti: list


def costruisci(asse, chiave, cat, mod, Kt, S, T, mot, g_scala=1.0, cuscinetti: dict | None = None) -> Asse:
    """``cuscinetti``: k_0, k_1, a della linearizzazione armonica nel modello di sintesi; ``None``: modello senza cuscinetti.
    La simulazione usa sempre l'attrito di Dahl."""
    v = dict(S[chiave]["uniform_lnf"]["hybrid"]["parameters"])
    v["g_Kt"] *= g_scala
    m = mod(cat, v, Kt, {}).struttura
    akd, wb = aac.DRIVE[asse]
    g = T["gains"][chiave]["uniform_lnf"]
    p = dataclasses.replace(DriveParams.from_files(chiave, aac.VCL / akd, aac.VCL / wb, pl_kp_unit="hz"), **g)
    att = dt.Attrito(v["sigma0"], v["Tc"])
    d_c = Drive(p)
    d_v = Drive(dataclasses.replace(p, Kpp=0.0, Kip_hz=0.0, kvff=1.0))
    cascata = aac.AnelloChiuso(m, d_c, mot, att, PASSO_SIM)
    anello = aac.AnelloChiuso(m, DriveVelocita(d_v), mot, att, PASSO_SIM)
    As, Bs, Cs = anello.A, anello.B[:, 0], anello.uscite["p"]
    if cuscinetti is not None:
        vx = np.r_[m.v, np.zeros(anello.N - anello.n)]
        vv = np.r_[np.zeros(anello.n), m.v, np.zeros(anello.N - 2 * anello.n)]
        As, Bs, Cs = lqg.con_cuscinetti(As, Bs, anello.B[:, 1], Cs, vx, vv, cuscinetti["k_0"], cuscinetti["k_1"], cuscinetti["a"])
    Am, Bm, Cm, modi = lqg.forma_modale(As, Bs[:, None], Cs[None, :])
    Ar, Br, Cr, scelti = lqg.tronca(Am, Bm, Cm, modi, lambda x: x.f_hz <= F_MAX_MODELLO and x.hinf >= HINF_MIN)
    Ap, Bp, Cp = lqg.posizione_primo_stato(Ar, Br, Cr)
    return Asse(asse, chiave, m, v, att, cascata, anello, d_c, Ap, Bp, Cp, scelti, modi)


def corri_lqg(a: Asse, c: lqg.ControlloreLQG, pc: np.ndarray, I0: float, ritardo: int = RITARDO, coppia: np.ndarray | None = None, riga: np.ndarray | None = None) -> dict:
    """Stessa integrazione di ``AnelloChiuso.corri``, con il comando di velocita' a passi di ``PASSO_MASTER``. ``coppia``: disturbo sulla porta
    dei cuscinetti, passo per passo; ``riga``: uscita in piu', relativa all'inizio (la linea di vista)."""
    an, att = a.anello_v, a.att
    s0, n, sv = att.sigma0, an.n, a.m.v
    Tb0 = float(np.clip(an.m.Kt * I0, -att.Tc, att.Tc))
    z = Tb0 / s0
    s = a.cascata.equilibrio(Tb0)
    p0 = an.uscite["p"] @ s
    vv = np.r_[np.zeros(n), sv, np.zeros(an.N - 2 * n)]
    sotto = int(round(PASSO_MASTER / PASSO_SIM))
    st = c.stato_iniziale()
    coda = [0.0] * ritardo
    P, I, IC, VP, VI, U = (np.empty(len(pc)) for _ in range(6))
    Ad, Bd, rp, ri, ric, rin, rvp, vpin, rvi = (an.Ad, an.Bd, an.uscite["p"], an.uscite["i"], an.ic_row, an.ic_in, an.vp_row,
                                                an.vp_in, an.vint_row)
    u_app = y_s = 0.0
    L = np.empty(len(pc))
    l0 = riga @ s if riga is not None else 0.0
    for k, r in enumerate(pc):
        y = rp @ s - p0
        if k % sotto == 0:
            y_s = y
            if hasattr(c, "misura"):  # aggancio: un controllore con ``misura`` (PI sull'errore misurato, Gawronski Fig. 8.4b) riceve y(k) prima di ``uscita``
                c.misura(st, y)
            coda.append(c.uscita(st, float(r)))
            u_app = coda.pop(0)
        ing = np.array([u_app, s0 * z + (coppia[k] if coppia is not None else 0.0)])
        P[k], I[k], U[k] = y, ri @ s, u_app
        L[k] = riga @ s - l0 if riga is not None else y
        IC[k], VP[k], VI[k] = ric @ s + rin @ ing, rvp @ s + vpin @ ing, rvi @ s
        s_new = Ad @ s + Bd @ ing
        vd = 0.5 * (vv @ s + vv @ s_new)
        al = s0 * abs(vd) / att.g(vd)
        if al > 0:
            e = np.exp(-al * PASSO_SIM)
            z = z * e + vd / al * (1 - e)
        s = s_new
        if k % sotto == sotto - 1:
            c.stima(st, y_s, u_app)
    return {"p": P, "i": I, "i_cmd": IC, "v_phase": VP, "i_int": VI, "u": U, "los": L}


def maschere(w):
    vc = np.abs(np.gradient(w.cmd, dt.DT)) > 1e-6
    i0, i1 = np.flatnonzero(vc)[[0, -1]]
    return i0, i1


def metriche(w, y: np.ndarray, im: np.ndarray) -> dict[str, float]:
    i0, i1 = maschere(w)
    mv = slice(i0, min(i1 + int(1.0 / dt.DT), len(y)))
    e = np.rad2deg((w.cmd - w.cmd[0]) - y) * 1e3
    fuori = np.flatnonzero(np.abs(e[i1:]) > BANDA_MODULO)
    finale = np.rad2deg(w.cmd[i1] - w.cmd[0]) * 1e3
    pos = np.rad2deg(y[i1:]) * 1e3
    corsa = np.rad2deg(w.cmd[i1] - w.cmd[i0]) * 1e3
    over = float(max(0.0, np.max((pos - finale) * np.sign(corsa))))
    return {"stroke_mdeg": abs(corsa), "err_max_mdeg": float(np.abs(e[mv]).max()), "err_rms_mdeg": float(np.sqrt(np.mean(e[mv] ** 2))),
            "t_set_s": float((fuori[-1] + 1) * dt.DT) if fuori.size else 0.0, "over_mdeg": over,
            "i_max_A": float(np.abs(im[mv]).max()), "i_rms_A": float(np.sqrt(np.mean(im[mv] ** 2)))}


def carica():
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    Kt = par.sperimentale.Kt
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    S = json.loads((RIC / "summary.json").read_text())
    T = json.loads((RIC / "time_domain_summary.json").read_text())
    mot = MotorParams.from_yaml(RADICE / "closedLoopModel" / "data" / "motor_TBM7646A.yaml", "cold")
    return par, Kt, S, T, mot, (("azimut", "azimuth", cat_az, fl.azimut), ("polare", "elevation", cat_pol, fl.polare))


def simula(a: Asse, c, w, t_max: float | None = None):
    """Riferimento misurato sul modello con LQG (``c``) o con la cascata (``c = None``); uscite sulla griglia di ``w``."""
    k = w.t <= (t_max if t_max else w.t[-1])
    t = np.arange(0.0, w.t[k][-1], PASSO_SIM)
    pc, I0 = np.interp(t, w.t, w.cmd - w.cmd[0]), float(w.I[:int(0.4 / dt.DT)].mean())
    r = a.cascata.corri(pc, I0) if c is None else corri_lqg(a, c, pc, I0)
    return {n: np.interp(w.t[k], t, r[n]) for n in ("p", "i")} | {f"max_{n}": float(np.abs(r[n]).max()) for n in ("i_cmd", "v_phase", "i_int")}


# ------------------------------------------------------------------------- analisi lineare


def varianti_cuscinetto(a: Asse) -> dict[str, np.ndarray]:
    """Matrice A della maglia di velocita' con i cuscinetti liberi (strisciamento) e bloccati (pre-scorrimento, sigma0)."""
    n = a.anello_v.n
    vv = np.r_[a.m.v, np.zeros(a.anello_v.N - n)]
    return {"sliding": a.anello_v.A, "stuck": a.anello_v.A + a.att.sigma0 * np.outer(a.anello_v.B[:, 1], vv)}


def frf(A, B, C, f):
    return np.array([C @ np.linalg.solve(2j * np.pi * x * np.eye(len(A)) - A, B) for x in f])


def zoh(f, T=PASSO_MASTER):
    return np.sinc(f * T) * np.exp(-1j * np.pi * f * T)


def anello_lqg(c, P, f, ritardo=RITARDO):
    z = np.exp(2j * np.pi * f * PASSO_MASTER)
    Fr, Fy = c.funzione_trasferimento(z, ritardo)
    return Fy * zoh(f) * P, Fr * zoh(f) * P


def margini(L: np.ndarray, f: np.ndarray) -> dict[str, float]:
    """Modulo minimo di 1 + L, margine di guadagno e di fase ai passaggi per -180 gradi e per 0 dB."""
    fase = np.unwrap(np.angle(L))
    mod = np.abs(L)
    pm, gm, fc = np.inf, np.inf, np.nan
    for i in np.flatnonzero((mod[:-1] - 1) * (mod[1:] - 1) < 0):
        x = abs(np.rad2deg(np.angle(-L[i])))
        if x < pm:
            pm, fc = x, f[i]
    k = np.round((fase + np.pi) / (2 * np.pi))
    for i in np.flatnonzero((fase[:-1] + np.pi - 2 * np.pi * k[:-1]) * (fase[1:] + np.pi - 2 * np.pi * k[:-1]) < 0):
        if mod[i] < 1:
            gm = min(gm, -20 * np.log10(mod[i]))
    return {"modulus": float(np.abs(1 + L).min()), "f_modulus_Hz": float(f[np.argmin(np.abs(1 + L))]), "pm_deg": float(pm),
            "gm_dB": float(gm), "f_c_Hz": float(fc)}


def autovalori_chiusi(a: Asse, c, A) -> np.ndarray:
    """Poli discreti dell'anello chiuso con il plant completo: stati [x, x_hat, e_i, u in coda]."""
    N, n = len(A), len(c.Ad)
    B, C = a.anello_v.B[:, 0], a.anello_v.uscite["p"]
    E = expm(np.block([[A, B[:, None]], [np.zeros((1, N + 1))]]) * PASSO_MASTER)
    Phi, Gam = E[:N, :N], E[:N, N]
    tot = N + n + 1 + RITARDO
    M = np.zeros((tot, tot))
    ix, ih, ie, iq = slice(0, N), slice(N, N + n), N + n, N + n + 1
    Ku = np.zeros(tot)
    Ku[ie] = c.k_i
    Ku[N] = -c.k_p
    Ku[N + 1:N + n] = -c.K_f
    pi_mis = isinstance(c, lqg.ControlloreLQGPIMisurato)  # Fig. 8.4b: il PI legge y = C x invece di y_hat
    if pi_mis:
        Ku[N] = 0.0
        Ku[ix] = -c.k_p * C
    M[ix, ix] = Phi
    M[ix, iq] = Gam
    M[ih, ix] = c.K_e[:, [0]] @ C[None, :]
    M[ih, ih] = c.Ad - c.K_e[:, [0]] @ np.eye(n)[:1]
    M[ih, iq] = c.Bd[:, 0]
    M[ie, ie] = 1.0
    if pi_mis:
        M[ie, ix] = -PASSO_MASTER * C
    else:
        M[ie, N] = -PASSO_MASTER
    M[iq, :] = Ku
    return np.linalg.eigvals(M)


def plant_misurato(a: Asse, f: np.ndarray, G: np.ndarray) -> np.ndarray:
    """y/v_cmd con il plant misurato (velocita'/corrente) e il drive del modello; la G di WorkBench comprende il loop di corrente."""
    d = a.drive
    s = 2j * np.pi * f
    Cv = d.C_vel(s) * d.F_fwd(s)
    return Cv * G / (1 + Cv * G * d.F_fb(s)) / s


F_GRID = np.geomspace(0.005, 480.0, 800)


def valuta_lineare(a: Asse, c, pm_cache: dict) -> dict:
    out = {}
    for nome, P in pm_cache.items():
        L, _ = anello_lqg(c, P[1], P[0])
        out[nome] = margini(L, P[0])
    return out


# ------------------------------------------------------------------------- taratura


KP = {"azimut": (10, 12, 14, 16, 18), "polare": (20, 30, 40, 60, 80, 100, 130, 160)}
RAPPORTI_KI = {"azimut": (0.03, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0), "polare": (0.002, 0.005, 0.01, 0.02, 0.03, 0.1)}
"""k_i come frazione del guadagno integrale critico 0.25 k_p^2 (7.9)."""
Q_F = (0.0, 1e3, 3e3, 1e4, 3e4)
RHO = (1.0, 1e-3)
INDICE_TOLL = 1.05
"""Peggioramento dell'indice ammesso per aumentare lo smorzamento dei modi strutturali."""
I_MAX = 15.0
MOD_MIN = 0.5


def pesi_da_guadagni(kp: float, rapporto: float) -> tuple[float, float]:
    """q_p, q_i per k_p, k_i dati sul modello rigido: k_i = sqrt(q_i), k_p^2 = q_p + 2 k_i."""
    ki = rapporto * 0.25 * kp ** 2
    return max(kp ** 2 - 2 * ki, 0.0), ki ** 2


def indice(ms: list[dict], rif: list[dict]) -> float:
    """Indice (8.27) ridotto alle metriche dei doppietti: errore rms rispetto alla cascata, assestamento rispetto a 1 s,
    sovraelongazione rispetto al 10 % della corsa. Con k_i > 0 l'anello e' di tipo 2 e sul gradino vale int e dt = 0: la
    sovraelongazione non si puo' annullare, quindi non si normalizza con quella quasi nulla della cascata."""
    return float(np.mean([m["err_rms_mdeg"] / r["err_rms_mdeg"] + m["t_set_s"] + m["over_mdeg"] / (0.1 * m["stroke_mdeg"])
                          for m, r in zip(ms, rif)]) / 3)


def valuta(a, c, prove, rif, pm_cache):
    ms = []
    for w in prove:
        r = simula(a, c, w)
        ms.append(metriche(w, r["p"], r["i"]))
    lin = valuta_lineare(a, c, pm_cache)
    stabile = all(np.abs(autovalori_chiusi(a, c, A)).max() < 1 + 1e-9 for A in varianti_cuscinetto(a).values())
    return {"indice": indice(ms, rif), "metriche": ms, "margini": lin, "stabile": bool(stabile),
            "mod_min": min(x["modulus"] for x in lin.values()), "i_max": max(m["i_max_A"] for m in ms)}


def taratura(a: Asse, mis, finestre) -> tuple[dict, list]:
    """Passi 8.2: pesi flessibili piccoli, q_p e q_i con i vincoli, poi correzione dei pesi flessibili e dello stimatore."""
    prove = []
    for t, mossa in FIGURA:
        w = next(x for x in finestre[t] if x.nome[-2:] == mossa)
        prove.append(w)
    rif = [metriche(w, (r := simula(a, None, w))["p"], r["i"]) for w in prove]
    pm_cache = {n: (F_GRID, frf(A, a.anello_v.B[:, 0], a.anello_v.uscite["p"], F_GRID))
                for n, A in varianti_cuscinetto(a).items()}
    k = (mis.f >= 1.0) & (mis.f <= 150.0)
    fm = mis.f[k]
    pm_cache["measured"] = (fm, plant_misurato(a, fm, mis.G[k]))
    nf = len(a.Ap) - 1
    registro, migliore = [], None
    for kp in KP[a.nome]:
        for ra in RAPPORTI_KI[a.nome]:
            qp, qi = pesi_da_guadagni(kp, ra)
            c = lqg.sintesi(a.Ap, a.Bp, a.Cp, PASSO_MASTER, qi, qp, np.zeros(nf))
            x = valuta(a, c, prove, rif, pm_cache)
            ok = x["stabile"] and x["i_max"] <= I_MAX and x["mod_min"] >= MOD_MIN and c.k_i <= 0.25 * c.k_p ** 2
            registro.append({"kp_target": kp, "ratio": ra, "q_p": qp, "q_i": qi, "k_p": c.k_p, "k_i": c.k_i, "ok": ok, "indice": x["indice"],
                             "i_max": x["i_max"], "mod_min": x["mod_min"], "stabile": x["stabile"]})
            print(a.nome, {k_: (round(v, 3) if isinstance(v, float) else v) for k_, v in registro[-1].items()}, flush=True)
            if ok and (migliore is None or x["indice"] < migliore["indice"]):
                migliore = {"q_p": qp, "q_i": qi, "indice": x["indice"]}
    if migliore is None:
        raise RuntimeError(f"{a.nome}: nessun candidato soddisfa i vincoli")
    # correzione dei pesi flessibili (solo modi strutturali) e dello stimatore
    fis = [m for m in a.modi if len(m.indici) == 2 and np.min(np.abs(a.m.omega / (2 * np.pi) - m.f_hz)) < 3.0]
    idx = [i - 1 for m in fis for i in m.indici]
    registro_f, best, zeta_best = [], {**migliore, "q_f": 0.0, "rho": 1.0}, 0.0
    for qf in Q_F:
        for rho in RHO:
            q = np.zeros(nf); q[idx] = qf
            c = lqg.sintesi(a.Ap, a.Bp, a.Cp, PASSO_MASTER, migliore["q_i"], migliore["q_p"], q, rho)
            x = valuta(a, c, prove, rif, pm_cache)
            ok = x["stabile"] and x["i_max"] <= I_MAX and x["mod_min"] >= MOD_MIN
            z = zeta_flessibili(a, c)
            registro_f.append({"q_f": qf, "rho": rho, "ok": ok, "indice": x["indice"], "i_max": x["i_max"], "mod_min": x["mod_min"], "zeta": z,
                               "K_f_flex": float(np.linalg.norm(c.K_f[idx]))})
            print(a.nome, "flex", {k_: (round(v, 4) if isinstance(v, float) else v) for k_, v in registro_f[-1].items()}, flush=True)
            zm = float(np.mean(list(z.values()))) if z else 0.0
            if ok and x["indice"] <= migliore["indice"] * INDICE_TOLL and zm > zeta_best + 1e-6:
                zeta_best = zm
                best = {"q_p": migliore["q_p"], "q_i": migliore["q_i"], "indice": x["indice"], "q_f": qf, "rho": rho}
    return best, registro, registro_f, idx, rif


def zeta_flessibili(a, c=None) -> dict[str, float]:
    """Smorzamento in anello chiuso dei poli vicini ai modi strutturali in banda, con i cuscinetti liberi."""
    if c is None:
        lam = np.linalg.eigvals(a.cascata.A)
    else:
        lam = np.log(autovalori_chiusi(a, c, a.anello_v.A).astype(complex)) / PASSO_MASTER
    out = {}
    for f0 in a.m.omega / (2 * np.pi):
        if 5.0 < f0 < 150.0:
            k = lam[lam.imag > 0]
            d = k[np.argmin(np.abs(np.abs(k) / (2 * np.pi) - f0) + 50 * (k.real / np.abs(k) < -0.12))]
            out[f"{f0:.1f}"] = float(-d.real / abs(d))
    return out


def misura_asse(cfg, asse):
    fn, et = afm.MISURA_AZIMUT if asse == "azimut" else afm.MISURE_POLARI[0]
    return afm.misura(cfg, fn, asse, et)


def finestre_asse(asse):
    return {t: dt.finestre(aac.VCL / "V2" / (aac.PROVE[asse] % t), asse) for t in aac.AMPIEZZE}


def controllore(a: Asse, b: dict) -> lqg.ControlloreLQG:
    nf = len(a.Ap) - 1
    q = np.zeros(nf)
    fis = [i - 1 for m in a.modi if len(m.indici) == 2 and np.min(np.abs(a.m.omega / (2 * np.pi) - m.f_hz)) < 3.0 for i in m.indici]
    q[fis] = b["q_f"]
    return lqg.sintesi(a.Ap, a.Bp, a.Cp, PASSO_MASTER, b["q_i"], b["q_p"], q, b["rho"])


def funzioni_anello(a: Asse, c, mis) -> dict:
    """T_ry e T_wy (w: coppia sui cuscinetti) di cascata e LQG, margini su plant liberi, bloccati e misurati."""
    f = F_GRID
    C, Bu, Bw = a.anello_v.uscite["p"], a.anello_v.B[:, 0], a.anello_v.B[:, 1]
    P = {n: frf(A, Bu, C, f) for n, A in varianti_cuscinetto(a).items()}
    Pw = frf(a.anello_v.A, Bw, C, f)
    L, Tr = anello_lqg(c, P["sliding"], f)
    k = (mis.f >= 1.0) & (mis.f <= 150.0)
    fm, Pm = mis.f[k], plant_misurato(a, mis.f[k], mis.G[k])
    Kpp = a.drive.p.Kpp
    marg = {"lqg": {n: margini(anello_lqg(c, x, f)[0], f) for n, x in P.items()} | {"measured": margini(anello_lqg(c, Pm, fm)[0], fm)},
            "cascade": {n: margini(Kpp * x, f) for n, x in P.items()} | {"measured": margini(Kpp * Pm, fm)}}
    return {"f": f, "T_lqg": Tr / (1 + L), "Tw_lqg": Pw / (1 + L), "T_cas": frf(a.cascata.A, a.cascata.B[:, 0], C, f),
            "Tw_cas": frf(a.cascata.A, a.cascata.B[:, 1], C, f), "P": P, "margini": marg, "fm": fm, "L_mis": anello_lqg(c, Pm, fm)[0],
            "L_mod": anello_lqg(c, frf(a.anello_v.A, Bu, C, fm), fm)[0], "Lc_mis": Kpp * Pm, "Lc_mod": Kpp * frf(a.anello_v.A, Bu, C, fm)}


def banda(f, T) -> float:
    k = np.flatnonzero(np.abs(T) < 10 ** (-3 / 20) * abs(T[0]))
    return float(f[k[0]]) if k.size else float(f[-1])


def robustezza(asse, chiave, cat, mod, Kt, S, T, mot, c, a) -> dict:
    """Poli e margini con g_Kt +-10 % e senza compensazione del ritardo (controllore fisso, plant variato)."""
    out = {}
    for g in (0.9, 1.0, 1.1):
        b = costruisci(asse, chiave, cat, mod, Kt, S, T, mot, g)
        for n, A in varianti_cuscinetto(b).items():
            P = frf(A, b.anello_v.B[:, 0], b.anello_v.uscite["p"], F_GRID)
            out[f"g{g:.1f}_{n}"] = {"rho_max": float(np.abs(autovalori_chiusi(b, c, A)).max()),
                                    "modulus": margini(anello_lqg(c, P, F_GRID)[0], F_GRID)["modulus"]}
    return out


def fig_modello(a: Asse, nome: str, seg):
    f = np.geomspace(0.1, 480, 900)
    C, B = a.anello_v.uscite["p"], a.anello_v.B[:, 0]
    G, Gr = frf(a.anello_v.A, B, C, f), np.array([(a.Cp @ np.linalg.solve(2j * np.pi * x * np.eye(len(a.Ap)) - a.Ap, a.Bp))[0, 0] for x in f])
    fig, ax = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True)
    ax[0].plot(f, pt._db(G), lw=1.5, color=pt.INCHIOSTRO, label="full velocity loop")
    ax[0].plot(f, pt._db(Gr), lw=1.3, label="reduced synthesis model", **pt._stile(1))
    ax[1].plot(f, np.rad2deg(np.angle(G)), lw=1.5, color=pt.INCHIOSTRO)
    ax[1].plot(f, np.rad2deg(np.angle(Gr)), lw=1.3, **pt._stile(1))
    for ax_ in ax:
        ax_.set_xscale("log")
    ax[0].set_ylabel("magnitude [dB rad/(rad/s)]"); ax[1].set_ylabel("phase [deg]"); ax[1].set_xlabel("frequency [Hz]")
    ax[1].set_ylim(-180, 180); ax[1].set_yticks(np.arange(-180, 181, 90)); ax[0].legend(fontsize=8)
    pt._titolo(fig, f"{nome}: velocity loop from velocity command to encoder, full and reduced", seg)
    return fig


def indice_cascata(rif: list[dict]) -> float:
    return indice(rif, rif)


def fig_scan(tun, seg):
    fig, assi = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, ch, t in zip(assi, ("azimuth", "elevation"), ("Azimuth", "Elevation")):
        reg = tun[ch]["scan"]
        for i, ra in enumerate(sorted({r["ratio"] for r in reg})[::2], 1):
            x = [r for r in reg if r["ratio"] == ra]
            ax.plot([r["k_p"] for r in x], [r["indice"] for r in x], lw=1.2, label=f"k_i = {ra:g} x critical", **pt._stile(i))
            okx = [r for r in x if r["ok"]]
            ax.plot([r["k_p"] for r in okx], [r["indice"] for r in okx], "o", ms=4, color=pt._stile(i)["color"])
        ax.set_title(t, fontsize=9); ax.set_xlabel("k_p [1/s]"); ax.set_xscale("log"); ax.axhline(indice_cascata(tun[ch]["reference"]), color=pt.INCHIOSTRO, lw=0.9, ls=":", label="cascade")
        ax.legend(fontsize=7)
    assi[0].set_ylabel("index f")
    pt._titolo(fig, "Tuning scan: index against k_p; markers: all constraints met", seg)
    return fig


def fig_tf(ris, seg):
    fig, ax = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    for j, (t, x) in enumerate(ris.items()):
        ax[0, j].plot(x["f"], pt._db(x["T_cas"]), lw=1.4, label="cascade", **pt._stile(1))
        ax[0, j].plot(x["f"], pt._db(x["T_lqg"]), lw=1.4, label="LQG", **pt._stile(2))
        ax[1, j].plot(x["f"], pt._db(x["Tw_cas"]), lw=1.4, **pt._stile(1))
        ax[1, j].plot(x["f"], pt._db(x["Tw_lqg"]), lw=1.4, **pt._stile(2))
        ax[0, j].set_title(t, fontsize=9); ax[1, j].set_xlabel("frequency [Hz]"); ax[0, j].set_ylim(-60, 12)
        for q in ax[:, j]:
            q.set_xscale("log")
    ax[0, 0].set_ylabel("T_ry, command to encoder [dB]"); ax[1, 0].set_ylabel("T_wy, bearing torque to encoder [dB rad/(N m)]")
    ax[0, 0].legend(fontsize=8)
    pt._titolo(fig, "Closed-loop transfer functions, bearings sliding (no friction in the linear model)", seg)
    return fig


def fig_anello(ris, seg):
    fig, ax = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    for j, (t, x) in enumerate(ris.items()):
        for L, et, st, lw in ((x["L_mod"], "LQG, model", pt._stile(2), 1.3), (x["L_mis"], "LQG, measured plant", pt._stile(1), 1.0),
                              (x["Lc_mis"], "cascade, measured plant", {"color": pt.INCHIOSTRO, "linestyle": ":"}, 1.0)):
            ax[0, j].plot(x["fm"], pt._db(L), lw=lw, label=et, **st)
            ax[1, j].plot(x["fm"], np.abs(1 + L), lw=lw, **st)
        ax[0, j].axhline(0, color=pt.INCHIOSTRO, lw=0.6); ax[1, j].axhline(MOD_MIN, color=pt.INCHIOSTRO, lw=0.6, ls="--")
        ax[0, j].set_title(t, fontsize=9); ax[1, j].set_xlabel("frequency [Hz]"); ax[1, j].set_yscale("log")
        for q in ax[:, j]:
            q.set_xscale("log")
    ax[0, 0].set_ylabel("|L| [dB]"); ax[1, 0].set_ylabel("|1 + L|; dashed: required minimum"); ax[0, 0].legend(fontsize=8)
    pt._titolo(fig, "Position loop opened at the encoder, with the model and with the measured frequency response", seg)
    return fig


def fig_doppietti(nome, curve, seg):
    fig, ax = plt.subplots(3, 3, figsize=(15, 9.5), sharex="col", gridspec_kw={"height_ratios": [2, 1.3, 1.6]})
    for j, key in enumerate(FIGURA):
        w, yc, ic, yl, il = curve[key]
        i0 = maschere(w)[0]
        k = (w.t >= w.t[i0] - 0.2) & (w.t <= w.t[i0] + 2.5)
        t = w.t[k]
        md = lambda q: np.rad2deg(q[k]) * 1e3  # noqa: E731
        rif = w.cmd - w.cmd[0]
        ax[0, j].plot(t, md(rif), color=pt.INCHIOSTRO, lw=0.8, ls=":", label="position command")
        ax[0, j].plot(t, md(w.pos), color=pt.INCHIOSTRO, lw=1.6, label="measured encoder (cascade on the real axis)")
        ax[0, j].plot(t, md(yc), lw=1.3, label="model, cascade", **pt._stile(1))
        ax[0, j].plot(t, md(yl), lw=1.3, label="model, LQG", **pt._stile(2))
        ax[1, j].plot(t, md(rif - w.pos), color=pt.INCHIOSTRO, lw=1.3)
        ax[1, j].plot(t, md(rif - yc), lw=1.2, **pt._stile(1))
        ax[1, j].plot(t, md(rif - yl), lw=1.2, **pt._stile(2))
        ax[2, j].plot(t, w.I[k], color=pt.INCHIOSTRO, lw=1.0)
        ax[2, j].plot(t, ic[k], lw=1.1, **pt._stile(1))
        ax[2, j].plot(t, il[k], lw=1.1, **pt._stile(2))
        ax[0, j].set_title(f"{key[0].split('_')[0].replace('01deg', '0.1deg').replace('05deg', '0.5deg')} doublet, move {key[1]}", fontsize=9)
        ax[2, j].set_xlabel("time from window start [s]")
    ax[0, 0].set_ylabel("position [mdeg]"); ax[1, 0].set_ylabel("following error [mdeg]"); ax[2, 0].set_ylabel("motor current [A rms]")
    ax[0, 1].legend(loc="lower right", fontsize=7)
    pt._titolo(fig, f"{nome}: V2 doublets, cascade and LQG on the model, same reference and same velocity loop", seg)
    return fig


def main() -> None:
    par, Kt, S, T, mot, assi = carica()
    cfg, seg = par.sperimentale, par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    fileT = CARTELLA / "tuning.json"
    tun = json.loads(fileT.read_text()) if "--riusa-taratura" in sys.argv and fileT.exists() else {}
    righe, ris, finale = [], {}, {}
    for asse, chiave, cat, mod in assi:
        nome = "Azimuth" if asse == "azimut" else "Elevation"
        a = costruisci(asse, chiave, cat, mod, Kt, S, T, mot)
        fin, mis = finestre_asse(asse), misura_asse(cfg, asse)
        if "best" not in tun.get(chiave, {}):
            best, reg, regf, idx, rif = taratura(a, mis, fin)
            tun[chiave] = {"best": best, "scan": reg, "scan_flex": regf, "flex_index": idx, "reference": rif,
                           "modes": [dataclasses.asdict(m) for m in a.modi], "modes_all": [dataclasses.asdict(m) for m in a.modi_tutti],
                           "zeta_cascade": zeta_flessibili(a)}
            fileT.write_text(json.dumps(tun, indent=1, default=float))
        c = controllore(a, tun[chiave]["best"])
        print(asse, "guadagni", c.k_p, c.k_i, c.K_f, flush=True)
        x = funzioni_anello(a, c, mis)
        ris[nome] = x
        curve, picchi = {}, {"cascade": {"i_cmd": 0.0, "v_phase": 0.0, "i_int": 0.0}, "lqg": {"i_cmd": 0.0, "v_phase": 0.0, "i_int": 0.0}}
        for t, ww in fin.items():
            for w in ww:
                rc, rl = simula(a, None, w), simula(a, c, w)
                for n in picchi["lqg"]:
                    picchi["cascade"][n] = max(picchi["cascade"][n], rc[f"max_{n}"])
                    picchi["lqg"][n] = max(picchi["lqg"][n], rl[f"max_{n}"])
                ris_w = {"measured": metriche(w, w.pos, w.I), "cascade": metriche(w, rc["p"], rc["i"]), "lqg": metriche(w, rl["p"], rl["i"])}
                for sch, mm in ris_w.items():
                    righe.append({"axis": asse, "test": t, "move": w.nome[-2:], "scheme": sch} | mm)
                curve[(t, w.nome[-2:])] = (w, rc["p"], rc["i"], rl["p"], rl["i"])
        pt.salva(fig_doppietti(nome, curve, seg), CARTELLA / "figures", f"{chiave}_doublets")
        pt.salva(fig_modello(a, nome, seg), CARTELLA / "figures", f"{chiave}_velocity_loop")
        plt.close("all")
        rob = robustezza(asse, chiave, cat, mod, Kt, S, T, mot, c, a)
        finale[chiave] = {"k_p": c.k_p, "k_i": c.k_i, "K_f": c.K_f, "K_e": c.K_e[:, 0], "weights": tun[chiave]["best"],
                          "regulator_poles_hz_zeta": [[float(abs(np.log(z)) / (2 * np.pi * PASSO_MASTER)), float(-np.log(z).real / abs(np.log(z)))]
                                                      for z in c.poli_regolatore.astype(complex) if z.imag >= 0],
                          "estimator_poles_hz_zeta": [[float(abs(np.log(z)) / (2 * np.pi * PASSO_MASTER)), float(-np.log(z).real / abs(np.log(z)))]
                                                      for z in c.poli_stimatore.astype(complex) if z.imag >= 0],
                          "margins": x["margini"], "bandwidth_Hz": {"cascade": banda(x["f"], x["T_cas"]), "lqg": banda(x["f"], x["T_lqg"])},
                          "peak_T_dB": {"cascade": float(pt._db(x["T_cas"]).max()), "lqg": float(pt._db(x["T_lqg"]).max())},
                          "T_wy_at_1Hz_dB": {"cascade": float(np.interp(1.0, x["f"], pt._db(x["Tw_cas"]))), "lqg": float(np.interp(1.0, x["f"], pt._db(x["Tw_lqg"])))},
                          "zeta_flex": {"cascade": zeta_flessibili(a), "lqg": zeta_flessibili(a, c), "open_velocity_loop": zeta_flessibili_v(a)},
                          "limits": picchi | {"i_limit": a.drive.p.i_limit, "vel_int_limit": a.drive.p.vel_int_limit, "V_phase_max": a.drive.V_phase_max},
                          "robustness": rob, "Kpp_cascade": a.drive.p.Kpp, "Kpv": a.drive.p.Kpv, "Kiv_hz": a.drive.p.Kiv_hz,
                          "modes": [dataclasses.asdict(m) for m in a.modi], "states_full": a.anello_v.N, "states_reduced": len(a.Ap)}
        print(asse, json.dumps({k: finale[chiave][k] for k in ("bandwidth_Hz", "peak_T_dB", "T_wy_at_1Hz_dB", "zeta_flex", "limits", "margins")}, indent=0, default=float), flush=True)
        (CARTELLA / "summary.json").write_text(json.dumps(finale, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))
    pt.salva(fig_tf(ris, seg), CARTELLA / "figures", "closed_loop_tf")
    pt.salva(fig_anello(ris, seg), CARTELLA / "figures", "position_loop")
    pt.salva(fig_scan(tun, seg), CARTELLA / "figures", "tuning_scan")
    plt.close("all")
    with open(CARTELLA / "doublet_metrics.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(righe[0]))
        wr.writeheader()
        wr.writerows([{k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()} for r in righe])
    mediane = {}
    for asse in ("azimut", "polare"):
        for sch in ("measured", "cascade", "lqg"):
            sel = [r for r in righe if r["axis"] == asse and r["scheme"] == sch and (r["axis"], r["test"], r["move"]) != ANOMALA]
            mediane.setdefault(asse, {})[sch] = {k: float(np.median([r[k] for r in sel])) for k in righe[0] if k not in ("axis", "test", "move", "scheme")}
    print(json.dumps(mediane, indent=1))
    caccia = {}
    for asse, chiave, cat, mod in assi:
        b = costruisci(asse, chiave, cat, mod, Kt, S, T, mot)
        caccia[chiave] = verifica_hunting(b, finestre_asse(asse), controllore(b, tun[chiave]["best"]))
    (CARTELLA / "hunting.json").write_text(json.dumps(caccia, indent=1))
    (CARTELLA / "medians.json").write_text(json.dumps(mediane, indent=1))


def verifica_hunting(a: Asse, fin, finale: lqg.ControlloreLQG) -> list[dict]:
    """Doppietto da 0.1 gradi, k_p circa 18, k_i crescente: errore residuo e frazione di tempo in moto dopo il comando."""
    w = next(x for x in fin["01deg_minus01deg"] if x.nome[-2:] == "#1")
    i1 = maschere(w)[1]
    Pb = frf(varianti_cuscinetto(a)["stuck"], a.anello_v.B[:, 0], a.anello_v.uscite["p"], F_GRID)
    out = []
    candidati = [lqg.sintesi(a.Ap, a.Bp, a.Cp, PASSO_MASTER, pesi_da_guadagni(18, ra)[1], pesi_da_guadagni(18, ra)[0], np.zeros(len(a.Ap) - 1))
                 for ra in (0.03, 0.1, 0.35, 1.0)] + [finale]
    for c in candidati:
        r = simula(a, c, w)
        e = np.rad2deg((w.cmd - w.cmd[0]) - r["p"])[i1:] * 1e3
        v = np.gradient(r["p"], dt.DT)[i1:]
        out.append({"final": c is finale, "k_p": c.k_p, "k_i": c.k_i, "e_max_last_3s_mdeg": float(np.abs(e[-6000:]).max()), "moving_fraction": float((np.abs(v) > 1e-5).mean()),
                    "sign_changes": int(np.sum(np.diff(np.sign(e)) != 0)), "modulus_stuck": margini(anello_lqg(c, Pb, F_GRID)[0], F_GRID)["modulus"]})
    return out


def zeta_flessibili_v(a) -> dict[str, float]:
    """Smorzamento dei modi nella maglia di velocita' chiusa, senza anello di posizione."""
    lam = np.linalg.eigvals(a.anello_v.A)
    out = {}
    for f0 in a.m.omega / (2 * np.pi):
        if 5.0 < f0 < 150.0:
            k = lam[lam.imag > 0]
            d = k[np.argmin(np.abs(np.abs(k) / (2 * np.pi) - f0) + 50 * (k.real / np.abs(k) < -0.12))]
            out[f"{f0:.1f}"] = float(-d.real / abs(d))
    return out


if __name__ == "__main__":
    main()
