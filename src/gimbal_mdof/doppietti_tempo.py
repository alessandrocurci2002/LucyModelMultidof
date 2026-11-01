"""Simulazione nel tempo dei modelli modali della forcella con la corrente misurata dei doppietti V2 come ingresso.

Prove: ``LucyIdentification/VClosedLoop/V2/MPCLucy_doublets_*.csv`` (tempo in ms, posizioni in mdeg, corrente
``AXIS*_AMPS``; AXIS1 = azimut, AXIS2 = polare). La colonna AMPS vale 2 I_rms (convenzione di ``closedLoopModel``,
``amps_per_A`` = 2). Ogni prova ha tre movimenti (+A, -2A, +A); ogni movimento e' una finestra.

Modello: la struttura modale senza cuscinetti (``ModelloModale`` con k_g = 0) piu' l'attrito dei cuscinetti sulla
deformazione v^T x (Dahl con curva di Stribeck; Dahl puro con ``Ts_ratio`` = 1). La parte lineare e' discretizzata
esattamente (ZOH, passo 0.5 ms), lo stato d'attrito e' aggiornato esattamente a velocita' costante sul passo. Nella
condizione iniziale l'attrito (pre-scorrimento) regge la coppia a riposo. Ingresso filtrato a 150 Hz (come lo studio
``current_input_study``): il modo dell'albero a ~950 Hz ha zeta al limite inferiore e il rumore della corrente lo
eccita.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import expm, lstsq
from scipy.signal import butter, filtfilt

from . import forcella_modale as fm

DT = 5e-4
AMPS_PER_A = 2.0
COLONNE = {"azimut": ("AzimuthAxis/ActualPosition", "AzimuthAxis/CommandPosition", "VM/VARIABLESTOMONITOR.AXIS1_AMPS"),
           "polare": ("PolarAxis/ActualPosition", "PolarAxis/CommandPosition", "VM/VARIABLESTOMONITOR.AXIS2_AMPS")}
_FILTRO = butter(4, 150 * 2 * DT)


@dataclass
class Finestra:
    nome: str
    t: np.ndarray
    I: np.ndarray
    """Corrente, A rms."""
    pos: np.ndarray
    """Encoder, rad, relativo all'inizio della finestra."""
    cmd: np.ndarray


def carica(percorso: Path, asse: str) -> dict[str, np.ndarray]:
    with open(percorso) as fh:
        cols = [c.strip() for c in fh.readline().split(",")]
    A = np.genfromtxt(percorso, delimiter=",", skip_header=1, usecols=range(len(cols)), filling_values=np.nan)
    ip, ic, ii = (cols.index(c) for c in COLONNE[asse])
    A = A[~np.isnan(A[:, [ip, ic, ii]]).any(axis=1)]
    return {"t": (A[:, 0] - A[0, 0]) * 1e-3, "pos": np.deg2rad(A[:, ip] * 1e-3), "cmd": np.deg2rad(A[:, ic] * 1e-3),
            "I": A[:, ii] / AMPS_PER_A}


def finestre(percorso: Path, asse: str, prima: float = 0.5, dopo: float = 12.0) -> list[Finestra]:
    """Una finestra per movimento del comando: da ``prima`` s prima a ``dopo`` s dopo (o fino al movimento seguente)."""
    d = carica(percorso, asse)
    t = d["t"]
    idx = np.flatnonzero(np.abs(np.diff(d["cmd"])) > 1e-12)
    inizi = idx[np.r_[True, np.diff(idx) > 200]]
    out = []
    for k, s in enumerate(inizi):
        t0 = t[s] - prima
        t1 = min(t[s] + dopo, t[inizi[k + 1]] - prima if k + 1 < len(inizi) else t[-1])
        m = (t >= t0) & (t <= t1)
        tt = np.arange(t0, t1, DT)
        I = filtfilt(*_FILTRO, np.interp(tt, t[m], d["I"][m]))
        p, c = np.interp(tt, t[m], d["pos"][m]), np.interp(tt, t[m], d["cmd"][m])
        out.append(Finestra(f"{percorso.stem.removeprefix('MPCLucy_doublets_')} #{k + 1}", tt - tt[0], I, p - p[0],
                            c - p[0]))
    return out


@dataclass
class Attrito:
    sigma0: float
    Tc: float
    Ts_ratio: float = 1.0
    vs: float = 1.0

    def g(self, v: float) -> float:
        return self.Tc * (1 + (self.Ts_ratio - 1) * np.exp(-(v / self.vs) ** 2))


class Simulatore:
    """Struttura modale (senza cuscinetti) + attrito sulla deformazione dei cuscinetti, coppia motore in ingresso."""

    def __init__(self, m: fm.ModelloModale, dt: float = DT):
        n = len(m.M)
        C = m.M @ m.Phi @ np.diag(2 * m.zeta * m.omega) @ m.Phi.T @ m.M
        K0 = m.K - m.k_g * np.outer(m.v, m.v)
        Mi = np.linalg.inv(m.M)
        A = np.block([[np.zeros((n, n)), np.eye(n)], [-Mi @ K0, -Mi @ C]])
        Bf = np.vstack([np.zeros((n, n)), Mi]) @ np.column_stack([m.e, -m.v])
        E = expm(np.block([[A, Bf], [np.zeros((2, 2 * n + 2))]]) * dt)
        self.Ad, self.Bd = E[:2 * n, :2 * n], E[:2 * n, 2 * n:]
        self.n, self.m, self.dt, self.K0 = n, m, dt, K0

    def corri(self, I: np.ndarray, att: Attrito) -> np.ndarray:
        """Posizione dell'encoder (rad, relativa all'inizio) per la corrente I (A rms)."""
        n, m, dt = self.n, self.m, self.dt
        tau = m.Kt * I
        tau0 = tau[:int(0.4 / dt)].mean()
        z = np.clip(tau0 / att.sigma0, -att.g(0) / att.sigma0, att.g(0) / att.sigma0)
        x0 = lstsq(np.vstack([self.K0, m.v[None, :]]), np.r_[m.e * tau0 - m.v * att.sigma0 * z, 0.0])[0]
        s = np.r_[x0, np.zeros(n)]
        out = np.empty(len(tau))
        Ad, Bd, v, c, s0 = self.Ad, self.Bd, m.v, m.c, att.sigma0
        for k in range(len(tau)):
            out[k] = c @ s[:n]
            s_new = Ad @ s + Bd @ np.array([tau[k], s0 * z])
            vd = 0.5 * (v @ s[n:] + v @ s_new[n:])
            a = s0 * abs(vd) / att.g(vd)
            if a > 0:
                e = np.exp(-a * dt)
                z = z * e + vd / a * (1 - e)
            s = s_new
        return out - out[0]


def fit(y: np.ndarray, y_mis: np.ndarray) -> float:
    """NRMSE, %."""
    return float(100 * (1 - np.linalg.norm(y - y_mis) / np.linalg.norm(y_mis - y_mis.mean())))
