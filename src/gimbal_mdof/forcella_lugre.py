"""Modello modale della forcella con i cuscinetti ad attrito di LuGre (funzione descrittiva) al posto di k_g (1 + j eta_g).

Richiesta dell'utente (2026-10-31): stessa formulazione e stessi vincoli della stima salvata in
``saved_fork_modes_modal_20261031_cad20_kt10`` (zeta_k <= 0.05, CAD +-20 %, K_t +-10 %, minimi quadrati sulle parti
reale e immaginaria), con l'attrito di LuGre sulla deformazione dei cuscinetti v^T x:

    dz/dt = v_b - sigma0 |v_b| z / g(v_b),   g(v) = T_c + (T_s - T_c) exp(-(v/v_s)^2),
    T_b = sigma0 z + sigma1 dz/dt + sigma2 v_b.

L'attrito di Dahl (richiesta successiva, stesso giorno) e' il caso T_s = T_c, sigma1 = sigma2 = 0: restano liberi
sigma0 e T_c (:func:`parametri_dahl`, ``famiglia="dahl"``).

In frequenza l'attrito entra con la prima armonica (funzione descrittiva) K_b(A, f) = T_b / x_b sotto
x_b = A sin(omega t). L'ampiezza A del cuscinetto in ogni punto e' quella della prova:

    A(f) = |G_mis(f)| I_exc |S_mis(f)| / omega,   S_mis = 1 / (1 + L_mis),

con I_exc la corrente iniettata e L_mis l'anello aperto misurato (la prova e' in anello chiuso: nel motore circola
I_exc S). Il moto dello statore sotto i modi della forcella e' trascurabile, quindi la rotazione dell'asse letta
dall'encoder e' la deformazione del cuscinetto dove l'attrito conta.

La parte strutturale e' quella di :mod:`gimbal_mdof.forcella_modale` senza il vincolo dei cuscinetti (k_g = 0): il
primo modo diventa il moto rigido dell'asse e i cuscinetti entrano con la correzione di rango 1 (Sherman-Morrison) con
h = K_b(A, f), esatta come per la parte isteretica del modello lineare.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from . import forcella_modale as fm
from .modi_forcella import CatenaAzimut, CatenaPolare, DatoFRF, Parametro, _misurata

NOMI_LUGRE = ("sigma0", "Tc", "Ts_ratio", "vs", "sigma1", "sigma2")
DAHL = {"Ts_ratio": 1.0, "vs": 1.0, "sigma1": 0.0, "sigma2": 0.0}
"""Valori fissi che riducono LuGre a Dahl (con T_s = T_c la velocita' di Stribeck non conta)."""
PASSI_CICLO = 256
CICLI = 3


def funzione_descrittiva(p: Mapping[str, float], A: np.ndarray, f: np.ndarray, n: int = PASSI_CICLO,
                         cicli: int = CICLI) -> np.ndarray:
    """Prima armonica T_b / x_b di LuGre sotto x_b = A sin(omega t), per tutti i punti insieme (regime: ultimo ciclo).

    Stato z aggiornato esattamente a velocita' costante sul passo; n passi per ciclo."""
    A, f = np.asarray(A, float), np.asarray(f, float)
    w = 2 * np.pi * f
    dt = 1.0 / (n * f)
    s0, Tc, s1, s2 = p["sigma0"], p["Tc"], p["sigma1"], p["sigma2"]
    Ts, vs = Tc * p["Ts_ratio"], p["vs"]
    z = np.zeros_like(A)
    num = np.zeros(A.shape, complex)
    for k in range(cicli * n):
        ph = 2 * np.pi * (k + 0.5) / n
        v = A * w * np.cos(ph)
        g = Tc + (Ts - Tc) * np.exp(-(v / vs) ** 2)
        a = s0 * np.abs(v) / g
        e = np.exp(-a * dt)
        with np.errstate(divide="ignore", invalid="ignore"):
            z_new = np.where(a > 0, z * e + np.where(a > 0, v / a, 0.0) * (1 - e), z + v * dt)
        if k >= (cicli - 1) * n:
            zm = 0.5 * (z + z_new)  # z a meta' passo, come v e la fase
            T = s0 * zm + s1 * (v - a * zm) + s2 * v
            num += T * np.exp(-1j * ph)
        z = z_new
    # x = A sin(ph): sum x e^{-j ph} = -j A n / 2
    return num / (-0.5j * A * n)


def theta_tau(m: fm.ModelloModale, s: np.ndarray, h: np.ndarray) -> np.ndarray:
    """theta_enc / tau con la rigidezza complessa dei cuscinetti h(s) in rango 1 sulla base dei modi senza cuscinetti."""
    b, g, w, D = m._modale(s)
    G0 = (g * b / D).sum(axis=1)
    gw, wb, ww = (g * w / D).sum(axis=1), (w * b / D).sum(axis=1), (w * w / D).sum(axis=1)
    return G0 - h * gw * wb / (1 + h * ww)


@dataclass
class ModelloLuGre:
    """Struttura senza cuscinetti (``ModelloModale`` con k_g = 0) piu' i cuscinetti di LuGre; ampiezze per misura."""

    struttura: fm.ModelloModale
    lugre: dict[str, float]
    ampiezze: Mapping[str, np.ndarray]

    def K_b(self, d: DatoFRF) -> np.ndarray:
        return funzione_descrittiva(self.lugre, self.ampiezze[d.etichetta], d.f)

    def frf_dato(self, d: DatoFRF) -> np.ndarray:
        m = self.struttura
        s = 2j * np.pi * d.f
        return _misurata(theta_tau(m, s, self.K_b(d)), s, m.Kt, m.ritardo, m.f_corrente)

    def frf_linearizzata(self, f: np.ndarray, A: np.ndarray) -> np.ndarray:
        """Pianta con i cuscinetti alla rigidezza equivalente all'ampiezza A (per figure su griglie fitte)."""
        m = self.struttura
        s = 2j * np.pi * np.asarray(f, float)
        return _misurata(theta_tau(m, s, funzione_descrittiva(self.lugre, A, f)), s, m.Kt, m.ritardo, m.f_corrente)


def _strutturali(v: Mapping[str, float]) -> dict[str, float]:
    return {**v, "k_g": 0.0, "eta_g": 0.0, "zeta_1": fm.ZETA_MIN}


def _attrito(v: Mapping[str, float]) -> dict[str, float]:
    """Parametri d'attrito: quelli mancanti in ``v`` (modello di Dahl) prendono i valori di :data:`DAHL`."""
    return {k: v[k] if k in v else DAHL[k] for k in NOMI_LUGRE}


def polare(cat: CatenaPolare, v: Mapping[str, float], Kt: float, ampiezze: Mapping[str, np.ndarray]) -> ModelloLuGre:
    return ModelloLuGre(fm.polare(cat, _strutturali(v), Kt), _attrito(v), ampiezze)


def azimut(cat: CatenaAzimut, v: Mapping[str, float], Kt: float, ampiezze: Mapping[str, np.ndarray]) -> ModelloLuGre:
    return ModelloLuGre(fm.azimut(cat, _strutturali(v), Kt), _attrito(v), ampiezze)


def parametri_lugre() -> list[Parametro]:
    return [Parametro("sigma0", 10.0, 1e7, unita="N m/rad", nota="LuGre: bristle stiffness (free)"),
            Parametro("Tc", 1e-3, 10.0, unita="N m", nota="LuGre: Coulomb torque (free)"),
            Parametro("Ts_ratio", 1.0, 20.0, log=False, nota="LuGre: static / Coulomb torque (free)"),
            Parametro("vs", 1e-6, 1.0, unita="rad/s", nota="LuGre: Stribeck velocity (free)"),
            Parametro("sigma1", 1e-3, 1e5, unita="N m s/rad", nota="LuGre: bristle damping (free)"),
            Parametro("sigma2", 1e-5, 1e3, unita="N m s/rad", nota="LuGre: viscous friction (free)")]


def parametri_dahl() -> list[Parametro]:
    return [p for p in parametri_lugre() if p.nome in ("sigma0", "Tc")]


def _sostituisci(par: Sequence[Parametro], famiglia: str) -> list[Parametro]:
    """Parametri del modello lineare con l'attrito al posto di k_g, eta_g; zeta_1 (moto rigido) non c'e' piu'."""
    out = []
    for p in par:
        if p.nome == "k_g":
            out += parametri_lugre() if famiglia == "lugre" else parametri_dahl()
        elif p.nome not in ("eta_g", "zeta_1"):
            out.append(p)
    return out


def parametri_polari(famiglia: str = "lugre", **kw) -> list[Parametro]:
    return _sostituisci(fm.parametri_polari(**kw), famiglia)


def parametri_azimut(famiglia: str = "lugre", **kw) -> list[Parametro]:
    return _sostituisci(fm.parametri_azimut(**kw), famiglia)


def ampiezza_misurata(G_mis: np.ndarray, L_mis: np.ndarray, f: np.ndarray, I_exc: float) -> np.ndarray:
    """Ampiezza della rotazione dell'asse durante la prova, rad."""
    return np.abs(G_mis) * I_exc * np.abs(1 / (1 + L_mis)) / (2 * np.pi * f)
