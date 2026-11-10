"""LQG-PI con la parte proporzionale-integrale sull'errore misurato (Gawronski 2008, Fig. 8.4b).

Richiesta dell'utente (2026-11-05): tenere l'implementazione della Fig. 8.4b, con il PI sull'errore misurato e = r - y (baseline
robusta, che non dipende dal modello interno) e la retroazione degli stati flessibili stimati, invece della legge (8.15) sullo
stato stimato usata finora (``lucy_cl.lqg.ControlloreLQG``, e = r - y_hat):

    u = k_i e_i + k_p e - K_f x_f_hat,   e = r - y,   e_i = sum T e.

Guadagni e stimatore sono quelli di ``lqg.sintesi`` (stessa CARE e stessa FARE): cambia solo il segnale su cui agisce il PI.
Lo stimatore resta quello predittivo della (8.13), alimentato da y e dal comando applicato.

Le simulazioni passano y(k) con ``misura`` prima di ``uscita`` (aggancio in ``g.corri_lqg``, ``cc.simula_coppia`` e
``p1.Banco.gradino``). ``attiva()`` sostituisce, nel solo processo che la chiama, ``lqg.sintesi`` e le funzioni dei poli in
anello chiuso, in modo che taratura e simulazioni esistenti usino questo controllore senza essere modificate.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

QUI = Path(__file__).resolve().parent
RADICE = QUI.parent
sys.path[:0] = [str(RADICE / "src"), str(RADICE), str(RADICE / "closedLoopModel")]

from scipy.linalg import expm  # noqa: E402

from lucy_cl import lqg  # noqa: E402

_SINTESI = lqg.sintesi


ControlloreLQGPIMisurato = lqg.ControlloreLQGPIMisurato
"""Dal 2026-11-05 la classe e' in ``lucy_cl.lqg``, perche' e' il controllore della tesi."""
da_stimato = lqg.pi_misurato


def sintesi(*a, **k) -> ControlloreLQGPIMisurato:
    return da_stimato(_SINTESI(*a, **k))


def poli_chiusi(c, A, B, C, ritardo: int = 1) -> np.ndarray:
    """Poli discreti dell'anello chiuso, stati [x, x_hat, e_i, u in coda], come ``p1.poli_chiusi``; per questo controllore il PI
    legge y = C x invece di y_hat."""
    if not isinstance(c, ControlloreLQGPIMisurato):
        if _POLI is None:  # prima di ``attiva``: la funzione originale
            import analisi_lqg_passo1_flessibili as p1
            return p1.poli_chiusi(c, A, B, C, ritardo)
        return _POLI(c, A, B, C, ritardo)
    N, n = len(A), len(c.Ad)
    E = expm(np.block([[A, B[:, None]], [np.zeros((1, N + 1))]]) * c.passo)
    tot = N + n + 1 + ritardo
    M = np.zeros((tot, tot))
    ie, iq = N + n, N + n + 1
    M[:N, :N], M[:N, iq] = E[:N, :N], E[:N, N]
    M[N:N + n, :N] = c.K_e[:, [0]] @ C[None, :]
    M[N:N + n, N:N + n] = c.Ad - c.K_e[:, [0]] @ np.eye(n)[:1]
    M[N:N + n, iq] = c.Bd[:, 0]
    M[ie, ie], M[ie, :N] = 1.0, -c.passo * C
    M[iq, ie], M[iq, :N] = c.k_i, -c.k_p * C
    M[iq, N + 1:N + n] = -c.K_f
    return np.linalg.eigvals(M)


_POLI = None


def attiva() -> None:
    """Sostituisce, in questo processo, la sintesi e i poli in anello chiuso usati dagli script della taratura e delle simulazioni."""
    global _POLI
    import analisi_lqg_gawronski as g
    import analisi_lqg_passo1_flessibili as p1
    if _POLI is None:
        _POLI = p1.poli_chiusi
    lqg.sintesi = sintesi
    p1.poli_chiusi = poli_chiusi
    g.autovalori_chiusi = lambda a, c, A: poli_chiusi(c, A, a.anello_v.B[:, 0], a.anello_v.uscite["p"], g.RITARDO)
