"""Test della sezione 13, punto 20: scenari S1 e S3."""

from __future__ import annotations

import numpy as np
import pytest

from gimbal_mdof import coupled as cp
from gimbal_mdof import params as pr
from gimbal_mdof import simulate as sm


@pytest.fixture(scope="module")
def acc(par: pr.Parametri) -> cp.ModelloAccoppiato:
    return cp.costruisci_sistema(par).accoppiato


@pytest.mark.lento
def test_S1_confronto_lineare(acc: cp.ModelloAccoppiato) -> None:
    """Canali dell'asse forzato entro l'1% dal lineare.

    I canali polari, a Omega_0 = 0, sono risposte incrociate dovute ai soli
    prodotti d'inerzia: il termine quadratico (1/2) J' dpsi^2 vi e' dello stesso
    ordine, quindi non sono un banco di prova della linearizzazione.
    """
    r = sm.S1(acc, np.deg2rad(15.0), t_fine=0.1, n_punti=501)
    for u in ("psi_e", "psi_m", "psi_L", "err_psi"):
        assert r.errori[u] < 1e-2, f"{u}: errore relativo {r.errori[u]:.3e}"


@pytest.mark.lento
def test_S3_conservazione(acc: cp.ModelloAccoppiato, par: pr.Parametri) -> None:
    r = sm.S3(acc, np.deg2rad(15.0), par.omega_max, t_fine=0.1, n_punti=201)
    assert r.derive["energia"] < 1e-8
    assert r.derive["momento_psi"] < 1e-8
