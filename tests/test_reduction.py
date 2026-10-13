"""Test della sezione 13, punti 10-12: riduzione di Guyan."""

from __future__ import annotations

import numpy as np
import pytest

from gimbal_mdof import fem
from gimbal_mdof import params as pr
from gimbal_mdof import reduction as rd


def _modelli(par: pr.Parametri):
    return {nome: fem.costruisci_asse(par.asse(nome)) for nome in ("azimut", "polare")}


# --- test 10 ----------------------------------------------------------------


@pytest.mark.parametrize("nome", ["azimut", "polare"])
def test_forma_chiusa_catena(par: pr.Parametri, nome: str) -> None:
    mod = fem.costruisci_asse(par.asse(nome))
    assert all(b == a + 1 for a, b in mod.collegamenti)  # e' una catena
    rid = rd.riduci(mod, "consistent")
    chiusa = rd.forma_chiusa_catena([p.k for p in mod.elementi], mod.masters)
    assert np.allclose(rid.K_r, chiusa.K_r(), rtol=1e-11, atol=1e-9 * np.max(rid.K_r))
    assert np.allclose(rid.T, chiusa.T(mod.n, mod.masters), rtol=0, atol=1e-12)
    assert np.allclose(rid.k_catena, chiusa.k_serie, rtol=1e-11)
    # K_r e M_r tridiagonali
    for A in (rid.K_r, rid.M_r):
        assert np.allclose(np.triu(A, 2), 0.0, atol=1e-14 * np.max(np.abs(A)))


def test_rigidezze_catena_polare(par: pr.Parametri) -> None:
    """k_theta,1: el. 1, 2, 3a; k_theta,2: 3b, 4, 5, 6; k_theta,3: 7..10 (sezione 4.2)."""
    mod = fem.costruisci_asse(par.polare)
    rid = rd.riduci(mod, "lumped")
    k = {p.etichetta: p.k for p in mod.elementi}
    serie = lambda nomi: 1.0 / sum(1.0 / k[n] for n in nomi)  # noqa: E731
    attese = [serie(["1", "2", "3a"]), serie(["3b", "4", "5", "6"]),
              serie(["7", "8", "9", "10"])]
    assert rid.ruoli == ("e", "p", "m", "r")
    assert np.allclose(rid.k_catena, attese, rtol=1e-12)


# --- test 11 ----------------------------------------------------------------


@pytest.mark.parametrize("nome", ["azimut", "polare"])
@pytest.mark.parametrize("forma", ["guyan", "guyan_lumped"])
def test_cedevolezza_statica_e_inerzia(par: pr.Parametri, nome: str, forma: str) -> None:
    mod = fem.costruisci_asse(par.asse(nome))
    rid = rd.riduci(mod, "consistent", forma)
    # cedevolezza fra i master con il primo master vincolato
    liberi_full = [i for i in range(mod.n) if i != mod.masters[0]]
    F_full = np.linalg.inv(mod.K[np.ix_(liberi_full, liberi_full)])
    idx = [liberi_full.index(i) for i in mod.masters[1:]]
    F_full_m = F_full[np.ix_(idx, idx)]
    F_red = np.linalg.inv(rid.K_r[1:, 1:])
    assert np.allclose(F_red, F_full_m, rtol=1e-10)
    uno_r, uno = np.ones(rid.n), np.ones(mod.n)
    assert uno_r @ rid.M_r @ uno_r == pytest.approx(uno @ mod.M["consistent"] @ uno, rel=1e-12)


# --- test 12 ----------------------------------------------------------------


@pytest.mark.parametrize("nome", ["azimut", "polare"])
@pytest.mark.parametrize("modello", ["lumped", "consistent"])
@pytest.mark.parametrize("k_agg", [0.0, 50.0, -5.0])
def test_rayleigh_ritz(par: pr.Parametri, nome: str, modello: str, k_agg: float) -> None:
    """Con mass_form: guyan, omega_k^r >= omega_k^full (anche con k_g al nodo P)."""
    mod = fem.costruisci_asse(par.asse(nome))
    rid = rd.riduci(mod, modello, "guyan")
    k_agg = k_agg if nome == "polare" else 0.0
    M, K = fem.matrici_isolate(mod, modello, J_aggiunta=1.0, k_aggiunta=k_agg)
    Mr, Kr = rd.matrici_ridotte_isolate(rid, J_aggiunta=1.0, k_aggiunta=k_agg)
    lam_full = np.linalg.eigvalsh(np.linalg.solve(np.linalg.cholesky(M), np.linalg.solve(
        np.linalg.cholesky(M), K).T))
    lam_red = np.sort(np.real(np.linalg.eigvals(np.linalg.solve(Mr, Kr))))
    lam_full = np.sort(lam_full)[: rid.n]
    scala = np.max(np.abs(lam_red))
    assert np.all(lam_red >= lam_full - 1e-9 * scala)


@pytest.mark.parametrize("nome", ["azimut", "polare"])
def test_craig_bampton(par: pr.Parametri, nome: str) -> None:
    mod = fem.costruisci_asse(par.asse(nome))
    M, K = fem.matrici_isolate(mod, "consistent", J_aggiunta=1.0)
    rid = rd.riduci(mod, "consistent", "guyan")
    _, M0, K0 = rd.riduci_craig_bampton(M, K, mod.masters, 0)
    Mr, Kr = rd.matrici_ridotte_isolate(rid, J_aggiunta=1.0)
    assert np.allclose(M0, Mr, rtol=1e-12, atol=1e-18) and np.allclose(K0, Kr, rtol=1e-12)
    # convergenza: aggiungendo modi a interfaccia fissa le prime frequenze scendono
    # verso quelle del modello completo
    f_full = fem.frequenze_proprie(M, K)[: rid.n]
    errori = []
    for n in (0, 1, 2, 4):
        _, Mc, Kc = rd.riduci_craig_bampton(M, K, mod.masters, n)
        f = fem.frequenze_proprie(Mc, Kc)[: rid.n]
        assert np.all(f >= f_full - 1e-9 * f_full.max())
        errori.append(np.max(np.abs(f[1:] / f_full[1:] - 1)))
    assert all(b <= a + 1e-15 for a, b in zip(errori[:-1], errori[1:]))
