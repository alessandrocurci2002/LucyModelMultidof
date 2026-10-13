"""Test della sezione 13, punti 4-9: modello FEM di ciascun asse."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from gimbal_mdof import fem
from gimbal_mdof import geometry as geo
from gimbal_mdof import params as pr

from conftest import CONFIG, config_dati, parametri

ACCIAIO = pr.Materiale(nome="acciaio", rho=8000.0, G=77e9, segnaposto=False)


def _uniforme(L: float, D: float) -> pr.Elemento:
    seg = pr.Segmento(lunghezza=L, strati=(pr.Strato("solid", D, 0.0, ACCIAIO, True),))
    return pr.Elemento(id=1, modo="geometry", lunghezza=L, segmenti=(seg,), k_joint=None,
                       extra=(), J_manuale=None, k_manuale=None, m_manuale=None,
                       I_t_manuale=None, y_c_manuale=None, segnaposto=False)


# --- test 4 -----------------------------------------------------------------


def test_matrici_elemento_uniforme() -> None:
    p = geo.proprieta_elemento(_uniforme(0.03, 0.02))
    J = p.J
    assert np.allclose(fem.matrice_massa_consistente(p), J / 6 * np.array([[2, 1], [1, 2]]),
                       rtol=1e-14, atol=0)
    assert np.allclose(fem.matrice_massa_concentrata(p), J / 2 * np.eye(2), rtol=1e-14, atol=0)
    assert np.allclose(fem.matrice_massa_consistente(p, "compliance"),
                       J / 6 * np.array([[2, 1], [1, 2]]), rtol=1e-14, atol=0)


def test_consistente_a_tratti_integrale_numerico() -> None:
    """La forma chiusa segmento per segmento coincide con la quadratura."""
    par = pr.parametri_da_config(CONFIG)
    p = geo.proprieta_elemento(par.azimut.elemento(3))  # due tratti diversi
    L = p.lunghezza
    x = np.linspace(0.0, L, 200001)
    Jbar = np.zeros_like(x)
    for s in p.segmenti:
        Jbar[(x >= s.x_inizio) & (x <= s.x_inizio + s.lunghezza)] = s.J_bar
    N = np.vstack([1 - x / L, x / L])
    M_num = np.array([[np.trapezoid(N[a] * N[b] * Jbar, x) for b in range(2)] for a in range(2)])
    assert np.allclose(fem.matrice_massa_consistente(p), M_num, rtol=1e-4)
    assert fem.matrice_massa_consistente(p).sum() == pytest.approx(p.J, rel=1e-14)


# --- test 5 -----------------------------------------------------------------


@pytest.mark.parametrize("nome", ["azimut", "polare"])
@pytest.mark.parametrize("modello", ["lumped", "consistent"])
def test_proprieta_globali(par: pr.Parametri, nome: str, modello: str) -> None:
    mod = fem.costruisci_asse(par.asse(nome))
    M, K = mod.M[modello], mod.K
    assert np.allclose(M, M.T) and np.allclose(K, K.T)
    assert np.max(np.abs(K @ np.ones(mod.n))) < 1e-9 * np.max(np.abs(K))
    totale = sum(p.J for p in mod.elementi_originali) + sum(n.J for n in par.asse(nome).inerzie_nodali)
    assert np.ones(mod.n) @ M @ np.ones(mod.n) == pytest.approx(totale, rel=1e-13)
    assert fem.numero_modi_rigidi(K, M) == 1


# --- test 6 -----------------------------------------------------------------


def test_ordinamento_e_identita_per_elementi_uniformi() -> None:
    rng = np.random.default_rng(1)
    k = rng.uniform(1e3, 1e5, 9)
    J = rng.uniform(1e-6, 1e-4, 9)
    Ml, K = fem.catena_uniforme(k, J, "lumped")
    Mc, _ = fem.catena_uniforme(k, J, "consistent")
    assert np.min(np.linalg.eigvalsh(Ml - Mc)) > -1e-18
    for _ in range(5):
        v = rng.normal(size=10)
        atteso = sum(J[i] / 6 * (v[i] - v[i + 1]) ** 2 for i in range(9))
        assert v @ (Ml - Mc) @ v == pytest.approx(atteso, rel=1e-12)
    f_l = fem.frequenze_proprie(Ml, K)
    f_c = fem.frequenze_proprie(Mc, K)
    assert np.all(f_l <= f_c * (1 + 1e-12))


def test_asimmetria_elemento_determinante() -> None:
    """det(M^l - M^c) = -delta^2: l'ordinamento non vale per elementi asimmetrici."""
    par = pr.parametri_da_config(CONFIG)
    p = geo.proprieta_elemento(par.azimut.elemento(3))
    D = fem.matrice_massa_concentrata(p) - fem.matrice_massa_consistente(p)
    assert np.linalg.det(D) == pytest.approx(-fem.asimmetria_elemento(p) ** 2, rel=1e-9)
    assert fem.asimmetria_elemento(geo.proprieta_elemento(_uniforme(0.03, 0.02))) == pytest.approx(0.0, abs=1e-20)


# --- test 7 -----------------------------------------------------------------


def _polare(xi: float, clamp_mm: float = 0.0, *, manuale: bool = False):
    dati = config_dati()
    dati["polar"]["payload_xi"] = xi
    dati["polar"]["payload_clamp_rigid_mm"] = clamp_mm
    if manuale:
        dati["polar"]["elements"][2] = {"id": 3, "mode": "manual", "J_kgm2": 2.0e-5,
                                        "k_Nm_rad": 1.0e5, "length_mm": 60.0, "m_kg": 0.2,
                                        "I_t_kgm2": 6e-5, "y_c_mm": 30.0}
    par = parametri(dati)
    return par, fem.costruisci_asse(par.polare)


@pytest.mark.parametrize("manuale", [False, True])
@pytest.mark.parametrize("xi", [0.3, 0.5, 0.71])
def test_nodo_P_conserva_rigidezza_e_inerzia(xi: float, manuale: bool) -> None:
    par, mod = _polare(xi, manuale=manuale)
    p3 = geo.proprieta_elemento(par.polare.elemento(3))
    a, b = (p for p in mod.elementi if p.etichetta in ("3a", "3b"))
    assert 1 / a.k + 1 / b.k == pytest.approx(1 / p3.k, rel=1e-13)
    assert a.J + b.J == pytest.approx(p3.J, rel=1e-13)
    assert mod.etichette[mod.i_P] == "P" and mod.n == 12
    x_P_locale = mod.x[mod.i_P] - mod.x[mod.indice(3)]
    assert x_P_locale == pytest.approx(xi * p3.lunghezza, rel=1e-12)


@pytest.mark.parametrize("xi, alias", [(0.0, "3"), (1.0, "4")])
def test_nodo_P_agli_estremi(xi: float, alias: str) -> None:
    _, mod = _polare(xi)
    assert mod.n == 11 and "P" not in mod.etichette
    assert mod.etichette[mod.i_P] == alias


def test_nodo_P_massa_concentrata_meta() -> None:
    par, mod = _polare(0.5)
    J3 = geo.proprieta_elemento(par.polare.elemento(3)).J
    assert mod.M["lumped"][mod.i_P, mod.i_P] == pytest.approx(J3 / 2, rel=1e-13)
    M_iso, _ = fem.matrici_isolate(mod, "lumped", J_aggiunta=1.234)
    assert M_iso[mod.i_P, mod.i_P] == pytest.approx(J3 / 2 + 1.234, rel=1e-13)


def test_tratto_rigido_del_morsetto() -> None:
    par, mod0 = _polare(0.5, 0.0)
    _, mod = _polare(0.5, 10.0)
    k0 = [p.k for p in mod0.elementi if p.etichetta in ("3a", "3b")]
    k1 = [p.k for p in mod.elementi if p.etichetta in ("3a", "3b")]
    assert all(b > a for a, b in zip(k0, k1))
    assert sum(p.J for p in mod.elementi) == pytest.approx(sum(p.J for p in mod0.elementi), rel=1e-13)


# --- test 8 -----------------------------------------------------------------


def test_cuscinetti() -> None:
    dati = config_dati()
    dati["azimuth"]["bearing_b"] = {3: 0.02, 6: 0.04}
    par = parametri(dati)
    mod = fem.costruisci_asse(par.azimut)
    C_b = mod.C_b
    assert C_b[2, 2] == pytest.approx(0.01) and C_b[3, 3] == pytest.approx(0.01)
    assert C_b[5, 5] == pytest.approx(0.02) and C_b[6, 6] == pytest.approx(0.02)
    assert np.count_nonzero(C_b) == 4
    M = mod.M["consistent"]
    C = fem.matrice_smorzamento(mod.K, C_b, fem.beta_smorzamento(M, mod.K, 0.01))
    n = mod.n
    A = np.block([[np.zeros((n, n)), np.eye(n)],
                  [-np.linalg.solve(M, mod.K), -np.linalg.solve(M, C)]])
    poli = np.linalg.eigvals(A)
    reali = np.sort(np.real(poli[np.abs(poli) < 1.0]))
    atteso = -0.06 / mod.J_totale
    assert len(reali) == 2
    # il polo "nullo" e' tale a meno del rumore di eigvals su autovalori ~1e5 rad/s
    scala = np.max(np.abs(poli))
    assert abs(reali[1]) < 1e-11 * scala
    assert reali[0] == pytest.approx(atteso, rel=1e-6)


# --- test 9 -----------------------------------------------------------------


def test_convergenza_barra_libera() -> None:
    L, D, rho, G = 0.3, 0.02, ACCIAIO.rho, ACCIAIO.G
    omega_an = np.array([n * np.pi / L * np.sqrt(G / rho) for n in (1, 2, 3)])
    errori_l, errori_c = [], []
    for n_el in (6, 12, 24, 48):
        p = geo.proprieta_elemento(_uniforme(L / n_el, D))
        k, J = [p.k] * n_el, [p.J] * n_el
        for modello, errori in (("lumped", errori_l), ("consistent", errori_c)):
            M, K = fem.catena_uniforme(k, J, modello)
            w = 2 * np.pi * fem.frequenze_proprie(M, K)[1:4]
            if modello == "lumped":
                assert np.all(w < omega_an)
            else:
                assert np.all(w > omega_an)
            errori.append(np.max(np.abs(w / omega_an - 1)))
    assert all(np.diff(errori_l) < 0) and all(np.diff(errori_c) < 0)
    assert errori_l[-1] < 2e-3 and errori_c[-1] < 2e-3


# --- telescopio senza albero (mode: telescope) -------------------------------


def _polare_telescopio(**polar: object):
    dati = config_dati()
    dati["polar"]["elements"][2] = {"id": 3, "mode": "telescope", "length_mm": 60.0}
    dati["polar"].update(polar)
    par = parametri(dati)
    return par, fem.costruisci_asse(par.polare)


def test_telescopio_senza_albero_nodi_coincidenti() -> None:
    par, mod = _polare_telescopio(bearing_b={2: 0.02, 4: 0.04})
    assert mod.n == 10 and mod.telescopio_rigido
    assert mod.indice(3) == mod.indice(4) == mod.indice("P") == mod.i_P
    assert "3" not in mod.etichette and "4" not in mod.etichette
    assert all(p.etichetta != "3" for p in mod.elementi)
    # inerzia: solo alberi e nodi, nessun contributo dall'elemento 3
    J_alberi = sum(p.J for p in mod.elementi_originali)
    assert mod.J_totale == pytest.approx(J_alberi + sum(n.J for n in par.polare.inerzie_nodali),
                                         rel=1e-13)
    # P raccoglie l'inerzia concentrata di meta' elemento 2 e meta' elemento 4
    J2, J4 = (geo.proprieta_elemento(par.polare.elemento(e)).J for e in (2, 4))
    assert mod.M["lumped"][mod.i_P, mod.i_P] == pytest.approx((J2 + J4) / 2, rel=1e-13)
    # i cuscinetti degli elementi 2 e 4 agiscono entrambi su P
    assert mod.C_b[mod.i_P, mod.i_P] == pytest.approx(0.01 + 0.02)
    # posizione di P: centro del telescopio
    x3 = sum(par.polare.elemento(e).lunghezza for e in (1, 2))
    assert mod.x[mod.i_P] == pytest.approx(x3 + 0.030, rel=1e-12)


def test_telescopio_senza_albero_catena_ridotta() -> None:
    from gimbal_mdof import reduction as rd

    _, mod = _polare_telescopio()
    rid = rd.riduci(mod, "consistent")
    k = {p.etichetta: p.k for p in mod.elementi}
    serie = lambda nomi: 1.0 / sum(1.0 / k[n] for n in nomi)  # noqa: E731
    assert rid.ruoli == ("e", "p", "m", "r")
    assert np.allclose(rid.k_catena, [serie(["1", "2"]), serie(["4", "5", "6"]),
                                      serie(["7", "8", "9", "10"])], rtol=1e-12)


def test_telescopio_validazione() -> None:
    dati = config_dati()
    dati["polar"]["elements"][2] = {"id": 3, "mode": "telescope", "length_mm": 60.0,
                                    "J_kgm2": 1e-3}
    with pytest.raises(pr.ErroreConfigurazione, match="payload"):
        parametri(dati)
    dati = config_dati()
    dati["polar"]["elements"][4] = {"id": 5, "mode": "telescope", "length_mm": 30.0}
    with pytest.raises(pr.ErroreConfigurazione, match="telescope"):
        parametri(dati)
