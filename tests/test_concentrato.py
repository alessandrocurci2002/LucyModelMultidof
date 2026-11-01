"""Modello a parametri concentrati (3 + 4 GdL): elementi equivalenti, Newton-Eulero, FdT, simulatore."""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest
import sympy as sp

from gimbal_mdof import concentrato as cc
from gimbal_mdof import concentrato_analisi as ca
from gimbal_mdof import concentrato_output as co
from gimbal_mdof import coupled as cp
from gimbal_mdof import fem
from gimbal_mdof import linearize as lz
from gimbal_mdof import params as pr
from gimbal_mdof import payload as plm
from gimbal_mdof import simulatore_concentrato as scm
from gimbal_mdof import tf


@pytest.fixture(scope="module")
def assi(par: pr.Parametri) -> tuple[fem.ModelloAsse, fem.ModelloAsse]:
    return fem.costruisci_asse(par.azimut), fem.costruisci_asse(par.polare)


@pytest.fixture(scope="module")
def sistema(par: pr.Parametri, assi) -> cp.Sistema:
    return cc.costruisci_sistema_concentrato(par, plm.costruisci_payload(par, assi[1]), assi)


@pytest.fixture(scope="module")
def sim(sistema: cp.Sistema) -> scm.SimulatoreNE:
    return scm.costruisci_simulatore(sistema)


def test_elementi_equivalenti(assi) -> None:
    """J_eq = somma, 1/k_eq = somma delle cedevolezze; inerzia totale e statica conservate."""
    for mod in assi:
        ee = cc.elementi_equivalenti(mod)
        assert [e.nodo_a for e in ee] + [ee[-1].nodo_b] == [mod.etichetta(i) for i in mod.masters]
        tutti = [p for p in mod.elementi]
        assert sum(e.J_elementi for e in ee) == pytest.approx(sum(p.J for p in tutti), rel=1e-14)
        for e in ee:
            ks = [p.k for p in mod.elementi if p.etichetta in e.elementi]
            assert 1.0 / e.k == pytest.approx(sum(1.0 / k for k in ks), rel=1e-14)
        for forma in ("lumped", "consistent"):
            M, K = cc.matrici_concentrate(mod, forma)
            uno = np.ones(M.shape[0])
            assert float(uno @ M @ uno) == pytest.approx(mod.J_totale, rel=1e-12)
            assert np.max(np.abs(K @ uno)) < 1e-9 * np.max(np.abs(K))
        # la catena di molle equivalenti e' la condensazione statica esatta (K_r di Guyan)
        from gimbal_mdof import reduction as rd

        rg = rd.riduci(mod, "consistent")
        assert np.allclose(cc.matrici_concentrate(mod)[1], rg.K_r, rtol=1e-10, atol=1e-6)


def test_newton_eulero_uguale_a_lagrange() -> None:
    assert cc.verifica_contro_lagrange() == sp.zeros(7, 1)


def test_forma_compatta() -> None:
    co.forma_compatta()  # solleva AssertionError se diversa da Newton-Eulero


def test_simulatore_uguale_a_lagrange(sim: scm.SimulatoreNE, sistema: cp.Sistema) -> None:
    rng = np.random.default_rng(3)
    for _ in range(10):
        q, dq = rng.normal(size=(2, 7))
        tau = rng.normal(size=2)
        a1 = sim.accelerazioni(q, dq, tau)
        a2 = sistema.accoppiato.accelerazioni(q, dq, tau)
        assert np.max(np.abs(a1 - a2)) <= 1e-12 * np.max(np.abs(a2))


def test_codice_esportato(sim: scm.SimulatoreNE, tmp_path) -> None:
    file = scm.esporta_codice(sim, tmp_path)
    spec = importlib.util.spec_from_file_location("ne_concentrato", file[0])
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    rng = np.random.default_rng(5)
    for _ in range(5):
        x = rng.normal(size=14)
        tau = rng.normal(size=2)
        a = sim.accelerazioni(x[:7], x[7:], tau)
        dx = mod.rhs(0.0, x, lambda t: tau)
        assert np.max(np.abs(dx[7:] - a)) <= 1e-12 * np.max(np.abs(a))
    testo_m = file[1].read_text()
    assert "function dx = ne_concentrato" in testo_m and "sin(theta_p)" in testo_m


def test_zeri_forma_chiusa(sistema: cp.Sistema, par: pr.Parametri) -> None:
    """Antirisonanze delle catene concentrate = zeri invarianti del modello accoppiato a Omega_0 = 0
    con payload bilanciato e senza prodotti d'inerzia (catene disaccoppiate)."""
    from dataclasses import replace

    pl0 = replace(sistema.payload.prodotti_nulli(), d=np.zeros(3))
    s0 = cc.costruisci_sistema_concentrato(par, pl0, (sistema.asse_psi, sistema.asse_theta))
    for th in (-0.5, 0.2):
        ml = lz.linearizza(s0.accoppiato, th, 0.0)
        for nome, f in ca.zeri_numerici(s0, th).items():
            u, j = {"psi_e": ("psi_e", 0), "psi_m": ("psi_m", 0), "theta_e": ("theta_e", 1)}[nome.split("/")[0]]
            z = tf.zeri_invarianti(ml.A, ml.Bu[:, j], ml.Cy[ml.riga(u)])
            fz = np.abs(z[np.imag(z) > 0]) / (2 * np.pi)
            assert np.min(np.abs(fz / f - 1.0)) < 1e-6, nome


@pytest.mark.lento
def test_simulatore_S1_S3(sim: scm.SimulatoreNE, par: pr.Parametri) -> None:
    from gimbal_mdof import simulate as sm

    r = sm.S1(sim, 0.0, t_fine=0.05, n_punti=251)
    for u in ("psi_e", "psi_m", "psi_L"):
        assert r.errori[u] < 1e-2
    r3 = sm.S3(sim, 0.0, par.omega_max, t_fine=0.05, n_punti=101)
    assert r3.derive["energia"] < 1e-8 and r3.derive["momento_psi"] < 1e-8


@pytest.mark.lento
def test_S4_non_linearita(sim: scm.SimulatoreNE, par: pr.Parametri) -> None:
    """Lo scarto non lineare/lineare cresce con l'ampiezza (termini trigonometrici)."""
    th0 = np.deg2rad(-45.0)
    e = [scm.S4(sim, th0, par.omega_max, delta_theta=np.deg2rad(d), t_fine=0.8, n_punti=401).errori["theta_p"]
         for d in (0.5, 5.0)]
    assert e[0] < 1e-2 and e[1] > 3 * e[0]
