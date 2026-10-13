"""Test della sezione 13, punti 13-19 e 21: modello accoppiato e linearizzazione."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
import sympy as sp
from gimbal import linearize as lz2

from gimbal_mdof import coupled as cp
from gimbal_mdof import fem
from gimbal_mdof import linearize as lz
from gimbal_mdof import params as pr
from gimbal_mdof import payload as plm
from gimbal_mdof import tf


@pytest.fixture(scope="module")
def sistema(par: pr.Parametri) -> cp.Sistema:
    return cp.costruisci_sistema(par, "consistent", "guyan")


@pytest.fixture(scope="module")
def sistema_lumped(par: pr.Parametri, sistema: cp.Sistema) -> cp.Sistema:
    return cp.costruisci_sistema(par, "consistent", "guyan_lumped", payload=sistema.payload,
                                 assi=(sistema.asse_psi, sistema.asse_theta))


def _con_payload(s: cp.Sistema, pl: plm.Payload, *, k_scala: float = 1.0,
                 smorzato: bool = True) -> cp.ModelloAccoppiato:
    rp = replace(s.rid_psi, K_r=s.rid_psi.K_r * k_scala,
                 C_b_r=s.rid_psi.C_b_r * (1.0 if smorzato else 0.0))
    rt = replace(s.rid_theta, K_r=s.rid_theta.K_r * k_scala,
                 C_b_r=s.rid_theta.C_b_r * (1.0 if smorzato else 0.0))
    b = (s.beta_psi, s.beta_theta) if smorzato else (0.0, 0.0)
    return cp.costruisci_accoppiato(rp, rt, pl, beta_psi=b[0], beta_theta=b[1])


# --- test 13 ----------------------------------------------------------------


def test_poisson() -> None:
    psi, th, dpsi, dth = sp.symbols("psi theta dpsi dtheta", real=True)
    w = cp.velocita_angolare_poisson(psi, th, dpsi, dth)
    atteso = sp.Matrix([dpsi * sp.sin(th), -dth, dpsi * sp.cos(th)])  # terna P, theta alza
    assert sp.simplify(w - atteso) == sp.zeros(3, 1)


# --- test 14 ----------------------------------------------------------------


def test_equazioni_simboliche_target(sistema: cp.Sistema) -> None:
    mod = sistema.accoppiato.simbolico
    st = mod.struttura
    sost, k = cp.sostituzione_catena_lumped(mod)
    sost.update(cp.PRODOTTI_NULLI)
    eom = mod.eom.subs(sost)
    q, dq, ddq = mod.q, mod.dq, mod.ddq
    th = q[st.j_theta_p]
    Jpsi = cp.I_xx * sp.sin(th) ** 2 + cp.I_zz * sp.cos(th) ** 2
    Jpsi1 = sp.diff(Jpsi, th)
    Vg = cp.m_p * cp.g * (cp.d_x * sp.sin(th) + cp.d_z * sp.cos(th))
    M = mod.M_psi.subs(sost), mod.M_theta.subs(sost)
    J_L = M[0][st.i_psi_L, st.i_psi_L] + cp.J_ax
    atteso_L = ((J_L + Jpsi) * ddq[st.j_psi_L] + Jpsi1 * dq[st.j_theta_p] * dq[st.j_psi_L]
                - k["k_psi_2"] * (q[st.j_psi_m] - q[st.j_psi_L]))
    atteso_p = ((M[1][st.i_theta_p, st.i_theta_p] + cp.I_yy) * ddq[st.j_theta_p]
                - sp.Rational(1, 2) * Jpsi1 * dq[st.j_psi_L] ** 2 + sp.diff(Vg, th)
                - k["k_theta_1"] * (q[st.j_theta_e] - th) + k["k_theta_2"] * (th - q[st.j_theta_m]))
    assert sp.simplify(sp.expand(eom[st.j_psi_L] - atteso_L)) == 0
    assert sp.simplify(sp.expand(eom[st.j_theta_p] - atteso_p)) == 0
    # le altre equazioni restano lineari e identiche alle catene isolate
    Mr = sp.diag(M[0], M[1])
    Kr = sp.diag(mod.K_psi.subs(sost), mod.K_theta.subs(sost))
    for i in range(st.n):
        if i in (st.j_psi_L, st.j_theta_p):
            continue
        lineare = (Mr[i, :] * ddq + Kr[i, :] * q)[0] - mod.Q[i]
        assert sp.simplify(eom[i] - lineare) == 0, f"equazione {i}"


def test_antisimmetria(sistema: cp.Sistema) -> None:
    S = cp.verifica_antisimmetria(sistema.accoppiato.simbolico)
    assert S.shape == (7, 7)


def test_target_linearizzazione(sistema: cp.Sistema) -> None:
    """Sezione 10.2: rispetto alle catene isolate cambiano solo gamma, k_c, k_g."""
    pl0 = sistema.payload.prodotti_nulli()
    acc = _con_payload(sistema, pl0)
    st = acc.struttura
    th0, Om0 = np.deg2rad(-30.0), np.deg2rad(10.0)
    ml = lz.linearizza(acc, th0, Om0)
    rp, rt = sistema.rid_psi, sistema.rid_theta
    M_bl = np.block([[rp.M_r, np.zeros((3, 4))], [np.zeros((4, 3)), rt.M_r]])
    C_bl = np.block([[sistema.beta_psi * rp.K_r + rp.C_b_r, np.zeros((3, 4))],
                     [np.zeros((4, 3)), sistema.beta_theta * rt.K_r + rt.C_b_r]])
    K_bl = np.block([[rp.K_r, np.zeros((3, 4))], [np.zeros((4, 3)), rt.K_r]])
    L, P = st.j_psi_L, st.j_theta_p
    M_bl[L, L] += pl0.J_ax + pl0.J_psi(th0)
    M_bl[P, P] += pl0.I_yy
    gamma = (pl0.I_xx - pl0.I_zz) * Om0 * np.sin(2 * th0)
    C_bl[L, P] += gamma
    C_bl[P, L] -= gamma
    K_bl[P, P] += -(pl0.I_xx - pl0.I_zz) * Om0**2 * np.cos(2 * th0) + pl0.k_g(th0)
    for A, B in ((ml.M0, M_bl), (ml.C0, C_bl), (ml.K0, K_bl)):
        assert np.allclose(A, B, rtol=1e-12, atol=1e-12 * np.max(np.abs(B)))


# --- test 15 ----------------------------------------------------------------


@pytest.mark.parametrize("theta0_deg", [-30.0, 15.0])
def test_disaccoppiamento(sistema: cp.Sistema, theta0_deg: float) -> None:
    pl0 = sistema.payload.prodotti_nulli()
    acc = _con_payload(sistema, pl0)
    th0 = np.deg2rad(theta0_deg)
    ml = lz.linearizza(acc, th0, 0.0)
    w = 2 * np.pi * np.logspace(-1, 4, 300)
    G = ml.frf(w, ["psi_e", "theta_e"])
    scala = np.max(np.abs(G))
    assert np.max(np.abs(G[:, 0, 1])) < 1e-12 * scala
    assert np.max(np.abs(G[:, 1, 0])) < 1e-12 * scala
    s0 = replace(sistema, payload=pl0)
    for k, (asse, rid) in enumerate((("azimut", sistema.rid_psi), ("polare", sistema.rid_theta))):
        M, K, C = s0.isolato_ridotto(asse, th0)
        e_in = np.zeros((rid.n, 1))
        e_in[rid.i_motore] = 1.0
        c = np.zeros((1, rid.n))
        c[0, rid.i_encoder] = 1.0
        G_iso = tf.frf_secondo_ordine(M, C, K, e_in, c, w)[:, 0, 0]
        assert np.allclose(G[:, k, k], G_iso, rtol=1e-9, atol=1e-12 * np.max(np.abs(G_iso)))


# --- test 16 ----------------------------------------------------------------


@pytest.mark.parametrize("theta0_deg, omega_deg", [(-30.0, 10.0), (15.0, 0.0)])
def test_limite_rigido(sistema: cp.Sistema, theta0_deg: float, omega_deg: float) -> None:
    """Rigidezze x 1e6: sotto la prima risonanza si ritrova il modello a 2 GdL.

    Lo scarto fisico scala come omega^2 J / k e si annulla; quello numerico
    dell'inversione di P(j omega) cresce come k / (omega^2 J) (stessa perdita di
    cifre documentata in ``gimbal.flex4dof.verifica_limite_rigido``). Con il
    fattore 1e6 della specifica entrambi sono < 1e-5 nella banda 1-10 Hz, che e'
    dunque quella del confronto; sotto 0.1 Hz domina l'errore di arrotondamento.
    """
    pl = sistema.payload
    acc = _con_payload(sistema, pl, k_scala=1e6, smorzato=False)
    th0, Om0 = np.deg2rad(theta0_deg), np.deg2rad(omega_deg)
    ml = lz.linearizza(acc, th0, Om0)
    J_az = float(sistema.rid_psi.M_r.sum()) + pl.J_ax
    J_pol = float(sistema.rid_theta.M_r.sum())
    m2 = plm.modello_2dof_equivalente(pl, J_az, J_pol)
    m2 = replace(m2, par=replace(m2.par, b_psi=0.0, b_theta=0.0))
    # il 2 GdL usa theta che abbassa il puntamento: theta_2dof = -theta, e le FdT
    # incrociate cambiano segno (S G S con S = diag(1, -1))
    ml2 = lz2.linearizza(m2, -th0, Om0)
    w = 2 * np.pi * np.logspace(0, 1, 40)
    G = ml.frf(w, ["psi_e", "theta_e"])
    S = np.diag([1.0, -1.0])
    G2 = S @ ml2.frf(w) @ S
    scala = np.max(np.abs(G2), axis=(1, 2))
    assert np.all(np.abs(G - G2) <= 1e-5 * scala[:, None, None])
    # Poli a bassa frequenza: quelli del modello rigido. Con rigidezze x 1e6 lo
    # spettro copre ~14 decadi in omega^2 e la doppia precisione impone un
    # pavimento (lo stesso con eig(A) e con il pencil generalizzato senza M^-1):
    # ~1e-5 relativo sui poli non nulli e |lambda| ~ 1e-9 lambda_max sul doppio
    # polo nullo (autovalore difettivo). Con fattore 1e2 lo scarto e' ~1e-8.
    p = ml.poli
    lam_max = float(np.max(np.abs(p)))
    p_bassi = p[np.abs(p) < 1e2]
    p2 = ml2.poli
    assert len(p_bassi) == len(p2) == 4
    nulli = np.abs(p2) < 1e-9
    assert np.sum(np.abs(p_bassi) < 1e-8 * lam_max) == np.sum(nulli)
    for z in p2[~nulli]:
        assert np.min(np.abs(p_bassi - z)) < 1e-4 * abs(z)


# --- test 17 ----------------------------------------------------------------


@pytest.mark.parametrize("asse", ["azimut", "polare"])
def test_collocazione(sistema: cp.Sistema, asse: str) -> None:
    th0 = np.deg2rad(-45.0)  # k_g > 0 (theta < 0 abbassa il puntamento): polare stabile
    M, K, _ = sistema.isolato_completo(asse, th0)
    mod = sistema.asse_psi if asse == "azimut" else sistema.asse_theta
    n = mod.n
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-np.linalg.solve(M, K), np.zeros((n, n))]])
    b = np.concatenate([np.zeros(n), np.linalg.solve(M, mod.vettore_forze())])
    c = np.concatenate([mod.vettore_forze(), np.zeros(n)])
    zeri = tf.zeri_invarianti(A, b, c)
    wz = np.sort(np.abs(np.imag(zeri)))[::2]
    wp = np.sort(np.abs(2 * np.pi * fem.frequenze_proprie(M, K)))
    # gli zeri collocati sono le frequenze con il nodo motore bloccato
    liberi = [i for i in range(n) if i != mod.i_motore]
    w_bloc = np.sort(2 * np.pi * fem.frequenze_proprie(M[np.ix_(liberi, liberi)], K[np.ix_(liberi, liberi)]))
    assert len(wz) == len(w_bloc) == len(wp) - 1
    # il pencil di Rosenbrock non e' bilanciato: precisione relativa ~1e-7
    assert np.allclose(wz, w_bloc, rtol=1e-6)
    # interlacciamento in senso debole (teorema di Cauchy). Un modo locale con il
    # nodo motore fermo (es. elementi polari 6 e 7 uguali e simmetrici attorno al
    # nodo 7) da' una cancellazione polo-zero esatta.
    assert np.all(wp[:-1] <= w_bloc * (1 + 1e-10)) and np.all(w_bloc <= wp[1:] * (1 + 1e-10))


# --- test 18 ----------------------------------------------------------------


def test_coerenza_rappresentazioni(sistema: cp.Sistema) -> None:
    ml = lz.linearizza(sistema.accoppiato, np.deg2rad(15.0), np.deg2rad(10.0))
    w = 2 * np.pi * np.logspace(-1, 4, 60)
    G1 = ml.frf(w)
    G2 = ml.frf_spazio_stati(w)
    assert np.allclose(G1, G2, rtol=1e-8, atol=1e-10 * np.max(np.abs(G1)))
    for lam in ml.poli:
        sv = np.linalg.svd(ml.P(lam), compute_uv=False)
        assert sv[-1] < 1e-8 * sv[0]


# --- test 19 ----------------------------------------------------------------


def test_reazione_polare_senza_lavoro_su_psi_L(sistema: cp.Sistema) -> None:
    mod = sistema.accoppiato.simbolico
    st = mod.struttura
    assert mod.Q[st.j_psi_L] == 0
    assert sp.diff(mod.eom[st.j_psi_L], cp.tau_theta) == 0
    assert mod.Q[st.j_theta_m] == cp.tau_theta and mod.Q[st.j_psi_m] == cp.tau_psi


# --- test 21 ----------------------------------------------------------------


@pytest.mark.parametrize("theta0_deg, omega_deg", [(-45.0, 0.0), (0.0, 10.0), (45.0, 10.0)])
def test_equilibrio_statico_polare(sistema: cp.Sistema, theta0_deg: float, omega_deg: float) -> None:
    ml = lz.linearizza(sistema.accoppiato, np.deg2rad(theta0_deg), np.deg2rad(omega_deg))
    th_e, th_p, th_m, th_r = ml.punto.torsioni_polari(3)
    k2 = sistema.rid_theta.k_catena[1]
    tau_th = ml.punto.tau0[1]
    assert th_e == pytest.approx(th_p, rel=1e-12, abs=1e-14)
    assert th_r == pytest.approx(th_m, rel=1e-12, abs=1e-14)
    assert th_m - th_p == pytest.approx(tau_th / k2, rel=1e-10)
    # senza attrito ai cuscinetti l'azimut ruota rigidamente
    assert np.allclose(ml.punto.q0[:3], 0.0, atol=1e-14)
    # coppia di mantenimento: gravita' + termine centrifugo
    pl = sistema.payload
    th0, Om0 = np.deg2rad(theta0_deg), np.deg2rad(omega_deg)
    atteso = pl.dV_dtheta(th0) - 0.5 * pl.J_psi_prime(th0) * Om0**2
    assert tau_th == pytest.approx(atteso, rel=1e-10)
