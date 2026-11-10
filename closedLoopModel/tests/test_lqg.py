"""Sintesi LQG di Gawronski: forma modale, guadagni e funzione di trasferimento del controllore."""

import numpy as np
import pytest

ct = pytest.importorskip("control")

from lucy_cl import lqg  # noqa: E402

A = np.array([[0, 1, 0, 0], [0, -16.0, 16.0, 0], [0, 0, 0, 1], [0, 0, -(2 * np.pi * 50) ** 2, -2 * 0.02 * 2 * np.pi * 50]])
A[1, 2] = 0.0
A[1, 3] = 0.1
B = np.array([[0.0], [16.0], [0.0], [100.0]])
C = np.array([[1.0, 0.0, 1e-3, 0.0]])


def frf(A_, B_, C_, f):
    return (C_ @ np.linalg.solve(2j * np.pi * f * np.eye(len(A_)) - A_, B_))[0, 0]


def test_forma_modale_conserva_la_frf_e_mette_il_rigido_per_primo():
    Am, Bm, Cm, modi = lqg.forma_modale(A, B, C)
    assert not np.isfinite(modi[0].hinf)
    assert abs(modi[1].f_hz - 2.546) < 0.01 and abs(modi[-1].f_hz - 50.0) < 0.1
    for f in (0.1, 3.0, 50.0):
        assert frf(Am, Bm, Cm, f) == pytest.approx(frf(A, B, C, f), rel=1e-8)
    assert np.allclose(Am[0], 0.0, atol=1e-8)


def test_posizione_primo_stato_conserva_la_frf():
    Ap, Bp, Cp = lqg.posizione_primo_stato(*lqg.forma_modale(A, B, C)[:3])
    assert np.allclose(Cp, np.eye(4)[:1])
    assert frf(Ap, Bp, Cp, 3.0) == pytest.approx(frf(A, B, C, 3.0), rel=1e-8)


def test_integratore_puro_guadagni_del_libro():
    """Su G = 1/s: k_i = sqrt(q_i), k_p^2 = q_p + 2 k_i (a meno della discretizzazione)."""
    c = lqg.sintesi(np.zeros((1, 1)), np.ones((1, 1)), np.ones((1, 1)), 1e-3, 25.0, 90.0, np.zeros(0))
    assert c.k_i == pytest.approx(5.0, rel=0.02)
    assert c.k_p == pytest.approx(10.0, rel=0.02)
    assert np.abs(c.poli_regolatore).max() < 1 and np.abs(c.poli_stimatore).max() < 1


@pytest.mark.parametrize("ritardo", [0, 1])
def test_funzione_trasferimento_uguale_al_passo_nel_tempo(ritardo):
    Ap, Bp, Cp = lqg.posizione_primo_stato(*lqg.forma_modale(A, B, C)[:3])
    c = lqg.sintesi(Ap, Bp, Cp, 1e-3, 100.0, 200.0, np.zeros(3))
    w = 2 * np.pi * 7.0 * 1e-3
    n = 40000
    st, coda, u = c.stato_iniziale(), [0.0] * ritardo, np.empty(n)
    for k in range(n):
        y = np.cos(w * k)
        coda.append(c.uscita(st, 0.0))
        u[k] = coda.pop(0)
        c.stima(st, y, u[k])
    k = np.arange(n - 7000, n)  # 49 periodi esatti
    misurata = 2 * np.mean(u[k] * np.exp(-1j * w * k))
    _, Fy = c.funzione_trasferimento(np.array([np.exp(1j * w)]), ritardo)
    assert -misurata == pytest.approx(Fy[0], rel=1e-3)


def test_con_cuscinetti_ha_la_funzione_trasferimento_attesa():
    """K_b(s) = k_0 + k_1 s/(s + a): stessa risposta di A + K_b(s) Bw v' calcolata in frequenza."""
    A0 = np.array([[0.0, 1.0], [0.0, -2.0]])
    Bu, Bw, C = np.array([0.0, 3.0]), np.array([0.0, -0.5]), np.array([1.0, 0.0])
    vx, vv = np.array([1.0, 0.0]), np.array([0.0, 1.0])
    k0, k1, a = 40.0, 300.0, 50.0
    A2, B2, C2 = lqg.con_cuscinetti(A0, Bu, Bw, C, vx, vv, k0, k1, a)
    s = 2j * np.pi * 7.0
    h = k0 + k1 * s / (s + a)
    atteso = C @ np.linalg.solve(s * np.eye(2) - A0 - h * np.outer(Bw, vx), Bu)
    assert C2 @ np.linalg.solve(s * np.eye(3) - A2, B2) == pytest.approx(atteso, rel=1e-10)


def test_v_continua_da_i_poli_della_fare():
    """Con ``v_continua`` i poli dello stimatore discreto sono quelli della (8.8)."""
    A = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -3.0]])[1:, 1:]
    B, C = np.array([[1.0], [1.0]]), np.array([[1.0, 0.5]])
    q_p, q_f = 0.3, np.array([2.0])
    passo = 1e-3
    pc = ct.lqe(A, np.eye(2), C, np.diag(np.r_[q_p, q_f]), np.eye(1))[2]
    pd = lqg.sintesi(A, B, C, passo, 0.1, q_p, q_f, v_continua=True).poli_stimatore
    assert np.allclose(np.sort(np.abs(np.log(pd.astype(complex)) / passo)), np.sort(np.abs(pc)), rtol=2e-2)
