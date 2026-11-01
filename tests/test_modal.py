"""Troncamento modale alla Gawronski: test di accettazione T1-T9.

I test usano la configurazione congelata ``tests/config_riferimento.yaml`` (method:
guyan, senza blocco ``reduction.modal``: si usano i valori predefiniti, banda 2 kHz,
taglio 4 kHz, nessuna correzione statica). Le tolleranze sono quelle della
richiesta; i numeri del modello reale stanno nel report, non qui.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from gimbal_mdof import coupled as cp
from gimbal_mdof import linearize as lz
from gimbal_mdof import modal as md
from gimbal_mdof import params as pr


@pytest.fixture(scope="module")
def sistema(par: pr.Parametri) -> md.SistemaModale:
    """Sistema modale con ciclo di validazione (criteri T5) su tutti i punti di lavoro."""
    return md.costruisci_sistema_modale(par)


def _con_pm(s: md.SistemaModale, **kw) -> md.SistemaModale:
    return replace(s, pm=replace(s.pm, **kw), extra={"azimut": (), "polare": ()})


@pytest.fixture(scope="module")
def tutti(sistema: md.SistemaModale) -> md.SistemaModale:
    """Tutti i modi del FEM completo (base modale invertibile)."""
    return _con_pm(sistema, n_flessibili={"azimut": sistema.ast_psi.n - 1,
                                          "polare": sistema.ast_theta.n - 1})


def _rel(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-300))


# --- T1: struttura della base modale -------------------------------------------------


@pytest.mark.parametrize("asse", ["azimut", "polare"])
@pytest.mark.parametrize("theta_ref_deg", [0.0, -30.0])
def test_T1_base_modale(sistema: md.SistemaModale, asse: str, theta_ref_deg: float) -> None:
    ast = sistema.asse(asse)
    an = sistema.analisi(asse, np.deg2rad(theta_ref_deg))
    uno = np.ones(ast.n)
    assert np.max(np.abs(ast.K @ uno)) <= 1e-10 * np.max(np.abs(ast.K))
    assert np.array_equal(an.Phi[:, 0], uno)
    Mm = an.Phi.T @ an.M_bar @ an.Phi
    Km = an.Phi.T @ ast.K @ an.Phi
    attesa_M = np.diag([an.J_tot] + [1.0] * (ast.n - 1))
    assert _rel(Mm, attesa_M) < 1e-10
    assert np.max(np.abs(Km - np.diag(an.w2))) <= 1e-9 * np.max(an.w2)
    assert an.J_tot == pytest.approx(uno @ an.M_bar @ uno, rel=1e-14)
    for c in (an.c_e, an.c_m, an.c_agg):
        assert c[0] == 1.0
    # smorzamento: Phi^T beta K Phi = diag(0, 2 zeta omega)
    D = an.Phi.T @ (ast.beta * ast.K) @ an.Phi
    assert np.max(np.abs(D - np.diag(2 * an.zeta * an.omega))) <= 1e-9 * np.max(np.abs(D))
    # nel modello accoppiato: M_tilde diagonale e vettori d'aggancio con prima componente 1
    m = sistema.modello(np.deg2rad(theta_ref_deg))
    b = m.base_psi if asse == "azimut" else m.base_theta
    sl = slice(0, b.r) if asse == "azimut" else slice(m.base_psi.r, m.n)
    Mt = m.M_tilde[sl, sl]
    assert _rel(Mt, np.diag(np.diag(Mt))) < 1e-12
    assert Mt[0, 0] == pytest.approx(an.J_tot, rel=1e-12)
    for i in (ast.i_e, ast.i_m, ast.i_agg):
        assert b.T[i, 0] == 1.0


# --- T7: frequenze dei modi tenuti uguali a quelle del FEM completo --------------------


@pytest.mark.parametrize("asse", ["azimut", "polare"])
def test_T7_frequenze_esatte(sistema: md.SistemaModale, asse: str) -> None:
    import scipy.linalg as sla

    an = sistema.analisi(asse)
    b = sistema.base(asse)
    ast = sistema.asse(asse)
    Mb = md.massa_rif(ast, sistema.payload, an.theta_ref)
    w2_rid = sla.eigvalsh(b.T.T @ ast.K @ b.T, b.T.T @ Mb @ b.T)
    w2_full = an.w2[list(b.modi_fem)]
    assert np.allclose(np.sqrt(np.abs(w2_rid[1:])), np.sqrt(w2_full[1:]), rtol=1e-10, atol=0)
    assert abs(w2_rid[0]) < 1e-9 * w2_full[-1]


# --- controlli incrociati con il modello accoppiato SymPy di coupled.py -----------------


def test_base_guyan_uguale_a_coupled(par: pr.Parametri, sistema: md.SistemaModale) -> None:
    """Con la base di Guyan il modello su base generica coincide con coupled.py."""
    sg = cp.costruisci_sistema(par, sistema.modello_massa, "guyan", payload=sistema.payload,
                               assi=sistema.assi_fem)
    mg = sistema.guyan()
    rng = np.random.default_rng(1)
    for _ in range(5):
        q, dq, ddq = rng.normal(size=(3, mg.n)) * [[0.3], [0.5], [2.0]]
        tau = rng.normal(size=2)
        r_sym = (sg.accoppiato.M(q) @ ddq
                 + np.asarray(sg.accoppiato.funzioni.forze(q, dq, sg.accoppiato.parametri)).reshape(-1)
                 - sg.accoppiato.matrice_ingressi() @ tau)
        assert _rel(mg.residuo(q, dq, ddq, tau), r_sym) < 1e-12
        assert mg.energia(q, dq) == pytest.approx(sg.accoppiato.energia(q, dq), rel=1e-12)
    for th, Om in ((np.deg2rad(-30.0), np.deg2rad(10.0)), (0.2, 0.0)):
        a, b = mg.linearizza(th, Om), lz.linearizza(sg.accoppiato, th, Om)
        for X, Y in ((a.M0, b.M0), (a.C0, b.C0), (a.K0, b.K0)):
            assert _rel(X, Y) < 1e-12
        assert np.allclose(a.punto.tau0, b.punto.tau0, rtol=1e-10, atol=1e-12)


def test_payload_linearizzato(sistema: md.SistemaModale) -> None:
    """Jacobiani del payload (SymPy) contro le forme chiuse di Payload (doc. unificato §7.1)."""
    pl = sistema.payload
    fp = md.funzioni_payload()
    p = md.parametri_payload(pl)
    for th, Om in ((-0.4, 0.2), (0.3, -0.1)):
        G0 = np.asarray(fp.G0(th, Om, p), dtype=float)
        K0 = np.asarray(fp.K0(th, Om, p), dtype=float)
        M0 = np.asarray(fp.M0(th, Om, p), dtype=float)
        gamma = pl.gamma(th, Om)
        assert np.allclose(G0, [[0.0, gamma], [-gamma, 0.0]], rtol=1e-12, atol=1e-14)
        assert np.allclose(K0, [[0.0, 0.0], [0.0, pl.k_c(th, Om) + pl.k_g(th)]], rtol=1e-12, atol=1e-12)
        assert np.allclose(M0, [[pl.J_psi(th), pl.M_12(th)], [pl.M_12(th), pl.I_yy]], rtol=1e-12)


# --- T2: EOM non lineari modali = Phi^T x residuo del modello a tutti i GdL -------------


@pytest.mark.parametrize("theta_ref_deg", [0.0, 25.0])
def test_T2_residuo_non_lineare(tutti: md.SistemaModale, theta_ref_deg: float) -> None:
    mm = tutti.modello(np.deg2rad(theta_ref_deg))
    mf = tutti.completo()
    assert mm.n == mf.n
    Phi = np.zeros((mf.n, mm.n))
    rp = mm.base_psi.r
    Phi[: tutti.ast_psi.n, :rp] = mm.base_psi.T
    Phi[tutti.ast_psi.n:, rp:] = mm.base_theta.T
    rng = np.random.default_rng(7)
    for _ in range(20):
        eta, deta = rng.normal(size=(2, mm.n)) * [[0.4], [0.8]]
        ddeta = rng.normal(size=mm.n) * 3.0
        tau = rng.normal(size=2) * 5.0
        r_mod = mm.residuo(eta, deta, ddeta, tau)
        r_fis = mf.residuo(Phi @ eta, Phi @ deta, Phi @ ddeta, tau)
        scala = max(float(np.max(np.abs(Phi.T @ r_fis))), float(np.max(np.abs(mm.B @ tau))))
        assert np.max(np.abs(r_mod - Phi.T @ r_fis)) <= 1e-10 * scala
        assert mm.energia(eta, deta) == pytest.approx(mf.energia(Phi @ eta, Phi @ deta), rel=1e-10)


# --- T9: simulazioni S1-S3 nel modello modale e confronto con tutti i GdL ----------------


@pytest.mark.lento
def test_T9_S1_modale_contro_lineare(sistema: md.SistemaModale) -> None:
    from gimbal_mdof import simulate as sm

    r = sm.S1(sistema.modello(), np.deg2rad(15.0), t_fine=0.1, n_punti=501)
    for u in ("psi_e", "psi_m", "psi_L", "err_psi"):
        assert r.errori[u] < 1e-2, f"{u}: errore relativo {r.errori[u]:.3e}"


@pytest.mark.lento
def test_T9_S3_conservazione(sistema: md.SistemaModale, par: pr.Parametri) -> None:
    from gimbal_mdof import simulate as sm

    r = sm.S3(sistema.modello(), np.deg2rad(15.0), par.omega_max, t_fine=0.1, n_punti=201)
    assert r.derive["energia"] < 1e-8
    assert r.derive["momento_psi"] < 1e-8


@pytest.mark.lento
def test_T9_S1_modale_contro_completo(sistema: md.SistemaModale) -> None:
    """S1 modale (DOP853, rtol 1e-10) contro il modello a tutti i GdL (rtol 1e-8: i modi
    fino a ~260 kHz rendono l'integratore esplicito limitato dalla stabilita')."""
    from gimbal_mdof import simulate as sm

    th = np.deg2rad(15.0)
    kw = dict(durata=0.008, t_fine=0.03, n_punti=301)
    rm = sm.S1(sistema.modello(), th, **kw)
    rf = sm.S1(sistema.completo(), th, rtol=1e-8, atol=1e-10, **kw)
    for u in ("psi_e", "psi_m", "psi_L"):
        i = lz.USCITE.index(u)
        picco = float(np.max(np.abs(rf.y_nl[:, i])))
        assert np.max(np.abs(rm.y_nl[:, i] - rf.y_nl[:, i])) < 1e-3 * picco, u


# --- T3: modello lineare modale (forma 7.2 e stato a blocchi) = modello fisico ------------


def _frf_rel(G: np.ndarray, Gr: np.ndarray) -> float:
    scala = np.max(np.abs(Gr), axis=(1, 2))[:, None, None]
    return float(np.max(np.abs(G - Gr) / scala))


def _stessi_poli(a: np.ndarray, b: np.ndarray, *, rtol: float, atol: float) -> bool:
    """Ogni polo di ``b`` ha in ``a`` un polo entro ``atol + rtol |p|`` (e viceversa)."""
    if a.size != b.size:
        return False
    return all(np.min(np.abs(a - p)) <= atol + rtol * abs(p) for p in b) and all(
        np.min(np.abs(b - p)) <= atol + rtol * abs(p) for p in a)


def _lontano_dai_poli_non_smorzati(w: np.ndarray, poli: np.ndarray) -> np.ndarray:
    """Esclude +/-5 % attorno ai poli con zeta < 1e-3 (pendolo polare: modo rigido senza
    smorzamento strutturale, k_g ~ 1 N m/rad). Li' P = K - omega^2 M ha voci ~1e7 e un
    quasi annullamento nella direzione rigida: e' la FRF FISICA calcolata in doppia
    precisione a perdere cifre (verificato con mpmath a 40 cifre a 0.164 Hz: numpy sul
    modello fisico 3.3e-8, sul modale 7e-9), non il modello modale."""
    p = poli[(np.imag(poli) > 0) & (np.abs(np.real(poli)) < 1e-3 * np.abs(poli))]
    ok = np.ones(w.size, dtype=bool)
    for wp in np.abs(p):
        ok &= np.abs(w / wp - 1.0) > 0.05
    return ok


def test_T3_frf_tutti_i_modi(tutti: md.SistemaModale, par: pr.Parametri) -> None:
    """FRF modale con tutti i modi = FRF fisica, rtol 1e-8, da 1 Hz a 2 kHz.

    Sotto ~1 Hz la FRF fisica in doppia precisione ha un pavimento d'errore
    eps |K| / (omega^2 J) ~ 1e-9..3e-8 (voci di K ~ 1e7 N m/rad contro omega^2 J ~ 1):
    la dinamica lenta si verifica confrontando direttamente i poli sotto 10 Hz.
    """
    w_tutte = 2 * np.pi * np.logspace(0, np.log10(2000.0), 300)
    uscite = ["psi_e", "theta_e", "psi_m", "theta_m", "psi_L", "theta_p"]
    for th in par.theta0:
        for Om in par.omega0:
            ml = tutti.lineare(th, Om)
            mf = tutti.lineare_completo(th, Om)
            # i poli lenti (pendolo, moti rigidi) si confrontano direttamente; il doppio polo
            # nullo di azimut del modello fisico vale +/-6e-5 j rad/s per la "molla a terra"
            # da arrotondamento di K 1 (~1e-8 N m/rad): atol 1e-4 rad/s
            lenti = 2 * np.pi * 10.0
            assert _stessi_poli(ml.poli[np.abs(ml.poli) < lenti], mf.poli[np.abs(mf.poli) < lenti],
                                rtol=1e-8, atol=1e-4)
            w = w_tutte[_lontano_dai_poli_non_smorzati(w_tutte, mf.poli)]
            Gf = mf.frf(w, uscite)
            assert _frf_rel(ml.frf(w, uscite), Gf) < 1e-8
            # stato a blocchi di Gawronski (2.25)/(2.32)
            sg = md.stato_gawronski(ml)
            righe = [lz.USCITE.index(u) for u in uscite]
            Gg = np.stack([sg.C[righe] @ np.linalg.solve(1j * wk * np.eye(sg.A.shape[0]) - sg.A, sg.B)
                           for wk in w])
            assert _frf_rel(Gg, Gf) < 1e-8
            assert _stessi_poli(np.linalg.eigvals(sg.A), ml.poli, rtol=1e-8, atol=1e-4)


def test_T3_blocchi_di_gawronski(sistema: md.SistemaModale) -> None:
    """Blocchi diagonali (2.25) per i flessibili e (2.32) per i rigidi, a Omega_0 = 0 e
    payload senza prodotti d'inerzia ne' eccentricita' (modello disaccoppiato)."""
    pl0 = replace(sistema.payload.prodotti_nulli(), d=np.zeros(3))
    s0 = md.SistemaModale(par=sistema.par, modello_massa=sistema.modello_massa, payload=pl0,
                          assi_fem=sistema.assi_fem,
                          ast_psi=replace(sistema.ast_psi, C_b=0 * sistema.ast_psi.C_b),
                          ast_theta=replace(sistema.ast_theta, C_b=0 * sistema.ast_theta.C_b),
                          pm=sistema.pm, extra=sistema.extra)
    ml = s0.lineare(0.3, 0.0)
    sg = md.stato_gawronski(ml)
    b = ml.accoppiato
    om = np.concatenate([b.base_psi.omega, b.base_theta.omega])
    ze = np.concatenate([b.base_psi.zeta, b.base_theta.zeta])
    A_att = np.zeros_like(sg.A)
    for i, (w, z) in enumerate(zip(om, ze)):
        A_att[2 * i: 2 * i + 2, 2 * i: 2 * i + 2] = ([[0.0, 1.0], [0.0, 0.0]] if w == 0
                                                     else [[0.0, w], [-w, -2 * z * w]])
    assert np.max(np.abs(sg.A - A_att)) <= 1e-9 * np.max(np.abs(A_att))


def test_T3_jacobiano_numerico(tutti: md.SistemaModale) -> None:
    """La (7.2) e' lo jacobiano del residuo non lineare (differenze centrate)."""
    th, Om = np.deg2rad(-30.0), np.deg2rad(10.0)
    for m in (tutti.modello(th), tutti.completo()):
        M0, C0, K0 = m.matrici_lineari(th, Om)
        pt = m.equilibrio(th, Om)
        q0, dq0 = pt.q0, Om * m.direzione_rigida("azimut")
        z = np.zeros(m.n)
        for X, quale in ((M0, "ddq"), (C0, "dq"), (K0, "q")):
            J = np.zeros_like(X)
            for k in range(m.n):
                h = 1e-6
                e = np.zeros(m.n)
                e[k] = h
                args = {"q": q0, "dq": dq0, "ddq": z}
                ap, am = dict(args), dict(args)
                ap[quale] = args[quale] + e
                am[quale] = args[quale] - e
                J[:, k] = (m.residuo(ap["q"], ap["dq"], ap["ddq"], pt.tau0)
                           - m.residuo(am["q"], am["dq"], am["ddq"], pt.tau0)) / (2 * h)
            assert np.max(np.abs(J - X)) <= 1e-5 * np.max(np.abs(X)), quale


# --- T4: limite rigido = modello a 2 GdL / De Donatis linearizzato ----------------------


@pytest.mark.parametrize("theta0_deg, omega_deg", [(-30.0, 10.0), (15.0, 0.0), (40.0, 10.0)])
def test_T4_limite_rigido(sistema: md.SistemaModale, theta0_deg: float, omega_deg: float) -> None:
    from gimbal import linearize as lz2

    from gimbal_mdof import payload as plm

    s0 = _con_pm(sistema, n_flessibili={"azimut": 0, "polare": 0})
    th0, Om0 = np.deg2rad(theta0_deg), np.deg2rad(omega_deg)
    ml = s0.lineare(th0, Om0)
    assert ml.n == 2
    uno_p, uno_t = np.ones(sistema.ast_psi.n), np.ones(sistema.ast_theta.n)
    J_az = float(uno_p @ sistema.ast_psi.M @ uno_p)  # comprende J^ax
    J_pol = float(uno_t @ sistema.ast_theta.M @ uno_t)
    b_psi = float(uno_p @ sistema.ast_psi.C_b @ uno_p)
    b_th = float(uno_t @ sistema.ast_theta.C_b @ uno_t)
    m2 = plm.modello_2dof_equivalente(sistema.payload, J_az, J_pol)
    m2 = replace(m2, par=replace(m2.par, b_psi=b_psi, b_theta=b_th))
    ml2 = lz2.linearizza(m2, -th0, Om0)  # theta_2dof = -theta
    w = 2 * np.pi * np.logspace(-1, 3, 80)
    S = np.diag([1.0, -1.0])
    G2 = S @ ml2.frf(w) @ S
    assert _frf_rel(ml.frf(w, ["psi_e", "theta_e"]), G2) < 1e-10


# --- T5: modello troncato (default) contro il modello a tutti i GdL in banda --------------


def test_T5_troncato_in_banda(sistema: md.SistemaModale, par: pr.Parametri) -> None:
    """Criteri T5 in ogni punto di lavoro, dopo il ciclo di validazione."""
    assert len(sistema.validazione) == len(par.theta0) * len(par.omega0)
    for th in par.theta0:
        for Om in par.omega0:
            e = sistema.confronta(th, Om)
            assert max(e.err_db) < md.SOGLIA_DB, (th, Om, e.err_db)
            assert max(e.err_deg) < md.SOGLIA_DEG, (th, Om, e.err_deg)
            assert all(err <= md.TOLL_ZERO for *_, err in e.zeri), e.zeri
            assert e.err_incrociato <= md.SOGLIA_INCROCIATI
            assert e.passa


def test_T5_ciclo_di_validazione(sistema: md.SistemaModale) -> None:
    """I modi aggiunti sono nell'ordine di Hankel decrescente e il loro numero e' riportato."""
    assert sistema.n_aggiunti == len(sistema.storia_aggiunte)
    for asse in ("azimut", "polare"):
        an = sistema.analisi(asse)
        for k in sistema.extra[asse]:
            assert an.motivo[k] == "validazione"
            assert an.f[k] > sistema.pm.f_taglio


# --- T8: ordine dei modi e numero di stati ------------------------------------------------


def test_T8_ordine_dei_modi(sistema: md.SistemaModale) -> None:
    sp_ = sistema.spettro()
    assert [(a, k) for _, a, k, _ in sp_[:2]] == [("azimut", 0), ("polare", 0)]
    assert all(k > 0 for _, _, k, _ in sp_[2:])
    flessibili = [f for f, _, k, _ in sp_[2:]]
    assert flessibili == sorted(flessibili)
    m = sistema.modello()
    n_tenuti = sum(1 for *_, t in sp_ if t)
    assert m.n == n_tenuti
    assert m.base_psi.nomi[0] == "eta_psi0" and m.base_theta.nomi[0] == "eta_theta0"
    # stato a blocchi: 2 stati per modo, rigidi per primi in ciascun asse
    ml = sistema.lineare(np.deg2rad(-45.0), 0.0)
    assert md.stato_gawronski(ml).A.shape == (2 * n_tenuti, 2 * n_tenuti)


# --- T6: correzione statica con vettori residui ----------------------------------------------


def test_T6_vettori_residui(sistema: md.SistemaModale, par: pr.Parametri) -> None:
    """Modello predefinito (modi del ciclo di validazione) con i vettori statici residui:
    zeri dei canali coppia -> encoder entro 1e-3 e modi tenuti non spostati."""
    s_rv = replace(sistema, pm=replace(sistema.pm, correzione_statica="residual_vectors"))
    for asse in ("azimut", "polare"):
        b = s_rv.base(asse)
        assert b.n_residui >= 1
        assert b.spostamenti is not None and np.max(np.abs(b.spostamenti)) < 1e-3
        assert np.array_equal(b.T[:, 0], np.ones(b.T.shape[0]))
    # senza vettori residui lo zero "statico" di tau_psi -> psi_e e' spostato dal troncamento
    s_puro = _con_pm(sistema)
    zf = md.zeri_leggeri(s_puro.lineare_completo(0.0, 0.0), "psi_e", 0)
    zp = md.zeri_leggeri(s_puro.lineare(0.0, 0.0), "psi_e", 0)
    zr = md.zeri_leggeri(_con_pm(sistema, correzione_statica="residual_vectors").lineare(0.0, 0.0),
                         "psi_e", 0)
    k = int(np.argmin(np.abs(zf - 10.0)))  # zero carico-molla a ~10 Hz
    assert np.min(np.abs(zr / zf[k] - 1.0)) < 1e-6 < np.min(np.abs(zp / zf[k] - 1.0))
    for th in (np.deg2rad(-45.0), 0.0, np.deg2rad(45.0)):
        ml = s_rv.lineare(th, 0.0)
        mf = s_rv.lineare_completo(th, 0.0)
        for j, u in enumerate(("psi_e", "theta_e")):
            zf = md.zeri_leggeri(mf, u, j)
            zf = zf[(zf >= par.f_min) & (zf <= s_rv.pm.f_banda)]
            zm = md.zeri_leggeri(ml, u, j)
            assert zf.size > 0
            for fz in zf:
                assert np.min(np.abs(zm / fz - 1.0)) < 1e-3, (u, fz)
