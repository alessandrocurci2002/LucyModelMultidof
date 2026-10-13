"""Modello accoppiato dei due assi ridotti (sezione 9).

Coordinate :math:`\\mathbf q = [\\mathbf q_\\psi,\\ \\mathbf q_\\theta]`, con
:math:`\\mathbf q_\\psi = [\\psi_e, \\psi_m, \\psi_L]` e
:math:`\\mathbf q_\\theta = [\\theta_e, \\theta_p, \\theta_m, \\theta_r]` nella
topologia di riferimento. La forcella ruota con :math:`\\psi_L` (Q12); il
telescopio con :math:`\\theta_p`, rotazione del nodo P relativa alla forcella.

.. math::
    T = \\tfrac12\\dot{\\mathbf q}_\\psi^{\\mathsf T}\\mathbf M_{r,\\psi}\\dot{\\mathbf q}_\\psi
      + \\tfrac12 J^{ax}\\dot\\psi_L^2
      + \\tfrac12\\dot{\\mathbf q}_\\theta^{\\mathsf T}\\mathbf M_{r,\\theta}\\dot{\\mathbf q}_\\theta
      + \\tfrac12\\boldsymbol\\omega_p|_B^{\\mathsf T}\\mathbf I_{p,O}|_B\\boldsymbol\\omega_p|_B

.. math::
    V = \\tfrac12\\mathbf q_\\psi^{\\mathsf T}\\mathbf K_{r,\\psi}\\mathbf q_\\psi
      + \\tfrac12\\mathbf q_\\theta^{\\mathsf T}\\mathbf K_{r,\\theta}\\mathbf q_\\theta
      + m_p g(-d_x\\sin\\theta_p + d_z\\cos\\theta_p)

Le equazioni di Lagrange **non sono scritte a mano**: le deriva SymPy. Le
matrici ridotte entrano come simboli con la struttura tridiagonale di una
catena, cosi' la derivazione vale per qualunque geometria; i valori numerici si
sostituiscono solo alla fine.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable, Final

import numpy as np
import sympy as sp
from gimbal import kinematics as kin
from gimbal import lagrange as lg

from .payload import Payload
from .reduction import ModelloRidotto

# ---------------------------------------------------------------------------
# Simboli comuni
# ---------------------------------------------------------------------------

I_xx, I_yy, I_zz = lg.I_xx, lg.I_yy, lg.I_zz
I_xy, I_xz, I_yz = lg.I_xy, lg.I_xz, lg.I_yz
g, d_x, d_y, d_z = lg.g, lg.d_x, lg.d_y, lg.d_z
tau_psi, tau_theta = lg.tau_psi, lg.tau_theta
m_p = sp.Symbol("m_p", positive=True)
J_ax = sp.Symbol("J^{ax}", nonnegative=True)

PRODOTTI_NULLI: Final[dict[sp.Symbol, int]] = dict(lg.PRODOTTI_NULLI)


@dataclass(frozen=True)
class Struttura:
    """Struttura delle due catene ridotte: chiave della derivazione simbolica."""

    ruoli_psi: tuple[str, ...]
    ruoli_theta: tuple[str, ...]
    i_psi_e: int
    i_psi_m: int
    i_psi_L: int
    i_theta_e: int
    i_theta_p: int
    i_theta_m: int

    @property
    def n_psi(self) -> int:
        return len(self.ruoli_psi)

    @property
    def n_theta(self) -> int:
        return len(self.ruoli_theta)

    @property
    def n(self) -> int:
        return self.n_psi + self.n_theta

    # indici globali in q
    @property
    def j_psi_e(self) -> int:
        return self.i_psi_e

    @property
    def j_psi_m(self) -> int:
        return self.i_psi_m

    @property
    def j_psi_L(self) -> int:
        return self.i_psi_L

    @property
    def j_theta_e(self) -> int:
        return self.n_psi + self.i_theta_e

    @property
    def j_theta_p(self) -> int:
        return self.n_psi + self.i_theta_p

    @property
    def j_theta_m(self) -> int:
        return self.n_psi + self.i_theta_m

    @property
    def nomi(self) -> tuple[str, ...]:
        return tuple(f"psi_{r}" for r in self.ruoli_psi) + tuple(
            f"theta_{r}" for r in self.ruoli_theta)


def struttura_da(rid_psi: ModelloRidotto, rid_theta: ModelloRidotto) -> Struttura:
    assert rid_theta.i_P is not None
    return Struttura(
        ruoli_psi=rid_psi.ruoli, ruoli_theta=rid_theta.ruoli,
        i_psi_e=rid_psi.i_encoder, i_psi_m=rid_psi.i_motore, i_psi_L=rid_psi.i_estremo,
        i_theta_e=rid_theta.i_encoder, i_theta_p=rid_theta.i_P, i_theta_m=rid_theta.i_motore,
    )


def _matrice_tridiagonale(nome: str, apice: str, ruoli: tuple[str, ...]) -> tuple[sp.Matrix, list[sp.Symbol]]:
    """Matrice simmetrica tridiagonale di simboli ``nome^{apice}_{ij}``."""
    n = len(ruoli)
    A = sp.zeros(n, n)
    simboli: list[sp.Symbol] = []
    for i in range(n):
        for j in (i, i + 1):
            if j >= n:
                continue
            s = sp.Symbol(f"{nome}^{{{apice}}}_{{{ruoli[i]}{ruoli[j]}}}", real=True)
            A[i, j] = A[j, i] = s
            simboli.append(s)
    return A, simboli


def derivata_temporale(espr: sp.Expr, q: sp.Matrix, dq: sp.Matrix, ddq: sp.Matrix) -> sp.Expr:
    """:math:`\\frac{d}{dt}f(\\mathbf q, \\dot{\\mathbf q}) = \\sum_i
    \\frac{\\partial f}{\\partial q_i}\\dot q_i + \\frac{\\partial f}{\\partial\\dot q_i}\\ddot q_i`."""
    return sum((sp.diff(espr, q[i]) * dq[i] + sp.diff(espr, dq[i]) * ddq[i]
                for i in range(len(q))), sp.S.Zero)


# ---------------------------------------------------------------------------
# Cinematica del payload (sezione 9.1)
# ---------------------------------------------------------------------------


def velocita_angolare_poisson(psi_L: sp.Symbol, theta_p: sp.Symbol, dpsi_L: sp.Symbol,
                              dtheta_p: sp.Symbol) -> sp.Matrix:
    """:math:`\\boldsymbol\\omega_p|_B` dalla formula di Poisson, :math:`[\\omega]_\\times = R^{\\mathsf T}\\dot R`.

    ``R = R_z(psi_L) R_y(-theta_p)`` nella terna interna P (x_P puntamento, y_P = j_hat,
    z_P verso l'alto a theta = 0); theta > 0 alza il puntamento. La forcella ruota con :math:`\\psi_L`, il
    telescopio con :math:`\\theta_p` relativo alla forcella.
    """
    t = sp.Symbol("t", real=True)
    a = sp.Function("a", real=True)(t)
    b = sp.Function("b", real=True)(t)
    R = kin.R_z(a) * kin.R_y(-b)
    S = sp.simplify(R.T * sp.diff(R, t))
    w = sp.simplify(kin.vee(S))
    sost = {sp.Derivative(a, t): dpsi_L, sp.Derivative(b, t): dtheta_p}
    w = w.subs(sost).subs({a: psi_L, b: theta_p})
    return sp.simplify(w)


def velocita_angolare_forma_chiusa(theta_p: sp.Symbol, dpsi_L: sp.Symbol,
                                   dtheta_p: sp.Symbol) -> sp.Matrix:
    """Forma chiusa nella terna P: :math:`(\\dot\\psi_L\\sin\\theta_p, -\\dot\\theta_p, \\dot\\psi_L\\cos\\theta_p)`."""
    return sp.Matrix([dpsi_L * sp.sin(theta_p), -dtheta_p, dpsi_L * sp.cos(theta_p)])


# ---------------------------------------------------------------------------
# Modello simbolico
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModelloAccoppiatoSimbolico:
    """Equazioni di moto simboliche del sistema accoppiato ridotto.

    ``eom = 0`` e' il residuo delle equazioni di Lagrange; ``residuo`` e' la
    stessa quantita' senza le forze generalizzate, nella forma
    :math:`\\mathbf M\\ddot{\\mathbf q} + \\mathbf C\\dot{\\mathbf q}
    + \\mathbf C_d\\dot{\\mathbf q} + \\partial V/\\partial\\mathbf q`.
    """

    struttura: Struttura
    q: sp.Matrix
    dq: sp.Matrix
    ddq: sp.Matrix
    M_psi: sp.Matrix
    K_psi: sp.Matrix
    C_psi: sp.Matrix
    M_theta: sp.Matrix
    K_theta: sp.Matrix
    C_theta: sp.Matrix
    omega_p: sp.Matrix
    T: sp.Expr
    V: sp.Expr
    D: sp.Expr
    L: sp.Expr
    Q: sp.Matrix
    lavoro_virtuale: sp.Expr
    eom: sp.Matrix
    M: sp.Matrix
    C: sp.Matrix
    Cd: sp.Matrix
    gradV: sp.Matrix
    residuo: sp.Matrix
    parametri: tuple[sp.Symbol, ...]

    @property
    def tensore_payload(self) -> sp.Matrix:
        return lg.tensore_B()


def _christoffel(M: sp.Matrix, q: sp.Matrix, dq: sp.Matrix) -> sp.Matrix:
    """Matrice di Coriolis dai simboli di Christoffel di prima specie."""
    n = len(q)
    dM = [M.diff(q[i]) for i in range(n)]
    attivi = [i for i in range(n) if not dM[i].is_zero_matrix]
    C = sp.zeros(n, n)
    for k in range(n):
        for j in range(n):
            somma = sp.S.Zero
            for i in range(n):
                c = sp.S.Zero
                if i in attivi:
                    c += dM[i][k, j]
                if j in attivi:
                    c += dM[j][k, i]
                if k in attivi:
                    c -= dM[k][i, j]
                if c != 0:
                    somma += sp.Rational(1, 2) * c * dq[i]
            C[k, j] = sp.simplify(somma)
    return C


@lru_cache(maxsize=8)
def costruisci_simbolico(st: Struttura) -> ModelloAccoppiatoSimbolico:
    """Deriva con SymPy le equazioni di moto del sistema accoppiato (sezione 9.5)."""
    rp, rt = st.ruoli_psi, st.ruoli_theta
    q_psi = sp.Matrix([sp.Symbol(f"\\psi_{{{r}}}", real=True) for r in rp])
    q_th = sp.Matrix([sp.Symbol(f"\\theta_{{{r}}}", real=True) for r in rt])
    dq_psi = sp.Matrix([sp.Symbol(f"\\dot{{\\psi}}_{{{r}}}", real=True) for r in rp])
    dq_th = sp.Matrix([sp.Symbol(f"\\dot{{\\theta}}_{{{r}}}", real=True) for r in rt])
    ddq_psi = sp.Matrix([sp.Symbol(f"\\ddot{{\\psi}}_{{{r}}}", real=True) for r in rp])
    ddq_th = sp.Matrix([sp.Symbol(f"\\ddot{{\\theta}}_{{{r}}}", real=True) for r in rt])
    q = sp.Matrix.vstack(q_psi, q_th)
    dq = sp.Matrix.vstack(dq_psi, dq_th)
    ddq = sp.Matrix.vstack(ddq_psi, ddq_th)

    M_psi, sM_psi = _matrice_tridiagonale("M", "\\psi", rp)
    K_psi, sK_psi = _matrice_tridiagonale("K", "\\psi", rp)
    C_psi, sC_psi = _matrice_tridiagonale("C", "\\psi", rp)
    M_th, sM_th = _matrice_tridiagonale("M", "\\theta", rt)
    K_th, sK_th = _matrice_tridiagonale("K", "\\theta", rt)
    C_th, sC_th = _matrice_tridiagonale("C", "\\theta", rt)

    psi_L, dpsi_L = q[st.j_psi_L], dq[st.j_psi_L]
    theta_p, dtheta_p = q[st.j_theta_p], dq[st.j_theta_p]

    # --- cinematica del payload (Poisson, verificata con la forma chiusa) -----
    omega = velocita_angolare_poisson(psi_L, theta_p, dpsi_L, dtheta_p)
    scarto = sp.simplify(omega - velocita_angolare_forma_chiusa(theta_p, dpsi_L, dtheta_p))
    if scarto != sp.zeros(3, 1):
        raise AssertionError(f"omega_p da Poisson diversa dalla forma chiusa: {scarto.T}")

    # --- energie ---------------------------------------------------------------
    half = sp.Rational(1, 2)
    T_payload = sp.expand(sp.trigsimp(sp.expand(half * (omega.T * lg.tensore_B() * omega)[0, 0])))
    T = (half * (dq_psi.T * M_psi * dq_psi)[0, 0] + half * J_ax * dpsi_L**2
         + half * (dq_th.T * M_th * dq_th)[0, 0] + T_payload)
    Z_in_B = kin.R_y(-theta_p).T * sp.Matrix([0, 0, 1])
    V_g = sp.simplify(m_p * g * (sp.Matrix([d_x, d_y, d_z]).T * Z_in_B)[0, 0])
    V = half * (q_psi.T * K_psi * q_psi)[0, 0] + half * (q_th.T * K_th * q_th)[0, 0] + V_g
    D = half * (dq_psi.T * C_psi * dq_psi)[0, 0] + half * (dq_th.T * C_th * dq_th)[0, 0]
    L = T - V

    # --- forze generalizzate dal lavoro virtuale (sezione 9.4) -----------------
    # Rotazioni virtuali: rotore azimut dpsi_m k; statore azimut 0 (base fissa);
    # rotore polare dpsi_L k + dtheta_m J; statore polare (forcella) dpsi_L k.
    dlt = sp.Matrix([sp.Symbol(f"\\delta {s.name}", real=True) for s in q])
    # Versori in terna 1: k_hat = (0,0,1), asse polare J_hat = -j_hat = (0,-1,0).
    Zv, Jh = sp.Matrix([0, 0, 1]), sp.Matrix([0, -1, 0])
    phi_rot_psi = dlt[st.j_psi_m] * Zv
    phi_rot_th = dlt[st.j_psi_L] * Zv + dlt[st.j_theta_m] * Jh
    phi_stat_th = dlt[st.j_psi_L] * Zv
    dW = (tau_psi * Zv.T * phi_rot_psi)[0, 0] + (tau_theta * Jh.T * (phi_rot_th - phi_stat_th))[0, 0]
    dW = sp.expand(dW)
    Q = sp.Matrix([sp.diff(dW, dlt[i]) for i in range(st.n)])

    # --- equazioni di Lagrange ------------------------------------------------
    eom = sp.zeros(st.n, 1)
    for i in range(st.n):
        dL_ddq = sp.diff(L, dq[i])
        eom[i] = sp.expand(derivata_temporale(dL_ddq, q, dq, ddq) - sp.diff(L, q[i])
                           + sp.diff(D, dq[i]) - Q[i])

    # --- forma matriciale ------------------------------------------------------
    M = sp.hessian(T, list(dq)).applyfunc(sp.simplify)
    C = _christoffel(M, q, dq)
    Cd = sp.hessian(D, list(dq))
    gradV = sp.Matrix([sp.simplify(sp.diff(V, q[i])) for i in range(st.n)])
    residuo = M * ddq + C * dq + Cd * dq + gradV

    scarto = (eom - (residuo - Q)).applyfunc(lambda e: sp.simplify(sp.expand(e)))
    if scarto != sp.zeros(st.n, 1):
        raise AssertionError(f"forma matriciale e Lagrange non coincidono: {scarto.T}")

    parametri = tuple(sM_psi + sK_psi + sC_psi + sM_th + sK_th + sC_th) + (
        J_ax, I_xx, I_yy, I_zz, I_xy, I_xz, I_yz, m_p, g, d_x, d_y, d_z)

    return ModelloAccoppiatoSimbolico(
        struttura=st, q=q, dq=dq, ddq=ddq, M_psi=M_psi, K_psi=K_psi, C_psi=C_psi,
        M_theta=M_th, K_theta=K_th, C_theta=C_th, omega_p=omega, T=T, V=V, D=D, L=L,
        Q=Q, lavoro_virtuale=dW, eom=eom, M=M, C=C, Cd=Cd, gradV=gradV, residuo=residuo,
        parametri=parametri,
    )


# ---------------------------------------------------------------------------
# Verifiche strutturali
# ---------------------------------------------------------------------------


def verifica_antisimmetria(mod: ModelloAccoppiatoSimbolico) -> sp.Matrix:
    """:math:`\\dot{\\mathbf M} - 2\\mathbf C` antisimmetrica (passivita')."""
    Mdot = mod.M.applyfunc(lambda e: derivata_temporale(e, mod.q, mod.dq, mod.ddq))
    S = (Mdot - 2 * mod.C).applyfunc(sp.simplify)
    if (S + S.T).applyfunc(sp.simplify) != sp.zeros(*S.shape):
        raise AssertionError("dM/dt - 2C non e' antisimmetrica")
    return S


def sostituzione_catena_lumped(mod: ModelloAccoppiatoSimbolico) -> tuple[dict, dict[str, sp.Symbol]]:
    """Sostituzione che porta le matrici ridotte nella forma di catena pura.

    Massa diagonale (``guyan_lumped``), rigidezze :math:`K_{ii} = k_{i-1} + k_i`,
    :math:`K_{i,i+1} = -k_i`, smorzamento nullo. Restituisce anche i simboli
    :math:`k_{\\psi,i}`, :math:`k_{\\theta,i}`.
    """
    sost: dict = {}
    k: dict[str, sp.Symbol] = {}
    for nome, Mm, Km, Cm, n in (("psi", mod.M_psi, mod.K_psi, mod.C_psi, mod.struttura.n_psi),
                                ("theta", mod.M_theta, mod.K_theta, mod.C_theta, mod.struttura.n_theta)):
        ks = [sp.Symbol(f"k_{{\\{nome},{i + 1}}}", positive=True) for i in range(n - 1)]
        for i, s in enumerate(ks):
            k[f"k_{nome}_{i + 1}"] = s
        for i in range(n):
            for j in range(i, min(i + 2, n)):
                if j == i + 1:
                    sost[Mm[i, j]] = 0
                    sost[Km[i, j]] = -ks[i]
                else:
                    sost[Km[i, i]] = (ks[i - 1] if i > 0 else 0) + (ks[i] if i < n - 1 else 0)
                sost[Cm[i, j]] = 0
    return sost, k


# ---------------------------------------------------------------------------
# Modello numerico
# ---------------------------------------------------------------------------


def _verifica_tridiagonale(A: np.ndarray, nome: str) -> None:
    scala = max(float(np.max(np.abs(A))), 1e-300)
    if np.max(np.abs(np.triu(A, 2))) > 1e-10 * scala:
        raise AssertionError(f"{nome} non e' tridiagonale: i master non formano una catena")


def valori_parametri(mod: ModelloAccoppiatoSimbolico, rid_psi: ModelloRidotto,
                     rid_theta: ModelloRidotto, pl: Payload, C_psi: np.ndarray,
                     C_theta: np.ndarray) -> np.ndarray:
    """Vettore numerico allineato a ``mod.parametri``."""
    valori: dict[sp.Symbol, float] = {}
    for S, A, nome in ((mod.M_psi, rid_psi.M_r, "M_psi"), (mod.K_psi, rid_psi.K_r, "K_psi"),
                       (mod.C_psi, C_psi, "C_psi"), (mod.M_theta, rid_theta.M_r, "M_theta"),
                       (mod.K_theta, rid_theta.K_r, "K_theta"), (mod.C_theta, C_theta, "C_theta")):
        _verifica_tridiagonale(A, nome)
        n = A.shape[0]
        for i in range(n):
            for j in range(i, min(i + 2, n)):
                valori[S[i, j]] = float(A[i, j])
    valori.update({J_ax: pl.J_ax, I_xx: pl.I_xx, I_yy: pl.I_yy, I_zz: pl.I_zz, I_xy: pl.I_xy,
                   I_xz: pl.I_xz, I_yz: pl.I_yz, m_p: pl.m, g: pl.g, d_x: pl.d_x, d_y: pl.d_y,
                   d_z: pl.d_z})
    return np.array([valori[s] for s in mod.parametri], dtype=float)


@dataclass(frozen=True)
class FunzioniNumeriche:
    """Funzioni ``lambdify`` del modello, con i parametri come argomento."""

    M: Callable
    forze: Callable
    T: Callable
    V: Callable
    p_psi: Callable
    M0: Callable
    C0: Callable
    K0: Callable
    residuo_statico: Callable


@lru_cache(maxsize=8)
def funzioni_numeriche(st: Struttura) -> FunzioniNumeriche:
    """Compila le funzioni numeriche (una volta per struttura)."""
    mod = costruisci_simbolico(st)
    q, dq, ddq = list(mod.q), list(mod.dq), list(mod.ddq)
    p = list(mod.parametri)
    th0, Om0 = sp.symbols("theta_0 Omega_0", real=True)

    def lam(args: list, espr: sp.Basic) -> Callable:
        return sp.lambdify(args, espr, modules="numpy", cse=True, dummify=True)

    forze = mod.C * mod.dq + mod.Cd * mod.dq + mod.gradV
    p_psi = sum((sp.diff(mod.L, mod.dq[i]) for i in range(st.n_psi)), sp.S.Zero)

    # --- linearizzazione (sezione 10.2): jacobiani del residuo nel punto di lavoro
    punto = {mod.q[st.j_theta_p]: th0}
    for i in range(st.n):
        punto[mod.dq[i]] = Om0 if i < st.n_psi else 0
        punto[mod.ddq[i]] = 0
    r = mod.residuo
    M0 = r.jacobian(mod.ddq).subs(punto)
    C0 = r.jacobian(mod.dq).subs(punto)
    K0 = r.jacobian(mod.q).subs(punto)
    altre = [s for s in q if s != mod.q[st.j_theta_p]]
    for nome, A in (("M_0", M0), ("C_0", C0), ("K_0", K0)):
        if A.free_symbols & set(altre):
            raise AssertionError(f"{nome} dipende da coordinate diverse da theta_p")
    # residuo nel punto di lavoro come funzione delle coordinate (statica, 10.1)
    r_stat = r.subs({k: v for k, v in punto.items() if k != mod.q[st.j_theta_p]})

    return FunzioniNumeriche(
        M=lam([q, p], mod.M),
        forze=lam([q, dq, p], forze),
        T=lam([q, dq, p], mod.T),
        V=lam([q, p], mod.V),
        p_psi=lam([q, dq, p], p_psi),
        M0=lam([th0, Om0, p], M0),
        C0=lam([th0, Om0, p], C0),
        K0=lam([th0, Om0, p], K0),
        residuo_statico=lam([q, Om0, p], r_stat),
    )


@dataclass(frozen=True)
class ModelloAccoppiato:
    """Modello accoppiato numerico: parametri assegnati, pronto per linearizzare
    e integrare."""

    struttura: Struttura
    parametri: np.ndarray
    rid_psi: ModelloRidotto
    rid_theta: ModelloRidotto
    payload: Payload
    beta_psi: float
    beta_theta: float
    funzioni: FunzioniNumeriche

    @property
    def n(self) -> int:
        return self.struttura.n

    @property
    def simbolico(self) -> ModelloAccoppiatoSimbolico:
        return costruisci_simbolico(self.struttura)

    def matrice_ingressi(self) -> np.ndarray:
        """:math:`\\mathbf E` (n x 2): :math:`\\mathbf Q = \\mathbf E[\\tau_\\psi, \\tau_\\theta]^{\\mathsf T}`."""
        E = np.zeros((self.n, 2))
        E[self.struttura.j_psi_m, 0] = 1.0
        E[self.struttura.j_theta_m, 1] = 1.0
        return E

    def M(self, q: np.ndarray) -> np.ndarray:
        return np.asarray(self.funzioni.M(q, self.parametri), dtype=float)

    def accelerazioni(self, q: np.ndarray, dq: np.ndarray, tau: np.ndarray) -> np.ndarray:
        """:math:`\\ddot{\\mathbf q} = \\mathbf M^{-1}(\\mathbf E\\boldsymbol\\tau - \\mathbf C\\dot{\\mathbf q}
        - \\mathbf C_d\\dot{\\mathbf q} - \\partial V/\\partial\\mathbf q)`."""
        f = np.asarray(self.funzioni.forze(q, dq, self.parametri), dtype=float).reshape(-1)
        return np.linalg.solve(self.M(q), self.matrice_ingressi() @ tau - f)

    def energia(self, q: np.ndarray, dq: np.ndarray) -> float:
        return float(self.funzioni.T(q, dq, self.parametri) + self.funzioni.V(q, self.parametri))

    def momento_psi(self, q: np.ndarray, dq: np.ndarray) -> float:
        """:math:`p_\\psi = \\sum_i \\partial\\mathcal L/\\partial\\dot\\psi_i`."""
        return float(self.funzioni.p_psi(q, dq, self.parametri))

    def senza_dissipazione(self) -> "ModelloAccoppiato":
        """Copia con :math:`\\mathbf C = 0` (scenario S3)."""
        mod = self.simbolico
        p = self.parametri.copy()
        for S in (mod.C_psi, mod.C_theta):
            for s in S.free_symbols:
                p[mod.parametri.index(s)] = 0.0
        return ModelloAccoppiato(self.struttura, p, self.rid_psi, self.rid_theta, self.payload,
                                 0.0, 0.0, self.funzioni)


THETA_RIFERIMENTO_BETA: Final[float] = 0.0
"""Punto di lavoro di riferimento per tarare :math:`\\beta`: lo smorzamento
strutturale non deve dipendere dal puntamento, quindi si tara una volta sola."""


@dataclass(frozen=True)
class Sistema:
    """Catena completa dei quattro stadi per un modello di massa e una forma ridotta."""

    par: object
    modello_massa: str
    forma_massa: str
    asse_psi: object
    asse_theta: object
    rid_psi: ModelloRidotto
    rid_theta: ModelloRidotto
    payload: Payload
    beta_psi: float
    beta_theta: float
    accoppiato: ModelloAccoppiato

    def aggiunte_isolate(self, asse: str, theta0: float) -> tuple[float, float]:
        """``(J_aggiunta, k_aggiunta)`` del modello isolato (sezione 7.3)."""
        pl = self.payload
        if asse == "azimut":
            return pl.J_ax + float(pl.J_psi(theta0)), 0.0
        return pl.I_yy, float(pl.k_g(theta0))

    def isolato_completo(self, asse: str, theta0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(M, K, C)`` del modello FEM completo isolato in ``theta0``."""
        from . import fem

        mod = self.asse_psi if asse == "azimut" else self.asse_theta
        J, k = self.aggiunte_isolate(asse, theta0)
        M, K = fem.matrici_isolate(mod, self.modello_massa, J_aggiunta=J, k_aggiunta=k)  # type: ignore[arg-type]
        beta = self.beta_psi if asse == "azimut" else self.beta_theta
        C = fem.matrice_smorzamento(mod.K, mod.C_b, beta)  # type: ignore[attr-defined]
        return M, K, C

    def isolato_ridotto(self, asse: str, theta0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(M_r, K_r, C_r)`` della catena ridotta isolata in ``theta0``."""
        from . import reduction as rd

        rid = self.rid_psi if asse == "azimut" else self.rid_theta
        J, k = self.aggiunte_isolate(asse, theta0)
        M, K = rd.matrici_ridotte_isolate(rid, J_aggiunta=J, k_aggiunta=k)
        beta = self.beta_psi if asse == "azimut" else self.beta_theta
        return M, K, beta * rid.K_r + rid.C_b_r


def costruisci_sistema(par, modello_massa: str = "consistent", forma_massa: str = "guyan", *,
                       payload: Payload | None = None, assi: tuple | None = None) -> Sistema:
    """Esegue i quattro stadi della sezione 0 e restituisce il :class:`Sistema`.

    ``payload`` e ``assi`` permettono di riusare oggetti gia' costruiti (il
    preprocessing del 2 GdL e l'assemblaggio non dipendono dal modello di massa).
    """
    from . import fem
    from . import payload as pl_mod
    from . import reduction as rd

    if assi is None:
        assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
                fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    mod_psi, mod_theta = assi
    pl = payload if payload is not None else pl_mod.costruisci_payload(par, mod_theta)
    rid_psi = rd.riduci(mod_psi, modello_massa, forma_massa)  # type: ignore[arg-type]
    rid_theta = rd.riduci(mod_theta, modello_massa, forma_massa)  # type: ignore[arg-type]

    th = THETA_RIFERIMENTO_BETA
    M_psi, _ = fem.matrici_isolate(mod_psi, modello_massa,  # type: ignore[arg-type]
                                   J_aggiunta=pl.J_ax + float(pl.J_psi(th)))
    M_th, _ = fem.matrici_isolate(mod_theta, modello_massa, J_aggiunta=pl.I_yy)  # type: ignore[arg-type]
    beta_psi = fem.beta_smorzamento(M_psi, mod_psi.K, par.zeta)
    beta_theta = fem.beta_smorzamento(M_th, mod_theta.K, par.zeta)

    acc = costruisci_accoppiato(rid_psi, rid_theta, pl, beta_psi=beta_psi, beta_theta=beta_theta)
    return Sistema(par=par, modello_massa=modello_massa, forma_massa=forma_massa,
                   asse_psi=mod_psi, asse_theta=mod_theta, rid_psi=rid_psi, rid_theta=rid_theta,
                   payload=pl, beta_psi=beta_psi, beta_theta=beta_theta, accoppiato=acc)


def costruisci_accoppiato(rid_psi: ModelloRidotto, rid_theta: ModelloRidotto, pl: Payload,
                          *, beta_psi: float, beta_theta: float) -> ModelloAccoppiato:
    """Assegna i parametri numerici al modello simbolico.

    Lo smorzamento ridotto e' :math:`\\mathbf C_r = \\beta\\mathbf K_r + \\mathbf C_{b,r}`,
    con :math:`\\beta` tarato sull'asse isolato completo (sezione 7.3). Telescopio
    e :math:`J^{ax}` entrano solo attraverso :math:`T` e :math:`V`: le masse
    ridotte sono quelle degli alberi e delle inerzie nodali, senza aggiunte.
    """
    st = struttura_da(rid_psi, rid_theta)
    mod = costruisci_simbolico(st)
    C_psi = beta_psi * rid_psi.K_r + rid_psi.C_b_r
    C_theta = beta_theta * rid_theta.K_r + rid_theta.C_b_r
    p = valori_parametri(mod, rid_psi, rid_theta, pl, C_psi, C_theta)
    return ModelloAccoppiato(st, p, rid_psi, rid_theta, pl, beta_psi, beta_theta,
                             funzioni_numeriche(st))
