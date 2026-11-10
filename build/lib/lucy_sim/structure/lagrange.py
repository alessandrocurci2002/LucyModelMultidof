r"""Coupled model of the two lumped chains: Lagrange equations derived with SymPy, and linearisation.

Coordinates :math:`q = [q_\psi, q_\theta]`, with :math:`q_\psi = [\psi_e, \psi_m, \psi_L, \psi_p]` and
:math:`q_\theta = [\theta_e, \theta_h, \theta_m, \theta_r, \theta_p]`. The fork turns with :math:`\psi_L`, the
telescope about z_hat with :math:`\psi_p` (hanging from L) and about the polar axis with :math:`\theta_p`, relative
to the fork.

.. math::
    T = \tfrac12\dot q_\psi^T M_\psi\dot q_\psi + \tfrac12 J^{ax}\dot\psi_L^2 + \tfrac12\dot q_\theta^T M_\theta\dot q_\theta
      + \tfrac12\omega_p^T I_{p,O}\,\omega_p,\qquad
    V = \tfrac12 q_\psi^T K_\psi q_\psi + \tfrac12 q_\theta^T K_\theta q_\theta + m_p g(-d_x\sin\theta_p + d_z\cos\theta_p).

The Lagrange equations are not written by hand: SymPy derives them, with the reduced matrices as symbols with the
structure of the springs, so that the derivation holds for any geometry; the numbers are substituted at the end.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable

import numpy as np
import sympy as sp

from .lumped import LumpedChain

I_xx, I_yy, I_zz = sp.symbols("I_xx I_yy I_zz", positive=True)
I_xy, I_xz, I_yz = sp.symbols("I_xy I_xz I_yz", real=True)
g = sp.Symbol("g", positive=True)
d_x, d_y, d_z = sp.symbols("d_x d_y d_z", real=True)
tau_psi, tau_theta = sp.symbols("tau_psi tau_theta", real=True)
m_p = sp.Symbol("m_p", positive=True)
J_ax = sp.Symbol("J^{ax}", nonnegative=True)


def _R_z(a: sp.Expr) -> sp.Matrix:
    c, s = sp.cos(a), sp.sin(a)
    return sp.Matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def _R_y(a: sp.Expr) -> sp.Matrix:
    c, s = sp.cos(a), sp.sin(a)
    return sp.Matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def _vee(S: sp.Matrix) -> sp.Matrix:
    """Inverse of the skew operator, :math:`[w]_\\times v = w\\times v`."""
    return sp.Matrix([S[2, 1], S[0, 2], S[1, 0]])


def _tensor() -> sp.Matrix:
    return sp.Matrix([[I_xx, I_xy, I_xz], [I_xy, I_yy, I_yz], [I_xz, I_yz, I_zz]])


@dataclass(frozen=True)
class Structure:
    """Structure of the two chains: key of the symbolic derivation."""

    roles_psi: tuple[str, ...]
    roles_theta: tuple[str, ...]
    i_psi_e: int
    i_psi_m: int
    i_psi_L: int
    i_theta_e: int
    i_theta_p: int
    i_theta_m: int
    arcs_theta: tuple[tuple[int, int], ...]
    i_psi_p: int

    @property
    def arcs_psi(self) -> tuple[tuple[int, int], ...]:
        return tuple((i, i + 1) for i in range(self.n_psi - 1))

    @property
    def n_psi(self) -> int:
        return len(self.roles_psi)

    @property
    def n_theta(self) -> int:
        return len(self.roles_theta)

    @property
    def n(self) -> int:
        return self.n_psi + self.n_theta

    @property
    def j_theta_e(self) -> int:
        return self.n_psi + self.i_theta_e

    @property
    def j_theta_p(self) -> int:
        return self.n_psi + self.i_theta_p

    @property
    def j_theta_m(self) -> int:
        return self.n_psi + self.i_theta_m


def structure_of(az: LumpedChain, pol: LumpedChain) -> Structure:
    return Structure(az.roles, pol.roles, az.index("e"), az.index("m"), az.index("L"), pol.index("e"), pol.index("p"),
                     pol.index("m"), pol.arcs, az.index("p"))


def _spring_matrix(name: str, sup: str, roles: tuple[str, ...], arcs: tuple[tuple[int, int], ...]) -> tuple[sp.Matrix, list[sp.Symbol]]:
    """Symmetric matrix of symbols ``name^{sup}_{ij}`` with the structure of the springs (diagonal plus ``arcs``)."""
    n = len(roles)
    A = sp.zeros(n, n)
    symbols = []
    for i, j in sorted({(i, i) for i in range(n)} | {(min(a, b), max(a, b)) for a, b in arcs}):
        s = sp.Symbol(f"{name}^{{{sup}}}_{{{roles[i]}{roles[j]}}}", real=True)
        A[i, j] = A[j, i] = s
        symbols.append(s)
    return A, symbols


def _time_derivative(expr: sp.Expr, q: sp.Matrix, dq: sp.Matrix, ddq: sp.Matrix) -> sp.Expr:
    return sum((sp.diff(expr, q[i]) * dq[i] + sp.diff(expr, dq[i]) * ddq[i] for i in range(len(q))), sp.S.Zero)


def _angular_velocity(psi: sp.Symbol, theta: sp.Symbol, dpsi: sp.Symbol, dtheta: sp.Symbol) -> sp.Matrix:
    """:math:`\\omega_p` in frame P from Poisson's formula, :math:`[\\omega]_\\times = R^T\\dot R`, with
    :math:`R = R_z(\\psi) R_y(-\\theta)`; checked against the closed form
    :math:`(\\dot\\psi\\sin\\theta, -\\dot\\theta, \\dot\\psi\\cos\\theta)`."""
    t = sp.Symbol("t", real=True)
    a = sp.Function("a", real=True)(t)
    b = sp.Function("b", real=True)(t)
    R = _R_z(a) * _R_y(-b)
    w = sp.simplify(_vee(sp.simplify(R.T * sp.diff(R, t))))
    w = sp.simplify(w.subs({sp.Derivative(a, t): dpsi, sp.Derivative(b, t): dtheta}).subs({a: psi, b: theta}))
    if sp.simplify(w - sp.Matrix([dpsi * sp.sin(theta), -dtheta, dpsi * sp.cos(theta)])) != sp.zeros(3, 1):
        raise AssertionError("angular velocity from Poisson's formula differs from the closed form")
    return w


def _payload_energies(psi, theta, dpsi, dtheta):
    """Kinetic energy and gravity potential of the rigid telescope."""
    w = _angular_velocity(psi, theta, dpsi, dtheta)
    half = sp.Rational(1, 2)
    T = sp.expand(sp.trigsimp(sp.expand(half * (w.T * _tensor() * w)[0, 0])))
    z_in_P = _R_y(-theta).T * sp.Matrix([0, 0, 1])
    V = sp.simplify(m_p * g * (sp.Matrix([d_x, d_y, d_z]).T * z_in_P)[0, 0])
    return w, T, V


def _christoffel(M: sp.Matrix, q: sp.Matrix, dq: sp.Matrix) -> sp.Matrix:
    """Coriolis matrix from the Christoffel symbols of the first kind."""
    n = len(q)
    dM = [M.diff(q[i]) for i in range(n)]
    active = [i for i in range(n) if not dM[i].is_zero_matrix]
    C = sp.zeros(n, n)
    for k in range(n):
        for j in range(n):
            total = sp.S.Zero
            for i in range(n):
                c = sp.S.Zero
                if i in active:
                    c += dM[i][k, j]
                if j in active:
                    c += dM[j][k, i]
                if k in active:
                    c -= dM[k][i, j]
                if c != 0:
                    total += sp.Rational(1, 2) * c * dq[i]
            C[k, j] = sp.simplify(total)
    return C


@dataclass(frozen=True)
class SymbolicModel:
    """Equations of motion: ``residual = M ddq + C dq + Cd dq + dV/dq``, generalised forces ``Q = E tau``."""

    structure: Structure
    q: sp.Matrix
    dq: sp.Matrix
    ddq: sp.Matrix
    M_psi: sp.Matrix
    K_psi: sp.Matrix
    C_psi: sp.Matrix
    M_theta: sp.Matrix
    K_theta: sp.Matrix
    C_theta: sp.Matrix
    T: sp.Expr
    V: sp.Expr
    M: sp.Matrix
    C: sp.Matrix
    Cd: sp.Matrix
    gradV: sp.Matrix
    residual: sp.Matrix
    parameters: tuple[sp.Symbol, ...]


@lru_cache(maxsize=4)
def symbolic_model(st: Structure) -> SymbolicModel:
    """Derives the equations of motion with SymPy and checks the matrix form against the Lagrange equations."""
    rp, rt = st.roles_psi, st.roles_theta
    q_psi = sp.Matrix([sp.Symbol(f"\\psi_{{{r}}}", real=True) for r in rp])
    q_th = sp.Matrix([sp.Symbol(f"\\theta_{{{r}}}", real=True) for r in rt])
    dq_psi = sp.Matrix([sp.Symbol(f"\\dot{{\\psi}}_{{{r}}}", real=True) for r in rp])
    dq_th = sp.Matrix([sp.Symbol(f"\\dot{{\\theta}}_{{{r}}}", real=True) for r in rt])
    ddq_psi = sp.Matrix([sp.Symbol(f"\\ddot{{\\psi}}_{{{r}}}", real=True) for r in rp])
    ddq_th = sp.Matrix([sp.Symbol(f"\\ddot{{\\theta}}_{{{r}}}", real=True) for r in rt])
    q, dq, ddq = sp.Matrix.vstack(q_psi, q_th), sp.Matrix.vstack(dq_psi, dq_th), sp.Matrix.vstack(ddq_psi, ddq_th)
    M_psi, sM_psi = _spring_matrix("M", "\\psi", rp, st.arcs_psi)
    K_psi, sK_psi = _spring_matrix("K", "\\psi", rp, st.arcs_psi)
    C_psi, sC_psi = _spring_matrix("C", "\\psi", rp, st.arcs_psi)
    M_th, sM_th = _spring_matrix("M", "\\theta", rt, st.arcs_theta)
    K_th, sK_th = _spring_matrix("K", "\\theta", rt, st.arcs_theta)
    C_th, sC_th = _spring_matrix("C", "\\theta", rt, st.arcs_theta)
    dpsi_L = dq[st.i_psi_L]
    theta_p, dtheta_p = q[st.j_theta_p], dq[st.j_theta_p]
    # the telescope turns about z_hat with psi_p (hanging from L); J^ax (polar shafts, hub, ring) stays on the fork
    _, T_payload, V_g = _payload_energies(q[st.i_psi_p], theta_p, dq[st.i_psi_p], dtheta_p)
    half = sp.Rational(1, 2)
    T = (half * (dq_psi.T * M_psi * dq_psi)[0, 0] + half * J_ax * dpsi_L**2 + half * (dq_th.T * M_th * dq_th)[0, 0] + T_payload)
    V = half * (q_psi.T * K_psi * q_psi)[0, 0] + half * (q_th.T * K_th * q_th)[0, 0] + V_g
    D = half * (dq_psi.T * C_psi * dq_psi)[0, 0] + half * (dq_th.T * C_th * dq_th)[0, 0]
    L = T - V
    # generalised forces from the virtual work: azimuth rotor dpsi_m k (base fixed); polar rotor dpsi_L k + dtheta_m J,
    # polar stator (fork) dpsi_L k; in frame 1 k_hat = (0, 0, 1) and the polar axis J_hat = -j_hat = (0, -1, 0)
    dlt = sp.Matrix([sp.Symbol(f"\\delta {s.name}", real=True) for s in q])
    Zv, Jh = sp.Matrix([0, 0, 1]), sp.Matrix([0, -1, 0])
    rot_psi = dlt[st.i_psi_m] * Zv
    rot_th = dlt[st.i_psi_L] * Zv + dlt[st.j_theta_m] * Jh
    stat_th = dlt[st.i_psi_L] * Zv
    dW = sp.expand((tau_psi * Zv.T * rot_psi)[0, 0] + (tau_theta * Jh.T * (rot_th - stat_th))[0, 0])
    Q = sp.Matrix([sp.diff(dW, dlt[i]) for i in range(st.n)])
    eom = sp.zeros(st.n, 1)
    for i in range(st.n):
        eom[i] = sp.expand(_time_derivative(sp.diff(L, dq[i]), q, dq, ddq) - sp.diff(L, q[i]) + sp.diff(D, dq[i]) - Q[i])
    M = sp.hessian(T, list(dq)).applyfunc(sp.simplify)
    C = _christoffel(M, q, dq)
    Cd = sp.hessian(D, list(dq))
    gradV = sp.Matrix([sp.simplify(sp.diff(V, q[i])) for i in range(st.n)])
    residual = M * ddq + C * dq + Cd * dq + gradV
    if (eom - (residual - Q)).applyfunc(lambda e: sp.simplify(sp.expand(e))) != sp.zeros(st.n, 1):
        raise AssertionError("matrix form and Lagrange equations differ")
    parameters = tuple(sM_psi + sK_psi + sC_psi + sM_th + sK_th + sC_th) + (J_ax, I_xx, I_yy, I_zz, I_xy, I_xz, I_yz, m_p, g, d_x, d_y, d_z)
    return SymbolicModel(st, q, dq, ddq, M_psi, K_psi, C_psi, M_th, K_th, C_th, T, V, M, C, Cd, gradV, residual, parameters)


@dataclass(frozen=True)
class NumericFunctions:
    M: Callable
    forces: Callable
    """:math:`C\\dot q + C_d\\dot q + \\partial V/\\partial q`."""
    M0: Callable
    C0: Callable
    K0: Callable


@lru_cache(maxsize=4)
def numeric_functions(st: Structure) -> NumericFunctions:
    """Compiled functions of the model, the parameters as argument; M0, C0, K0 are the Jacobians of the residual at
    the operating point (theta_p = theta_0, all psi rates = Omega_0, theta rates and accelerations zero)."""
    mod = symbolic_model(st)
    q, dq, p = list(mod.q), list(mod.dq), list(mod.parameters)
    th0, Om0 = sp.symbols("theta_0 Omega_0", real=True)

    def lam(args, expr):
        return sp.lambdify(args, expr, modules="numpy", cse=True, dummify=True)

    point = {mod.q[st.j_theta_p]: th0}
    for i in range(st.n):
        point[mod.dq[i]] = Om0 if i < st.n_psi else 0
        point[mod.ddq[i]] = 0
    r = mod.residual
    M0, C0, K0 = r.jacobian(mod.ddq).subs(point), r.jacobian(mod.dq).subs(point), r.jacobian(mod.q).subs(point)
    others = [s for s in q if s != mod.q[st.j_theta_p]]
    for name, A in (("M_0", M0), ("C_0", C0), ("K_0", K0)):
        if A.free_symbols & set(others):
            raise AssertionError(f"{name} depends on coordinates other than theta_p")
    return NumericFunctions(M=lam([q, p], mod.M), forces=lam([q, dq, p], mod.C * mod.dq + mod.Cd * mod.dq + mod.gradV),
                            M0=lam([th0, Om0, p], M0), C0=lam([th0, Om0, p], C0), K0=lam([th0, Om0, p], K0))


def _check_structure(A: np.ndarray, name: str, arcs: tuple[tuple[int, int], ...]) -> None:
    n = A.shape[0]
    mask = np.eye(n, dtype=bool)
    for a, b in arcs:
        mask[a, b] = mask[b, a] = True
    scale = max(float(np.max(np.abs(A))), 1e-300)
    if np.max(np.abs(np.where(mask, 0.0, A))) > 1e-10 * scale:
        raise AssertionError(f"{name}: entries outside the structure of the springs")


@dataclass(frozen=True)
class CoupledModel:
    """Coupled model with numerical parameters."""

    structure: Structure
    parameters: np.ndarray
    azimuth: LumpedChain
    polar: LumpedChain
    payload: object
    functions: NumericFunctions

    @property
    def n(self) -> int:
        return self.structure.n

    def M(self, q: np.ndarray) -> np.ndarray:
        return np.asarray(self.functions.M(q, self.parameters), dtype=float)

    def forces(self, q: np.ndarray, dq: np.ndarray) -> np.ndarray:
        return np.asarray(self.functions.forces(q, dq, self.parameters), dtype=float).reshape(-1)

    def input_matrix(self) -> np.ndarray:
        """:math:`E` (n x 2): :math:`Q = E[\\tau_\\psi, \\tau_\\theta]^T`."""
        E = np.zeros((self.n, 2))
        E[self.structure.i_psi_m, 0] = 1.0
        E[self.structure.j_theta_m, 1] = 1.0
        return E

    def without_dissipation(self) -> "CoupledModel":
        mod = symbolic_model(self.structure)
        p = self.parameters.copy()
        for S in (mod.C_psi, mod.C_theta):
            for s in S.free_symbols:
                p[mod.parameters.index(s)] = 0.0
        return CoupledModel(self.structure, p, self.azimuth, self.polar, self.payload, self.functions)

    def linearize(self, theta0: float = 0.0, Omega0: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(M0, C0, K0)`` at the operating point."""
        p = self.parameters
        return tuple(np.asarray(f(theta0, Omega0, p), dtype=float) for f in (self.functions.M0, self.functions.C0, self.functions.K0))  # type: ignore[return-value]


def coupled_model(az: LumpedChain, pol: LumpedChain, payload, *, beta_psi: float, beta_theta: float) -> CoupledModel:
    """Assigns the numbers to the symbolic model; the reduced damping is :math:`\\beta K_r + C_{b,r}`."""
    st = structure_of(az, pol)
    mod = symbolic_model(st)
    C_psi, C_theta = beta_psi * az.K + az.C_b, beta_theta * pol.K + pol.C_b
    values: dict[sp.Symbol, float] = {}
    for S, A, name, arcs in ((mod.M_psi, az.M, "M_psi", st.arcs_psi), (mod.K_psi, az.K, "K_psi", st.arcs_psi),
                             (mod.C_psi, C_psi, "C_psi", st.arcs_psi), (mod.M_theta, pol.M, "M_theta", st.arcs_theta),
                             (mod.K_theta, pol.K, "K_theta", st.arcs_theta), (mod.C_theta, C_theta, "C_theta", st.arcs_theta)):
        _check_structure(A, name, arcs)
        for i in range(A.shape[0]):
            for j in range(i, A.shape[0]):
                if S[i, j] != 0:
                    values[S[i, j]] = float(A[i, j])
    pl = payload
    values.update({J_ax: pl.J_ax, I_xx: pl.I_xx, I_yy: pl.I_yy, I_zz: pl.I_zz, I_xy: pl.I_xy, I_xz: pl.I_xz, I_yz: pl.I_yz,
                   m_p: pl.m, g: pl.g, d_x: pl.d_x, d_y: pl.d_y, d_z: pl.d_z})
    return CoupledModel(st, np.array([values[s] for s in mod.parameters], dtype=float), az, pol, payload, numeric_functions(st))
