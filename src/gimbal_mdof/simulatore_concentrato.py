"""Simulatore non lineare del modello a parametri concentrati (3 + 4 GdL).

Integra le equazioni di Newton-Eulero di :func:`gimbal_mdof.concentrato.newton_eulero`,
con tutte le funzioni trigonometriche di :math:`\\theta_p`:

.. math::
    \\mathbf M(\\theta_p)\\,\\ddot{\\mathbf q} = \\mathbf E\\boldsymbol\\tau - \\mathbf h(\\mathbf q,
    \\dot{\\mathbf q}),\\qquad \\mathbf M = \\frac{\\partial\\mathbf r}{\\partial\\ddot{\\mathbf q}},\\quad
    \\mathbf h = \\mathbf r\\big|_{\\ddot{\\mathbf q} = 0,\\ \\boldsymbol\\tau = 0}

dove :math:`\\mathbf r` e' il residuo delle equazioni (``lhs - rhs``). Le espressioni sono
compilate con ``lambdify``. :class:`SimulatoreNE` offre l'interfaccia usata da
:mod:`gimbal_mdof.simulate` (``accelerazioni``, ``energia``, ``momento_psi``,
``direzione_rigida``, ``senza_dissipazione``, ``linearizza``), quindi gli scenari S1-S3 si
applicano senza modifiche; :func:`S4` e' uno scenario a grande ampiezza, in cui le
non linearita' trigonometriche contano.

Energia, momento :math:`p_\\psi` e linearizzazione sono quelli del modello di Lagrange
equivalente (:mod:`gimbal_mdof.coupled` con la catena concentrata), di cui le equazioni di
Newton-Eulero sono verificate identiche (:func:`gimbal_mdof.concentrato.verifica_contro_lagrange`).

:func:`esporta_codice` scrive il simulatore come codice autonomo (Python senza SymPy e
funzione MATLAB/Octave), con i parametri numerici e le equazioni esplicite.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Callable, Final

import numpy as np
import sympy as sp

from . import concentrato as cc
from . import coupled as cp
from . import linearize as lz
from .simulate import RTOL, ATOL, RisultatoScenario, integra


@dataclass(frozen=True)
class FunzioniNE:
    M: Callable
    """``M(q, p)``: matrice di massa (7 x 7), dipende solo da theta_p."""
    h: Callable
    """``h(q, dq, p)``: termini giroscopici, centrifughi, elastici, dissipativi e gravita'."""
    E: np.ndarray
    """Ingressi: tau_psi sul nodo motore di azimut, tau_theta su quello polare."""


@lru_cache(maxsize=2)
def funzioni_NE(variante: str = "catena") -> FunzioniNE:
    ne = cc.newton_eulero(variante)
    r = ne.residuo
    M = r.jacobian(ne.ddq)
    if any(s in M.free_symbols for s in list(ne.dq) + [ne.q[i] for i in range(ne.n) if i != ne.i_theta_p]):
        raise AssertionError("M dipende da coordinate diverse da theta_p")
    zero = {s: 0 for s in ne.ddq} | {cp.tau_psi: 0, cp.tau_theta: 0}
    h = r.subs(zero)
    E = -np.array(r.jacobian(sp.Matrix([cp.tau_psi, cp.tau_theta])), dtype=float)
    p = list(ne.simboli_parametri)
    lam = lambda args, e: sp.lambdify(args, e, modules="numpy", cse=True)  # noqa: E731
    return FunzioniNE(M=lam([list(ne.q), p], M), h=lam([list(ne.q), list(ne.dq), p], h), E=E)


def vettore_parametri(s: cp.Sistema) -> np.ndarray:
    ne = cc.newton_eulero(cc.variante_di(s))
    v = cc.valori_parametri_NE(s)
    return np.array([v[sym] for sym in ne.simboli_parametri], dtype=float)


@dataclass(frozen=True)
class SimulatoreNE:
    """Modello non lineare a parametri concentrati integrato dalle equazioni di Newton-Eulero."""

    sistema: cp.Sistema
    p: np.ndarray
    funzioni: FunzioniNE
    accoppiato: cp.ModelloAccoppiato
    """Modello di Lagrange equivalente (energia, p_psi, linearizzazione)."""

    @property
    def n(self) -> int:
        return self.accoppiato.n

    @property
    def variante(self) -> str:
        return cc.variante_di(self.sistema)

    @property
    def struttura(self) -> cp.Struttura:
        return self.accoppiato.struttura

    @property
    def segnaposto(self) -> bool:
        return self.sistema.par.segnaposto

    def matrice_ingressi(self) -> np.ndarray:
        return self.funzioni.E

    def matrice_uscite(self) -> np.ndarray:
        return lz.costruisci_uscite(self.accoppiato)

    def M(self, q: np.ndarray) -> np.ndarray:
        return np.asarray(self.funzioni.M(q, self.p), dtype=float)

    def accelerazioni(self, q: np.ndarray, dq: np.ndarray, tau: np.ndarray) -> np.ndarray:
        h = np.asarray(self.funzioni.h(q, dq, self.p), dtype=float).reshape(-1)
        return np.linalg.solve(self.M(q), self.funzioni.E @ tau - h)

    def campo(self, tau: Callable[[float], np.ndarray]) -> Callable[[float, np.ndarray], np.ndarray]:
        """Campo vettoriale :math:`\\dot{\\mathbf x} = \\mathbf f(t, \\mathbf x)`, :math:`\\mathbf x = [\\mathbf q, \\dot{\\mathbf q}]`."""
        def f(t: float, x: np.ndarray) -> np.ndarray:
            n = self.n
            return np.concatenate([x[n:], self.accelerazioni(x[:n], x[n:], tau(t))])
        return f

    def energia(self, q: np.ndarray, dq: np.ndarray) -> float:
        return self.accoppiato.energia(q, dq)

    def momento_psi(self, q: np.ndarray, dq: np.ndarray) -> float:
        return self.accoppiato.momento_psi(q, dq)

    def direzione_rigida(self, asse: str) -> np.ndarray:
        return self.accoppiato.direzione_rigida(asse)

    def linearizza(self, theta0: float, Omega0: float) -> lz.ModelloLineare:
        return lz.linearizza(self.accoppiato, theta0, Omega0)

    def senza_dissipazione(self) -> "SimulatoreNE":
        ne = cc.newton_eulero(self.variante)
        p = self.p.copy()
        for sym in list(ne.c.values()) + list(ne.b.values()):
            p[ne.simboli_parametri.index(sym)] = 0.0
        return replace(self, p=p, accoppiato=self.accoppiato.senza_dissipazione())


def costruisci_simulatore(s: cp.Sistema) -> SimulatoreNE:
    if s.forma_massa != cc.FORMA:
        raise ValueError("il simulatore Newton-Eulero richiede il sistema a parametri concentrati")
    return SimulatoreNE(sistema=s, p=vettore_parametri(s), funzioni=funzioni_NE(cc.variante_di(s)),
                        accoppiato=s.accoppiato)


# ---------------------------------------------------------------------------
# Scenario a grande ampiezza
# ---------------------------------------------------------------------------

AMPIEZZA_S4: Final[float] = float(np.deg2rad(10.0))


def S4(sim: SimulatoreNE, theta0: float, Omega0: float, *, delta_theta: float = AMPIEZZA_S4,
       t_fine: float = 4.0, n_punti: int = 4001) -> RisultatoScenario:
    """Oscillazione polare di grande ampiezza durante l'inseguimento.

    Stato iniziale: equilibrio stazionario in :math:`(\\theta_0, \\Omega_0)` con la coppia
    nominale :math:`\\boldsymbol\\tau_0` applicata, piu' una rotazione rigida
    ``delta_theta`` di tutto l'asse polare. Il telescopio oscilla come un pendolo attorno
    a :math:`\\theta_0` (gravita' :math:`M_{g,J}(\\theta_p)`), e l'azimut risponde con
    l'inerzia variabile :math:`J_{\\psi,p}(\\theta_p)` e i termini giroscopici. Il confronto
    con il modello linearizzato, integrato dallo stesso stato iniziale, mostra l'effetto
    delle funzioni trigonometriche: lo scarto cresce con l'ampiezza.
    """
    ml = sim.linearizza(theta0, Omega0)
    q0, tau0 = ml.punto.q0, ml.punto.tau0
    u_psi, u_th = sim.direzione_rigida("azimut"), sim.direzione_rigida("polare")
    t = np.linspace(0.0, t_fine, n_punti)
    x0 = np.concatenate([q0 + delta_theta * u_th, Omega0 * u_psi])
    X = integra(sim.campo(lambda tt: tau0), x0, t)
    n = sim.n
    dq_nl = X[:, :n] - q0[None] - Omega0 * t[:, None] * u_psi[None]

    def campo_lin(tt: float, x: np.ndarray) -> np.ndarray:
        return ml.A @ x

    X_lin = integra(campo_lin, np.concatenate([delta_theta * u_th, np.zeros(n)]), t)
    Cq = sim.matrice_uscite()
    y_nl, y_lin = dq_nl @ Cq.T, X_lin[:, :n] @ Cq.T
    errori = {}
    for i, u in enumerate(lz.USCITE):
        picco = float(np.max(np.abs(y_lin[:, i])))
        errori[u] = float(np.max(np.abs(y_nl[:, i] - y_lin[:, i])) / picco) if picco > 0 else float("nan")
    E = np.array([sim.energia(x[:n], x[n:]) for x in X[:: max(1, n_punti // 400)]])
    return RisultatoScenario(nome="S4", descrizione=f"polar oscillation of {np.rad2deg(delta_theta):g} deg "
                             "while tracking", t=t, y_nl=y_nl, y_lin=y_lin, theta0=theta0,
                             Omega0=Omega0, ampiezza=delta_theta, errori=errori, energia=E)


# ---------------------------------------------------------------------------
# Esportazione del simulatore come codice autonomo
# ---------------------------------------------------------------------------


def _nomi_semplici(variante: str = "catena") -> dict[sp.Symbol, sp.Symbol]:
    """Simboli LaTeX -> identificatori, con le componenti del telescopio nella terna 2
    (I_II, I_JJ, I_KK, I_IJ, I_IK, I_JK; d_I, d_J, d_K), come nel documento."""
    ne = cc.newton_eulero(variante)
    sost: dict[sp.Symbol, sp.Expr] = {}
    for vett, pre in ((ne.q, ""), (ne.dq, "d"), (ne.ddq, "dd")):
        for s, nome in zip(vett, ne.nomi):
            sost[s] = sp.Symbol(pre + nome, real=True)
    for d, pre in ((ne.J, "J_"), (ne.k, "k_"), (ne.c, "c_"), (ne.b, "b_")):
        for nome, s in d.items():
            sost[s] = sp.Symbol(pre + nome, real=True)
    I = {n: sp.Symbol(n, real=True) for n in ("I_II", "I_JJ", "I_KK", "I_IJ", "I_IK", "I_JK",
                                               "d_I", "d_J", "d_K", "m_p", "g", "J_ax")}
    sost |= {cp.I_xx: I["I_KK"], cp.I_yy: I["I_JJ"], cp.I_zz: I["I_II"], cp.I_xy: -I["I_JK"],
             cp.I_xz: I["I_IK"], cp.I_yz: -I["I_IJ"], cp.d_x: I["d_K"], cp.d_y: -I["d_J"],
             cp.d_z: I["d_I"], cp.m_p: I["m_p"], cp.g: I["g"], cp.J_ax: I["J_ax"],
             cp.tau_psi: sp.Symbol("tau_psi", real=True), cp.tau_theta: sp.Symbol("tau_theta", real=True)}
    return sost  # type: ignore[return-value]


def parametri_documento(sim: SimulatoreNE) -> dict[str, float]:
    """Valori numerici con i nomi del documento (terna 2 per il telescopio)."""
    ne = cc.newton_eulero(sim.variante)
    v = cc.valori_parametri_NE(sim.sistema)
    out: dict[str, float] = {}
    for d, pre in ((ne.J, "J_"), (ne.k, "k_"), (ne.c, "c_"), (ne.b, "b_")):
        for nome, s in d.items():
            out[pre + nome] = float(v[s])
    pl = sim.sistema.payload
    I2, d2 = pl.I_2, pl.d_2
    out |= {"I_II": I2[0, 0], "I_JJ": I2[1, 1], "I_KK": I2[2, 2], "I_IJ": I2[0, 1], "I_IK": I2[0, 2],
            "I_JK": I2[1, 2], "d_I": d2[0], "d_J": d2[1], "d_K": d2[2], "m_p": pl.m, "g": pl.g,
            "J_ax": pl.J_ax}
    return {k: float(x) for k, x in out.items()}


def esporta_codice(sim: SimulatoreNE, cartella: str | Path) -> list[Path]:
    """Scrive ``ne_concentrato.py`` (NumPy/SciPy, senza SymPy) e ``ne_concentrato.m`` (MATLAB/Octave)."""
    from sympy.printing.octave import octave_code
    from sympy.printing.pycode import pycode

    ne = cc.newton_eulero(sim.variante)
    sost = _nomi_semplici(sim.variante)
    n = ne.n
    q_txt = ", ".join(ne.nomi)
    gdl = "3 + 5 GdL, telescopio appeso al mozzo" if sim.variante == "mozzo" else "3 + 4 GdL"
    r = ne.residuo.xreplace(sost)
    ddq = [sost[s] for s in ne.ddq]
    M = r.jacobian(sp.Matrix(ddq)).applyfunc(sp.simplify)
    h = r.subs({s: 0 for s in ddq} | {sost[cp.tau_psi]: 0, sost[cp.tau_theta]: 0}).applyfunc(sp.expand)
    par = parametri_documento(sim)
    q = [sost[s] for s in ne.q]
    dq = [sost[s] for s in ne.dq]
    cartella = Path(cartella)
    cartella.mkdir(parents=True, exist_ok=True)
    avviso = "risultato da parametri segnaposto" if sim.segnaposto else ""

    # --- Python ---------------------------------------------------------------
    L = [f'"""Modello a parametri concentrati ({gdl}): equazioni di Newton-Eulero non lineari.',
         "",
         "Generato da gimbal_mdof.simulatore_concentrato.esporta_codice. Non modificare a mano.",
         f"q = [{q_txt}] (rad); x = [q, dq].",
         "M(theta_p) ddq = E tau - h(q, dq);  tau = [tau_psi, tau_theta] (N m).",
         "Telescopio con le componenti nella terna 2 (I_hat, J_hat, K_hat; K_hat = puntamento).",
         (f"AVVISO: {avviso}." if avviso else ""), '"""', "",
         "import numpy as np", "from numpy import sin, cos", "",
         "PARAMETRI = {"]
    L += [f"    {k!r}: {v!r}," for k, v in par.items()]
    L += ["}", "", "", "def massa(q, p=PARAMETRI):"]
    L += ["    " + ", ".join(str(s) for s in q) + " = q"]
    L += [f"    {k} = p[{k!r}]" for k in par]
    L += ["    return np.array(["]
    L += ["        [" + ", ".join(pycode(M[i, j]).replace("math.", "") for j in range(n)) + "]," for i in range(n)]
    L += ["    ])", "", "", "def termini(q, dq, p=PARAMETRI):",
          '    """h(q, dq): giroscopici, centrifughi, elastici, dissipativi, gravita\'."""']
    L += ["    " + ", ".join(str(s) for s in q) + " = q", "    " + ", ".join(str(s) for s in dq) + " = dq"]
    L += [f"    {k} = p[{k!r}]" for k in par]
    L += ["    return np.array(["] + ["        " + pycode(h[i]).replace("math.", "") + "," for i in range(n)]
    L += ["    ])", "", "",
          f"E = np.array({np.array2string(sim.funzioni.E, separator=', ')})", "", "",
          "def rhs(t, x, tau, p=PARAMETRI):",
          '    """dx/dt per scipy.integrate.solve_ivp; tau(t) restituisce [tau_psi, tau_theta]."""',
          f"    q, dq = x[:{n}], x[{n}:]",
          "    ddq = np.linalg.solve(massa(q, p), E @ np.asarray(tau(t)) - termini(q, dq, p))",
          "    return np.concatenate([dq, ddq])", ""]
    p_py = cartella / "ne_concentrato.py"
    p_py.write_text("\n".join(L), encoding="utf-8")

    # --- MATLAB / Octave ------------------------------------------------------
    m = ["function dx = ne_concentrato(t, x, tau, p)",
         f"% Modello a parametri concentrati ({gdl}): equazioni di Newton-Eulero non lineari.",
         "% Generato da gimbal_mdof.simulatore_concentrato.esporta_codice. Non modificare a mano.",
         f"% x = [q; dq], q = [{' '.join(ne.nomi)}]' (rad);",
         "% tau = @(t) [tau_psi; tau_theta] (N m); p = ne_concentrato_parametri().",
         "% Uso: [t, X] = ode45(@(t, x) ne_concentrato(t, x, tau, p), [0 T], x0);",
         (f"% AVVISO: {avviso}." if avviso else "%"),
         f"q = x(1:{n}); dq = x({n + 1}:{2 * n});"]
    m += [f"{s} = q({i + 1}); d{s} = dq({i + 1});" for i, s in enumerate(map(str, q))]
    m += [f"{k} = p.{k};" for k in par]
    m += [f"M = zeros({n}, {n});"]
    m += [f"M({i + 1},{j + 1}) = {octave_code(M[i, j])};" for i in range(n) for j in range(n) if M[i, j] != 0]
    m += [f"h = zeros({n}, 1);"] + [f"h({i + 1}) = {octave_code(h[i])};" for i in range(n)]
    m += [f"E = {octave_code(sp.Matrix(sim.funzioni.E))};",
          "u = tau(t); u = u(:);", "ddq = M \\ (E * u - h);", "dx = [dq; ddq];", "end", ""]
    p_m = cartella / "ne_concentrato.m"
    p_m.write_text("\n".join(m), encoding="utf-8")
    pm = ["function p = ne_concentrato_parametri()",
          "% Parametri numerici del modello a parametri concentrati (SI)."
          + (f" AVVISO: {avviso}." if avviso else "")]
    pm += [f"p.{k} = {v!r};" for k, v in par.items()] + ["end", ""]
    p_pm = cartella / "ne_concentrato_parametri.m"
    p_pm.write_text("\n".join(pm), encoding="utf-8")
    return [p_py, p_m, p_pm]


__all__ = ["SimulatoreNE", "costruisci_simulatore", "funzioni_NE", "S4", "esporta_codice",
           "parametri_documento", "RTOL", "ATOL"]
