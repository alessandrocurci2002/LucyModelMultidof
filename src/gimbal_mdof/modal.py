"""Troncamento modale alla Gawronski e modello accoppiato su base generica.

Riferimenti: Gawronski, *Modeling and Control of Antennas and Telescopes* (2008),
cap. 2 e 4; documento ``MODELLO_unificato_gawronski_dedonatis.md`` (citato come
"doc. unificato"); analisi ``ANALISI_guyan_vs_gawronski.md``.

**Idea.** Ogni riduzione usata nel pacchetto e' una base di Ritz dei GdL del FEM
completo di ciascun asse, :math:`\\mathbf q_a = \\mathbf T_a\\boldsymbol\\eta_a`:

* ``completa``: :math:`\\mathbf T = \\mathbb 1` (modello di riferimento "esatto", 9 + 10 GdL);
* ``guyan`` e ``craig_bampton``: le basi di :mod:`gimbal_mdof.reduction`;
* ``modale``: :math:`\\mathbf T = \\boldsymbol\\Phi = [\\mathbf 1, \\boldsymbol\\phi_1, \\ldots]`,
  autovettori del FEM COMPLETO con il payload congelato in :math:`\\theta_{ref}`
  (doc. unificato §3-4, costruiti pero' dal FEM completo e non dalle catene di
  Guyan, come raccomanda ``ANALISI_guyan_vs_gawronski.md`` §5.1).

Il payload e' il corpo rigido di De Donatis agganciato ai nodi L (azimut, nodo 9)
e P (polare) tramite le righe della base in quei nodi,
:math:`[\\psi_L, \\theta_p]^{\\mathsf T} = \\mathbf E_\\eta\\boldsymbol\\eta`. Le equazioni
del moto sono la (5.1) del doc. unificato:

.. math::
    \\tilde{\\mathbf M}\\ddot{\\boldsymbol\\eta} + \\tilde{\\mathbf D}\\dot{\\boldsymbol\\eta}
    + \\tilde{\\mathbf K}\\boldsymbol\\eta + \\mathbf E_\\eta^{\\mathsf T}\\big\\{[\\mathbf M_p(\\theta_p)
    - \\bar{\\mathbf M}_p]\\ddot{\\mathbf q}_p + \\mathbf C_p\\dot{\\mathbf q}_p + \\mathbf g_p\\big\\}
    = \\mathbf T^{\\mathsf T}\\mathbf B\\boldsymbol\\tau

con :math:`\\tilde{\\mathbf M} = \\mathbf T^{\\mathsf T}\\bar{\\mathbf M}\\mathbf T` (per la base
modale :math:`\\operatorname{diag}(J_{tot}, 1, \\ldots)`), :math:`\\tilde{\\mathbf K} =
\\mathbf T^{\\mathsf T}\\mathbf K\\mathbf T` (:math:`\\operatorname{diag}(0, \\omega_k^2)`) e
:math:`\\tilde{\\mathbf D} = \\mathbf T^{\\mathsf T}(\\beta\\mathbf K + \\mathbf C_b)\\mathbf T`
(:math:`\\operatorname{diag}(0, 2\\zeta_k\\omega_k) + \\boldsymbol\\Phi^{\\mathsf T}\\mathbf C_b\\boldsymbol\\Phi`).
:math:`\\bar{\\mathbf M}_p` e' l'inerzia del payload gia' contenuta nei modi (nulla per
le basi non modali). :math:`\\mathbf M_p, \\mathbf C_p, \\mathbf g_p` sono ricavati con
SymPy dalla stessa cinematica di :mod:`gimbal_mdof.coupled`
(:func:`gimbal_mdof.coupled.energie_payload`) e compilati in funzioni numeriche.

Convenzioni: indici dei nodi in base 1 in ingresso e in uscita (etichette del FEM);
:math:`\\theta > 0` alza il puntamento; tutte le grandezze in SI.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Callable, Final, Literal, Sequence

import numpy as np
import scipy.linalg as sla
import sympy as sp

from . import coupled as cp
from . import fem
from . import linearize as lz
from . import reduction as rd
from . import tf
from .params import Avviso, Parametri, ParametriModali, in_inglese
from .payload import Payload

TipoBase = Literal["completa", "guyan", "craig_bampton", "modale"]

SOGLIA_DB: Final[float] = 1.0
"""Scarto massimo di modulo sui canali diretti (T5), dB."""
SOGLIA_DEG: Final[float] = 5.0
"""Scarto massimo di fase sui canali diretti (T5), gradi."""
BANDA_ZERO: Final[float] = 0.02
"""Semiampiezza relativa della banda esclusa attorno a ogni zero di riferimento (T5)."""
TOLL_ZERO: Final[float] = 0.01
"""Scarto relativo ammesso sulla frequenza di uno zero escluso (T5)."""
SOGLIA_INCROCIATI: Final[float] = 0.12
"""Scarto relativo massimo della matrice 2x2 (norma 2) in banda (T5)."""
SOGLIA_MAGLIA: Final[float] = 0.02
"""Scarto massimo fra massa concentrata e consistente sui modi tenuti."""
N_GRIGLIA_VALIDAZIONE: Final[int] = 1500


# ---------------------------------------------------------------------------
# Payload: M_p, C_p, g_p da SymPy (doc. unificato §5)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FunzioniPayload:
    """Funzioni numeriche del payload nelle coordinate d'aggancio :math:`(\\psi_L, \\theta_p)`.

    ``M(theta, p)``: :math:`\\mathbf M_p` (2x2); ``forze(theta, dpsi, dtheta, p)``:
    :math:`\\mathbf C_p\\dot{\\mathbf q}_p + \\mathbf g_p`; ``T``, ``V``: energie;
    ``M0``, ``G0``, ``K0``: jacobiani del residuo del payload rispetto a
    :math:`\\ddot{\\mathbf q}_p, \\dot{\\mathbf q}_p, \\mathbf q_p` nel punto di lavoro
    :math:`(\\theta_0, \\dot\\psi_L = \\Omega_0, \\dot\\theta_p = 0)` (doc. unificato §7.1).
    """

    M: Callable
    forze: Callable
    T: Callable
    V: Callable
    M0: Callable
    G0: Callable
    K0: Callable
    simboli: tuple[sp.Symbol, ...]
    M_simbolica: sp.Matrix
    C_simbolica: sp.Matrix
    g_simbolica: sp.Matrix


PARAMETRI_PAYLOAD: Final[tuple[sp.Symbol, ...]] = (
    cp.I_xx, cp.I_yy, cp.I_zz, cp.I_xy, cp.I_xz, cp.I_yz, cp.m_p, cp.g, cp.d_x, cp.d_y, cp.d_z)


def parametri_payload(pl: Payload) -> np.ndarray:
    """Vettore numerico allineato a :data:`PARAMETRI_PAYLOAD`."""
    return np.array([pl.I_xx, pl.I_yy, pl.I_zz, pl.I_xy, pl.I_xz, pl.I_yz, pl.m, pl.g,
                     pl.d_x, pl.d_y, pl.d_z], dtype=float)


@lru_cache(maxsize=1)
def funzioni_payload() -> FunzioniPayload:
    """Lagrangiano del solo payload (cinematica di :mod:`gimbal_mdof.coupled`), compilato.

    Si verifica che :math:`\\mathbf M_p\\ddot{\\mathbf q} + \\mathbf C_p\\dot{\\mathbf q} +
    \\mathbf g_p` coincida con le equazioni di Lagrange di :math:`T_p - V_g`.
    """
    psi, th = sp.symbols(r"\psi_L \theta_p", real=True)
    dpsi, dth = sp.symbols(r"\dot{\psi}_L \dot{\theta}_p", real=True)
    ddpsi, ddth = sp.symbols(r"\ddot{\psi}_L \ddot{\theta}_p", real=True)
    q, dq, ddq = sp.Matrix([psi, th]), sp.Matrix([dpsi, dth]), sp.Matrix([ddpsi, ddth])
    _, T_p, V_g = cp.energie_payload(psi, th, dpsi, dth)
    M = sp.hessian(T_p, list(dq)).applyfunc(sp.simplify)
    C = cp._christoffel(M, q, dq)
    gV = sp.Matrix([sp.simplify(sp.diff(V_g, s)) for s in q])
    f = M * ddq + C * dq + gV
    lag = sp.Matrix([cp.derivata_temporale(sp.diff(T_p, dq[i]), q, dq, ddq) - sp.diff(T_p, q[i])
                     + sp.diff(V_g, q[i]) for i in range(2)])
    if (lag - f).applyfunc(lambda e: sp.simplify(sp.expand(e))) != sp.zeros(2, 1):
        raise AssertionError("payload: forma matriciale e Lagrange non coincidono")
    if any(psi in e.free_symbols for e in list(M) + list(C) + list(gV)):
        raise AssertionError("payload: le matrici dipendono da psi_L (coordinata ciclica)")
    th0, Om0 = sp.symbols("theta_0 Omega_0", real=True)
    punto = {th: th0, dpsi: Om0, dth: 0, ddpsi: 0, ddth: 0}
    M0 = f.jacobian(ddq).subs(punto)
    G0 = f.jacobian(dq).subs(punto)
    K0 = f.jacobian(q).subs(punto)
    p = list(PARAMETRI_PAYLOAD)

    def lam(args: list, e: sp.Basic) -> Callable:
        return sp.lambdify(args, e, modules="numpy", cse=True, dummify=True)

    return FunzioniPayload(
        M=lam([th, p], M), forze=lam([th, dpsi, dth, p], C * dq + gV),
        T=lam([th, dpsi, dth, p], T_p), V=lam([th, p], V_g),
        M0=lam([th0, Om0, p], M0), G0=lam([th0, Om0, p], G0), K0=lam([th0, Om0, p], K0),
        simboli=(psi, th, dpsi, dth), M_simbolica=M, C_simbolica=C, g_simbolica=gV)


# ---------------------------------------------------------------------------
# Matrici strutturali del FEM completo di ciascun asse
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AsseStrutturale:
    """FEM completo di un asse senza il telescopio (sezione 7 della specifica).

    ``M`` comprende alberi e inerzie nodali (forcella al nodo 9) e, per l'azimut,
    :math:`J^{ax}` al nodo L; ``K`` sono i soli alberi (un modo rigido); ``C =
    beta K + C_b``. ``i_agg`` e' il nodo d'aggancio del payload (L = nodo 9 per
    l'azimut, P per il polare). Indici interni (base 0).
    """

    nome: str
    modello: fem.ModelloAsse
    modello_massa: str
    M: np.ndarray
    K: np.ndarray
    C_b: np.ndarray
    beta: float
    i_e: int
    i_m: int
    i_agg: int
    i_r: int

    @property
    def n(self) -> int:
        return self.M.shape[0]

    @property
    def C(self) -> np.ndarray:
        return self.beta * self.K + self.C_b

    def etichetta(self, i: int) -> str:
        return self.modello.etichetta(i)


def assi_strutturali(par: Parametri, pl: Payload, assi: tuple[fem.ModelloAsse, fem.ModelloAsse],
                     modello_massa: str) -> tuple[AsseStrutturale, AsseStrutturale]:
    """Matrici dei due assi; :math:`\\beta` tarato come in :func:`gimbal_mdof.coupled.costruisci_sistema`."""
    mod_psi, mod_th = assi
    th = cp.THETA_RIFERIMENTO_BETA
    M_psi = mod_psi.M[modello_massa].copy()
    M_psi[mod_psi.i_estremo, mod_psi.i_estremo] += pl.J_ax
    M_th = mod_th.M[modello_massa].copy()
    Mi_psi, _ = fem.matrici_isolate(mod_psi, modello_massa,  # type: ignore[arg-type]
                                    J_aggiunta=pl.J_ax + float(pl.J_psi(th)))
    Mi_th, _ = fem.matrici_isolate(mod_th, modello_massa, J_aggiunta=pl.I_yy)  # type: ignore[arg-type]
    beta_psi = fem.beta_smorzamento(Mi_psi, mod_psi.K, par.zeta)
    beta_th = fem.beta_smorzamento(Mi_th, mod_th.K, par.zeta)
    assert mod_th.i_P is not None
    return (
        AsseStrutturale("azimut", mod_psi, modello_massa, M_psi, mod_psi.K.copy(), mod_psi.C_b.copy(),
                        beta_psi, mod_psi.i_encoder, mod_psi.i_motore, mod_psi.i_estremo,
                        mod_psi.i_estremo),
        AsseStrutturale("polare", mod_th, modello_massa, M_th, mod_th.K.copy(), mod_th.C_b.copy(),
                        beta_th, mod_th.i_encoder, mod_th.i_motore, mod_th.i_P, mod_th.i_estremo),
    )


def inerzia_payload_rif(asse: str, pl: Payload, theta_ref: float) -> float:
    """Inerzia del payload congelata nel nodo d'aggancio (doc. unificato §3).

    Azimut: :math:`J_{\\psi,p}(\\theta_{ref}) = M_{p,11}(\\theta_{ref})`; polare:
    :math:`I_{yy,p} = M_{p,22}`. :math:`J^{ax}` e' gia' nella massa strutturale.
    """
    return float(pl.J_psi(theta_ref)) if asse == "azimut" else float(pl.I_yy)


def massa_rif(ast: AsseStrutturale, pl: Payload, theta_ref: float) -> np.ndarray:
    """:math:`\\bar{\\mathbf M}_a`: massa strutturale con il payload congelato (k_g = 0)."""
    M = ast.M.copy()
    M[ast.i_agg, ast.i_agg] += inerzia_payload_rif(ast.nome, pl, theta_ref)
    return M


# ---------------------------------------------------------------------------
# Analisi modale di un asse (doc. unificato §4, Gawronski §2.2 e §4.3)
# ---------------------------------------------------------------------------


def modi_con_rigido(M: np.ndarray, K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """:math:`(\\omega^2, \\boldsymbol\\Phi)` con il modo rigido :math:`\\mathbf 1` NON normalizzato.

    ``eigh(K, M)`` (LAPACK dsygvd) normalizza in massa; la prima colonna, il modo
    rigido (:math:`\\mathbf K\\mathbf 1 = 0`, Gawronski 2.30), si sostituisce con
    :math:`\\mathbf 1`, cosi' :math:`\\eta_0` e' un angolo e la sua massa modale e'
    :math:`J_{tot} = \\mathbf 1^{\\mathsf T}\\bar{\\mathbf M}\\mathbf 1`. I flessibili hanno
    segno tale che la componente piu' grande sia positiva.
    """
    w2, Phi = sla.eigh(K, M)
    scala = float(np.max(np.abs(w2)))
    nulli = np.abs(w2) < 1e-9 * scala
    if int(nulli.sum()) != 1 or not nulli[0]:
        raise AssertionError(f"atteso un solo modo rigido, trovati {int(nulli.sum())}")
    v = Phi[:, 0] / Phi[0, 0]
    if np.max(np.abs(v - 1.0)) > 1e-6:
        raise AssertionError("il modo a frequenza nulla non e' il moto rigido 1")
    w2 = w2.copy()
    w2[0] = 0.0
    Phi = Phi.copy()
    Phi[:, 0] = 1.0
    for k in range(1, Phi.shape[1]):
        if Phi[np.argmax(np.abs(Phi[:, k])), k] < 0:
            Phi[:, k] *= -1.0
    return w2, Phi


def norme_hankel(Phi: np.ndarray, w2: np.ndarray, beta: float, i_m: int, i_e: int) -> np.ndarray:
    """Norma di Hankel di ogni modo per l'uscita di posizione encoder (Gawronski 4.4-4.6).

    :math:`h_i = |\\phi_i(m)\\,\\phi_i(e)|/(4\\zeta_i\\omega_i^2)`, con :math:`\\zeta_i =
    \\beta\\omega_i/2` e :math:`\\boldsymbol\\Phi` normalizzata in massa. Modo rigido:
    :math:`+\\infty` (si tiene sempre, Gawronski 4.11-4.12). Senza smorzamento
    (:math:`\\beta = 0`) tutte le norme sono infinite.
    """
    h = np.full(Phi.shape[1], np.inf)
    om = np.sqrt(w2[1:])
    z = 0.5 * beta * om
    with np.errstate(divide="ignore"):
        h[1:] = np.abs(Phi[i_m, 1:] * Phi[i_e, 1:]) / (4.0 * z * om**2)
    return h


@dataclass(frozen=True)
class AnalisiModaleAsse:
    """Tutti i modi del FEM completo di un asse con il payload congelato in ``theta_ref``.

    Le colonne di ``Phi`` sono il modo rigido (colonna 0, :math:`\\mathbf 1`) e i modi
    flessibili in ordine di frequenza (colonna k = modo flessibile k).
    """

    asse: str
    theta_ref: float
    w2: np.ndarray
    Phi: np.ndarray
    M_bar: np.ndarray
    beta: float
    J_tot: float
    hankel: np.ndarray
    tenuti: tuple[int, ...]
    """Indici (colonne di ``Phi``, >= 1) dei modi flessibili tenuti."""
    motivo: dict[int, str]
    """Per ogni modo tenuto: ``taglio``, ``hankel``, ``forzato`` o ``validazione``."""
    i_e: int
    i_m: int
    i_agg: int

    @property
    def omega(self) -> np.ndarray:
        return np.sqrt(np.maximum(self.w2, 0.0))

    @property
    def f(self) -> np.ndarray:
        return self.omega / (2 * np.pi)

    @property
    def zeta(self) -> np.ndarray:
        return 0.5 * self.beta * self.omega

    @property
    def masse_modali(self) -> np.ndarray:
        return np.einsum("ik,ij,jk->k", self.Phi, self.M_bar, self.Phi)

    @property
    def c_e(self) -> np.ndarray:
        """Riga di :math:`\\boldsymbol\\Phi` all'encoder (Gawronski 2.22)."""
        return self.Phi[self.i_e]

    @property
    def c_m(self) -> np.ndarray:
        """Riga al nodo motore: partecipazione modale della coppia (Gawronski 2.21)."""
        return self.Phi[self.i_m]

    @property
    def c_agg(self) -> np.ndarray:
        """Riga al nodo d'aggancio del payload: :math:`\\mathbf c_L` o :math:`\\mathbf c_P`."""
        return self.Phi[self.i_agg]

    @property
    def scartati(self) -> tuple[int, ...]:
        return tuple(k for k in range(1, self.Phi.shape[1]) if k not in self.tenuti)

    def flessibilita_residua(self, i_out: int, i_in: int) -> float:
        """:math:`\\sum_{i\\ \\text{scartati}}\\phi_i(out)\\phi_i(in)/\\omega_i^2` (rad/(N m))."""
        return float(sum(self.Phi[i_out, k] * self.Phi[i_in, k] / self.w2[k] for k in self.scartati))


def seleziona_modi(f: np.ndarray, hankel: np.ndarray, pm: ParametriModali, asse: str,
                   extra: Sequence[int] = ()) -> tuple[tuple[int, ...], dict[int, str]]:
    """Modi flessibili tenuti: taglio in frequenza, norma di Hankel, numero forzato, extra."""
    n = len(f)
    motivo: dict[int, str] = {}
    forzato = pm.n_flessibili.get(asse)
    if forzato is not None:
        for k in range(1, min(forzato, n - 1) + 1):
            motivo[k] = "forzato"
    else:
        for k in range(1, n):
            if f[k] <= pm.f_taglio:
                motivo[k] = "taglio"
        if pm.hankel_min is not None and n > 1:
            h_max = float(np.max(hankel[1:]))
            for k in range(1, n):
                if k not in motivo and hankel[k] >= pm.hankel_min * h_max:
                    motivo[k] = "hankel"
    for k in extra:
        motivo.setdefault(int(k), "validazione")
    return tuple(sorted(motivo)), motivo


def analizza_asse(ast: AsseStrutturale, pl: Payload, theta_ref: float, pm: ParametriModali,
                  extra: Sequence[int] = ()) -> AnalisiModaleAsse:
    """Base modale completa dell'asse in ``theta_ref`` e scelta dei modi."""
    Mb = massa_rif(ast, pl, theta_ref)
    w2, Phi = modi_con_rigido(Mb, ast.K)
    h = norme_hankel(Phi, w2, ast.beta, ast.i_m, ast.i_e)
    f = np.sqrt(w2) / (2 * np.pi)
    tenuti, motivo = seleziona_modi(f, h, pm, ast.nome, extra)
    J_tot = float(np.ones(ast.n) @ Mb @ np.ones(ast.n))
    return AnalisiModaleAsse(asse=ast.nome, theta_ref=float(theta_ref), w2=w2, Phi=Phi, M_bar=Mb,
                             beta=ast.beta, J_tot=J_tot, hankel=h, tenuti=tenuti, motivo=motivo,
                             i_e=ast.i_e, i_m=ast.i_m, i_agg=ast.i_agg)


# ---------------------------------------------------------------------------
# Basi di Ritz
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BaseAsse:
    """Base di Ritz di un asse, :math:`\\mathbf q_{FEM} = \\mathbf T\\boldsymbol\\eta`.

    ``M_bar_p`` e' l'inerzia del payload contenuta in :math:`\\tilde{\\mathbf M}` (solo
    per le basi modali). ``omega`` e ``zeta`` sono definiti per le basi modali
    (:math:`\\omega_0 = 0` per il modo rigido) e servono allo stato a blocchi di
    Gawronski. ``spostamenti`` (solo correzione statica): spostamento relativo
    delle frequenze dei modi tenuti.
    """

    asse: str
    tipo: TipoBase
    T: np.ndarray
    nomi: tuple[str, ...]
    M_bar_p: float = 0.0
    omega: np.ndarray | None = None
    zeta: np.ndarray | None = None
    theta_ref: float | None = None
    modi_fem: tuple[int, ...] = ()
    """Per le basi modali: modi del FEM completo tenuti (0 = rigido)."""
    n_residui: int = 0
    spostamenti: np.ndarray | None = None

    @property
    def r(self) -> int:
        return self.T.shape[1]


def _simbolo(asse: str) -> str:
    return "psi" if asse == "azimut" else "theta"


def base_completa(ast: AsseStrutturale) -> BaseAsse:
    s = _simbolo(ast.nome)
    return BaseAsse(ast.nome, "completa", np.eye(ast.n),
                    tuple(f"{s}_{ast.etichetta(i)}" for i in range(ast.n)))


def base_guyan(ast: AsseStrutturale) -> BaseAsse:
    mod = ast.modello
    T = rd.trasformazione_guyan(ast.K, mod.masters)
    s = _simbolo(ast.nome)
    return BaseAsse(ast.nome, "guyan", T, tuple(f"{s}_{r}" for r in rd._ruoli(mod)))


def base_craig_bampton(ast: AsseStrutturale, n_modi: int) -> BaseAsse:
    """Guyan piu' ``n_modi`` modi a interfaccia fissa (sezione 8.3)."""
    mod = ast.modello
    T, _, _ = rd.riduci_craig_bampton(ast.M, ast.K, mod.masters, n_modi)
    s = _simbolo(ast.nome)
    nomi = tuple(f"{s}_{r}" for r in rd._ruoli(mod)) + tuple(
        f"{s}_cb{k + 1}" for k in range(T.shape[1] - len(mod.masters)))
    return BaseAsse(ast.nome, "craig_bampton", T, nomi)


def base_modale(an: AnalisiModaleAsse, pl: Payload, ast: AsseStrutturale,
                correzione: str = "none") -> BaseAsse:
    """Modo rigido piu' modi flessibili tenuti; opzionalmente vettori statici residui.

    ``residual_vectors``: per ogni nodo caricato (motore e aggancio del payload) il
    vettore :math:`\\mathbf r_j = \\sum_{i\\ \\text{scartati}}\\boldsymbol\\phi_i\\phi_i(j)/\\omega_i^2`,
    flessibilita' residua con "inertia relief" (l'unica definita per una catena
    libera). E' gia' :math:`\\bar{\\mathbf M}`- e :math:`\\mathbf K`-ortogonale ai modi
    tenuti; si ortogonalizza comunque numericamente, si scartano le colonne
    dipendenti e si risolve il problema agli autovalori nello spazio aumentato:
    la base resta modale (diagonale) e il modello resta un'ODE valida anche nel
    non lineare.
    """
    colonne = [0, *an.tenuti]
    T = an.Phi[:, colonne]
    s = _simbolo(an.asse)
    nomi = [f"eta_{s}{k}" for k in colonne]
    w2 = an.w2[colonne]
    n_res = 0
    spost = None
    if correzione == "residual_vectors" and an.scartati:
        Mb = an.M_bar
        R = np.column_stack([sum(an.Phi[:, k] * an.Phi[j, k] / an.w2[k] for k in an.scartati)
                             for j in sorted({an.i_m, an.i_agg})])
        # Gram-Schmidt in metrica M rispetto ai modi tenuti e fra i residui
        Q = T.copy()
        for c in R.T:
            v = c.copy()
            for _ in range(2):
                for q in Q.T:
                    v = v - (q @ Mb @ v) / (q @ Mb @ q) * q
            nv = float(np.sqrt(v @ Mb @ v))
            if nv > 1e-8 * float(np.sqrt(c @ Mb @ c)):
                Q = np.column_stack([Q, v / nv])
                n_res += 1
        Kr, Mr = Q.T @ ast.K @ Q, Q.T @ Mb @ Q
        w2n, V = sla.eigh(0.5 * (Kr + Kr.T), 0.5 * (Mr + Mr.T))
        Tn = Q @ V
        w2n_rid, Tn = _rigido_in_testa(w2n, Tn)
        spost = np.sqrt(w2n_rid[1:len(colonne)] / w2[1:]) - 1.0
        T, w2 = Tn, w2n_rid
        # modi tenuti corretti ("c") e pseudo-modi dei vettori residui ("r")
        nomi = ([f"eta_{s}0"] + [f"eta_{s}{k}c" for k in colonne[1:]]
                + [f"eta_{s}_r{k + 1}" for k in range(n_res)])
    om = np.sqrt(np.maximum(w2, 0.0))
    return BaseAsse(an.asse, "modale", T, tuple(nomi),
                    M_bar_p=inerzia_payload_rif(an.asse, pl, an.theta_ref), omega=om,
                    zeta=0.5 * an.beta * om, theta_ref=an.theta_ref, modi_fem=tuple(colonne),
                    n_residui=n_res, spostamenti=spost)


def _rigido_in_testa(w2: np.ndarray, V: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rimette il modo rigido come :math:`\\mathbf 1` non normalizzato in prima colonna."""
    scala = float(np.max(np.abs(w2)))
    w2 = w2.copy()
    if abs(w2[0]) > 1e-9 * scala:
        raise AssertionError("base aumentata: manca il modo rigido")
    v = V[:, 0] / V[0, 0]
    if np.max(np.abs(v - 1.0)) > 1e-6:
        raise AssertionError("base aumentata: il primo modo non e' il moto rigido")
    V = V.copy()
    V[:, 0] = 1.0
    w2[0] = 0.0
    for k in range(1, V.shape[1]):
        if V[np.argmax(np.abs(V[:, k])), k] < 0:
            V[:, k] *= -1.0
    return w2, V


# ---------------------------------------------------------------------------
# Modello accoppiato su base generica (doc. unificato §5.1)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StrutturaRitz:
    """Dimensioni e nomi delle coordinate (analogo di :class:`gimbal_mdof.coupled.Struttura`)."""

    n_psi: int
    n_theta: int
    nomi: tuple[str, ...]

    @property
    def n(self) -> int:
        return self.n_psi + self.n_theta


@dataclass(frozen=True)
class ModelloRitz:
    """Modello accoppiato non lineare su una base di Ritz dei due assi.

    Offre la stessa interfaccia numerica di :class:`gimbal_mdof.coupled.ModelloAccoppiato`
    usata da :mod:`gimbal_mdof.simulate` e :mod:`gimbal_mdof.linearize`.
    """

    asse_psi: AsseStrutturale
    asse_theta: AsseStrutturale
    base_psi: BaseAsse
    base_theta: BaseAsse
    payload: Payload
    M_tilde: np.ndarray
    K_tilde: np.ndarray
    D_tilde: np.ndarray
    E_eta: np.ndarray
    """(2 x n): :math:`[\\psi_L, \\theta_p]^{\\mathsf T} = \\mathbf E_\\eta\\boldsymbol\\eta`."""
    B: np.ndarray
    """(n x 2): :math:`\\mathbf T^{\\mathsf T}\\mathbf B_o`, coppie ai nodi motore."""
    M_bar_p: np.ndarray
    p_payload: np.ndarray
    segnaposto: bool
    funzioni: FunzioniPayload

    # -- dimensioni ------------------------------------------------------------
    @property
    def n(self) -> int:
        return self.M_tilde.shape[0]

    @property
    def struttura(self) -> StrutturaRitz:
        return StrutturaRitz(self.base_psi.r, self.base_theta.r,
                             self.base_psi.nomi + self.base_theta.nomi)

    @property
    def tipo(self) -> str:
        return self.base_psi.tipo

    # -- matrici --------------------------------------------------------------
    def matrice_ingressi(self) -> np.ndarray:
        return self.B

    def matrice_uscite(self) -> np.ndarray:
        """Righe :data:`gimbal_mdof.linearize.USCITE` in funzione di :math:`\\boldsymbol\\eta`."""
        bp, bt, ap, at = self.base_psi, self.base_theta, self.asse_psi, self.asse_theta
        n_psi = bp.r

        def riga(asse: str, i: int) -> np.ndarray:
            c = np.zeros(self.n)
            if asse == "azimut":
                c[:n_psi] = bp.T[i]
            else:
                c[n_psi:] = bt.T[i]
            return c

        r = {"psi_e": riga("azimut", ap.i_e), "theta_e": riga("polare", at.i_e),
             "psi_m": riga("azimut", ap.i_m), "theta_m": riga("polare", at.i_m),
             "psi_L": riga("azimut", ap.i_agg), "theta_p": riga("polare", at.i_agg),
             "theta_r": riga("polare", at.i_r)}
        r["err_psi"] = r["psi_e"] - r["psi_L"]
        r["err_theta"] = r["theta_e"] - r["theta_p"]
        return np.array([r[u] for u in lz.USCITE])

    def theta_p(self, q: np.ndarray) -> float:
        return float(self.E_eta[1] @ q)

    def M(self, q: np.ndarray) -> np.ndarray:
        """:math:`\\tilde{\\mathbf M} + \\mathbf E_\\eta^{\\mathsf T}[\\mathbf M_p(\\theta_p) - \\bar{\\mathbf M}_p]\\mathbf E_\\eta`."""
        Mp = np.asarray(self.funzioni.M(self.theta_p(q), self.p_payload), dtype=float)
        return self.M_tilde + self.E_eta.T @ (Mp - self.M_bar_p) @ self.E_eta

    def forze(self, q: np.ndarray, dq: np.ndarray) -> np.ndarray:
        """:math:`\\tilde{\\mathbf D}\\dot{\\boldsymbol\\eta} + \\tilde{\\mathbf K}\\boldsymbol\\eta
        + \\mathbf E_\\eta^{\\mathsf T}(\\mathbf C_p\\dot{\\mathbf q}_p + \\mathbf g_p)`."""
        dqp = self.E_eta @ dq
        fp = np.asarray(self.funzioni.forze(self.theta_p(q), dqp[0], dqp[1], self.p_payload),
                        dtype=float).reshape(2)
        return self.D_tilde @ dq + self.K_tilde @ q + self.E_eta.T @ fp

    def residuo(self, q: np.ndarray, dq: np.ndarray, ddq: np.ndarray, tau: np.ndarray) -> np.ndarray:
        """Residuo della (5.1): :math:`\\mathbf M\\ddot{\\boldsymbol\\eta} + \\ldots - \\mathbf B\\boldsymbol\\tau`."""
        return self.M(q) @ ddq + self.forze(q, dq) - self.B @ tau

    def accelerazioni(self, q: np.ndarray, dq: np.ndarray, tau: np.ndarray) -> np.ndarray:
        return np.linalg.solve(self.M(q), self.B @ tau - self.forze(q, dq))

    def energia(self, q: np.ndarray, dq: np.ndarray) -> float:
        """:math:`E = \\tfrac12\\dot{\\boldsymbol\\eta}^{\\mathsf T}\\mathbf M\\dot{\\boldsymbol\\eta}
        + \\tfrac12\\boldsymbol\\eta^{\\mathsf T}\\tilde{\\mathbf K}\\boldsymbol\\eta + V_g(\\theta_p)`."""
        return float(0.5 * dq @ self.M(q) @ dq + 0.5 * q @ self.K_tilde @ q
                     + self.funzioni.V(self.theta_p(q), self.p_payload))

    def direzione_rigida(self, asse: str) -> np.ndarray:
        """:math:`\\mathbf u` con :math:`\\mathbf T_a\\mathbf u = \\mathbf 1` (moto rigido dell'asse)."""
        b = self.base_psi if asse == "azimut" else self.base_theta
        u_a, *_ = np.linalg.lstsq(b.T, np.ones(b.T.shape[0]), rcond=None)
        if np.max(np.abs(b.T @ u_a - 1.0)) > 1e-8:
            raise AssertionError(f"il moto rigido dell'asse {asse} non e' nello span della base")
        u = np.zeros(self.n)
        if asse == "azimut":
            u[: self.base_psi.r] = u_a
        else:
            u[self.base_psi.r:] = u_a
        return u

    def momento_psi(self, q: np.ndarray, dq: np.ndarray) -> float:
        """:math:`p_\\psi = \\mathbf u_\\psi^{\\mathsf T}\\mathbf M\\dot{\\boldsymbol\\eta}`: la
        traslazione rigida di azimut e' ciclica (:math:`\\mathbf K\\mathbf 1 = 0`, payload
        indipendente da :math:`\\psi_L`)."""
        return float(self.direzione_rigida("azimut") @ self.M(q) @ dq)

    def senza_dissipazione(self) -> "ModelloRitz":
        return replace(self, D_tilde=np.zeros_like(self.D_tilde))

    # -- punto di lavoro e linearizzazione (doc. unificato §7) ------------------
    def equilibrio(self, theta0: float, Omega0: float) -> lz.PuntoDiLavoro:
        """Moto stazionario :math:`\\theta_p = \\theta_0`, :math:`\\dot\\psi = \\Omega_0`, con
        :math:`\\psi_e = 0`. Con :math:`\\theta_p` fissato il residuo e' affine: sistema lineare."""
        n = self.n
        u = self.direzione_rigida("azimut")
        dq0 = Omega0 * u
        fp = np.asarray(self.funzioni.forze(theta0, Omega0 * (self.E_eta[0] @ u), 0.0, self.p_payload),
                        dtype=float).reshape(2)
        r0 = self.D_tilde @ dq0 + self.E_eta.T @ fp
        C = self.matrice_uscite()
        A = np.zeros((n + 2, n + 2))
        A[:n, :n] = self.K_tilde
        A[:n, n:] = -self.B
        A[n, :n] = self.E_eta[1]
        A[n + 1, :n] = C[lz.USCITE.index("psi_e")]
        b = np.concatenate([-r0, [theta0, 0.0]])
        z = np.linalg.solve(A, b)
        q0, tau0 = z[:n], z[n:]
        res = self.forze(q0, dq0) - self.B @ tau0
        scala = max(float(np.max(np.abs(r0))), float(np.max(np.abs(self.B @ tau0))), 1e-12)
        if np.max(np.abs(res)) > 1e-8 * scala or abs(self.theta_p(q0) - theta0) > 1e-12:
            raise AssertionError(f"equilibrio non soddisfatto: residuo {np.max(np.abs(res)):.3e}")
        return lz.PuntoDiLavoro(theta0=float(theta0), Omega0=float(Omega0), q0=q0, tau0=tau0,
                                residuo=float(np.max(np.abs(res))))

    def matrici_lineari(self, theta0: float, Omega0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """:math:`(\\mathbf M_0, \\mathbf C_0, \\mathbf K_0)` (doc. unificato 7.2).

        :math:`\\mathbf M_0 = \\tilde{\\mathbf M} + \\mathbf E_\\eta^{\\mathsf T}[\\mathbf M_p(\\theta_0) -
        \\bar{\\mathbf M}_p]\\mathbf E_\\eta` (con :math:`\\theta_{ref} = \\theta_0` resta il solo
        :math:`M_{12}`), :math:`\\mathbf C_0 = \\tilde{\\mathbf D} + \\mathbf E_\\eta^{\\mathsf T}\\mathbf G_{p0}\\mathbf E_\\eta`
        (:math:`\\gamma`), :math:`\\mathbf K_0 = \\tilde{\\mathbf K} + (k_c + k_g)\\hat{\\mathbf c}_P\\hat{\\mathbf c}_P^{\\mathsf T}`.
        """
        p = self.p_payload
        M0p = np.asarray(self.funzioni.M0(theta0, Omega0, p), dtype=float)
        G0p = np.asarray(self.funzioni.G0(theta0, Omega0, p), dtype=float)
        K0p = np.asarray(self.funzioni.K0(theta0, Omega0, p), dtype=float)
        E = self.E_eta
        return (self.M_tilde + E.T @ (M0p - self.M_bar_p) @ E, self.D_tilde + E.T @ G0p @ E,
                self.K_tilde + E.T @ K0p @ E)

    def linearizza(self, theta0: float, Omega0: float) -> lz.ModelloLineare:
        M0, C0, K0 = self.matrici_lineari(theta0, Omega0)
        n = self.n
        Minv = np.linalg.inv(M0)
        A = np.block([[np.zeros((n, n)), np.eye(n)], [-Minv @ K0, -Minv @ C0]])
        Bu = np.vstack([np.zeros((n, 2)), Minv @ self.B])
        return lz.ModelloLineare(theta0=float(theta0), Omega0=float(Omega0), M0=M0, C0=C0, K0=K0,
                                 E=self.B, A=A, Bu=Bu, Cq=self.matrice_uscite(),
                                 punto=self.equilibrio(theta0, Omega0), accoppiato=self)  # type: ignore[arg-type]


def costruisci_ritz(ast_psi: AsseStrutturale, ast_th: AsseStrutturale, b_psi: BaseAsse,
                    b_th: BaseAsse, pl: Payload, *, segnaposto: bool = False) -> ModelloRitz:
    """Proietta le matrici strutturali e i vettori d'aggancio sulle basi dei due assi."""
    rp, rt = b_psi.r, b_th.r
    n = rp + rt
    Mt, Kt, Dt = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
    for ast, b, sl in ((ast_psi, b_psi, slice(0, rp)), (ast_th, b_th, slice(rp, n))):
        Mbar = ast.M.copy()
        Mbar[ast.i_agg, ast.i_agg] += b.M_bar_p
        for A, S in ((Mt, Mbar), (Kt, ast.K), (Dt, ast.C)):
            P = b.T.T @ S @ b.T
            A[sl, sl] = 0.5 * (P + P.T)
    E = np.zeros((2, n))
    E[0, :rp] = b_psi.T[ast_psi.i_agg]
    E[1, rp:] = b_th.T[ast_th.i_agg]
    B = np.zeros((n, 2))
    B[:rp, 0] = b_psi.T[ast_psi.i_m]
    B[rp:, 1] = b_th.T[ast_th.i_m]
    return ModelloRitz(asse_psi=ast_psi, asse_theta=ast_th, base_psi=b_psi, base_theta=b_th,
                       payload=pl, M_tilde=Mt, K_tilde=Kt, D_tilde=Dt, E_eta=E, B=B,
                       M_bar_p=np.diag([b_psi.M_bar_p, b_th.M_bar_p]), p_payload=parametri_payload(pl),
                       segnaposto=segnaposto, funzioni=funzioni_payload())


# ---------------------------------------------------------------------------
# Stato a blocchi di Gawronski (eq. 2.25 e 2.32)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StatoGawronski:
    """Spazio degli stati con :math:`\\mathbf x_i = [\\omega_i\\eta_i,\\ \\dot\\eta_i]` (flessibili)
    e :math:`[\\eta_0,\\ \\dot\\eta_0]` (rigidi): blocchi diagonali (2.25)/(2.32) piu' le
    correzioni di rango basso fuori diagonale."""

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    S: np.ndarray
    """:math:`\\mathbf x = \\mathbf S[\\boldsymbol\\eta;\\ \\dot{\\boldsymbol\\eta}]`."""
    nomi: tuple[str, ...]


def stato_gawronski(ml: lz.ModelloLineare) -> StatoGawronski:
    acc = ml.accoppiato
    assert isinstance(acc, ModelloRitz) and acc.base_psi.omega is not None and acc.base_theta.omega is not None
    om = np.concatenate([acc.base_psi.omega, acc.base_theta.omega])
    n = om.size
    S = np.zeros((2 * n, 2 * n))
    nomi = []
    for i, (w, nm) in enumerate(zip(om, acc.struttura.nomi)):
        S[2 * i, i] = w if w > 0 else 1.0
        S[2 * i + 1, n + i] = 1.0
        nomi += [f"{'w*' if w > 0 else ''}{nm}", f"d{nm}"]
    Si = np.linalg.inv(S)
    return StatoGawronski(A=S @ ml.A @ Si, B=S @ ml.Bu, C=ml.Cy @ Si, S=S, nomi=tuple(nomi))


# ---------------------------------------------------------------------------
# Confronto in banda con il modello a tutti i GdL (test T5)
# ---------------------------------------------------------------------------


def zeri_leggeri(ml: lz.ModelloLineare, uscita: str, j: int) -> np.ndarray:
    """Frequenze (Hz) degli zeri poco smorzati (|Im| > |Re|, Im > 0) del canale SISO."""
    z = tf.zeri_invarianti(ml.A, ml.Bu[:, j], ml.Cy[ml.riga(uscita)])
    z = z[(np.imag(z) > 0) & (np.abs(np.imag(z)) > np.abs(np.real(z)))]
    return np.sort(np.abs(z)) / (2 * np.pi)


@dataclass(frozen=True)
class EsitoConfronto:
    """Scarti in banda del modello ridotto dal modello a tutti i GdL in un punto di lavoro."""

    theta0: float
    Omega0: float
    err_db: tuple[float, float]
    """Massimo scarto di modulo (dB) su tau_psi -> psi_e e tau_theta -> theta_e, fuori dalle bande degli zeri."""
    err_deg: tuple[float, float]
    zeri: tuple[tuple[str, float, float, float], ...]
    """(canale, f_riferimento, f_ridotto, scarto relativo) per ogni zero in banda."""
    err_incrociato: float
    """Massimo di :math:`\\|\\mathbf G - \\mathbf G_{ref}\\|_2/\\|\\mathbf G_{ref}\\|_2` in banda."""

    @property
    def passa(self) -> bool:
        ok_zeri = all(np.isfinite(e) and e <= TOLL_ZERO for *_, e in self.zeri)
        return (max(self.err_db) < SOGLIA_DB and max(self.err_deg) < SOGLIA_DEG and ok_zeri
                and self.err_incrociato <= SOGLIA_INCROCIATI)


def confronta(ml: lz.ModelloLineare, ml_ref: lz.ModelloLineare, f_min: float, f_banda: float,
              n: int = N_GRIGLIA_VALIDAZIONE) -> EsitoConfronto:
    """Criteri T5: canali diretti entro 1 dB e 5 gradi fuori da +/-2 % attorno agli zeri di
    riferimento, zeri esclusi coincidenti entro 1 %, matrice 2x2 entro il 12 %."""
    f = np.logspace(np.log10(f_min), np.log10(f_banda), n)
    w = 2 * np.pi * f
    G = ml.frf(w, ["psi_e", "theta_e"])
    Gr = ml_ref.frf(w, ["psi_e", "theta_e"])
    escludi = np.zeros(f.size, dtype=bool)
    err_db, err_deg, zeri = [], [], []
    for j, (u, nome) in enumerate((("psi_e", "tau_psi -> psi_e"), ("theta_e", "tau_theta -> theta_e"))):
        zr = zeri_leggeri(ml_ref, u, j)
        zr = zr[(zr >= f_min) & (zr <= f_banda)]
        zm = zeri_leggeri(ml, u, j)
        maschera = np.zeros(f.size, dtype=bool)
        for fz in zr:
            maschera |= np.abs(f / fz - 1.0) <= BANDA_ZERO
            if zm.size:
                k = int(np.argmin(np.abs(zm - fz)))
                zeri.append((nome, float(fz), float(zm[k]), float(abs(zm[k] / fz - 1.0))))
            else:
                zeri.append((nome, float(fz), float("nan"), float("inf")))
        escludi |= maschera
        r = G[:, j, j] / Gr[:, j, j]
        ok = ~maschera
        err_db.append(float(np.max(np.abs(20 * np.log10(np.abs(r[ok]))))) if ok.any() else 0.0)
        err_deg.append(float(np.max(np.abs(np.degrees(np.angle(r[ok]))))) if ok.any() else 0.0)
    ok = ~escludi
    num = np.linalg.norm(G[ok] - Gr[ok], ord=2, axis=(1, 2))
    den = np.linalg.norm(Gr[ok], ord=2, axis=(1, 2))
    return EsitoConfronto(theta0=float(ml.theta0), Omega0=float(ml.Omega0),
                          err_db=(err_db[0], err_db[1]), err_deg=(err_deg[0], err_deg[1]),
                          zeri=tuple(zeri), err_incrociato=float(np.max(num / den)))


# ---------------------------------------------------------------------------
# Zeri dell'asse isolato e flessibilita' residua (dal prototipo confronto_guyan_modale.py)
# ---------------------------------------------------------------------------


def zeri_da_modi(x: np.ndarray, r: np.ndarray, D: float = 0.0, n_scan: int = 4001) -> np.ndarray:
    """Zeri non smorzati (Hz) di :math:`g(x) = \\sum_i r_i/(x_i - x) + D`, :math:`x = \\omega^2`.

    ``x[0] = 0`` e' il modo rigido; ``D`` e' il termine statico della correzione dei
    modi scartati. Un cambio di segno per intervallo fra poli consecutivi, raffinato
    con ``brentq`` (funzione ripresa da ``tests/confronto_guyan_modale.py``).
    """
    from scipy.optimize import brentq

    ordine = np.argsort(x)
    x, r = np.asarray(x)[ordine], np.asarray(r)[ordine]

    def g(xx: float) -> float:
        return float(np.sum(r / (x - xx)) + D)

    bordi = list(x) + [10.0 * x[-1] if x[-1] > 0 else 1.0]
    t = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0.0, 1.0, n_scan)[1:-1])
    out = []
    for a, b in zip(bordi[:-1], bordi[1:]):
        pts = a + (b - a) * t
        val = np.array([g(p) for p in pts])
        for k in np.nonzero(np.sign(val[:-1]) * np.sign(val[1:]) < 0)[0]:
            out.append(np.sqrt(brentq(g, pts[k], pts[k + 1], xtol=1e-14 * max(pts[k], 1e-30),
                                      rtol=1e-13)) / (2 * np.pi))
    return np.array(sorted(out))


@dataclass(frozen=True)
class ResiduoAsse:
    """Flessibilita' residua dei modi scartati e zeri che sposta (asse isolato, k_g = 0)."""

    asse: str
    F_collocata: float
    """:math:`\\sum_{scartati}\\phi_i(m)^2/\\omega_i^2` (rad/(N m))."""
    F_non_collocata: float
    """:math:`\\sum_{scartati}\\phi_i(e)\\phi_i(m)/\\omega_i^2`."""
    F_collocata_tot: float
    F_non_collocata_tot: float
    zeri: tuple[tuple[str, float, float, float], ...]
    """(FdT, f completo, f troncato puro, f troncato + correzione statica), primi zeri in banda."""


def residuo_asse(an: AnalisiModaleAsse, f_banda: float) -> ResiduoAsse:
    tutti = range(1, an.Phi.shape[1])
    zeri = []
    mm = an.masse_modali  # J_tot per il rigido (colonna 1 non normalizzata), 1 per i flessibili
    for nome, i_out in (("tau -> encoder", an.i_e), ("tau -> motore", an.i_m)):
        res = an.Phi[i_out] * an.Phi[an.i_m] / mm
        tenuti = [0, *an.tenuti]
        D = an.flessibilita_residua(i_out, an.i_m)
        z_full = zeri_da_modi(an.w2, res)
        z_mod = zeri_da_modi(an.w2[tenuti], res[tenuti])
        z_cor = zeri_da_modi(an.w2[tenuti], res[tenuti], D)
        for k, fz in enumerate(z_full):
            if fz > f_banda:
                break
            zeri.append((nome, float(fz), float(z_mod[k]) if k < z_mod.size else float("nan"),
                         float(z_cor[k]) if k < z_cor.size else float("nan")))
    F = lambda o, i, modi: float(sum(an.Phi[o, k] * an.Phi[i, k] / an.w2[k] for k in modi))  # noqa: E731
    return ResiduoAsse(asse=an.asse, F_collocata=an.flessibilita_residua(an.i_m, an.i_m),
                       F_non_collocata=an.flessibilita_residua(an.i_e, an.i_m),
                       F_collocata_tot=F(an.i_m, an.i_m, tutti), F_non_collocata_tot=F(an.i_e, an.i_m, tutti),
                       zeri=tuple(zeri))


# ---------------------------------------------------------------------------
# Sistema modale completo: scelta dei modi, validazione, controllo di maglia
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ControlloMaglia:
    asse: str
    modo: int
    f_consistente: float
    f_concentrata: float

    @property
    def scarto(self) -> float:
        return self.f_concentrata / self.f_consistente - 1.0

    @property
    def avviso(self) -> bool:
        return abs(self.scarto) > SOGLIA_MAGLIA


@dataclass(frozen=True)
class SistemaModale:
    """Modello modale dei due assi, con il riferimento a tutti i GdL per il confronto."""

    par: Parametri
    modello_massa: str
    payload: Payload
    assi_fem: tuple[fem.ModelloAsse, fem.ModelloAsse]
    ast_psi: AsseStrutturale
    ast_theta: AsseStrutturale
    pm: ParametriModali
    extra: dict[str, tuple[int, ...]] = field(default_factory=lambda: {"azimut": (), "polare": ()})
    validazione: tuple[EsitoConfronto, ...] = ()
    n_aggiunti: int = 0
    """Modi aggiunti dal ciclo di validazione (in ordine decrescente di norma di Hankel)."""
    storia_aggiunte: tuple[tuple[str, int, float], ...] = ()
    """(asse, modo, f [Hz]) dei modi aggiunti dal ciclo di validazione."""

    @property
    def segnaposto(self) -> bool:
        return self.par.segnaposto

    def asse(self, nome: str) -> AsseStrutturale:
        return self.ast_psi if nome == "azimut" else self.ast_theta

    def analisi(self, asse: str, theta_ref: float | None = None) -> AnalisiModaleAsse:
        th = self.pm.theta_ref if theta_ref is None else theta_ref
        return analizza_asse(self.asse(asse), self.payload, th, self.pm, self.extra.get(asse, ()))

    def base(self, asse: str, theta_ref: float | None = None) -> BaseAsse:
        an = self.analisi(asse, theta_ref)
        # i modi tenuti si scelgono in theta_ref di config.yaml e restano gli stessi (per
        # indice) a ogni theta_ref: le frequenze dipendono poco da theta (payload congelato)
        an_cfg = self.analisi(asse)
        an = replace(an, tenuti=an_cfg.tenuti, motivo=an_cfg.motivo)
        return base_modale(an, self.payload, self.asse(asse), self.pm.correzione_statica)

    def modello(self, theta_ref: float | None = None) -> ModelloRitz:
        """Modello modale non lineare (un solo :math:`\\theta_{ref}`: quello di config.yaml)."""
        return costruisci_ritz(self.ast_psi, self.ast_theta, self.base("azimut", theta_ref),
                               self.base("polare", theta_ref), self.payload, segnaposto=self.segnaposto)

    def completo(self) -> ModelloRitz:
        """Modello di riferimento a tutti i GdL del FEM (9 + 10), senza riduzione."""
        return costruisci_ritz(self.ast_psi, self.ast_theta, base_completa(self.ast_psi),
                               base_completa(self.ast_theta), self.payload, segnaposto=self.segnaposto)

    def guyan(self) -> ModelloRitz:
        return costruisci_ritz(self.ast_psi, self.ast_theta, base_guyan(self.ast_psi),
                               base_guyan(self.ast_theta), self.payload, segnaposto=self.segnaposto)

    def craig_bampton(self, n_modi: int) -> ModelloRitz:
        return costruisci_ritz(self.ast_psi, self.ast_theta, base_craig_bampton(self.ast_psi, n_modi),
                               base_craig_bampton(self.ast_theta, n_modi), self.payload,
                               segnaposto=self.segnaposto)

    def lineare(self, theta0: float, Omega0: float) -> lz.ModelloLineare:
        """Modello lineare modale con :math:`\\theta_{ref} = \\theta_0` (doc. unificato 7.2)."""
        return self.modello(theta0).linearizza(theta0, Omega0)

    def lineare_completo(self, theta0: float, Omega0: float) -> lz.ModelloLineare:
        return self.completo().linearizza(theta0, Omega0)

    def confronta(self, theta0: float, Omega0: float) -> EsitoConfronto:
        return confronta(self.lineare(theta0, Omega0), self.lineare_completo(theta0, Omega0),
                         self.par.f_min, self.pm.f_banda)

    def controllo_maglia(self) -> tuple[ControlloMaglia, ...]:
        """Frequenze dei modi tenuti con massa consistente e concentrata (avviso oltre il 2 %)."""
        out = []
        for nome in ("azimut", "polare"):
            an = self.analisi(nome)
            fs = {}
            for mm in ("consistent", "lumped"):
                ast_mm = assi_strutturali(self.par, self.payload, self.assi_fem, mm)[0 if nome == "azimut" else 1]
                fs[mm] = analizza_asse(ast_mm, self.payload, an.theta_ref, self.pm).f
            for k in an.tenuti:
                out.append(ControlloMaglia(nome, k, float(fs["consistent"][k]), float(fs["lumped"][k])))
        return tuple(out)

    def residui(self) -> tuple[ResiduoAsse, ResiduoAsse]:
        return (residuo_asse(self.analisi("azimut"), self.pm.f_banda),
                residuo_asse(self.analisi("polare"), self.pm.f_banda))

    def spettro(self) -> list[tuple[float, str, int, bool]]:
        """Modi dei due assi insieme in :math:`\\theta_{ref}`: rigidi per primi, poi i
        flessibili per frequenza (asse disaccoppiato, :math:`\\Omega_0 = 0`, k_g = 0)."""
        righe = []
        for nome in ("azimut", "polare"):
            an = self.analisi(nome)
            for k in range(an.Phi.shape[1]):
                righe.append((float(an.f[k]), nome, k, k == 0 or k in an.tenuti))
        return sorted(righe, key=lambda r: (r[2] != 0, r[0], r[1]))

    def avvisi(self) -> list[Avviso]:
        out: list[Avviso] = []
        for c in self.controllo_maglia():
            if c.avviso:
                out.append(Avviso(
                    f"MAGLIA_{c.asse.upper()}_{c.modo}", "incoerenza",
                    f"Flexible mode {c.modo} of the {in_inglese(c.asse)} axis: lumped mass "
                    f"{c.f_concentrata:.1f} Hz against consistent {c.f_consistente:.1f} Hz "
                    f"({c.scarto * 100:+.1f} %, threshold {SOGLIA_MAGLIA * 100:g} %): the mesh "
                    "has not converged on this mode."))
        if self.validazione and not all(e.passa for e in self.validazione):
            out.append(Avviso("VALIDAZIONE_MODALE", "incoerenza",
                              "The modal model does not meet the T5 criteria at every operating "
                              "point even with all the available modes."))
        return out


def costruisci_sistema_modale(par: Parametri, modello_massa: str = "consistent", *,
                              payload: Payload | None = None,
                              assi: tuple[fem.ModelloAsse, fem.ModelloAsse] | None = None,
                              punti: Sequence[tuple[float, float]] | None = None,
                              valida: bool = True) -> SistemaModale:
    """Costruisce il sistema modale e, se ``valida``, esegue il ciclo di validazione.

    Ciclo: se in un punto di lavoro i criteri T5 non sono rispettati, si aggiunge il
    modo scartato (fra i due assi) con la norma di Hankel piu' alta e si ripete,
    finche' tutti i punti passano o non restano modi.
    """
    from . import payload as pl_mod

    if assi is None:
        assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
                fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    pl = payload if payload is not None else pl_mod.costruisci_payload(par, assi[1])
    a_psi, a_th = assi_strutturali(par, pl, assi, modello_massa)
    s = SistemaModale(par=par, modello_massa=modello_massa, payload=pl, assi_fem=assi,
                      ast_psi=a_psi, ast_theta=a_th, pm=par.riduzione.modale)
    if not valida:
        return s
    punti = list(punti) if punti is not None else [(t, o) for t in par.theta0 for o in par.omega0]
    storia: list[tuple[str, int, float]] = []
    while True:
        esiti = tuple(s.confronta(t, o) for t, o in punti)
        if all(e.passa for e in esiti):
            break
        candidati = []
        for nome in ("azimut", "polare"):
            an = s.analisi(nome)
            candidati += [(float(an.hankel[k]), nome, k, float(an.f[k])) for k in an.scartati]
        if not candidati:
            break
        h, nome, k, fk = max(candidati)
        storia.append((nome, k, fk))
        extra = dict(s.extra)
        extra[nome] = tuple(sorted(set(extra.get(nome, ())) | {k}))
        s = replace(s, extra=extra)
    return replace(s, validazione=esiti, n_aggiunti=len(storia), storia_aggiunte=tuple(storia))
