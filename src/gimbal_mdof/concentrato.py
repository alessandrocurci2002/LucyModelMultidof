"""Modello a parametri concentrati (3 + 4 GdL) e sue equazioni di Newton-Eulero.

Topologia (richiesta dell'utente, 2026-10-15)::

    azimut:  nodo 1 --[e_psi,1]-- nodo 5 --[e_psi,2]-- nodo 9
    polare:  nodo 1 --[e_th,1]--- nodo P --[e_th,2]--- nodo 7 --[e_th,3]-- nodo 11

Ogni elemento equivalente condensa gli elementi del FEM compresi fra due nodi
principali consecutivi:

.. math::
    J_{eq} = \\sum_i J_i, \\qquad \\frac{1}{k_{eq}} = \\sum_i\\frac{1}{k_i}

(rigidezze in serie). Le eventuali inerzie concentrate nei nodi intermedi si
sommano all'elemento equivalente che le contiene. Il modello e' a parametri
concentrati: ogni nodo principale porta la propria inerzia concentrata piu' meta'
dell'inerzia degli elementi equivalenti adiacenti (``forma="lumped"``); la forma
``"consistent"`` (elemento equivalente uniforme, :math:`\\frac{J}{6}[2\\ 1; 1\\ 2]`) e'
disponibile per confronto.

L'accoppiamento dinamico fra gli assi e' quello del modello gia' sviluppato: il
telescopio e' il corpo rigido agganciato al nodo P, la forcella ruota con il nodo 9
(:math:`\\psi_L`), e le equazioni del moto sono quelle di :mod:`gimbal_mdof.coupled`
(Lagrange con SymPy), di cui questo modulo riusa la catena di calcolo. In piu'
ricava **indipendentemente** le equazioni di Newton-Eulero di ogni GdL (bilancio
dei momenti su ogni nodo; equazioni di Eulero del telescopio attorno al punto
fisso :math:`o_2`) e verifica che coincidano con quelle di Lagrange.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Final

import numpy as np
import sympy as sp

from . import coupled as cp
from . import fem
from . import reduction as rd
from .payload import Payload

FORMA: Final[str] = "concentrato"
"""Etichetta della forma ridotta (``ModelloRidotto.forma_massa``) del modello a parametri concentrati."""


# ---------------------------------------------------------------------------
# Elementi equivalenti
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ElementoEquivalente:
    """Elemento fra due nodi principali consecutivi (etichette in base 1, come nel FEM)."""

    asse: str
    nome: str
    """``1``, ``2`` (azimut) o ``1``, ``2``, ``3`` (polare)."""
    nodo_a: str
    nodo_b: str
    elementi: tuple[str, ...]
    """Etichette degli elementi del FEM condensati."""
    J_elementi: float
    J_nodi_interni: float
    """Inerzie concentrate dei nodi intermedi, sommate all'elemento."""
    k: float
    """Rigidezza in serie, :math:`1/k = \\sum 1/k_i`."""
    ia: int = 0
    ib: int = 1
    """Posizioni dei due nodi fra i master (indici della catena ridotta)."""

    @property
    def J(self) -> float:
        return self.J_elementi + self.J_nodi_interni


def elementi_equivalenti(mod: fem.ModelloAsse) -> tuple[ElementoEquivalente, ...]:
    """Condensa gli elementi del FEM fra nodi principali consecutivi (``mod.masters``).

    Modello ``hub``: i nodi principali della catena sono i master tranne il telescopio P; il ramo
    H-P (molla k_b, senza massa) e' un elemento equivalente a se', l'ultimo. Lo stesso per il ramo L-P
    dell'azimut (``azimuth.telescope_attachment``)."""
    out = []
    masters = list(mod.masters)
    catena = [m for m in masters if not (mod.ramo and m == mod.i_P)]
    for n, (a, b) in enumerate(zip(catena[:-1], catena[1:]), 1):
        gruppo = [(p, (i, j)) for p, (i, j) in zip(mod.elementi, mod.collegamenti) if a <= i and j <= b]
        if not gruppo:
            raise ValueError(f"asse {mod.nome}: nessun elemento fra i nodi {mod.etichetta(a)} e {mod.etichetta(b)}")
        J_int = float(sum(mod.J_nodali[i] for i in range(a + 1, b)))
        out.append(ElementoEquivalente(
            asse=mod.nome, nome=str(n), nodo_a=mod.etichetta(a), nodo_b=mod.etichetta(b),
            elementi=tuple(p.etichetta for p, _ in gruppo),
            J_elementi=float(sum(p.J for p, _ in gruppo)), J_nodi_interni=J_int,
            k=1.0 / float(sum(1.0 / p.k for p, _ in gruppo)), ia=masters.index(a), ib=masters.index(b)))
    if mod.ramo:
        assert mod.i_H is not None and mod.i_P is not None
        gruppo = [p for p, c in zip(mod.elementi, mod.collegamenti) if c == (mod.i_H, mod.i_P)]
        out.append(ElementoEquivalente(
            asse=mod.nome, nome=str(len(out) + 1), nodo_a=mod.etichetta(mod.i_H), nodo_b="P",
            elementi=tuple(p.etichetta for p in gruppo), J_elementi=float(sum(p.J for p in gruppo)),
            J_nodi_interni=0.0, k=1.0 / float(sum(1.0 / p.k for p in gruppo)),
            ia=masters.index(mod.i_H), ib=masters.index(mod.i_P)))
    return tuple(out)


def matrici_concentrate(mod: fem.ModelloAsse, forma: str = "lumped") -> tuple[np.ndarray, np.ndarray]:
    """``(M, K)`` della catena a parametri concentrati sui nodi principali.

    ``lumped``: :math:`M_{jj} = J^n_j + \\tfrac12\\sum_{e\\ni j}J_e`; ``consistent``: elemento
    equivalente uniforme. :math:`K` e' la catena delle molle equivalenti (uguale al
    :math:`\\mathbf K_r` di Guyan, che per una catena e' esattamente in serie).
    """
    ee = elementi_equivalenti(mod)
    n = len(mod.masters)
    M = np.diag([float(mod.J_nodali[i]) for i in mod.masters])
    K = np.zeros((n, n))
    for e in ee:
        ids = np.ix_([e.ia, e.ib], [e.ia, e.ib])
        K[ids] += fem.matrice_rigidezza_elemento(e.k)
        if forma == "lumped":
            M[ids] += 0.5 * e.J * np.eye(2)
        elif forma == "consistent":
            M[ids] += e.J / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
        else:
            raise ValueError(f"forma '{forma}' non ammessa (lumped | consistent)")
    return M, K


def riduci_concentrato(mod: fem.ModelloAsse, forma: str = "lumped") -> rd.ModelloRidotto:
    """Catena a parametri concentrati nella struttura di :class:`gimbal_mdof.reduction.ModelloRidotto`.

    Cosi' :mod:`gimbal_mdof.coupled` (Lagrange con SymPy), :mod:`gimbal_mdof.linearize`,
    :mod:`gimbal_mdof.tf` e :mod:`gimbal_mdof.simulate` si applicano senza modifiche. ``T`` e'
    l'interpolazione statica dei nodi intermedi (pesi sulla cedevolezza), usata solo per
    riportare sui nodi principali l'attrito dei cuscinetti (``C_b``) e per ricostruire le
    deformate; non entra nelle masse.
    """
    M, K = matrici_concentrate(mod, forma)
    T = rd.trasformazione_guyan(mod.K, mod.masters)
    C_b = T.T @ mod.C_b @ T
    return rd.ModelloRidotto(asse=mod.nome, modello_massa="consistent", forma_massa=FORMA,  # type: ignore[arg-type]
                             etichette=mod.etichette_master, ruoli=rd._ruoli(mod), T=T, M_r=M, M_r_guyan=M,
                             K_r=K, C_b_r=0.5 * (C_b + C_b.T), completo=mod)


def costruisci_sistema_concentrato(par, payload: Payload, assi: tuple[fem.ModelloAsse, fem.ModelloAsse],
                                   forma: str = "lumped") -> cp.Sistema:
    """Sistema accoppiato a 7 GdL a parametri concentrati, con lo stesso accoppiamento del modello
    accoppiato (telescopio in P, forcella al nodo 9). :math:`\\beta` e' tarato sul primo modo
    flessibile di ciascun asse concentrato isolato (telescopio e forcella compresi, theta = 0)."""
    mod_psi, mod_th = assi
    rp, rt = riduci_concentrato(mod_psi, forma), riduci_concentrato(mod_th, forma)
    th = cp.THETA_RIFERIMENTO_BETA
    Mp, Kp = rd.matrici_ridotte_isolate(rp, J_aggiunta=payload.J_ax, J_telescopio=float(payload.J_psi(th)))
    Mt, Kt = rd.matrici_ridotte_isolate(rt, J_aggiunta=payload.I_yy)
    beta_psi = fem.beta_smorzamento(Mp, rp.K_r, par.zeta)
    beta_th = fem.beta_smorzamento(Mt, rt.K_r, par.zeta)
    acc = cp.costruisci_accoppiato(rp, rt, payload, beta_psi=beta_psi, beta_theta=beta_th)
    return cp.Sistema(par=par, modello_massa="consistent", forma_massa=FORMA, asse_psi=mod_psi,
                      asse_theta=mod_th, rid_psi=rp, rid_theta=rt, payload=payload, beta_psi=beta_psi,
                      beta_theta=beta_th, accoppiato=acc)


# ---------------------------------------------------------------------------
# Equazioni di Newton-Eulero
# ---------------------------------------------------------------------------

RUOLI_PSI: Final[tuple[str, ...]] = ("e", "m", "L")
RUOLI_PSI_TELESCOPIO: Final[tuple[str, ...]] = ("e", "m", "L", "p")
"""Telescopio appeso al load di azimut (variante ``mozzo_azimut``): encoder, motore, forcella + mozzo +
anello, telescopio."""
RUOLI_THETA: Final[tuple[str, ...]] = ("e", "p", "m", "r")
RUOLI_THETA_MOZZO: Final[tuple[str, ...]] = ("e", "h", "m", "r", "p")
"""Modello ``hub``: encoder, centro del mozzo, motore, estremita', telescopio (appeso a h)."""
ARCHI_THETA: Final[dict[str, tuple[tuple[int, int], ...]]] = {
    "catena": ((0, 1), (1, 2), (2, 3)),
    "mozzo": ((0, 1), (1, 2), (2, 3), (1, 4)),
    "mozzo_azimut": ((0, 1), (1, 2), (2, 3), (1, 4)),
}
"""Molle polari fra i nodi (indici dei ruoli), nell'ordine di k_theta,1, k_theta,2, ..."""
ARCHI_PSI: Final[tuple[tuple[int, int], ...]] = ((0, 1), (1, 2))


def ruoli_theta(variante: str) -> tuple[str, ...]:
    return RUOLI_THETA_MOZZO if variante in ("mozzo", "mozzo_azimut") else RUOLI_THETA


def ruoli_psi(variante: str) -> tuple[str, ...]:
    return RUOLI_PSI_TELESCOPIO if variante == "mozzo_azimut" else RUOLI_PSI


def archi_psi(variante: str) -> tuple[tuple[int, int], ...]:
    """Molle di azimut fra i nodi (indici dei ruoli): catena e - m - L (- p)."""
    n = len(ruoli_psi(variante))
    return tuple((i, i + 1) for i in range(n - 1))


def variante_di(s: cp.Sistema) -> str:
    """``"mozzo"`` se il telescopio e' appeso al mozzo (modello ``hub``), ``"mozzo_azimut"`` se e' appeso
    anche al load di azimut (``azimuth.telescope_attachment``), altrimenti ``"catena"``."""
    if s.rid_psi.completo.ramo:
        return "mozzo_azimut"
    return "mozzo" if s.rid_theta.completo.ramo else "catena"


@dataclass(frozen=True)
class EquazioniNE:
    """Equazioni di Newton-Eulero del modello concentrato, una per GdL.

    ``lhs[i] = rhs[i]``: a sinistra la derivata del momento angolare (inerzia, termini
    giroscopici e centrifughi), a destra i momenti esterni (molle, smorzatori, cuscinetti,
    motori, gravita'). ``residuo = lhs - rhs`` in funzione dei simboli ``q, dq, ddq``.
    """

    q: sp.Matrix
    dq: sp.Matrix
    ddq: sp.Matrix
    nomi: tuple[str, ...]
    lhs: tuple[sp.Expr, ...]
    rhs: tuple[sp.Expr, ...]
    residuo: sp.Matrix
    H_tel: sp.Matrix
    """Momento angolare del telescopio rispetto a o_2, componenti nella terna interna P."""
    H_z: sp.Expr
    """Componente verticale (lungo z_hat) del momento angolare del telescopio."""
    M_g: sp.Matrix
    """Momento della gravita' sul telescopio rispetto a o_2, terna P."""
    J: dict[str, sp.Symbol]
    k: dict[str, sp.Symbol]
    c: dict[str, sp.Symbol]
    b: dict[str, sp.Symbol]
    simboli_parametri: tuple[sp.Symbol, ...]
    variante: str = "catena"

    @property
    def n(self) -> int:
        return len(self.q)

    @property
    def i_theta_p(self) -> int:
        """Indice di theta_p (telescopio) in q."""
        return self.nomi.index("theta_p")


def _simbolo(asse: str, ruolo: str, ordine: int) -> sp.Symbol:
    base = r"\psi" if asse == "psi" else r"\theta"
    pre = {0: "", 1: r"\dot", 2: r"\ddot"}[ordine]
    nome = f"{pre}{{{base}}}_{{{ruolo}}}" if ordine else f"{base}_{{{ruolo}}}"
    return sp.Symbol(nome, real=True)


@lru_cache(maxsize=2)
def newton_eulero(variante: str = "catena") -> EquazioniNE:
    """Deriva con SymPy le equazioni di Newton-Eulero del modello concentrato.

    * Nodi di catena (azimut e, m; polare e, m, r e, nel modello ``hub``, il centro del mozzo h):
      :math:`J_i\\ddot q_i = \\sum` momenti assiali (molle :math:`k(q_a - q_b)`, smorzatori
      :math:`c(\\dot q_a - \\dot q_b)`, cuscinetti :math:`b\\dot q_i`, motore). Per i nodi polari,
      assialsimmetrici attorno a :math:`\\hat J`, la componente assiale dell'equazione di Eulero e'
      :math:`J\\ddot\\theta_i` perche' :math:`\\dot{\\hat J} = \\dot\\psi_L\\hat k\\times\\hat J\\perp\\hat J`.
    * Nodo L (forcella + alberi polari + telescopio), asse fisso :math:`\\hat z`:
      :math:`\\frac{d}{dt}\\bigl[(J_L + J^{ax})\\dot\\psi_L + \\hat z\\cdot\\mathbf H_{tel}\\bigr] = M_z`.
      Variante ``mozzo_azimut``: il telescopio ruota con :math:`\\psi_p`, appeso a L dalla molla
      :math:`k_{\\psi,3}`; su L resta :math:`\\frac{d}{dt}[(J_L + J^{ax})\\dot\\psi_L] = M_z` e su p
      :math:`\\frac{d}{dt}[J_p\\dot\\psi_p + \\hat z\\cdot\\mathbf H_{tel}] = M_z`.
    * Nodo P (telescopio), equazione di Eulero del telescopio rispetto al punto fisso :math:`o_2`
      proiettata su :math:`\\hat J`:
      :math:`J_p\\ddot\\theta_p + \\hat J\\cdot\\bigl(\\dot{\\mathbf H}_{tel}|_P + \\boldsymbol\\omega\\times
      \\mathbf H_{tel}\\bigr) = M_J`.

    ``variante = "catena"``: polare e - p - m - r (7 GdL); ``"mozzo"``: polare e - h - m - r con il
    telescopio p appeso a h dalla molla k_theta,4 (8 GdL); ``"mozzo_azimut"``: come ``"mozzo"`` e in piu'
    azimut e - m - L - p con il telescopio appeso a L dalla molla k_psi,3 (9 GdL).
    """
    rp = ruoli_psi(variante)
    rt = ruoli_theta(variante)
    archi_th = ARCHI_THETA[variante]
    t = sp.Symbol("t", real=True)
    f = {(a, r): sp.Function(f"{a}_{r}", real=True)(t)
         for a, ruoli in (("psi", rp), ("theta", rt)) for r in ruoli}
    J = {f"psi_{r}": sp.Symbol(f"J_{{\\psi,{r}}}", positive=True) for r in rp}
    J |= {f"theta_{r}": sp.Symbol(f"J_{{\\theta,{r}}}", positive=True) for r in rt}
    k = {f"psi_{i}": sp.Symbol(f"k_{{\\psi,{i}}}", positive=True) for i in range(1, len(archi_psi(variante)) + 1)}
    k |= {f"theta_{i}": sp.Symbol(f"k_{{\\theta,{i}}}", positive=True) for i in range(1, len(archi_th) + 1)}
    c = {n: sp.Symbol(s.name.replace("k_", "c_"), nonnegative=True) for n, s in k.items()}
    b = {f"psi_{r}": sp.Symbol(f"b_{{\\psi,{r}}}", nonnegative=True) for r in rp}
    b |= {f"theta_{r}": sp.Symbol(f"b_{{\\theta,{r}}}", nonnegative=True) for r in rt}
    tau_psi, tau_th = cp.tau_psi, cp.tau_theta
    D = lambda x: sp.diff(x, t)  # noqa: E731

    def molla(asse: str, i: int, a: str, b_: str) -> sp.Expr:
        """Momento trasmesso dalla molla i (positivo sul nodo b_): k (q_a - q_b) + c (dq_a - dq_b)."""
        qa, qb = f[(asse, a)], f[(asse, b_)]
        return k[f"{asse}_{i}"] * (qa - qb) + c[f"{asse}_{i}"] * (D(qa) - D(qb))

    def bilancio(asse: str, ruoli: tuple[str, ...], archi) -> dict[str, sp.Expr]:
        """Somma dei momenti delle molle su ciascun nodo (+ sul secondo nodo, - sul primo)."""
        out = {r: sp.S.Zero for r in ruoli}
        for i, (a, b_) in enumerate(archi, 1):
            T = molla(asse, i, ruoli[a], ruoli[b_])
            out[ruoli[a]] -= T
            out[ruoli[b_]] += T
        return out

    psi_e, psi_m, psi_L = (f[("psi", r)] for r in RUOLI_PSI)
    psi_tel = f[("psi", "p")] if "p" in rp else psi_L   # coordinata di azimut del telescopio
    th_p = f[("theta", "p")]
    mom_psi = bilancio("psi", rp, archi_psi(variante))
    mom_th = bilancio("theta", rt, archi_th)

    # --- telescopio: cinematica di coupled.py (terna interna P) --------------------
    omega = cp.velocita_angolare_forma_chiusa(th_p, D(psi_tel), D(th_p))
    I_P = cp.lg.tensore_B()
    H = I_P * omega
    z_P = cp.kin.R_y(-th_p).T * sp.Matrix([0, 0, 1])  # z_hat in terna P
    J_P = sp.Matrix([0, -1, 0])                          # J_hat = -y_P in terna P
    H_z = sp.expand((z_P.T * H)[0, 0])
    d = sp.Matrix([cp.d_x, cp.d_y, cp.d_z])
    M_g = d.cross(-cp.m_p * cp.g * z_P)
    dH_N = H.diff(t) + omega.cross(H)                    # derivata in N, componenti in P

    lhs = [J["psi_e"] * D(D(psi_e)), J["psi_m"] * D(D(psi_m)),
           D((J["psi_L"] + cp.J_ax) * D(psi_L) + (H_z if "p" not in rp else 0))]
    rhs = [mom_psi["e"] - b["psi_e"] * D(psi_e),
           tau_psi + mom_psi["m"] - b["psi_m"] * D(psi_m),
           mom_psi["L"] - b["psi_L"] * D(psi_L)]
    if "p" in rp:
        lhs.append(D(J["psi_p"] * D(psi_tel) + H_z))
        rhs.append(mom_psi["p"] - b["psi_p"] * D(psi_tel))
    for r in rt:
        q_r = f[("theta", r)]
        sx = J[f"theta_{r}"] * D(D(q_r))
        dx = mom_th[r] - b[f"theta_{r}"] * D(q_r)
        if r == "p":
            sx += (J_P.T * dH_N)[0, 0]
            dx += (J_P.T * M_g)[0, 0]
        if r == "m":
            dx += tau_th
        lhs.append(sx)
        rhs.append(dx)

    # --- dalle funzioni del tempo ai simboli q, dq, ddq ---------------------------
    ordine = [("psi", r) for r in rp] + [("theta", r) for r in rt]
    q = sp.Matrix([_simbolo(a, r, 0) for a, r in ordine])
    dq = sp.Matrix([_simbolo(a, r, 1) for a, r in ordine])
    ddq = sp.Matrix([_simbolo(a, r, 2) for a, r in ordine])
    sost = {}
    for (key, s0, s1, s2) in zip(ordine, q, dq, ddq):
        fn = f[key]
        sost[D(D(fn))] = s2
        sost[D(fn)] = s1
    sost_f = {f[key]: s0 for key, s0 in zip(ordine, q)}

    def a_simboli(e: sp.Expr) -> sp.Expr:
        return sp.expand(sp.trigsimp(sp.expand(e.subs(sost).subs(sost_f))))

    lhs_s = tuple(a_simboli(e) for e in lhs)
    rhs_s = tuple(a_simboli(e) for e in rhs)
    residuo = sp.Matrix([sp.expand(a - b_) for a, b_ in zip(lhs_s, rhs_s)])
    nomi = tuple(f"psi_{r}" for r in rp) + tuple(f"theta_{r}" for r in rt)
    parametri = (tuple(J.values()) + tuple(k.values()) + tuple(c.values()) + tuple(b.values())
                 + (cp.J_ax, cp.I_xx, cp.I_yy, cp.I_zz, cp.I_xy, cp.I_xz, cp.I_yz, cp.m_p, cp.g,
                    cp.d_x, cp.d_y, cp.d_z))
    return EquazioniNE(q=q, dq=dq, ddq=ddq, nomi=nomi, lhs=lhs_s, rhs=rhs_s, residuo=residuo,
                       H_tel=H.subs(sost).subs(sost_f), H_z=a_simboli(H_z), M_g=M_g.subs(sost_f),
                       J=J, k=k, c=c, b=b, simboli_parametri=parametri, variante=variante)


def struttura_NE(variante: str = "catena") -> cp.Struttura:
    """Struttura del modello di Lagrange equivalente alla catena concentrata ``variante``."""
    rt = ruoli_theta(variante)
    rp = ruoli_psi(variante)
    return cp.Struttura(ruoli_psi=rp, ruoli_theta=rt, i_psi_e=0, i_psi_m=1, i_psi_L=2,
                        i_theta_e=rt.index("e"), i_theta_p=rt.index("p"), i_theta_m=rt.index("m"),
                        archi_theta=ARCHI_THETA[variante] if variante != "catena" else (),
                        i_psi_p=rp.index("p") if "p" in rp else None)


def verifica_contro_lagrange(variante: str = "catena") -> sp.Matrix:
    """Differenza fra Newton-Eulero e Lagrange (:func:`gimbal_mdof.coupled.costruisci_simbolico`)
    per la catena concentrata: masse diagonali, molle e smorzatori sugli archi, cuscinetti a terra.
    Restituisce il vettore (nullo) delle differenze; solleva ``AssertionError`` se non lo e'."""
    ne = newton_eulero(variante)
    st = struttura_NE(variante)
    mod = cp.costruisci_simbolico(st)
    sost: dict = {}
    for nome, Mm, Km, Cm, ruoli, archi in (("psi", mod.M_psi, mod.K_psi, mod.C_psi, ruoli_psi(variante),
                                            archi_psi(variante)),
                                           ("theta", mod.M_theta, mod.K_theta, mod.C_theta,
                                            ruoli_theta(variante), ARCHI_THETA[variante])):
        n = len(ruoli)
        K_d = [sp.S.Zero] * n
        C_d = [sp.S.Zero] * n
        for i, (a, b_) in enumerate(archi, 1):
            ki, ci = ne.k[f"{nome}_{i}"], ne.c[f"{nome}_{i}"]
            lo, hi = min(a, b_), max(a, b_)
            sost[Mm[lo, hi]] = 0
            sost[Km[lo, hi]] = -ki
            sost[Cm[lo, hi]] = -ci
            K_d[a] += ki
            K_d[b_] += ki
            C_d[a] += ci
            C_d[b_] += ci
        for i in range(n):
            sost[Mm[i, i]] = ne.J[f"{nome}_{ruoli[i]}"]
            sost[Km[i, i]] = K_d[i]
            sost[Cm[i, i]] = C_d[i] + ne.b[f"{nome}_{ruoli[i]}"]
    eom = mod.eom.subs(sost)
    eom = eom.subs(dict(zip(mod.ddq, ne.ddq))).subs(dict(zip(mod.dq, ne.dq))).subs(dict(zip(mod.q, ne.q)))
    diff = (eom - ne.residuo).applyfunc(lambda e: sp.simplify(sp.expand(e)))
    if diff != sp.zeros(len(ne.q), 1):
        raise AssertionError(f"Newton-Eulero e Lagrange non coincidono: {diff.T}")
    return diff


def valori_parametri_NE(s: cp.Sistema) -> dict[sp.Symbol, float]:
    """Valori numerici dei simboli di :func:`newton_eulero` per il sistema concentrato ``s``."""
    variante = variante_di(s)
    ne = newton_eulero(variante)
    out: dict[sp.Symbol, float] = {}
    for nome, rid, beta, ruoli in (("psi", s.rid_psi, s.beta_psi, ruoli_psi(variante)),
                                   ("theta", s.rid_theta, s.beta_theta, ruoli_theta(variante))):
        kc = rid.k_catena
        for i, r in enumerate(ruoli):
            out[ne.J[f"{nome}_{r}"]] = float(rid.M_r[i, i])
            out[ne.b[f"{nome}_{r}"]] = float(rid.C_b_r[i, i])
        for i, kv in enumerate(kc, 1):
            out[ne.k[f"{nome}_{i}"]] = float(kv)
            out[ne.c[f"{nome}_{i}"]] = float(beta * kv)
    pl = s.payload
    out.update({cp.J_ax: pl.J_ax, cp.I_xx: pl.I_xx, cp.I_yy: pl.I_yy, cp.I_zz: pl.I_zz, cp.I_xy: pl.I_xy,
                cp.I_xz: pl.I_xz, cp.I_yz: pl.I_yz, cp.m_p: pl.m, cp.g: pl.g, cp.d_x: pl.d_x,
                cp.d_y: pl.d_y, cp.d_z: pl.d_z})
    return out
