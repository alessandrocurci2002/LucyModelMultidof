"""Modello FEM torsionale di ciascun asse, isolato (sezione 7).

Elemento a due nodi con rigidezza :math:`\\mathbf K_e = k_e[1, -1; -1, 1]` e due
modelli di massa:

* **concentrata** (7.1): :math:`\\mathbf M_e^\\ell = \\frac{J_e}{2}\\mathbb 1`;
* **consistente** (7.2): :math:`(\\mathbf M_e^c)_{ab} = \\sum_s \\bar J_s
  \\int_{x_s}^{x_s+L_s} N_a N_b\\,dx`, con inerzia per unita' di lunghezza
  costante a tratti.

Poiche' su ogni segmento le funzioni di forma sono lineari in :math:`x` (sia con
``consistent_shape: linear`` sia con ``compliance``, dove sono lineari nella
cedevolezza accumulata, che a sua volta e' lineare a tratti), l'integrale e'
calcolato in forma chiusa con

.. math:: \\int_{x_s}^{x_s+L_s} f\\,g\\,dx = \\frac{L_s}{6}
          \\left(2 f_0 g_0 + f_0 g_1 + f_1 g_0 + 2 f_1 g_1\\right),

esatta per il prodotto di due funzioni lineari (:math:`f_0, f_1` sono i valori
agli estremi del segmento).

Il nodo ausiliario P (sezione 7.5) divide l'elemento del telescopio in 3a e 3b.
Internamente l'ordine dei nodi polari e' 1, 2, 3, P, 4, ..., 11.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Final, Literal, Sequence

import numpy as np
import scipy.linalg as sla

from . import geometry as geo
from .geometry import Pezzo, ProprietaElemento
from .params import Avviso, Elemento, NomeAsse, ParametriAsse, in_inglese

ModelloMassa = Literal["lumped", "consistent"]
MODELLI_MASSA: Final[tuple[ModelloMassa, ...]] = ("lumped", "consistent")


# ---------------------------------------------------------------------------
# Matrici di elemento
# ---------------------------------------------------------------------------


def matrice_rigidezza_elemento(k: float) -> np.ndarray:
    """:math:`\\mathbf K_e = k_e\\begin{bmatrix}1&-1\\\\-1&1\\end{bmatrix}`."""
    return k * np.array([[1.0, -1.0], [-1.0, 1.0]])


def matrice_massa_concentrata(p: ProprietaElemento) -> np.ndarray:
    """Massa concentrata: meta' dell'inerzia (distribuita + extra) su ciascun nodo."""
    return 0.5 * p.J * np.eye(2)


def _valori_forma(p: ProprietaElemento, forma: str) -> list[tuple[float, float]]:
    """Valori di :math:`N_2` agli estremi di ciascun segmento.

    ``linear``: :math:`N_2 = x/L`. ``compliance``: :math:`N_2 = C(x)/C_{tot}`,
    con :math:`C(x) = \\int_0^x \\bar c\\,d\\xi` la cedevolezza accumulata dei
    soli tratti d'albero (la cedevolezza di giunzione, di posizione ignota,
    non entra nella forma). Se la cedevolezza dei tratti e' nulla si ricade
    nella forma lineare.
    """
    L = p.lunghezza
    if forma == "compliance":
        C_tot = sum(s.c_bar * s.lunghezza for s in p.segmenti)
        if C_tot > 0.0:
            valori: list[tuple[float, float]] = []
            C = 0.0
            for s in p.segmenti:
                C1 = C + s.c_bar * s.lunghezza
                valori.append((C / C_tot, C1 / C_tot))
                C = C1
            return valori
    return [(s.x_inizio / L, (s.x_inizio + s.lunghezza) / L) for s in p.segmenti]


def matrice_massa_consistente(p: ProprietaElemento, forma: str = "linear") -> np.ndarray:
    """Massa consistente con inerzia costante a tratti, integrali in forma chiusa.

    Le inerzie extra non distribuite si dividono a meta' sui due nodi.
    """
    M = np.zeros((2, 2))
    for s, (n2_0, n2_1) in zip(p.segmenti, _valori_forma(p, forma)):
        N0 = (1.0 - n2_0, n2_0)
        N1 = (1.0 - n2_1, n2_1)
        for a in range(2):
            for b in range(2):
                M[a, b] += s.J_bar * s.lunghezza / 6.0 * (
                    2.0 * N0[a] * N0[b] + N0[a] * N1[b] + N1[a] * N0[b] + 2.0 * N1[a] * N1[b]
                )
    return M + 0.5 * p.J_extra * np.eye(2)


def matrice_massa(p: ProprietaElemento, modello: ModelloMassa, forma: str = "linear") -> np.ndarray:
    if modello == "lumped":
        return matrice_massa_concentrata(p)
    if modello == "consistent":
        return matrice_massa_consistente(p, forma)
    raise ValueError(f"modello di massa '{modello}' non riconosciuto")


# ---------------------------------------------------------------------------
# Nodo P (sezione 7.5)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DivisioneP:
    """Esito della divisione dell'elemento del telescopio nel nodo P."""

    diviso: bool
    elemento_a: Elemento | None
    elemento_b: Elemento | None
    elemento_intero: Elemento
    """Elemento usato quando non si divide (xi = 0 o 1), con l'eventuale tratto rigido."""
    alias: int | None
    """Nodo utente (base 1) con cui P coincide quando non si divide."""
    x_P: float
    avvisi: tuple[Avviso, ...]


def dividi_elemento_payload(el: Elemento, xi: float, tratto_rigido: float) -> DivisioneP:
    """Divide l'elemento del telescopio nel punto :math:`x_P = \\xi_p L`.

    * ``geometry``: i segmenti vengono tagliati in :math:`x_P` (e ai bordi del
      tratto rigido, centrato in P, se ``tratto_rigido > 0``); inerzie e
      rigidezze di 3a e 3b si ricalcolano con la sezione 6.
    * ``manual``: elemento uniforme, quindi :math:`k_{3a} = k_3 L/(x_P - c/2)`,
      :math:`k_{3b} = k_3 L/(L - x_P - c/2)`, :math:`J_{3a} = \\xi_p J_3`,
      :math:`J_{3b} = (1-\\xi_p) J_3`, che per :math:`c = 0` sono le formule
      della specifica.

    Una cedevolezza di giunzione dell'elemento si ripartisce in proporzione
    alla lunghezza dei due sotto-elementi (ipotesi). Le inerzie extra
    dell'elemento restano a meta' sui nodi 3 e 4 (non su P), come nel modello
    senza divisione: le aggiunge :func:`costruisci_asse`.
    """
    L = el.lunghezza
    x_P = xi * L
    c = tratto_rigido
    a_r, b_r = max(0.0, x_P - c / 2.0), min(L, x_P + c / 2.0)
    avvisi: list[Avviso] = []
    diviso = 0.0 < xi < 1.0

    if diviso and c / 2.0 >= min(x_P, L - x_P):
        raise ValueError(
            f"il tratto rigido del morsetto ({c * 1e3:g} mm) copre per intero uno dei due "
            "sotto-elementi 3a/3b: la rigidezza sarebbe infinita"
        )
    if not diviso and c >= L:
        raise ValueError("il tratto rigido del morsetto copre tutto l'elemento del telescopio")

    def _joint(frazione: float) -> float | None:
        return None if el.k_joint is None else el.k_joint / frazione

    if el.modo == "geometry":
        tagli = [x_P] + ([a_r, b_r] if c > 0.0 else [])
        segs = geo.taglia_segmenti(el.segmenti, tagli)
        if c > 0.0:
            segs = geo.marca_rigidi(segs, a_r, b_r)
        if not diviso:
            intero = replace(el, segmenti=tuple(segs), extra=())
            return DivisioneP(False, None, None, intero, el.id if xi == 0.0 else el.id + 1,
                              x_P, tuple(avvisi))
        segs_a, segs_b, x = [], [], 0.0
        for s in segs:
            (segs_a if x + s.lunghezza / 2.0 < x_P else segs_b).append(s)
            x += s.lunghezza
        el_a = replace(el, lunghezza=x_P, segmenti=tuple(segs_a), extra=(),
                       k_joint=_joint(xi))
        el_b = replace(el, lunghezza=L - x_P, segmenti=tuple(segs_b), extra=(),
                       k_joint=_joint(1.0 - xi))
        return DivisioneP(True, el_a, el_b, el, None, x_P, tuple(avvisi))

    # modalita' manual: elemento uniforme
    assert el.J_manuale is not None and el.k_manuale is not None
    avvisi.append(Avviso(
        "P_MANUALE_UNIFORME", "informazione",
        f"Polar element {el.id} in manual mode: the split at node P assumes "
        "a uniform element (k and J shared in proportion to the length)."))
    cedevolezza_unitaria = 1.0 / (el.k_manuale * L)
    if not diviso:
        k_int = 1.0 / (cedevolezza_unitaria * (L - (b_r - a_r)))
        intero = replace(el, k_manuale=k_int, extra=())
        return DivisioneP(False, None, None, intero, el.id if xi == 0.0 else el.id + 1,
                          x_P, tuple(avvisi))
    k_a = 1.0 / (cedevolezza_unitaria * (x_P - c / 2.0))
    k_b = 1.0 / (cedevolezza_unitaria * (L - x_P - c / 2.0))
    el_a = replace(el, lunghezza=x_P, J_manuale=xi * el.J_manuale, k_manuale=k_a, extra=(),
                   k_joint=_joint(xi), m_manuale=None, I_t_manuale=None, y_c_manuale=None)
    el_b = replace(el, lunghezza=L - x_P, J_manuale=(1.0 - xi) * el.J_manuale, k_manuale=k_b,
                   extra=(), k_joint=_joint(1.0 - xi), m_manuale=None, I_t_manuale=None,
                   y_c_manuale=None)
    return DivisioneP(True, el_a, el_b, el, None, x_P, tuple(avvisi))


# ---------------------------------------------------------------------------
# Modello di un asse
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PezzoPosizionato:
    """Pezzo rigido dell'asse con la posizione assiale del baricentro dal nodo 1."""

    pezzo: Pezzo
    x: float


@dataclass(frozen=True)
class ModelloAsse:
    """Modello FEM completo di un asse, senza telescopio ne' termini del load.

    Attributes
    ----------
    etichette
        Etichette dei nodi nell'ordine interno (``"1"``, ..., ``"P"``, ...).
    x
        Posizioni dei nodi dal nodo 1, m.
    elementi, collegamenti
        Elementi FEM effettivi (con 3a e 3b) e coppie di indici interni.
    elementi_originali
        Elementi dell'utente, non divisi: servono alle tabelle e a
        :math:`J^{ax}`.
    K, M
        Rigidezza degli alberi e masse (alberi + inerzie nodali) per modello.
    C_b
        Attrito viscoso dei cuscinetti (sezione 7.3).
    masters
        Indici interni dei nodi master, ordinati lungo la catena.
    pezzi
        Pezzi rigidi dell'asse posizionati, per :math:`J^{ax}` e sottrazioni.
    """

    nome: NomeAsse
    etichette: tuple[str, ...]
    x: np.ndarray
    elementi: tuple[ProprietaElemento, ...]
    collegamenti: tuple[tuple[int, int], ...]
    elementi_originali: tuple[ProprietaElemento, ...]
    K: np.ndarray
    M: dict[str, np.ndarray]
    C_b: np.ndarray
    J_nodali: np.ndarray
    i_encoder: int
    i_motore: int
    i_estremo: int
    i_P: int | None
    masters: tuple[int, ...]
    pezzi: tuple[PezzoPosizionato, ...]
    forma_consistente: str
    divisione: DivisioneP | None
    segnaposto: bool
    avvisi: tuple[Avviso, ...] = field(default_factory=tuple)
    alias: dict[str, int] = field(default_factory=dict)
    """Nodi utente che coincidono con P quando il telescopio non ha albero."""
    i_H: int | None = None
    """Centro del mozzo (modello ``hub``): nodo della catena a cui e' appeso il telescopio P."""

    @property
    def ramo(self) -> bool:
        """``True`` se il telescopio P e' un nodo fuori dalla catena, appeso a H tramite k_b."""
        return self.i_H is not None and self.i_P is not None and self.i_H != self.i_P

    @property
    def n(self) -> int:
        return len(self.etichette)

    @property
    def J_totale(self) -> float:
        """Inerzia assiale totale: elementi effettivi piu' inerzie nodali.

        Le inerzie extra dell'elemento del telescopio sono gia' in
        ``J_nodali`` (meta' sui nodi 3 e 4), non negli elementi 3a e 3b.
        """
        return float(sum(p.J for p in self.elementi) + self.J_nodali.sum())

    def indice(self, nodo: int | str) -> int:
        """Indice interno del nodo utente (base 1, oppure ``"P"``)."""
        if str(nodo) == "P":
            if self.i_P is None:
                raise KeyError(f"il nodo P non esiste sull'asse {self.nome}")
            return self.i_P
        if str(nodo) in self.alias:
            return self.alias[str(nodo)]
        return self.etichette.index(str(nodo))

    @property
    def telescopio_rigido(self) -> bool:
        """``True`` se l'elemento del telescopio non ha albero (nodi 3, 4 e P coincidenti)."""
        return bool(self.alias)

    def etichetta(self, i: int) -> str:
        """Etichetta del nodo di indice interno ``i`` (``"P"`` anche se coincide con 3 o 4)."""
        return "P" if i == self.i_P else self.etichette[i]

    @property
    def etichette_master(self) -> tuple[str, ...]:
        return tuple(self.etichetta(i) for i in self.masters)

    def vettore_forze(self) -> np.ndarray:
        """Vettore :math:`\\mathbf e_p` del nodo motore: :math:`\\mathbf Q = \\tau\\mathbf e_p`."""
        e = np.zeros(self.n)
        e[self.i_motore] = 1.0
        return e


def _proprieta(el: Elemento, etichetta: str | None = None) -> ProprietaElemento:
    return geo.proprieta_elemento(el, etichetta)


def costruisci_asse(asse: ParametriAsse, forma_consistente: str = "linear", *,
                   con_telescopio: bool = True) -> ModelloAsse:
    """Assembla il modello FEM completo di un asse (sezioni 7.1-7.3 e 7.5).

    Modello ``hub`` (``asse.mozzo``): le due meta' del mozzo sostituiscono gli elementi
    ``elementi_lato`` (molle senza massa), il centro H prende i nodi dell'elemento del telescopio
    con l'inerzia del mozzo, e il telescopio e' il nodo P fuori dalla catena, collegato a H dalla
    molla k_b. Con ``con_telescopio = False`` (prova senza telescopio) il ramo non c'e' e H fa da P.
    """
    from . import mozzo as mz

    avvisi: list[Avviso] = []
    pm = asse.mozzo
    originali = tuple(
        mz.proprieta_lato(pm, str(e.id), e.segnaposto or pm.geometria_segnaposto)
        if pm is not None and e.id in pm.elementi_lato else _proprieta(e) for e in asse.elementi)
    ramo = pm is not None and con_telescopio

    # --- elementi effettivi, con l'eventuale divisione nel nodo P -----------
    elementi: list[ProprietaElemento] = []
    etichette: list[str] = ["1"]
    x_nodi: list[float] = [0.0]
    collegamenti: list[tuple[int, int]] = []
    extra_nodali: dict[str, float] = {}
    divisione: DivisioneP | None = None
    i_P: int | None = None
    alias: dict[str, int] = {}
    xu = 0.0  # posizione reale dell'ultimo nodo utente

    for el, p_orig in zip(asse.elementi, originali):
        i = len(etichette) - 1
        if asse.nome == "polare" and el.id == asse.elemento_payload and el.modo == "telescope":
            # Nessun albero: il telescopio rigido rende i nodi i e i+1 un unico nodo, P.
            # Lo si disegna nella posizione del centro del telescopio.
            assert asse.xi_payload is not None
            i_P = i
            alias[str(el.id)] = alias[str(el.id + 1)] = i
            etichette[i] = "H" if ramo else "P"
            x_nodi[i] = xu + asse.xi_payload * el.lunghezza
            xu += el.lunghezza
            continue
        if asse.nome == "polare" and el.id == asse.elemento_payload:
            assert asse.xi_payload is not None
            divisione = dividi_elemento_payload(el, asse.xi_payload, asse.tratto_rigido)
            avvisi.extend(divisione.avvisi)
            # le inerzie extra restano a meta' sui nodi utente i e i+1
            J_ex = p_orig.J_extra
            if J_ex > 0.0:
                extra_nodali[str(el.id)] = extra_nodali.get(str(el.id), 0.0) + J_ex / 2
                extra_nodali[str(el.id + 1)] = extra_nodali.get(str(el.id + 1), 0.0) + J_ex / 2
            if divisione.diviso:
                assert divisione.elemento_a is not None and divisione.elemento_b is not None
                p_a = _proprieta(divisione.elemento_a, f"{el.id}a")
                p_b = _proprieta(divisione.elemento_b, f"{el.id}b")
                etichette.append("P")
                x_nodi.append(xu + p_a.lunghezza)
                i_P = len(etichette) - 1
                collegamenti.append((i, i_P))
                elementi.append(p_a)
                etichette.append(str(el.id + 1))
                xu += p_a.lunghezza + p_b.lunghezza
                x_nodi.append(xu)
                collegamenti.append((i_P, i_P + 1))
                elementi.append(p_b)
                continue
            p_int = _proprieta(divisione.elemento_intero, str(el.id))
            etichette.append(str(el.id + 1))
            xu += p_int.lunghezza
            x_nodi.append(xu)
            collegamenti.append((i, i + 1))
            elementi.append(p_int)
            i_P = etichette.index(str(divisione.alias))
            continue
        etichette.append(str(el.id + 1))
        xu += p_orig.lunghezza
        x_nodi.append(xu)
        collegamenti.append((i, i + 1))
        elementi.append(p_orig)

    # --- modello hub: telescopio P appeso al centro del mozzo H tramite k_b ----------
    i_H: int | None = None
    if pm is not None:
        i_H = i_P
        assert i_H is not None
        J_H = pm.J_polare
        if not con_telescopio and pm.tratto_telescopio is not None:
            # prova senza telescopio: il tratto di tubo stretto nell'anello se ne va con il telescopio
            J_H -= mz.inerzia_tratto(pm)[2]  # type: ignore[index]
        extra_nodali[str(pm.elemento_centro)] = extra_nodali.get(str(pm.elemento_centro), 0.0) + J_H
        if ramo:
            etichette.append("P")
            x_nodi.append(x_nodi[i_H])
            i_P = len(etichette) - 1
            collegamenti.append((i_H, i_P))
            elementi.append(mz.proprieta_attacco(pm))

    # --- azimut: telescopio P appeso al load L (forcella + mozzo + anello) tramite k ---------
    if asse.attacco_telescopio is not None and con_telescopio:
        ramo = True
        i_H = etichette.index(str(asse.nodo_estremo))
        etichette.append("P")
        x_nodi.append(x_nodi[i_H])
        i_P = len(etichette) - 1
        collegamenti.append((i_H, i_P))
        elementi.append(mz.proprieta_attacco_azimut(asse.attacco_telescopio))

    n = len(etichette)
    idx = {e: k for k, e in enumerate(etichette)}
    idx.update(alias)

    # --- assemblaggio --------------------------------------------------------
    K = np.zeros((n, n))
    M = {m: np.zeros((n, n)) for m in MODELLI_MASSA}
    for p, (a, b) in zip(elementi, collegamenti):
        ids = np.ix_([a, b], [a, b])
        K[ids] += matrice_rigidezza_elemento(p.k)
        for m in MODELLI_MASSA:
            M[m][ids] += matrice_massa(p, m, forma_consistente)

    J_nodali = np.zeros(n)
    for nodale in asse.inerzie_nodali:
        J_nodali[idx[str(nodale.nodo)]] += nodale.J
    for nodo, J in extra_nodali.items():
        J_nodali[idx[nodo]] += J
    for m in MODELLI_MASSA:
        M[m] += np.diag(J_nodali)

    # --- attrito dei cuscinetti ------------------------------------------------
    C_b = np.zeros((n, n))
    for el_id, b in asse.b_cuscinetti.items():
        a_i, b_i = idx[str(el_id)], idx[str(el_id + 1)]
        C_b[a_i, a_i] += b / 2.0
        C_b[b_i, b_i] += b / 2.0

    # --- master ------------------------------------------------------------
    # con il ramo H-P (L-P sull'azimut) sono master sia il nodo d'aggancio sia il telescopio
    masters = sorted({(i_P if str(m) == "P" else idx[str(m)]) for m in asse.masters}  # type: ignore[misc]
                     | ({i_H, i_P} if ramo and i_H is not None and i_P is not None else set()))
    assert i_P is None or i_P in masters

    # --- pezzi posizionati (per J^ax e sottrazioni) --------------------------
    x_utente = {str(k + 1): 0.0 for k in range(asse.n_nodi)}
    xc = 0.0
    for k, p in enumerate(originali):
        x_utente[str(k + 1)] = xc
        xc += p.lunghezza
    x_utente[str(asse.n_nodi)] = xc
    pezzi: list[PezzoPosizionato] = []
    for k, p in enumerate(originali):
        for pz in p.pezzi:
            pezzi.append(PezzoPosizionato(pezzo=pz, x=x_utente[str(k + 1)] + pz.x_c))
    for nodale in asse.inerzie_nodali:
        pezzi.append(PezzoPosizionato(
            pezzo=Pezzo(nome=f"nodo {nodale.nodo} '{nodale.nome}'", m=nodale.m or 0.0,
                        J=nodale.J, I_t=nodale.I_t or 0.0, x_c=0.0),
            x=x_utente[str(nodale.nodo)]))

    modello = ModelloAsse(
        nome=asse.nome,
        etichette=tuple(etichette),
        x=np.array(x_nodi),
        elementi=tuple(elementi),
        collegamenti=tuple(collegamenti),
        elementi_originali=originali,
        K=K,
        M=M,
        C_b=C_b,
        J_nodali=J_nodali,
        i_encoder=idx[str(asse.nodo_encoder)],
        i_motore=idx[str(asse.nodo_coppia)],
        i_estremo=idx[str(asse.nodo_estremo)],
        i_P=i_P,
        masters=tuple(masters),
        pezzi=tuple(pezzi),
        forma_consistente=forma_consistente,
        divisione=divisione,
        segnaposto=asse.segnaposto,
        avvisi=tuple(avvisi),
        alias=alias,
        i_H=i_H if ramo else None,
    )
    for m in MODELLI_MASSA:
        verifica_matrici(modello.M[m], modello.K, modello.J_totale, etichetta=f"{asse.nome}/{m}",
                         senza_massa=(i_P,) if ramo and i_P is not None else ())

    # M^l - M^c e' semidefinita solo se ogni elemento ha A_e = C_e (sezione 7.4):
    # con inerzia non simmetrica rispetto alla mezzeria l'ordinamento
    # omega^l <= omega^c non e' garantito, e va detto invece di supporlo.
    asimmetrici = [p.etichetta for p in elementi if abs(asimmetria_elemento(p)) > 1e-9 * p.J]
    if forma_consistente == "linear" and asimmetrici:
        lam = float(autovalori_differenza_masse(modello).min())
        avvisi.append(Avviso(
            f"MASSA_NON_ORDINATA_{asse.nome.upper()}", "informazione",
            f"{in_inglese(asse.nome).capitalize()} axis: elements {asimmetrici} have an inertia that is not "
            "symmetric about the midpoint, so M^l - M^c is not semidefinite (smallest eigenvalue "
            f"{lam:.3e} kg m^2) and omega^l <= omega^c is not guaranteed for every mode "
            "(section 7.4: the property holds for uniform elements)."))
        modello = replace(modello, avvisi=tuple(avvisi))
    return modello


def asimmetria_elemento(p: ProprietaElemento) -> float:
    """:math:`\\delta_e = (C_e - A_e)/2`, con :math:`A_e = \\int N_1^2\\bar J`,
    :math:`C_e = \\int N_2^2\\bar J` (funzioni di forma lineari).

    Per l'elemento vale :math:`\\det(\\mathbf M_e^\\ell - \\mathbf M_e^c) = -\\delta_e^2`:
    la differenza e' semidefinita se e solo se :math:`\\delta_e = 0`.
    """
    Mc = matrice_massa_consistente(p, "linear")
    return 0.5 * float(Mc[1, 1] - Mc[0, 0])


# ---------------------------------------------------------------------------
# Matrici del modello isolato
# ---------------------------------------------------------------------------


def matrici_isolate(
    mod: ModelloAsse,
    modello_massa: ModelloMassa,
    *,
    J_aggiunta: float = 0.0,
    k_aggiunta: float = 0.0,
    J_telescopio: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Matrici dell'asse isolato con le inerzie dominanti (sezione 7.3).

    * polare: ``J_aggiunta`` = :math:`I_{yy,p}` e ``k_aggiunta`` =
      :math:`k_g(\\theta_0)` al nodo P;
    * azimut: ``J_aggiunta`` = :math:`J^{ax}` e ``J_telescopio`` = :math:`J_{\\psi,p}(\\theta_0)` al
      nodo 9; con il telescopio appeso al load (``azimuth.telescope_attachment``) ``J_telescopio`` va
      al nodo P.

    Nel modello accoppiato queste grandezze entrano solo attraverso :math:`T` e
    :math:`V` e **non** vanno sommate qui (flag ``isolated`` della specifica:
    si chiama questa funzione solo per l'asse isolato).
    """
    M = mod.M[modello_massa].copy()
    K = mod.K.copy()
    i = mod.i_P if mod.nome == "polare" else mod.i_estremo
    i_tel = mod.i_P if mod.nome == "azimut" and mod.ramo else i
    assert i is not None and i_tel is not None
    M[i, i] += J_aggiunta
    M[i_tel, i_tel] += J_telescopio
    K[i, i] += k_aggiunta
    return M, K


def beta_smorzamento(M: np.ndarray, K_alberi: np.ndarray, zeta: float) -> float:
    """:math:`\\beta = 2\\zeta/\\omega_1`, con :math:`\\omega_1` il primo modo flessibile.

    ``K_alberi`` e' la sola rigidezza degli alberi, che ha esattamente un modo
    rigido: la rigidezza gravitazionale non e' strutturale e non entra ne' in
    :math:`\\beta` ne' in :math:`\\mathbf C = \\beta\\mathbf K`.
    """
    w = frequenze_proprie(M, K_alberi) * 2.0 * np.pi
    flessibili = w[w > 1e-6 * max(w.max(), 1.0)]
    if flessibili.size == 0:
        return 0.0
    return 2.0 * zeta / float(flessibili[0])


def matrice_smorzamento(K_alberi: np.ndarray, C_b: np.ndarray, beta: float) -> np.ndarray:
    """:math:`\\mathbf C = \\beta\\mathbf K + \\mathbf C_b` (sezione 7.3)."""
    return beta * K_alberi + C_b


# ---------------------------------------------------------------------------
# Analisi modale
# ---------------------------------------------------------------------------


def modi(M: np.ndarray, K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Autovalori e forme modali normalizzate in massa di :math:`(\\mathbf K, \\mathbf M)`.

    Returns
    -------
    (omega2, Phi)
        :math:`\\omega^2` crescenti (negativi se la rigidezza gravitazionale
        rende instabile il moto rigido) e forme modali come colonne.
    """
    omega2, Phi = sla.eigh(K, M)
    return omega2, Phi


def frequenze_proprie(M: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Frequenze proprie in Hz, con segno: negative per autovalori negativi.

    Un autovalore negativo (rigidezza gravitazionale negativa) e' riportato
    come :math:`-\\sqrt{|\\omega^2|}/2\\pi`, per non confonderlo con un modo.
    """
    omega2, _ = modi(M, K)
    scala = max(float(np.max(np.abs(omega2))), 1.0)
    omega2 = np.where(np.abs(omega2) < 1e-12 * scala, 0.0, omega2)
    return np.sign(omega2) * np.sqrt(np.abs(omega2)) / (2.0 * np.pi)


# ---------------------------------------------------------------------------
# Controlli (sezione 7.4)
# ---------------------------------------------------------------------------


def verifica_matrici(M: np.ndarray, K: np.ndarray, J_totale: float, *, etichetta: str = "",
                     senza_massa: tuple[int, ...] = ()) -> None:
    """Simmetria, definita positivita' di M, :math:`K\\mathbf 1 = 0`, inerzia totale.

    ``senza_massa``: nodi ammessi senza inerzia (il telescopio P del modello ``hub``, la cui
    inerzia entra dal payload); M deve essere definita positiva sugli altri nodi.

    Raises
    ------
    AssertionError
        Se uno dei controlli fallisce.
    """
    scala_K = max(float(np.max(np.abs(K))), 1e-30)
    assert np.allclose(M, M.T, atol=0.0, rtol=1e-14), f"{etichetta}: M non simmetrica"
    assert np.allclose(K, K.T, atol=1e-14 * scala_K), f"{etichetta}: K non simmetrica"
    altri = [i for i in range(M.shape[0]) if i not in senza_massa]
    assert np.all(np.linalg.eigvalsh(M[np.ix_(altri, altri)]) > 0.0), f"{etichetta}: M non definita positiva"
    uno = np.ones(K.shape[0])
    assert np.max(np.abs(K @ uno)) <= 1e-10 * scala_K, f"{etichetta}: K 1 != 0"
    somma = float(uno @ M @ uno)
    assert abs(somma - J_totale) <= 1e-12 * max(J_totale, 1e-30), (
        f"{etichetta}: 1^T M 1 = {somma:.12e} diverso dall'inerzia totale {J_totale:.12e}"
    )


def autovalori_differenza_masse(mod: ModelloAsse) -> np.ndarray:
    """Autovalori di :math:`\\mathbf M^\\ell - \\mathbf M^c`: tutti :math:`\\ge 0`
    con funzioni di forma lineari (sezione 7.4)."""
    return np.linalg.eigvalsh(mod.M["lumped"] - mod.M["consistent"])


def numero_modi_rigidi(K: np.ndarray, M: np.ndarray) -> int:
    """Numero di autovalori nulli del problema :math:`(\\mathbf K, \\mathbf M)`."""
    omega2, _ = modi(M, K)
    return int(np.sum(np.abs(omega2) < 1e-9 * max(float(np.max(np.abs(omega2))), 1.0)))


def catena_uniforme(k: Sequence[float], J: Sequence[float], modello: ModelloMassa) -> tuple[np.ndarray, np.ndarray]:
    """Matrici di una catena di elementi uniformi (usata nei test di convergenza)."""
    n = len(k) + 1
    K = np.zeros((n, n))
    M = np.zeros((n, n))
    for i, (ki, Ji) in enumerate(zip(k, J)):
        ids = np.ix_([i, i + 1], [i, i + 1])
        K[ids] += matrice_rigidezza_elemento(ki)
        if modello == "lumped":
            M[ids] += 0.5 * Ji * np.eye(2)
        else:
            M[ids] += Ji / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
    return M, K
