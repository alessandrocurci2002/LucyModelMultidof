"""Proprieta' degli elementi d'albero per composizione di cilindri (sezione 6).

Per uno strato anulare di lunghezza :math:`L`, diametri :math:`D_o`, :math:`D_i`
(:math:`D_i = 0` per il pieno), densita' :math:`\\rho` e modulo :math:`G`:

.. math::
    m = \\rho\\frac{\\pi}{4}(D_o^2 - D_i^2)L, \\qquad
    J = \\frac{\\pi\\rho L}{32}(D_o^4 - D_i^4) = \\frac{m}{8}(D_o^2 + D_i^2),

.. math::
    I_t = m\\left[\\frac{D_o^2 + D_i^2}{16} + \\frac{L^2}{12}\\right], \\qquad
    J_p = \\frac{\\pi}{32}(D_o^4 - D_i^4).

Composizione: gli strati coassiali di un segmento sommano inerzie e masse, e le
rigidezze dei soli strati portanti in parallelo,
:math:`k_s = \\frac{1}{L_s}\\sum_\\ell G_\\ell J_{p,\\ell}`; i segmenti di un
elemento sommano inerzie e masse, e le cedevolezze in serie,
:math:`1/k_e = \\sum_s 1/k_s + 1/k_{joint}`.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .params import MM, Elemento, Segmento, Strato

# ---------------------------------------------------------------------------
# Strato
# ---------------------------------------------------------------------------


def massa_strato(rho: float, d_out: float, d_in: float, L: float) -> float:
    """Massa di un cilindro (cavo), kg."""
    return rho * np.pi / 4.0 * (d_out**2 - d_in**2) * L


def momento_polare_area(d_out: float, d_in: float) -> float:
    """Momento polare d'area :math:`J_p = \\pi(D_o^4 - D_i^4)/32`, m^4."""
    return np.pi / 32.0 * (d_out**4 - d_in**4)


def inerzia_assiale(rho: float, d_out: float, d_in: float, L: float) -> float:
    """Inerzia assiale :math:`J = \\pi\\rho L (D_o^4 - D_i^4)/32`, kg m^2."""
    return rho * L * momento_polare_area(d_out, d_in)


def inerzia_trasversale(m: float, d_out: float, d_in: float, L: float) -> float:
    """Inerzia trasversale nel baricentro :math:`m[(D_o^2+D_i^2)/16 + L^2/12]`."""
    return m * ((d_out**2 + d_in**2) / 16.0 + L**2 / 12.0)


@dataclass(frozen=True)
class ProprietaStrato:
    """Contributo di uno strato, per la tabella di scomposizione."""

    elemento: str
    segmento: int
    strato: int
    forma: str
    d_out: float
    d_in: float
    materiale: str
    portante: bool
    lunghezza: float
    x_c: float
    """Baricentro dello strato dal nodo ``i`` dell'elemento, m."""
    m: float
    J: float
    I_t: float
    J_p: float
    GJ_p: float
    """Rigidezza torsionale della sezione :math:`G J_p` (0 se non portante)."""
    rigido: bool


def proprieta_strato(
    strato: Strato, L: float, *, elemento: str = "", segmento: int = 0, indice: int = 0,
    x_inizio: float = 0.0, rigido: bool = False,
) -> ProprietaStrato:
    """Massa, inerzie e rigidezza di sezione di uno strato di lunghezza ``L``."""
    rho, G = strato.materiale.rho, strato.materiale.G
    m = massa_strato(rho, strato.d_out, strato.d_in, L)
    J_p = momento_polare_area(strato.d_out, strato.d_in)
    return ProprietaStrato(
        elemento=elemento, segmento=segmento, strato=indice, forma=strato.forma,
        d_out=strato.d_out, d_in=strato.d_in, materiale=strato.materiale.nome,
        portante=strato.portante, lunghezza=L, x_c=x_inizio + L / 2.0, m=m,
        J=rho * L * J_p, I_t=inerzia_trasversale(m, strato.d_out, strato.d_in, L),
        J_p=J_p, GJ_p=G * J_p if strato.portante else 0.0, rigido=rigido,
    )


# ---------------------------------------------------------------------------
# Elemento
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Pezzo:
    """Pezzo assialsimmetrico rigido, per :math:`J^{ax}` e per le sottrazioni.

    ``x_c`` e' la posizione del baricentro dal nodo ``i`` dell'elemento.
    """

    nome: str
    m: float
    J: float
    I_t: float
    x_c: float


@dataclass(frozen=True)
class SegmentoFEM:
    """Tratto a sezione costante, in coordinate locali dell'elemento.

    Serve alla massa consistente (sezione 7.2): ``J_bar`` e' l'inerzia per
    unita' di lunghezza, ``c_bar`` la cedevolezza per unita' di lunghezza
    :math:`1/\\sum G J_p` (0 se il tratto e' rigido).
    """

    x_inizio: float
    lunghezza: float
    J_bar: float
    c_bar: float


@dataclass(frozen=True)
class ProprietaElemento:
    """Proprieta' di un elemento, in SI.

    Attributes
    ----------
    J_distribuita
        Inerzia assiale distribuita lungo l'elemento (strati, o J manuale).
    J_extra
        Inerzie concentrate dell'elemento: in entrambi i modelli di massa si
        dividono a meta' sui due nodi.
    k
        Rigidezza torsionale equivalente, N m/rad.
    m, I_t, x_G
        Massa totale, inerzia trasversale nel baricentro dell'elemento e
        posizione del baricentro dal nodo ``i``.
    """

    etichetta: str
    modo: str
    lunghezza: float
    J_distribuita: float
    J_extra: float
    k: float
    k_segmenti: tuple[float, ...]
    k_joint: float | None
    m: float
    I_t: float
    x_G: float
    segmenti: tuple[SegmentoFEM, ...]
    strati: tuple[ProprietaStrato, ...]
    pezzi: tuple[Pezzo, ...]
    segnaposto: bool

    @property
    def J(self) -> float:
        """Inerzia assiale totale dell'elemento."""
        return self.J_distribuita + self.J_extra


def _rigidezza_serie(cedevolezze: Iterable[float]) -> float:
    c = float(sum(cedevolezze))
    if c <= 0.0:
        raise ValueError("cedevolezza totale nulla: l'elemento sarebbe infinitamente rigido")
    return 1.0 / c


def proprieta_elemento(el: Elemento, etichetta: str | None = None) -> ProprietaElemento:
    """Calcola inerzia, rigidezza, massa e scomposizione di un elemento.

    Modalita' ``geometry``: composizione degli strati e dei segmenti come in
    sezione 6. Se ``k_manuale`` e' dato, sostituisce la rigidezza dei segmenti
    (l'inerzia resta quella geometrica). Modalita' ``manual``: elemento
    uniforme con ``J_manuale`` e ``k_manuale``; il contributo trasversale viene
    da ``m_manuale``, ``I_t_manuale`` e ``y_c_manuale`` se dati, altrimenti e'
    nullo (con avviso in :mod:`params`).

    In entrambe le modalita' una ``k_joint`` si somma in serie. Le inerzie
    extra non distribuite sono collocate a meta' elemento per il calcolo di
    :math:`J^{ax}`.
    """
    nome = etichetta or str(el.id)
    if el.modo == "telescope":
        # nessun albero: inerzia nulla (il telescopio e' nel payload), collegamento rigido
        return ProprietaElemento(
            etichetta=nome, modo="telescope", lunghezza=el.lunghezza, J_distribuita=0.0,
            J_extra=0.0, k=float("inf"), k_segmenti=(), k_joint=None, m=0.0, I_t=0.0,
            x_G=el.lunghezza / 2.0, segmenti=(), strati=(), pezzi=(), segnaposto=el.segnaposto)
    strati: list[ProprietaStrato] = []
    segmenti: list[SegmentoFEM] = []
    pezzi: list[Pezzo] = []
    k_seg: list[float] = []

    if el.modo == "geometry":
        x = 0.0
        for s, seg in enumerate(el.segmenti, start=1):
            props = [
                proprieta_strato(st, seg.lunghezza, elemento=nome, segmento=s, indice=j,
                                 x_inizio=x, rigido=seg.rigido)
                for j, st in enumerate(seg.strati, start=1)
            ]
            strati.extend(props)
            J_s = sum(p.J for p in props)
            GJ = sum(p.GJ_p for p in props)
            if seg.rigido:
                c_bar = 0.0
            else:
                if GJ <= 0.0 and el.k_manuale is None:
                    raise ValueError(f"elemento {nome}, segmento {s}: nessuno strato portante")
                c_bar = 1.0 / GJ if GJ > 0.0 else 0.0
                if GJ > 0.0:
                    k_seg.append(GJ / seg.lunghezza)
            segmenti.append(SegmentoFEM(x_inizio=x, lunghezza=seg.lunghezza,
                                        J_bar=J_s / seg.lunghezza, c_bar=c_bar))
            pezzi.extend(
                Pezzo(nome=f"el.{nome} seg.{s} str.{p.strato}", m=p.m, J=p.J, I_t=p.I_t, x_c=p.x_c)
                for p in props
            )
            x += seg.lunghezza
        J_dist = float(sum(p.J for p in strati))
        if el.k_manuale is not None:
            cedevolezze = [1.0 / el.k_manuale]
        else:
            cedevolezze = [1.0 / k for k in k_seg]
    else:
        assert el.J_manuale is not None and el.k_manuale is not None
        J_dist = float(el.J_manuale)
        L = el.lunghezza
        segmenti.append(SegmentoFEM(x_inizio=0.0, lunghezza=L, J_bar=J_dist / L,
                                    c_bar=1.0 / (el.k_manuale * L)))
        k_seg.append(el.k_manuale)
        cedevolezze = [1.0 / el.k_manuale]
        if None not in (el.m_manuale, el.I_t_manuale, el.y_c_manuale):
            pezzi.append(Pezzo(nome=f"el.{nome} (manuale)", m=float(el.m_manuale),  # type: ignore[arg-type]
                               J=J_dist, I_t=float(el.I_t_manuale),  # type: ignore[arg-type]
                               x_c=float(el.y_c_manuale)))  # type: ignore[arg-type]
        else:
            pezzi.append(Pezzo(nome=f"el.{nome} (manuale, senza m/I_t)", m=0.0, J=J_dist,
                               I_t=0.0, x_c=L / 2.0))

    if el.k_joint is not None:
        cedevolezze.append(1.0 / el.k_joint)
    k = _rigidezza_serie(cedevolezze)

    J_extra = float(sum(e.J for e in el.extra))
    for e in el.extra:
        pezzi.append(Pezzo(nome=f"el.{nome} extra '{e.nome}'", m=e.m or 0.0, J=e.J,
                           I_t=e.I_t or 0.0, x_c=el.lunghezza / 2.0))

    m_tot = float(sum(p.m for p in pezzi))
    x_G = float(sum(p.m * p.x_c for p in pezzi) / m_tot) if m_tot > 0.0 else el.lunghezza / 2.0
    I_t = float(sum(p.I_t + p.m * (p.x_c - x_G) ** 2 for p in pezzi))

    return ProprietaElemento(
        etichetta=nome, modo=el.modo, lunghezza=el.lunghezza, J_distribuita=J_dist,
        J_extra=J_extra, k=k, k_segmenti=tuple(k_seg), k_joint=el.k_joint, m=m_tot,
        I_t=I_t, x_G=x_G, segmenti=tuple(segmenti), strati=tuple(strati),
        pezzi=tuple(pezzi), segnaposto=el.segnaposto,
    )


# ---------------------------------------------------------------------------
# Taglio dei segmenti (serve al nodo P, sezione 7.5)
# ---------------------------------------------------------------------------


def taglia_segmenti(segmenti: Sequence[Segmento], tagli: Iterable[float]) -> list[Segmento]:
    """Taglia i segmenti nelle ascisse ``tagli`` (dal nodo ``i``).

    Un segmento attraversato da un taglio diventa due segmenti con gli stessi
    strati: massa e inerzia si conservano esattamente, perche' entrambe sono
    lineari nella lunghezza.
    """
    punti = sorted({float(t) for t in tagli})
    out: list[Segmento] = []
    x = 0.0
    for seg in segmenti:
        a, b = x, x + seg.lunghezza
        interni = [t for t in punti if a + 1e-12 < t < b - 1e-12]
        bordi = [a, *interni, b]
        for x0, x1 in zip(bordi[:-1], bordi[1:]):
            out.append(Segmento(lunghezza=x1 - x0, strati=seg.strati, rigido=seg.rigido))
        x = b
    return out


def marca_rigidi(segmenti: Sequence[Segmento], a: float, b: float) -> list[Segmento]:
    """Rende rigidi (cedevolezza nulla) i segmenti contenuti in ``[a, b]``.

    I segmenti devono essere gia' tagliati in ``a`` e ``b``.
    """
    out: list[Segmento] = []
    x = 0.0
    for seg in segmenti:
        centro = x + seg.lunghezza / 2.0
        rigido = seg.rigido or (a - 1e-12 <= centro <= b + 1e-12 and b > a)
        out.append(Segmento(lunghezza=seg.lunghezza, strati=seg.strati, rigido=rigido))
        x += seg.lunghezza
    return out


# ---------------------------------------------------------------------------
# Tabella di scomposizione
# ---------------------------------------------------------------------------

INTESTAZIONE_SCOMPOSIZIONE: tuple[str, ...] = (
    "elemento", "segmento", "strato", "forma", "d_out_mm", "d_in_mm", "L_mm", "materiale",
    "portante", "rigido", "x_c_mm", "m_kg", "J_kgm2", "I_t_kgm2", "GJp_Nm2", "k_Nm_rad",
    "segnaposto",
)


def righe_scomposizione(props: Sequence[ProprietaElemento]) -> list[list[object]]:
    """Righe della tabella ``element_breakdown_<asse>.csv``.

    Per ogni elemento: una riga per strato e una riga di totale, cosi' che ogni
    numero del modello sia tracciabile fino al singolo cilindro.
    """
    righe: list[list[object]] = []
    for p in props:
        for s in p.strati:
            righe.append([
                p.etichetta, s.segmento, s.strato, s.forma, s.d_out / MM, s.d_in / MM,
                s.lunghezza / MM, s.materiale, s.portante, s.rigido, s.x_c / MM,
                f"{s.m:.6e}", f"{s.J:.6e}", f"{s.I_t:.6e}", f"{s.GJ_p:.6e}",
                f"{s.GJ_p / s.lunghezza:.6e}" if s.GJ_p > 0 else "", p.segnaposto,
            ])
        for pz in p.pezzi:
            if "extra" in pz.nome or "manuale" in pz.nome:
                righe.append([p.etichetta, "", "", pz.nome, "", "", p.lunghezza / MM, "", "",
                              "", pz.x_c / MM, f"{pz.m:.6e}", f"{pz.J:.6e}", f"{pz.I_t:.6e}",
                              "", "", p.segnaposto])
        righe.append([
            p.etichetta, "TOTALE", "", p.modo, "", "", p.lunghezza / MM, "", "", "",
            p.x_G / MM, f"{p.m:.6e}", f"{p.J:.6e}", f"{p.I_t:.6e}", "", f"{p.k:.6e}",
            p.segnaposto,
        ])
    return righe


def scrivi_scomposizione(props: Sequence[ProprietaElemento], percorso: str | Path) -> Path:
    """Scrive la tabella di scomposizione per strato (sezione 6)."""
    percorso = Path(percorso)
    percorso.parent.mkdir(parents=True, exist_ok=True)
    with percorso.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(INTESTAZIONE_SCOMPOSIZIONE)
        w.writerows(righe_scomposizione(props))
    return percorso
