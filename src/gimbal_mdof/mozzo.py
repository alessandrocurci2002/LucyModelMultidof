"""Mozzo dell'asse polare e attacco cedevole del telescopio (modello ``hub``).

Conformazione (richiesta dell'utente, 2026-10-17)::

    (tubo + flangia) - (anello di supporto del telescopio) - (flangia + tubo)

* mozzo laterale, corpo unico: flangia piena coassiale all'asse polare, unita all'anello, e tubo
  sede dei cuscinetti verso l'esterno; due mozzi simmetrici;
* anello di supporto: il telescopio vi e' infilato, quindi il suo asse e' quello del telescopio
  (:math:`\\hat K`), perpendicolare all'asse polare; i mozzi sono fissati in due punti
  diametralmente opposti.

Le molle stanno in serie. Tubo e flangia sono alberi in torsione, :math:`k = G J_p / L`.
L'anello, caricato da coppie attorno a :math:`\\hat J` nei due punti d'attacco, lavora come due
semianelli in parallelo; lungo ciascuno il momento e' costante (in direzione :math:`\\hat J`) e si
scompone in flessione fuori piano e torsione della parete. Con l'energia complementare
(Castigliano)

.. math::
    \\theta = \\frac{\\pi R}{4}\\left(\\frac{1}{E I_r} + \\frac{1}{G J_t}\\right) M,
    \\qquad k_{anello} = \\frac{4}{\\pi R}\\left(\\frac{1}{E I_r} + \\frac{1}{G J_t}\\right)^{-1},

con :math:`R` raggio medio, :math:`I_r = t b^3/12` (flessione attorno al raggio),
:math:`J_t = \\beta b t^3` (parete sottile in torsione), :math:`t` spessore radiale, :math:`b`
larghezza. L'anello e' considerato da solo, senza il contributo del tubo del telescopio
(scelta dell'utente, 2026-10-17).

Il centro H del mozzo porta l'inerzia del mozzo attorno all'asse polare (esportazione Onshape).
Le due meta' del mozzo, ``(tubo + flangia + meta' anello)``, sono molle senza massa fra i nodi dei
cuscinetti e H. Il telescopio e' un nodo P fuori dalla catena, collegato a H da :math:`k_b`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from gimbal import inertia as inz
from .geometry import ProprietaElemento
from .params import MM, AttaccoAzimut, ParametriMozzo, ParametriTerne, TrattoMozzo


def _J_polare(d_out: float, d_in: float) -> float:
    return np.pi * (d_out**4 - d_in**4) / 32.0


def rigidezza_tratto(t: TrattoMozzo) -> float:
    """Torsione di un tratto coassiale all'asse polare, :math:`k = G J_p / L`."""
    return t.materiale.G * _J_polare(t.d_out, t.d_in) / t.lunghezza


@dataclass(frozen=True)
class RigidezzeMozzo:
    """Rigidezze delle parti del mozzo e loro serie, N m/rad."""

    k_tubo: float
    k_flangia: float
    k_anello: float
    """Anello intero fra i due punti d'attacco (due semianelli in parallelo)."""
    EI_r: float
    GJ_t: float
    E_anello: float
    R_medio: float
    spessore: float
    k_lato: float
    """Meta' mozzo: tubo + flangia + meta' anello (cedevolezza dell'anello divisa a meta')."""
    k_totale: float
    """Serie completa fra i due cuscinetti."""
    k_semianelli: float = 0.0
    """Anello da solo (modello a semianelli), prima del fattore di irrigidimento."""
    fattore: float = 1.0


def rigidezze(pm: ParametriMozzo) -> RigidezzeMozzo:
    an = pm.anello
    k_tubo = rigidezza_tratto(pm.tubo)
    k_flangia = rigidezza_tratto(pm.flangia)
    t = (an.d_out - an.d_in) / 2.0
    R = (an.d_out + an.d_in) / 4.0
    b = an.larghezza
    E = 2.0 * an.materiale.G * (1.0 + an.poisson)
    EI_r = E * t * b**3 / 12.0
    beta = (1.0 - 0.63 * t / b + 0.052 * (t / b) ** 5) / 3.0
    GJ_t = an.materiale.G * beta * b * t**3
    k_semianelli = 4.0 / (np.pi * R * (1.0 / EI_r + 1.0 / GJ_t))
    k_anello = an.fattore_rigidezza * k_semianelli
    k_lato = 1.0 / (1.0 / k_tubo + 1.0 / k_flangia + 1.0 / (2.0 * k_anello))
    return RigidezzeMozzo(k_tubo=k_tubo, k_flangia=k_flangia, k_anello=k_anello, EI_r=EI_r, GJ_t=GJ_t,
                          E_anello=E, R_medio=R, spessore=t, k_lato=k_lato, k_totale=k_lato / 2.0,
                          k_semianelli=k_semianelli, fattore=an.fattore_rigidezza)


@dataclass(frozen=True)
class InerziaMozzo:
    """Inerzia del mozzo rispetto a o_2 nella terna interna P (x puntamento, y = j_hat, z in alto)."""

    I_2: np.ndarray
    m: float
    J_polare: float
    """Attorno all'asse polare (al nodo H)."""
    J_azimut: float
    """Attorno a z_hat a theta = 0 (si somma a J^ax se il mozzo non e' nel load di azimut)."""
    scarto_baricentro: float
    """Distanza del baricentro da o_2, m."""
    m_geometria: float
    """Massa stimata dalla geometria e dalle densita' dei materiali, kg (controllo)."""
    V_geometria: float


def volume_geometria(pm: ParametriMozzo) -> tuple[float, float]:
    """``(volume, massa)`` dalla geometria nominale: anello + 2 (flangia + tubo)."""
    an = pm.anello
    V_an = np.pi * (an.d_out**2 - an.d_in**2) / 4.0 * an.larghezza
    V_lato = 0.0
    m_lato = 0.0
    for t in (pm.flangia, pm.tubo):
        V = np.pi * (t.d_out**2 - t.d_in**2) / 4.0 * t.lunghezza
        V_lato += V
        m_lato += V * t.materiale.rho
    return V_an + 2.0 * V_lato, V_an * an.materiale.rho + 2.0 * m_lato


def inerzia(pm: ParametriMozzo, tr: ParametriTerne) -> InerziaMozzo:
    """Trasporta il tensore Onshape (baricentro, assi di N) in o_2 e lo ruota nella terna P."""
    r = pm.r_G_N - tr.r_o2
    I_2 = inz.ruota(inz.trasporta(pm.I_G_N, pm.m, r), tr.R_NP_0)
    V, m_geo = volume_geometria(pm)
    return InerziaMozzo(I_2=I_2, m=pm.m, J_polare=float(I_2[1, 1]), J_azimut=float(I_2[2, 2]),
                        scarto_baricentro=float(np.linalg.norm(r)), m_geometria=m_geo, V_geometria=V)


def inerzia_tratto(pm: ParametriMozzo) -> tuple[float, float, float] | None:
    """``(m, J_assiale, J_trasversale)`` del tratto di telescopio solidale all'anello (``None`` se assente).

    Tubo con l'asse lungo il telescopio, centrato in o_2: :math:`J_{ax} = m(r_o^2 + r_i^2)/2` attorno al
    proprio asse (puntamento), :math:`J_t = m[(r_o^2 + r_i^2)/4 + L^2/12]` attorno a un diametro, cioe'
    attorno all'asse polare e a z_hat a theta = 0."""
    t = pm.tratto_telescopio
    if t is None:
        return None
    ro, ri = t.d_out / 2.0, t.d_in / 2.0
    m = t.materiale.rho * np.pi * (ro**2 - ri**2) * t.lunghezza
    return m, m * (ro**2 + ri**2) / 2.0, m * ((ro**2 + ri**2) / 4.0 + t.lunghezza**2 / 12.0)


def lunghezza_lato(pm: ParametriMozzo) -> float:
    """Lunghezza lungo l'asse polare di meta' mozzo (tubo + flangia + raggio esterno dell'anello)."""
    return pm.tubo.lunghezza + pm.flangia.lunghezza + pm.anello.d_out / 2.0


def proprieta_lato(pm: ParametriMozzo, etichetta: str, segnaposto: bool) -> ProprietaElemento:
    """Meta' mozzo come molla senza massa (l'inerzia del mozzo sta tutta nel nodo H)."""
    k = rigidezze(pm).k_lato
    return ProprietaElemento(
        etichetta=etichetta, modo="hub", lunghezza=lunghezza_lato(pm), J_distribuita=0.0, J_extra=0.0,
        k=k, k_segmenti=(k,), k_joint=None, m=0.0, I_t=0.0, x_G=0.0, segmenti=(), strati=(), pezzi=(),
        segnaposto=segnaposto)


def proprieta_attacco(pm: ParametriMozzo) -> ProprietaElemento:
    """Molla k_b fra il centro del mozzo H e il telescopio P (senza massa, lunghezza nulla)."""
    return ProprietaElemento(
        etichetta="b", modo="attachment", lunghezza=0.0, J_distribuita=0.0, J_extra=0.0, k=pm.k_b,
        k_segmenti=(pm.k_b,), k_joint=None, m=0.0, I_t=0.0, x_G=0.0, segmenti=(), strati=(), pezzi=(),
        segnaposto=pm.k_b_segnaposto)


def proprieta_attacco_azimut(at: AttaccoAzimut) -> ProprietaElemento:
    """Molla fra il load di azimut L (forcella + mozzo + anello) e il telescopio P (senza massa)."""
    return ProprietaElemento(
        etichetta="b", modo="attachment", lunghezza=0.0, J_distribuita=0.0, J_extra=0.0, k=at.k,
        k_segmenti=(at.k,), k_joint=None, m=0.0, I_t=0.0, x_G=0.0, segmenti=(), strati=(), pezzi=(),
        segnaposto=at.segnaposto)


def righe_tabella(pm: ParametriMozzo) -> list[tuple[str, str, float]]:
    """``(parte, geometria, k)`` per il report; la geometria e' in mm."""
    rg = rigidezze(pm)
    tu, fl, an = pm.tubo, pm.flangia, pm.anello
    return [
        ("side hub: bearing-seat tube", f"{tu.d_out / MM:g}/{tu.d_in / MM:g} x {tu.lunghezza / MM:g}, "
         f"{tu.materiale.nome}", rg.k_tubo),
        ("side hub: flange", f"{fl.d_out / MM:g} x {fl.lunghezza / MM:g} (solid), {fl.materiale.nome}", rg.k_flangia),
        ("support ring (two half rings)", f"{an.d_out / MM:g}/{an.d_in / MM:g} x {an.larghezza / MM:g}, "
         f"{an.materiale.nome}", rg.k_anello),
    ]
