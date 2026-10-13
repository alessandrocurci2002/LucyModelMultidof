"""Telescopio (payload, corpo 2) e termine :math:`J^{ax}`.

Terne (convenzione della tesi De Donatis, :class:`gimbal_mdof.params.ParametriTerne`):
inerziale N (origine :math:`O_N`, versori :math:`\\hat x, \\hat y, \\hat z`), corpo 1
(azimut, origine :math:`o_1`, versori :math:`\\hat\\imath, \\hat\\jmath, \\hat k`), corpo 2
(polare, origine :math:`o_2`, versori :math:`\\hat I, \\hat J, \\hat K`, con
:math:`\\hat K` lungo il puntamento verso le ottiche, :math:`\\hat J` asse polare,
:math:`\\hat I = \\hat J\\times\\hat K`).

La dinamica usa internamente la terna di puntamento P del modello a 2 GdL
(:math:`x_P = \\hat K`, :math:`y_P = \\hat J`, :math:`z_P = -\\hat I`): ``Payload.I`` e
``Payload.d`` sono in P, e gli attributi ``I_xx, I_yy, I_zz, ...`` sono le sue
componenti. ``Payload.I_2`` e ``Payload.d_2`` danno le stesse grandezze nella terna 2
del documento (:math:`I_{KK} = I_{xx}`, :math:`I_{II} = I_{zz}`, :math:`d_K = d_x`,
:math:`d_I = -d_z`, ...), tramite :data:`gimbal_mdof.params.R_P2`.

Due sorgenti (``payload.source``):

* ``manual``: dati da ``config.yaml``. Con ``frame: cad_cg`` il tensore e' dato
  rispetto al baricentro con assi paralleli a N (l'esportazione di Onshape): lo si
  porta in convenzione tensoriale, lo si trasporta in :math:`o_2` (Huygens-Steiner)
  e lo si ruota nella terna 2 con la posa CAD di ``frames``, usando le funzioni del
  modello a 2 GdL (**riusate, non duplicate**). Con ``frame: body2_o2`` il tensore
  e' gia' quello richiesto.
* ``gimbal_2dof``: corpo B del CAD del modello a 2 GdL.

Se il tensore comprende anche parti dell'asse polare gia' modellate negli elementi
FEM, ``subtract_shaft_parts: true`` le sottrae: ogni pezzo, assialsimmetrico
attorno a :math:`\\hat J`, ha tensore rispetto a :math:`o_2` in terna 2
:math:`\\operatorname{diag}(I_t + m y_c^2,\\ J,\\ I_t + m y_c^2)`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Final

import numpy as np
from gimbal import inertia as inz
from gimbal import params as pr2
from gimbal.inertia import ModelloInerziale

from .fem import ModelloAsse
from .params import R_P2, Avviso, Parametri

AVVISI_2DOF_ESCLUSI: Final[frozenset[str]] = frozenset(
    {"SMORZAMENTI_NULLI", "OMEGA0_SEGNAPOSTO", "FLEX_PARAMETRI_MANCANTI", "AZ_AXIS_XY"}
)
"""Avvisi del 2 GdL che non riguardano il telescopio (hanno un analogo qui)."""



@dataclass(frozen=True)
class Payload:
    """Telescopio rigido collegato al nodo P, e contributo trasversale polare.

    Attributes
    ----------
    I
        Tensore :math:`\\mathbf I_{p,o_2}|_2` del telescopio, kg m^2.
    m, d
        Massa ed eccentricita' del baricentro rispetto a :math:`o_2`, in terna 2.
    J_ax
        :math:`J^{ax} = \\sum (I_t + m y_c^2)` sui pezzi polari, kg m^2.
    sottratto
        Totali sottratti dal tensore inserito (massa, :math:`\\sum J`,
        :math:`\\sum(I_t + m y_c^2)`), per la tracciabilita'.
    par_2dof
        Parametri del modello a 2 GdL, usati come contenitore nel confronto del
        limite rigido.
    """

    abilitato: bool
    I: np.ndarray
    m: float
    d: np.ndarray
    g: float
    J_ax: float
    sottratto: dict[str, float]
    par_2dof: pr2.Parametri | None
    avvisi: tuple[Avviso, ...] = field(default_factory=tuple)

    @property
    def I_2(self) -> np.ndarray:
        """Tensore nella terna 2 del documento (I_hat, J_hat, K_hat), K_hat = puntamento."""
        return R_P2.T @ self.I @ R_P2

    @property
    def d_2(self) -> np.ndarray:
        """Eccentricita' nella terna 2 del documento: (d_I, d_J, d_K)."""
        return R_P2.T @ self.d

    @property
    def I_xx(self) -> float:
        return float(self.I[0, 0])

    @property
    def I_yy(self) -> float:
        return float(self.I[1, 1])

    @property
    def I_zz(self) -> float:
        return float(self.I[2, 2])

    @property
    def I_xy(self) -> float:
        return float(self.I[0, 1])

    @property
    def I_xz(self) -> float:
        return float(self.I[0, 2])

    @property
    def I_yz(self) -> float:
        return float(self.I[1, 2])

    @property
    def d_x(self) -> float:
        return float(self.d[0])

    @property
    def d_y(self) -> float:
        return float(self.d[1])

    @property
    def d_z(self) -> float:
        return float(self.d[2])

    # -- grandezze del punto di lavoro (sezioni 9.5 e 10.2) ------------------

    # Convenzione: theta > 0 alza il puntamento. Componenti nella terna interna P
    # (x_P puntamento, y_P = j_hat = -J_hat, z_P verso l'alto a theta = 0).

    def J_psi(self, theta: float | np.ndarray) -> float | np.ndarray:
        """:math:`J_{\\psi,p}(\\theta) = I_{xx}\\sin^2\\theta + I_{zz}\\cos^2\\theta + I_{xz}\\sin 2\\theta` (terna P)."""
        return (self.I_xx * np.sin(theta) ** 2 + self.I_zz * np.cos(theta) ** 2
                + self.I_xz * np.sin(2 * theta))

    def J_psi_prime(self, theta: float | np.ndarray) -> float | np.ndarray:
        return (self.I_xx - self.I_zz) * np.sin(2 * theta) + 2 * self.I_xz * np.cos(2 * theta)

    def J_psi_second(self, theta: float | np.ndarray) -> float | np.ndarray:
        return 2 * (self.I_xx - self.I_zz) * np.cos(2 * theta) - 4 * self.I_xz * np.sin(2 * theta)

    def M_12(self, theta: float | np.ndarray) -> float | np.ndarray:
        return -(self.I_xy * np.sin(theta) + self.I_yz * np.cos(theta))

    def V_g(self, theta: float | np.ndarray) -> float | np.ndarray:
        """:math:`m_p g\\,(d_x\\sin\\theta + d_z\\cos\\theta)`: quota del baricentro sopra :math:`o_2`."""
        return self.m * self.g * (self.d_x * np.sin(theta) + self.d_z * np.cos(theta))

    def dV_dtheta(self, theta: float | np.ndarray) -> float | np.ndarray:
        return self.m * self.g * (self.d_x * np.cos(theta) - self.d_z * np.sin(theta))

    def k_g(self, theta: float | np.ndarray) -> float | np.ndarray:
        """Rigidezza gravitazionale :math:`-m_p g (d_x\\sin\\theta + d_z\\cos\\theta)`."""
        return -self.m * self.g * (self.d_x * np.sin(theta) + self.d_z * np.cos(theta))

    def k_c(self, theta: float, Omega: float) -> float:
        """Rigidezza centrifuga :math:`-(I_{xx}-I_{zz})\\Omega_0^2\\cos 2\\theta_0`."""
        return float(-0.5 * self.J_psi_second(theta) * Omega**2)

    def gamma(self, theta: float, Omega: float) -> float:
        """Coefficiente giroscopico :math:`\\gamma = J_{\\psi,p}'(\\theta_0)\\,\\Omega_0`."""
        return float(self.J_psi_prime(theta) * Omega)

    def prodotti_nulli(self) -> "Payload":
        """Copia con i prodotti d'inerzia nulli (target di regressione)."""
        return replace(self, I=np.diag(np.diag(self.I)))


def posizioni_polari(mod: ModelloAsse, y_nodo1: float | None,
                     verso: float = 1.0) -> list[tuple[str, float, float, float, float]]:
    """Pezzi polari ``(nome, m, J, I_t, y_c)``, con :math:`y_c` da O lungo :math:`y_A`.

    :math:`y_c = y_1 + s\\,x`, con :math:`x` l'ascissa del baricentro del pezzo dal
    nodo 1 e :math:`s = \\pm 1` il verso in cui cresce la numerazione dei nodi
    (``nodes_increase_along``). Se ``y_nodo1`` e' ``None`` (Q7) si pone
    :math:`y_c = 0`: nessuno Steiner.
    """
    return [
        (pp.pezzo.nome, pp.pezzo.m, pp.pezzo.J, pp.pezzo.I_t,
         0.0 if y_nodo1 is None else y_nodo1 + verso * pp.x)
        for pp in mod.pezzi
    ]


def J_ax(mod: ModelloAsse, y_nodo1: float | None, verso: float = 1.0) -> float:
    """:math:`J^{ax} = \\sum (I_t + m y_c^2)` sui pezzi dell'asse polare (sezione 9.2)."""
    return float(sum(I_t + m * y**2 for _, m, _, I_t, y in posizioni_polari(mod, y_nodo1, verso)))


def tensore_da_config(par: Parametri) -> tuple[np.ndarray, float, np.ndarray, list[Avviso]]:
    """:math:`(\\mathbf I_{p,o_2}|_2, m_p, \\mathbf d_2)` dai dati di ``config.yaml``.

    ``frame: cad_cg`` (come Onshape): :math:`\\mathbf I_{p,G}|_N` rispetto al
    baricentro, assi paralleli a N, baricentro :math:`\\mathbf r_G|_N` misurato da
    :math:`O_N`. Con :math:`\\mathbf r = \\mathbf r_G - \\mathbf r_{o_2}` e l'assetto
    :math:`\\mathbf R_{N2}^{(0)}` del corpo 2 nella posa CAD:

    .. math:: \\mathbf I_{p,o_2}|_N = \\mathbf I_{p,G}|_N + m(\\|\\mathbf r\\|^2\\mathbb 1 - \\mathbf r\\mathbf r^{\\mathsf T}),
              \\quad \\mathbf I_{p,o_2}|_2 = \\mathbf R_{N2}^{(0)\\mathsf T}\\mathbf I_{p,o_2}|_N\\mathbf R_{N2}^{(0)},
              \\quad \\mathbf d_2 = \\mathbf R_{N2}^{(0)\\mathsf T}\\mathbf r.

    ``frame: body2_o2``: tensore ed eccentricita' gia' rispetto a :math:`o_2` in terna 2.
    """
    pp = par.payload
    tr = par.terne
    assert pp.I is not None and pp.m is not None
    avvisi: list[Avviso] = []
    m = float(pp.m)
    I_in = inz.normalizza_prodotti(pp.I, pp.segno_prodotti)
    if pp.terna == "cad_cg":
        assert pp.r_G_N is not None
        inz.controlli_fisici(I_in, "I_p,G|N (payload in config.yaml)")
        r = pp.r_G_N - tr.r_o2
        I_2 = inz.ruota(inz.trasporta(I_in, m, r), tr.R_NP_0)
        d = tr.R_NP_0.T @ r
        I_G_2 = inz.ruota(I_in, tr.R_NP_0)
        avvisi.append(Avviso(
            "PAYLOAD_SEGNO_PRODOTTI", "segnaposto",
            f"payload.products_sign = '{pp.segno_prodotti}': la convenzione di segno dei prodotti "
            "d'inerzia di Onshape non e' verificata con un corpo di prova; influenza i prodotti "
            "I_IJ, I_IK, I_JK, non i momenti diagonali."))
    else:
        assert pp.d is not None
        # input nella terna 2 (K_hat puntamento): si passa alla terna interna P
        I_2 = R_P2 @ I_in @ R_P2.T
        d = R_P2 @ np.asarray(pp.d, dtype=float)
        I_G_2 = inz.trasporta(I_2, -m, d)
    inz.controlli_fisici(I_2, "I_p,o2|2 (payload in config.yaml)")

    # Il tubo ottico ruota attorno al proprio asse: l'asse di puntamento I_hat deve
    # essere quello di inerzia baricentrica minima. Se non lo e', la terna di input
    # o la posa CAD sono probabilmente sbagliate: lo si segnala, non si corregge.
    diag = np.diag(I_G_2)
    if int(np.argmin(diag)) != 0:
        avvisi.append(Avviso(
            "PAYLOAD_ASSE_SOSPETTO", "incoerenza",
            f"Momenti baricentrici del telescopio in terna 2: {np.array2string(diag, precision=5)} "
            "kg m^2 (lungo K_hat, J_hat, -I_hat). Il minimo non e' attorno a K_hat (puntamento): "
            "verificare payload.frame, gli "
            "assi del tensore inserito e la posa CAD (frames.pointing_cad, theta_cad_deg)."))
    return I_2, m, np.asarray(d, dtype=float), avvisi


def _sottrai_pezzi(I: np.ndarray, m: float, d: np.ndarray,
                   pezzi: list[tuple[str, float, float, float, float]]) -> tuple[np.ndarray, float, np.ndarray, dict[str, float]]:
    """Sottrae i pezzi assialsimmetrici dell'asse polare dal tensore rispetto a o_2 in terna 2."""
    S = np.zeros((3, 3))
    primo = np.zeros(3)
    sottratto = {"m": 0.0, "J": 0.0, "trasversale": 0.0}
    for _, m_i, J_i, It_i, y_i in pezzi:
        S += np.diag([It_i + m_i * y_i**2, J_i, It_i + m_i * y_i**2])
        primo += m_i * np.array([0.0, y_i, 0.0])
        sottratto["m"] += m_i
        sottratto["J"] += J_i
        sottratto["trasversale"] += It_i + m_i * y_i**2
    m_p = m - sottratto["m"]
    if m_p <= 0.0:
        raise ValueError(f"massa del telescopio non positiva dopo la sottrazione degli alberi "
                         f"polari ({m_p:.4g} kg)")
    I_p = I - S
    try:
        inz.controlli_fisici(I_p, "I_p,o2|2 (tensore inserito meno alberi polari)")
    except inz.TensoreNonFisico as e:
        raise ValueError(
            f"{e}. La sottrazione degli elementi polari da' un tensore non fisico: verificare la "
            "geometria polare o porre payload.subtract_shaft_parts: false.") from e
    return I_p, m_p, (m * d - primo) / m_p, sottratto


def costruisci_payload(par: Parametri, mod_polare: ModelloAsse) -> Payload:
    """Telescopio collegato al nodo P, dalla sorgente scelta in ``config.yaml``."""
    avvisi: list[Avviso] = []
    pp = par.payload
    pezzi = posizioni_polari(mod_polare, par.polare.y_nodo1_da_O, par.polare.verso_numerazione)
    Jax = float(sum(I_t + m * y**2 for _, m, _, I_t, y in pezzi))
    par2 = pr2.parametri_da_config(pp.config_2dof) if pp.config_2dof.exists() else None

    if not pp.abilitato:
        avvisi.append(Avviso("PAYLOAD_DISABILITATO", "informazione",
                             "payload.enabled = false: nessun telescopio al nodo P."))
        return Payload(abilitato=False, I=np.zeros((3, 3)), m=0.0, d=np.zeros(3), g=par.g,
                       J_ax=Jax, sottratto={}, par_2dof=par2, avvisi=tuple(avvisi))

    if pp.sorgente == "manual":
        I, m, d, avv = tensore_da_config(par)
        avvisi.extend(avv)
        origine = (f"tensore inserito in config.yaml (frame: {pp.terna}"
                   + (f", {pp.nota}" if pp.nota else "") + ")")
    else:
        if par2 is None:
            raise ValueError(f"payload.source = gimbal_2dof ma {pp.config_2dof} non esiste")
        mi = inz.costruisci_modello_inerziale(par2)
        avvisi.extend(Avviso(f"2DOF_{a.codice}", a.severita, f"(dal modello a 2 GdL) {a.messaggio}")
                      for a in mi.avvisi if a.codice not in AVVISI_2DOF_ESCLUSI)
        I, m, d = mi.I_B_O_B.copy(), float(mi.m_B), mi.d.copy()
        origine = "corpo B del CAD del modello a 2 GdL"
        scarto = abs(float(mi.r_O_W[2]) - par.terne.z_o2)
        if scarto > 1e-6:
            avvisi.append(Avviso(
                "PAYLOAD_O2_2DOF", "incoerenza",
                f"Con source: gimbal_2dof il tensore e' trasportato nel punto O del 2 GdL "
                f"(z = {mi.r_O_W[2] * 1e3:.3f} mm), non in o_2 di frames (z = "
                f"{par.terne.z_o2 * 1e3:.3f} mm): scarto {scarto * 1e3:.3f} mm."))

    sottratto = {"m": 0.0, "J": 0.0, "trasversale": 0.0}
    if pp.sottrai_alberi:
        m_in = m
        I, m, d, sottratto = _sottrai_pezzi(I, m, d, pezzi)
        avvisi.append(Avviso(
            "Q3_SOTTRAZIONE_ALBERI", "segnaposto",
            f"Q3: dal {origine} (m = {m_in:.3f} kg) sono stati sottratti gli elementi e le "
            f"inerzie nodali polari del FEM: massa {sottratto['m']:.4f} kg, sum J = "
            f"{sottratto['J']:.3e} kg m^2, sum(I_t + m y_c^2) = {sottratto['trasversale']:.3e} "
            "kg m^2. Si assume che il tensore inserito contenga esattamente questi pezzi."))
    else:
        avvisi.append(Avviso(
            "Q3_TELESCOPIO", "informazione",
            f"Q3: telescopio dal {origine}, senza sottrazioni (subtract_shaft_parts: false): "
            "il tensore deve comprendere il solo telescopio, non le parti dell'asse polare "
            "gia' modellate negli elementi."))

    return Payload(abilitato=True, I=I, m=m, d=np.asarray(d, dtype=float), g=par.g, J_ax=Jax,
                   sottratto=sottratto, par_2dof=par2, avvisi=tuple(avvisi))


def modello_2dof_equivalente(pl: Payload, J_azimut: float, J_polare: float) -> ModelloInerziale:
    """Modello rigido a 2 GdL con le inerzie totali del modello multi-GdL (test 16).

    :math:`J_A = J_{azimut}` (alberi, nodi, forcella e :math:`J^{ax}`), tensore del
    telescopio con :math:`I_{yy} + J_{polare}` (alberi polari), stessi
    prodotti, massa e offset del telescopio. Il preprocessing del 2 GdL serve
    solo da contenitore, i cui campi inerziali vengono sostituiti.
    """
    assert pl.par_2dof is not None
    base = inz.costruisci_modello_inerziale(pl.par_2dof, verifica=False)
    I = pl.I.copy()
    I[1, 1] += J_polare
    return replace(base, J_A=float(J_azimut), I_B_O_B=I, d=pl.d.copy(), m_B=pl.m, g=pl.g)
