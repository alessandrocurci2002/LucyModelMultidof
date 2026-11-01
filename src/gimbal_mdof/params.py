"""Lettura e validazione di ``config.yaml``, conversione in SI (sezione 5).

E' l'unico punto del pacchetto in cui entrano i dati configurabili e in cui
avviene la conversione delle unita' dei file di input (mm, kg m^2, N m/rad)
in SI. Tutto il resto del pacchetto riceve un oggetto :class:`Parametri`
immutabile, gia' validato.

Ogni dato segnaposto o ipotesi di default produce un
:class:`gimbal.params.Avviso` (riusato dal modello a 2 GdL), che viene
propagato fino a ``report.md``, alle figure e al documento LaTeX: la sezione 0
della specifica vieta di presentare un risultato da parametri segnaposto senza
avvertenza.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Final, Literal, Mapping, Sequence

import numpy as np
import yaml
from gimbal.params import Avviso, Severita

__all__ = [
    "Avviso",
    "ErroreConfigurazione",
    "Materiale",
    "Strato",
    "Segmento",
    "InerziaExtra",
    "Elemento",
    "InerziaNodale",
    "ParametriAsse",
    "ParametriPayload",
    "ParametriModali",
    "ParametriRiduzione",
    "ParametriSperimentali",
    "VoceMisura",
    "ParametriFEM",
    "Parametri",
    "parametri_da_config",
    "RADICE",
]

MM: Final[float] = 1e-3
"""Da mm a m: l'unico fattore di conversione delle lunghezze del pacchetto."""

RADICE: Final[Path] = Path(__file__).resolve().parents[2]
"""Radice del progetto ``gimbal_multidof/``."""

NOTA_REALISMO: Final[str] = "result from placeholder parameters"
"""Avviso da riportare su ogni figura, tabella o frequenza da dati segnaposto (in inglese,
come tutti i testi del report e delle figure)."""

NOME_ASSE: Final[dict[str, str]] = {"azimut": "azimuth", "polare": "polar"}
"""Nome inglese degli assi per report e figure (le chiavi restano quelle del codice)."""

NOME_MATERIALE: Final[dict[str, str]] = {"acciaio": "steel", "alluminio": "aluminium"}
"""Nome inglese dei materiali di ``config.yaml`` per report e figure."""


def in_inglese(nome: str) -> str:
    """Nome di asse o materiale da mostrare nel report (invariato se non mappato)."""
    return NOME_ASSE.get(nome, NOME_MATERIALE.get(nome, nome))

Forma = Literal["solid", "hollow"]
ModoElemento = Literal["geometry", "manual", "telescope"]
NomeAsse = Literal["azimut", "polare"]


class ErroreConfigurazione(ValueError):
    """Il file di configurazione non rispetta lo schema della sezione 5."""


# ---------------------------------------------------------------------------
# Dati degli elementi (SI)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Materiale:
    """Materiale omogeneo e isotropo."""

    nome: str
    rho: float
    """Densita', kg/m^3."""
    G: float
    """Modulo di taglio, Pa."""
    segnaposto: bool


@dataclass(frozen=True)
class Strato:
    """Strato anulare coassiale di un segmento (sezione 6).

    Per il cilindro pieno ``d_in = 0``.
    """

    forma: Forma
    d_out: float
    """Diametro esterno, m."""
    d_in: float
    """Diametro interno, m (0 per il pieno)."""
    materiale: Materiale
    portante: bool
    """``True`` se lo strato porta coppia e contribuisce alla rigidezza."""


@dataclass(frozen=True)
class Segmento:
    """Tratto di elemento con sezione costante; i segmenti sono in serie.

    ``rigido = True`` marca un tratto a cedevolezza nulla (morsetto del
    telescopio, Q2): non viene mai letto da ``config.yaml``, lo imposta
    :mod:`gimbal_mdof.fem` quando divide l'elemento del telescopio.
    """

    lunghezza: float
    strati: tuple[Strato, ...]
    rigido: bool = False


@dataclass(frozen=True)
class InerziaExtra:
    """Inerzia non distribuita associata a un elemento (es. disco encoder).

    Si divide a meta' sui due nodi dell'elemento in entrambi i modelli di massa.
    """

    nome: str
    J: float
    m: float | None
    I_t: float | None
    segnaposto: bool


@dataclass(frozen=True)
class Elemento:
    """Elemento ``i`` dell'albero, dal nodo ``i`` al nodo ``i+1``.

    In modalita' ``geometry`` inerzia e rigidezza vengono da
    :mod:`gimbal_mdof.geometry`; ``k_manuale``, se dato, sostituisce la sola
    rigidezza. In modalita' ``manual`` valgono ``J_manuale`` e ``k_manuale``.
    """

    id: int
    modo: ModoElemento
    lunghezza: float
    segmenti: tuple[Segmento, ...]
    k_joint: float | None
    extra: tuple[InerziaExtra, ...]
    J_manuale: float | None
    k_manuale: float | None
    m_manuale: float | None
    I_t_manuale: float | None
    y_c_manuale: float | None
    """Baricentro dal nodo ``i``, m (solo modalita' manual)."""
    segnaposto: bool


@dataclass(frozen=True)
class InerziaNodale:
    """Inerzia concentrata a un nodo (rotori, forcella).

    Con ``frame: cad_cg`` (esportazione Onshape: tensore rispetto al baricentro,
    assi paralleli a N, baricentro da :math:`O_N`) ``J`` e' il momento attorno
    all'asse di azimut, :math:`J = I_{zz,G} + m\\,[(x_G - x_a)^2 + (y_G - y_a)^2]`
    (Huygens-Steiner), e ``J_G`` il momento baricentrico di partenza.
    """

    nodo: int
    nome: str
    J: float
    m: float | None
    I_t: float | None
    fonte: str
    segnaposto: bool
    J_G: float | None = None
    """Momento baricentrico :math:`I_{zz,G}` dell'esportazione CAD (solo ``frame: cad_cg``)."""
    distanza_asse: float | None = None
    """Distanza del baricentro dall'asse di azimut, m (solo ``frame: cad_cg``)."""
    contiene_alberi_polari: bool = False
    """``True`` se l'esportazione comprende gli alberi polari (non il telescopio): il loro
    contributo attorno all'asse di azimut, :math:`J^{ax}`, e' gia' qui e non va aggiunto."""


@dataclass(frozen=True)
class TrattoMozzo:
    """Tratto assialsimmetrico di un mozzo laterale, coassiale all'asse polare (SI)."""

    d_out: float
    d_in: float
    lunghezza: float
    materiale: Materiale


@dataclass(frozen=True)
class AnelloMozzo:
    """Anello di supporto del telescopio: asse lungo il telescopio (K_hat), perpendicolare
    all'asse polare; i due mozzi laterali vi sono fissati in punti diametralmente opposti (SI)."""

    d_out: float
    d_in: float
    larghezza: float
    """Estensione lungo l'asse dell'anello (``width_mm``)."""
    materiale: Materiale
    poisson: float
    fattore_rigidezza: float = 1.0
    """``stiffness_factor``: moltiplica la rigidezza dei semianelli (1 = anello da solo, limite inferiore)."""


@dataclass(frozen=True)
class ParametriMozzo:
    """Mozzo dell'asse polare e attacco del telescopio (modello ``hub``, richiesta 2026-10-17).

    Il mozzo (due mozzi laterali tubo + flangia, anello di supporto) sta sul percorso fra i
    cuscinetti: le sue due meta' sostituiscono gli elementi ``elementi_lato`` e il suo centro H
    e' il nodo dell'elemento ``elemento_centro`` (nodi coincidenti). Il telescopio non e' piu'
    rigido in H: e' un nodo P fuori dalla catena, collegato a H dalla molla ``k_b``.
    """

    elementi_lato: tuple[int, int]
    elemento_centro: int
    m: float
    r_G_N: np.ndarray
    """Baricentro in N da O_N, m (esportazione Onshape senza mate connector)."""
    I_G_N: np.ndarray
    """Tensore rispetto al baricentro, assi paralleli a N, kg m^2."""
    nel_load_azimut: bool
    """``True`` se l'inerzia del mozzo attorno a z_hat e' gia' nel load di azimut."""
    fonte_inerzia: str
    inerzia_segnaposto: bool
    tubo: TrattoMozzo
    flangia: TrattoMozzo
    anello: AnelloMozzo
    geometria_segnaposto: bool
    k_b: float
    fonte_k_b: str
    k_b_segnaposto: bool
    tratto_telescopio: TrattoMozzo | None = None
    """``telescope_on_ring``: tratto del tubo del telescopio stretto nell'anello, solidale al mozzo
    (asse lungo il telescopio, centrato in o_2); la sua inerzia passa dal telescopio al nodo H."""
    tratto_segnaposto: bool = False
    J_polare: float = 0.0
    """Inerzia del mozzo attorno all'asse polare, kg m^2 (da ``I_G_N`` e dalle terne; comprende il
    tratto di telescopio solidale all'anello, se presente)."""
    J_azimut: float = 0.0
    """Inerzia del mozzo attorno a z_hat a theta = 0, kg m^2."""


@dataclass(frozen=True)
class AttaccoAzimut:
    """Attacco cedevole del telescopio sull'asse di azimut (``azimuth.telescope_attachment``, 2026-10-25).

    Il load di azimut non e' piu' un corpo unico: (forcella + mozzo + anello) ruota con il nodo L e il
    telescopio e' un nodo P fuori dalla catena, collegato a L dalla molla ``k`` attorno a z_hat (a theta = 0
    e' il ribaltamento del tubo nell'anello attorno a I_hat, cioe' lo stesso giunto a viti di k_b).
    """

    k: float
    fonte: str
    segnaposto: bool


@dataclass(frozen=True)
class ParametriAsse:
    """Dati di un asse, con nodi ed elementi numerati da 1."""

    nome: NomeAsse
    n_nodi: int
    nodo_encoder: int
    nodo_coppia: int
    nodo_estremo: int
    """``load_node`` per l'azimut, ``end_node`` per il polare."""
    masters: tuple[int | str, ...]
    elementi: tuple[Elemento, ...]
    inerzie_nodali: tuple[InerziaNodale, ...]
    attrito_cuscinetti: tuple[tuple[int, float], ...]
    """Coppie ``(elemento, b)`` in N m s/rad."""
    elemento_payload: int | None = None
    xi_payload: float | None = None
    tratto_rigido: float = 0.0
    elementi_cuscinetto: tuple[int, ...] = ()
    y_nodo1_da_O: float | None = None
    verso_numerazione: float = 1.0
    """+1 se l'indice dei nodi polari cresce lungo +y_A, -1 se lungo -y_A."""
    mozzo: ParametriMozzo | None = None
    """Mozzo e attacco cedevole del telescopio (solo polare, ``polar.hub.enabled: true``)."""
    attacco_telescopio: AttaccoAzimut | None = None
    """Telescopio appeso al load con una molla (solo azimut, ``azimuth.telescope_attachment.enabled: true``)."""

    @property
    def n_elementi(self) -> int:
        return self.n_nodi - 1

    def elemento(self, id_: int) -> Elemento:
        """Elemento con identificativo ``id_`` (base 1)."""
        return self.elementi[id_ - 1]

    @property
    def b_cuscinetti(self) -> dict[int, float]:
        return dict(self.attrito_cuscinetti)

    @property
    def segnaposto(self) -> bool:
        return any(e.segnaposto for e in self.elementi) or any(
            n.segnaposto for n in self.inerzie_nodali
        )


@dataclass(frozen=True)
class ParametriPayload:
    """Telescopio rigido collegato al nodo P (SI).

    ``sorgente = "manual"``: tensore, massa e baricentro da ``config.yaml``,
    in una delle due terne di input:

    * ``"cad_cg"``: tensore rispetto al baricentro con assi paralleli alla
      terna inerziale N (x_hat, y_hat, z_hat), che coincide con la terna
      dell'assieme Onshape (esportazione *Mass properties* senza mate
      connector); baricentro ``r_G_N`` in N, misurato da O_N;
    * ``"body2_o2"``: tensore rispetto a o_2 nella terna del corpo 2
      (I_hat, J_hat, K_hat), eccentricita' del baricentro ``d`` in terna 2.

    ``sorgente = "gimbal_2dof"``: corpo B del CAD del modello a 2 GdL.
    """

    abilitato: bool
    sorgente: Literal["manual", "gimbal_2dof"]
    config_2dof: Path
    sottrai_alberi: bool
    terna: Literal["cad_cg", "body2_o2"] | None = None
    m: float | None = None
    I: np.ndarray | None = None
    """Tensore inserito, kg m^2, con i fuori diagonale come scritti nel file."""
    r_G_N: np.ndarray | None = None
    d: np.ndarray | None = None
    segno_prodotti: Literal["tensor", "positive"] = "tensor"
    nota: str = ""
    segnaposto: bool = False


R_P2: Final[np.ndarray] = np.array([[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]])
"""Colonne: I_hat, J_hat, K_hat nella terna interna di puntamento P (x_P puntamento,
y_P = j_hat, z_P = x_P x y_P verso l'alto a theta = 0): I_hat = z_P, J_hat = -y_P,
K_hat = x_P. Rotazione propria (det = +1)."""


@dataclass(frozen=True)
class ParametriTerne:
    """Terne di riferimento (convenzione della tesi De Donatis), in SI.

    * N, inerziale: origine O_N, versori x_hat, y_hat, z_hat; coincide con la
      terna dell'assieme Onshape, z_hat verticale verso l'alto lungo l'asse
      di azimut.
    * Corpo 1 (azimut, forcella): origine o_1 sull'asse di azimut, versori
      i_hat, j_hat, k_hat = z_hat; ruota con psi_L attorno a k_hat. j_hat e'
      l'asse polare.
    * Corpo 2 (polare, telescopio): origine o_2 = centro di rotazione dell'albero
      polare sull'asse di azimut, versori K_hat (puntamento, verso le ottiche),
      I_hat (verso l'alto con il telescopio orizzontale), J_hat = K_hat x I_hat
      (asse polare, J_hat = -j_hat); ruota con theta_p attorno a J_hat, quindi
      theta > 0 ALZA il puntamento (theta = elevazione).

    Internamente la dinamica usa la terna "di puntamento" P del modello a 2 GdL
    (x_P = K_hat, y_P = j_hat = -J_hat, z_P = I_hat), che ruota di -theta attorno a
    y_P: :data:`R_P2` converte le componenti fra P e la terna 2.

    O_N, o_1, o_2 stanno su z_hat e differiscono solo per una traslazione lungo
    z_hat. Le origini dei corpi sono sugli assi di rotazione, non nei
    baricentri: l'eccentricita' del baricentro e' un dato a parte.
    """

    asse_az_xy: np.ndarray
    """Posizione dell'asse di azimut nel piano (x_hat, y_hat), m."""
    z_o1: float
    """dist(O_N, o_1) lungo z_hat, m."""
    z_o2: float
    """dist(O_N, o_2) lungo z_hat, m."""
    i_hat_0: np.ndarray
    """Direzione orizzontale del puntamento nella posa CAD, in N (i_hat a psi = 0)."""
    J_hat_0: np.ndarray
    """Asse polare J_hat della terna 2 nella posa CAD, in N (J_hat = -j_hat)."""
    theta_cad: float
    """Elevazione della posa CAD (theta > 0 alza il puntamento), rad."""
    R_N1_0: np.ndarray
    """Assetto del corpo 1 nella posa CAD: colonne i_hat, j_hat, k_hat in N."""
    R_N2_0: np.ndarray
    """Assetto del corpo 2 nella posa CAD: colonne I_hat, J_hat, K_hat in N."""
    R_NP_0: np.ndarray
    """Assetto della terna interna di puntamento nella posa CAD (x_P = K_hat)."""

    @property
    def r_o1(self) -> np.ndarray:
        return np.array([self.asse_az_xy[0], self.asse_az_xy[1], self.z_o1])

    @property
    def r_o2(self) -> np.ndarray:
        return np.array([self.asse_az_xy[0], self.asse_az_xy[1], self.z_o2])


MetodoRiduzione = Literal["guyan", "craig_bampton", "modal"]
CorrezioneStatica = Literal["none", "residual_vectors"]


@dataclass(frozen=True)
class ParametriModali:
    """Troncamento modale alla Gawronski (``reduction.modal``), unita' SI.

    Attributes
    ----------
    f_banda
        Banda di interesse, Hz.
    fattore_taglio
        Si tengono i modi flessibili con :math:`f_i \\le` ``fattore_taglio * f_banda``.
    hankel_min
        Se non ``None``: si tengono anche i modi sopra soglia con norma di Hankel
        :math:`\\ge` ``hankel_min`` volte la massima dell'asse.
    n_flessibili
        Numero di modi flessibili forzato per asse (``None`` = criterio automatico).
    theta_ref
        :math:`\\theta_{ref}` della base delle simulazioni non lineari, rad.
    correzione_statica
        ``none`` oppure ``residual_vectors`` (vettori statici residui).
    """

    f_banda: float = 2000.0
    fattore_taglio: float = 2.0
    hankel_min: float | None = None
    n_flessibili: Mapping[str, int | None] = field(
        default_factory=lambda: {"azimut": None, "polare": None})
    theta_ref: float = 0.0
    correzione_statica: CorrezioneStatica = "none"

    @property
    def f_taglio(self) -> float:
        return self.fattore_taglio * self.f_banda


@dataclass(frozen=True)
class ParametriRiduzione:
    metodo: MetodoRiduzione
    n_modi_interfaccia: int
    forma_massa: Literal["guyan", "guyan_lumped"]
    modale: ParametriModali = field(default_factory=ParametriModali)


@dataclass(frozen=True)
class VoceMisura:
    """Una misura sperimentale da confrontare con il modello (sezione ``experimental``)."""

    file: str
    asse: NomeAsse
    theta0: float
    """Elevazione durante la prova, rad (theta > 0 alza il puntamento)."""
    telescopio: bool
    etichetta: str
    theta0_segnaposto: bool = False


@dataclass(frozen=True)
class ParametriSperimentali:
    """Misure di risposta in frequenza (Kollmorgen WorkBench) per il confronto con il modello."""

    cartella: Path
    Kt: float
    """Costante di coppia del motore, N m/A (corrente delle misure)."""
    traccia: str
    f_min: float
    f_max: float
    banda_livello: tuple[float, float]
    misure: tuple[VoceMisura, ...]


@dataclass(frozen=True)
class VincoloTerra:
    """Rigidezza equivalente verso terra di un asse (sezione ``ground_constraint``).

    Equivalente lineare, nel dominio della frequenza, dell'attrito di pre-scorrimento dei
    cuscinetti sotto pochi hertz: coppia :math:`k(1 + j\eta)\,\theta + c\,\dot\theta` sul nodo
    dei cuscinetti (azimut: verso la base; polare: verso la forcella). Vale per l'ampiezza
    della prova su cui e' stimato (l'attrito non e' lineare) e non entra nel simulatore nel tempo.
    """

    k: float
    """Rigidezza, N m/rad."""
    eta: float = 0.0
    """Fattore di perdita isteretico (perdita indipendente dalla frequenza)."""
    c: float = 0.0
    """Smorzamento viscoso, N m s/rad."""
    fonte: str = ""
    segnaposto: bool = True


@dataclass(frozen=True)
class ParametriVincoloTerra:
    """Vincoli verso terra per asse; per il polare un valore distinto per la prova senza telescopio
    (a parita' di corrente l'asse leggero si muove molto meno: l'attrito e' piu' rigido)."""

    azimut: VincoloTerra | None = None
    polare: VincoloTerra | None = None
    polare_senza_telescopio: VincoloTerra | None = None

    def per(self, asse: str, telescopio: bool = True) -> VincoloTerra | None:
        if asse == "azimut":
            return self.azimut
        return self.polare if telescopio else self.polare_senza_telescopio


@dataclass(frozen=True)
class ParametriFEM:
    modelli_massa: tuple[str, ...]
    forma_consistente: Literal["linear", "compliance"]


@dataclass(frozen=True)
class Parametri:
    """Parametri del modello multi-GdL, in unita' SI."""

    materiali: Mapping[str, Materiale]
    azimut: ParametriAsse
    polare: ParametriAsse
    payload: ParametriPayload
    riduzione: ParametriRiduzione
    fem: ParametriFEM
    terne: ParametriTerne
    zeta: float
    theta0: tuple[float, ...]
    """Punti di lavoro in theta, rad."""
    omega0: tuple[float, ...]
    """Velocita' di inseguimento, rad/s."""
    theta0_sweep: np.ndarray
    g: float
    f_min: float
    f_max: float
    n_punti: int
    latex_lingua: str
    latex_compila: bool
    sorgente: Path | None
    avvisi: tuple[Avviso, ...] = field(default_factory=tuple)
    latex_guyan: bool = True
    """``latex.include_guyan``: nel PDF anche la riduzione di Guyan e i confronti con essa."""
    sperimentale: ParametriSperimentali | None = None
    """Sezione ``experimental`` (``None`` se assente o disabilitata)."""
    cartella_risultati: str = "results"
    """``results_dir``: cartella dei risultati, relativa alla radice del progetto."""
    vincolo_terra: ParametriVincoloTerra | None = None
    """Sezione ``ground_constraint`` (``None`` se assente o disabilitata)."""

    @property
    def mozzo(self) -> ParametriMozzo | None:
        """Modello ``hub`` attivo (``polar.hub.enabled: true``), altrimenti ``None``."""
        return self.polare.mozzo

    @property
    def omega_max(self) -> float:
        return max(self.omega0)

    @property
    def segnaposto(self) -> bool:
        """``True`` se almeno un dato e' segnaposto: le uscite vanno marcate."""
        return any(a.severita == "segnaposto" for a in self.avvisi)

    def asse(self, nome: NomeAsse) -> ParametriAsse:
        return self.azimut if nome == "azimut" else self.polare

    def avvisi_per_severita(self, severita: Severita) -> tuple[Avviso, ...]:
        return tuple(a for a in self.avvisi if a.severita == severita)

    def stampa_avvisi(self) -> None:
        stampa_avvisi(self.avvisi)


def stampa_avvisi(avvisi: Sequence[Avviso]) -> None:
    """Stampa gli avvisi in modo visibile (sezione 0: nessun risultato senza avvertenza)."""
    if not avvisi:
        return
    print("=" * 78)
    print("AVVISI SUI DATI DI INGRESSO".center(78))
    print("=" * 78)
    for a in avvisi:
        print(f"  {a}")
    print("=" * 78)


# ---------------------------------------------------------------------------
# Lettura con validazione
# ---------------------------------------------------------------------------


def _richiedi(d: Mapping[str, Any], chiave: str, contesto: str) -> Any:
    if not isinstance(d, Mapping) or chiave not in d or d[chiave] is None:
        raise ErroreConfigurazione(f"{contesto}: manca il campo obbligatorio '{chiave}'")
    return d[chiave]


def _positivo(valore: Any, nome: str, contesto: str) -> float:
    try:
        v = float(valore)
    except (TypeError, ValueError) as e:
        raise ErroreConfigurazione(f"{contesto}: '{nome}' = {valore!r} non e' un numero") from e
    if not np.isfinite(v) or v <= 0.0:
        raise ErroreConfigurazione(f"{contesto}: '{nome}' = {v} deve essere positivo")
    return v


def _opzionale(valore: Any, nome: str, contesto: str, *, scala: float = 1.0,
               positivo: bool = True) -> float | None:
    if valore is None:
        return None
    try:
        v = float(valore)
    except (TypeError, ValueError) as e:
        raise ErroreConfigurazione(f"{contesto}: '{nome}' = {valore!r} non e' un numero") from e
    if positivo and v < 0.0:
        raise ErroreConfigurazione(f"{contesto}: '{nome}' = {v} non puo' essere negativo")
    return v * scala


def _materiali(dati: Mapping[str, Any]) -> dict[str, Materiale]:
    grezzi = _richiedi(dati, "materials", "config")
    if not isinstance(grezzi, Mapping) or not grezzi:
        raise ErroreConfigurazione("materials: serve almeno un materiale")
    out: dict[str, Materiale] = {}
    for nome, m in grezzi.items():
        ctx = f"materials.{nome}"
        out[str(nome)] = Materiale(
            nome=str(nome),
            rho=_positivo(_richiedi(m, "rho", ctx), "rho", ctx),
            G=_positivo(_richiedi(m, "G", ctx), "G", ctx),
            segnaposto=bool(m.get("placeholder", False)),
        )
    return out


def leggi_strato(d: Mapping[str, Any], materiali: Mapping[str, Materiale], ctx: str) -> Strato:
    """Legge e valida uno strato (regole della sezione 5)."""
    forma = str(_richiedi(d, "shape", ctx)).strip().lower()
    d_out = _positivo(_richiedi(d, "d_out_mm", ctx), "d_out_mm", ctx) * MM
    d_in_raw = d.get("d_in_mm", None)
    if forma == "solid":
        if d_in_raw is not None and float(d_in_raw) != 0.0:
            raise ErroreConfigurazione(
                f"{ctx}: cilindro pieno (shape: solid) con d_in_mm = {d_in_raw}: "
                "per un pieno d_in_mm deve essere assente o nullo; usare shape: hollow"
            )
        d_in = 0.0
    elif forma == "hollow":
        if d_in_raw is None:
            raise ErroreConfigurazione(f"{ctx}: cilindro cavo (shape: hollow) senza d_in_mm")
        d_in = float(d_in_raw) * MM
        if not (0.0 < d_in < d_out):
            raise ErroreConfigurazione(
                f"{ctx}: cilindro cavo con d_in = {d_in / MM:g} mm e d_out = "
                f"{d_out / MM:g} mm: serve 0 < d_in < d_out"
            )
    else:
        raise ErroreConfigurazione(f"{ctx}: shape = '{forma}', ammessi 'solid' e 'hollow'")
    nome_mat = str(_richiedi(d, "material", ctx))
    if nome_mat not in materiali:
        raise ErroreConfigurazione(
            f"{ctx}: materiale '{nome_mat}' non definito in materials ({sorted(materiali)})"
        )
    return Strato(
        forma=forma,  # type: ignore[arg-type]
        d_out=d_out,
        d_in=d_in,
        materiale=materiali[nome_mat],
        portante=bool(d.get("carries_torque", True)),
    )


def _extra(lista: Any, ctx: str) -> tuple[InerziaExtra, ...]:
    out: list[InerziaExtra] = []
    for j, e in enumerate(lista or []):
        c = f"{ctx}.extra_inertias[{j}]"
        out.append(
            InerziaExtra(
                nome=str(e.get("name", f"extra {j + 1}")),
                J=float(_opzionale(_richiedi(e, "J_kgm2", c), "J_kgm2", c) or 0.0),
                m=_opzionale(e.get("m_kg"), "m_kg", c),
                I_t=_opzionale(e.get("I_t_kgm2"), "I_t_kgm2", c),
                segnaposto=bool(e.get("placeholder", False)),
            )
        )
    return tuple(out)


def leggi_elemento(
    d: Mapping[str, Any], materiali: Mapping[str, Materiale], ctx: str
) -> Elemento:
    """Legge e valida un elemento in modalita' ``geometry`` o ``manual``."""
    id_ = int(_richiedi(d, "id", ctx))
    ctx = f"{ctx} (elemento {id_})"
    modo = str(d.get("mode", "geometry")).strip().lower()
    k_joint = _opzionale(d.get("k_joint_Nm_rad"), "k_joint_Nm_rad", ctx)
    if k_joint is not None and k_joint <= 0.0:
        raise ErroreConfigurazione(f"{ctx}: k_joint_Nm_rad deve essere positivo o null")
    extra = _extra(d.get("extra_inertias"), ctx)
    k_man = _opzionale(d.get("k_Nm_rad"), "k_Nm_rad", ctx)
    if k_man is not None and k_man <= 0.0:
        raise ErroreConfigurazione(f"{ctx}: k_Nm_rad deve essere positivo")

    if modo == "geometry":
        segs_raw = _richiedi(d, "segments", ctx)
        if not segs_raw:
            raise ErroreConfigurazione(f"{ctx}: modalita' geometry senza segmenti")
        segmenti: list[Segmento] = []
        for s, seg in enumerate(segs_raw):
            cs = f"{ctx}.segments[{s}]"
            L = _positivo(_richiedi(seg, "length_mm", cs), "length_mm", cs) * MM
            strati_raw = _richiedi(seg, "layers", cs)
            if not strati_raw:
                raise ErroreConfigurazione(f"{cs}: segmento senza strati")
            strati = tuple(
                leggi_strato(st, materiali, f"{cs}.layers[{j}]")
                for j, st in enumerate(strati_raw)
            )
            if k_man is None and not any(st.portante for st in strati):
                raise ErroreConfigurazione(
                    f"{cs}: nessuno strato con carries_torque: true e nessuna k_Nm_rad "
                    "manuale: la rigidezza del segmento sarebbe nulla"
                )
            segmenti.append(Segmento(lunghezza=L, strati=strati))
        lunghezza = float(sum(s.lunghezza for s in segmenti))
        if d.get("length_mm") is not None:
            L_dich = float(d["length_mm"]) * MM
            if abs(L_dich - lunghezza) > 1e-9:
                raise ErroreConfigurazione(
                    f"{ctx}: length_mm = {L_dich / MM:g} mm diverso dalla somma dei "
                    f"segmenti ({lunghezza / MM:g} mm)"
                )
        return Elemento(
            id=id_, modo="geometry", lunghezza=lunghezza, segmenti=tuple(segmenti),
            k_joint=k_joint, extra=extra, J_manuale=None, k_manuale=k_man,
            m_manuale=None, I_t_manuale=None, y_c_manuale=None,
            segnaposto=bool(d.get("placeholder", False)),
        )

    if modo == "manual":
        J = _positivo(_richiedi(d, "J_kgm2", ctx), "J_kgm2", ctx)
        if k_man is None:
            raise ErroreConfigurazione(f"{ctx}: modalita' manual senza k_Nm_rad")
        L = _positivo(_richiedi(d, "length_mm", ctx), "length_mm", ctx) * MM
        y_c = _opzionale(d.get("y_c_mm"), "y_c_mm", ctx, scala=MM)
        if y_c is not None and y_c > L:
            raise ErroreConfigurazione(f"{ctx}: y_c_mm oltre la lunghezza dell'elemento")
        return Elemento(
            id=id_, modo="manual", lunghezza=L, segmenti=(), k_joint=k_joint,
            extra=extra, J_manuale=J, k_manuale=k_man,
            m_manuale=_opzionale(d.get("m_kg"), "m_kg", ctx),
            I_t_manuale=_opzionale(d.get("I_t_kgm2"), "I_t_kgm2", ctx),
            y_c_manuale=y_c,
            segnaposto=bool(d.get("placeholder", False)),
        )

    if modo == "telescope":
        # Nessun albero: il telescopio (rigido) collega direttamente i nodi i e i+1,
        # che coincidono con P. length_mm e' la distanza fra le due flange, solo per
        # le posizioni e il disegno.
        if d.get("segments") or d.get("J_kgm2") is not None or k_man is not None:
            raise ErroreConfigurazione(
                f"{ctx}: mode: telescope non ha albero: niente segments, J_kgm2 o k_Nm_rad "
                "(l'inerzia del telescopio va nella sezione payload)")
        L = _positivo(_richiedi(d, "length_mm", ctx), "length_mm", ctx) * MM
        return Elemento(
            id=id_, modo="telescope", lunghezza=L, segmenti=(), k_joint=None, extra=(),
            J_manuale=None, k_manuale=None, m_manuale=None, I_t_manuale=None,
            y_c_manuale=None, segnaposto=bool(d.get("placeholder", False)),
        )

    raise ErroreConfigurazione(f"{ctx}: mode = '{modo}', ammessi 'geometry', 'manual' e 'telescope'")


def _leggi_asse(
    d: Mapping[str, Any],
    nome: NomeAsse,
    materiali: Mapping[str, Materiale],
    avvisi: list[Avviso],
    asse_xy: np.ndarray | None = None,
) -> ParametriAsse:
    """``asse_xy``: posizione dell'asse di azimut nel piano (x, y) di N, m; serve alle
    inerzie nodali date come esportazione CAD (``frame: cad_cg``)."""
    chiave_cfg = "azimuth" if nome == "azimut" else "polar"
    ctx = chiave_cfg
    n_nodi = int(_richiedi(d, "n_nodes", ctx))
    if n_nodi < 2:
        raise ErroreConfigurazione(f"{ctx}: servono almeno 2 nodi")
    nodo_encoder = int(_richiedi(d, "encoder_node", ctx))
    nodo_coppia = int(_richiedi(d, "torque_node", ctx))
    nodo_estremo = int(_richiedi(d, "load_node" if nome == "azimut" else "end_node", ctx))
    for etichetta, n in (("encoder_node", nodo_encoder), ("torque_node", nodo_coppia),
                         ("load_node/end_node", nodo_estremo)):
        if not 1 <= n <= n_nodi:
            raise ErroreConfigurazione(f"{ctx}: {etichetta} = {n} fuori da 1..{n_nodi}")

    # --- elementi ---------------------------------------------------------
    elementi_raw = _richiedi(d, "elements", ctx)
    elementi = tuple(
        leggi_elemento(e, materiali, f"{ctx}.elements[{j}]") for j, e in enumerate(elementi_raw)
    )
    if len(elementi) != n_nodi - 1:
        raise ErroreConfigurazione(
            f"{ctx}: {len(elementi)} elementi ma n_nodes - 1 = {n_nodi - 1}"
        )
    ids = [e.id for e in elementi]
    if ids != list(range(1, n_nodi)):
        raise ErroreConfigurazione(f"{ctx}: gli id degli elementi devono essere 1..{n_nodi - 1} "
                                   f"contigui e ordinati, trovati {ids}")

    # --- inerzie nodali -----------------------------------------------------
    nodali: list[InerziaNodale] = []
    for j, n in enumerate(d.get("node_inertias") or []):
        c = f"{ctx}.node_inertias[{j}]"
        nodo = int(_richiedi(n, "node", c))
        if not 1 <= nodo <= n_nodi:
            raise ErroreConfigurazione(f"{c}: nodo {nodo} fuori da 1..{n_nodi}")
        alberi_polari = bool(n.get("includes_polar_shafts", False))
        if alberi_polari and not (nome == "azimut" and nodo == nodo_estremo):
            raise ErroreConfigurazione(
                f"{c}: includes_polar_shafts ha senso solo per il load di azimut (nodo {nodo_estremo})")
        if str(n.get("frame", "axis")).strip() == "cad_cg":
            if nome != "azimut":
                raise ErroreConfigurazione(f"{c}: frame: cad_cg e' previsto solo per l'asse di azimut")
            if n.get("J_kgm2") is not None:
                raise ErroreConfigurazione(f"{c}: con frame: cad_cg non si da' J_kgm2 (lo si calcola)")
            segno = str(n.get("products_sign", "tensor")).strip().lower()
            if segno not in ("tensor", "positive"):
                raise ErroreConfigurazione(f"{c}.products_sign = '{segno}': ammessi 'tensor' e 'positive'")
            m_c = _positivo(_richiedi(n, "m_kg", c), "m_kg", c)
            I_G = _tensore(_richiedi(n, "I_kgm2", c), c)
            r_G = _vettore3(_richiedi(n, "r_G_mm", c), "r_G_mm", c) * MM
            r_ax = asse_xy if asse_xy is not None else np.zeros(2)
            dist = float(np.hypot(r_G[0] - r_ax[0], r_G[1] - r_ax[1]))
            J_G = float(I_G[2, 2])  # momento diagonale: indipendente dal segno dei prodotti
            if J_G <= 0.0:
                raise ErroreConfigurazione(f"{c}: I_zz del baricentro non positivo")
            nodali.append(InerziaNodale(
                nodo=nodo, nome=str(n.get("name", f"nodo {nodo}")), J=J_G + m_c * dist**2,
                m=m_c, I_t=None, fonte=str(n.get("source", "")),
                segnaposto=bool(n.get("placeholder", False)), J_G=J_G, distanza_asse=dist,
                contiene_alberi_polari=alberi_polari))
            continue
        J_raw = n.get("J_kgm2")
        if J_raw is None:
            if nome == "azimut" and nodo == nodo_estremo:
                raise ErroreConfigurazione(
                    f"{c}: l'inerzia del load (forcella) al nodo {nodo} e' obbligatoria "
                    "(Q9): inserire il valore da un'esportazione Onshape della sola "
                    "forcella, oppure un segnaposto con placeholder: true"
                )
            raise ErroreConfigurazione(f"{c}: J_kgm2 mancante")
        nodali.append(
            InerziaNodale(
                nodo=nodo,
                nome=str(n.get("name", f"nodo {nodo}")),
                J=float(_opzionale(J_raw, "J_kgm2", c) or 0.0),
                m=_opzionale(n.get("m_kg"), "m_kg", c),
                I_t=_opzionale(n.get("I_t_kgm2"), "I_t_kgm2", c),
                fonte=str(n.get("source", "")),
                segnaposto=bool(n.get("placeholder", False)),
                contiene_alberi_polari=alberi_polari,
            )
        )
    if nome == "azimut" and not any(n.nodo == nodo_estremo for n in nodali):
        raise ErroreConfigurazione(
            f"{ctx}: manca l'inerzia del load (forcella) al nodo {nodo_estremo} (Q9)"
        )

    # --- attrito dei cuscinetti ---------------------------------------------
    attrito: list[tuple[int, float]] = []
    for k, b in (d.get("bearing_b") or {}).items():
        el = int(k)
        if not 1 <= el <= n_nodi - 1:
            raise ErroreConfigurazione(f"{ctx}.bearing_b: elemento {el} inesistente")
        b_val = float(b)
        if b_val < 0.0:
            raise ErroreConfigurazione(f"{ctx}.bearing_b: b = {b_val} negativo")
        attrito.append((el, b_val))

    # --- specifici del polare -----------------------------------------------
    el_payload: int | None = None
    xi: float | None = None
    tratto_rigido = 0.0
    el_cusc: tuple[int, ...] = ()
    y1: float | None = None
    verso = 1.0
    if nome == "polare":
        verso_raw = str(d.get("nodes_increase_along", "+J_hat")).strip()
        if verso_raw not in ("+J_hat", "-J_hat"):
            raise ErroreConfigurazione(
                f"{ctx}: nodes_increase_along = '{verso_raw}', ammessi '+J_hat' e '-J_hat'")
        verso = 1.0 if verso_raw == "+J_hat" else -1.0
        el_payload = int(_richiedi(d, "payload_element", ctx))
        if not 1 <= el_payload <= n_nodi - 1:
            raise ErroreConfigurazione(f"{ctx}: payload_element = {el_payload} inesistente")
        for e in elementi:
            if e.modo == "telescope" and e.id != el_payload:
                raise ErroreConfigurazione(
                    f"{ctx}: solo l'elemento del telescopio ({el_payload}) puo' avere "
                    f"mode: telescope, non l'elemento {e.id}")
        telescopio = elementi[el_payload - 1].modo == "telescope"
        L_p = elementi[el_payload - 1].lunghezza
        centro = d.get("payload_center_from_node3_mm")
        if telescopio:
            # P coincide con i nodi 3 e 4 (collegamento rigido): xi serve solo al disegno
            xi = 0.5 if centro is None else float(centro) * MM / L_p
            avvisi.append(Avviso(
                "TELESCOPIO_SENZA_ALBERO", "informazione",
                f"Polar element {el_payload} in mode: telescope: no shaft under the "
                "telescope. The rigid telescope directly joins nodes "
                f"{el_payload} and {el_payload + 1}, which coincide with node P; its "
                "inertia enters only through the payload section."))
        elif centro is not None:
            x_c = float(centro) * MM
            if not 0.0 <= x_c <= L_p:
                raise ErroreConfigurazione(
                    f"{ctx}: payload_center_from_node3_mm = {float(centro):g} mm fuori "
                    f"dall'elemento {el_payload} (lungo {L_p / MM:g} mm)"
                )
            xi = x_c / L_p
        else:
            xi = float(d.get("payload_xi", 0.5))
            if not 0.0 <= xi <= 1.0:
                raise ErroreConfigurazione(f"{ctx}: payload_xi = {xi} fuori da [0, 1]")
            if xi == 0.5:
                avvisi.append(Avviso(
                    "Q1_XI_DEFAULT", "segnaposto",
                    "Q1: telescope node P at the midpoint of polar element "
                    f"{el_payload} (xi_p = 0.5), i.e. telescope spread over the whole "
                    "element. Default assumption: give "
                    "payload_center_from_node3_mm if the centre of the telescope is known."))
        tratto_rigido = float(d.get("payload_clamp_rigid_mm", 0.0) or 0.0) * MM
        if tratto_rigido < 0.0:
            raise ErroreConfigurazione(f"{ctx}: payload_clamp_rigid_mm negativo")
        if telescopio:
            tratto_rigido = 0.0
        elif tratto_rigido == 0.0:
            avvisi.append(Avviso(
                "Q2_MORSETTO", "segnaposto",
                "Q2: the telescope clamp does not stiffen the shaft "
                "(payload_clamp_rigid_mm = 0, conservative assumption)."))
        el_cusc = tuple(int(e) for e in d.get("bearing_elements", []) or [])
        for e in el_cusc:
            if not 1 <= e <= n_nodi - 1:
                raise ErroreConfigurazione(f"{ctx}: bearing_elements contiene {e}, inesistente")
        if el_payload in el_cusc:
            raise ErroreConfigurazione(
                f"{ctx}: l'elemento del telescopio ({el_payload}) non puo' essere su cuscinetto"
            )
        y1 = _opzionale(d.get("y_node1_from_o2_mm"), "y_node1_from_o2_mm", ctx,
                        scala=MM, positivo=False)

    # --- nodi master ---------------------------------------------------------
    masters_raw = tuple(_richiedi(d, "masters", ctx))
    masters: list[int | str] = []
    for m in masters_raw:
        if isinstance(m, str) and m.strip().upper() == "P":
            if nome != "polare":
                raise ErroreConfigurazione(f"{ctx}: il nodo P esiste solo sul polare")
            masters.append("P")
        else:
            mi = int(m)
            if not 1 <= mi <= n_nodi:
                raise ErroreConfigurazione(f"{ctx}: master {mi} fuori da 1..{n_nodi}")
            masters.append(mi)
    if len(set(masters)) != len(masters):
        raise ErroreConfigurazione(f"{ctx}: nodi master duplicati {masters_raw}")
    obbligatori: list[int | str] = [nodo_encoder, nodo_coppia, nodo_estremo]
    if nome == "polare":
        obbligatori.append("P")
    assert all(o in masters for o in obbligatori), (
        f"{ctx}: encoder, motore, load/estremita' (e P per il polare) devono essere "
        f"nodi master (sezione 4.2): masters = {masters}, obbligatori = {obbligatori}"
    )

    return ParametriAsse(
        nome=nome, n_nodi=n_nodi, nodo_encoder=nodo_encoder, nodo_coppia=nodo_coppia,
        nodo_estremo=nodo_estremo, masters=tuple(masters), elementi=elementi,
        inerzie_nodali=tuple(nodali), attrito_cuscinetti=tuple(attrito),
        elemento_payload=el_payload, xi_payload=xi, tratto_rigido=tratto_rigido,
        elementi_cuscinetto=el_cusc, y_nodo1_da_O=y1, verso_numerazione=verso,
        mozzo=_leggi_mozzo(d.get("hub"), materiali, avvisi, ctx, elementi, el_payload)
        if nome == "polare" else None,
        attacco_telescopio=_leggi_attacco_azimut(d.get("telescope_attachment"), avvisi, ctx)
        if nome == "azimut" else None,
    )


def _leggi_attacco_azimut(d: Any, avvisi: list[Avviso], ctx_asse: str) -> AttaccoAzimut | None:
    """Sezione ``azimuth.telescope_attachment``; ``None`` se assente o con ``enabled: false``."""
    if not d or not bool(d.get("enabled", False)):
        return None
    ctx = f"{ctx_asse}.telescope_attachment"
    at = AttaccoAzimut(k=_positivo(_richiedi(d, "k_Nm_rad", ctx), "k_Nm_rad", ctx), fonte=str(d.get("source", "")),
                       segnaposto=bool(d.get("placeholder", True)))
    avvisi.append(Avviso(
        "AZIMUT_TELESCOPIO", "segnaposto" if at.segnaposto else "informazione",
        f"Azimuth telescope node (azimuth.telescope_attachment): the azimuth load (fork + hub + ring) turns with "
        f"node L and the telescope is a node P off the chain, attached to L by k = {at.k:.4g} N m/rad about z_hat"
        + (f" ({at.fonte})" if at.fonte else "") + ". The spring does not depend on theta (tilt of the tube in "
        "the ring at theta = 0)."))
    return at



def _tratto(d: Mapping[str, Any], materiali: Mapping[str, Materiale], ctx: str,
            pieno: bool = False) -> TrattoMozzo:
    mat = str(_richiedi(d, "material", ctx))
    if mat not in materiali:
        raise ErroreConfigurazione(f"{ctx}: materiale '{mat}' non definito in materials")
    d_out = _positivo(_richiedi(d, "d_out_mm", ctx), "d_out_mm", ctx) * MM
    d_in = 0.0 if pieno else float(_richiedi(d, "d_in_mm", ctx)) * MM
    if not 0.0 <= d_in < d_out:
        raise ErroreConfigurazione(f"{ctx}: serve 0 <= d_in_mm < d_out_mm")
    return TrattoMozzo(d_out=d_out, d_in=d_in,
                       lunghezza=_positivo(_richiedi(d, "length_mm", ctx), "length_mm", ctx) * MM,
                       materiale=materiali[mat])


def _leggi_mozzo(d: Any, materiali: Mapping[str, Materiale], avvisi: list[Avviso], ctx_asse: str,
                 elementi: Sequence[Elemento], el_payload: int | None) -> ParametriMozzo | None:
    """Sezione ``polar.hub``; ``None`` se assente o con ``enabled: false`` (modello originale)."""
    if not d or not bool(d.get("enabled", False)):
        return None
    ctx = f"{ctx_asse}.hub"
    lato = tuple(int(e) for e in _richiedi(d, "side_elements", ctx))
    centro = int(_richiedi(d, "center_element", ctx))
    if len(lato) != 2:
        raise ErroreConfigurazione(f"{ctx}: side_elements deve avere due elementi")
    if centro != el_payload or elementi[centro - 1].modo != "telescope":
        raise ErroreConfigurazione(f"{ctx}: center_element deve essere l'elemento del telescopio "
                                   f"(payload_element = {el_payload}) in mode: telescope")
    if sorted(lato) != [centro - 1, centro + 1]:
        raise ErroreConfigurazione(f"{ctx}: side_elements devono essere gli elementi adiacenti a {centro}")
    inz_d = _richiedi(d, "inertia", ctx)
    ci = f"{ctx}.inertia"
    if str(inz_d.get("frame", "cad_cg")) != "cad_cg":
        raise ErroreConfigurazione(f"{ci}: e' ammessa solo frame: cad_cg (esportazione Onshape)")
    I_G = np.array(_richiedi(inz_d, "I_kgm2", ci), dtype=float)
    if I_G.shape != (3, 3) or not np.allclose(I_G, I_G.T, atol=1e-9):
        raise ErroreConfigurazione(f"{ci}: I_kgm2 deve essere una matrice 3x3 simmetrica")
    r_G = np.array(_richiedi(inz_d, "r_G_mm", ci), dtype=float) * MM
    an = _richiedi(d, "ring", ctx)
    cr = f"{ctx}.ring"
    mat_an = str(_richiedi(an, "material", cr))
    if mat_an not in materiali:
        raise ErroreConfigurazione(f"{cr}: materiale '{mat_an}' non definito in materials")
    anello = AnelloMozzo(
        d_out=_positivo(_richiedi(an, "d_out_mm", cr), "d_out_mm", cr) * MM,
        d_in=_positivo(_richiedi(an, "d_in_mm", cr), "d_in_mm", cr) * MM,
        larghezza=_positivo(_richiedi(an, "width_mm", cr), "width_mm", cr) * MM,
        materiale=materiali[mat_an], poisson=float(an.get("poisson", 0.33)),
        fattore_rigidezza=_positivo(an.get("stiffness_factor", 1.0), "stiffness_factor", cr))
    if anello.d_in >= anello.d_out:
        raise ErroreConfigurazione(f"{cr}: serve d_in_mm < d_out_mm")
    lat = _richiedi(d, "side", ctx)
    att = _richiedi(d, "telescope_attachment", ctx)
    tr_d = d.get("telescope_on_ring")
    tratto = None
    if tr_d and bool(tr_d.get("enabled", False)):
        tratto = _tratto(tr_d, materiali, f"{ctx}.telescope_on_ring")
    ca = f"{ctx}.telescope_attachment"
    mz = ParametriMozzo(
        elementi_lato=(lato[0], lato[1]), elemento_centro=centro,
        m=_positivo(_richiedi(inz_d, "m_kg", ci), "m_kg", ci), r_G_N=r_G, I_G_N=I_G,
        nel_load_azimut=bool(inz_d.get("in_azimuth_load", False)),
        fonte_inerzia=str(inz_d.get("source", "")), inerzia_segnaposto=bool(inz_d.get("placeholder", False)),
        tubo=_tratto(_richiedi(lat, "tube", f"{ctx}.side"), materiali, f"{ctx}.side.tube"),
        flangia=_tratto(_richiedi(lat, "flange", f"{ctx}.side"), materiali, f"{ctx}.side.flange", pieno=True),
        anello=anello, geometria_segnaposto=bool(d.get("geometry_placeholder", False)),
        k_b=_positivo(_richiedi(att, "k_b_Nm_rad", ca), "k_b_Nm_rad", ca),
        fonte_k_b=str(att.get("source", "")), k_b_segnaposto=bool(att.get("placeholder", True)),
        tratto_telescopio=tratto, tratto_segnaposto=bool(tr_d.get("placeholder", True)) if tratto else False)
    if tratto is not None:
        avvisi.append(Avviso(
            "MOZZO_TRATTO_TELESCOPIO", "segnaposto" if mz.tratto_segnaposto else "informazione",
            f"hub.telescope_on_ring: the section of the telescope tube clamped in the support ring "
            f"({tratto.d_out / MM:g}/{tratto.d_in / MM:g} x {tratto.lunghezza / MM:g} mm, {tratto.materiale.nome}) "
            "is assumed rigid with the ring: its inertia is moved from the telescope to the hub centre H."))
    avvisi[:] = [Avviso(a.codice, a.severita,
                        f"Polar element {centro} in mode: telescope: no shaft under the telescope. In the hub "
                        f"model its nodes {centro} and {centro + 1} coincide with the hub centre H; the telescope "
                        "is the node P attached to H by k_b, and its inertia enters only through the payload "
                        "section.") if a.codice == "TELESCOPIO_SENZA_ALBERO" else a for a in avvisi]
    avvisi.append(Avviso(
        "MOZZO_TELESCOPIO", "informazione",
        f"Hub model (polar.hub.enabled: true): the hub (two side hubs, tube + flange, and the support ring) "
        f"replaces polar elements {list(lato)} and its centre H takes the nodes of element {centro}; the "
        f"telescope is a node P off the chain, attached to H by k_b = {mz.k_b:.4g} N m/rad."))
    if anello.fattore_rigidezza != 1.0:
        avvisi.append(Avviso(
            "MOZZO_ANELLO_IRRIGIDITO", "segnaposto",
            f"hub.ring.stiffness_factor = {anello.fattore_rigidezza:g}: the stiffness of the support ring from the "
            "half-ring model (ring alone, lower bound) is multiplied by this factor, chosen from the test without "
            "the telescope (no structural mode below ~400 Hz)."))
    if mz.k_b_segnaposto:
        avvisi.append(Avviso(
            "MOZZO_K_B", "segnaposto",
            f"k_b = {mz.k_b:.4g} N m/rad (telescope to hub) is a placeholder"
            + (f" ({mz.fonte_k_b})" if mz.fonte_k_b else "") + ": it sets the antiresonance at the encoder."))
    if not mz.nel_load_azimut:
        avvisi.append(Avviso(
            "MOZZO_INERZIA_AZIMUT", "segnaposto",
            "hub.inertia.in_azimuth_load = false: the inertia of the hub about z_hat is added to J^ax. "
            "Check whether the azimuth load exported from Onshape already contains the hub parts."))
    return mz

def _avvisi_asse(asse: ParametriAsse, avvisi: list[Avviso]) -> None:
    """Avvisi sui segnaposto e sui dati incompleti di un asse."""
    seg = [e.id for e in asse.elementi if e.segnaposto]
    if seg:
        avvisi.append(Avviso(
            f"GEOMETRIA_{asse.nome.upper()}", "segnaposto",
            f"{in_inglese(asse.nome).capitalize()} axis: elements {seg} marked placeholder. Shaft "
            "inertias and stiffnesses, and therefore every frequency and transfer function, are a "
            f"{NOTA_REALISMO}."))
    nodi_con_inerzia = {n.nodo for n in asse.inerzie_nodali}
    senza = [f"{ruolo} (node {nodo})" for ruolo, nodo in (("encoder", asse.nodo_encoder),
                                                          ("motor", asse.nodo_coppia))
             if nodo not in nodi_con_inerzia]
    if senza:
        avvisi.append(Avviso(
            f"NODI_SENZA_INERZIA_{asse.nome.upper()}", "segnaposto",
            f"{in_inglese(asse.nome).capitalize()} axis: no concentrated inertia at {', '.join(senza)} "
            f"({'node_inertias entirely commented out' if not asse.inerzie_nodali else 'rotor/disc missing'}): "
            "those nodes carry only the shaft mass and are almost inertia-free. This affects "
            "the zeros and the high modes (typically strongly) and the choice of the truncation modes."))
    for n in asse.inerzie_nodali:
        if n.segnaposto:
            avvisi.append(Avviso(
                f"INERZIA_NODALE_{asse.nome.upper()}_{n.nodo}", "segnaposto",
                f"{in_inglese(asse.nome).capitalize()} axis, node {n.nodo} ({n.nome}): J = {n.J:.4g} kg m^2 "
                f"placeholder ({n.fonte or 'source not given'})."))
    if all(e.k_joint is None for e in asse.elementi):
        avvisi.append(Avviso(
            f"Q5_GIUNZIONI_{asse.nome.upper()}", "segnaposto",
            f"Q5: {in_inglese(asse.nome)} axis without joint compliances (shrink fits, "
            "keys, flanges): the stiffnesses are those of the shafts alone and "
            "probably overestimate the real one."))
    if asse.nome == "polare":
        for e in asse.elementi:
            if e.modo == "manual" and None in (e.m_manuale, e.I_t_manuale, e.y_c_manuale):
                avvisi.append(Avviso(
                    f"MANUALE_POLARE_{e.id}", "incoerenza",
                    f"Polar element {e.id} in manual mode without m_kg, I_t_kgm2 "
                    "or y_c_mm: its contribution to the azimuth inertia is set to zero."))
        for n in asse.inerzie_nodali:
            if n.m is None or n.I_t is None:
                avvisi.append(Avviso(
                    f"NODALE_POLARE_{n.nodo}", "incoerenza",
                    f"Polar node inertia at node {n.nodo} ({n.nome}) without m_kg or "
                    "I_t_kgm2: its contribution to the azimuth inertia is set to zero."))
        if asse.y_nodo1_da_O is None:
            avvisi.append(Avviso(
                "Q7_POSIZIONE_POLARE", "segnaposto",
                "Q7: y_node1_from_o2_mm = null. The contribution of the polar shaft "
                "to the azimuth inertia (J^ax) is computed from the transverse moments "
                "of the parts about their own centres only, without Steiner terms about o_2."))
    b = asse.b_cuscinetti
    if not b or all(v == 0.0 for v in b.values()):
        avvisi.append(Avviso(
            f"Q8_CUSCINETTI_{asse.nome.upper()}", "segnaposto",
            f"Q8: zero viscous bearing friction on the {in_inglese(asse.nome)} axis."))


def _vettore3(v: Any, nome: str, ctx: str) -> np.ndarray:
    try:
        a = np.array([float(x) for x in v], dtype=float)
    except (TypeError, ValueError) as e:
        raise ErroreConfigurazione(f"{ctx}: '{nome}' deve essere una lista di 3 numeri") from e
    if a.shape != (3,):
        raise ErroreConfigurazione(f"{ctx}: '{nome}' deve avere 3 componenti, ne ha {a.size}")
    return a


def _tensore(v: Any, ctx: str) -> np.ndarray:
    """Tensore 3x3 simmetrico, come scritto nel file (kg m^2)."""
    try:
        I = np.array([[float(x) for x in riga] for riga in v], dtype=float)
    except (TypeError, ValueError) as e:
        raise ErroreConfigurazione(f"{ctx}: 'I_kgm2' deve essere una matrice 3x3 di numeri") from e
    if I.shape != (3, 3):
        raise ErroreConfigurazione(f"{ctx}: 'I_kgm2' deve essere 3x3, trovato {I.shape}")
    scala = max(float(np.max(np.abs(I))), 1e-300)
    if np.max(np.abs(I - I.T)) > 1e-9 * scala:
        raise ErroreConfigurazione(
            f"{ctx}: 'I_kgm2' non e' simmetrico (I_ij deve coincidere con I_ji)")
    return 0.5 * (I + I.T)


def _leggi_payload(pl: Mapping[str, Any], base: Path, avvisi: list[Avviso]) -> ParametriPayload:
    """Sezione ``payload`` di ``config.yaml``."""
    ctx = "payload"
    sorgente = str(pl.get("source", "manual")).strip().lower()
    if sorgente not in ("manual", "gimbal_2dof"):
        raise ErroreConfigurazione(f"{ctx}.source = '{sorgente}': ammessi 'manual' e 'gimbal_2dof'")
    cfg2 = Path(str(pl.get("config_2dof", "../gimbal_2dof/config.yaml")))
    if not cfg2.is_absolute():
        cfg2 = (base / cfg2).resolve()
    abilitato = bool(pl.get("enabled", True))
    if sorgente == "gimbal_2dof":
        return ParametriPayload(abilitato=abilitato, sorgente="gimbal_2dof", config_2dof=cfg2,
                                sottrai_alberi=bool(pl.get("subtract_shaft_parts", True)))

    terna = str(_richiedi(pl, "frame", ctx)).strip()
    if terna not in ("cad_cg", "body2_o2"):
        raise ErroreConfigurazione(f"{ctx}.frame = '{terna}': ammessi 'cad_cg' e 'body2_o2'")
    segno = str(pl.get("products_sign", "tensor")).strip().lower()
    if segno not in ("tensor", "positive"):
        raise ErroreConfigurazione(f"{ctx}.products_sign = '{segno}': ammessi 'tensor' e 'positive'")
    m = _positivo(_richiedi(pl, "m_kg", ctx), "m_kg", ctx)
    I = _tensore(_richiedi(pl, "I_kgm2", ctx), ctx)
    r_G_N = d = None
    if terna == "cad_cg":
        r_G_N = _vettore3(_richiedi(pl, "r_G_mm", ctx), "r_G_mm", ctx) * MM
    else:
        d = _vettore3(_richiedi(pl, "d_mm", ctx), "d_mm", ctx) * MM
    segnaposto = bool(pl.get("placeholder", False))
    nota = str(pl.get("source_note", ""))
    if segnaposto:
        avvisi.append(Avviso(
            "PAYLOAD_SEGNAPOSTO", "segnaposto",
            f"Telescope tensor, mass and centre of mass in config.yaml marked placeholder"
            + (f" ({nota})" if nota else "") + "."))
    return ParametriPayload(
        abilitato=abilitato, sorgente="manual", config_2dof=cfg2,
        sottrai_alberi=bool(pl.get("subtract_shaft_parts", False)), terna=terna,  # type: ignore[arg-type]
        m=m, I=I, r_G_N=r_G_N, d=d, segno_prodotti=segno, nota=nota,  # type: ignore[arg-type]
        segnaposto=segnaposto,
    )


def _leggi_sperimentale(d: Mapping[str, Any] | None, base: Path,
                        avvisi: list[Avviso]) -> ParametriSperimentali | None:
    """Sezione ``experimental``: misure da confrontare con il modello (cartella ``Validation`` esclusa)."""
    if not d or not bool(d.get("enabled", True)):
        return None
    ctx = "experimental"
    cartella = Path(str(_richiedi(d, "data_dir", ctx)))
    if not cartella.is_absolute():
        cartella = (base / cartella).resolve()
    Kt = _positivo(_richiedi(d, "torque_constant_Nm_per_A", ctx), "torque_constant_Nm_per_A", ctx)
    banda = [float(v) for v in d.get("level_band_Hz", [15.0, 300.0])]
    if len(banda) != 2 or not 0 < banda[0] < banda[1]:
        raise ErroreConfigurazione(f"{ctx}.level_band_Hz deve essere [f_a, f_b] con 0 < f_a < f_b")
    voci = []
    for k, m in enumerate(_richiedi(d, "measurements", ctx)):
        c = f"{ctx}.measurements[{k}]"
        file = str(_richiedi(m, "file", c))
        percorso = (cartella / file).resolve()
        if "Validation" in percorso.relative_to(cartella).parts:
            raise ErroreConfigurazione(f"{c}: i dati della cartella Validation non vanno usati")
        if not percorso.exists():
            raise ErroreConfigurazione(f"{c}: file {percorso} inesistente")
        asse = str(_richiedi(m, "axis", c))
        if asse not in ("azimut", "polare"):
            raise ErroreConfigurazione(f"{c}.axis = '{asse}': ammessi azimut, polare")
        th_seg = bool(m.get("theta0_placeholder", False))
        voci.append(VoceMisura(file=file, asse=asse, theta0=float(np.deg2rad(float(m.get("theta0_deg", 0.0)))),  # type: ignore[arg-type]
                               telescopio=bool(m.get("telescope", True)), etichetta=str(m.get("label", file)),
                               theta0_segnaposto=th_seg))
        if th_seg:
            avvisi.append(Avviso(
                f"MISURA_THETA0_{k}", "segnaposto",
                f"{file}: elevation during the test not recorded in the file; theta0 = "
                f"{float(m.get('theta0_deg', 0.0)):g} deg assumed."))
    avvisi.append(Avviso(
        "MISURE_SPERIMENTALI", "informazione",
        f"Comparison with {len(voci)} plant measurements (Kollmorgen WorkBench, {cartella.name}); K_t = {Kt:g} N m/A. "
        "The Validation folder is not used."))
    return ParametriSperimentali(cartella=cartella, Kt=Kt, traccia=str(d.get("trace", "Plant")),
                                 f_min=_positivo(d.get("f_min_Hz", 0.1), "f_min_Hz", ctx),
                                 f_max=_positivo(d.get("f_max_Hz", 2000.0), "f_max_Hz", ctx),
                                 banda_livello=(banda[0], banda[1]), misure=tuple(voci))


def _leggi_vincolo_terra(d: Mapping[str, Any] | None, avvisi: list[Avviso]) -> ParametriVincoloTerra | None:
    """Sezione ``ground_constraint``: ``azimuth``, ``polar``, ``polar_without_telescope``, ciascuna con
    ``k_Nm_rad`` e facoltativi ``eta``, ``c_Nms_rad``, ``source``, ``placeholder``."""
    if not d or not d.get("enabled", True):
        return None
    out: dict[str, VincoloTerra | None] = {}
    for chiave, nome in (("azimuth", "azimut"), ("polar", "polare"),
                         ("polar_without_telescope", "polare_senza_telescopio")):
        v = d.get(chiave)
        if not v:
            out[nome] = None
            continue
        ctx = f"ground_constraint.{chiave}"
        k = _positivo(_richiedi(v, "k_Nm_rad", ctx), "k_Nm_rad", ctx)
        eta, c = float(v.get("eta", 0.0)), float(v.get("c_Nms_rad", 0.0))
        if eta < 0 or c < 0:
            raise ErroreConfigurazione(f"{ctx}: eta e c_Nms_rad devono essere non negativi")
        out[nome] = VincoloTerra(k=k, eta=eta, c=c, fonte=str(v.get("source", "")),
                                 segnaposto=bool(v.get("placeholder", True)))
    avvisi.append(Avviso(
        "VINCOLO_TERRA", "segnaposto",
        "ground_constraint: linear equivalent of the bearing friction below a few hertz, "
        + ", ".join(f"{in_inglese(n) if n != 'polare_senza_telescopio' else 'polar without telescope'} "
                    f"k = {v.k:.4g} N m/rad, eta = {v.eta:g}" for n, v in out.items() if v)
        + ". It holds for the amplitude of the test it was fitted on, it acts in the frequency "
          "responses only (hysteretic loss) and it is not in the time-domain simulator."))
    return ParametriVincoloTerra(**out)


def _leggi_modale(d: Mapping[str, Any]) -> ParametriModali:
    """Sezione ``reduction.modal`` (troncamento modale alla Gawronski)."""
    ctx = "reduction.modal"
    f_banda = _positivo(d.get("f_band_Hz", 2000.0), "f_band_Hz", ctx)
    fattore = float(d.get("cutoff_factor", 2.0))
    if fattore < 1.0:
        raise ErroreConfigurazione(f"{ctx}.cutoff_factor = {fattore} < 1: il taglio sarebbe dentro la banda")
    h = d.get("hankel_min_ratio")
    hankel = None if h is None else float(h)
    if hankel is not None and not 0.0 < hankel <= 1.0:
        raise ErroreConfigurazione(f"{ctx}.hankel_min_ratio = {hankel} fuori da (0, 1]")
    nf_raw = d.get("n_flex", {}) or {}
    nf: dict[str, int | None] = {}
    for chiave, nome in (("azimut", "azimut"), ("polare", "polare")):
        v = nf_raw.get(chiave, nf_raw.get({"azimut": "azimuth", "polare": "polar"}[chiave]))
        if v is not None and int(v) < 0:
            raise ErroreConfigurazione(f"{ctx}.n_flex.{chiave} negativo")
        nf[nome] = None if v is None else int(v)
    corr = str(d.get("static_correction", "none")).strip()
    if corr not in ("none", "residual_vectors"):
        raise ErroreConfigurazione(f"{ctx}.static_correction = '{corr}': ammessi none, residual_vectors")
    return ParametriModali(f_banda=f_banda, fattore_taglio=fattore, hankel_min=hankel,
                           n_flessibili=nf, theta_ref=float(np.deg2rad(float(d.get("theta_ref_deg", 0.0)))),
                           correzione_statica=corr)  # type: ignore[arg-type]


def _leggi_terne(f: Mapping[str, Any], avvisi: list[Avviso]) -> ParametriTerne:
    """Sezione ``frames``: quote delle origini e posa CAD (convenzione De Donatis)."""
    from gimbal.params import costruisci_terne, direzione

    ctx = "frames"
    xy = np.array([float(v) for v in f.get("azimuth_axis_xy_mm", [0.0, 0.0])]) * MM
    if xy.shape != (2,):
        raise ErroreConfigurazione(f"{ctx}.azimuth_axis_xy_mm deve avere 2 componenti")
    z_o1 = float(_richiedi(f, "o1_z_mm", ctx)) * MM
    z_o2 = float(_richiedi(f, "o2_z_mm", ctx)) * MM
    try:
        I0 = direzione(f.get("pointing_cad", "-Y"))
        J0 = direzione(f.get("polar_axis_cad", "+X"))
    except ValueError as e:
        raise ErroreConfigurazione(f"{ctx}: {e}") from e
    theta_cad = float(np.deg2rad(float(f.get("theta_cad_deg", 0.0))))
    try:
        # costruisci_terne (2 GdL) vuole j_hat = z_hat x i_hat e theta che abbassa il
        # puntamento: J_hat = -j_hat e theta_P = -theta.
        j0_calc, R1, R2 = costruisci_terne(I0, -J0, -theta_cad)
    except (AssertionError, ValueError) as e:
        raise ErroreConfigurazione(
            f"{ctx}: pointing_cad e polar_axis_cad incoerenti (serve J_hat = K_hat x I_hat = "
            f"i_hat x z_hat, terna destra): {e}") from e
    J0_calc = -j0_calc
    avvisi.append(Avviso(
        "TERNE_ORIGINI", "informazione",
        f"Origins on z_hat: dist(O_N, o_1) = {z_o1 / MM:g} mm, dist(O_N, o_2) = {z_o2 / MM:g} mm "
        "(o_2 = centre of rotation of the polar shaft). O_N is the origin of the Onshape assembly "
        "(centre of the base of the azimuth shaft) and x_hat, y_hat, z_hat are its axes."))
    return ParametriTerne(asse_az_xy=xy, z_o1=z_o1, z_o2=z_o2, i_hat_0=I0, J_hat_0=J0_calc,
                          theta_cad=theta_cad, R_N1_0=R1, R_N2_0=R2 @ R_P2, R_NP_0=R2)


def carica_config(percorso: str | Path) -> dict[str, Any]:
    percorso = Path(percorso)
    with percorso.open("r", encoding="utf-8") as f:
        dati = yaml.safe_load(f)
    if not isinstance(dati, dict):
        raise ErroreConfigurazione(f"{percorso}: il file non contiene una mappa YAML")
    return dati


def _griglia(spec: Any, nome: str) -> np.ndarray:
    if isinstance(spec, Mapping):
        start, stop, step = (float(spec[k]) for k in ("start", "stop", "step"))
        if step == 0.0:
            raise ErroreConfigurazione(f"{nome}: step nullo")
        n = int(round(abs(stop - start) / abs(step))) + 1
        return np.linspace(start, stop, n)
    if isinstance(spec, Sequence):
        return np.array([float(v) for v in spec])
    raise ErroreConfigurazione(f"{nome}: formato non riconosciuto")


def parametri_da_config(
    percorso: str | Path | None = None, *, dati: dict[str, Any] | None = None
) -> Parametri:
    """Costruisce :class:`Parametri` da ``config.yaml`` o da un dizionario.

    Parameters
    ----------
    percorso
        File YAML. Se ``None`` (e ``dati`` e' ``None``) si usa
        ``config.yaml`` nella radice del progetto. Con ``dati`` serve solo a
        risolvere i percorsi relativi (``payload.config_2dof``).
    dati
        Dizionario gia' caricato, per costruire configurazioni nei test.
    """
    if dati is None:
        percorso = Path(percorso) if percorso is not None else RADICE / "config.yaml"
        dati = carica_config(percorso)
    base = Path(percorso).resolve().parent if percorso is not None else RADICE

    avvisi: list[Avviso] = []
    materiali = _materiali(dati)
    mat_seg = [m.nome for m in materiali.values() if m.segnaposto]
    if mat_seg:
        avvisi.append(Avviso(
            "Q6_MATERIALI", "segnaposto",
            f"Q6: placeholder materials {[in_inglese(n) for n in mat_seg]} (" + ", ".join(
                f"{in_inglese(m.nome)}: rho = {m.rho:g} kg/m^3, G = {m.G / 1e9:g} GPa"
                for m in materiali.values() if m.segnaposto) + ")."))

    asse_xy = np.array([float(v) for v in (dati.get("frames", {}) or {}).get(
        "azimuth_axis_xy_mm", [0.0, 0.0])]) * MM
    azimut = _leggi_asse(_richiedi(dati, "azimuth", "config"), "azimut", materiali, avvisi, asse_xy)
    polare = _leggi_asse(_richiedi(dati, "polar", "config"), "polare", materiali, avvisi)
    _avvisi_asse(azimut, avvisi)
    _avvisi_asse(polare, avvisi)
    if azimut.attacco_telescopio is not None and polare.mozzo is None:
        raise ErroreConfigurazione("azimuth.telescope_attachment richiede il modello hub (polar.hub.enabled: true): "
                                   "il telescopio e' appeso al mozzo su entrambi gli assi")
    if any(n.contiene_alberi_polari for n in azimut.inerzie_nodali):
        # J^ax non si aggiunge (e' gia' nel load di azimut): la posizione del nodo 1 polare non conta
        avvisi = [Avviso(a.codice, "informazione",
                         a.messaggio + " Irrelevant: the azimuth load includes the polar shafts "
                         "(includes_polar_shafts: true), so J^ax is not added.")
                  if a.codice == "Q7_POSIZIONE_POLARE" else a for a in avvisi]

    # --- payload ---------------------------------------------------------------
    terne = _leggi_terne(_richiedi(dati, "frames", "config"), avvisi)
    if polare.mozzo is not None:
        from .mozzo import inerzia as _inerzia_mozzo

        from .mozzo import inerzia_tratto as _inerzia_tratto

        im = _inerzia_mozzo(polare.mozzo, terne)
        it = _inerzia_tratto(polare.mozzo)
        dJ = it[2] if it is not None else 0.0
        polare = replace(polare, mozzo=replace(polare.mozzo, J_polare=im.J_polare + dJ,
                                                J_azimut=im.J_azimut + dJ))
    payload = _leggi_payload(dati.get("payload", {}) or {}, base, avvisi)

    # --- riduzione e FEM --------------------------------------------------------
    rd = dati.get("reduction", {}) or {}
    metodo = str(rd.get("method", "modal"))
    forma_massa = str(rd.get("mass_form", "guyan"))
    if metodo not in ("guyan", "craig_bampton", "modal"):
        raise ErroreConfigurazione(
            f"reduction.method = '{metodo}' non ammesso (guyan | craig_bampton | modal)")
    if forma_massa not in ("guyan", "guyan_lumped"):
        raise ErroreConfigurazione(f"reduction.mass_form = '{forma_massa}' non ammesso")
    n_modi = int(rd.get("n_fixed_interface_modes", 0))
    if n_modi < 0:
        raise ErroreConfigurazione("reduction.n_fixed_interface_modes negativo")
    riduzione = ParametriRiduzione(metodo=metodo, n_modi_interfaccia=n_modi,  # type: ignore[arg-type]
                                   forma_massa=forma_massa,  # type: ignore[arg-type]
                                   modale=_leggi_modale(rd.get("modal", {}) or {}))
    if metodo == "modal":
        avvisi.append(Avviso(
            "RIDUZIONE_MODALE", "informazione",
            "Reduction: Gawronski modal truncation of the full FEM of each axis "
            f"(band {riduzione.modale.f_banda:g} Hz, cut-off {riduzione.modale.f_taglio:g} Hz, "
            f"static correction: {riduzione.modale.correzione_statica})."))
    else:
        avvisi.append(Avviso(
            "Q4_RIDUZIONE", "informazione",
            f"Reduction: {metodo} on the master nodes (user decision 2026-10-12, Q4); "
            "modal truncation (method: modal) is the default method since 2026-10-13."))

    fm = dati.get("fem", {}) or {}
    modelli = tuple(str(m) for m in fm.get("mass_models", ["lumped", "consistent"]))
    for m in modelli:
        if m not in ("lumped", "consistent"):
            raise ErroreConfigurazione(f"fem.mass_models: '{m}' non ammesso")
    forma_c = str(fm.get("consistent_shape", "linear"))
    if forma_c not in ("linear", "compliance"):
        raise ErroreConfigurazione(f"fem.consistent_shape = '{forma_c}' non ammesso")
    fem = ParametriFEM(modelli_massa=modelli, forma_consistente=forma_c)  # type: ignore[arg-type]

    # --- smorzamento ----------------------------------------------------------
    dm = dati.get("damping", {}) or {}
    if str(dm.get("model", "stiffness_proportional")) != "stiffness_proportional":
        raise ErroreConfigurazione("damping.model: ammesso solo 'stiffness_proportional'")
    zeta = float(dm.get("zeta_first_flexible", 0.01))
    if not 0.0 <= zeta < 1.0:
        raise ErroreConfigurazione(f"damping.zeta_first_flexible = {zeta} fuori da [0, 1)")
    avvisi.append(Avviso(
        "Q8_SMORZAMENTO", "segnaposto",
        f"Q8: stiffness-proportional damping, C = beta K, with zeta = {zeta:g} "
        "on the first flexible mode of each axis. Placeholder: the height of the FRF peaks and "
        "notches is not representative."))

    # --- punti di lavoro --------------------------------------------------------
    op = dati.get("operating_points", {}) or {}
    theta0 = tuple(float(np.deg2rad(v)) for v in op.get("theta0_deg", [0.0]))
    omega0_deg = tuple(float(v) for v in op.get("omega0_deg_s", [0.0, 10.0]))
    omega0 = tuple(float(np.deg2rad(v)) for v in omega0_deg)
    avvisi.append(Avviso(
        "Q11_OMEGA_MAX", "segnaposto",
        f"Q11: tracking rates {list(omega0_deg)} deg/s; Omega_max = "
        f"{max(omega0_deg):g} deg/s placeholder. It sets the gyroscopic coupling."))
    sweep = np.deg2rad(_griglia(op.get("theta0_sweep_deg",
                                       {"start": -45.0, "stop": 45.0, "step": 5.0}),
                                "theta0_sweep_deg"))
    avvisi.append(Avviso(
        "Q12_FORCELLA_RIGIDA", "informazione",
        "Q12: rigid fork between azimuth node 9 and the polar bearings. The "
        "reaction of tau_theta does no work on psi_L (section 9.4)."))

    fr = dati.get("frequency", {}) or {}
    lt = dati.get("latex", {}) or {}
    vincolo_terra = _leggi_vincolo_terra(dati.get("ground_constraint"), avvisi)

    return Parametri(
        materiali=materiali,
        azimut=azimut,
        polare=polare,
        payload=payload,
        riduzione=riduzione,
        fem=fem,
        terne=terne,
        zeta=zeta,
        theta0=theta0,
        omega0=omega0,
        theta0_sweep=sweep,
        g=float(dati.get("g", 9.80665)),
        f_min=_positivo(fr.get("f_min_Hz", 0.1), "f_min_Hz", "frequency"),
        f_max=_positivo(fr.get("f_max_Hz", 5000.0), "f_max_Hz", "frequency"),
        n_punti=int(fr.get("n_points", 4000)),
        latex_lingua=str(lt.get("language", "en")),
        latex_compila=bool(lt.get("compile", True)),
        latex_guyan=bool(lt.get("include_guyan", True)),
        sperimentale=_leggi_sperimentale(dati.get("experimental"), base, avvisi),
        sorgente=Path(percorso) if percorso is not None else None,
        avvisi=tuple(avvisi),
        cartella_risultati=str(dati.get("results_dir", "results")),
        vincolo_terra=vincolo_terra,
    )
