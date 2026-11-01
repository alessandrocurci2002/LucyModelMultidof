"""Confronto del modello a 7 GdL (parametri concentrati) con le risposte in frequenza misurate.

Misure: esportazioni CSV di Kollmorgen WorkBench ("Manual Tuning Measurement", tipo
``Plant``), cartella ``sim/LucyIdentification`` (la sottocartella ``Validation`` non si usa,
per richiesta dell'utente). Il drive inietta una corrente sinusoidale (seno a passi,
``CurrentExcitation`` A) con gli anelli chiusi e stima la pianta come rapporto fra la
velocita' di retroazione e la corrente, in (rad/s)/A: la verifica e' che il moto rigido del
modello, :math:`K_t/(J\\omega)`, coincide con la misura di elevazione a 10 Hz entro 0.2 dB.

Il modello e' convertito nelle stesse unita':

.. math:: G_{meas}(j\\omega) = \\frac{\\dot\\theta(j\\omega)}{I(j\\omega)} = j\\omega\\,K_t\\,
          \\frac{\\theta(j\\omega)}{\\tau(j\\omega)},

con :math:`\\theta/\\tau` la FdT coppia -> posizione del modello linearizzato (Omega_0 = 0) al
nodo encoder o al nodo motore. Il lettore dei CSV riprende ``load_bode`` di
``sim/LucyIdentification/plot_bode.py``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from . import concentrato as cc
from . import coupled as cp
from . import fem
from . import linearize as lz
from .params import Parametri, ParametriSperimentali
from .payload import Payload


# ---------------------------------------------------------------------------
# Lettura dei CSV di Kollmorgen WorkBench
# ---------------------------------------------------------------------------


def carica_bode(percorso: str | Path) -> tuple[np.ndarray, dict[str, tuple[np.ndarray, np.ndarray]],
                                               dict[str, str]]:
    """``(f [Hz], {traccia: (modulo [dB], fase [deg])}, metadati)`` da un'esportazione WorkBench.

    Le prime colonne sono la misura; le colonne ``#ParentName, #Name, #Type, #Value`` in coda
    sono metadati (``"Padre/Nome" -> valore``).
    """
    with open(percorso, newline="", encoding="utf-8-sig") as f:
        righe = list(csv.reader(f))
    intest = righe[0]
    dati = [r for r in righe[1:] if r and r[0].strip()]
    n_col = next(i for i, h in enumerate(intest) if not h.strip())
    tab = np.array([[float(x) for x in r[:n_col]] for r in dati])
    tracce: dict[str, list] = {}
    for i, h in enumerate(intest[:n_col]):
        if h.startswith("Frequency"):
            continue
        nome, _, tipo = h.partition("[")[0].strip().rpartition(" ")
        if tipo == "Gain":
            tracce.setdefault(nome, [None, None])[0] = tab[:, i]
        elif tipo == "Phase":
            tracce.setdefault(nome, [None, None])[1] = tab[:, i]
    meta: dict[str, str] = {}
    if "#ParentName" in intest:
        j = intest.index("#ParentName")
        for r in righe[1:]:
            if len(r) > j + 3 and r[j + 1]:
                meta.setdefault(f"{r[j]}/{r[j + 1]}", r[j + 3])
    return tab[:, 0], {k: (v[0], v[1]) for k, v in tracce.items()}, meta


@dataclass(frozen=True)
class Misura:
    """Una misura di pianta: dati, condizioni e parte del controllore attivo durante la prova."""

    file: str
    asse: str
    theta0: float
    telescopio: bool
    etichetta: str
    f: np.ndarray
    modulo_db: np.ndarray
    fase_deg: np.ndarray
    eccitazione_A: float
    f_fine: float
    filtri: tuple[tuple[float, float, float, float], ...]
    """Filtri anti-risonanza in avanti (f_num, Q_num, f_den, Q_den) della configurazione di misura."""
    meta: Mapping[str, str] = field(default_factory=dict)

    @property
    def G(self) -> np.ndarray:
        return 10 ** (self.modulo_db / 20) * np.exp(1j * np.deg2rad(self.fase_deg))


def carica_misura(cfg: ParametriSperimentali, voce) -> Misura:
    f, tracce, meta = carica_bode(cfg.cartella / voce.file)
    if cfg.traccia not in tracce:
        raise ValueError(f"{voce.file}: traccia '{cfg.traccia}' assente (disponibili: {list(tracce)})")
    g, p = tracce[cfg.traccia]
    filtri = []
    for k in range(8):
        base = f"AntiResonanceFilters[{k}]/"
        if meta.get(base + "ArfGroup") == "ArfGroupForward":
            filtri.append(tuple(float(meta[base + c]) for c in ("NumeratorF", "NumeratorQ", "DenominatorF",
                                                                    "DenominatorQ")))
    return Misura(file=voce.file, asse=voce.asse, theta0=voce.theta0, telescopio=voce.telescopio,
                  etichetta=voce.etichetta, f=f, modulo_db=g, fase_deg=p,
                  eccitazione_A=float(meta.get("MeasurementConfiguration/CurrentExcitation", "nan")),
                  f_fine=float(meta.get("MeasurementConfiguration/SineExcitationEndFrequency", f[-1])),
                  filtri=tuple(filtri), meta=meta)


# ---------------------------------------------------------------------------
# FdT del modello nelle unita' della misura
# ---------------------------------------------------------------------------


def payload_senza_telescopio(pl: Payload) -> Payload:
    """Payload disabilitato (prova senza telescopio): restano gli alberi polari e J^ax."""
    return replace(pl, abilitato=False, I=np.zeros((3, 3)), m=0.0, d=np.zeros(3))


def canali(asse: str) -> tuple[str, str]:
    return ("psi_e", "psi_m") if asse == "azimut" else ("theta_e", "theta_m")


def frf_velocita(ml: lz.ModelloLineare, asse: str, f: np.ndarray, Kt: float) -> dict[str, np.ndarray]:
    """:math:`j\\omega K_t\\,\\theta/\\tau` ai nodi encoder e motore dell'asse eccitato, (rad/s)/A."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    j = 0 if asse == "azimut" else 1
    u = canali(asse)
    G = ml.frf(w, list(u))[:, :, j]
    return {c: 1j * w * Kt * G[:, i] for i, c in enumerate(u)}


@dataclass(frozen=True)
class Caratteristica:
    """Picco o valle del modulo misurato e corrispondente polo/zero del modello."""

    tipo: str
    """``risonanza`` o ``antirisonanza``."""
    f_misura: float
    modulo_db: float
    f_modello: float | None
    canale: str | None


def estremi(f: np.ndarray, g_db: np.ndarray, *, prominenza_db: float = 8.0,
            finestra_dec: float = 0.25) -> list[tuple[str, float, float]]:
    """Picchi e valli del modulo con prominenza locale almeno ``prominenza_db`` (dB).

    La prominenza si misura in una finestra di +/- ``finestra_dec`` decadi: altezza del
    picco (o profondita' della valle) rispetto al piu' alto dei due minimi (massimi) ai lati.
    Esclude gli estremi ai bordi della misura.
    """
    out = []
    lf = np.log10(f)
    for i in range(1, len(g_db) - 1):
        for tipo, segno in (("risonanza", 1.0), ("antirisonanza", -1.0)):
            x = segno * g_db
            if not (x[i] > x[i - 1] and x[i] >= x[i + 1]):
                continue
            sx = x[(lf >= lf[i] - finestra_dec) & (lf < lf[i])]
            dx = x[(lf > lf[i]) & (lf <= lf[i] + finestra_dec)]
            if sx.size == 0 or dx.size == 0 or lf[i] - finestra_dec < lf[0] or lf[i] + finestra_dec > lf[-1]:
                continue
            if x[i] - max(sx.min(), dx.min()) >= prominenza_db:
                out.append((tipo, float(f[i]), float(g_db[i])))
    return out


def poli_zeri_modello(ml: lz.ModelloLineare, uscita: str, j: int, f_max: float) -> tuple[np.ndarray, np.ndarray]:
    """Frequenze (Hz) dei poli e degli zeri poco smorzati (Im > |Re|) del canale, fino a ``f_max``."""
    from . import tf

    p = ml.poli
    z = tf.zeri_invarianti(ml.A, ml.Bu[:, j], ml.Cy[ml.riga(uscita)])
    sel = lambda r: np.sort(np.abs(r[(np.imag(r) > 0) & (np.imag(r) > np.abs(np.real(r)))])) / (2 * np.pi)  # noqa: E731
    fp, fz = sel(p), sel(z)
    return fp[fp <= f_max], fz[fz <= f_max]


RAPPORTO_CORRISPONDENZA = 2.0
"""Una caratteristica misurata corrisponde a un polo/zero del modello se le frequenze
differiscono meno di questo fattore."""


BANDE_FIT: tuple[tuple[str, float, float], ...] = (("all", 0.0, np.inf), ("1-5 Hz", 1.0, 5.0),
                                                   ("5-100 Hz", 5.0, 100.0), (">100 Hz", 100.0, np.inf))
"""Bande del fit: intera misura, vincolo a terra, dinamica del carico/telescopio, modi strutturali."""


def fit_nrmse(y: np.ndarray, y_modello: np.ndarray) -> float:
    r"""Fit in per cento come in System Identification Toolbox (``compare``, ``goodnessOfFit`` 'NRMSE', report
    di ``greyest``): :math:`100\,(1 - \lVert y - \hat y\rVert / \lVert y - \bar y\rVert)`, con
    :math:`\bar y` la media della misura. 100 = modello perfetto, 0 = non meglio della media, puo' essere
    negativo. Per dati in frequenza :math:`y` e :math:`\hat y` sono le risposte complesse sui punti di misura."""
    y, y_modello = np.asarray(y), np.asarray(y_modello)
    return float(100.0 * (1.0 - np.linalg.norm(y - y_modello) / np.linalg.norm(y - np.mean(y))))


def fit_bande(f: np.ndarray, y: np.ndarray, y_modello: np.ndarray) -> dict[str, float]:
    r""":func:`fit_nrmse` sull'intera misura e nelle bande :data:`BANDE_FIT` (``nan`` con meno di 3 punti)."""
    out = {}
    for nome, a, b in BANDE_FIT:
        m = (f >= a) & (f <= b)
        out[nome] = fit_nrmse(y[m], y_modello[m]) if m.sum() >= 3 else float("nan")
    return out


@dataclass(frozen=True)
class Confronto:
    """Misura contro modello a 7 GdL concentrato."""

    misura: Misura
    f_modello: np.ndarray
    G_modello: dict[str, np.ndarray]
    caratteristiche: tuple[Caratteristica, ...]
    scarto_livello_db: dict[str, float]
    """Mediana di (misura - modello) in dB nella banda di livello."""
    banda_livello: tuple[float, float]
    fit: dict[str, float] = field(default_factory=dict)
    """Fit NRMSE (per cento) del canale misurato (encoder) sui punti di misura, per banda (:data:`BANDE_FIT`)."""


def sistemi_per_misure(par: Parametri) -> dict[bool, cp.Sistema]:
    """Modello a 7 GdL concentrato con il telescopio (``True``) e senza (``False``)."""
    from . import payload as plm

    assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
            fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    pl = plm.costruisci_payload(par, assi[1])
    return {True: cc.costruisci_sistema_concentrato(par, pl, assi),
            False: cc.costruisci_sistema_concentrato(par, payload_senza_telescopio(pl), assi_senza_telescopio(par, assi))}


def assi_senza_telescopio(par: Parametri, assi: tuple) -> tuple:
    """Assi per la prova senza telescopio: nel modello ``hub`` il ramo del telescopio sparisce (il nodo
    P non avrebbe inerzia) e il centro del mozzo H fa da nodo P; altrimenti gli assi invariati. Lo stesso
    per il ramo L-P dell'azimut (``azimuth.telescope_attachment``)."""
    if par.mozzo is None:
        return assi
    psi = assi[0]
    if par.azimut.attacco_telescopio is not None:
        psi = fem.costruisci_asse(par.azimut, par.fem.forma_consistente, con_telescopio=False)
    return (psi, fem.costruisci_asse(par.polare, par.fem.forma_consistente, con_telescopio=False))


def confronta(par: Parametri, mis: Misura, sistemi: Mapping[bool, cp.Sistema] | None = None,
              f: np.ndarray | None = None) -> Confronto:
    cfg = par.sperimentale
    assert cfg is not None
    s = (sistemi or sistemi_per_misure(par))[mis.telescopio]
    f = np.logspace(np.log10(cfg.f_min), np.log10(cfg.f_max), 2000) if f is None else f
    ml = lz.linearizza(s.accoppiato, mis.theta0, 0.0)
    # vincolo verso terra (ground_constraint): solo nella risposta in frequenza; poli e zeri sotto
    # restano quelli strutturali del modello senza vincolo
    from . import vincolo_terra as vt

    ml_v = vt.applica(ml, par, mis.telescopio)
    G = frf_velocita(ml_v, mis.asse, f, cfg.Kt)
    # fit sui punti di misura (modello valutato alle frequenze misurate, canale encoder)
    G_mis_mod = frf_velocita(ml_v, mis.asse, mis.f, cfg.Kt)[canali(mis.asse)[0]]
    fit = fit_bande(mis.f, mis.G, G_mis_mod)
    j = 0 if mis.asse == "azimut" else 1
    car = []
    # risonanze: modi flessibili dell'asse eccitato (catena concentrata isolata in theta0, con k_g);
    # antirisonanze: zeri invarianti del canale encoder del modello accoppiato. Il drive chiude gli
    # anelli sull'encoder (indicazione dell'utente, 2026-10-16): la pianta misurata e' encoder/corrente,
    # il canale motore si mostra solo come riferimento.
    f_asse = fem.frequenze_proprie(*s.isolato_ridotto(mis.asse, mis.theta0)[:2])
    f_asse = f_asse[f_asse > 1.0]
    for tipo, fm, gm in estremi(mis.f, mis.modulo_db):
        migliore: tuple[float, str] | None = None
        for c in canali(mis.asse)[:1]:
            _, fz = poli_zeri_modello(ml, c, j, 2.0 * cfg.f_max)
            cand = f_asse if tipo == "risonanza" else fz
            if cand.size:
                k = int(np.argmin(np.abs(np.log(cand / fm))))
                if migliore is None or abs(np.log(cand[k] / fm)) < abs(np.log(migliore[0] / fm)):
                    migliore = (float(cand[k]), c)
        if migliore is not None and abs(np.log(migliore[0] / fm)) > np.log(RAPPORTO_CORRISPONDENZA):
            migliore = None  # nessun polo/zero del modello vicino: caratteristica non modellata
        car.append(Caratteristica(tipo, fm, gm, migliore[0] if migliore else None,
                                  migliore[1] if migliore else None))
    a, b = cfg.banda_livello
    m = (mis.f >= a) & (mis.f <= b)
    scarti = {c: float(np.median(mis.modulo_db[m] - np.interp(mis.f[m], f, 20 * np.log10(np.abs(Gc)))))
              for c, Gc in G.items()}
    return Confronto(misura=mis, f_modello=f, G_modello=G, caratteristiche=tuple(car),
                     scarto_livello_db=scarti, banda_livello=(a, b), fit=fit)


# ---------------------------------------------------------------------------
# Parametri equivalenti suggeriti dalle caratteristiche misurate
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StimaEquivalente:
    """Stima di primo ordine di un parametro che riprodurrebbe una caratteristica misurata."""

    misura: str
    grandezza: str
    valore: float
    unita: str
    modello: float | None
    base: str


def stime_equivalenti(conf: Confronto, s: cp.Sistema) -> list[StimaEquivalente]:
    """Stime di primo ordine (non un'identificazione):

    * risonanza a bassa frequenza (sotto l'antirisonanza carico-molla): rigidezza verso terra
      :math:`k_0 = J_{tot}(2\\pi f_0)^2`, assente nel modello (moto rigido libero);
    * antirisonanza carico-molla dell'azimut: :math:`f_a = \\sqrt{k_{\\psi,2}/J_L^*}/2\\pi` ->
      rapporto :math:`k_{\\psi,2}/J_L^*` misurato;
    * livello fra antirisonanza e prima risonanza (cedevolezza): scarto in dB -> fattore sulla
      cedevolezza fra motore e carico.
    """
    mis = conf.misura
    pl = s.payload
    out: list[StimaEquivalente] = []
    if mis.asse == "azimut":
        n = s.rid_psi.n  # 4 con il telescopio appeso a L: il carico su k_psi,2 e' L + p
        J = float(np.ones(n) @ s.rid_psi.M_r @ np.ones(n)) + pl.J_ax + float(pl.J_psi(mis.theta0))
        J_L = float(s.rid_psi.M_r[2:, 2:].sum()) + pl.J_ax + float(pl.J_psi(mis.theta0))
        k2 = float(s.rid_psi.k_catena[1])
    else:
        J = float(np.ones(4) @ s.rid_theta.M_r @ np.ones(4)) + (pl.I_yy if mis.telescopio else 0.0)
        J_L = float(s.rid_theta.M_r[1, 1]) + (pl.I_yy if mis.telescopio else 0.0)
        k2 = float(s.rid_theta.k_catena[1])
    antiris = [c for c in conf.caratteristiche if c.tipo == "antirisonanza"]
    f_a = min((c.f_misura for c in antiris if c.f_misura > 5.0), default=None)
    bassi = [c for c in conf.caratteristiche if c.tipo == "risonanza" and c.f_misura < (f_a or 20.0)
             and c.f_misura < 20.0]
    if bassi:
        f0 = max(bassi, key=lambda c: c.modulo_db).f_misura  # la risonanza bassa dominante
        out.append(StimaEquivalente(mis.etichetta, "rigidezza verso terra k_0", J * (2 * np.pi * f0) ** 2,
                                    "N m/rad", 0.0, f"risonanza misurata a {f0:.3g} Hz, J_tot = {J:.4g} kg m^2"))
    if f_a is not None and mis.telescopio:
        out.append(StimaEquivalente(mis.etichetta, "k_2 / J_L* (molla motore-carico / inerzia del carico)",
                                    (2 * np.pi * f_a) ** 2, "1/s^2", k2 / J_L,
                                    f"antirisonanza misurata a {f_a:.3g} Hz"))
    c = canali(mis.asse)[1 if mis.asse == "polare" else 0]
    out.append(StimaEquivalente(mis.etichetta, "fattore di livello (misura/modello)",
                                10 ** (conf.scarto_livello_db[c] / 20), "-", 1.0,
                                f"mediana {conf.scarto_livello_db[c]:+.1f} dB in {conf.banda_livello[0]:g}-"
                                f"{conf.banda_livello[1]:g} Hz, canale {c}"))
    return out
