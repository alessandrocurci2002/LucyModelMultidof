r"""Vibrazioni della forcella nelle piante misurate: modello data-driven (richiesta dell'utente, 2026-10-29).

Le prove del 28-29/10 indicano che le risonanze polari sotto i 150 Hz sono vibrazioni della forcella: lo statore
del motore polare scarica la reazione della coppia sulla struttura su cui e' montato, la struttura vibra e la testina
dell'encoder, fissata sull'altro braccio, ne legge il moto. Per questo le stesse frequenze compaiono sull'azimut.

**Polare** (:class:`ModelloPolare`). Catena concentrata H7 (nodi e, h, m, r, p) con, in piu', :math:`N` masse
modali agganciate allo statore. Ogni massa e' un modo della struttura normalizzato sullo statore: inerzia
:math:`J_r`, molla :math:`k_r = J_r\omega_r^2` e smorzatore :math:`c_r = 2\zeta_r J_r\omega_r` verso terra,
rotazione :math:`q_r`. La reazione :math:`-\tau` dello statore agisce su tutte le masse; la testina dell'encoder
ruota di :math:`\sum_r \kappa_r q_r`, con :math:`\kappa_r = \pm 1` (in fase o in opposizione con lo statore nel
modo :math:`r`), e l'encoder legge la rotazione relativa

.. math:: \theta_{enc} = \theta_e - \sum_r \kappa_r\,q_r .

Dalla sola misura encoder/corrente si identifica il residuo :math:`a_r = \kappa_r / J_r` (segno e modulo), non
:math:`J_r` e :math:`\kappa_r` separati: :math:`J_r = 1/|a_r|` e' l'inerzia modale riferita allo statore con la
testina che si muove quanto lo statore. Con :math:`a_r < 0` e :math:`|a_r| J \gtrsim 1` il modo produce zeri a
fase non minima: la fase scende di 180 gradi senza risalire, come nelle misure fra 45 e 120 Hz.

I cuscinetti polari stanno sulla stessa struttura dello statore: il vincolo di pre-scorrimento
:math:`k_g(1 + j\eta)` (:mod:`gimbal_mdof.vincolo_terra`) collega il centro del mozzo H alla rotazione dello
statore :math:`\sum_r q_r`, non a terra.

**Azimut** (:class:`ModelloAzimut`). Il load L di H7 (forcella + alberi polari + mozzo + anello) non e' piu' rigido:
si divide in una base L (piatto e parte bassa della forcella, sui cuscinetti) e una sommita' A (bracci, mozzo,
anello) collegate da :math:`k_{LA}`, con :math:`J_L + J_A` pari all'inerzia CAD del load; il telescopio P e' appeso
alla sommita' con :math:`k_{\psi,3}`. Catena e - m - L - A - P, vincolo a terra sulla base.

Le risposte sono nelle unita' di WorkBench (velocita' dell'encoder su corrente, (rad/s)/A), con l'anello di corrente
del primo ordine (``IL.BW``) e un ritardo puro equivalente :math:`e^{-s\tau_d}`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares

from . import linearize as lz
from . import tf

F_CORRENTE = 1500.0
"""Banda dell'anello di corrente del drive (``IL.BW`` = 1.5 kHz, file .akd), Hz."""


# ---------------------------------------------------------------------------
# Catene di H7
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CatenaPolare:
    """Blocco polare di H7 (coordinate e, h, m, r, p) linearizzato in theta0 = 0."""

    M: np.ndarray
    C: np.ndarray
    K: np.ndarray
    ruoli: tuple[str, ...]

    def j(self, ruolo: str) -> int:
        return self.ruoli.index(ruolo)

    @property
    def k_b(self) -> float:
        """Rigidezza telescopio - mozzo (fra h e p), N m/rad."""
        return float(-self.K[self.j("h"), self.j("p")])

    @property
    def c_b(self) -> float:
        return float(-self.C[self.j("h"), self.j("p")])

    @property
    def J(self) -> float:
        """Inerzia rigida dell'asse, kg m^2."""
        return float(np.sum(self.M))


@dataclass(frozen=True)
class CatenaAzimut:
    """Catena di azimut di H7 (e - m - L - p): inerzie, molle e smorzatori degli elementi."""

    J_e: float
    J_m: float
    J_L: float
    J_p: float
    k: tuple[float, float, float]
    """:math:`k_{\\psi,1}` (e-m), :math:`k_{\\psi,2}` (m-L), :math:`k_{\\psi,3}` (L-p), N m/rad."""
    c: tuple[float, float, float]


def catene_h7(ml: lz.ModelloLineare) -> tuple[CatenaAzimut, CatenaPolare]:
    """Catene dei due assi dal modello concentrato H7 linearizzato (prodotti d'inerzia fra gli assi trascurati,
    ~1e-6)."""
    st = ml.accoppiato.struttura
    n_psi = st.n_psi
    ip = np.arange(n_psi, st.n)
    pol = CatenaPolare(*(X[np.ix_(ip, ip)].real.copy() for X in (ml.M0, ml.C0, ml.K0)), tuple(st.ruoli_theta))
    r = list(st.ruoli_psi)
    if r != ["e", "m", "L", "p"]:
        raise ValueError(f"catena di azimut inattesa: {r} (serve H7, telescopio appeso al load)")
    M, C, K = ml.M0, ml.C0, ml.K0
    kk = tuple(float(-K[i, i + 1]) for i in range(3))
    cc = tuple(float(-C[i, i + 1]) for i in range(3))
    az = CatenaAzimut(*(float(M[i, i]) for i in range(4)), kk, cc)
    return az, pol


def _uscita(M: np.ndarray, C: np.ndarray, K: np.ndarray, e: np.ndarray, c: np.ndarray, s: np.ndarray) -> np.ndarray:
    """:math:`\\mathbf c^T(\\mathbf M s^2 + \\mathbf C s + \\mathbf K)^{-1}\\mathbf e` su tutte le ``s``."""
    P = M[None] * s[:, None, None] ** 2 + C[None] * s[:, None, None] + K[None]
    x = np.linalg.solve(P, np.broadcast_to(e.astype(complex), (s.size, e.size))[..., None])[..., 0]
    return x @ c


def _misurata(theta_tau: np.ndarray, s: np.ndarray, Kt: float, ritardo: float, f_corrente: float) -> np.ndarray:
    """Pianta WorkBench: :math:`s K_t\\,\\theta/\\tau` con anello di corrente e ritardo."""
    wi = 2 * np.pi * f_corrente
    return s * Kt * theta_tau * wi / (s + wi) * np.exp(-s * ritardo)


def _radici(M: np.ndarray, C: np.ndarray, K: np.ndarray, e: np.ndarray, c: np.ndarray,
            f_max: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Poli e zeri (f [Hz], zeta) con Im > 0 fino a ``f_max`` del canale ``c^T x / tau``, rigidezza reale."""
    n = M.shape[0]
    Mi = np.linalg.inv(M)
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-Mi @ K.real, -Mi @ C]])
    b = np.r_[np.zeros(n), Mi @ e]
    cc = np.r_[c, np.zeros(n)]

    def sel(r: np.ndarray) -> list[tuple[float, float]]:
        r = r[(np.imag(r) > 1e-9) & (np.abs(r) <= 2 * np.pi * f_max)]
        return sorted((float(abs(x) / (2 * np.pi)), float(-x.real / abs(x))) for x in r)
    return sel(np.linalg.eigvals(A)), sel(tf.zeri_invarianti(A, b, cc))


# ---------------------------------------------------------------------------
# Polare: masse modali sullo statore
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModoStatore:
    """Massa modale agganciata allo statore polare."""

    f: float
    """Frequenza propria della massa sulla sua molla, Hz."""
    zeta: float
    a: float
    """Residuo :math:`\\kappa/J`, 1/(kg m^2): il segno e' quello di :math:`\\kappa`."""

    @property
    def J(self) -> float:
        return 1.0 / abs(self.a)

    @property
    def kappa(self) -> float:
        return float(np.sign(self.a))

    @property
    def k(self) -> float:
        return self.J * (2 * np.pi * self.f) ** 2

    @property
    def c(self) -> float:
        return 2 * self.zeta * self.J * 2 * np.pi * self.f


@dataclass(frozen=True)
class ModelloPolare:
    catena: CatenaPolare
    k_b: float
    c_b: float
    k_g: float
    eta_g: float
    ritardo: float
    modi: tuple[ModoStatore, ...]
    Kt: float
    f_corrente: float = F_CORRENTE
    cuscinetti_su_statore: bool = True
    """Vincolo di pre-scorrimento fra il mozzo e lo statore (True) o fra il mozzo e terra (False, come in H7)."""

    def matrici(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """``(M, C, K, e, c)``: coordinate [catena H7, q_1..q_N], ``K`` complessa (vincolo isteretico),
        ingresso :math:`\\tau` al motore e reazione sulle masse, uscita encoder relativa alla testina."""
        cat, N = self.catena, len(self.modi)
        n0 = cat.M.shape[0]
        n = n0 + N
        M, C, K = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n), complex)
        M[:n0, :n0], C[:n0, :n0], K[:n0, :n0] = cat.M, cat.C, cat.K
        h, p = cat.j("h"), cat.j("p")
        for X, nuovo, vecchio in ((K, self.k_b, cat.k_b), (C, self.c_b, cat.c_b)):
            d = nuovo - vecchio
            X[h, h] += d
            X[p, p] += d
            X[h, p] -= d
            X[p, h] -= d
        e, c = np.zeros(n), np.zeros(n)
        e[cat.j("m")], c[cat.j("e")] = 1.0, 1.0
        v = np.zeros(n)  # allungamento del vincolo dei cuscinetti
        v[h] = 1.0
        for r, m in enumerate(self.modi):
            j = n0 + r
            M[j, j], K[j, j], C[j, j] = m.J, m.k, m.c
            e[j] = -1.0
            c[j] = -m.kappa
            if self.cuscinetti_su_statore:
                v[j] = -1.0
        K += self.k_g * (1 + 1j * self.eta_g) * np.outer(v, v)
        return M, C, K, e, c

    def frf(self, f: np.ndarray) -> np.ndarray:
        """Velocita' dell'encoder su corrente, (rad/s)/A."""
        s = 2j * np.pi * np.asarray(f, dtype=float)
        return _misurata(_uscita(*self.matrici(), s), s, self.Kt, self.ritardo, self.f_corrente)

    def contributi(self, f: np.ndarray) -> dict[str, np.ndarray]:
        """Pianta scomposta nei termini della lettura dell'encoder (la somma e' :meth:`frf`): rotazione del disco
        :math:`\\theta_e` e, per ogni massa, :math:`-\\kappa_r q_r`."""
        M, C, K, e, c = self.matrici()
        n0 = self.catena.M.shape[0]
        s = 2j * np.pi * np.asarray(f, dtype=float)
        out = {}
        for nome, j in [("rotor (encoder disc)", self.catena.j("e"))] + [
                (f"mass {r + 1} ({m.f:.1f} Hz)", n0 + r) for r, m in enumerate(self.modi)]:
            cj = np.zeros_like(c)
            cj[j] = c[j]
            out[nome] = _misurata(_uscita(M, C, K, e, cj, s), s, self.Kt, self.ritardo, self.f_corrente)
        return out

    def radici(self, f_max: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        """Poli e zeri (f, zeta) del canale encoder fino a ``f_max``, con la sola rigidezza del vincolo."""
        return _radici(*self.matrici(), f_max)


def modello_polare(catena: CatenaPolare, v: Mapping[str, float], n_modi: int, Kt: float,
                   cuscinetti_su_statore: bool = True) -> ModelloPolare:
    """:class:`ModelloPolare` dai valori dei parametri di :func:`parametri_polari`."""
    modi = tuple(ModoStatore(v[f"f_{r}"], v[f"zeta_{r}"], v[f"a_{r}"]) for r in range(1, n_modi + 1))
    return ModelloPolare(catena, v["k_b"], v["c_b"], v["k_g"], v["eta_g"], v["tau_d"], modi, Kt,
                         cuscinetti_su_statore=cuscinetti_su_statore)


# ---------------------------------------------------------------------------
# Azimut: load flessibile
# ---------------------------------------------------------------------------


NODI_AZIMUT = ("e", "m", "L", "A", "P")


@dataclass(frozen=True)
class ModelloAzimut:
    catena: CatenaAzimut
    alfa: float
    """Frazione dell'inerzia CAD del load nella sommita' della forcella, :math:`J_A / (J_L + J_A)`."""
    k_LA: float
    c_LA: float
    k_3: float
    c_3: float
    k_2: float
    c_2: float
    k_g: float
    eta_g: float
    ritardo: float
    Kt: float
    f_corrente: float = F_CORRENTE

    @property
    def J_A(self) -> float:
        return self.alfa * self.catena.J_L

    @property
    def J_L(self) -> float:
        return (1.0 - self.alfa) * self.catena.J_L

    def matrici(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        cat = self.catena
        M = np.diag([cat.J_e, cat.J_m, self.J_L, self.J_A, cat.J_p])
        K, C = np.zeros((5, 5), complex), np.zeros((5, 5))
        for X, valori in ((K, (cat.k[0], self.k_2, self.k_LA, self.k_3)), (C, (cat.c[0], self.c_2, self.c_LA, self.c_3))):
            for i, x in enumerate(valori):
                X[i, i] += x
                X[i + 1, i + 1] += x
                X[i, i + 1] -= x
                X[i + 1, i] -= x
        K[2, 2] += self.k_g * (1 + 1j * self.eta_g)
        e, c = np.zeros(5), np.zeros(5)
        e[1], c[0] = 1.0, 1.0
        return M, C, K, e, c

    def frf(self, f: np.ndarray) -> np.ndarray:
        s = 2j * np.pi * np.asarray(f, dtype=float)
        return _misurata(_uscita(*self.matrici(), s), s, self.Kt, self.ritardo, self.f_corrente)

    def radici(self, f_max: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        return _radici(*self.matrici(), f_max)

    def modi(self, f_max: float) -> list[tuple[float, float, dict[str, float]]]:
        """``(f, zeta, forma)`` dei modi fino a ``f_max``; forma = ampiezze reali dei nodi normalizzate sulla
        massima in modulo (rigidezza reale del vincolo)."""
        M, C, K, _, _ = self.matrici()
        n = 5
        Mi = np.linalg.inv(M)
        A = np.block([[np.zeros((n, n)), np.eye(n)], [-Mi @ K.real, -Mi @ C]])
        w, V = np.linalg.eig(A)
        out = []
        for i in np.flatnonzero((w.imag > 1e-9) & (np.abs(w) <= 2 * np.pi * f_max)):
            x = V[:n, i]
            x = x / x[np.argmax(np.abs(x))]
            out.append((float(abs(w[i]) / (2 * np.pi)), float(-w[i].real / abs(w[i])),
                        {nodo: float(x[j].real) for j, nodo in enumerate(NODI_AZIMUT)}))
        return sorted(out, key=lambda t: t[0])


def modello_azimut(catena: CatenaAzimut, v: Mapping[str, float], Kt: float) -> ModelloAzimut:
    return ModelloAzimut(catena, v["alpha"], v["k_LA"], v["c_LA"], v["k_3"], v["c_3"], v["k_2"], v["c_2"],
                         v["k_g"], v["eta_g"], v["tau_d"], Kt)


# ---------------------------------------------------------------------------
# Stima ai minimi quadrati
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Parametro:
    nome: str
    lo: float
    hi: float
    log: bool = True
    """Ricerca in scala logaritmica (parametri positivi) o lineare (residui con segno)."""
    unita: str = ""
    nota: str = ""

    def x(self, v: float) -> float:
        v = float(np.clip(v, self.lo, self.hi))
        return float(np.log(v)) if self.log else v

    def v(self, x: float) -> float:
        return float(np.exp(x)) if self.log else float(x)

    @property
    def limiti_x(self) -> tuple[float, float]:
        return (np.log(self.lo), np.log(self.hi)) if self.log else (self.lo, self.hi)


@dataclass(frozen=True)
class DatoFRF:
    """Punti di una misura nella banda del fit, con il peso per punto."""

    etichetta: str
    f: np.ndarray
    G: np.ndarray
    peso: np.ndarray


def dato(etichetta: str, f: np.ndarray, G: np.ndarray, banda: tuple[float, float]) -> DatoFRF:
    """Punti in ``banda``; il peso :math:`\\sqrt{\\Delta\\ln f / 0.01}` integra l'errore in log f, cosi' le griglie
    piu' fitte non contano di piu' (1 per una griglia all'1 %)."""
    k = (f >= banda[0]) & (f <= banda[1])
    lf = np.log(f[k])
    return DatoFRF(etichetta, f[k], G[k], np.sqrt(np.gradient(lf) / 0.01))


@dataclass
class Esito:
    parametri: list[Parametro]
    valori: dict[str, float]
    std_rel: dict[str, float]
    """Deviazione standard relativa dalla jacobiana (:math:`\\sigma_x` per i parametri in log, :math:`\\sigma/|v|`
    per quelli lineari); ``inf`` se il parametro non e' determinato."""
    costo: float
    nfev: int
    partenze: list[tuple[dict[str, float], float]] = field(default_factory=list)
    """Valori iniziali e costo finale di ogni partenza."""

    def al_limite(self, nome: str, tol: float = 1e-3) -> bool:
        p = next(q for q in self.parametri if q.nome == nome)
        lo, hi = p.limiti_x
        x = p.x(self.valori[nome])
        return bool(min(x - lo, hi - x) < tol * max(1.0, abs(hi - lo)))


PERDITA = "soft_l1"
SCALA_PERDITA = 0.1
"""Perdita robusta sull'errore logaritmico complesso: lineare oltre ~0.1 (circa 0.9 dB o 6 gradi), cosi' le valli
profonde con poca coerenza non dominano il fit."""


def residui(modello: Callable[[Mapping[str, float]], object], parametri: Sequence[Parametro],
            dati: Sequence[DatoFRF]) -> Callable[[np.ndarray], np.ndarray]:
    """Errore :math:`\\ln(G_{mod}/G_{mis})` (parte reale: modulo in neper; immaginaria: fase in rad) pesato."""

    def r(x: np.ndarray) -> np.ndarray:
        m = modello({p.nome: p.v(xi) for p, xi in zip(parametri, x)})
        out = []
        for d in dati:
            e = np.log(m.frf(d.f) / d.G)
            out += [e.real * d.peso, e.imag * d.peso]
        return np.concatenate(out)
    return r


def stima(modello: Callable[[Mapping[str, float]], object], parametri: Sequence[Parametro],
          dati: Sequence[DatoFRF], partenze: Sequence[Mapping[str, float]], max_nfev: int = 800) -> Esito:
    """Minimi quadrati da piu' partenze (``least_squares``, perdita :data:`PERDITA`); tiene la migliore."""
    parametri = list(parametri)
    fun = residui(modello, parametri, dati)
    lo, hi = (np.array(v) for v in zip(*(p.limiti_x for p in parametri)))
    migliore, storia = None, []
    for p0 in partenze:
        x0 = np.clip([p.x(p0[p.nome]) for p in parametri], lo + 1e-9 * (hi - lo), hi - 1e-9 * (hi - lo))
        r = least_squares(fun, x0, bounds=(lo, hi), loss=PERDITA, f_scale=SCALA_PERDITA, max_nfev=max_nfev,
                          x_scale="jac")
        storia.append((dict(p0), float(r.cost)))
        if migliore is None or r.cost < migliore.cost:
            migliore = r
    assert migliore is not None
    J = migliore.jac
    m, n = J.shape
    s2 = 2 * migliore.cost / max(m - n, 1)
    cov = np.linalg.pinv(J.T @ J) * s2
    sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
    valori = {p.nome: p.v(x) for p, x in zip(parametri, migliore.x)}
    std = {p.nome: float(sd[i] if p.log else sd[i] / max(abs(migliore.x[i]), 1e-300))
           for i, p in enumerate(parametri)}
    return Esito(parametri, valori, std, float(migliore.cost), int(migliore.nfev), storia)


def parametri_polari(f_modi: Sequence[float], catena: CatenaPolare) -> list[Parametro]:
    """Parametri del modello polare: attacco del telescopio, vincolo dei cuscinetti, ritardo e, per ogni massa,
    frequenza (entro il 15 % del valore iniziale sotto i 20 Hz, 10 % sopra), smorzamento e residuo."""
    out = [Parametro("k_b", 1e4, 1e9, unita="N m/rad", nota="telescope - hub"),
           Parametro("c_b", 1e-3, 1e4, unita="N m s/rad", nota="telescope - hub"),
           Parametro("k_g", 1.0, 1e5, unita="N m/rad", nota="bearings (pre-sliding) to the stator structure"),
           Parametro("eta_g", 1e-3, 30.0, nota="loss factor of the bearing constraint"),
           Parametro("tau_d", 1e-6, 1e-3, unita="s", nota="equivalent drive delay")]
    for r, f0 in enumerate(f_modi, start=1):
        df = 0.15 if f0 < 20.0 else 0.10
        out += [Parametro(f"f_{r}", f0 * (1 - df), f0 * (1 + df), unita="Hz", nota=f"stator mass {r}"),
                Parametro(f"zeta_{r}", 1e-3, 0.5, nota=f"stator mass {r}"),
                Parametro(f"a_{r}", -50.0, 50.0, log=False, unita="1/(kg m^2)",
                          nota=f"stator mass {r}: kappa/J")]
    return out


def parametri_azimut(catena: CatenaAzimut) -> list[Parametro]:
    k2 = catena.k[1]
    return [Parametro("alpha", 0.02, 0.95, nota="share of the CAD load inertia in the fork top"),
            Parametro("k_LA", 1e3, 1e8, unita="N m/rad", nota="fork base - fork top"),
            Parametro("c_LA", 1e-3, 1e3, unita="N m s/rad", nota="fork base - fork top"),
            Parametro("k_3", 1e3, 1e8, unita="N m/rad", nota="fork top - telescope"),
            Parametro("c_3", 1e-3, 1e3, unita="N m s/rad", nota="fork top - telescope"),
            Parametro("k_2", 0.5 * k2, 3.0 * k2, unita="N m/rad", nota="motor - fork base (shaft)"),
            Parametro("c_2", 1e-3, 1e2, unita="N m s/rad", nota="motor - fork base (shaft)"),
            Parametro("k_g", 10.0, 1e5, unita="N m/rad", nota="bearings to ground"),
            Parametro("eta_g", 1e-3, 10.0, nota="loss factor of the bearing constraint"),
            Parametro("tau_d", 1e-6, 2e-3, unita="s", nota="equivalent drive delay")]


def con_modi(valori: Mapping[str, float], nuovi: Sequence[tuple[float, float, float]],
             primo: int = 1) -> dict[str, float]:
    """Copia di ``valori`` con le masse ``(f, zeta, a)`` dalla numero ``primo``."""
    out = dict(valori)
    for r, (f, z, a) in enumerate(nuovi, start=primo):
        out[f"f_{r}"], out[f"zeta_{r}"], out[f"a_{r}"] = f, z, a
    return out
