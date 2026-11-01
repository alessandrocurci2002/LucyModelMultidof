"""Modelli della forcella in forma modale, stimati ai minimi quadrati sulle parti reale e immaginaria.

Richiesta dell'utente (2026-10-31), su entrambi gli assi:

* telescopio collegato rigidamente all'anello (polare: fuso con il mozzo; azimut: fuso con la sommita' della forcella);
* cuscinetti con il vincolo lineare isteretico k_g (1 + j eta_g) (polare: fra mozzo e struttura dello statore; azimut:
  fra base della forcella e terra);
* formulazione modale: modi non smorzati di (M, K) con la parte reale del vincolo, smorzamento viscoso modale zeta_k
  per ogni modo (zeta_k <= 0.05; prima stima 0.01, vedi ``ZETA_MAX``) e la parte isteretica dei cuscinetti proiettata sulle coordinate modali (correzione di
  rango 1, esatta, con la formula di Sherman-Morrison);
* vincoli: inerzie e rigidezze del CAD entro +-20 % (prima +-15 %), K_t entro la tolleranza di catalogo +-10 %, rigidezze a terra, rigidezze introdotte dall'utente (base -
  sommita' della forcella) e parametri modali delle masse sullo statore liberi;
* errore: differenza fra le parti reali e fra le parti immaginarie della pianta (minimi quadrati, nessun logaritmo e
  nessuna perdita robusta).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares

from .modi_forcella import F_CORRENTE, CatenaAzimut, CatenaPolare, DatoFRF, Parametro, _misurata

ZETA_MAX = 0.05
"""Limite dello smorzamento viscoso modale di ogni modo (2026-10-31: alzato da 0.01, che limitava i modi delle masse
sullo statore)."""
ZETA_MIN = 1e-4
VARIAZIONE_CAD = 0.20
"""Variazione relativa massima delle rigidezze del CAD (2026-10-31: alzata da 0.15)."""
VARIAZIONE_CAD_J = 0.20
"""Variazione relativa massima delle inerzie del CAD (2026-10-31: alzata da 0.15)."""
VARIAZIONE_KT = 0.10
"""Tolleranza di K_t del datasheet Kollmorgen TBM(S)-7646-A (+-10 %): fattore g_Kt per asse (2026-10-31)."""
"""Variazione massima relativa delle inerzie e delle rigidezze del CAD."""
NODI_POLARE = ("e", "hp", "m", "r", "q1", "q2")
"""Encoder, mozzo con telescopio, rotore del motore, estremita' dell'albero, masse sullo statore."""
NODI_AZIMUT = ("e", "m", "L", "AP")
"""Encoder, rotore del motore, base della forcella, sommita' della forcella con telescopio."""


def molle(K: np.ndarray) -> tuple[list[tuple[int, int, float]], np.ndarray]:
    """Molle ``(i, j, k)`` (i < j) e rigidezze a terra per nodo di una matrice di rigidezza a catena."""
    n = K.shape[0]
    archi = [(i, j, float(-K[i, j])) for i in range(n) for j in range(i + 1, n) if K[i, j] != 0.0]
    return archi, K.sum(axis=1)


def _molla(K: np.ndarray, i: int, j: int | None, k: float) -> None:
    K[i, i] += k
    if j is not None:
        K[j, j] += k
        K[i, j] -= k
        K[j, i] -= k


@dataclass(frozen=True)
class Modo:
    f: float
    zeta: float
    residuo: float
    """R_k = (c^T phi_k)(phi_k^T e), 1/(kg m^2), con phi_k normalizzato sulla massa."""
    eta_cuscinetti: float
    """Fattore di perdita modale equivalente della parte isteretica dei cuscinetti, h w_k^2 / omega_k^2."""
    forma: dict[str, float]


@dataclass(frozen=True)
class ModelloModale:
    """Pianta WorkBench di un asse in forma modale (vedi il docstring del modulo)."""

    M: np.ndarray
    K: np.ndarray
    """Rigidezza reale, con la parte reale k_g v v^T del vincolo dei cuscinetti."""
    v: np.ndarray
    k_g: float
    eta_g: float
    e: np.ndarray
    c: np.ndarray
    zeta: np.ndarray
    """Smorzamento viscoso dei modi in ordine di frequenza crescente."""
    ritardo: float
    Kt: float
    nodi: tuple[str, ...]
    f_corrente: float = F_CORRENTE

    def __post_init__(self):
        w2, Phi = eigh(self.K, self.M)
        object.__setattr__(self, "omega", np.sqrt(np.maximum(w2, 0.0)))
        object.__setattr__(self, "Phi", Phi)

    def _modale(self, s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        b, g, w = self.Phi.T @ self.e, self.Phi.T @ self.c, self.Phi.T @ self.v
        D = s[:, None] ** 2 + 2 * self.zeta * self.omega * s[:, None] + self.omega ** 2
        return b, g, w, D

    def theta_tau(self, s: np.ndarray) -> np.ndarray:
        """theta_enc / tau: somma modale piu' la correzione di rango 1 della parte isteretica dei cuscinetti."""
        b, g, w, D = self._modale(s)
        h = 1j * self.k_g * self.eta_g
        G0 = (g * b / D).sum(axis=1)
        gw, wb, ww = (g * w / D).sum(axis=1), (w * b / D).sum(axis=1), (w * w / D).sum(axis=1)
        return G0 - h * gw * wb / (1 + h * ww)

    def frf(self, f: np.ndarray) -> np.ndarray:
        s = 2j * np.pi * np.asarray(f, dtype=float)
        return _misurata(self.theta_tau(s), s, self.Kt, self.ritardo, self.f_corrente)

    def contributi(self, f: np.ndarray) -> dict[str, np.ndarray]:
        """Termini modali R_k / D_k e correzione dei cuscinetti, ciascuno con i fattori della pianta WorkBench."""
        s = 2j * np.pi * np.asarray(f, dtype=float)
        b, g, w, D = self._modale(s)
        fatt = lambda x: _misurata(x, s, self.Kt, self.ritardo, self.f_corrente)  # noqa: E731
        out = {f"mode {k + 1} ({self.omega[k] / (2 * np.pi):.1f} Hz)": fatt(g[k] * b[k] / D[:, k])
               for k in range(len(self.omega))}
        out["bearing hysteresis (rank-1 term)"] = fatt(self.theta_tau(s) - (g * b / D).sum(axis=1))
        return out

    def modi(self) -> list[Modo]:
        b, g, w = self.Phi.T @ self.e, self.Phi.T @ self.c, self.Phi.T @ self.v
        out = []
        for k, om in enumerate(self.omega):
            x = self.Phi[:, k] / self.Phi[np.argmax(np.abs(self.Phi[:, k])), k]
            out.append(Modo(float(om / (2 * np.pi)), float(self.zeta[k]), float(g[k] * b[k]),
                            float(self.k_g * self.eta_g * w[k] ** 2 / om ** 2),
                            {n: float(xi) for n, xi in zip(self.nodi, x)}))
        return out


# ---------------------------------------------------------------------------
# Assi
# ---------------------------------------------------------------------------


def polare(cat: CatenaPolare, v: Mapping[str, float], Kt: float) -> ModelloModale:
    """H7 polare con il telescopio fuso con il mozzo, due masse sullo statore e cuscinetti fra mozzo e statore."""
    ruoli = list(cat.ruoli)
    ie, ih, im, ir, ip = (ruoli.index(x) for x in ("e", "h", "m", "r", "p"))
    J = np.diag(cat.M)
    archi, terra = molle(cat.K)
    rig = {frozenset((i, j)): k for i, j, k in archi}
    n = len(NODI_POLARE)
    M, K = np.zeros((n, n)), np.zeros((n, n))
    M[0, 0] = v["s_J_e"] * J[ie]
    M[1, 1] = v["s_J_hp"] * (J[ih] + J[ip])
    M[2, 2] = v["s_J_m"] * J[im]
    M[3, 3] = v["s_J_r"] * J[ir]
    _molla(K, 0, 1, v["s_k_1"] * rig[frozenset((ie, ih))])
    _molla(K, 1, 2, v["s_k_2"] * rig[frozenset((ih, im))])
    _molla(K, 2, 3, v["s_k_3"] * rig[frozenset((im, ir))])
    _molla(K, 1, None, terra[ih] + terra[ip])  # residuo a terra di H7 (~0.8 N m/rad), sul nodo fuso
    e, c, vv = np.zeros(n), np.zeros(n), np.zeros(n)
    e[2], c[0], vv[1] = 1.0, 1.0, 1.0
    for r in (1, 2):
        j = 3 + r
        Jr = 1.0 / abs(v[f"a_{r}"])
        M[j, j] = Jr
        _molla(K, j, None, Jr * (2 * np.pi * v[f"f_{r}"]) ** 2)
        e[j], c[j], vv[j] = -1.0, -np.sign(v[f"a_{r}"]), -1.0
    K += v["k_g"] * np.outer(vv, vv)
    zeta = np.array([v[f"zeta_{k}"] for k in range(1, n + 1)])
    return ModelloModale(M, K, vv, v["k_g"], v["eta_g"], e, c, zeta, v["tau_d"], Kt * v.get("g_Kt", 1.0), NODI_POLARE)


def azimut(cat: CatenaAzimut, v: Mapping[str, float], Kt: float) -> ModelloModale:
    """H7 azimut con il load diviso in base e sommita' della forcella e il telescopio fuso con la sommita'.

    Un solo fattore CAD scala il load e il telescopio (solo J_L e J_AP sono identificabili)."""
    n = len(NODI_AZIMUT)
    s = v["s_J_L"]
    M = np.diag([v["s_J_e"] * cat.J_e, v["s_J_m"] * cat.J_m, (1 - v["alpha"]) * s * cat.J_L,
                 v["alpha"] * s * cat.J_L + s * cat.J_p])
    K = np.zeros((n, n))
    _molla(K, 0, 1, v["s_k_1"] * cat.k[0])
    _molla(K, 1, 2, v["s_k_2"] * cat.k[1])
    _molla(K, 2, 3, v["k_LA"])
    vv, e, c = np.zeros(n), np.zeros(n), np.zeros(n)
    vv[2], e[1], c[0] = 1.0, 1.0, 1.0
    K += v["k_g"] * np.outer(vv, vv)
    zeta = np.array([v[f"zeta_{k}"] for k in range(1, n + 1)])
    return ModelloModale(M, K, vv, v["k_g"], v["eta_g"], e, c, zeta, v["tau_d"], Kt * v.get("g_Kt", 1.0), NODI_AZIMUT)


def _fattore_cad(cad: float, cad_J: float | None) -> Callable[[str, str], Parametro]:
    """Fattore CAD lineare: le inerzie (``s_J_*``) entro +-cad_J, le rigidezze entro +-cad."""
    def f(nome: str, nota: str) -> Parametro:
        v = cad if cad_J is None or not nome.startswith("s_J_") else cad_J
        return Parametro(nome, 1 - v, 1 + v, log=False, nota=nota)
    return f


def parametri_polari(zeta_max: float = ZETA_MAX, cad: float = VARIAZIONE_CAD,
                     cad_J: float | None = VARIAZIONE_CAD_J, kt: float = VARIAZIONE_KT) -> list[Parametro]:
    _cad = _fattore_cad(cad, cad_J)
    return [_cad("s_J_e", "CAD factor: encoder disc inertia"),
            _cad("s_J_hp", "CAD factor: hub + telescope inertia"),
            _cad("s_J_m", "CAD factor: motor rotor inertia"),
            _cad("s_J_r", "CAD factor: rear-end inertia"),
            _cad("s_k_1", "CAD factor: shaft encoder - hub"),
            _cad("s_k_2", "CAD factor: shaft hub - motor"),
            _cad("s_k_3", "CAD factor: shaft motor - rear end"),
            Parametro("k_g", 1.0, 1e6, unita="N m/rad", nota="bearings to the stator structure (free)"),
            Parametro("eta_g", 1e-3, 30.0, nota="hysteretic loss factor of the bearings (free)"),
            Parametro("tau_d", 1e-6, 2e-3, unita="s", nota="equivalent drive delay"),
            Parametro("f_1", 40.0, 80.0, unita="Hz", nota="stator mass 1"),
            Parametro("a_1", -50.0, 50.0, log=False, unita="1/(kg m^2)", nota="stator mass 1: kappa/J"),
            Parametro("f_2", 80.0, 140.0, unita="Hz", nota="stator mass 2"),
            Parametro("a_2", -50.0, 50.0, log=False, unita="1/(kg m^2)", nota="stator mass 2: kappa/J")
            ] + _kt(kt) + _zeta(len(NODI_POLARE), zeta_max)


def parametri_azimut(zeta_max: float = ZETA_MAX, cad: float = VARIAZIONE_CAD,
                     cad_J: float | None = VARIAZIONE_CAD_J, kt: float = VARIAZIONE_KT) -> list[Parametro]:
    _cad = _fattore_cad(cad, cad_J)
    return [_cad("s_J_e", "CAD factor: encoder inertia"),
            _cad("s_J_m", "CAD factor: motor rotor inertia"),
            _cad("s_J_L", "CAD factor: load and telescope inertia"),
            _cad("s_k_1", "CAD factor: shaft encoder - motor"),
            _cad("s_k_2", "CAD factor: shaft motor - fork base"),
            Parametro("alpha", 0.02, 0.98, log=False, nota="share of the CAD load inertia in the fork top"),
            Parametro("k_LA", 1e3, 1e8, unita="N m/rad", nota="fork base - fork top (free)"),
            Parametro("k_g", 1.0, 1e6, unita="N m/rad", nota="bearings to ground (free)"),
            Parametro("eta_g", 1e-3, 30.0, nota="hysteretic loss factor of the bearings (free)"),
            Parametro("tau_d", 1e-6, 2e-3, unita="s", nota="equivalent drive delay"),
            ] + _kt(kt) + _zeta(len(NODI_AZIMUT), zeta_max)


def _kt(kt: float) -> list[Parametro]:
    """Fattore su K_t entro la tolleranza di catalogo; nessun parametro se kt = 0."""
    return [Parametro("g_Kt", 1 - kt, 1 + kt, log=False, nota="factor on K_t (catalogue tolerance)")] if kt > 0 else []


def _zeta(n: int, zeta_max: float) -> list[Parametro]:
    return [Parametro(f"zeta_{k}", ZETA_MIN, zeta_max, nota=f"viscous modal damping of mode {k}")
            for k in range(1, n + 1)]


# ---------------------------------------------------------------------------
# Stima
# ---------------------------------------------------------------------------


@dataclass
class Esito:
    parametri: list[Parametro]
    valori: dict[str, float]
    std_rel: dict[str, float]
    costo: float
    nfev: int
    partenze: list[tuple[dict[str, float], float]] = field(default_factory=list)

    def al_limite(self, nome: str, tol: float = 1e-3) -> bool:
        p = next(q for q in self.parametri if q.nome == nome)
        lo, hi = p.limiti_x
        x = p.x(self.valori[nome])
        return bool(min(x - lo, hi - x) < tol * (hi - lo))


def residui(modello: Callable[[Mapping[str, float]], ModelloModale], parametri: Sequence[Parametro],
            dati: Sequence[DatoFRF]) -> Callable[[np.ndarray], np.ndarray]:
    """Differenze pesate delle parti reali e immaginarie, w_i [Re(G - G_mis), Im(G - G_mis)]."""

    def r(x: np.ndarray) -> np.ndarray:
        m = modello({p.nome: p.v(xi) for p, xi in zip(parametri, x)})
        out = []
        for d in dati:
            # i modelli con cuscinetti non lineari (funzione descrittiva) dipendono dalla misura, non solo da f
            e = (m.frf_dato(d) if hasattr(m, "frf_dato") else m.frf(d.f)) - d.G
            out += [e.real * d.peso, e.imag * d.peso]
        return np.concatenate(out)
    return r


def stima(modello: Callable[[Mapping[str, float]], ModelloModale], parametri: Sequence[Parametro],
          dati: Sequence[DatoFRF], partenze: Sequence[Mapping[str, float]], max_nfev: int = 3000) -> Esito:
    """Minimi quadrati con limiti da piu' partenze (``least_squares``, perdita quadratica); tiene la migliore."""
    parametri = list(parametri)
    fun = residui(modello, parametri, dati)
    lo, hi = (np.array(v) for v in zip(*(p.limiti_x for p in parametri)))
    migliore, storia = None, []
    for p0 in partenze:
        x0 = np.clip([p.x(p0[p.nome]) for p in parametri], lo + 1e-6 * (hi - lo), hi - 1e-6 * (hi - lo))
        r = least_squares(fun, x0, bounds=(lo, hi), loss="linear", max_nfev=max_nfev, x_scale="jac")
        storia.append((dict(p0), float(r.cost)))
        if migliore is None or r.cost < migliore.cost:
            migliore = r
    assert migliore is not None
    # covarianza sui soli parametri liberi: quelli al limite sono fissati dal vincolo (std = nan)
    libero = np.minimum(migliore.x - lo, hi - migliore.x) > 1e-3 * (hi - lo)
    J = migliore.jac[:, libero]
    m, n = J.shape
    s2 = 2 * migliore.cost / max(m - n, 1)
    sd = np.full(len(parametri), np.nan)
    sd[libero] = np.sqrt(np.maximum(np.diag(np.linalg.pinv(J.T @ J) * s2), 0.0))
    nulla = np.linalg.norm(J, axis=0) < 1e-9 * np.linalg.norm(J, axis=0).max()  # i dati non lo vedono
    sd[np.flatnonzero(libero)[nulla]] = np.inf
    valori = {p.nome: p.v(x) for p, x in zip(parametri, migliore.x)}
    std = {p.nome: float(sd[i] if p.log else sd[i] / max(abs(migliore.x[i]), 1e-300))
           for i, p in enumerate(parametri)}
    return Esito(parametri, valori, std, float(migliore.cost), int(migliore.nfev), storia)
