"""Tabelle e analisi parametriche (sezioni 8.3 e 11).

* validazione della riduzione: frequenze proprie del modello completo contro il
  ridotto, per ogni combinazione (concentrata/consistente) x
  (guyan/guyan_lumped), e frequenza oltre la quale il ridotto non e' piu'
  rappresentativo;
* poli e zeri del modello accoppiato per ogni punto di lavoro;
* frequenze in funzione di :math:`\\theta_0` e FRF incrociate in funzione di
  :math:`\\Omega_0`.

Ogni CSV porta, se i dati sono segnaposto, la colonna/nota "risultato da
parametri segnaposto".
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from . import fem
from . import linearize as lz
from . import reduction as rd
from . import tf
from .coupled import Sistema
from .params import NOTA_REALISMO

SOGLIA_DB: float = 3.0
"""Scarto di modulo oltre il quale il ridotto non e' piu' rappresentativo."""


def salva_csv(percorso: str | Path, intestazione: Sequence[str], righe: Iterable[Sequence[object]],
              *, segnaposto: bool = False) -> Path:
    percorso = Path(percorso)
    percorso.parent.mkdir(parents=True, exist_ok=True)
    with percorso.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        if segnaposto:
            w.writerow([f"# AVVISO: {NOTA_REALISMO}"])
        w.writerow(intestazione)
        w.writerows(righe)
    return percorso


# ---------------------------------------------------------------------------
# Frequenze proprie dei modelli isolati e validazione della riduzione
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConfrontoRiduzione:
    asse: str
    modello_massa: str
    forma_massa: str
    f_full: np.ndarray
    f_red: np.ndarray
    f_cb: dict[int, np.ndarray]
    f_limite: float
    """Frequenza oltre la quale la FRF coppia -> encoder del ridotto si scosta di
    piu' di :data:`SOGLIA_DB` dB dal completo (``nan`` se mai nella banda)."""

    @property
    def errore_relativo(self) -> np.ndarray:
        n = len(self.f_red)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(np.abs(self.f_full[:n]) > 0, self.f_red / self.f_full[:n] - 1.0, 0.0)


def frf_encoder_isolato(M: np.ndarray, C: np.ndarray, K: np.ndarray, i_motore: int,
                        i_uscita: int, w: np.ndarray) -> np.ndarray:
    n = M.shape[0]
    b = np.zeros((n, 1))
    b[i_motore] = 1.0
    c = np.zeros((1, n))
    c[0, i_uscita] = 1.0
    return tf.frf_secondo_ordine(M, C, K, b, c, w)[:, 0, 0]


def frequenza_limite(w: np.ndarray, G_full: np.ndarray, G_red: np.ndarray,
                     soglia_db: float = SOGLIA_DB) -> float:
    """Prima frequenza (Hz) in cui il modulo del ridotto si scosta oltre la soglia."""
    d = np.abs(20 * np.log10(np.abs(G_red) / np.abs(G_full)))
    fuori = np.nonzero(d > soglia_db)[0]
    return float(w[fuori[0]] / (2 * np.pi)) if fuori.size else float("nan")


def confronta_riduzione(s: Sistema, asse: str, theta0: float, w: np.ndarray,
                        n_cb: Sequence[int] = (1, 2, 4)) -> ConfrontoRiduzione:
    """Frequenze e FRF coppia -> encoder: completo contro ridotto, asse isolato."""
    M, K, C = s.isolato_completo(asse, theta0)
    Mr, Kr, Cr = s.isolato_ridotto(asse, theta0)
    mod = s.asse_psi if asse == "azimut" else s.asse_theta
    rid = s.rid_psi if asse == "azimut" else s.rid_theta
    f_full = fem.frequenze_proprie(M, K)
    f_red = fem.frequenze_proprie(Mr, Kr)
    f_cb = {}
    for n in n_cb:
        _, Mc, Kc = rd.riduci_craig_bampton(M, K, mod.masters, n)
        f_cb[n] = fem.frequenze_proprie(Mc, Kc)
    G_full = frf_encoder_isolato(M, C, K, mod.i_motore, mod.i_encoder, w)
    G_red = frf_encoder_isolato(Mr, Cr, Kr, rid.i_motore, rid.i_encoder, w)
    return ConfrontoRiduzione(asse=asse, modello_massa=s.modello_massa, forma_massa=s.forma_massa,
                              f_full=f_full, f_red=f_red, f_cb=f_cb,
                              f_limite=frequenza_limite(w, G_full, G_red))


def righe_validazione(confronti: Sequence[ConfrontoRiduzione]) -> list[list[object]]:
    righe = []
    for c in confronti:
        for k in range(len(c.f_red)):
            righe.append([c.asse, c.modello_massa, c.forma_massa, k + 1, f"{c.f_full[k]:.6g}",
                          f"{c.f_red[k]:.6g}", f"{c.errore_relativo[k]:+.3e}",
                          *[f"{c.f_cb[n][k]:.6g}" for n in sorted(c.f_cb)],
                          f"{c.f_limite:.6g}"])
    return righe


def intestazione_validazione(confronti: Sequence[ConfrontoRiduzione]) -> list[str]:
    n_cb = sorted(confronti[0].f_cb) if confronti else []
    return ["asse", "modello_massa", "forma_massa", "modo", "f_completo_Hz", "f_ridotto_Hz",
            "errore_relativo", *[f"f_CB{n}_Hz" for n in n_cb], f"f_limite_{SOGLIA_DB:g}dB_Hz"]


# ---------------------------------------------------------------------------
# Poli e zeri del modello accoppiato
# ---------------------------------------------------------------------------

INTESTAZIONE_RADICI = ("modello_massa", "forma_massa", "theta0_deg", "Omega0_deg_s", "grandezza",
                       "Re", "Im", "f_Hz", "zeta", "tipo")


def righe_poli_zeri(ml: lz.ModelloLineare, modello_massa: str, forma_massa: str,
                    uscite: Sequence[str] = ("psi_e", "theta_e")) -> list[list[object]]:
    """Poli del sistema e zeri invarianti della matrice coppia -> encoder."""
    base = [modello_massa, forma_massa, f"{np.rad2deg(ml.theta0):g}", f"{np.rad2deg(ml.Omega0):g}"]
    righe = []
    for r in tf.classifica(ml.poli):
        if np.imag(r.valore) < 0:
            continue
        righe.append(base + ["polo", f"{np.real(r.valore):.6g}", f"{np.imag(r.valore):.6g}",
                             f"{r.f_hz:.6g}", "" if r.zeta is None else f"{r.zeta:.4g}", r.tipo])
    ingressi = ("tau_psi", "tau_theta")
    for u in uscite:
        i = ml.riga(u)
        for j, ing in enumerate(ingressi):
            z = tf.zeri_invarianti(ml.A, ml.Bu[:, j], ml.Cy[i])
            for r in tf.classifica(z):
                if np.imag(r.valore) < 0:
                    continue
                righe.append(base + [f"zero {u}/{ing}", f"{np.real(r.valore):.6g}",
                                     f"{np.imag(r.valore):.6g}", f"{r.f_hz:.6g}",
                                     "" if r.zeta is None else f"{r.zeta:.4g}", r.tipo])
    return righe


def frequenze_modi(ml: lz.ModelloLineare) -> np.ndarray:
    """Frequenze (Hz) dei poli complessi del modello accoppiato, crescenti."""
    p = ml.poli
    complessi = p[np.imag(p) > 1e-6 * max(1.0, float(np.max(np.abs(p))))]
    return np.sort(np.abs(complessi)) / (2 * np.pi)


# ---------------------------------------------------------------------------
# Sensibilita' (sezione 11, punto 4)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SweepTheta:
    theta0: np.ndarray
    omega0: tuple[float, ...]
    frequenze: dict[float, np.ndarray]
    """Per ogni Omega_0: matrice (len(theta0), n_modi) delle frequenze dei modi
    complessi (Hz); NaN dove il numero di modi complessi cambia."""
    reali: dict[float, np.ndarray]
    """Massimo polo reale (rad/s) per ogni theta_0: > 0 se instabile."""


def sweep_theta0(s: Sistema, theta_grid: np.ndarray, omega_list: Sequence[float]) -> SweepTheta:
    freq: dict[float, np.ndarray] = {}
    reali: dict[float, np.ndarray] = {}
    for Om in omega_list:
        righe, rmax = [], []
        for th in theta_grid:
            ml = lz.linearizza(s.accoppiato, float(th), float(Om))
            righe.append(frequenze_modi(ml))
            p = ml.poli
            rmax.append(float(np.max(np.real(p[np.abs(np.imag(p)) < 1e-6 * max(1.0, np.max(np.abs(p)))]),
                                     initial=-np.inf)))
        n = max(len(r) for r in righe)
        F = np.full((len(theta_grid), n), np.nan)
        for i, r in enumerate(righe):
            F[i, n - len(r):] = r
        freq[float(Om)] = F
        reali[float(Om)] = np.array(rmax)
    return SweepTheta(theta0=np.asarray(theta_grid), omega0=tuple(float(o) for o in omega_list),
                      frequenze=freq, reali=reali)


def frf_incrociate(s: Sistema, theta0: float, omega_list: Sequence[float], w: np.ndarray) -> dict[float, np.ndarray]:
    """FRF coppia -> encoder (2x2) per ciascun :math:`\\Omega_0`."""
    return {float(Om): lz.linearizza(s.accoppiato, theta0, float(Om)).frf(w, ["psi_e", "theta_e"])
            for Om in omega_list}
