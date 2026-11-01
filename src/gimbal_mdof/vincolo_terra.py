"""Rigidezza equivalente verso terra (sezione ``ground_constraint``, richiesta 2026-10-19).

Sotto pochi hertz le misure di pianta salgono di +20 dB/decade con fase positiva: l'asse non e'
un'inerzia libera ma un'inerzia su una molla dissipativa, l'attrito di pre-scorrimento dei
cuscinetti. Sotto il primo modo flessibile l'asse e' rigido e la rigidezza dinamica misurata vale

.. math::
    K(j\\omega) = \\frac{\\tau}{\\theta_e} = \\frac{j\\omega K_t}{G(j\\omega)}
               \\simeq k\\,(1 + j\\eta) + j\\omega c - J\\omega^2 ,

con :math:`G` la pianta misurata (velocita' dell'encoder su corrente) e :math:`J` l'inerzia
rigida dell'asse nel modello. :func:`stima` ricava :math:`k` e la perdita (isteretica
:math:`\\eta` o viscosa :math:`c`) ai minimi quadrati su modulo e fase; :func:`applica` aggiunge il
vincolo al modello lineare sul nodo dei cuscinetti (azimut: la forcella, nodo L / 9; polare: il
centro del mozzo H, o P se manca). Il termine isteretico rende complessa la rigidezza: il modello
vincolato serve solo per le risposte in frequenza, poli e zeri strutturali restano quelli del
modello senza vincolo.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal

import numpy as np

from . import linearize as lz
from .params import Parametri, VincoloTerra

Perdita = Literal["isteretica", "viscosa"]


@dataclass(frozen=True)
class Stima:
    """Esito di :func:`stima` nella banda :math:`[f_1, f_2]`."""

    k: float
    eta: float
    c: float
    J: float
    banda: tuple[float, float]
    rms_db: float
    rms_deg: float
    max_db: float
    max_deg: float

    @property
    def f0(self) -> float:
        """Frequenza propria dell'asse rigido sulla molla, Hz."""
        return float(np.sqrt(self.k / self.J) / (2 * np.pi))


def risposta(k: float, eta: float, c: float, J: float, f: np.ndarray, Kt: float) -> np.ndarray:
    """Pianta dell'asse rigido vincolato, :math:`j\\omega K_t / (k(1 + j\\eta) + j\\omega c - J\\omega^2)`."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    return 1j * w * Kt / (k * (1 + 1j * eta) + 1j * w * c - J * w**2)


def stima(f: np.ndarray, modulo_db: np.ndarray, fase_deg: np.ndarray, J: float, Kt: float, *,
          banda: tuple[float, float] = (1.0, 5.0), perdita: Perdita = "isteretica") -> Stima:
    """:math:`k` e perdita del vincolo ai minimi quadrati su modulo (dB) e fase (gradi / 10) nella banda,
    con l'inerzia rigida :math:`J` del modello fissata."""
    from scipy.optimize import least_squares

    m = (f >= banda[0]) & (f <= banda[1])
    fb, gd, pd = f[m], modulo_db[m], np.deg2rad(fase_deg[m])

    def residui(x: np.ndarray) -> np.ndarray:
        k, d = np.exp(x)
        G = risposta(k, d, 0.0, J, fb, Kt) if perdita == "isteretica" else risposta(k, 0.0, d, J, fb, Kt)
        return np.r_[20 * np.log10(np.abs(G)) - gd, np.angle(G * np.exp(-1j * pd), deg=True) / 10.0]

    d0 = 1.0 if perdita == "isteretica" else 20.0
    migliore = min((least_squares(residui, [np.log(k0), np.log(d0)]) for k0 in (50.0, 300.0, 1000.0, 3000.0)),
                   key=lambda r: r.cost)
    k, d = np.exp(migliore.x)
    r = residui(migliore.x)
    n = fb.size
    e_db, e_deg = r[:n], 10.0 * r[n:]
    return Stima(k=float(k), eta=float(d) if perdita == "isteretica" else 0.0,
                 c=float(d) if perdita == "viscosa" else 0.0, J=float(J), banda=banda,
                 rms_db=float(np.sqrt(np.mean(e_db**2))), rms_deg=float(np.sqrt(np.mean(e_deg**2))),
                 max_db=float(np.max(np.abs(e_db))), max_deg=float(np.max(np.abs(e_deg))))


def inerzia_rigida(ml: lz.ModelloLineare, asse: str, Kt: float, f: float = 1.0) -> float:
    """Inerzia rigida dell'asse vista dall'encoder, :math:`-\\mathrm{Re}\\,K(j\\omega)/\\omega^2` a bassa
    frequenza (molto sotto il primo modo, con :math:`k_g` trascurabile)."""
    from . import sperimentale as spm

    c = spm.canali(asse)[0]
    G = spm.frf_velocita(ml, asse, np.array([f]), Kt)[c][0]
    w = 2 * np.pi * f
    return float(-(1j * w * Kt / G).real / w**2)


def coordinata(ml: lz.ModelloLineare, asse: str) -> int:
    """Indice della coordinata del modello su cui agisce il vincolo: azimut, il nodo della forcella
    (ruolo L del modello ridotto, nodo 9 del FEM); polare, il centro del mozzo H fra i cuscinetti
    (P nella configurazione di base e nella prova senza telescopio)."""
    acc = ml.accoppiato
    if hasattr(acc, "base_psi"):  # FEM (base completa) o base modale: nomi dei gradi di liberta'
        nomi = list(acc.base_psi.nomi) + list(acc.base_theta.nomi)
        if asse == "azimut":
            return len(acc.base_psi.nomi) - 1
        for nome in ("theta_H", "theta_P"):
            if nome in nomi:
                return nomi.index(nome)
        raise ValueError("vincolo a terra: nessun nodo H o P nella base polare")
    st = acc.struttura
    if asse == "azimut":
        return st.j_psi_L
    ruoli = st.ruoli_theta
    return st.n_psi + (ruoli.index("h") if "h" in ruoli else ruoli.index("p"))


def applica(ml: lz.ModelloLineare, par: Parametri, telescopio: bool = True) -> lz.ModelloLineare:
    """Modello lineare con i vincoli verso terra di ``par.vincolo_terra`` (rigidezza complessa: solo per
    le risposte in frequenza). Senza la sezione, il modello invariato."""
    vt = par.vincolo_terra
    if vt is None:
        return ml
    K0, C0 = ml.K0.astype(complex), ml.C0.copy()
    for asse in ("azimut", "polare"):
        v: VincoloTerra | None = vt.per(asse, telescopio)
        if v is None:
            continue
        j = coordinata(ml, asse)
        K0[j, j] += v.k * (1 + 1j * v.eta)
        C0[j, j] += v.c
    return replace(ml, K0=K0, C0=C0)
