"""Riduzione di ciascun asse ai nodi master (sezione 8).

Condensazione statica di Guyan: partizionati i GdL in master (m) e slave (s),

.. math::
    \\mathbf T = \\begin{bmatrix}\\mathbb 1\\\\ -\\mathbf K_{ss}^{-1}\\mathbf K_{sm}\\end{bmatrix},
    \\qquad \\mathbf M_r = \\mathbf T^{\\mathsf T}\\mathbf M\\mathbf T,\\quad
    \\mathbf K_r = \\mathbf T^{\\mathsf T}\\mathbf K\\mathbf T,\\quad
    \\mathbf C_r = \\mathbf T^{\\mathsf T}\\mathbf C\\mathbf T.

Per una catena la condensazione ha forma chiusa (rigidezze in serie fra master
consecutivi, slave interpolati con pesi sulla cedevolezza), che il modulo
calcola anche indipendentemente come target di regressione.

``T`` dipende solo dalla rigidezza degli alberi: la rigidezza gravitazionale e
le inerzie del telescopio e del load stanno su nodi master e passano invariate
attraverso la condensazione. Per questo la stessa ``T`` serve sia al modello
isolato sia a quello accoppiato.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal, Sequence

import numpy as np
import scipy.linalg as sla

from .fem import ModelloAsse, ModelloMassa

FormaMassa = Literal["guyan", "guyan_lumped"]
FORME_MASSA: Final[tuple[FormaMassa, ...]] = ("guyan", "guyan_lumped")

RUOLI: Final[dict[str, tuple[str, ...]]] = {
    "azimut": ("e", "m", "L"),
    "polare": ("e", "p", "m", "r"),
}
"""Ruoli dei master nella topologia di riferimento (sezione 4.2)."""


# ---------------------------------------------------------------------------
# Guyan
# ---------------------------------------------------------------------------


def trasformazione_guyan(K: np.ndarray, masters: Sequence[int]) -> np.ndarray:
    """Matrice di condensazione :math:`\\mathbf T` (n x m), righe nell'ordine globale."""
    n = K.shape[0]
    m = list(masters)
    s = [i for i in range(n) if i not in m]
    T = np.zeros((n, len(m)))
    T[m, np.arange(len(m))] = 1.0
    if s:
        T[s, :] = -np.linalg.solve(K[np.ix_(s, s)], K[np.ix_(s, m)])
    return T


@dataclass(frozen=True)
class ModelloRidotto:
    """Catena ridotta di un asse.

    Attributes
    ----------
    etichette
        Etichette dei master (``"1"``, ``"P"``, ...), nell'ordine della catena.
    ruoli
        Ruolo di ciascun master (``e``, ``m``, ``L`` / ``e``, ``p``, ``m``, ``r``)
        se i master sono quelli di riferimento, altrimenti le etichette.
    T
        Matrice di condensazione (n x m).
    M_r, K_r, C_b_r
        Massa ridotta (nella forma richiesta), rigidezza ridotta e attrito dei
        cuscinetti ridotto. Lo smorzamento strutturale e' :math:`\\beta\\mathbf K_r`.
    k_catena
        Rigidezze equivalenti delle molle della catena, :math:`-K_{r,i,i+1}`.
    """

    asse: str
    modello_massa: ModelloMassa
    forma_massa: FormaMassa
    etichette: tuple[str, ...]
    ruoli: tuple[str, ...]
    T: np.ndarray
    M_r: np.ndarray
    M_r_guyan: np.ndarray
    K_r: np.ndarray
    C_b_r: np.ndarray
    completo: ModelloAsse

    @property
    def n(self) -> int:
        return len(self.etichette)

    @property
    def k_catena(self) -> np.ndarray:
        return np.array([-self.K_r[i, i + 1] for i in range(self.n - 1)])

    def indice_ruolo(self, ruolo: str) -> int:
        return self.ruoli.index(ruolo)

    @property
    def i_encoder(self) -> int:
        return self.completo.masters.index(self.completo.i_encoder)

    @property
    def i_motore(self) -> int:
        return self.completo.masters.index(self.completo.i_motore)

    @property
    def i_estremo(self) -> int:
        return self.completo.masters.index(self.completo.i_estremo)

    @property
    def i_P(self) -> int | None:
        if self.completo.i_P is None:
            return None
        return self.completo.masters.index(self.completo.i_P)

    def vettore_forze(self) -> np.ndarray:
        e = np.zeros(self.n)
        e[self.i_motore] = 1.0
        return e


def _ruoli(mod: ModelloAsse) -> tuple[str, ...]:
    """Ruoli dei master; se i master non sono quelli di riferimento, etichette."""
    attesi = RUOLI[mod.nome]
    if len(mod.masters) != len(attesi):
        return tuple(mod.etichetta(i) for i in mod.masters)
    if mod.nome == "azimut":
        ordine = (mod.i_encoder, mod.i_motore, mod.i_estremo)
    else:
        ordine = (mod.i_encoder, mod.i_P, mod.i_motore, mod.i_estremo)  # type: ignore[assignment]
    if tuple(mod.masters) != tuple(ordine):
        return tuple(mod.etichetta(i) for i in mod.masters)
    return attesi


def massa_lumped(M_r: np.ndarray) -> np.ndarray:
    """``guyan_lumped``: somma per righe sulla diagonale (conserva inerzia e moto rigido)."""
    return np.diag(M_r.sum(axis=1))


def riduci(mod: ModelloAsse, modello_massa: ModelloMassa, forma_massa: FormaMassa = "guyan") -> ModelloRidotto:
    """Condensazione di Guyan dell'asse sui suoi master (sezione 8.1-8.2)."""
    T = trasformazione_guyan(mod.K, mod.masters)
    M_g = T.T @ mod.M[modello_massa] @ T
    M_g = 0.5 * (M_g + M_g.T)
    K_r = T.T @ mod.K @ T
    K_r = 0.5 * (K_r + K_r.T)
    C_b_r = T.T @ mod.C_b @ T
    M_r = M_g if forma_massa == "guyan" else massa_lumped(M_g)
    return ModelloRidotto(
        asse=mod.nome, modello_massa=modello_massa, forma_massa=forma_massa,
        etichette=mod.etichette_master, ruoli=_ruoli(mod), T=T, M_r=M_r, M_r_guyan=M_g,
        K_r=K_r, C_b_r=C_b_r, completo=mod,
    )


def matrici_ridotte_isolate(
    rid: ModelloRidotto, *, J_aggiunta: float = 0.0, k_aggiunta: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """Catena ridotta isolata, con telescopio (polare) o load completo (azimut).

    Le aggiunte stanno su un master, quindi coincidono con quelle del modello
    completo isolato (:func:`gimbal_mdof.fem.matrici_isolate`).
    """
    M, K = rid.M_r.copy(), rid.K_r.copy()
    i = rid.i_P if rid.asse == "polare" else rid.i_estremo
    assert i is not None
    M[i, i] += J_aggiunta
    K[i, i] += k_aggiunta
    return M, K


# ---------------------------------------------------------------------------
# Forma chiusa per una catena (target di regressione, sezione 8.1)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FormaChiusaCatena:
    """Rigidezze in serie fra master consecutivi e pesi degli slave."""

    k_serie: np.ndarray
    pesi: dict[int, tuple[int, int, float]]
    """Per ogni slave ``j``: ``(a, b, w_j)`` con
    :math:`\\varphi_j = (1-w_j)\\varphi_a + w_j\\varphi_b` (indici dei master
    nella catena); gli slave esterni seguono rigidamente il master adiacente."""

    def K_r(self) -> np.ndarray:
        n = len(self.k_serie) + 1
        K = np.zeros((n, n))
        for i, k in enumerate(self.k_serie):
            K[np.ix_([i, i + 1], [i, i + 1])] += k * np.array([[1.0, -1.0], [-1.0, 1.0]])
        return K

    def T(self, n: int, masters: Sequence[int]) -> np.ndarray:
        T = np.zeros((n, len(masters)))
        for c, i in enumerate(masters):
            T[i, c] = 1.0
        for j, (a, b, w) in self.pesi.items():
            T[j, a] += 1.0 - w
            T[j, b] += w
        return T


def forma_chiusa_catena(k_elementi: Sequence[float], masters: Sequence[int]) -> FormaChiusaCatena:
    """Condensazione di Guyan di una catena in forma chiusa.

    .. math::
        \\frac{1}{k_{ab}} = \\sum_{i=a}^{b-1}\\frac{1}{k_i}, \\qquad
        w_j = \\frac{\\sum_{i=a}^{j-1} 1/k_i}{\\sum_{i=a}^{b-1} 1/k_i}

    ``k_elementi[i]`` collega il nodo ``i`` al nodo ``i+1`` (indici interni).
    """
    n = len(k_elementi) + 1
    m = sorted(masters)
    c = np.concatenate([[0.0], np.cumsum(1.0 / np.asarray(k_elementi, dtype=float))])
    k_serie = np.array([1.0 / (c[b] - c[a]) for a, b in zip(m[:-1], m[1:])])
    pesi: dict[int, tuple[int, int, float]] = {}
    for j in range(n):
        if j in m:
            continue
        if j < m[0]:
            pesi[j] = (0, 0, 0.0)
        elif j > m[-1]:
            pesi[j] = (len(m) - 1, len(m) - 1, 0.0)
        else:
            q = int(np.searchsorted(m, j)) - 1
            a, b = m[q], m[q + 1]
            pesi[j] = (q, q + 1, float((c[j] - c[a]) / (c[b] - c[a])))
    return FormaChiusaCatena(k_serie=k_serie, pesi=pesi)


# ---------------------------------------------------------------------------
# Craig-Bampton (solo verifica)
# ---------------------------------------------------------------------------


def riduci_craig_bampton(
    M: np.ndarray, K: np.ndarray, masters: Sequence[int], n_modi: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Craig-Bampton con ``n_modi`` modi a interfaccia fissa.

    .. math:: \\mathbf T_{CB} = \\left[\\mathbf T_{Guyan}\\;\\middle|\\;
              \\begin{smallmatrix}\\mathbf 0\\\\ \\boldsymbol\\Phi_{s}\\end{smallmatrix}\\right]

    Con ``n_modi = 0`` coincide con Guyan (sezione 8.3).

    Returns
    -------
    (T_CB, M_CB, K_CB)
    """
    T_g = trasformazione_guyan(K, masters)
    n = K.shape[0]
    s = [i for i in range(n) if i not in masters]
    n_modi = min(n_modi, len(s))
    if n_modi == 0:
        T = T_g
    else:
        _, Phi = sla.eigh(K[np.ix_(s, s)], M[np.ix_(s, s)])
        Psi = np.zeros((n, n_modi))
        Psi[s, :] = Phi[:, :n_modi]
        T = np.hstack([T_g, Psi])
    M_cb = T.T @ M @ T
    K_cb = T.T @ K @ T
    return T, 0.5 * (M_cb + M_cb.T), 0.5 * (K_cb + K_cb.T)
