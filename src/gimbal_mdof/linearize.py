"""Punto di lavoro e linearizzazione del modello accoppiato (sezioni 10.1-10.2).

Punto di lavoro: :math:`\\theta_p = \\theta_0`, tutte le :math:`\\dot\\psi =
\\Omega_0`, tutte le :math:`\\dot\\theta = 0`, :math:`\\ddot{\\mathbf q} = 0`. Le
altre coordinate (torsioni statiche) e le coppie nominali
:math:`\\boldsymbol\\tau_0` si ottengono dall'equilibrio stazionario, con
:math:`\\psi_e = 0` (coordinata ciclica). Con :math:`\\theta_p` fissato il
residuo e' affine nelle coordinate rimanenti (i termini del payload dipendono
solo da :math:`\\theta_p`), quindi l'equilibrio e' un sistema lineare risolto in
forma chiusa.

Le matrici linearizzate sono jacobiani simbolici del residuo, valutati nel
punto di lavoro (:func:`gimbal_mdof.coupled.funzioni_numeriche`):

.. math:: \\mathbf M_0\\delta\\ddot{\\mathbf q} + \\mathbf C_0\\delta\\dot{\\mathbf q}
          + \\mathbf K_0\\delta\\mathbf q = \\mathbf E\\,\\delta\\boldsymbol\\tau.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Sequence

import numpy as np

from .coupled import ModelloAccoppiato

USCITE: Final[tuple[str, ...]] = (
    "psi_e", "theta_e", "psi_m", "theta_m", "psi_L", "theta_p", "theta_r",
    "err_psi", "err_theta",
)
"""Uscite della sezione 10.3; ``err_*`` = encoder meno angolo vero."""

ETICHETTE_USCITE: Final[dict[str, str]] = {
    "psi_e": r"$\psi_e$ (encoder az.)",
    "theta_e": r"$\theta_e$ (encoder pol.)",
    "psi_m": r"$\psi_m$ (motore az.)",
    "theta_m": r"$\theta_m$ (motore pol.)",
    "psi_L": r"$\psi_L$ (forcella)",
    "theta_p": r"$\theta_p$ (telescopio)",
    "theta_r": r"$\theta_r$ (estremita' pol.)",
    "err_psi": r"$\psi_e - \psi_L$",
    "err_theta": r"$\theta_e - \theta_p$",
}


@dataclass(frozen=True)
class PuntoDiLavoro:
    """Equilibrio stazionario (sezione 10.1)."""

    theta0: float
    Omega0: float
    q0: np.ndarray
    tau0: np.ndarray
    """Coppie nominali :math:`[\\tau_{\\psi 0}, \\tau_{\\theta 0}]`."""
    residuo: float
    """Norma del residuo dell'equilibrio, per controllo."""

    def torsioni_polari(self, n_psi: int) -> np.ndarray:
        """Coordinate polari all'equilibrio, :math:`[\\theta_e, \\theta_p, \\theta_m, \\theta_r]`."""
        return self.q0[n_psi:]


def equilibrio(acc: ModelloAccoppiato, theta0: float, Omega0: float) -> PuntoDiLavoro:
    """Risolve l'equilibrio stazionario: torsioni statiche e coppie nominali."""
    st, p = acc.struttura, acc.parametri
    n = st.n
    q_base = np.zeros(n)
    q_base[st.j_theta_p] = theta0
    r0 = np.asarray(acc.funzioni.residuo_statico(q_base, Omega0, p), dtype=float).reshape(n)
    K0 = np.asarray(acc.funzioni.K0(theta0, Omega0, p), dtype=float)
    libere = [j for j in range(n) if j not in (st.j_psi_e, st.j_theta_p)]
    E = acc.matrice_ingressi()
    A = np.hstack([K0[:, libere], -E])
    z = np.linalg.solve(A, -r0)
    q0 = q_base.copy()
    q0[libere] = z[: len(libere)]
    tau0 = z[len(libere):]
    res = np.asarray(acc.funzioni.residuo_statico(q0, Omega0, p), dtype=float).reshape(n) - E @ tau0
    scala = max(float(np.max(np.abs(r0))), float(np.max(np.abs(E @ tau0))), 1e-12)
    if np.max(np.abs(res)) > 1e-8 * scala:
        raise AssertionError(f"equilibrio non soddisfatto: residuo {np.max(np.abs(res)):.3e}")
    return PuntoDiLavoro(theta0=float(theta0), Omega0=float(Omega0), q0=q0, tau0=tau0,
                         residuo=float(np.max(np.abs(res))))


@dataclass(frozen=True)
class ModelloLineare:
    """Modello linearizzato nel punto di lavoro, con spazio degli stati.

    :math:`\\mathbf x = [\\delta\\mathbf q, \\delta\\dot{\\mathbf q}]`,
    :math:`\\mathbf u = \\delta\\boldsymbol\\tau`, uscite :data:`USCITE`.
    """

    theta0: float
    Omega0: float
    M0: np.ndarray
    C0: np.ndarray
    K0: np.ndarray
    E: np.ndarray
    A: np.ndarray
    Bu: np.ndarray
    Cq: np.ndarray
    """Matrice delle uscite sulle coordinate (n_uscite x n)."""
    punto: PuntoDiLavoro
    accoppiato: ModelloAccoppiato

    @property
    def n(self) -> int:
        return self.M0.shape[0]

    @property
    def Cy(self) -> np.ndarray:
        """Uscite sullo stato: :math:`[\\mathbf C_q, \\mathbf 0]`."""
        return np.hstack([self.Cq, np.zeros_like(self.Cq)])

    def riga(self, uscita: str) -> int:
        return USCITE.index(uscita)

    def P(self, s: complex) -> np.ndarray:
        return self.M0 * s**2 + self.C0 * s + self.K0

    def frf(self, omega: np.ndarray, uscite: Sequence[str] = USCITE) -> np.ndarray:
        """:math:`\\mathbf G(j\\omega) = \\mathbf C_q\\mathbf P(j\\omega)^{-1}\\mathbf E`.

        Returns
        -------
        np.ndarray
            Forma ``(len(omega), len(uscite), 2)``.
        """
        from .tf import frf_secondo_ordine

        righe = [self.riga(u) for u in uscite]
        return frf_secondo_ordine(self.M0, self.C0, self.K0, self.E, self.Cq[righe], omega)

    def frf_spazio_stati(self, omega: np.ndarray) -> np.ndarray:
        """FRF da :math:`\\mathbf C(j\\omega\\mathbb 1 - \\mathbf A)^{-1}\\mathbf B`, per controllo."""
        omega = np.atleast_1d(omega)
        I = np.eye(self.A.shape[0])
        return np.stack([self.Cy @ np.linalg.solve(1j * w * I - self.A, self.Bu) for w in omega])

    @property
    def poli(self) -> np.ndarray:
        return np.linalg.eigvals(self.A)


def costruisci_uscite(acc: ModelloAccoppiato) -> np.ndarray:
    st = acc.struttura
    n = st.n
    righe = []
    j_theta_r = st.n_psi + (st.n_theta - 1)
    indici = {"psi_e": st.j_psi_e, "theta_e": st.j_theta_e, "psi_m": st.j_psi_m,
              "theta_m": st.j_theta_m, "psi_L": st.j_psi_L, "theta_p": st.j_theta_p,
              "theta_r": j_theta_r}
    for u in USCITE:
        c = np.zeros(n)
        if u == "err_psi":
            c[st.j_psi_e] += 1.0
            c[st.j_psi_L] -= 1.0
        elif u == "err_theta":
            c[st.j_theta_e] += 1.0
            c[st.j_theta_p] -= 1.0
        else:
            c[indici[u]] = 1.0
        righe.append(c)
    return np.array(righe)


def linearizza(acc: ModelloAccoppiato, theta0: float, Omega0: float) -> ModelloLineare:
    """Matrici linearizzate e spazio degli stati in :math:`(\\theta_0, \\Omega_0)`."""
    p = acc.parametri
    M0 = np.asarray(acc.funzioni.M0(theta0, Omega0, p), dtype=float)
    C0 = np.asarray(acc.funzioni.C0(theta0, Omega0, p), dtype=float)
    K0 = np.asarray(acc.funzioni.K0(theta0, Omega0, p), dtype=float)
    E = acc.matrice_ingressi()
    n = M0.shape[0]
    Minv = np.linalg.inv(M0)
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-Minv @ K0, -Minv @ C0]])
    Bu = np.vstack([np.zeros((n, 2)), Minv @ E])
    return ModelloLineare(theta0=float(theta0), Omega0=float(Omega0), M0=M0, C0=C0, K0=K0,
                          E=E, A=A, Bu=Bu, Cq=costruisci_uscite(acc),
                          punto=equilibrio(acc, theta0, Omega0), accoppiato=acc)
