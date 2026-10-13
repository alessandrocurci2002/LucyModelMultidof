"""Simulazioni non lineari del sistema accoppiato ridotto (sezione 11, punto 5).

Integrazione con :func:`scipy.integrate.solve_ivp`, metodo ``DOP853``,
``rtol = 1e-10``, ``atol = 1e-12``. Il doppietto di coppia (riusato da
:func:`gimbal.simulate.doppietto`) e' discontinuo: l'integrazione e' spezzata
nei suoi istanti di commutazione, cosi' l'integratore non attraversa mai un
salto della forzante.

In ogni scenario con punto di lavoro le coppie nominali :math:`\\boldsymbol\\tau_0`
sono applicate in avanti e lo stato iniziale e' l'equilibrio stazionario della
sezione 10.1: la perturbazione e' il solo doppietto.

* **S1**: doppietto su :math:`\\tau_\\psi` a :math:`\\Omega_0 = 0`, confronto con il
  modello lineare (errore relativo < 1%).
* **S2**: inseguimento a :math:`\\Omega_{max}` con doppietto su :math:`\\tau_\\psi`;
  risposta incrociata su :math:`\\theta_e` e :math:`\\theta_p`.
* **S3**: :math:`\\mathbf C = 0`, :math:`\\boldsymbol\\tau = 0`: conservazione di
  :math:`E = T + V` e di :math:`p_\\psi`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Final, Sequence

import numpy as np
from gimbal.simulate import doppietto
from scipy.integrate import solve_ivp

from . import linearize as lz
from .coupled import ModelloAccoppiato

RTOL: Final[float] = 1e-10
ATOL: Final[float] = 1e-12
METODO: Final[str] = "DOP853"


def integra(campo: Callable[[float, np.ndarray], np.ndarray], x0: np.ndarray,
            t_eval: np.ndarray, commutazioni: Sequence[float] = ()) -> np.ndarray:
    """Integra a tratti fra le commutazioni della forzante; restituisce ``x(t_eval)``."""
    bordi = sorted({float(t_eval[0]), float(t_eval[-1]),
                    *[c for c in commutazioni if t_eval[0] < c < t_eval[-1]]})
    X = np.empty((t_eval.size, x0.size))
    x = np.asarray(x0, dtype=float)
    for a, b in zip(bordi[:-1], bordi[1:]):
        maschera = (t_eval >= a) & ((t_eval < b) | ((b == bordi[-1]) & (t_eval <= b)))
        sol = solve_ivp(campo, (a, b), x, method=METODO, rtol=RTOL, atol=ATOL,
                        t_eval=t_eval[maschera] if maschera.any() else None, dense_output=False)
        if not sol.success:
            raise RuntimeError(f"integrazione non riuscita su [{a}, {b}]: {sol.message}")
        if maschera.any():
            X[maschera] = sol.y.T
        sol_fine = solve_ivp(campo, (a, b), x, method=METODO, rtol=RTOL, atol=ATOL) \
            if not maschera.any() or not np.isclose(sol.t[-1], b) else sol
        x = sol_fine.y[:, -1]
    return X


@dataclass(frozen=True)
class RisultatoScenario:
    """Storia temporale di uno scenario e confronto con il lineare."""

    nome: str
    descrizione: str
    t: np.ndarray
    y_nl: np.ndarray
    """Uscite :data:`gimbal_mdof.linearize.USCITE` del non lineare, come scarto dall'equilibrio."""
    y_lin: np.ndarray | None
    theta0: float
    Omega0: float
    ampiezza: float
    errori: dict[str, float] = field(default_factory=dict)
    energia: np.ndarray | None = None
    momento: np.ndarray | None = None
    derive: dict[str, float] = field(default_factory=dict)


def _uscite(ml: lz.ModelloLineare, dq: np.ndarray) -> np.ndarray:
    return dq @ ml.Cq.T


def scenario_doppietto(acc: ModelloAccoppiato, theta0: float, Omega0: float, *, ampiezza: float,
                       durata: float, t_fine: float, n_punti: int = 2001, nome: str = "S1",
                       descrizione: str = "") -> RisultatoScenario:
    """Doppietto su :math:`\\tau_\\psi` intorno al punto di lavoro, non lineare e lineare."""
    ml = lz.linearizza(acc, theta0, Omega0)
    q0, tau0 = ml.punto.q0, ml.punto.tau0
    n = acc.n
    t = np.linspace(0.0, t_fine, n_punti)
    t_ini = 0.1 * durata
    comm = (t_ini, t_ini + durata, t_ini + 2 * durata)

    def tau(tt: float) -> np.ndarray:
        return tau0 + np.array([doppietto(tt, t_inizio=t_ini, durata=durata, ampiezza=ampiezza), 0.0])

    def campo_nl(tt: float, x: np.ndarray) -> np.ndarray:
        return np.concatenate([x[n:], acc.accelerazioni(x[:n], x[n:], tau(tt))])

    dq0 = np.zeros(n)
    dq0[: acc.struttura.n_psi] = Omega0
    X = integra(campo_nl, np.concatenate([q0, dq0]), t, comm)
    deriva = np.zeros((t.size, n))
    deriva[:, : acc.struttura.n_psi] = Omega0 * t[:, None]
    dq_nl = X[:, :n] - q0[None] - deriva

    def campo_lin(tt: float, x: np.ndarray) -> np.ndarray:
        du = np.array([doppietto(tt, t_inizio=t_ini, durata=durata, ampiezza=ampiezza), 0.0])
        return ml.A @ x + ml.Bu @ du

    X_lin = integra(campo_lin, np.zeros(2 * n), t, comm)
    y_nl = _uscite(ml, dq_nl)
    y_lin = _uscite(ml, X_lin[:, :n])
    errori = {}
    for i, u in enumerate(lz.USCITE):
        picco = float(np.max(np.abs(y_lin[:, i])))
        errori[u] = float(np.max(np.abs(y_nl[:, i] - y_lin[:, i])) / picco) if picco > 0 else np.nan
    return RisultatoScenario(nome=nome, descrizione=descrizione, t=t, y_nl=y_nl, y_lin=y_lin,
                             theta0=theta0, Omega0=Omega0, ampiezza=ampiezza, errori=errori)


def S1(acc: ModelloAccoppiato, theta0: float, *, ampiezza: float = 0.05, durata: float = 0.02,
       t_fine: float = 0.3, n_punti: int = 1501) -> RisultatoScenario:
    return scenario_doppietto(acc, theta0, 0.0, ampiezza=ampiezza, durata=durata, t_fine=t_fine,
                              n_punti=n_punti, nome="S1",
                              descrizione="doppietto su tau_psi a Omega_0 = 0")


def S2(acc: ModelloAccoppiato, theta0: float, Omega_max: float, *, ampiezza: float = 0.05,
       durata: float = 0.02, t_fine: float = 0.3, n_punti: int = 1501) -> RisultatoScenario:
    return scenario_doppietto(acc, theta0, Omega_max, ampiezza=ampiezza, durata=durata,
                              t_fine=t_fine, n_punti=n_punti, nome="S2",
                              descrizione="inseguimento a Omega_max con doppietto su tau_psi")


def S3(acc: ModelloAccoppiato, theta0: float, Omega_max: float, *, t_fine: float = 0.2,
       n_punti: int = 1001, dtheta0: float = 0.02) -> RisultatoScenario:
    """Conservazione con :math:`\\mathbf C = 0` e :math:`\\boldsymbol\\tau = 0`.

    Stato iniziale di moto rigido: tutte le :math:`\\dot\\psi = \\Omega_{max}`,
    tutte le :math:`\\dot\\theta = ` ``dtheta0``, catene indeformate. I modi
    torsionali vengono eccitati solo dalle forze giroscopiche e gravitazionali,
    come in un moto di puntamento reale.
    """
    cons = acc.senza_dissipazione()
    st = cons.struttura
    n = cons.n
    q0 = np.zeros(n)
    q0[st.n_psi:] = theta0
    dq0 = np.zeros(n)
    dq0[: st.n_psi] = Omega_max
    dq0[st.n_psi:] = dtheta0
    t = np.linspace(0.0, t_fine, n_punti)
    zero = np.zeros(2)

    def campo(tt: float, x: np.ndarray) -> np.ndarray:
        return np.concatenate([x[n:], cons.accelerazioni(x[:n], x[n:], zero)])

    X = integra(campo, np.concatenate([q0, dq0]), t)
    E = np.array([cons.energia(x[:n], x[n:]) for x in X])
    p = np.array([cons.momento_psi(x[:n], x[n:]) for x in X])
    derive = {"energia": float(np.max(np.abs(E - E[0])) / abs(E[0])),
              "momento_psi": float(np.max(np.abs(p - p[0])) / abs(p[0]))}
    y = (X[:, :n] - q0[None]) @ lz.costruisci_uscite(cons).T
    return RisultatoScenario(nome="S3", descrizione="conservazione con C = 0 e tau = 0", t=t,
                             y_nl=y, y_lin=None, theta0=theta0, Omega0=Omega_max, ampiezza=0.0,
                             energia=E, momento=p, derive=derive)
