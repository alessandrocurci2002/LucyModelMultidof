r"""LQG di Gawronski (*Modeling and Control of Antennas and Telescopes*, 2008, cap. 6 e 8) su un asse.

Plant: la maglia di velocita' chiusa, ingresso comando di velocita' :math:`u`, uscita encoder :math:`y` (Fig. 8.1).
Passi, con i numeri delle equazioni del libro:

1. coordinate modali a blocchi reali, con norma :math:`H_\infty` di ogni modo (4.5), e troncamento;
2. posizione come primo stato: :math:`x_p = [y;\ x_f]` (6.1-6.3);
3. discretizzazione ZOH (2.34-2.37) e aumento con l'integrale della posizione (6.4-6.7);
4. :math:`K_c = [k_i\ k_p\ K_f]` con ``control.dlqr`` e :math:`K_e` con ``control.dlqe`` (riquadro di p. 98);
5. legge d'inseguimento :math:`u = k_i\hat e_i + k_p\hat e - K_f\hat x_f` (8.12, 8.15).

Lo stimatore lavora su :math:`[y;\ x_f]`, senza l'integrale: :math:`y_i` non e' osservabile da :math:`y` e la Riccati
dello stimatore sul modello aumentato non ha soluzione. :math:`\hat e_i` e' l'integrale di :math:`\hat e`, calcolato
nel controllore.
"""

from __future__ import annotations

from dataclasses import dataclass

import control as ct
import numpy as np
from scipy.linalg import eig, solve


@dataclass(frozen=True)
class Modo:
    """Un modo reale (una riga di stato) o una coppia complessa (due righe) delle coordinate modali."""

    indici: tuple[int, ...]
    f_hz: float
    zeta: float
    hinf: float
    """Altezza del picco, :math:`|c\\,b|/(2\\zeta\\omega)` (4.5); per un polo reale, :math:`|c\\,b|/|\\lambda|`."""


def forma_modale(A: np.ndarray, B: np.ndarray, C: np.ndarray, tol_zero: float = 1e-7) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Modo]]:
    """Coordinate modali a blocchi reali: ogni modo ha :math:`\\|B_i\\| = \\|C_i\\|` (la scala rende confrontabili i
    pesi). Il modo rigido (polo nell'origine) e' il primo. ``B`` a una colonna, ``C`` a una riga."""
    lam, V = eig(A)
    W = solve(V, np.eye(len(A), dtype=complex))
    sc = np.abs(lam).max()
    fatti = np.zeros(len(lam), bool)
    colonne, modi = [], []
    k = 0
    for i in np.argsort(np.abs(lam)):
        if fatti[i]:
            continue
        fatti[i] = True
        cb, bb = (C @ V[:, i])[0], (W[i, :] @ B)[0]
        res = abs(cb * bb)
        scala = np.sqrt(abs(bb / cb)) if res > 0 else 1.0
        if abs(lam[i].imag) <= 1e-9 * max(abs(lam[i]), 1.0):
            v = V[:, i] / np.exp(1j * np.angle(V[np.argmax(np.abs(V[:, i])), i]))
            colonne.append(scala * v.real[:, None])
            rigido = abs(lam[i]) < tol_zero * sc
            modi.append(Modo((k,), 0.0 if rigido else abs(lam[i]) / (2 * np.pi), 1.0, np.inf if rigido else float(res / abs(lam[i]))))
            k += 1
        else:
            j = int(np.argmin(np.abs(lam - lam[i].conj()) + fatti * 1e30))
            fatti[j] = True
            om, z = abs(lam[i]), -lam[i].real / abs(lam[i])
            colonne += [scala * V[:, i].real[:, None], scala * V[:, i].imag[:, None]]
            modi.append(Modo((k, k + 1), om / (2 * np.pi), z, float(res / max(z * om, 1e-30))))
            k += 2
    T = np.hstack(colonne)
    Ti = solve(T, np.eye(len(A)))
    return Ti @ A @ T, Ti @ B, C @ T, modi


def tronca(Am, Bm, Cm, modi: list[Modo], tieni) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Modo]]:
    """Coordinate modali troncate ai modi ``tieni`` (predicato su ``Modo``); il modo rigido si tiene sempre (4.11)."""
    scelti = [m for m in modi if not np.isfinite(m.hinf) or tieni(m)]
    idx = [i for m in scelti for i in m.indici]
    ridotti, k = [], 0
    for m in scelti:
        ridotti.append(Modo(tuple(range(k, k + len(m.indici))), m.f_hz, m.zeta, m.hinf))
        k += len(m.indici)
    return Am[np.ix_(idx, idx)], Bm[idx], Cm[:, idx], ridotti


def con_cuscinetti(A, Bu, Bw, C, v_x, v_v, k_0: float, k_1: float, a: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Plant con la linearizzazione armonica dei cuscinetti, :math:`K_b(s) = k_0 + k_1 s/(s + a)`: una molla in parallelo a
    un elemento di Maxwell, con un solo stato in piu' :math:`z`: :math:`\\dot z = \\dot\\delta - a z`, :math:`T_b = k_0\\delta + k_1 z`.
    ``Bw`` e' la colonna di :math:`T_b`; ``v_x``, ``v_v`` danno :math:`\\delta` e :math:`\\dot\\delta` dallo stato."""
    n = len(A)
    A2 = np.zeros((n + 1, n + 1))
    A2[:n, :n] = A + k_0 * np.outer(Bw, v_x)
    A2[:n, n] = k_1 * Bw
    A2[n, :n] = v_v
    A2[n, n] = -a
    return A2, np.r_[Bu, 0.0], np.r_[C, 0.0]


def posizione_primo_stato(Am, Bm, Cm) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """:math:`x_p = [y;\\ x_f] = P x_m`, :math:`P = [C_m;\\ 0\\ I]` (6.1-6.3); richiede il modo rigido per primo."""
    n = len(Am)
    P = np.vstack([Cm, np.eye(n)[1:]])
    Pi = solve(P, np.eye(n))
    return P @ Am @ Pi, P @ Bm, np.eye(n)[:1]


@dataclass
class ControlloreLQG:
    """Guadagni e matrici discrete del controllore; ``passo`` e' il periodo del master."""

    Ad: np.ndarray
    Bd: np.ndarray
    k_i: float
    k_p: float
    K_f: np.ndarray
    K_e: np.ndarray
    passo: float
    poli_regolatore: np.ndarray
    poli_stimatore: np.ndarray

    def stato_iniziale(self) -> dict:
        return {"x": np.zeros(len(self.Ad)), "e_i": 0.0}

    def uscita(self, st: dict, r: float) -> float:
        """Legge (8.15) con :math:`\\hat e = r - \\hat y`; aggiorna :math:`\\hat e_i`."""
        e = r - st["x"][0]
        u = self.k_i * st["e_i"] + self.k_p * e - self.K_f @ st["x"][1:]
        st["e_i"] += self.passo * e
        return float(u)

    def stima(self, st: dict, y: float, u_applicato: float) -> None:
        """Stimatore (8.13) discreto: ``u_applicato`` e' il comando che il drive ha ricevuto in questo periodo."""
        st["x"] = self.Ad @ st["x"] + self.Bd[:, 0] * u_applicato + self.K_e[:, 0] * (y - st["x"][0])

    def funzione_trasferimento(self, z: np.ndarray, ritardo: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """:math:`u = F_r(z) r - F_y(z) y`, con il ritardo di ``ritardo`` periodi fra il controllore e il drive,
        di cui lo stimatore tiene conto."""
        n = len(self.Ad)
        Fy, Fr = np.empty(len(z), complex), np.empty(len(z), complex)
        Kx = np.r_[self.k_p, self.K_f]
        for i, zi in enumerate(z):
            q = zi ** -ritardo
            Ii = self.passo / (zi - 1.0)
            # [zI - A + Ke e1' + B q (kp e1' + Kf) + B q ki Ii e1'] xh = Ke y + B q (kp + ki Ii) r
            D = np.eye(n) * zi - self.Ad + np.outer(self.K_e[:, 0], np.eye(n)[0]) + q * np.outer(self.Bd[:, 0], Kx) \
                + q * self.k_i * Ii * np.outer(self.Bd[:, 0], np.eye(n)[0])
            xy = np.linalg.solve(D, self.K_e[:, 0])
            xr = np.linalg.solve(D, q * (self.k_p + self.k_i * Ii) * self.Bd[:, 0])
            Fy[i] = q * (Kx @ xy + self.k_i * Ii * xy[0])
            Fr[i] = q * (self.k_p + self.k_i * Ii - Kx @ xr - self.k_i * Ii * xr[0])
        return Fr, Fy


class ControlloreLQGPIMisurato(ControlloreLQG):
    """Gawronski 2008, Fig. 8.4b: PI sull'errore misurato e = r - y, retroazione degli stati flessibili stimati; stessi guadagni e
    stesso stimatore di :class:`ControlloreLQG`. La simulazione passa y(k) con :meth:`misura` prima di :meth:`uscita`."""

    def stato_iniziale(self) -> dict:
        return {"x": np.zeros(len(self.Ad)), "e_i": 0.0, "y": 0.0}

    def misura(self, st: dict, y: float) -> None:
        st["y"] = float(y)

    def uscita(self, st: dict, r: float) -> float:
        e = r - st["y"]
        u = self.k_i * st["e_i"] + self.k_p * e - self.K_f @ st["x"][1:]
        st["e_i"] += self.passo * e
        return float(u)

    def funzione_trasferimento(self, z: np.ndarray, ritardo: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """:math:`u = F_r(z) r - F_y(z) y`, con g = k_p + k_i T/(z-1) e q = z^-ritardo:
        :math:`[zI - A + K_e e_1^T + B q K_f'] \\hat x = K_e y + B q g (r - y)`, :math:`u = g (r - y) - K_f' \\hat x`."""
        n = len(self.Ad)
        Fy, Fr = np.empty(len(z), complex), np.empty(len(z), complex)
        Kf = np.r_[0.0, self.K_f]
        for i, zi in enumerate(z):
            q = zi ** -ritardo
            gi = self.k_p + self.k_i * self.passo / (zi - 1.0)
            D = np.eye(n) * zi - self.Ad + np.outer(self.K_e[:, 0], np.eye(n)[0]) + q * np.outer(self.Bd[:, 0], Kf)
            xy = np.linalg.solve(D, self.K_e[:, 0] - q * gi * self.Bd[:, 0])
            xr = np.linalg.solve(D, q * gi * self.Bd[:, 0])
            Fy[i] = q * (gi + Kf @ xy)
            Fr[i] = q * (gi - Kf @ xr)
        return Fr, Fy


def pi_misurato(c: ControlloreLQG) -> ControlloreLQGPIMisurato:
    """Lo stesso controllore con il PI sull'errore misurato (Fig. 8.4b)."""
    from dataclasses import fields
    return ControlloreLQGPIMisurato(**{f.name: getattr(c, f.name) for f in fields(c)})


def sintesi(Ap, Bp, Cp, passo: float, q_i: float, q_p: float, q_f, rho_e: float = 1.0, v_f=None, v_continua: bool = False) -> ControlloreLQG:
    """Pesi :math:`Q = \\mathrm{diag}(q_i, q_p, q_f)`, :math:`R = 1` (8.23-8.25); :math:`V = \\rho_e\\,\\mathrm{diag}(q_p,
    q_f)`, :math:`W = 1` (8.26 per :math:`\\rho_e = 1`). ``q_f`` per stato di :math:`x_f`. ``v_f``: pesi di :math:`x_f` in
    :math:`V` al posto di :math:`\\rho_e q_f`, per tarare :math:`Q` e :math:`V` separatamente (p. 119).
    ``v_continua``: :math:`V` e' l'intensita' della FARE continua (8.8). La discretizzazione da' :math:`V_d = V T`, :math:`W_d = W/T`, cioe' a
    :math:`W_d = 1` il peso :math:`V T^2`: solo cosi' i poli dello stimatore sono quelli della (8.8). Per ``dlqr`` :math:`Q` e :math:`R` scalano
    insieme e questo non serve. Con ``False`` :math:`V` entra per campione (riquadro di p. 98)."""
    d = ct.c2d(ct.ss(Ap, Bp, Cp, 0.0), passo, "zoh")
    Ad, Bd, Cd = np.asarray(d.A), np.asarray(d.B), np.asarray(d.C)
    n = len(Ad)
    Ao = np.block([[np.ones((1, 1)), passo * Cd], [np.zeros((n, 1)), Ad]])
    Bo = np.vstack([np.zeros((1, 1)), Bd])
    q_f = np.asarray(q_f, float)
    Kc, _, pc = ct.dlqr(Ao, Bo, np.diag(np.r_[q_i, q_p, q_f]), np.eye(1))
    V = rho_e * np.diag(np.r_[q_p, q_f]) if v_f is None else np.diag(np.r_[rho_e * q_p, np.asarray(v_f, float)])
    Ke, _, pe = ct.dlqe(Ad, np.eye(n), Cd, passo ** 2 * V if v_continua else V, np.eye(1))
    Kc = np.asarray(Kc)
    return ControlloreLQG(Ad, Bd, float(Kc[0, 0]), float(Kc[0, 1]), Kc[0, 2:], np.asarray(Ke), passo,
                          np.asarray(pc), np.asarray(pe))
