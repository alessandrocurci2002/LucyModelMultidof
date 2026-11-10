r"""LQG-PI controller of Gawronski (*Modeling and Control of Antennas and Telescopes*, 2008, ch. 6 and 8) for one axis.

Plant: the closed velocity loop of the drive, input velocity command :math:`u`, output encoder :math:`y` (Fig. 8.1). The
LQG replaces the position loop in the master at 1 kHz; current and velocity loops stay in the drive. Steps, with the
equation numbers of the book:

1. real block-modal coordinates, with the :math:`H_\infty` norm of every mode (4.5), and truncation;
2. position as the first state, :math:`x_p = [y;\ x_f]` (6.1-6.3);
3. ZOH discretisation (2.34-2.37) and augmentation with the integral of the position (6.4-6.7);
4. :math:`K_c = [k_i\ k_p\ K_f]` with ``control.dlqr`` and :math:`K_e` with ``control.dlqe`` (box on p. 98).

Two tracking laws with the same gains and the same estimator:

* :class:`LQGPIMeasured` (Fig. 8.4b, the controller of the thesis): :math:`u = k_i e_i + k_p e - K_f\hat x_f` with the
  measured error :math:`e = r - y`;
* :class:`LQGController` (8.12, 8.15): the same law on the estimated error :math:`\hat e = r - \hat y`.

The estimator works on :math:`[y;\ x_f]` without the integral (:math:`y_i` is not observable from :math:`y`); the
integral of the error is computed in the controller.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import control as ct
import numpy as np
from scipy.linalg import eig, expm, solve


@dataclass(frozen=True)
class Mode:
    """A real mode (one state row) or a complex pair (two rows) of the modal coordinates."""

    indices: tuple[int, ...]
    f_hz: float
    zeta: float
    hinf: float
    """Peak height, :math:`|c\\,b|/(2\\zeta\\omega)` (4.5); for a real pole :math:`|c\\,b|/|\\lambda|`."""


def modal_form(A: np.ndarray, B: np.ndarray, C: np.ndarray, tol_zero: float = 1e-7) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Mode]]:
    """Real block-modal coordinates, every mode scaled to :math:`\\|B_i\\| = \\|C_i\\|`; the rigid mode (pole at the
    origin) first. ``B`` with one column, ``C`` with one row."""
    lam, V = eig(A)
    W = solve(V, np.eye(len(A), dtype=complex))
    sc = np.abs(lam).max()
    done = np.zeros(len(lam), bool)
    columns, modes = [], []
    k = 0
    for i in np.argsort(np.abs(lam)):
        if done[i]:
            continue
        done[i] = True
        cb, bb = (C @ V[:, i])[0], (W[i, :] @ B)[0]
        res = abs(cb * bb)
        scale = np.sqrt(abs(bb / cb)) if res > 0 else 1.0
        if abs(lam[i].imag) <= 1e-9 * max(abs(lam[i]), 1.0):
            v = V[:, i] / np.exp(1j * np.angle(V[np.argmax(np.abs(V[:, i])), i]))
            columns.append(scale * v.real[:, None])
            rigid = abs(lam[i]) < tol_zero * sc
            modes.append(Mode((k,), 0.0 if rigid else abs(lam[i]) / (2 * np.pi), 1.0, np.inf if rigid else float(res / abs(lam[i]))))
            k += 1
        else:
            j = int(np.argmin(np.abs(lam - lam[i].conj()) + done * 1e30))
            done[j] = True
            om, z = abs(lam[i]), -lam[i].real / abs(lam[i])
            columns += [scale * V[:, i].real[:, None], scale * V[:, i].imag[:, None]]
            modes.append(Mode((k, k + 1), om / (2 * np.pi), z, float(res / max(z * om, 1e-30))))
            k += 2
    T = np.hstack(columns)
    Ti = solve(T, np.eye(len(A)))
    return Ti @ A @ T, Ti @ B, C @ T, modes


def truncate(Am, Bm, Cm, modes: list[Mode], keep) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Mode]]:
    """Modal coordinates truncated to the modes ``keep`` (predicate on :class:`Mode`); the rigid mode is always kept (4.11)."""
    chosen = [m for m in modes if not np.isfinite(m.hinf) or keep(m)]
    idx = [i for m in chosen for i in m.indices]
    reduced, k = [], 0
    for m in chosen:
        reduced.append(Mode(tuple(range(k, k + len(m.indices))), m.f_hz, m.zeta, m.hinf))
        k += len(m.indices)
    return Am[np.ix_(idx, idx)], Bm[idx], Cm[:, idx], reduced


def position_first(Am, Bm, Cm) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """:math:`x_p = [y;\\ x_f] = P x_m`, :math:`P = [C_m;\\ 0\\ I]` (6.1-6.3); needs the rigid mode first."""
    n = len(Am)
    P = np.vstack([Cm, np.eye(n)[1:]])
    Pi = solve(P, np.eye(n))
    return P @ Am @ Pi, P @ Bm, np.eye(n)[:1]


@dataclass
class LQGController:
    """Gains and discrete matrices of the controller (PI on the estimated error, 8.15); ``dt`` is the master period.

    Interface used by the simulators, once per master period: :meth:`measure` (only :class:`LQGPIMeasured`),
    :meth:`output` (command for the reference ``r``), then, after the drive received the command, :meth:`update`."""

    Ad: np.ndarray
    Bd: np.ndarray
    k_i: float
    k_p: float
    K_f: np.ndarray
    K_e: np.ndarray
    dt: float
    regulator_poles: np.ndarray
    estimator_poles: np.ndarray

    def initial_state(self) -> dict:
        return {"x": np.zeros(len(self.Ad)), "e_i": 0.0}

    def output(self, st: dict, r: float) -> float:
        """Law (8.15) with :math:`\\hat e = r - \\hat y`; updates :math:`\\hat e_i`."""
        e = r - st["x"][0]
        u = self.k_i * st["e_i"] + self.k_p * e - self.K_f @ st["x"][1:]
        st["e_i"] += self.dt * e
        return float(u)

    def update(self, st: dict, y: float, u_applied: float) -> None:
        """Estimator (8.13): ``u_applied`` is the command the drive received in this period."""
        st["x"] = self.Ad @ st["x"] + self.Bd[:, 0] * u_applied + self.K_e[:, 0] * (y - st["x"][0])

    def transfer_functions(self, z: np.ndarray, delay: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """:math:`u = F_r(z) r - F_y(z) y`, with ``delay`` periods between controller and drive."""
        n = len(self.Ad)
        Fy, Fr = np.empty(len(z), complex), np.empty(len(z), complex)
        Kx = np.r_[self.k_p, self.K_f]
        for i, zi in enumerate(z):
            q = zi ** -delay
            Ii = self.dt / (zi - 1.0)
            D = np.eye(n) * zi - self.Ad + np.outer(self.K_e[:, 0], np.eye(n)[0]) + q * np.outer(self.Bd[:, 0], Kx) \
                + q * self.k_i * Ii * np.outer(self.Bd[:, 0], np.eye(n)[0])
            xy = np.linalg.solve(D, self.K_e[:, 0])
            xr = np.linalg.solve(D, q * (self.k_p + self.k_i * Ii) * self.Bd[:, 0])
            Fy[i] = q * (Kx @ xy + self.k_i * Ii * xy[0])
            Fr[i] = q * (self.k_p + self.k_i * Ii - Kx @ xr - self.k_i * Ii * xr[0])
        return Fr, Fy


class LQGPIMeasured(LQGController):
    """Gawronski 2008, Fig. 8.4b: PI on the measured error :math:`e = r - y` and feedback of the estimated flexible
    states; same gains and estimator as :class:`LQGController`. The simulators pass y(k) with :meth:`measure` before
    :meth:`output`."""

    def initial_state(self) -> dict:
        return {"x": np.zeros(len(self.Ad)), "e_i": 0.0, "y": 0.0}

    def measure(self, st: dict, y: float) -> None:
        st["y"] = float(y)

    def output(self, st: dict, r: float) -> float:
        e = r - st["y"]
        u = self.k_i * st["e_i"] + self.k_p * e - self.K_f @ st["x"][1:]
        st["e_i"] += self.dt * e
        return float(u)

    def transfer_functions(self, z: np.ndarray, delay: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """With :math:`g = k_p + k_i T/(z-1)` and :math:`q = z^{-d}`:
        :math:`[zI - A + K_e e_1^T + B q K_f'] \\hat x = K_e y + B q g (r - y)`, :math:`u = g (r - y) - K_f' \\hat x`."""
        n = len(self.Ad)
        Fy, Fr = np.empty(len(z), complex), np.empty(len(z), complex)
        Kf = np.r_[0.0, self.K_f]
        for i, zi in enumerate(z):
            q = zi ** -delay
            gi = self.k_p + self.k_i * self.dt / (zi - 1.0)
            D = np.eye(n) * zi - self.Ad + np.outer(self.K_e[:, 0], np.eye(n)[0]) + q * np.outer(self.Bd[:, 0], Kf)
            xy = np.linalg.solve(D, self.K_e[:, 0] - q * gi * self.Bd[:, 0])
            xr = np.linalg.solve(D, q * gi * self.Bd[:, 0])
            Fy[i] = q * (gi + Kf @ xy)
            Fr[i] = q * (gi - Kf @ xr)
        return Fr, Fy


def with_measured_pi(c: LQGController) -> LQGPIMeasured:
    """The same controller with the PI on the measured error (Fig. 8.4b)."""
    return LQGPIMeasured(**{f.name: getattr(c, f.name) for f in fields(c)})


def synthesize(Ap, Bp, Cp, dt: float, q_i: float, q_p: float, q_f, rho_e: float = 1.0, v_f=None, v_continuous: bool = False,
               measured_pi: bool = True) -> LQGController:
    """Weights :math:`Q = \\mathrm{diag}(q_i, q_p, q_f)`, :math:`R = 1` (8.23-8.25); :math:`V = \\rho_e\\,\\mathrm{diag}(q_p,
    q_f)`, :math:`W = 1` (8.26 for :math:`\\rho_e = 1`). ``q_f`` per state of :math:`x_f`; ``v_f``: weights of :math:`x_f`
    in :math:`V` instead of :math:`\\rho_e q_f`.

    ``v_continuous``: :math:`V` is the intensity of the continuous FARE (8.8); the discretisation gives :math:`V_d = V T`,
    :math:`W_d = W/T`, i.e. the weight :math:`V T^2` at :math:`W_d = 1`, so that the estimator poles are those of (8.8).
    With ``False`` :math:`V` enters per sample (box on p. 98). ``measured_pi``: :class:`LQGPIMeasured` (default, the
    thesis) or :class:`LQGController`."""
    d = ct.c2d(ct.ss(Ap, Bp, Cp, 0.0), dt, "zoh")
    Ad, Bd, Cd = np.asarray(d.A), np.asarray(d.B), np.asarray(d.C)
    n = len(Ad)
    Ao = np.block([[np.ones((1, 1)), dt * Cd], [np.zeros((n, 1)), Ad]])
    Bo = np.vstack([np.zeros((1, 1)), Bd])
    q_f = np.asarray(q_f, float)
    Kc, _, pc = ct.dlqr(Ao, Bo, np.diag(np.r_[q_i, q_p, q_f]), np.eye(1))
    V = rho_e * np.diag(np.r_[q_p, q_f]) if v_f is None else np.diag(np.r_[rho_e * q_p, np.asarray(v_f, float)])
    Ke, _, pe = ct.dlqe(Ad, np.eye(n), Cd, dt ** 2 * V if v_continuous else V, np.eye(1))
    Kc = np.asarray(Kc)
    c = LQGController(Ad, Bd, float(Kc[0, 0]), float(Kc[0, 1]), Kc[0, 2:], np.asarray(Ke), dt, np.asarray(pc), np.asarray(pe))
    return with_measured_pi(c) if measured_pi else c


def closed_loop_poles(c: LQGController, A: np.ndarray, B: np.ndarray, C: np.ndarray, delay: int = 1) -> np.ndarray:
    """Discrete closed-loop poles with any continuous plant ``(A, B, C)``: states [x, x_hat, e_i, u in the queue]."""
    N, n = len(A), len(c.Ad)
    E = expm(np.block([[A, B[:, None]], [np.zeros((1, N + 1))]]) * c.dt)
    tot = N + n + 1 + delay
    M = np.zeros((tot, tot))
    ie, iq = N + n, N + n + 1
    M[:N, :N], M[:N, iq] = E[:N, :N], E[:N, N]
    M[N:N + n, :N] = c.K_e[:, [0]] @ C[None, :]
    M[N:N + n, N:N + n] = c.Ad - c.K_e[:, [0]] @ np.eye(n)[:1]
    M[N:N + n, iq] = c.Bd[:, 0]
    M[ie, ie], M[ie, N] = 1.0, -c.dt
    M[iq, ie], M[iq, N] = c.k_i, -c.k_p
    if isinstance(c, LQGPIMeasured):  # Fig. 8.4b: the PI reads y = C x instead of y_hat
        M[ie, N], M[iq, N] = 0.0, 0.0
        M[ie, :N], M[iq, :N] = -c.dt * C, -c.k_p * C
    M[iq, N + 1:N + n] = -c.K_f
    return np.linalg.eigvals(M)


def worst_pole(z: np.ndarray, dt: float) -> dict[str, float]:
    """Largest modulus among the discrete poles ``z`` (sample time ``dt``) and its frequency."""
    i = int(np.argmax(np.abs(z)))
    return {"abs_z": float(abs(z[i])), "f_Hz": float(abs(np.log(complex(z[i])) / dt) / (2 * np.pi))}
