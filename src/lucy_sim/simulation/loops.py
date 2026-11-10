"""Time-domain simulators: one axis with the LQG in the master, and the two coupled axes (the pointer of the thesis).

Integration as :meth:`lucy_sim.plant.loop.ClosedLoopAxis.run`: exact ZOH of the linear part at the simulation step
(0.1 ms), Dahl state updated exactly at constant velocity over the step. The master runs every ``dt_master``: it reads
the encoder, computes the command, which reaches the drive ``delay`` periods later (EtherCAT), and updates the
estimator with the command the drive received. No current limit in the simulation (the current is an output).
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
from scipy.linalg import expm

from ..plant.axis import Axis
from .coupling import LagrangeCoupling
from .disturbances import Disturbance, NoisyMeasurement, disturbance_torque

THETA_STEP = np.deg2rad(0.1)
"""Elevation step at which the discrete azimuth model with J(theta) is recomputed."""


def run_lqg(a: Axis, c, p_cmd: np.ndarray, I0: float, dt_master: float, delay: int = 1, torque: np.ndarray | None = None,
            los: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """One axis with the drive in velocity mode and the controller ``c`` in the master; reference ``p_cmd`` sampled at
    the simulation step, relative to the start. Same outputs as :meth:`ClosedLoopAxis.run`, plus the command ``u``."""
    an, fr = a.velocity_loop, a.friction
    s0, n = fr.sigma0, an.n
    Tb0 = float(np.clip(an.m.Kt * I0, -fr.Tc, fr.Tc))
    z = Tb0 / s0
    s = a.cascade.equilibrium(Tb0)
    p0 = an.out_position @ s
    vv = an.bearing_rate_row
    sub = int(round(dt_master / an.dt))
    st = c.initial_state()
    queue = [0.0] * delay
    P, I, IC, VP, VI, U, L = (np.empty(len(p_cmd)) for _ in range(7))
    Ad, Bd, rp, ri = an.Ad, an.Bd, an.out_position, an.out_current
    row = None if los is None else np.r_[los[:n], np.zeros(an.N - n)]
    l0 = row @ s if row is not None else 0.0
    u_app = y_s = 0.0
    measure = hasattr(c, "measure")
    for k, r in enumerate(p_cmd):
        y = rp @ s - p0
        if k % sub == 0:
            y_s = y
            if measure:
                c.measure(st, y)
            queue.append(c.output(st, float(r)))
            u_app = queue.pop(0)
        inp = np.array([u_app, s0 * z + (torque[k] if torque is not None else 0.0)])
        P[k], I[k], U[k] = y, ri @ s, u_app
        L[k] = row @ s - l0 if row is not None else y
        IC[k], VP[k], VI[k] = an.row_i_cmd @ s + an.in_i_cmd @ inp, an.row_v_phase @ s + an.in_v_phase @ inp, an.row_i_int @ s
        s_new = Ad @ s + Bd @ inp
        vd = 0.5 * (vv @ s + vv @ s_new)
        al = s0 * abs(vd) / fr.g(vd)
        if al > 0:
            e = np.exp(-al * an.dt)
            z = z * e + vd / al * (1 - e)
        s = s_new
        if k % sub == sub - 1:
            c.update(st, y_s, u_app)
    return {"p": P, "i": I, "i_cmd": IC, "v_phase": VP, "i_int": VI, "u": U, "los": L}


def simulate_two_axes(axes: dict[str, Axis], controllers: dict, reference: dict[str, np.ndarray], t: np.ndarray,
                      disturbance: Disturbance, seed: int, coupling: LagrangeCoupling | None, theta0: float,
                      dt_master: float, delay: int = 1, inertia: bool = True, with_H: bool = True) -> dict[str, dict[str, np.ndarray]]:
    """Co-simulation of the two axes (keys ``azimuth``, ``polar``), each in cascade (controller ``None``) or with its
    controller in the master, exchanging at every step the coupling torques of the telescope computed at the start
    of the step from the rates of the telescope nodes; the torques enter those nodes, as in the Lagrange equations.

    ``reference``: absolute azimuth and elevation of the line of sight [rad] at the times ``t`` (simulation step);
    each loop follows it relative to its first value, ``theta0`` is the elevation at the start. ``seed``: seed of the
    disturbance (``seed + 1`` on the polar axis). ``inertia``: the inertia of the telescope node about the azimuth
    follows J(theta) (only the row of its acceleration changes; the discrete model is recomputed every
    :data:`THETA_STEP`; the modal damping stays the one at theta = 0). ``with_H``: the inertial terms -H th'' and
    -H psi'' are added with the acceleration of the other node over the previous step (H^2 / (J I_yy) ~ 1e-12).
    ``coupling = None``: independent axes.

    Outputs per axis: encoder ``p``, current ``i``, line of sight ``los`` (relative to the start) and ``coupling_torque``."""
    names = ("azimuth", "polar")
    n_t = len(t)
    D = {}
    for i, q in enumerate(names):
        a, c = axes[q], controllers[q]
        an = a.cascade if c is None else a.velocity_loop
        Mi = np.linalg.inv(a.model.M)
        bc = np.zeros(an.N)
        bc[an.n:2 * an.n] = Mi @ a.los
        B3 = np.c_[an.B, bc]
        E = expm(np.block([[an.A, B3], [np.zeros((3, an.N + 3))]]) * an.dt)
        r = reference[q] - reference[q][0]
        sigma = disturbance.noise
        if c is None and sigma:
            noise = np.random.default_rng(seed + i).standard_normal(int(n_t * an.dt / dt_master) + 2) * sigma
            r = r - noise[np.floor(t / dt_master).astype(int)]
        x = None if c is None else (NoisyMeasurement(c, sigma, np.random.default_rng(seed + i)) if sigma else c)
        # model error: Coulomb level of the plant scaled, sigma0 and controllers unchanged
        fr = dataclasses.replace(a.friction, Tc=a.friction.Tc * disturbance.friction_factor) if disturbance.friction_factor != 1.0 else a.friction
        D[q] = {"an": an, "Ad": E[:an.N, :an.N], "Bd": E[:an.N, an.N:], "s": a.cascade.equilibrium(disturbance.mean), "z": 0.0,
                "s0": a.friction.sigma0, "fr": fr, "r": r, "dist": disturbance_torque(disturbance, t, seed + i, an.dt), "c": x,
                "st": x.initial_state() if x else None, "queue": [0.0] * delay, "u": 0.0, "ys": 0.0,
                "measure": x is not None and hasattr(x, "measure"), "vv": an.bearing_rate_row, "rp": an.out_position,
                "ri": an.out_current, "rl": np.r_[a.los, np.zeros(an.N - an.n)], "rv": np.r_[np.zeros(an.n), a.los, np.zeros(an.N - 2 * an.n)],
                "P": np.empty(n_t), "I": np.empty(n_t), "L": np.empty(n_t), "T": np.empty(n_t)}
    an0 = D["azimuth"]["an"]
    sub = int(round(dt_master / an0.dt))
    Da, De = D["azimuth"], D["polar"]
    a_az = axes["azimuth"]
    i_v = an0.n + int(np.argmax(a_az.los))
    J0 = float(a_az.model.M[i_v - an0.n, i_v - an0.n])
    bc0 = np.zeros(an0.N)
    bc0[an0.n:2 * an0.n] = np.linalg.inv(a_az.model.M) @ a_az.los
    cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    def azimuth_model(th: float) -> tuple[np.ndarray, np.ndarray]:
        key = int(round(th / THETA_STEP))
        if key not in cache:
            f = J0 / (J0 + coupling.delta_J(key * THETA_STEP))
            A_ = an0.A.copy()
            B_ = np.c_[an0.B, bc0]
            A_[i_v, :] *= f
            B_[i_v, :] *= f
            E_ = expm(np.block([[A_, B_], [np.zeros((3, an0.N + 3))]]) * an0.dt)
            cache[key] = (E_[:an0.N, :an0.N], E_[:an0.N, an0.N:])
        return cache[key]

    v_a0 = v_e0 = 0.0
    for k in range(n_t):
        if coupling is not None:
            th_k = theta0 + float(De["rl"] @ De["s"])
            v_a, v_e = float(Da["rv"] @ Da["s"]), float(De["rv"] @ De["s"])
            ta, te = coupling.torques(th_k, v_a, v_e)
            if with_H:
                h = coupling.H(th_k)
                ta -= h * (v_e - v_e0) / an0.dt
                te -= h * (v_a - v_a0) / an0.dt
                v_a0, v_e0 = v_a, v_e
            if inertia:
                Da["Ad"], Da["Bd"] = azimuth_model(th_k)
        else:
            ta, te = 0.0, 0.0
        Da["T"][k], De["T"][k] = ta, te
        for X, tq in ((Da, ta), (De, te)):
            s = X["s"]
            y = float(X["rp"] @ s)
            if X["c"] is None:
                u0 = X["r"][k]
            else:
                if k % sub == 0:
                    X["ys"] = y
                    if X["measure"]:
                        X["c"].measure(X["st"], y)
                    X["queue"].append(X["c"].output(X["st"], float(X["r"][k])))
                    X["u"] = X["queue"].pop(0)
                u0 = X["u"]
            inp = np.array([u0, X["s0"] * X["z"] + X["dist"][k], tq])
            X["P"][k], X["I"][k], X["L"][k] = y, float(X["ri"] @ s), float(X["rl"] @ s)
            s_new = X["Ad"] @ s + X["Bd"] @ inp
            vd = 0.5 * float(X["vv"] @ s + X["vv"] @ s_new)
            al = X["s0"] * abs(vd) / X["fr"].g(vd)
            if al > 0:
                e = math.exp(-al * an0.dt)
                X["z"] = X["z"] * e + vd / al * (1 - e)
            X["s"] = s_new
            if X["c"] is not None and k % sub == sub - 1:
                X["c"].update(X["st"], X["ys"], X["u"])
    return {q: {"p": D[q]["P"], "i": D[q]["I"], "los": D[q]["L"], "coupling_torque": D[q]["T"]} for q in names}
