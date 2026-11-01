r"""Nonlinear mechanical plant: the lumped model of ``gimbal_mdof`` (by default the 8-DOF H6 model).

The plant is built from any configuration file of the ``gimbal_mdof`` project, so the model can be changed by
pointing ``plant.config`` of ``config.yaml`` to another file (the 7-DOF baseline ``config.yaml`` or any hub
version). The equations are the Newton--Euler equations of ``gimbal_mdof.concentrato`` with all the
trigonometric terms of the telescope (``gimbal_mdof.simulatore_concentrato``):

.. math::
    \mathbf M(\theta_p)\,\ddot{\mathbf q} = \mathbf E\,\boldsymbol\tau_m + \mathbf Q_f - \mathbf h(\mathbf q, \dot{\mathbf q}),

with :math:`\boldsymbol\tau_m = [\tau_\psi, \tau_\theta]` the motor torques on the motor nodes and
:math:`\mathbf Q_f` the bearing friction. For the H6 model
:math:`\mathbf q = [\psi_e, \psi_m, \psi_L, \theta_e, \theta_h, \theta_m, \theta_r, \theta_p]`.

Bearing friction. H6 represents it below a few hertz as a linear equivalent :math:`k(1 + j\eta)` on the bearing
node (section ``ground_constraint``), valid at the amplitude of the test and only in the frequency domain. In
the time domain it is replaced by its physical counterpart, a Dahl friction element,

.. math::
    \dot T_f = \sigma_0\Bigl(\dot q - |\dot q|\,\frac{T_f}{T_c}\Bigr),\qquad Q_f = -T_f,

whose describing function at the test amplitude :math:`A` is :math:`k(1 + j\eta)` (:func:`dahl_from_hysteretic`):
stiff with little loss for small motions, Coulomb friction :math:`T_c` in gross sliding. In the integration
:math:`|\dot q|` is replaced by :math:`\sqrt{\dot q^2 + v_\varepsilon^2} - v_\varepsilon`, with :math:`v_\varepsilon`
far below the speeds of interest (default :math:`10^{-6}` rad/s, 1/70 of the sidereal rate): the kink of
:math:`|\dot q|` at rest otherwise stalls the implicit solver when the bearing stops after sliding. The slope of
the smoothed function stays below 1, so the friction never decreases with the speed (a smoothing with
:math:`\dot q\tanh(\dot q/v_\varepsilon)`, slope up to 1.2, would create a falling friction branch and chatter).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from . import GIMBAL_NAME


# --------------------------------------------------------------------------- Dahl friction


def _dahl_loop(xi: float, n: int = 4000) -> tuple[float, float]:
    """First harmonic of the steady Dahl loop (exponent 1) for :math:`x = A\\sin\\phi`, in units of
    :math:`T_c` and with :math:`\\xi = \\sigma_0 A / T_c`: returns :math:`(k_{eq} A/T_c,\\ \\eta)`.

    Rising branch from :math:`x = -A`: :math:`T = T_c - (T_c + T_1)e^{-\\sigma_0 (x + A)/T_c}`,
    :math:`T_1 = T_c\\tanh\\xi`; the falling branch is its odd mirror."""
    phi = np.linspace(-np.pi / 2, 3 * np.pi / 2, n, endpoint=False)
    u = np.sin(phi)  # x / A
    t1 = np.tanh(xi)
    rising = np.cos(phi) >= 0
    T = np.where(rising, 1 - (1 + t1) * np.exp(-xi * (u + 1)), -(1 - (1 + t1) * np.exp(-xi * (1 - u))))
    a = 2 * np.mean(T * np.sin(phi))
    b = 2 * np.mean(T * np.cos(phi))
    return float(a), float(b / a)


def dahl_describing(sigma0: float, Tc: float, A: float) -> complex:
    """Equivalent complex stiffness :math:`k(1 + j\\eta)` of the Dahl element at amplitude ``A``."""
    a, eta = _dahl_loop(sigma0 * A / Tc)
    k = a * Tc / A
    return complex(k, k * eta)


def dahl_from_hysteretic(k: float, eta: float, A: float) -> tuple[float, float]:
    """:math:`(\\sigma_0, T_c)` of the Dahl element whose describing function at amplitude ``A`` is
    :math:`k(1 + j\\eta)`. :math:`\\eta` depends on :math:`\\xi = \\sigma_0 A/T_c` only and grows
    monotonically from 0 (elastic) to infinity (Coulomb)."""
    xi = brentq(lambda x: _dahl_loop(x)[1] - eta, 1e-4, 1e4)
    a, _ = _dahl_loop(xi)
    Tc = k * A / a
    return xi * Tc / A, Tc


@dataclass(frozen=True)
class DahlElement:
    axis: str
    index: int
    sigma0: float
    Tc: float
    k_eq: float = 0.0
    eta_eq: float = 0.0
    amplitude: float = 0.0
    v_eps: float = 1e-6
    """Velocity scale of the smoothing of :math:`|\\dot q|`, rad/s."""
    Ts: float | None = None
    """Break-away (static) level, N m: the friction can hold up to :math:`T_s` at rest and falls to
    :math:`T_c` in sliding, :math:`g(v) = T_c + (T_s - T_c)e^{-(v/v_s)^2}` in place of :math:`T_c` (Stribeck
    curve of the LuGre model). ``None``: pure Dahl element."""
    vs: float = 1e-3
    """Stribeck velocity, rad/s."""

    def level(self, v: float) -> float:
        if self.Ts is None:
            return self.Tc
        return self.Tc + (self.Ts - self.Tc) * float(np.exp(-(v / self.vs) ** 2))


# --------------------------------------------------------------------------- plant


@dataclass
class GimbalPlant:
    """Nonlinear lumped plant of a ``gimbal_mdof`` configuration."""

    config: Path
    theta0: float = 0.0
    balanced_telescope: bool = True
    zeta_first_flexible: float | None = None
    """Damping ratio of the first flexible mode of each axis; ``None`` keeps the value of the configuration."""
    yaml_overrides: dict = field(default_factory=dict)
    """Values replacing those of the configuration file, by dotted path (e.g.
    ``polar.hub.ring.stiffness_factor``), applied before the model is built."""
    parameter_overrides: dict = field(default_factory=dict)
    """Lumped parameters replaced after the model is built, by symbol (e.g. ``k_{\\theta,4}``)."""
    ground_override: dict = field(default_factory=dict)
    """Ground constraint per axis, ``{axis: {"k": .., "eta": .., "c": ..}}``, replacing the configuration."""
    friction: list[DahlElement] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        from gimbal_mdof import concentrato as cc
        from gimbal_mdof import params as pr
        from gimbal_mdof import simulatore_concentrato as scm
        from gimbal_mdof import sperimentale as spm

        self.config = Path(self.config)
        self._added: set[tuple[str, str, str]] = set()
        if self.yaml_overrides:
            data = pr.carica_config(self.config)
            for path, value in self.yaml_overrides.items():
                node = data
                keys = path.split(".")
                for k in keys[:-1]:
                    node = node[k]
                node[keys[-1]] = value
            self.par = pr.parametri_da_config(self.config, dati=data)
        else:
            self.par = pr.parametri_da_config(self.config)
        self.system = spm.sistemi_per_misure(self.par)[True]
        self.sim = scm.costruisci_simulatore(self.system)
        self.ne = cc.newton_eulero(self.sim.variante)
        self.p = np.array(self.sim.p, dtype=float)
        self._sym = [str(s) for s in self.ne.simboli_parametri]
        self.names = [str(x).replace("\\", "").replace("{", "").replace("}", "") for x in self.ne.q]
        self.n = self.sim.n
        self.E = np.asarray(self.sim.funzioni.E, dtype=float)
        st = self.sim.struttura
        rt = st.ruoli_theta
        self._idx = {
            ("azimuth", "encoder"): st.j_psi_e, ("azimuth", "motor"): st.j_psi_m, ("azimuth", "bearing"): st.j_psi_L,
            ("elevation", "encoder"): st.j_theta_e, ("elevation", "motor"): st.j_theta_m,
            ("elevation", "bearing"): st.n_psi + (rt.index("h") if "h" in rt else rt.index("p")),
            ("elevation", "telescope"): st.j_theta_p,
        }
        if self.balanced_telescope:
            # g = 0 removes the gravity moment and keeps the inertia of the offset CG (m d^2 and products)
            m, g = self.value("m_p"), self.value("g")
            d = np.hypot(self.value("d_x"), self.value("d_z"))
            self.notes.append(f"telescope balanced (g = 0): the static gravity moment m g d = {m * g * d:.2f} N m "
                              f"(CG {1e3 * d:.1f} mm from the polar axis in the CAD) is removed; inertia unchanged")
            self.set_value("g", 0.0)
        self.f_ground = self.par.vincolo_terra
        self._c0 = {name: self.value(name) for name in self._sym if name.startswith("c_")}
        self.zeta_config = float(self.par.zeta)
        self.zeta = self.zeta_config
        if self.zeta_first_flexible is not None:
            self.set_zeta(self.zeta_first_flexible)
            self.notes.append(f"damping of the first flexible mode {self.zeta_config:g} (configuration) "
                              f"replaced by {self.zeta:g}")
        for name, v in self.parameter_overrides.items():
            self.notes.append(f"{name} = {self.value(name):.4g} replaced by {float(v):.4g}")
            self.set_value(name, float(v))
            if name.startswith("c_"):
                self._c0[name] = float(v) * self.zeta_config / self.zeta
        for ax, g in self.ground_override.items():
            self.notes.append(f"ground constraint, {ax}: k = {g['k']:.4g} N m/rad, eta = {g.get('eta', 0.0):.4g}")

    # ---- parameters ------------------------------------------------------------------------------------

    def value(self, name: str) -> float:
        return float(self.p[self._sym.index(name)])

    def set_value(self, name: str, v: float) -> None:
        self.p[self._sym.index(name)] = v

    def set_zeta(self, zeta: float) -> None:
        r"""Stiffness-proportional damping :math:`\mathbf C = \beta\mathbf K`, :math:`\beta = 2\zeta_1/\omega_1`
        per axis: changing :math:`\zeta_1` scales all element dampers by :math:`\zeta_1/\zeta_{1,config}`."""
        for name, c in self._c0.items():
            self.set_value(name, c * zeta / self.zeta_config)
        self.zeta = float(zeta)

    def ground(self, axis: str) -> tuple[float, float, float] | None:
        r"""Ground constraint :math:`k(1 + j\eta) + j\omega c` of ``axis``: override or configuration."""
        if axis in self.ground_override:
            g = self.ground_override[axis]
            return float(g["k"]), float(g.get("eta", 0.0)), float(g.get("c", 0.0))
        if self.f_ground is None:
            return None
        v = self.f_ground.per(GIMBAL_NAME[axis], True)
        return None if v is None else (v.k, v.eta, v.c)

    def index(self, axis: str, node: str) -> int:
        """Coordinate of a node: ``encoder``, ``motor``, ``bearing`` (fork for the azimuth, hub centre or
        telescope node for the polar axis), ``telescope``."""
        return self._idx[(axis, node)]

    def add_inertia(self, axis: str, node: str, J: float, tag: str) -> None:
        """Adds a concentrated inertia (e.g. the motor rotor) to the lumped inertia of a node, once per ``tag``."""
        if (axis, node, tag) in self._added:
            return
        role = {"encoder": "e", "motor": "m"}[node]
        sym = rf"J_{{\{'psi' if axis == 'azimuth' else 'theta'},{role}}}"
        J0 = self.value(sym)
        self.p[self._sym.index(sym)] += J
        self._added.add((axis, node, tag))
        self.notes.append(f"{axis} {node} node: {tag} inertia {J:.3g} kg m^2 added to {J0:.3g} kg m^2")

    def add_dahl(self, axis: str, k: float, eta: float, amplitude: float, v_eps: float = 1e-6) -> DahlElement:
        """Dahl friction on the bearing node of ``axis``, equivalent to :math:`k(1 + j\\eta)` at ``amplitude``."""
        s0, Tc = dahl_from_hysteretic(k, eta, amplitude)
        el = DahlElement(axis, self.index(axis, "bearing"), s0, Tc, k, eta, amplitude, v_eps)
        self.friction.append(el)
        return el

    # ---- dynamics --------------------------------------------------------------------------------------

    @property
    def n_friction(self) -> int:
        return len(self.friction)

    def M(self, q: np.ndarray) -> np.ndarray:
        return np.asarray(self.sim.funzioni.M(q, self.p), dtype=float)

    def h(self, q: np.ndarray, dq: np.ndarray) -> np.ndarray:
        return np.asarray(self.sim.funzioni.h(q, dq, self.p), dtype=float).reshape(-1)

    def accelerations(self, q: np.ndarray, dq: np.ndarray, tau_m: np.ndarray, Tf: np.ndarray | None = None,
                      Q_extra: np.ndarray | None = None) -> np.ndarray:
        Q = self.E @ tau_m - self.h(q, dq)
        if Q_extra is not None:
            Q = Q + Q_extra
        if Tf is not None:
            for el, t in zip(self.friction, Tf):
                Q[el.index] -= t
        return np.linalg.solve(self.M(q), Q)

    def friction_derivatives(self, dq: np.ndarray, Tf: np.ndarray) -> np.ndarray:
        out = np.empty(self.n_friction)
        for k, el in enumerate(self.friction):
            v = dq[el.index]
            out[k] = el.sigma0 * (v - (np.sqrt(v * v + el.v_eps**2) - el.v_eps) * Tf[k] / el.level(v))
        return out

    def q_equilibrium(self) -> np.ndarray:
        """Rest position at the elevation ``theta0`` (all polar nodes at theta0)."""
        q = np.zeros(self.n)
        q[self.sim.struttura.n_psi:] = self.theta0
        return q

    # ---- linearisation ---------------------------------------------------------------------------------

    def linear_matrices(self, eps: float = 1e-7) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """:math:`\\mathbf M`, :math:`\\mathbf C = \\partial\\mathbf h/\\partial\\dot{\\mathbf q}` and
        :math:`\\mathbf K = \\partial\\mathbf h/\\partial\\mathbf q` at rest at ``theta0`` (central differences)."""
        q0, dq0 = self.q_equilibrium(), np.zeros(self.n)
        K, C = np.zeros((self.n, self.n)), np.zeros((self.n, self.n))
        for j in range(self.n):
            e = np.zeros(self.n)
            e[j] = eps
            K[:, j] = (self.h(q0 + e, dq0) - self.h(q0 - e, dq0)) / (2 * eps)
            C[:, j] = (self.h(q0, dq0 + e) - self.h(q0, dq0 - e)) / (2 * eps)
        return self.M(q0), C, K

    def frf(self, axis: str, f: np.ndarray, *, ground: str = "hysteretic") -> dict[str, np.ndarray]:
        """Responses to the motor torque of ``axis`` (the other motor torque is zero): rotation of the encoder
        and of the motor node per N m. ``ground``: ``hysteretic`` adds the H6 constraint :math:`k(1 + j\\eta)` of
        the configuration on the bearing nodes (frequency domain, as in H6), ``dahl`` the small-amplitude
        stiffness :math:`\\sigma_0` of the Dahl elements, ``none`` nothing."""
        M, C, K = self.linear_matrices()
        K = K.astype(complex)
        if ground == "hysteretic":
            for ax in ("azimuth", "elevation"):
                g = self.ground(ax)
                if g is not None:
                    j = self.index(ax, "bearing")
                    K[j, j] += g[0] * (1 + 1j * g[1])
                    C[j, j] += g[2]
        elif ground == "dahl":
            for el in self.friction:
                K[el.index, el.index] += el.sigma0
        j_in = 0 if axis == "azimuth" else 1
        out = {"encoder": [], "motor": []}
        ie, im = self.index(axis, "encoder"), self.index(axis, "motor")
        for fk in np.atleast_1d(f):
            s = 2j * np.pi * fk
            x = np.linalg.solve(M * s**2 + C * s + K, self.E[:, j_in])
            out["encoder"].append(x[ie])
            out["motor"].append(x[im])
        return {k: np.array(v) for k, v in out.items()}
