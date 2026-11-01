r"""Fit of plant parameters on a measured plant frequency response (encoder velocity over current command).

Parameters (names used in the configuration and in the results):

* ``k_b``      -- spring between the hub centre H and the telescope P, :math:`k_{\theta,4}` [N m/rad];
* ``ring``     -- stiffness factor :math:`\kappa_r` of the support ring (``polar.hub.ring.stiffness_factor`` of
  gimbal_mdof): the ring is in series in the two hub halves, so that :math:`1/k_{\theta,i} = a_i + b_i/\kappa_r`
  for :math:`k_{\theta,1}` (encoder - hub) and :math:`k_{\theta,2}` (hub - motor), with :math:`a_i, b_i` from two
  builds of the model (:class:`RingMap`; checked exact on other values);
* ``kg_az``, ``eta_az``, ``kg_el``, ``eta_el`` -- ground constraint :math:`k(1 + j\eta)` on the bearing node
  of each axis [N m/rad, -];
* ``zeta_az``, ``zeta_el`` -- damping of each axis as in H6: stiffness proportional,
  :math:`c_i = \beta\,k_i` with :math:`\beta = 2\zeta_1/\omega_1`, :math:`\omega_1` the first flexible mode of
  the axis in the base model (so :math:`\zeta_1` is the damping ratio of that mode); the dampers follow the
  fitted stiffnesses;
* ``c_*``      -- single element dampers of the lumped model [N m s/rad], by symbol (``c_{\psi,1}``, ...).

The residual of one frequency point is :math:`\ln(\hat G/G) = \ln|\hat G/G| + j\,\angle(\hat G/G)`, real and
imaginary parts, on the measurement points in a band; the parameters are fitted on their logarithms.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .closed_loop import ClosedLoop, loop_frf

K1, K2, K4 = r"k_{\theta,1}", r"k_{\theta,2}", r"k_{\theta,4}"
GROUND = {"kg_az": ("azimuth", "k"), "eta_az": ("azimuth", "eta"), "kg_el": ("elevation", "k"),
          "eta_el": ("elevation", "eta")}


@dataclass(frozen=True)
class RingMap:
    """:math:`1/k_{\\theta,i} = a_i + b_i/\\kappa_r` for the two hub springs."""

    a: dict[str, float]
    b: dict[str, float]

    @classmethod
    def from_config(cls, config) -> "RingMap":
        from .plant import GimbalPlant

        k = {}
        for f in (1000.0, 1.0):
            pl = GimbalPlant(config, yaml_overrides={"polar.hub.ring.stiffness_factor": f})
            k[f] = {s: pl.value(s) for s in (K1, K2)}
        a, b = {}, {}
        for s in (K1, K2):
            sol = np.linalg.solve([[1.0, 1.0 / 1000.0], [1.0, 1.0]], [1.0 / k[1000.0][s], 1.0 / k[1.0][s]])
            a[s], b[s] = float(sol[0]), float(sol[1])
        return cls(a, b)

    def k(self, s: str, ring: float) -> float:
        return 1.0 / (self.a[s] + self.b[s] / ring)


@dataclass(frozen=True)
class DampingMap:
    r""":math:`\beta_1 = c_i/(\zeta_1 k_i)` per axis of the base model (stiffness-proportional damping)."""

    beta_per_zeta: dict[str, float]
    zeta: dict[str, float]

    @classmethod
    def from_plant(cls, pl) -> "DampingMap":
        b, z = {}, {}
        for ax, tag in (("azimuth", "psi"), ("elevation", "theta")):
            names = [n for n in pl._sym if n.startswith(f"c_{{\\{tag}")]
            ratios = [pl.value(n) / pl.value("k" + n[1:]) for n in names]
            b[ax] = float(np.mean(ratios)) / pl.zeta
            z[ax] = pl.zeta
        return cls(b, z)

    def set(self, pl, axis: str, zeta: float) -> None:
        tag = "psi" if axis == "azimuth" else "theta"
        for n in [n for n in pl._sym if n.startswith(f"c_{{\\{tag}")]:
            pl.set_value(n, self.beta_per_zeta[axis] * zeta * pl.value("k" + n[1:]))
        self.zeta[axis] = float(zeta)


def apply(cl: ClosedLoop, values: dict[str, float], ring_map: RingMap | None,
          damping: DampingMap | None = None) -> None:
    """Writes the parameters into the plant of the loop (lumped values and ground constraint); with ``damping``
    the stiffness-proportional dampers are updated after the stiffnesses."""
    pl = cl.plant
    for name, v in values.items():
        if name == "k_b":
            pl.set_value(K4, v)
        elif name == "ring":
            for s in (K1, K2):
                pl.set_value(s, ring_map.k(s, v))
        elif name in ("zeta_az", "zeta_el"):
            continue
        elif name in GROUND:
            ax, key = GROUND[name]
            g = dict(pl.ground_override.get(ax) or {})
            if not g:
                base = pl.ground(ax) or (0.0, 0.0, 0.0)
                g = {"k": base[0], "eta": base[1], "c": base[2]}
            g[key] = float(v)
            pl.ground_override[ax] = g
        else:
            pl.set_value(name, v)
    if damping is not None:
        for ax, key in (("azimuth", "zeta_az"), ("elevation", "zeta_el")):
            damping.set(pl, ax, values.get(key, damping.zeta[ax]))


def current_values(cl: ClosedLoop, names: list[str], ring_map: RingMap | None,
                   damping: DampingMap | None = None) -> dict[str, float]:
    pl = cl.plant
    out = {}
    for n in names:
        if n == "k_b":
            out[n] = pl.value(K4)
        elif n == "ring":
            # invert the map on k_theta,1
            out[n] = float(ring_map.b[K1] / (1.0 / pl.value(K1) - ring_map.a[K1]))
        elif n in GROUND:
            ax, key = GROUND[n]
            g = pl.ground(ax)
            out[n] = g[0] if key == "k" else g[1]
        elif n in ("zeta_az", "zeta_el"):
            out[n] = damping.zeta["azimuth" if n == "zeta_az" else "elevation"]
        else:
            out[n] = pl.value(n)
    return out


def residual(G_model: np.ndarray, G_meas: np.ndarray) -> np.ndarray:
    r = np.log(G_model / G_meas)
    return np.concatenate([r.real, np.angle(G_model / G_meas)])


@dataclass
class FitResult:
    axis: str
    names: list[str]
    start: dict[str, float]
    values: dict[str, float]
    band: tuple[float, float]
    cost_start: float
    cost: float
    n_eval: int
    history: list = field(default_factory=list)


def fit(cl: ClosedLoop, axis: str, f: np.ndarray, G: np.ndarray, names: list[str], bounds: dict[str, tuple],
        band: tuple[float, float], ring_map: RingMap | None, damping: DampingMap | None = None, *,
        loss: str = "soft_l1", f_scale: float = 0.5, max_nfev: int = 200,
        prior: dict[str, tuple[float, float]] | None = None) -> FitResult:
    """Least squares on the logarithms of ``names`` against the measured plant ``G`` in ``band``. ``prior``:
    ``{name: (value, weight)}`` adds the residual ``weight * ln(x / value)`` (regularisation of the parameters
    that the band does not determine, e.g. dampers of modes outside it)."""
    from scipy.optimize import least_squares

    m = (f >= band[0]) & (f <= band[1])
    fb, Gb = f[m], G[m]
    start = current_values(cl, names, ring_map, damping)
    lo = np.log([bounds[n][0] for n in names])
    hi = np.log([bounds[n][1] for n in names])
    x0 = np.clip(np.log([start[n] for n in names]), lo + 1e-9, hi - 1e-9)
    hist = []

    def res(x: np.ndarray) -> np.ndarray:
        apply(cl, dict(zip(names, np.exp(x))), ring_map, damping)
        r = residual(loop_frf(cl, axis, fb).plant, Gb)
        if prior:
            r = np.concatenate([r, [w * (xi - np.log(prior[n][0])) for n, xi in zip(names, x)
                                    for w in [prior.get(n, (1.0, 0.0))[1]] if n in prior]])
        hist.append(float(np.sum(r**2)))
        return r

    r0 = res(x0)
    sol = least_squares(res, x0, bounds=(lo, hi), loss=loss, f_scale=f_scale, diff_step=1e-3, max_nfev=max_nfev,
                        x_scale=1.0)
    values = dict(zip(names, np.exp(sol.x)))
    apply(cl, values, ring_map, damping)
    return FitResult(axis, names, start, values, band, float(np.sum(r0**2)), float(np.sum(sol.fun**2)),
                     int(sol.nfev), hist)
