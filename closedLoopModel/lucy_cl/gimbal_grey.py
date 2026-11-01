r"""Adapter between the lumped plant of the closed loop and the grey-box estimator (:mod:`lucy_cl.greyest`).

For given parameter values the adapter writes the plant (lumped inertias, springs, dampers, telescope inertia,
ground constraint), builds the linear matrix problem about rest

.. math:: \bigl(-\omega^2\mathbf M + j\omega\mathbf C + \mathbf K + \mathbf K_g(1 + j\eta)\bigr)\,\mathbf x = \mathbf e\,\tau

and returns the plant as WorkBench measures it, encoder velocity over current command (:func:`loop_frf`: motor
current loop with back EMF, motor viscous friction, drive delay). The ground constraint is hysteretic, so the model
is defined in frequency only (``GreyModel.frfun``).

Parameters of each axis (``axis_parameters``), with the kind that sets the bounds of the estimate:

* ``inertia``    -- lumped inertias of the chain (rotor included at the motor node), :math:`J^{ax}` (hub and ring
  about :math:`\hat z`), telescope :math:`I_{zz}` (azimuth) and :math:`I_{yy}` (elevation): CAD or geometry;
* ``geometric``  -- springs computed from the shaft geometry; for the two hub springs of the elevation
  (:math:`k_{\theta,1}`, :math:`k_{\theta,2}`) only the hub and shaft part: the support ring is in series and is a
  parameter of its own, :math:`1/k_{\theta,i} = 1/k^{geo}_{\theta,i} + b_i/\kappa_r` (:class:`RingMap`);
* ``free``       -- springs not determined by the geometry: telescope-hub :math:`k_{\theta,4}`, telescope-fork
  :math:`k_{\psi,3}`, ring factor :math:`\kappa_r`, ground constraint :math:`k_g`, :math:`\eta`;
* ``damping``    -- element dampers and loss factor of the ground constraint.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Mapping

import numpy as np

from .closed_loop import ClosedLoop, loop_frf
from .driver import Drive
from .greyest import GreyModel, GreyParameter
from .plant_fit import K1, K2, RingMap

Kind = Literal["inertia", "geometric", "free", "damping"]

GROUND = {"azimuth": ("kg_az", "eta_az"), "elevation": ("kg_el", "eta_el")}


@dataclass(frozen=True)
class PlantParameter:
    name: str
    kind: Kind
    nominal: float
    """CAD or geometric value (``inertia``, ``geometric``) or starting value (``free``, ``damping``)."""
    note: str = ""


def axis_parameters(cl: ClosedLoop, axis: str, ring_map: RingMap | None) -> list[PlantParameter]:
    """Parameters of ``axis`` with their kind and current value. Nodes without inertia (telescope nodes of the
    chains) are left out; their inertia is that of the telescope (``I_zz``, ``I_yy``)."""
    pl = cl.plant
    v = pl.value
    tag = "psi" if axis == "azimuth" else "theta"
    roles = [r for r in pl.sim.struttura.ruoli_psi] if axis == "azimuth" else list(pl.sim.struttura.ruoli_theta)
    out = []
    for r in roles:
        s = rf"J_{{\{tag},{r}}}"
        if v(s) > 0:
            out.append(PlantParameter(s, "inertia", v(s), "rotor included" if r == "m" else ""))
    if axis == "azimuth":
        out += [PlantParameter("J^{ax}", "inertia", v("J^{ax}"), "hub and ring about z"),
                PlantParameter("I_zz", "inertia", v("I_zz"), "telescope about z at theta = 0")]
    else:
        out.append(PlantParameter("I_yy", "inertia", v("I_yy"), "telescope about the polar axis"))
    n_springs = len([s for s in pl._sym if s.startswith(rf"k_{{\{tag},")])
    for i in range(1, n_springs + 1):
        s = rf"k_{{\{tag},{i}}}"
        if axis == "azimuth":
            kind: Kind = "free" if i == 3 else "geometric"
            note = "telescope - fork" if i == 3 else "shafts"
            out.append(PlantParameter(s, kind, v(s), note))
        elif s in (K1, K2):
            assert ring_map is not None
            out.append(PlantParameter(s + "^{geo}", "geometric", 1.0 / ring_map.a[s], "hub and shafts, ring excluded"))
        else:
            kind = "free" if i == 4 else "geometric"
            out.append(PlantParameter(s, kind, v(s), "telescope - hub" if i == 4 else "shafts"))
    if axis == "elevation":
        assert ring_map is not None
        ring = float(ring_map.b[K1] / (1.0 / v(K1) - ring_map.a[K1]))
        out.append(PlantParameter("ring", "free", ring, "stiffness factor of the support ring"))
    for i in range(1, n_springs + 1):
        s = rf"c_{{\{tag},{i}}}"
        out.append(PlantParameter(s, "damping", v(s)))
    g = pl.ground(axis) or (1.0, 0.1, 0.0)
    out += [PlantParameter(GROUND[axis][0], "free", g[0], "ground constraint stiffness"),
            PlantParameter(GROUND[axis][1], "damping", g[1], "ground constraint loss factor")]
    return out


@dataclass(frozen=True)
class Bounds:
    """Bounds of the estimate per kind: relative band around the nominal value for ``inertia`` and ``geometric``,
    absolute ranges for the free ones."""

    relative: float = 0.2
    free_stiffness: tuple[float, float] = (1.0, 1.0e9)
    ring: tuple[float, float] = (1.0e-3, 1.0e6)
    damper: tuple[float, float] = (1.0e-6, 1.0e4)
    eta: tuple[float, float] = (1.0e-4, 1.0e2)

    def of(self, p: PlantParameter) -> tuple[float, float]:
        if p.kind in ("inertia", "geometric"):
            return p.nominal * (1 - self.relative), p.nominal * (1 + self.relative)
        if p.name == "ring":
            return self.ring
        if p.name.startswith("eta"):
            return self.eta
        if p.kind == "damping":
            return self.damper
        return self.free_stiffness


@dataclass
class PlantGreyModel:
    """Writes parameter values into the plant of ``cl`` and evaluates the plant response of ``axis``."""

    cl: ClosedLoop
    axis: str
    drive: Drive
    ring_map: RingMap | None
    params: list[PlantParameter]
    evaluations: int = field(default=0)

    def apply(self, values: Mapping[str, float]) -> None:
        pl = self.cl.plant
        ring = values.get("ring")
        for name, x in values.items():
            if name == "ring":
                continue
            if name.endswith("^{geo}"):
                s = name[: -len("^{geo}")]
                assert self.ring_map is not None and ring is not None
                pl.set_value(s, 1.0 / (1.0 / x + self.ring_map.b[s] / ring))
            elif name in GROUND[self.axis]:
                g = dict(pl.ground_override.get(self.axis) or {})
                if not g:
                    base = pl.ground(self.axis) or (0.0, 0.0, 0.0)
                    g = {"k": base[0], "eta": base[1], "c": base[2]}
                g["k" if name.startswith("kg") else "eta"] = float(x)
                pl.ground_override[self.axis] = g
            else:
                pl.set_value(name, float(x))

    def frfun(self, values: Mapping[str, float], w: np.ndarray) -> np.ndarray:
        self.apply(values)
        self.evaluations += 1
        return loop_frf(self.cl, self.axis, np.asarray(w) / (2 * np.pi), drive=self.drive).plant

    def matrices(self, values: Mapping[str, float] | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(M, C, K)`` of the linear problem (``K`` complex, ground constraint included), for inspection."""
        if values:
            self.apply(values)
        pl = self.cl.plant
        M, C, K = pl.linear_matrices()
        K = K.astype(complex)
        for ax in ("azimuth", "elevation"):
            g = pl.ground(ax)
            if g is not None:
                j = pl.index(ax, "bearing")
                K[j, j] += g[0] * (1 + 1j * g[1])
                C[j, j] += g[2]
        return M, C, K

    def modes(self, values: Mapping[str, float] | None = None, f_max: float = np.inf) -> list[tuple[float, float]]:
        r"""``(f [Hz], zeta)`` of the modes of the axis up to ``f_max``: eigenvalues of
        :math:`\mathbf M\ddot{\mathbf q} + \mathbf C\dot{\mathbf q} + \mathbf K\mathbf q = 0` restricted to the
        coordinates of the axis (the axes are uncoupled at rest, products of inertia ~1e-6), with the stiffness of
        the ground constraint (not its hysteretic loss, which is bearing friction). A real root (overdamped motion)
        counts as :math:`\zeta = 1` at its own magnitude; roots below 0.01 Hz (rigid motion) are left out."""
        if values:
            self.apply(values)
        pl = self.cl.plant
        M, C, K = pl.linear_matrices()
        for ax in ("azimuth", "elevation"):
            g = pl.ground(ax)
            if g is not None:
                j = pl.index(ax, "bearing")
                K[j, j] += g[0]
                C[j, j] += g[2]
        n_psi = pl.sim.struttura.n_psi
        idx = np.arange(n_psi) if self.axis == "azimuth" else np.arange(n_psi, M.shape[0])
        M, C, K = (X[np.ix_(idx, idx)] for X in (M, C, K))
        n = M.shape[0]
        Mi = np.linalg.inv(M)
        lam = np.linalg.eigvals(np.block([[np.zeros((n, n)), np.eye(n)], [-Mi @ K, -Mi @ C]]))
        out = []
        for r in lam[np.imag(lam) >= 0]:
            f = abs(r) / (2 * np.pi)
            if 0.01 < f <= f_max:
                out.append((float(f), float(-r.real / abs(r)) if abs(r.imag) > 1e-9 * abs(r) else 1.0))
        return sorted(out)

    def damping_penalty(self, zeta_max: float, f_max: float, weight: float = 1.0e3):
        """Soft constraint for :class:`GreyModel`: ``weight * max(0, zeta - zeta_max)`` for every mode up to
        ``f_max`` (the modes above do not shape the response in the band and are left free)."""
        n_roots = 2 * self.cl.plant.n  # upper bound of the roots with Im >= 0 (real roots count one by one)

        def penalty(values: Mapping[str, float]) -> np.ndarray:
            # fixed length: the number of modes below f_max changes during the search
            v = np.sort([max(zeta - zeta_max, 0.0) for _, zeta in self.modes(values, f_max)])[::-1]
            out = np.zeros(n_roots)
            out[: v.size] = v[:n_roots]
            return weight * out
        return penalty

    def grey_model(self, bounds: Bounds, start: Mapping[str, float] | None = None,
                   fixed: tuple[str, ...] = (), penalty=None) -> GreyModel:
        """Grey-box model with the parameters of the axis; ``start`` replaces the starting values of the free
        ones (clipped into the bounds), ``fixed`` keeps parameters at their nominal value, ``penalty`` is a soft
        constraint (e.g. :meth:`damping_penalty`)."""
        gp = []
        for p in self.params:
            lo, hi = bounds.of(p)
            x = float(np.clip((start or {}).get(p.name, p.nominal), lo, hi))
            gp.append(GreyParameter(p.name, x, lo, hi, free=p.name not in fixed, scale="log"))
        return GreyModel(odefun=None, parameters=tuple(gp), frfun=self.frfun, penalty_fun=penalty,
                         name=f"LUCY {self.axis} plant")
