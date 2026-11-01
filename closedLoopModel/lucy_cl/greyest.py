r"""Grey-box estimation of the parameters of a linear ODE model, in the style of MATLAB ``greyest``.

The model is a user function that returns the state-space matrices for given parameter values, as ``idgrey``:

.. math::
    \dot x = A(\theta)\,x + B(\theta)\,u,\qquad y = C(\theta)\,x + D(\theta)\,u

(continuous, ``Ts = 0``) or :math:`x_{k+1} = A x_k + B u_k` (discrete, ``Ts > 0``). :func:`greyest` estimates the
free parameters :math:`\theta` by nonlinear least squares on

* frequency-response data (:class:`FrequencyData`, as ``idfrd``): error :math:`\hat G(j\omega) - G(j\omega)`
  (``frequency_error = "complex"``, as MATLAB) or :math:`\ln(\hat G/G)` (``"log"``: dB and phase, for data with a
  wide dynamic range such as the plant measurements of WorkBench);
* time-domain data (:class:`TimeData`, as ``iddata``): simulation error :math:`\hat y(t) - y(t)` with the input
  held constant between samples (ZOH), optionally with the initial state estimated.

Several experiments can be passed together (list of data of the same kind). Correspondence with MATLAB:

=========================== =================================================================================
MATLAB                      here
=========================== =================================================================================
``idgrey(odefun, p, ...)``  :class:`GreyModel` (``odefun(p, Ts, *aux) -> (A, B, C, D)``)
``Structure.Parameters``    :class:`GreyParameter` (``value``, ``minimum``, ``maximum``, ``free``)
``greyestOptions``          :class:`GreyestOptions`
``greyest(data, init_sys)`` :func:`greyest` -> :class:`GreyEstimate` (``.model``, ``.report``)
``sys.Report.Fit``          :class:`GreyestReport` (``fit_percent``, ``mse``, ``fpe``, ``aic``, ``aicc``, ``bic``)
``getcov`` / ``getpvec``    ``report.covariance``, ``report.std``, ``GreyEstimate.parameters``
``compare(data, sys)``      :func:`compare`
=========================== =================================================================================

Additions for the models of this project: a frequency-only model (``GreyModel.frfun``) for responses that are not
rational in :math:`s`, e.g. the hysteretic ground constraint :math:`k(1 + j\eta)`; a robust loss (``loss``) for
isolated bad points of the measurements; parameters searched on their logarithm (``scale = "log"``) when they span
decades.

Parameter statistics: with :math:`J` the Jacobian of the weighted residuals with respect to the free parameters
at the estimate (central differences, without the robust loss) and :math:`s^2 = r^Tr/(N_r - d)`, either the
classic :math:`s^2 (J^TJ)^{-1}` of MATLAB (residuals of equal variance) or, by default, the sandwich estimate that
stays valid when the noise changes along the data (``GreyestOptions.covariance``). Both assume uncorrelated
residuals: for time data whose residuals are correlated they are optimistic.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Callable, Literal, Mapping, Sequence

import numpy as np
import scipy.linalg as sla

Matrices = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


# --------------------------------------------------------------------------- model


@dataclass(frozen=True)
class GreyParameter:
    """One parameter of the grey-box model (``Structure.Parameters`` of ``idgrey``)."""

    name: str
    value: float
    minimum: float = -np.inf
    maximum: float = np.inf
    free: bool = True
    scale: Literal["auto", "lin", "log"] = "auto"
    """Coordinate of the search: ``log`` (only for positive bounds) searches :math:`\\ln\\theta`, suited to
    stiffnesses and dampers spanning decades; ``auto`` = ``log`` when ``minimum > 0``."""

    def __post_init__(self) -> None:
        if not self.minimum <= self.value <= self.maximum:
            raise ValueError(f"parameter {self.name}: value {self.value} outside [{self.minimum}, {self.maximum}]")
        if self.coordinate == "log" and self.minimum <= 0:
            raise ValueError(f"parameter {self.name}: scale 'log' needs minimum > 0")

    @property
    def coordinate(self) -> str:
        if self.scale == "auto":
            return "log" if self.minimum > 0 else "lin"
        return self.scale


@dataclass(frozen=True)
class GreyModel:
    """Linear grey-box model (``idgrey``).

    ``odefun(p, Ts, *aux)`` receives the parameter values by name and returns ``(A, B, C, D)``; ``Ts = 0`` is a
    continuous model. ``frfun(p, w, *aux)``, optional, returns the frequency response ``(nw, ny, nu)`` directly
    (``w`` in rad/s): it is used for frequency data in place of the matrices, so that models that are not rational
    in ``s`` can be estimated (time data still need ``odefun``).
    """

    odefun: Callable[..., Matrices] | None
    parameters: tuple[GreyParameter, ...]
    Ts: float = 0.0
    aux: tuple = ()
    frfun: Callable[..., np.ndarray] | None = None
    penalty_fun: Callable[..., np.ndarray] | None = None
    """Soft constraints (not in MATLAB greyest): ``penalty_fun(p, *aux)`` returns residuals appended to the error,
    zero when the constraints hold and growing with the violation (e.g. ``w * max(0, zeta - zeta_max)``)."""
    input_names: tuple[str, ...] = ()
    output_names: tuple[str, ...] = ()
    name: str = ""

    def __post_init__(self) -> None:
        if self.odefun is None and self.frfun is None:
            raise ValueError("GreyModel needs odefun or frfun")
        names = [p.name for p in self.parameters]
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate parameter names: {names}")

    @classmethod
    def create(cls, odefun: Callable[..., Matrices] | None, parameters: Mapping[str, Any], Ts: float = 0.0,
               **kw) -> "GreyModel":
        """``parameters``: ``{name: value}`` or ``{name: GreyParameter | dict(value=..., minimum=..., ...)}``."""
        out = []
        for n, v in parameters.items():
            if isinstance(v, GreyParameter):
                out.append(replace(v, name=n))
            elif isinstance(v, Mapping):
                out.append(GreyParameter(n, **v))
            else:
                out.append(GreyParameter(n, float(v)))
        return cls(odefun=odefun, parameters=tuple(out), Ts=Ts, **kw)

    # ---- parameters

    @property
    def names(self) -> list[str]:
        return [p.name for p in self.parameters]

    @property
    def values(self) -> dict[str, float]:
        return {p.name: p.value for p in self.parameters}

    @property
    def free_names(self) -> list[str]:
        return [p.name for p in self.parameters if p.free]

    def parameter(self, name: str) -> GreyParameter:
        return self.parameters[self.names.index(name)]

    def with_values(self, values: Mapping[str, float]) -> "GreyModel":
        return replace(self, parameters=tuple(replace(p, value=float(values.get(p.name, p.value)))
                                              for p in self.parameters))

    def with_free(self, free: Sequence[str]) -> "GreyModel":
        """Copy with exactly the parameters ``free`` free (the others fixed)."""
        unknown = set(free) - set(self.names)
        if unknown:
            raise KeyError(f"unknown parameters {sorted(unknown)}")
        return replace(self, parameters=tuple(replace(p, free=p.name in free) for p in self.parameters))

    # ---- responses

    def ssdata(self, values: Mapping[str, float] | None = None) -> Matrices:
        if self.odefun is None:
            raise ValueError("model defined only in frequency (frfun): no state-space matrices")
        p = self.values if values is None else {**self.values, **values}
        A, B, C, D = (np.atleast_2d(np.asarray(X, dtype=float)) for X in self.odefun(p, self.Ts, *self.aux))
        n = A.shape[0]
        if A.shape != (n, n) or B.shape[0] != n or C.shape[1] != n or D.shape != (C.shape[0], B.shape[1]):
            raise ValueError(f"inconsistent matrices: A {A.shape}, B {B.shape}, C {C.shape}, D {D.shape}")
        return A, B, C, D

    def frequency_response(self, w: np.ndarray, values: Mapping[str, float] | None = None) -> np.ndarray:
        """``(nw, ny, nu)`` complex response at ``w`` rad/s."""
        w = np.asarray(w, dtype=float)
        if self.frfun is not None:
            p = self.values if values is None else {**self.values, **values}
            G = np.asarray(self.frfun(p, w, *self.aux), dtype=complex)
            return G.reshape(len(w), 1, 1) if G.ndim == 1 else G
        return state_space_response(*self.ssdata(values), w, self.Ts)

    def simulate(self, t: np.ndarray, u: np.ndarray, x0: np.ndarray | None = None,
                 values: Mapping[str, float] | None = None) -> np.ndarray:
        """Output ``(N, ny)`` for the input ``u`` ``(N, nu)`` on the uniform grid ``t`` (ZOH between samples)."""
        A, B, C, D = self.ssdata(values)
        dt = _uniform_step(t)
        if self.Ts > 0 and not np.isclose(dt, self.Ts, rtol=1e-6):
            raise ValueError(f"discrete model with Ts = {self.Ts} s, data sampled at {dt} s")
        Ad, Bd = (A, B) if self.Ts > 0 else zoh(A, B, dt)
        return simulate_discrete(Ad, Bd, C, D, np.asarray(u, dtype=float).reshape(len(t), -1),
                                 np.zeros(A.shape[0]) if x0 is None else np.asarray(x0, dtype=float))


def state_space_response(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray, w: np.ndarray,
                         Ts: float = 0.0) -> np.ndarray:
    r""":math:`C(sI - A)^{-1}B + D`, :math:`s = j\omega` (continuous) or :math:`e^{j\omega T_s}` (discrete);
    ``(nw, ny, nu)``."""
    s = np.exp(1j * w * Ts) if Ts > 0 else 1j * w
    n = A.shape[0]
    if n == 0:
        return np.broadcast_to(D.astype(complex), (len(w),) + D.shape).copy()
    X = np.linalg.solve(s[:, None, None] * np.eye(n) - A, np.broadcast_to(B.astype(complex), (len(w),) + B.shape))
    return C @ X + D


def zoh(A: np.ndarray, B: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Zero-order-hold discretisation, from the exponential of the augmented matrix."""
    n, m = B.shape
    E = sla.expm(np.block([[A, B], [np.zeros((m, n + m))]]) * dt)
    return E[:n, :n], E[:n, n:]


def simulate_discrete(Ad: np.ndarray, Bd: np.ndarray, C: np.ndarray, D: np.ndarray, u: np.ndarray,
                      x0: np.ndarray) -> np.ndarray:
    N = u.shape[0]
    x = np.empty((N, Ad.shape[0]))
    xk = x0.astype(float)
    Bu = u @ Bd.T
    for k in range(N):
        x[k] = xk
        xk = Ad @ xk + Bu[k]
    return x @ C.T + u @ D.T


def _uniform_step(t: np.ndarray) -> float:
    d = np.diff(np.asarray(t, dtype=float))
    if d.size == 0 or np.ptp(d) > 1e-6 * d.mean():
        raise ValueError("time data must be uniformly sampled")
    return float(d.mean())


# --------------------------------------------------------------------------- data


@dataclass(frozen=True)
class FrequencyData:
    """Frequency-response data (``idfrd``): ``f`` in Hz, ``G`` complex ``(nf,)`` (SISO) or ``(nf, ny, nu)``.
    ``weights`` ``(nf,)``, optional, multiplies the error of each frequency (a measure of confidence, e.g. the
    coherence)."""

    f: np.ndarray
    G: np.ndarray
    weights: np.ndarray | None = None
    name: str = ""

    def __post_init__(self) -> None:
        G = np.asarray(self.G, dtype=complex)
        object.__setattr__(self, "f", np.asarray(self.f, dtype=float))
        object.__setattr__(self, "G", G.reshape(len(self.f), 1, 1) if G.ndim == 1 else G)
        if self.G.shape[0] != len(self.f):
            raise ValueError("FrequencyData: G and f have different lengths")

    @property
    def w(self) -> np.ndarray:
        return 2 * np.pi * self.f


@dataclass(frozen=True)
class TimeData:
    """Time-domain data (``iddata``): ``t`` uniform, ``u`` ``(N, nu)``, ``y`` ``(N, ny)``."""

    t: np.ndarray
    u: np.ndarray
    y: np.ndarray
    name: str = ""

    def __post_init__(self) -> None:
        t = np.asarray(self.t, dtype=float)
        object.__setattr__(self, "t", t)
        object.__setattr__(self, "u", np.asarray(self.u, dtype=float).reshape(len(t), -1))
        object.__setattr__(self, "y", np.asarray(self.y, dtype=float).reshape(len(t), -1))
        _uniform_step(t)


Data = FrequencyData | TimeData


# --------------------------------------------------------------------------- options and results


@dataclass(frozen=True)
class GreyestOptions:
    """Options of :func:`greyest` (``greyestOptions``)."""

    search_method: Literal["auto", "lsqnonlin", "lm", "dogbox"] = "auto"
    """``lsqnonlin``: trust-region reflective with bounds (scipy ``trf``); ``lm``: Levenberg-Marquardt, only without
    bounds and without robust loss; ``auto``: ``lsqnonlin``."""
    weighting_filter: Any = None
    """Frequency data: ``None``; a band ``[f1, f2]`` Hz or a list of bands ``[[f1, f2], ...]`` (passbands, as
    MATLAB); a vector of weights as long as the data (same frequencies in every experiment; otherwise use
    ``FrequencyData.weights``); ``"inv"`` (:math:`1/|G|`) or ``"invsqrt"`` (:math:`1/\\sqrt{|G|}`)."""
    frequency_error: Literal["complex", "log"] = "complex"
    output_weight: Any = None
    """``None`` (identity), ``"normalize"`` (each output divided by its spread in the data) or a vector ``(ny,)``."""
    initial_state: Literal["zero", "estimate"] = "zero"
    """Time data: ``estimate`` adds the initial state of each experiment to the unknowns."""
    enforce_stability: bool = False
    """Penalises estimates with unstable poles (requires ``odefun``)."""
    loss: Literal["linear", "soft_l1", "huber", "cauchy", "arctan"] = "linear"
    f_scale: float = 1.0
    """Residual size at which the robust loss starts to act."""
    regularization_lambda: float = 0.0
    regularization_R: np.ndarray | None = None
    """Weight matrix of the free parameters (identity if ``None``); the penalty is
    :math:`\\lambda(\\theta - \\theta^*)^TR(\\theta - \\theta^*)` in the search coordinates (logarithm for the
    ``log`` parameters), as ``Regularization`` of MATLAB."""
    regularization_nominal: Literal["zero", "model"] | Mapping[str, float] = "model"
    max_iterations: int = 200
    """Maximum number of model evaluations of the search (scipy ``max_nfev``)."""
    tolerance: float = 1e-8
    """Relative tolerance on cost, step and gradient."""
    diff_step: float = 1e-4
    """Relative step of the finite differences."""
    covariance: Literal["robust", "classic"] = "robust"
    """``classic``: :math:`s^2(J^TJ)^{-1}` as MATLAB, valid when all residuals have the same variance (e.g. the
    ``log`` error with noise proportional to the response); ``robust``: heteroscedasticity-consistent sandwich
    :math:`(J^TJ)^{-1}J^T\\mathrm{diag}(r^2)J(J^TJ)^{-1}\\,N/(N - d)`, valid also when the noise changes along the
    data (e.g. the ``complex`` error near a resonance). Neither accounts for correlated residuals."""
    display: bool = False


@dataclass(frozen=True)
class ParameterEstimate:
    name: str
    value: float
    initial: float
    free: bool
    std: float
    """Standard deviation, nan for fixed parameters or when the information matrix is singular."""
    at_bound: bool

    @property
    def ci95(self) -> tuple[float, float]:
        return (self.value - 1.96 * self.std, self.value + 1.96 * self.std)


@dataclass(frozen=True)
class GreyestReport:
    """Quality of the estimate (``Report`` of the estimated ``idgrey``)."""

    method: str
    status: str
    termination: str
    n_evaluations: int
    fit_percent: np.ndarray
    """NRMSE fit per experiment and output, ``(n_experiments, ny)``: :math:`100(1 - \\lVert y - \\hat y\\rVert /
    \\lVert y - \\bar y\\rVert)`, on the data used (inside the weighting bands), complex for frequency data."""
    mse: float
    """Mean squared error, like the following criteria computed on the weighted errors (weighting filter and
    output weight), without the regularization."""
    fpe: float
    aic: float
    aicc: float
    bic: float
    loss_start: float
    loss: float
    """Weighted sum of squared residuals (with regularization) at the start and at the estimate."""
    n_data: int
    n_free: int
    parameters: tuple[ParameterEstimate, ...]
    covariance: np.ndarray
    """Covariance of the free parameters, in their units, order of ``free_names``."""
    correlation: np.ndarray
    free_names: tuple[str, ...]
    jacobian_condition: float
    """Condition number of the column-scaled Jacobian: large values mean parameters the data do not separate."""
    residuals: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))

    @property
    def std(self) -> dict[str, float]:
        return {p.name: p.std for p in self.parameters if p.free}

    def summary(self) -> str:
        lines = [f"greyest: {self.method}, {self.status} ({self.termination}), {self.n_evaluations} evaluations",
                 "fit [%]: " + "; ".join(", ".join(f"{v:.1f}" for v in row) for row in self.fit_percent),
                 f"MSE {self.mse:.4g}  FPE {self.fpe:.4g}  AIC {self.aic:.4g}  AICc {self.aicc:.4g}  BIC {self.bic:.4g}",
                 f"loss {self.loss_start:.4g} -> {self.loss:.4g}, data {self.n_data}, free parameters {self.n_free}, "
                 f"Jacobian condition {self.jacobian_condition:.3g}",
                 f"{'parameter':>16s} {'value':>12s} {'std':>10s} {'std/|value|':>11s} {'initial':>12s}  free"]
        for p in self.parameters:
            if p.free:
                rel = f"{p.std / abs(p.value):11.2%}" if p.value != 0 else f"{'-':>11s}"
                std, flag = f"{p.std:10.3g}", "yes" + (" (at bound)" if p.at_bound else "")
            else:
                rel, std, flag = f"{'-':>11s}", f"{'-':>10s}", "no"
            lines.append(f"{p.name:>16s} {p.value:12.5g} {std} {rel} {p.initial:12.5g}  {flag}")
        return "\n".join(lines)


@dataclass(frozen=True)
class GreyEstimate:
    model: GreyModel
    report: GreyestReport
    x0: tuple[np.ndarray, ...] = ()
    """Estimated initial states, one per time experiment (``initial_state = "estimate"``)."""

    @property
    def parameters(self) -> dict[str, float]:
        return self.model.values


# --------------------------------------------------------------------------- estimation


def frequency_selection(d: FrequencyData, weighting_filter: Any) -> tuple[np.ndarray, np.ndarray]:
    """``(mask, weights)`` of the frequency points of ``d`` after the weighting filter (see :class:`GreyestOptions`)."""
    m = np.ones(len(d.f), dtype=bool)
    w = np.ones(len(d.f)) if d.weights is None else np.asarray(d.weights, dtype=float).copy()
    wf = weighting_filter
    if isinstance(wf, str):
        if wf not in ("inv", "invsqrt"):
            raise ValueError(f"weighting_filter '{wf}': use 'inv' or 'invsqrt'")
        mag = np.abs(d.G).mean(axis=(1, 2))
        w *= 1.0 / mag if wf == "inv" else 1.0 / np.sqrt(mag)
    elif wf is not None:
        arr = np.asarray(wf, dtype=float)
        if arr.ndim == 1 and arr.size == len(d.f) and arr.size != 2:
            w *= arr
        else:
            bands = np.atleast_2d(arr)
            if bands.shape[1] != 2:
                raise ValueError("weighting_filter: bands must be [f1, f2] pairs in Hz")
            m[:] = False
            for f1, f2 in bands:
                m |= (d.f >= f1) & (d.f <= f2)
    m &= w > 0
    if not m.any():
        raise ValueError(f"data set '{d.name}': no frequency point left after the weighting filter")
    return m, w


def fit_percent(y: np.ndarray, y_model: np.ndarray) -> float:
    """NRMSE fit in per cent (as ``compare`` and ``goodnessOfFit``): 100 = perfect, 0 = as good as the mean."""
    y, y_model = np.asarray(y), np.asarray(y_model)
    den = np.linalg.norm(y - np.mean(y))
    return float(100.0 * (1.0 - np.linalg.norm(y - y_model) / den)) if den > 0 else float("nan")


@dataclass
class _Problem:
    """Residual function of one estimation, in the search coordinates ``x`` (free parameters, then initial states)."""

    model: GreyModel
    data: list
    opt: GreyestOptions
    free: list[GreyParameter]
    n_x0: int
    masks: list = field(default_factory=list)
    weights: list = field(default_factory=list)
    out_w: list = field(default_factory=list)
    nominal: np.ndarray = field(default_factory=lambda: np.zeros(0))
    L_reg: np.ndarray | None = None

    def __post_init__(self) -> None:
        for d in self.data:
            if isinstance(d, FrequencyData):
                m, w = frequency_selection(d, self.opt.weighting_filter)
                self.masks.append(m)
                self.weights.append(w)
                spread = np.std(d.G[m], axis=0).mean(axis=1) if m.any() else np.ones(d.G.shape[1])
            else:
                self.masks.append(np.ones(len(d.t), dtype=bool))
                self.weights.append(None)
                spread = np.std(d.y, axis=0)
            self.out_w.append(self._output_weight(spread))
        ow = self.opt.regularization_nominal
        if isinstance(ow, Mapping):
            nom = [ow.get(p.name, p.value) for p in self.free]
        elif ow == "zero":
            nom = [0.0 if p.coordinate == "lin" else 1.0 for p in self.free]
        else:
            nom = [p.value for p in self.free]
        self.nominal = np.array([np.log(v) if p.coordinate == "log" else v for p, v in zip(self.free, nom)])
        if self.opt.regularization_lambda > 0:
            R = np.eye(len(self.free)) if self.opt.regularization_R is None else np.asarray(self.opt.regularization_R)
            # symmetric square root: R may be only semidefinite (zero weight on some parameters)
            ev, V = np.linalg.eigh(0.5 * (R + R.T))
            self.L_reg = np.sqrt(self.opt.regularization_lambda) * (V * np.sqrt(np.clip(ev, 0, None))) @ V.T

    def _output_weight(self, spread: np.ndarray) -> np.ndarray:
        ow = self.opt.output_weight
        if ow is None:
            return np.ones_like(spread, dtype=float)
        if isinstance(ow, str):
            if ow != "normalize":
                raise ValueError(f"output_weight '{ow}': use None, 'normalize' or a vector")
            return 1.0 / np.where(spread > 0, spread, 1.0)
        return np.asarray(ow, dtype=float)

    # ---- coordinates

    def to_x(self, values: Mapping[str, float], x0s: Sequence[np.ndarray] = ()) -> np.ndarray:
        x = [np.log(values[p.name]) if p.coordinate == "log" else values[p.name] for p in self.free]
        return np.concatenate([np.array(x, dtype=float)] + [np.asarray(v, dtype=float) for v in x0s])

    def values(self, x: np.ndarray) -> dict[str, float]:
        out = {}
        for p, xi in zip(self.free, x[: len(self.free)]):
            out[p.name] = float(np.exp(xi)) if p.coordinate == "log" else float(xi)
        return out

    def x0s(self, x: np.ndarray) -> list[np.ndarray]:
        rest = x[len(self.free):]
        n = self.n_x0
        return [rest[i * n:(i + 1) * n] for i in range(len(rest) // n)] if n else []

    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        lo = [np.log(p.minimum) if p.coordinate == "log" else p.minimum for p in self.free]
        hi = [np.log(p.maximum) if p.coordinate == "log" else p.maximum for p in self.free]
        n_extra = self.n_x0 * sum(isinstance(d, TimeData) for d in self.data)
        return np.r_[lo, [-np.inf] * n_extra], np.r_[hi, [np.inf] * n_extra]

    # ---- residuals

    def predictions(self, x: np.ndarray) -> list[np.ndarray]:
        """Model output on each data set: ``(nf, ny, nu)`` complex or ``(N, ny)``."""
        vals = self.values(x)
        x0s = self.x0s(x)
        out, i_time = [], 0
        for d in self.data:
            if isinstance(d, FrequencyData):
                out.append(self.model.frequency_response(d.w, vals))
            else:
                x0 = x0s[i_time] if x0s else None
                out.append(self.model.simulate(d.t, d.u, x0, vals))
                i_time += 1
        return out

    def errors(self, x: np.ndarray) -> list[np.ndarray]:
        """Weighted errors per data set: ``(n_used, ny*nu)`` complex (frequency) or ``(N, ny)`` (time)."""
        out = []
        for d, yhat, m, w, ow in zip(self.data, self.predictions(x), self.masks, self.weights, self.out_w):
            if isinstance(d, FrequencyData):
                Gm, Gd = yhat[m], d.G[m]
                if Gm.shape != Gd.shape:
                    raise ValueError(f"model response {Gm.shape[1:]} and data {Gd.shape[1:]} have different sizes")
                e = np.log(Gm / Gd) if self.opt.frequency_error == "log" else Gm - Gd
                if self.opt.frequency_error == "log":
                    e = e.real + 1j * np.angle(Gm / Gd)
                e = e * w[m][:, None, None] * ow[None, :, None]
                out.append(e.reshape(e.shape[0], -1))
            else:
                if yhat.shape != d.y.shape:
                    raise ValueError(f"model output {yhat.shape} and data {d.y.shape} have different sizes")
                out.append((yhat - d.y) * ow[None, :])
        return out

    def penalty(self, x: np.ndarray) -> np.ndarray:
        r = []
        if self.L_reg is not None:
            r.append(self.L_reg @ (x[: len(self.free)] - self.nominal))
        if self.model.penalty_fun is not None:
            r.append(np.atleast_1d(np.asarray(self.model.penalty_fun({**self.model.values, **self.values(x)},
                                                                     *self.model.aux), dtype=float)))
        if self.opt.enforce_stability and self.model.odefun is not None:
            A = self.model.ssdata(self.values(x))[0]
            ev = np.linalg.eigvals(A)
            margin = (np.abs(ev).max() - 1.0) if self.model.Ts > 0 else ev.real.max()
            r.append(np.array([1e3 * max(margin, 0.0)]))
        return np.concatenate(r) if r else np.zeros(0)

    def residual(self, x: np.ndarray) -> np.ndarray:
        parts = []
        for e in self.errors(x):
            parts += [e.real.ravel(), e.imag.ravel()] if np.iscomplexobj(e) else [e.ravel()]
        return np.concatenate(parts + [self.penalty(x)])


def greyest(data: Data | Sequence[Data], init_sys: GreyModel, options: GreyestOptions | None = None) -> GreyEstimate:
    """Estimates the free parameters of ``init_sys`` on ``data`` (one or more experiments of the same kind).

    The initial values of the parameters are the starting point; the parameters with ``free = False`` stay fixed.
    Returns the model with the estimated values and the report (fit, information criteria, covariance).
    """
    from scipy.optimize import least_squares

    opt = options or GreyestOptions()
    data_list = list(data) if isinstance(data, (list, tuple)) else [data]
    kinds = {type(d) for d in data_list}
    if len(kinds) != 1:
        raise ValueError("all the experiments must be of the same kind (frequency or time)")
    is_time = kinds == {TimeData}
    if is_time and init_sys.odefun is None:
        raise ValueError("time data need odefun (state-space matrices)")
    free = [p for p in init_sys.parameters if p.free]
    if not free:
        raise ValueError("no free parameter")
    n_x0 = init_sys.ssdata()[0].shape[0] if is_time and opt.initial_state == "estimate" else 0
    pb = _Problem(init_sys, data_list, opt, free, n_x0)

    method = {"auto": "trf", "lsqnonlin": "trf", "lm": "lm", "dogbox": "dogbox"}[opt.search_method]
    lo, hi = pb.bounds()
    if method == "lm" and (np.isfinite(lo).any() or np.isfinite(hi).any() or opt.loss != "linear"):
        raise ValueError("search_method 'lm' does not support bounds or a robust loss")
    x_start = pb.to_x(init_sys.values, [np.zeros(n_x0)] * (len(data_list) if n_x0 else 0))
    margin = lambda b: np.where(np.isfinite(b), 1e-12 * (1 + np.abs(b)), 0.0)  # noqa: E731
    x_start = np.clip(x_start, lo + margin(lo), hi - margin(hi))
    if n_x0:  # initial state of each experiment by linear least squares at the starting parameters
        x_start = _initial_states(pb, x_start)
    r0 = pb.residual(x_start)
    kw: dict[str, Any] = dict(method=method, max_nfev=opt.max_iterations, ftol=opt.tolerance, xtol=opt.tolerance,
                              gtol=opt.tolerance, diff_step=opt.diff_step, x_scale="jac", verbose=2 if opt.display else 0)
    if method != "lm":
        kw.update(bounds=(lo, hi), loss=opt.loss, f_scale=opt.f_scale)
    sol = least_squares(pb.residual, x_start, **kw)
    x = sol.x
    est = init_sys.with_values(pb.values(x))
    report = _report(pb, x, x_start, r0, sol, method, lo, hi)
    return GreyEstimate(model=est, report=report, x0=tuple(pb.x0s(x)))


def _initial_states(pb: _Problem, x: np.ndarray) -> np.ndarray:
    """Least-squares initial states (the output is linear in ``x0``) at the parameters of ``x``."""
    vals = pb.values(x)
    A, B, C, D = pb.model.ssdata(vals)
    n = A.shape[0]
    x = x.copy()
    i = 0
    for d, ow in zip(pb.data, pb.out_w):
        dt = _uniform_step(d.t)
        Ad = A if pb.model.Ts > 0 else zoh(A, B, dt)[0]
        y0 = pb.model.simulate(d.t, d.u, None, vals)
        # free response C Ad^k e_j
        Phi = np.empty((len(d.t), C.shape[0], n))
        M = np.eye(n)
        for k in range(len(d.t)):
            Phi[k] = C @ M
            M = Ad @ M
        Phi = Phi * ow[None, :, None]
        x0 = np.linalg.lstsq(Phi.reshape(-1, n), ((d.y - y0) * ow[None, :]).ravel(), rcond=None)[0]
        x[len(pb.free) + i * n: len(pb.free) + (i + 1) * n] = x0
        i += 1
    return x


def _report(pb: _Problem, x: np.ndarray, x_start: np.ndarray, r0: np.ndarray, sol, method: str,
            lo: np.ndarray, hi: np.ndarray) -> GreyestReport:
    d_free = len(pb.free)
    errs = pb.errors(x)
    # fit per experiment and output, on the used points, without weights
    fits = []
    for dset, yhat, m in zip(pb.data, pb.predictions(x), pb.masks):
        if isinstance(dset, FrequencyData):
            Gm, Gd = yhat[m], dset.G[m]
            fits.append([fit_percent(Gd[:, i, :].ravel(), Gm[:, i, :].ravel()) for i in range(Gd.shape[1])])
        else:
            fits.append([fit_percent(dset.y[:, i], yhat[:, i]) for i in range(dset.y.shape[1])])
    E = np.vstack(errs)
    N, ny = E.shape
    V = np.linalg.det((E.conj().T @ E).real / N) if ny > 1 else float(np.mean(np.abs(E) ** 2))
    V = max(float(V), np.finfo(float).tiny)
    d = len(x)  # free parameters and estimated initial states
    mse = float(np.mean(np.abs(E) ** 2) * ny)
    fpe = V * (1 + d / N) / (1 - d / N) if N > d else float("inf")
    const = N * ny * (np.log(2 * np.pi) + 1)
    aic = N * np.log(V) + 2 * d + const
    aicc = aic + 2 * d * (d + 1) / (N - d - 1) if N - d - 1 > 0 else float("inf")
    bic = N * np.log(V) + d * np.log(N) + const

    # covariance from the Jacobian of the plain residuals (no robust loss), central differences
    r = pb.residual(x)
    n_data = r.size - pb.penalty(x).size
    J = np.empty((r.size, len(x)))
    for j in range(len(x)):
        h = pb.opt.diff_step * max(1.0, abs(x[j]))
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        J[:, j] = (pb.residual(xp) - pb.residual(xm)) / (2 * h)
    dof = max(n_data - len(x), 1)
    s2 = float(np.sum(r[:n_data] ** 2)) / dof
    Jf = J[:, :d_free]
    col = np.linalg.norm(Jf, axis=0)
    col[col == 0] = 1.0
    sv = np.linalg.svd(Jf / col, compute_uv=False)
    cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf")
    bread = np.linalg.pinv(J.T @ J)
    if pb.opt.covariance == "classic":
        cov_x = s2 * bread
    else:
        Jd, rd = J[:n_data], r[:n_data]
        cov_x = bread @ ((Jd * rd[:, None] ** 2).T @ Jd) @ bread * n_data / dof
    cov_x = cov_x[:d_free, :d_free]
    # from the search coordinates to the parameters: d theta / d x = theta for log parameters
    jac = np.array([np.exp(xi) if p.coordinate == "log" else 1.0 for p, xi in zip(pb.free, x[:d_free])])
    cov = cov_x * np.outer(jac, jac)
    std = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(std, std)
    vals = pb.values(x)
    tol = 1e-6
    params = []
    for p in pb.model.parameters:
        if p.free:
            i = [q.name for q in pb.free].index(p.name)
            at = bool(x[i] - lo[i] <= tol * (1 + abs(lo[i])) or hi[i] - x[i] <= tol * (1 + abs(hi[i])))
            params.append(ParameterEstimate(p.name, vals[p.name], p.value, True,
                                            float(std[i]) if np.isfinite(cond) else float("nan"), at))
        else:
            params.append(ParameterEstimate(p.name, p.value, p.value, False, float("nan"), False))
    status = "converged" if sol.status > 0 else "maximum evaluations reached" if sol.status == 0 else "failed"
    return GreyestReport(
        method=method, status=status, termination=str(sol.message), n_evaluations=int(sol.nfev),
        fit_percent=np.array(fits), mse=mse, fpe=float(fpe), aic=float(aic), aicc=float(aicc), bic=float(bic),
        loss_start=float(np.sum(r0**2)), loss=float(np.sum(r**2)), n_data=int(N), n_free=d_free,
        parameters=tuple(params), covariance=cov, correlation=corr, free_names=tuple(p.name for p in pb.free),
        jacobian_condition=cond, residuals=r)


def compare(data: Data | Sequence[Data], model: GreyModel, x0: Sequence[np.ndarray] | None = None,
            options: GreyestOptions | None = None) -> tuple[list[np.ndarray], np.ndarray]:
    """Model response on each data set and NRMSE fit ``(n_experiments, ny)`` (``compare``), on the points selected
    by the weighting filter of ``options`` (all points by default)."""
    data_list = list(data) if isinstance(data, (list, tuple)) else [data]
    opt = options or GreyestOptions()
    out, fits = [], []
    for k, d in enumerate(data_list):
        if isinstance(d, FrequencyData):
            m = frequency_selection(d, opt.weighting_filter)[0]
            G = model.frequency_response(d.w)
            out.append(G)
            fits.append([fit_percent(d.G[m][:, i, :].ravel(), G[m][:, i, :].ravel()) for i in range(d.G.shape[1])])
        else:
            y = model.simulate(d.t, d.u, None if x0 is None else x0[k])
            out.append(y)
            fits.append([fit_percent(d.y[:, i], y[:, i]) for i in range(d.y.shape[1])])
    return out, np.array(fits)
