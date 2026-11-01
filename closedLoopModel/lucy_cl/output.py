"""Figures and tables of ``main.py`` (English labels, PDF and PNG)."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .closed_loop import LoopFRF, SimResult  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "legend.fontsize": 8,
                     "figure.dpi": 110, "savefig.bbox": "tight"})

MDEG = 1e3 * 180.0 / np.pi
"""rad -> millidegrees."""


def save(fig: plt.Figure, folder: Path, name: str) -> list[Path]:
    folder.mkdir(parents=True, exist_ok=True)
    out = [folder / f"{name}.pdf", folder / f"{name}.png"]
    for p in out:
        fig.savefig(p)
    plt.close(fig)
    return out


def _db(x: np.ndarray) -> np.ndarray:
    return 20 * np.log10(np.abs(x))


def _deg(x: np.ndarray) -> np.ndarray:
    return np.rad2deg(np.angle(x))


def bode(curves: Sequence[tuple[str, np.ndarray, np.ndarray, dict]], title: str, ylabel: str,
         zero_db: bool = False) -> plt.Figure:
    """Magnitude and phase of several responses: (label, f, H, plot kwargs)."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.6, 5.0), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    for lab, f, H, kw in curves:
        a1.semilogx(f, _db(H), label=lab, **kw)
        a2.semilogx(f, _deg(H), **kw)
    if zero_db:
        a1.axhline(0.0, color="k", lw=0.6)
    a1.set_ylabel(f"|{ylabel}| [dB]")
    a2.set_ylabel("phase [deg]")
    a2.set_yticks(np.arange(-180, 181, 90))
    a2.set_xlabel("frequency [Hz]")
    a1.set_title(title)
    a1.legend(loc="best")
    fig.align_ylabels()
    return fig


def loop_comparison(axis: str, f_meas: np.ndarray, meas: Mapping[str, np.ndarray], model: LoopFRF,
                    model_akd: LoopFRF | None, label_meas: str) -> plt.Figure:
    """2 x 3 panels: Plant, Open Loop, Closed Loop (magnitude over phase), measurement against model."""
    fig, ax = plt.subplots(2, 3, figsize=(11.5, 5.4), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    traces = (("Plant", model.plant, None, "(rad/s)/A"), ("Open Loop", model.open_loop,
              None if model_akd is None else model_akd.open_loop, "-"),
              ("Closed Loop", model.closed_loop, None if model_akd is None else model_akd.closed_loop, "-"))
    for j, (name, H, Ha, unit) in enumerate(traces):
        a1, a2 = ax[0, j], ax[1, j]
        a1.semilogx(f_meas, _db(meas[name]), "k.", ms=3, label=label_meas)
        a2.semilogx(f_meas, _deg(meas[name]), "k.", ms=3)
        a1.semilogx(model.f, _db(H), "C0", lw=1.2, label="model, tuning of the test")
        a2.semilogx(model.f, _deg(H), "C0", lw=1.2)
        if Ha is not None:
            a1.semilogx(model_akd.f, _db(Ha), "C3--", lw=1.0, label="model, .akd parameters")
            a2.semilogx(model_akd.f, _deg(Ha), "C3--", lw=1.0)
        if name == "Open Loop":
            a1.axhline(0.0, color="k", lw=0.6)
        a1.set_title(f"{axis}: {name}" + (f" [{unit}]" if unit != "-" else ""))
        a2.set_yticks(np.arange(-180, 181, 90))
        a2.set_xlabel("frequency [Hz]")
    ax[0, 0].set_ylabel("magnitude [dB]")
    ax[1, 0].set_ylabel("phase [deg]")
    ax[0, 1].legend(loc="lower left")
    fig.tight_layout()
    return fig


def zeta_sweep(f: np.ndarray, loops: Sequence[tuple[float, np.ndarray]], axis: str, title: str) -> plt.Figure:
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.6, 5.0), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    for k, (z, L) in enumerate(loops):
        a1.semilogx(f, _db(L), color=plt.cm.viridis(k / max(1, len(loops) - 1)), label=f"zeta_1 = {z:g}")
        a2.semilogx(f, np.rad2deg(np.unwrap(np.angle(L))), color=plt.cm.viridis(k / max(1, len(loops) - 1)))
    a1.axhline(0.0, color="k", lw=0.6)
    for y in (-180, -540):
        a2.axhline(y, color="k", lw=0.6, ls=":")
    a1.set_ylabel("|L| [dB]")
    a2.set_ylabel("phase (unwrapped) [deg]")
    a2.set_xlabel("frequency [Hz]")
    a1.set_title(title)
    a1.legend()
    return fig


def dahl_loops(elements: Sequence, n: int = 4000) -> plt.Figure:
    """Steady Dahl loop at the test amplitude against the ellipse of the equivalent k(1 + j eta)."""
    fig, axs = plt.subplots(1, len(elements), figsize=(4.0 * len(elements), 3.3), squeeze=False)
    for a, el in zip(axs[0], elements):
        A, Tc, xi = el.amplitude, el.Tc, el.sigma0 * el.amplitude / el.Tc
        phi = np.linspace(-np.pi / 2, 3 * np.pi / 2, n, endpoint=False)
        u = np.sin(phi)
        t1 = np.tanh(xi)
        rising = np.cos(phi) >= 0
        T = Tc * np.where(rising, 1 - (1 + t1) * np.exp(-xi * (u + 1)), -(1 - (1 + t1) * np.exp(-xi * (1 - u))))
        a.plot(1e3 * A * u, T, "C0", label="Dahl element")
        a.plot(1e3 * A * u, el.k_eq * A * (np.sin(phi) + el.eta_eq * np.cos(phi)), "C1--",
               label=f"k(1 + j eta), k = {el.k_eq:g}, eta = {el.eta_eq:g}")
        a.set_xlabel("bearing rotation [mrad]")
        a.set_ylabel("friction torque [N m]")
        a.set_title(f"{el.axis}: A = {1e3 * A:.2f} mrad, sigma0 = {el.sigma0:.0f} N m/rad, Tc = {Tc:.2f} N m",
                    fontsize=8)
        a.legend(fontsize=7)
    fig.tight_layout()
    return fig


def scenario(r: SimResult, axes: Sequence[str], linear: Mapping[str, np.ndarray] | None = None) -> plt.Figure:
    """Time histories of the commanded axes: position, error, current, torques and phase voltage."""
    fig, ax = plt.subplots(4, len(axes), figsize=(5.4 * len(axes), 8.4), sharex=True, squeeze=False)
    for j, a in enumerate(axes):
        s, t = r.signals[a], r.t
        load = "fork" if a == "azimuth" else "telescope"
        ax[0, j].plot(t, MDEG * s["p_cmd"], "k--", lw=0.9, label="command")
        ax[0, j].plot(t, MDEG * s["p_enc"], "C0", label="encoder")
        ax[0, j].plot(t, MDEG * s["p_load"], "C2", lw=0.9, label=f"{load} node")
        if linear is not None and a in linear:
            ax[0, j].plot(t, MDEG * linear[a], "C3:", lw=1.4, label="linearised model (Dahl -> sigma0)")
        ax[0, j].set_ylabel("position [mdeg]")
        ax[0, j].set_title(f"{r.scenario.name}: {a}")
        ax[0, j].legend(fontsize=7)
        ax[1, j].plot(t, MDEG * s["error"], "C0", label="command - encoder")
        ax[1, j].plot(t, MDEG * (s["p_cmd"] - s["p_load"]), "C2", lw=0.9, label=f"command - {load}")
        ax[1, j].set_ylabel("error [mdeg]")
        ax[1, j].legend(fontsize=7)
        ax[2, j].plot(t, s["i_cmd"], "C1", label="current command")
        ax[2, j].plot(t, s["i"], "C0", lw=0.9, label="motor current")
        ax[2, j].set_ylabel("current [A rms]")
        ax[2, j].legend(fontsize=7)
        ax[3, j].plot(t, s["tau_m"], "C0", label="motor torque")
        ax[3, j].plot(t, s["T_f"], "C3", label="bearing friction (Dahl)")
        ax[3, j].set_ylabel("torque [N m]")
        ax[3, j].set_xlabel("time [s]")
        ax[3, j].legend(fontsize=7)
    fig.tight_layout()
    return fig


def markdown_table(header: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(out)


def fmt(x: float, nd: int = 3) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "-"
    if isinstance(x, float) and np.isinf(x):
        return "inf"
    return f"{x:.{nd}g}"
