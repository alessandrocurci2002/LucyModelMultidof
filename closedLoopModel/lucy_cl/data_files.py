"""Readers for the drive data: WorkBench parameter files (.akd) and WorkBench frequency-response exports (.csv)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


def read_akd(path: str | Path) -> dict[str, str]:
    """Parameters of a WorkBench ``.akd`` file (lines ``NAME value``; header lines start with ``#``)."""
    out: dict[str, str] = {}
    for line in Path(path).read_text(errors="ignore").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        out[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return out


@dataclass(frozen=True)
class WorkBenchExport:
    """A WorkBench frequency-response export: traces (magnitude dB, phase deg) and metadata."""

    path: Path
    f: np.ndarray
    traces: dict[str, tuple[np.ndarray, np.ndarray]]
    meta: dict[str, str]

    def response(self, trace: str) -> np.ndarray:
        """Complex response of a trace (``Plant``, ``Controller``, ``Open Loop``, ``Closed Loop``,
        ``Feedback Biquad``)."""
        g, p = self.traces[trace]
        return 10 ** (g / 20) * np.exp(1j * np.deg2rad(p))

    @property
    def excitation_A(self) -> float:
        return float(self.meta.get("MeasurementConfiguration/CurrentExcitation", "nan"))

    def tuning(self, key: str) -> float:
        """A value of the ``TuningConfiguration`` block (gains active during the test)."""
        return float(self.meta[f"TuningConfiguration/{key}"])

    def filters(self) -> list[tuple[str, float, float, float, float]]:
        """Velocity-loop biquads recorded during the test: (group, fz, qz, fp, qp), forward group first."""
        out = []
        for k in range(8):
            b = f"AntiResonanceFilters[{k}]/"
            if b + "NumeratorF" not in self.meta:
                break
            out.append((self.meta[b + "ArfGroup"], float(self.meta[b + "NumeratorF"]),
                        float(self.meta[b + "NumeratorQ"]), float(self.meta[b + "DenominatorF"]),
                        float(self.meta[b + "DenominatorQ"])))
        return out

    def amplitude_rad(self, f_ref: float) -> float:
        """Amplitude of the axis rotation at ``f_ref`` during the test, :math:`|G(f)|\\,I / (2\\pi f)`, with
        :math:`G` the measured plant (rad/s)/A and :math:`I` the current excitation."""
        G = self.response("Plant")
        g = np.interp(np.log(f_ref), np.log(self.f), np.abs(G))
        return float(g * self.excitation_A / (2 * np.pi * f_ref))


def read_workbench(path: str | Path) -> WorkBenchExport:
    """Reads a WorkBench Bode export with the parser of ``gimbal_mdof.sperimentale`` (columns
    ``<trace> Gain[dB]``, ``<trace> Phase[deg]`` and the ``#ParentName, #Name, #Type, #Value`` metadata)."""
    from gimbal_mdof import sperimentale as spm

    f, traces, meta = spm.carica_bode(path)
    return WorkBenchExport(path=Path(path), f=f, traces=traces, meta=meta)
