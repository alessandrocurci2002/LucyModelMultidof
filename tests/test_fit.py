"""Fit NRMSE come in System Identification Toolbox: 100 (1 - ||y - y_hat|| / ||y - mean(y)||)."""

from __future__ import annotations

import numpy as np
import pytest

from gimbal_mdof import sperimentale as spm


def test_fit_nrmse_casi_limite() -> None:
    rng = np.random.default_rng(0)
    y = rng.normal(size=50) + 1j * rng.normal(size=50)
    assert spm.fit_nrmse(y, y) == pytest.approx(100.0)
    assert spm.fit_nrmse(y, np.full_like(y, y.mean())) == pytest.approx(0.0, abs=1e-12)
    assert spm.fit_nrmse(y, -y) < 0.0
    # coerente con la definizione esplicita
    yh = y + 0.1 * rng.normal(size=50)
    atteso = 100 * (1 - np.linalg.norm(y - yh) / np.linalg.norm(y - y.mean()))
    assert spm.fit_nrmse(y, yh) == pytest.approx(atteso)


def test_fit_bande() -> None:
    f = np.logspace(0, 3, 60)
    y = 1.0 / (1j * f)
    fb = spm.fit_bande(f, y, y)
    assert set(fb) == {b for b, _, _ in spm.BANDE_FIT}
    assert all(v == pytest.approx(100.0) for v in fb.values())
