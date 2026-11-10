import json
from pathlib import Path

import numpy as np
import pytest

from lucy_sim import GimbalModel

REF = Path(__file__).parent / "reference"


@pytest.fixture(scope="session")
def model() -> GimbalModel:
    return GimbalModel.load()


@pytest.fixture(scope="session")
def ref_model() -> dict:
    return json.loads((REF / "model.json").read_text())


def rel(a, b) -> float:
    """Largest difference relative to the largest magnitude of the reference."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-300))
