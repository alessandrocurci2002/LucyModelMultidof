"""Fixture condivise dalla suite di test.

I test verificano il **codice**, non i dati di ``config.yaml``: le proprieta'
controllate (simmetrie, conservazioni, forme chiuse, limiti) valgono qualunque
siano i dati. Per questo la suite usa una configurazione di riferimento
congelata, ``tests/config_riferimento.yaml``, indipendente dal file del
progetto che l'utente modifica; dove serve un caso controllato se ne modifica
una copia con :func:`config_dati`.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import pytest

from gimbal_mdof import params as pr

RADICE = Path(__file__).resolve().parents[1]
CONFIG = Path(__file__).resolve().parent / "config_riferimento.yaml"

_DATI_BASE = pr.carica_config(CONFIG)


def config_dati() -> dict[str, Any]:
    """Copia profonda del dizionario di ``config.yaml``, modificabile nei test."""
    return copy.deepcopy(_DATI_BASE)


def parametri(dati: dict[str, Any]) -> pr.Parametri:
    """Parametri da un dizionario, con i percorsi relativi risolti come il file."""
    return pr.parametri_da_config(CONFIG, dati=dati)


@pytest.fixture(scope="session")
def par() -> pr.Parametri:
    """Configurazione di riferimento."""
    return pr.parametri_da_config(CONFIG)
