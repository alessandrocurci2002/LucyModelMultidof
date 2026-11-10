"""Simulazioni della tesi (circonferenza, cap. 6, e traiettorie registrate del drone, cap. 8) con l'LQG-PI sull'errore misurato.

Stesse simulazioni, scenari, semi e metriche di ``analisi_tesi_hil.py``; cambiano solo il controllore (``controllore.py``) e il
tuning (``config/lqg_tuning.json``, da ``taratura.py``). Risultati, tabelle e figure restano in questa cartella: la tesi non cambia.

    LAVORI=4 .venv/bin/python -I lqg_pi_errore_misurato/simulazioni.py [--figure]
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

QUI = Path(__file__).resolve().parent
SOLO_FIGURE = "--figure" in sys.argv
sys.argv = [sys.argv[0]] + (["--figure"] if SOLO_FIGURE else [])
sys.path.insert(0, str(QUI))
os.environ.setdefault("LAVORI", "4")

import controllore  # noqa: E402

controllore.attiva()

import analisi_lqg_circonferenza as cc  # noqa: E402

cc.TUNING = QUI / "config" / "lqg_tuning.json"

import analisi_tesi_hil as th  # noqa: E402  (importa analisi_tesi_puntamento, che costruisce assi e controllori con il tuning sopra)

th.CARTELLA = QUI / "results" / "hil"
th.FIGURE = QUI / "figures" / "hil"
th.FIGURE_LQG = QUI / "figures" / "lqg"

if __name__ == "__main__":
    th.main()
