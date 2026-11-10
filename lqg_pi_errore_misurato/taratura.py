"""Taratura dell'LQG-PI con il PI sull'errore misurato (Fig. 8.4b), con la stessa procedura e gli stessi vincoli della tesi.

Si esegue ``analisi_lqg_passi_2_4.main`` (passi 2-4 di Gawronski, taratura fine, rialzo di q_p, prova su V; modello interno senza
cuscinetti) con il controllore di ``controllore.py``. Vincoli invariati: corrente sotto 18 A sui doppietti (plant Dahl),
min |1 + L| >= 0.5, anello stabile con i cuscinetti in strisciamento e bloccati, banda sotto la fondamentale.

In piu', richiesta dell'utente (2026-11-05): tracking system di bassa autorita' (Gawronski, par. 8.1.4), definita sui pesi,
Q_t >> Q_f. Scelta mia per "molto maggiore": un ordine di grandezza, q_p / q_f >= ``RAPPORTO_MIN`` per ogni modo il cui peso
viene alzato dal passo 4 (il peso di base 0.1 degli stati flessibili resta quello del libro).

    .venv/bin/python -I lqg_pi_errore_misurato/taratura.py

Risultati in lqg_pi_errore_misurato/results/steps_2_4/ (summary.json e figure); la configurazione e' copiata in
config/lqg_tuning.json.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

QUI = Path(__file__).resolve().parent
sys.argv = [sys.argv[0], "--senza-cuscinetti"]
sys.path.insert(0, str(QUI))

import controllore  # noqa: E402

controllore.attiva()

import analisi_lqg_passi_2_4 as p24  # noqa: E402

RAPPORTO_MIN = 10.0
_LIMITI = p24.limiti


def limiti(x: dict, f_fond: float, integrale: bool = False) -> list[str]:
    v = _LIMITI(x, f_fond, integrale)
    if x.get("q_f") and x["q_p"] / max(x["q_f"].values()) < RAPPORTO_MIN:
        v.append("low authority")
    return v


def main() -> None:
    p24.CARTELLA = QUI / "results" / "steps_2_4"
    p24.limiti = limiti
    p24.main()
    (QUI / "config").mkdir(exist_ok=True)
    shutil.copy(p24.CARTELLA / "summary.json", QUI / "config" / "lqg_tuning.json")


if __name__ == "__main__":
    main()
