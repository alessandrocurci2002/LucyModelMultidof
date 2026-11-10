"""LQG con modello interno senza cuscinetti, a parita' di errore d'inseguimento con la cascata, sul plant Dahl.

Senza disturbi l'LQG dovrebbe dare prestazioni almeno pari a quelle della cascata. q_p si cerca per bisezione finche' l'errore rms mediano sui
doppietti e' quello della cascata, che non ha integrale di posizione; q_i quasi nullo la riproduce, q_i = 0.1 mostra l'effetto dell'integrale.

    .venv/bin/python -I analisi_lqg_equivalenza_cascata.py

Risultati in results_lqg_tuning/comparison_measured/cascade_equivalence.json.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]

import analisi_lqg_confronto_misure as cm  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_passi_2_4 as p24  # noqa: E402
from lucy_cl import lqg  # noqa: E402

Q_I_MIN = 1e-6
PROVE = ((Q_I_MIN, 1.0), (Q_I_MIN, 1e4), (0.1, 1.0))
"""q_i, rho."""


def main() -> None:
    par, Kt, S, T, mot, assi = g.carica()
    out = {}
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        W, Kpp = p24.comandi_tipici(asse), a.drive.p.Kpp
        sintesi = lambda qp, qi, rho=1.0: lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, qi, qp, np.zeros(len(a.Ap) - 1), rho_e=rho, v_continua=True)
        casc = cm.prova(a, None, W)
        r = {"K_pp": Kpp, "cascade": casc, "P_encoder": cm.prova(a, cm.PEncoder(Kpp), W), "same_k_p": {}, "same_error": []}
        for qi in (0.1, Q_I_MIN):
            c = sintesi(Kpp ** 2, qi)
            r["same_k_p"][str(qi)] = {"k_p": c.k_p, "k_i": c.k_i} | cm.prova(a, c, W)
        for qi, rho in PROVE:
            lo, hi = Kpp ** 2, 16 * Kpp ** 2
            for _ in range(6):
                m = float(np.sqrt(lo * hi))
                if cm.prova(a, sintesi(m, qi, rho), W)["err_rms_median_mdeg"] > casc["err_rms_median_mdeg"]:
                    lo = m
                else:
                    hi = m
            c = sintesi(hi, qi, rho)
            r["same_error"].append({"q_i": qi, "rho": rho, "q_p": hi, "k_p": c.k_p, "k_i": c.k_i, "lqg": cm.prova(a, c, W), "without_K_f": cm.prova(a, dataclasses.replace(c, K_f=0 * c.K_f), W)})
            print(ch, qi, rho, round(hi, 1), flush=True)
        out[ch] = r
    (cm.CARTELLA / "cascade_equivalence.json").write_text(json.dumps(p24.senza_array(out), indent=1))


if __name__ == "__main__":
    main()
