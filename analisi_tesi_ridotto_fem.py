"""Figura della tesi (cap. Multibody, sottosezione "Reduced model - 18-element FEM Comparison"): risposta in frequenza del
modello ridotto a 7 gradi di liberta' e del FEM a 18 elementi da cui deriva.

Richiesta dell'utente (2026-11-05): grafici di confronto fra modello ridotto e FEM a 18 elementi, senza numeri sull'asse delle
frequenze (policy aziendale). Stessi modelli della Tab. ``tab:mb_lump_fem7``, cioe' della configurazione di base ``config.yaml``
(``development_report/build_development_report.py``, tag "7"): theta_0 = 0, assi accoppiati, con e senza telescopio.
Grandezza: velocita' dell'encoder per ampere, (rad/s)/A, come le misure (eq. ``eq:mb_plant``), canale dell'encoder, che e' quello
misurato. Stile delle figure della tesi (``style.py``, come ``analisi_tesi_hil.py``).

    .venv/bin/python -I analisi_tesi_ridotto_fem.py

Figura in Thesis-Writing/plots/multibody/reduced_vs_fem.pdf.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE / "src"), str(RADICE / "development_report")]

import build_development_report as br  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402

STILE = Path("/home/alessandro/Desktop/Thesis/OpenVins/OpenVincenzo/colcon_ws/src/open_vins/ov_data/sc_cantina01_thesis/output/plots/scripts")
FIGURA = Path("/home/alessandro/Desktop/Thesis/Thesis-Writing/plots/multibody")
F = np.geomspace(br.F_MIN, br.F_MAX, 4000)
MODELLI = {"fem": dict(label="18-element FEM", color="#000000", ls="-", lw=1.2),
           "red": dict(label="Reduced 7-DOF model", color="#0072B2", ls="--", lw=1.4)}
CASI = (("azimut", True, "(a) Azimuth"), ("polare", True, "(b) Polar axis, with telescope"), ("polare", False, "(c) Polar axis, without telescope"))


def main() -> None:
    sys.path.insert(0, str(STILE))
    import style  # noqa: E402
    style.setup()
    par, sis = br.build_models(br.CFG_7)
    fem_t, fem_n = br.full_fem(par)
    Kt = par.sperimentale.Kt
    fig, ax = style.figure(aspect=1.05, nrows=3, ncols=2, sharex=True)
    for r, (asse, tel, titolo) in enumerate(CASI):
        ce = spm.canali(asse)[0]
        G = {"fem": spm.frf_velocita((fem_t if tel else fem_n).lineare_completo(0.0, 0.0), asse, F, Kt)[ce],
             "red": spm.frf_velocita(lz.linearizza(sis[tel].accoppiato, 0.0, 0.0), asse, F, Kt)[ce]}
        for k, s in MODELLI.items():
            ax[r, 0].semilogx(F, 20 * np.log10(np.abs(G[k])), **s)
            ax[r, 1].semilogx(F, np.rad2deg(np.angle(G[k])), **{q: v for q, v in s.items() if q != "label"})
        ax[r, 0].set_title(titolo, loc="left", fontsize=style.FONT_SIZE)
        ax[r, 0].set_ylabel("Magnitude [dB (rad/s)/A]")
        ax[r, 1].set_ylabel("Phase [deg]")
        ax[r, 1].set_ylim(-180, 180)
        ax[r, 1].set_yticks(np.arange(-180, 181, 90))
    for x in ax.ravel():
        x.set_xlim(F[0], F[-1])
        # frequenze non divulgabili: griglia logaritmica senza numeri
        x.tick_params(axis="x", which="both", labelbottom=False)
        x.grid(True, which="minor", axis="x", color="#EEEEEE", lw=0.4)
    for x in ax[-1]:
        x.set_xlabel("Frequency (logarithmic scale)")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=True, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    FIGURA.mkdir(parents=True, exist_ok=True)
    style.save(fig, FIGURA / "reduced_vs_fem.pdf")


if __name__ == "__main__":
    main()
