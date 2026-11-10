"""Figura della tesi (cap. Multibody, sezione dell'anti-risonanza del polare): risposta misurata del polare con il primo montaggio
del telescopio e dopo l'inserimento delle viti mancanti del collegamento anello-telescopio.

Richiesta dell'utente (2026-11-05): confronto fra la misura polare usata fino a quel punto del capitolo (la stessa di
``cmp7_pol``, ``Elevation_BodeMeasurement_bode_New_16_10_26``, 4 viti, 3.5 A) e quella con 14 viti
(``VClosedLoop/VTesi/BodeMeasurement_Polar_14viti.csv``, 23/10/2026, 3 A), usata poi per l'identificazione. Asse delle frequenze
senza numeri (policy aziendale). Stile delle figure della tesi (``style.py``).

    .venv/bin/python -I analisi_tesi_viti.py

Figura in Thesis-Writing/plots/multibody/screws_14_polar.png.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE / "src")]

from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402

STILE = Path("/home/alessandro/Desktop/Thesis/OpenVins/OpenVincenzo/colcon_ws/src/open_vins/ov_data/sc_cantina01_thesis/output/plots/scripts")
FIGURA = Path("/home/alessandro/Desktop/Thesis/Thesis-Writing/plots/multibody")
MISURE = ((RADICE / "config.yaml", "Elevation_BodeMeasurement_bode_New_16_10_26.csv",
           dict(label="First assembly, 4 screws", color="#D55E00", ls="--", lw=1.2)),
          (RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml", "VClosedLoop/VTesi/BodeMeasurement_Polar_14viti.csv",
           dict(label="Complete assembly, 14 screws", color="#000000", ls="-", lw=1.2)))


def misura(cfg_path: Path, nome: str):
    cfg = pr.parametri_da_config(cfg_path).sperimentale
    return spm.carica_misura(cfg, next(v for v in cfg.misure if v.file == nome))


def main() -> None:
    sys.path.insert(0, str(STILE))
    import style  # noqa: E402
    style.setup()
    fig, ax = style.figure(aspect=0.62, nrows=2, ncols=1, sharex=True)
    f_min, f_max = np.inf, 0.0
    for cfg_path, nome, s in MISURE:
        m = misura(cfg_path, nome)
        ax[0].semilogx(m.f, m.modulo_db, **s)
        ax[1].semilogx(m.f, m.fase_deg, **{k: v for k, v in s.items() if k != "label"})
        f_min, f_max = min(f_min, m.f[0]), max(f_max, m.f[-1])
    for x in ax:
        x.set_xlim(f_min, f_max)
        # frequenze non divulgabili: griglia logaritmica senza numeri
        x.tick_params(axis="x", which="both", labelbottom=False)
        x.grid(True, which="minor", axis="x", color="#EEEEEE", lw=0.4)
    ax[0].set_ylabel("Magnitude [dB (rad/s)/A]")
    ax[1].set_ylabel("Phase [deg]")
    ax[1].set_ylim(-180, 180)
    ax[1].set_yticks(np.arange(-180, 181, 90))
    ax[1].set_xlabel("Frequency (logarithmic scale)")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=True, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    FIGURA.mkdir(parents=True, exist_ok=True)
    style.save(fig, FIGURA / "screws_14_polar.png")


if __name__ == "__main__":
    main()
