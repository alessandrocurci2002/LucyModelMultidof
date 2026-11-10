"""Figure della tesi (cap. Multibody, sezione "Experimental Frequency Responses"): risposta misurata e modello a 7 gradi di liberta'.

Richiesta dell'utente (2026-11-05): le figure ``cmp7_az``, ``cmp7_pol`` e ``cmp7_nt`` (Fig. ``fig:mb_cmp7``) senza i numeri
sull'asse delle frequenze (policy aziendale), come ``analisi_tesi_ridotto_fem.py``. Stessi dati e stesso modello di
``development_report/build_development_report.py`` (``bode_comparison``, configurazione di base ``config.yaml``): misura (nera)
e modello a 7 gradi di liberta' al canale dell'encoder, che e' quello misurato; area in grigio: fuori dalla banda misurata.
Stile delle figure della tesi (``style.py``).

    .venv/bin/python -I analisi_tesi_cmp7.py

Figure in Thesis-Writing/plots/multibody/cmp7_{az,pol,nt}.png.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE / "src"), str(RADICE / "development_report")]

import build_development_report as br  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402

STILE = Path("/home/alessandro/Desktop/Thesis/OpenVins/OpenVincenzo/colcon_ws/src/open_vins/ov_data/sc_cantina01_thesis/output/plots/scripts")
FIGURA = Path("/home/alessandro/Desktop/Thesis/Thesis-Writing/plots/multibody")
SEL = {"az": "Azimuth_BodeMeasurement_Po_0", "pol": "Elevation_BodeMeasurement_bode_New_16_10_26", "nt": "BodeMeasurement_NoTelescope"}


def main() -> None:
    sys.path.insert(0, str(STILE))
    import style  # noqa: E402
    style.setup()
    par, sis = br.build_models(br.CFG_7)
    cfg = par.sperimentale
    misure = {Path(v.file).stem: spm.carica_misura(cfg, v) for v in cfg.misure}
    for k, stem in SEL.items():
        mis = misure[stem]
        conf = spm.confronta(par, mis, sis, br.F_GRID)
        ce = spm.canali(mis.asse)[0]
        G = conf.G_modello[ce]
        fig, ax = style.figure(aspect=0.62, nrows=2, ncols=1, sharex=True)
        ax[0].semilogx(conf.f_modello, 20 * np.log10(np.abs(G)), color="#0072B2", ls="-", lw=1.4, label="7-DOF model")
        ax[1].semilogx(conf.f_modello, np.rad2deg(np.angle(G)), color="#0072B2", ls="-", lw=1.4)
        ax[0].semilogx(mis.f, mis.modulo_db, color="#000000", marker="o", ms=2.2, lw=0.8, label="Measurement")
        ax[1].semilogx(mis.f, mis.fase_deg, color="#000000", marker="o", ms=2.2, lw=0.8)
        for x in ax:
            x.set_xlim(br.F_MIN, br.F_MAX)
            x.axvspan(br.F_MIN, mis.f[0], color="#EEEEEE", zorder=0)
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
        style.save(fig, FIGURA / f"cmp7_{k}.png")


if __name__ == "__main__":
    main()
