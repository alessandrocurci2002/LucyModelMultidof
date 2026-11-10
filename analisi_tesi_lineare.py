"""Figura della tesi (cap. LQG, sezione del modello lineare): pianta corrente -> velocita' dell'encoder del modello senza
cuscinetti, quello interno dell'LQG-PI, contro una misura per asse.

Richiesta dell'utente (2026-11-03): un solo confronto sperimentale per asse; asse delle frequenze senza numeri (policy
aziendale); stile delle altre figure della tesi (``style.py`` della campagna sperimentale, come ``analisi_tesi_hil.py``).
Misure: quelle della stima ibrida (``analisi_lqg_gawronski.misura_asse``). Il modello e' la struttura di
``analisi_lqg_gawronski.costruisci`` (``Asse.m``), con ritardo e filtro della maglia di corrente come la misura di WorkBench.
Banda 1-150 Hz: la banda della stima e del troncamento del modello interno.

    .venv/bin/python -I analisi_tesi_lineare.py

Figura in Thesis-Writing/plots/lqg/lqg_plant_frf.png.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]

import analisi_lqg_gawronski as g  # noqa: E402
import analisi_tesi_hil as th  # noqa: E402

FIGURA = Path("/home/alessandro/Desktop/Thesis/Thesis-Writing/plots/lqg")
BANDA = (1.0, 150.0)
MODELLO = dict(label="Model without bearings", color="#2C5F99", ls="-", lw=1.3)
MISURA = dict(label="Measurement", color="#000000", ls="-", lw=0.4, marker="o", ms=1.0)


def main() -> None:
    st = th.stile()
    par, Kt, S, T, mot, assi = g.carica()
    f = np.geomspace(*BANDA, 2000)
    fig, ax = st.figure(aspect=0.62, nrows=2, ncols=2, sharex=True)
    for j, (asse, ch, cat, mod) in enumerate(assi):
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        ms = g.misura_asse(par.sperimentale, asse)
        k = (ms.f >= BANDA[0]) & (ms.f <= BANDA[1])
        G = a.m.frf(f)
        ax[0, j].semilogx(f, 20 * np.log10(np.abs(G)), **MODELLO)
        ax[0, j].semilogx(ms.f[k], ms.modulo_db[k], **MISURA)
        ax[1, j].semilogx(f, np.angle(G, deg=True), **{**MODELLO, "label": None})
        ax[1, j].semilogx(ms.f[k], (ms.fase_deg[k] + 180) % 360 - 180, **{**MISURA, "label": None})
        ax[1, j].set_ylim(-180, 180)
        ax[1, j].set_yticks(np.arange(-180, 181, 90))
        ax[1, j].set_xlim(*BANDA)
        ax[1, j].set_xlabel("Frequency (logarithmic scale)")
        ax[0, j].set_title("(%s) %s" % ("ab"[j], ("Azimuth", "Polar axis")[j]), loc="left", fontsize=9)
        print(asse, ms.file, "%.2f A" % ms.eccitazione_A, flush=True)
    for x in ax.ravel():
        # frequenze non divulgabili: griglia logaritmica senza numeri
        x.tick_params(axis="x", which="both", labelbottom=False)
        x.grid(True, which="minor", axis="x", color="#EEEEEE", lw=0.4)
    ax[0, 0].set_ylabel("Magnitude [dB (rad/s)/A]")
    ax[1, 0].set_ylabel("Phase [deg]")
    h, l = ax[0, 0].get_legend_handles_labels()
    th.legenda(fig, h, l, 2)
    fig.tight_layout()
    FIGURA.mkdir(parents=True, exist_ok=True)
    st.save(fig, FIGURA / "lqg_plant_frf.png")


if __name__ == "__main__":
    main()
