"""Figura della tesi (cap. LQG, sezione della taratura): funzione di trasferimento in anello chiuso dal comando alla posizione
dell'encoder, :math:`T_{ry}`, della cascata e dell'LQG-PI tarato, per mostrare la banda dei due controllori.

Richiesta dell'utente (2026-11-03): banda del nuovo controllore rispetto alla cascata; asse delle frequenze senza numeri (policy
aziendale); stile delle altre figure della tesi (``style.py``, come ``analisi_tesi_hil.py``).
Entrambi gli anelli sono chiusi sullo stesso plant lineare, il modello completo dell'asse con i cuscinetti in strisciamento
(``varianti_cuscinetto(a)["sliding"]``), come la banda della cascata nel report di taratura. L'LQG-PI e' quello del capitolo HIL
(``analisi_lqg_circonferenza.casi``, "LQG-PI tuned"), con campionamento a 1 ms e un periodo di ritardo (``anello_lqg``).
Banda: primo attraversamento di -3 dB (scelta A5 del report).

    .venv/bin/python -I analisi_tesi_banda.py

Figura in Thesis-Writing/plots/lqg/lqg_bandwidth.png.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]

import analisi_lqg_circonferenza as cc  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_tesi_hil as th  # noqa: E402

FIGURA = Path("/home/alessandro/Desktop/Thesis/Thesis-Writing/plots/lqg")
BANDA = (0.1, 150.0)
"""Da sotto la banda della cascata alla banda delle misure usate nell'identificazione."""
SOGLIA = dict(color="#777777", ls=":", lw=0.8)


def main() -> None:
    st = th.stile()
    par, Kt, S, T, mot, assi = g.carica()
    Fi = json.loads(cc.TUNING.read_text())
    f = np.geomspace(*BANDA, 4000)
    fig, ax = st.figure(aspect=0.42, nrows=1, ncols=2, sharey=True)
    bande = {}
    for j, (asse, ch, cat, mod) in enumerate(assi):
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        m0 = next(m for m in a.modi if len(m.indici) == 2 and m.zeta < 0.1)
        c = cc.casi(a, Fi[ch]["final"], m0.f_hz)["LQG-PI tuned"]
        an = a.anello_v
        P = g.frf(g.varianti_cuscinetto(a)["sliding"], an.B[:, 0], an.uscite["p"], f)
        L, N = g.anello_lqg(c, P, f)
        T_ry = {"cascade": g.frf(a.cascata.A, a.cascata.B[:, 0], a.cascata.uscite["p"], f), "lqg": N / (1 + L)}
        bande[asse] = {}
        for k in ("cascade", "lqg"):
            stl = th.CTRL[k]
            ax[j].semilogx(f, 20 * np.log10(np.abs(T_ry[k])), color=stl["color"], ls=stl["ls"], lw=stl["lw"], label=stl["label"])
            fb = g.banda(f, T_ry[k])
            ax[j].plot(fb, -3.0, marker=stl["marker"], color=stl["color"], ms=4, ls="none")
            bande[asse][k] = fb
        ax[j].axhline(-3.0, **SOGLIA, label="$-3$ dB")
        ax[j].set_xlim(*BANDA)
        ax[j].set_ylim(-40, 10)
        ax[j].set_xlabel("Frequency (logarithmic scale)")
        ax[j].set_title("(%s) %s" % ("ab"[j], ("Azimuth", "Polar axis")[j]), loc="left", fontsize=9)
        r = bande[asse]["lqg"] / bande[asse]["cascade"]
        print(asse, "k_p %.2f, K_pp %.2f, bandwidth ratio LQG-PI / cascade %.2f" % (c.k_p, a.drive.p.Kpp, r), flush=True)
    for x in ax:
        # frequenze non divulgabili: griglia logaritmica senza numeri
        x.tick_params(axis="x", which="both", labelbottom=False)
        x.grid(True, which="minor", axis="x", color="#EEEEEE", lw=0.4)
    ax[0].set_ylabel("$|T_{ry}|$ [dB]")
    h, l = ax[0].get_legend_handles_labels()
    th.legenda(fig, h, l, 3)
    fig.tight_layout()
    FIGURA.mkdir(parents=True, exist_ok=True)
    st.save(fig, FIGURA / "lqg_bandwidth.png")
    (RADICE / "results_lqg_tuning" / "no_bearings" / "thesis_bandwidth.json").write_text(json.dumps(bande, indent=2))


if __name__ == "__main__":
    main()
