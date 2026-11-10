"""T_ry, T_wy, T_ra e T_wa (Gawronski 2008, Fig. 8.6-8.8) della cascata, dell'LQG-PI sull'errore stimato (tesi) e dell'LQG-PI
sull'errore misurato (Fig. 8.4b, tuning di ``taratura.py``).

r: comando di posizione; w: coppia di disturbo sulla porta dei cuscinetti, come lo scenario 1 della tesi (``cc.t_wy``); y: posizione
dell'encoder; a: accelerazione dell'encoder, a = s^2 y. Plant lineare completo con i cuscinetti in strisciamento, come la Fig. 6.5.
Asse delle frequenze senza numeri (policy aziendale); stile delle figure della tesi.

    .venv/bin/python -I lqg_pi_errore_misurato/funzioni_trasferimento.py

Figure in figures/transfer_functions_{azimuth,polar}.png, margini e bande in results/transfer_functions.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

QUI = Path(__file__).resolve().parent
sys.argv = [sys.argv[0], "--figure"]
sys.path.insert(0, str(QUI))

import controllore as ct  # noqa: E402
import analisi_lqg_circonferenza as cc  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_tesi_hil as th  # noqa: E402

F = np.geomspace(0.05, 300.0, 3000)
STILI = {"cascade": th.CTRL["cascade"],
         "estimated": th.CTRL["lqg"] | {"label": "LQG-PI, estimated error (thesis)"},
         "measured": dict(label="LQG-PI, measured error (Fig. 8.4b)", color="#0072B2", ls="-.", marker="^", lw=1.3)}


def controllori(a, ch: str) -> dict:
    m0 = next(m for m in a.modi if len(m.indici) == 2 and m.zeta < 0.1)
    vecchio = json.loads(cc.TUNING.read_text())[ch]["final"]
    nuovo = json.loads((QUI / "config" / "lqg_tuning.json").read_text())[ch]["final"]
    return {"cascade": None, "estimated": cc.casi(a, vecchio, m0.f_hz)["LQG-PI tuned"],
            "measured": ct.da_stimato(cc.casi(a, nuovo, m0.f_hz)["LQG-PI tuned"])}


def risposte(a, c, A) -> dict:
    an = a.anello_v
    s2 = (2j * np.pi * F) ** 2
    if c is None:
        Ac = a.cascata.A if A is an.A else a.cascata.A + a.att.sigma0 * np.outer(a.cascata.B[:, 1], np.r_[a.m.v, np.zeros(a.cascata.N - a.cascata.n)])
        T_ry = g.frf(Ac, a.cascata.B[:, 0], a.cascata.uscite["p"], F)
        T_wy = g.frf(Ac, a.cascata.B[:, 1], a.cascata.uscite["p"], F)
        L = a.drive.p.Kpp * g.frf(A, an.B[:, 0], an.uscite["p"], F)
    else:
        P = g.frf(A, an.B[:, 0], an.uscite["p"], F)
        L, N = g.anello_lqg(c, P, F)
        T_ry, T_wy = N / (1 + L), g.frf(A, an.B[:, 1], an.uscite["p"], F) / (1 + L)
    return {"T_ry": T_ry, "T_wy": T_wy, "T_ra": s2 * T_ry, "T_wa": s2 * T_wy, "L": L}


def main() -> None:
    st = th.stile()
    par, Kt, S, T, mot, assi = g.carica()
    out = {}
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        C = controllori(a, ch)
        fig, ax = st.figure(aspect=0.78, nrows=2, ncols=2, sharex=True)
        out[ch] = {}
        for k, c in C.items():
            out[ch][k] = {}
            for n, A in g.varianti_cuscinetto(a).items():
                R = risposte(a, c, A)
                m = g.margini(R["L"], F)
                out[ch][k][n] = {"bandwidth_Hz": g.banda(F, R["T_ry"]), "T_wy_max": float(np.abs(R["T_wy"]).max()), "T_wa_max": float(np.abs(R["T_wa"]).max())} | m
                if n != "sliding":
                    continue
                s = STILI[k]
                for x, q in zip(ax.ravel(), ("T_ry", "T_wy", "T_ra", "T_wa")):
                    x.loglog(F, np.abs(R[q]), color=s["color"], ls=s["ls"], lw=s["lw"], label=s["label"])
        for x, q, u in zip(ax.ravel(), ("T_ry", "T_wy", "T_ra", "T_wa"), ("[-]", "[rad/(N m)]", "[1/s$^2$]", "[1/(kg m$^2$)]")):
            x.set_title("$|%s_{%s}|$ %s" % (q[0], q[2:], u), loc="left", fontsize=st.FONT_SIZE)
            x.set_xlim(F[0], F[-1])
            # frequenze non divulgabili: griglia logaritmica senza numeri
            x.tick_params(axis="x", which="both", labelbottom=False)
            x.grid(True, which="minor", axis="x", color="#EEEEEE", lw=0.4)
        for x in ax[-1]:
            x.set_xlabel("Frequency (logarithmic scale)")
        h, l = ax[0, 0].get_legend_handles_labels()
        th.legenda(fig, h, l, 3)
        fig.tight_layout()
        (QUI / "figures").mkdir(exist_ok=True)
        st.save(fig, QUI / "figures" / f"transfer_functions_{'azimuth' if asse == 'azimut' else 'polar'}.png")
    (QUI / "results").mkdir(exist_ok=True)
    (QUI / "results" / "transfer_functions.json").write_text(json.dumps(out, indent=1))
    for ch, d in out.items():
        for k, v in d.items():
            print(ch, k, {n: (round(x["bandwidth_Hz"], 2), round(x["modulus"], 3), round(x["pm_deg"], 1), "%.2e" % x["T_wy_max"]) for n, x in v.items()})


if __name__ == "__main__":
    main()
