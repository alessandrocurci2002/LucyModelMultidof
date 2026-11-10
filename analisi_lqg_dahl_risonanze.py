"""Dahl linearizzato armonicamente con priorita' alle risonanze invece che alla bassa frequenza.

Richiesta dell'utente (2026-11-01). Sul plant u -> y della maglia di velocita', a confronto con la misura:

* Dahl (funzione descrittiva all'ampiezza della prova): riferimento, non e' un modello di stato;
* Dahl linearizzato all'ampiezza della risonanza principale: K_b costante, k(1 + j eta), a una sola ampiezza;
* Maxwell, k_0 + k_1 s/(s + a), fittato sul riferimento con peso quasi nullo fuori dalle bande di risonanza: modello di stato.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_lqg_dahl_risonanze.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE / "src"), str(RADICE), str(RADICE / "closedLoopModel")]

import matplotlib.pyplot as plt  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

import analisi_forcella_lugre as afl  # noqa: E402
from analisi_lqg_gawronski import carica, costruisci, misura_asse, plant_misurato  # noqa: E402
from analisi_lqg_punto_lavoro import dahl, errore  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning" / "resonance_priority"
BANDE = {"azimut": ((48.0, 62.0), (62.0, 72.0)), "polare": ((48.0, 64.0), (90.0, 110.0))}
PESO_FUORI = 0.03
"""Peso del fit fuori dalle bande di risonanza."""


def fit_risonanze(f, plant, Pd, br) -> dict:
    s = 2j * np.pi * f
    w = np.where(br, 1.0, PESO_FUORI)

    def K(x):
        q = np.exp(x)
        return q[0] + q[1] * s / (s + q[2])

    def res(x):
        e = w * np.log(plant(K(x)) / Pd)
        return np.r_[e.real, e.imag]
    best = None
    for k0 in (10.0, 300.0, 3000.0):
        for k1 in (300.0, 3000.0, 20000.0):
            for fa in (1.0, 10.0, 50.0, 200.0):
                r = least_squares(res, np.log([k0, k1, 2 * np.pi * fa]), x_scale="jac", max_nfev=150)
                if best is None or r.cost < best.cost:
                    best = r
    q = np.exp(best.x)
    return {"k_0": float(q[0]), "k_1": float(q[1]), "a": float(q[2]), "f_a_Hz": float(q[2] / (2 * np.pi)), "h": K(best.x)}


def picchi(f, P, bande) -> list[list[float]]:
    out = []
    for lo, hi in bande:
        m = (f >= lo) & (f <= hi)
        db = pt._db(P[m])
        out.append([float(f[m][db.argmax()]), float(db.max())])
    return out


def in_banda(P, Pm, br) -> list[float]:
    e = np.log(P[br] / Pm[br])
    return [float(np.sqrt(np.mean((20 / np.log(10) * e.real) ** 2))), float(np.sqrt(np.mean(np.rad2deg(e.imag) ** 2)))]


def main() -> None:
    par, Kt, S, T, mot, assi = carica()
    cfg, seg = par.sperimentale, par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    O = json.loads((RADICE / "results_lqg_tuning" / "operating_point.json").read_text())
    out = {}
    for asse, chiave, cat, mod in assi:
        nome = "Azimuth" if asse == "azimut" else "Elevation"
        a = costruisci(asse, chiave, cat, mod, Kt, S, T, mot)
        ms = misura_asse(cfg, asse)
        _, amp = afl.dati_con_ampiezza(cfg, [ms])
        k = (ms.f >= 1) & (ms.f <= 150)
        f, A, Pm = ms.f[k], amp[ms.etichetta], plant_misurato(a, ms.f[k], ms.G[k])
        an = a.anello_v
        Bu, Bw, C = an.B[:, 0], an.B[:, 1], an.uscite["p"]
        vv = np.r_[a.m.v, np.zeros(an.N - an.n)]
        vd = np.r_[np.zeros(an.n), a.m.v, np.zeros(an.N - 2 * an.n)]
        R = [np.linalg.solve(2j * np.pi * x * np.eye(an.N) - an.A, np.c_[Bu, Bw]) for x in f]
        g0, gw = np.array([C @ r[:, 0] for r in R]), np.array([C @ r[:, 1] for r in R])
        vu, vw = np.array([vv @ r[:, 0] for r in R]), np.array([vv @ r[:, 1] for r in R])

        def plant(h):
            return g0 + h * gw * vu / (1 - h * vw)
        d = dahl(a)
        Kb = fl.funzione_descrittiva(d, A, f)
        Pd = plant(Kb)
        s = 2j * np.pi * f
        br = np.zeros_like(f, bool)
        for lo, hi in BANDE[asse]:
            br |= (f >= lo) & (f <= hi)
        m0 = (f >= BANDE[asse][0][0]) & (f <= BANDE[asse][0][1])
        A_r = float(np.median(A[m0]))
        K_r = complex(fl.funzione_descrittiva(d, np.array([A_r]), np.array([56.0]))[0])
        mx = fit_risonanze(f, plant, Pd, br)
        Pmx = plant(mx.pop("h"))
        Pk = plant(K_r * np.ones_like(f))
        mp = O[chiave]["maxwell"]
        Pp = plant(mp["k_0"] + mp["k_1"] * s / (s + mp["a"]))
        P0 = plant(0.0 * s)
        A2 = np.block([[an.A + mx["k_0"] * np.outer(Bw, vv), mx["k_1"] * Bw[:, None]], [vd[None, :], -mx["a"] * np.ones((1, 1))]])
        lam = np.linalg.eigvals(A2)
        lam = lam[(lam.imag > 0) & (np.abs(lam) / (2 * np.pi) > 45) & (np.abs(lam) / (2 * np.pi) < 110) & (-lam.real / np.abs(lam) < 0.1)]
        curve = {"no bearings": P0, "true Dahl, describing function at the test amplitude (reference)": Pd,
                 f"Dahl linearised at the resonance amplitude, A = {1e6 * A_r:.0f} urad": Pk, "Maxwell, previous fit (low frequency priority)": Pp,
                 "Maxwell, resonance priority (state space)": Pmx}
        out[chiave] = {"A_resonance_urad": 1e6 * A_r, "K_b_at_resonance_amplitude": [K_r.real, K_r.imag], "eta": K_r.imag / K_r.real, "maxwell_resonance": mx,
                       "damping": [[float(abs(x) / (2 * np.pi)), float(-x.real / abs(x))] for x in sorted(lam, key=abs)],
                       "curves": {n: {"errors": errore(P, Pm, f), "peaks_Hz_dB": picchi(f, P, BANDE[asse]), "band_vs_measured_dB_deg": in_banda(P, Pm, br)} for n, P in curve.items()},
                       "measured_peaks_Hz_dB": picchi(f, Pm, BANDE[asse])}
        fig, ax = plt.subplots(2, 2, figsize=(14, 7.6), gridspec_kw={"width_ratios": [1.15, 1]})
        for j, (lo, hi) in enumerate(((1.0, 150.0), (40.0, 110.0))):
            z = (f >= lo) & (f <= hi)
            for i, (n, P) in enumerate(curve.items(), 1):
                st = {"color": pt.INCHIOSTRO_SECONDARIO, "linestyle": ":"} if n == "no bearings" else pt._stile(i)
                ax[0, j].plot(f[z], pt._db(P[z]), lw=1.3 if "resonance priority" in n else 1.0, label=n, **st)
                ax[1, j].plot(f[z], np.rad2deg(np.angle(P[z])), lw=1.3 if "resonance priority" in n else 1.0, **st)
            ax[0, j].plot(f[z], pt._db(Pm[z]), color=pt.INCHIOSTRO, lw=1.2, marker="o", ms=1.8, label="from measured FRF")
            ax[1, j].plot(f[z], np.rad2deg(np.angle(Pm[z])), color=pt.INCHIOSTRO, lw=1.2, marker="o", ms=1.8)
            for q in ax[:, j]:
                q.set_xscale("log")
            ax[1, j].set_xlabel("frequency [Hz]"); ax[1, j].set_ylim(-180, 180); ax[1, j].set_yticks(np.arange(-180, 181, 90))
        ax[0, 1].set_title("zoom on the resonances", fontsize=9); ax[0, 0].set_title("whole band", fontsize=9)
        ax[0, 0].set_ylabel("|y/u| [dB s]"); ax[1, 0].set_ylabel("phase [deg]"); ax[0, 0].legend(fontsize=7, loc="lower left")
        pt._titolo(fig, f"{nome}: velocity-loop plant, Dahl linearised with priority on the resonances", seg)
        pt.salva(fig, CARTELLA / "figures", f"{chiave}_resonance_priority")
        plt.close(fig)
    (CARTELLA / "summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
