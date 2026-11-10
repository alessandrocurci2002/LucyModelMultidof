"""Punto di lavoro della linearizzazione dei cuscinetti nel modello di sintesi dell'LQG-PI.

Richiesta dell'utente (2026-11-01): includere i cuscinetti nel modello lineare dell'LQG con la linearizzazione armonica,
senza sovrastimare la dissipazione ai modi, con priorita' alle risonanze. Lo script confronta il viscoso k_eq + c_eq con
l'elemento di Maxwell, una molla in serie a uno smorzatore, K_b(s) = k_1 s/(s + a), rispetto al plant u -> y ottenuto dalla FRF misurata, costruisce il modello di sintesi e misura il moto
dei cuscinetti sulla circonferenza di Attorrese con la cascata, senza rumore. Non stima i parametri della struttura.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_lqg_punto_lavoro.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE / "src"), str(RADICE), str(RADICE / "closedLoopModel")]

import matplotlib.pyplot as plt  # noqa: E402

import analisi_forcella_lugre as afl  # noqa: E402
from analisi_lqg_gawronski import PASSO_SIM, carica, costruisci, frf, misura_asse, plant_misurato  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402
from lucy_cl import lqg  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning"
BANDE = ((1.0, 5.0), (5.0, 30.0), (30.0, 150.0))
AMPIEZZE = (1e-5, 1e-4, 1e-3, 1e-2)
CENTRO, RAGGIO, QUOTA, OMEGA = 100.0, 10.0, 5.0, 0.2
"""Circonferenza di Attorrese, come in ``closedLoopModel/config_dob_tracking.yaml``."""
BANDA_PUNTO = (1.0, 10.0)
BANDE_PICCO = ((50.0, 62.0), (90.0, 110.0))
BANDE_RISONANZA = {"azimut": ((48.0, 62.0), (62.0, 72.0)), "polare": ((48.0, 64.0), (90.0, 110.0))}
"""Bande con piu' peso nel fit; in azimut la seconda e' il risonatore del drive a 65 Hz."""
PESO_RISONANZE = 3.0
PESI_SCANSIONE = (1.0, 3.0, 10.0, 30.0)
T_AVVIO = 5.0
"""Avvio da fermo sulla circonferenza, escluso dalle metriche come in ``config_dob_tracking.yaml``."""


def dahl(a) -> dict[str, float]:
    return {"sigma0": a.att.sigma0, "Tc": a.att.Tc, "Ts_ratio": 1.0, "vs": 1.0, "sigma1": 0.0, "sigma2": 0.0}


def fit_maxwell(f: np.ndarray, plant, Kb: np.ndarray, br: np.ndarray, peso: float = PESO_RISONANZE, rami: int = 1, k0_libero: bool = False) -> dict:
    """Molla k_i in serie a uno smorzatore c_i = k_i/a_i, K_b(s) = sum k_i s/(s + a_i), con k_0 = 0 salvo ``k0_libero``.

    Il plant u -> y segue quello con la funzione descrittiva all'ampiezza della prova: errore su ln(y/u), 1-150 Hz, con ``peso``
    sui punti ``br`` delle bande di risonanza. Griglia misurata geometrica: ogni banda di ln f pesa uguale. Il riferimento e' il
    modello, non la misura, che resta il controllo."""
    s = 2j * np.pi * f
    Pd = plant(Kb)
    w = np.where(br, float(peso), 1.0)
    n0 = 1 if k0_libero else 0

    def K(x):
        q = np.exp(x)
        return (q[0] if k0_libero else 0.0) + sum(q[n0 + 2 * i] * s / (s + q[n0 + 1 + 2 * i]) for i in range(rami))

    def res(x):
        e = w * np.log(plant(K(x)) / Pd)
        return np.r_[e.real, e.imag]
    partenze = ((0.5,), (3.0,), (10.0,), (50.0,)) if rami == 1 else ((1.0, 30.0), (1.0, 70.0), (5.0, 70.0))
    best = None
    for k1 in (300.0, 3000.0, 20000.0):
        for fa in partenze:
            x0 = [np.log(100.0)] if k0_libero else []
            for x in fa:
                x0 += [np.log(k1), np.log(2 * np.pi * x)]
            r = least_squares(res, x0, x_scale="jac", max_nfev=200)
            if best is None or r.cost < best.cost:
                best = r
    q = np.exp(best.x)
    out = {"k_0": float(q[0]) if k0_libero else 0.0, "h": K(best.x)}
    if rami == 1:
        out |= {"k_1": float(q[n0]), "a": float(q[n0 + 1]), "f_a_Hz": float(q[n0 + 1] / (2 * np.pi)), "c": float(q[n0] / q[n0 + 1])}
    else:
        out["branches"] = [[float(q[n0 + 2 * i]), float(q[n0 + 1 + 2 * i] / (2 * np.pi))] for i in range(rami)]
    return out


def picchi(f, P) -> list[float]:
    """Altezza del picco vicino ai modi strutturali in banda, in dB."""
    return [float(pt._db(P[(f >= a) & (f <= b)]).max()) for a, b in BANDE_PICCO]


def errore(P, Pm, f) -> list[dict[str, float]]:
    out = []
    for lo, hi in BANDE:
        k = (f >= lo) & (f <= hi)
        e = np.log(P[k] / Pm[k])
        out.append({"rms_dB": float(np.sqrt(np.mean((20 / np.log(10) * e.real) ** 2))),
                    "rms_deg": float(np.sqrt(np.mean(np.rad2deg(e.imag) ** 2)))})
    return out


def circonferenza(t: np.ndarray) -> dict[str, np.ndarray]:
    x, y = CENTRO + RAGGIO * np.cos(OMEGA * t), RAGGIO * np.sin(OMEGA * t)
    return {"azimut": np.unwrap(np.arctan2(y, x)), "polare": np.arctan2(QUOTA, np.hypot(x, y))}


def moto_cuscinetto(a, pc: np.ndarray) -> dict:
    """Cascata sulla circonferenza; il moto del cuscinetto e' letto dall'encoder, come in report_dahl_resampled."""
    r = a.cascata.corri(pc - pc[0], 0.0)
    p = r["p"]
    v = np.gradient(p, PASSO_SIM)
    # stato di Dahl normalizzato s = sigma0 z / T_c, aggiornato esattamente a velocita' costante sul passo
    s0, Tc = a.att.sigma0, a.att.Tc
    z, s = 0.0, np.empty(len(p))
    for k, vk in enumerate(v):
        al = s0 * abs(vk) / Tc
        if al > 0:
            e = np.exp(-al * PASSO_SIM)
            z = z * e + vk / al * (1 - e)
        s[k] = s0 * z / Tc
    t = np.arange(len(p)) * PASSO_SIM
    k = t >= T_AVVIO
    pre = np.abs(s[k]) < 0.95
    vk = v[k]
    return {"p": p, "v": v, "s": s, "t": t,
            "velocity_rms_rad_s": float(np.sqrt(np.mean(vk ** 2))), "velocity_max_rad_s": float(np.abs(vk).max()),
            "stroke_urad": float(1e6 * (pc.max() - pc.min())), "presliding_fraction": float(pre.mean()),
            "reversals": int(np.sum(np.diff(np.sign(vk[np.abs(vk) > 1e-5])) != 0)), "i_rms_A": float(np.sqrt(np.mean(r["i"][k] ** 2)))}


def main() -> None:
    par, Kt, S, T, mot, assi = carica()
    cfg, seg = par.sperimentale, par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    out = {}
    fig, ax = plt.subplots(2, 2, figsize=(12, 7.4), sharex=True)
    fig2, ax2 = plt.subplots(2, 2, figsize=(12, 6.2), sharex=True)
    t_c = np.arange(0.0, T_AVVIO + 2 * np.pi / OMEGA, PASSO_SIM)
    ref = circonferenza(t_c)
    for j, (asse, chiave, cat, mod) in enumerate(assi):
        nome = "Azimuth" if asse == "azimut" else "Elevation"
        a = costruisci(asse, chiave, cat, mod, Kt, S, T, mot)
        mis = misura_asse(cfg, asse)
        _, amp = afl.dati_con_ampiezza(cfg, [mis])
        k = (mis.f >= 1) & (mis.f <= 150)
        f, A_m, Pm = mis.f[k], amp[mis.etichetta], plant_misurato(a, mis.f[k], mis.G[k])
        an = a.anello_v
        Bu, Bw, C = an.B[:, 0], an.B[:, 1], an.uscite["p"]
        vv = np.r_[a.m.v, np.zeros(an.N - an.n)]

        R = [np.linalg.solve(2j * np.pi * x * np.eye(an.N) - an.A, np.c_[Bu, Bw]) for x in f]
        g0, gw = np.array([C @ r[:, 0] for r in R]), np.array([C @ r[:, 1] for r in R])
        vu, vw = np.array([vv @ r[:, 0] for r in R]), np.array([vv @ r[:, 1] for r in R])

        def plant(h) -> np.ndarray:
            """y/u con la rigidezza complessa h(f) dei cuscinetti, in rango 1."""
            return g0 + h * gw * vu / (1 - h * vw)

        d = dahl(a)
        Kb = fl.funzione_descrittiva(d, A_m, f)
        b = (f >= BANDA_PUNTO[0]) & (f <= BANDA_PUNTO[1])
        keq, ceq = float(np.median(Kb[b].real)), float(np.median(Kb[b].imag / (2 * np.pi * f[b])))
        br = np.zeros_like(f, bool)
        for lo, hi in BANDE_RISONANZA[asse]:
            br |= (f >= lo) & (f <= hi)
        mx = fit_maxwell(f, plant, Kb, br)
        due = fit_maxwell(f, plant, Kb, br, rami=2)
        libero = fit_maxwell(f, plant, Kb, br, k0_libero=True)
        sf = 2j * np.pi * f
        Kmx = mx.pop("h")
        vx0, vd0 = np.r_[a.m.v, np.zeros(an.N - an.n)], np.r_[np.zeros(an.n), a.m.v, np.zeros(an.N - 2 * an.n)]
        scansione = []
        lo0, hi0 = BANDE_RISONANZA[asse][0]
        m0 = (f >= lo0) & (f <= hi0)
        for peso in PESI_SCANSIONE:
            q = fit_maxwell(f, plant, Kb, br, peso)
            Pq = plant(q.pop("h"))
            Aq = lqg.con_cuscinetti(an.A, Bu, Bw, C, vx0, vd0, 0.0, q["k_1"], q["a"])[0]
            lam = np.linalg.eigvals(Aq)
            lam = lam[(lam.imag > 0) & (np.abs(lam) / (2 * np.pi) > lo0) & (np.abs(lam) / (2 * np.pi) < hi0) & (-lam.real / np.abs(lam) < 0.1)]
            db = pt._db(Pq[m0])
            scansione.append(q | {"weight": peso, "errors": errore(Pq, Pm, f), "peak_Hz": float(f[m0][db.argmax()]), "peak_dB": float(db.max()),
                                  "zeta_mode": float(-lam[0].real / abs(lam[0])) if lam.size else None})
        db = pt._db(Pm[m0])
        picco_mis = [float(f[m0][db.argmax()]), float(db.max())]
        cand = {"none (sliding)": plant(0.0), "stuck (sigma0)": plant(a.att.sigma0), "describing function, measured amplitude": plant(Kb),
                "viscous k_eq + c_eq": plant(keq + sf * ceq), "Maxwell (spring-damper)": plant(Kmx), "two Maxwell branches": plant(due.pop("h")),
                "Maxwell with free k_0": plant(libero.pop("h"))}
        tab = {n: errore(P, Pm, f) for n, P in cand.items()}
        pk = {n: picchi(f, P) for n, P in cand.items()} | {"measured": picchi(f, Pm)}
        for A0 in AMPIEZZE:
            tab[f"describing function, A = {1e6 * A0:.0f} urad"] = errore(plant(fl.funzione_descrittiva(d, np.full_like(f, A0), f)), Pm, f)
        vx, vv = np.r_[a.m.v, np.zeros(an.N - an.n)], np.r_[np.zeros(an.n), a.m.v, np.zeros(an.N - 2 * an.n)]
        Avisc = an.A + keq * np.outer(Bw, vx) + ceq * np.outer(Bw, vv)
        A2, B2, C2 = lqg.con_cuscinetti(an.A, Bu, Bw, C, vx, vv, mx["k_0"], mx["k_1"], mx["a"])
        _, _, _, modi = lqg.forma_modale(A2, B2[:, None], C2[None, :])
        _, _, _, modi_v = lqg.forma_modale(Avisc, Bu[:, None], C[None, :])
        _, _, _, modi_0 = lqg.forma_modale(an.A, Bu[:, None], C[None, :])
        f_low = 1e-4
        residuo = float(abs(C2 @ np.linalg.solve(2j * np.pi * f_low * np.eye(an.N + 1) - A2, B2)) * 2 * np.pi * f_low)
        ris = costruisci(asse, chiave, cat, mod, Kt, S, T, mot, cuscinetti=mx)
        Pr = np.array([(ris.Cp @ np.linalg.solve(2j * np.pi * x * np.eye(len(ris.Ap)) - ris.Ap, ris.Bp))[0, 0] for x in f])
        tab["reduced synthesis model"] = errore(Pr, Pm, f)
        c = lqg.sintesi(ris.Ap, ris.Bp, ris.Cp, 1e-3, 1.0, 100.0, np.zeros(len(ris.Ap) - 1))
        cand["reduced synthesis model"] = Pr
        m = moto_cuscinetto(a, ref[asse])
        k1 = np.argmin(np.abs(f - 1.0))
        scarto = {n: float(np.sqrt(np.mean((20 / np.log(10) * np.log(np.abs(cand[n][br] / Pm[br]))) ** 2))) for n in cand}
        fase = {n: float(np.sqrt(np.mean(np.rad2deg(np.angle(cand[n][br] / Pm[br])) ** 2))) for n in cand}
        out[chiave] = {"gain_at_1Hz_dB": {"maxwell_over_measured": float(pt._db(cand["Maxwell (spring-damper)"][k1] / Pm[k1])),
                                          "no_bearings_over_maxwell": float(pt._db(cand["none (sliding)"][k1] / cand["Maxwell (spring-damper)"][k1]))},
                       "resonance_band_errors": {n: {"rms_dB": scarto[n], "rms_deg": fase[n]} for n in cand},"presliding_scale_urad": 1e6 * a.att.Tc / a.att.sigma0, "frf_amplitude_urad_1_10Hz": [float(1e6 * A_m[b].min()), float(1e6 * np.median(A_m[b])), float(1e6 * A_m[b].max())],
                       "k_eq": keq, "c_eq": ceq, "c_eq_band_Hz": list(BANDA_PUNTO), "errors": tab, "origin_residue": residuo, "maxwell": mx, "two_branches": due, "free_k_0": libero, "peaks_dB": pk,
                       "weight_scan": scansione, "measured_peak_Hz_dB": picco_mis, "resonance_bands_Hz": BANDE_RISONANZA[asse],
                       "modes_with_bearings": [{"f_hz": x.f_hz, "zeta": x.zeta, "hinf": x.hinf} for x in modi],
                       "modes_viscous": [{"f_hz": x.f_hz, "zeta": x.zeta, "hinf": x.hinf} for x in modi_v],
                       "modes_no_bearings": [{"f_hz": x.f_hz, "zeta": x.zeta, "hinf": x.hinf} for x in modi_0],
                       "synthesis": {"order": len(ris.Ap), "modes": [{"f_hz": x.f_hz, "zeta": x.zeta, "hinf": x.hinf} for x in ris.modi],
                                     "test_kp_ki_at_qp100_qi1": [c.k_p, c.k_i], "test_regulator_max_abs": float(np.abs(c.poli_regolatore).max()), "test_estimator_max_abs": float(np.abs(c.poli_stimatore).max())},
                       "circle": {x: m[x] for x in ("velocity_rms_rad_s", "velocity_max_rad_s", "stroke_urad", "presliding_fraction", "reversals", "i_rms_A")}}
        print(nome, json.dumps({x: out[chiave][x] for x in ("presliding_scale_urad", "frf_amplitude_urad_1_10Hz", "k_eq", "c_eq", "maxwell", "peaks_dB", "origin_residue", "synthesis", "circle")}), flush=True)
        print(" poli", [(round(x.f_hz, 2), round(x.zeta, 3)) for x in modi[:7]], flush=True)
        for st_i, n in enumerate(("none (sliding)", "stuck (sigma0)", "Maxwell (spring-damper)"), 1):
            ax[0, j].plot(f, pt._db(cand[n]), lw=1.3, label=n, **pt._stile(st_i))
            ax[1, j].plot(f, np.rad2deg(np.angle(cand[n])), lw=1.1, **pt._stile(st_i))
        ax[0, j].plot(f, pt._db(cand["viscous k_eq + c_eq"]), lw=1.0, color=pt.INCHIOSTRO_SECONDARIO, ls=":", label="viscous k_eq + c_eq")
        ax[1, j].plot(f, np.rad2deg(np.angle(cand["viscous k_eq + c_eq"])), lw=0.9, color=pt.INCHIOSTRO_SECONDARIO, ls=":")
        ax[0, j].plot(f, pt._db(Pm), color=pt.INCHIOSTRO, lw=1.1, marker="o", ms=1.5, label="from measured FRF")
        ax[1, j].plot(f, np.rad2deg(np.angle(Pm)), color=pt.INCHIOSTRO, lw=1.1, marker="o", ms=1.5)
        ax[0, j].set_title(nome, fontsize=9); ax[1, j].set_xlabel("frequency [Hz]"); ax[1, j].set_ylim(-180, 180); ax[1, j].set_yticks(np.arange(-180, 181, 90))
        ax[0, j].set_xscale("log"); ax[1, j].set_xscale("log"); ax[0, j].legend(fontsize=7)
        dec = max(int(0.01 / PASSO_SIM), 1)
        ax2[0, j].plot(m["t"][::dec], m["v"][::dec], color=pt.INCHIOSTRO, lw=1.0); ax2[0, j].set_title(nome, fontsize=9)
        ax2[1, j].plot(m["t"][::dec], m["s"][::dec], color=pt.INCHIOSTRO, lw=1.0); ax2[1, j].set_xlabel("time [s]")
        ax2[1, j].axhline(0.95, color=pt.INCHIOSTRO, lw=0.6, ls=":"); ax2[1, j].axhline(-0.95, color=pt.INCHIOSTRO, lw=0.6, ls=":")
    ax[0, 0].set_ylabel("|y/u| [dB s]"); ax[1, 0].set_ylabel("phase [deg]")
    pt._titolo(fig, "Velocity-loop plant with the bearings at different operating points, against the measured frequency response", seg)
    ax2[0, 0].set_ylabel("bearing velocity [rad/s]"); ax2[1, 0].set_ylabel("sigma0 z / T_c")
    pt._titolo(fig2, "Bearing state on the Attorrese circle with the cascade, no noise: |sigma0 z / T_c| near 1 is sliding", seg)
    pt.salva(fig, CARTELLA / "figures", "plant_operating_points")
    pt.salva(fig2, CARTELLA / "figures", "circle_bearing_state")
    plt.close("all")
    (CARTELLA / "operating_point.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
