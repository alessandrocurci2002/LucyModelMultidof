r"""Cascata ritarata (P di posizione, PI di velocita', filtri dell'anello di velocita') con gli stessi vincoli dell'LQG-PI.

Richiesta dell'utente (2026-11-03): verificare che l'LQG-PI dia piu' banda anche rispetto a una cascata ritarata, sullo stesso
modello e con gli stessi vincoli. Plant: il modello dell'asse di ``analisi_lqg_gawronski.costruisci`` (struttura, motore, drive
AKD2G), con i cuscinetti di Dahl nelle simulazioni dei doppietti e in strisciamento o bloccati nelle verifiche lineari.

Parametri liberi: K_pp, K_pv, K_iv del drive e i due biquad sul ramo in avanti dell'anello di velocita'. Restano quelli del drive
l'anello di corrente (banda circa 1.5 kHz, non limita) e i filtri in retroazione. Nessun integrale di posizione (K_ip = 0), come
sulla cascata di LUCY; nessun feedforward.

Vincoli, gli stessi dell'LQG-PI (scelte A1-A3 del report di taratura):
    - picco del comando di corrente sotto 18 A sui doppietti, uno per ampiezza, plant Dahl, senza limite di corrente;
    - min |1 + L| >= 0.5 sull'anello di posizione (aperto al comando di velocita'), plant in strisciamento;
    - anello stabile con i cuscinetti in strisciamento e bloccati.
In piu', scelta mia: min |1 + L_v| dell'anello di velocita' (aperto al comando di corrente) non sotto
min(0.5, valore della cascata attuale), perche' l'LQG lascia l'anello di velocita' com'e' e una cascata ritarata non deve
ottenere banda rendendolo meno robusto di oggi.

Filtri provati sul ramo in avanti (le frequenze sono quelle esatte del modello: e' un vantaggio per la cascata):
    ``drive``  i filtri attuali del drive;
    ``notch``  un notch di -20 dB su ogni modo strutturale poco smorzato sotto 150 Hz (al piu' due);
    ``notch+lp`` notch sul modo principale e passa-basso del secondo ordine a 200 Hz.

Ricerca: griglia di K_pv (x0.5-x40 dell'attuale) e K_iv (x1, x3, x10); per ciascun punto il K_pp piu' alto che rispetta i vincoli
lineari (griglia logaritmica, si ferma al primo K_pp che li viola); sui candidati migliori si simulano i doppietti e, se la corrente
supera 18 A, K_pp si abbassa per bisezione. Banda: primo attraversamento di -3 dB di T_ry (scelta A5), sul plant in strisciamento.

    .venv/bin/python -I analisi_cascata_ritarata.py

Risultati in results_lqg_tuning/no_bearings/cascade_retuned/summary.json.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]

import analisi_anello_chiuso_ibrido as aac  # noqa: E402
import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_passi_2_4 as p24  # noqa: E402
from lucy_cl.driver import Biquad, Drive  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning" / "no_bearings" / "cascade_retuned"
I_MAX, MODULO = 18.0, 0.5
F = np.geomspace(0.05, 300.0, 3000)
S_KPV = np.geomspace(0.5, 40.0, 17)
S_KIV = (1.0, 3.0, 10.0)
KPP = np.geomspace(2.0, 3000.0, 60)
N_CANDIDATI = 6
UNITA = Biquad(500.0, 0.5, 500.0, 0.5, 0)


def notch(f_hz: float) -> Biquad:
    """-20 dB al centro (Q_p / Q_z = 0.1), poli ben smorzati."""
    return Biquad(f_hz, 7.0, f_hz, 0.7, 4)


def filtri(a: g.Asse) -> dict[str, tuple[Biquad, Biquad]]:
    m = a.m
    f = m.omega / (2 * np.pi)
    modi = [x for x, z in zip(f, m.zeta) if 1.0 < x <= 150.0 and z < 0.1][:2]
    nn = [notch(x) for x in modi] + [UNITA] * (2 - len(modi))
    return {"drive": a.drive.p.forward, "notch": (nn[0], nn[1]), "notch+lp": (notch(modi[0]), Biquad(5000.0, 0.707107, 200.0, 0.707107, 4))}


def cascata(a: g.Asse, p) -> aac.AnelloChiuso:
    return aac.AnelloChiuso(a.m, Drive(p), MOT, a.att, g.PASSO_SIM)


def stuck(an: aac.AnelloChiuso, a: g.Asse) -> np.ndarray:
    vx = np.r_[a.m.v, np.zeros(an.N - an.n)]
    return an.A + a.att.sigma0 * np.outer(an.B[:, 1], vx)


def stabile(A: np.ndarray, integratore: bool = False) -> bool:
    """Tutti i poli nel semipiano sinistro; con ``integratore`` si ammette un polo nell'origine: nell'anello di velocita'
    da solo la posizione non e' chiusa e, a cuscinetti bloccati, l'integrale della velocita' e' legato alla posizione."""
    lam = np.linalg.eigvals(A)
    if integratore:
        lam = np.delete(lam, np.argmin(np.abs(lam)))
    return bool(lam.real.max() < -1e-9)


def anello_velocita(a: g.Asse, p) -> dict:
    """Anello di velocita' chiuso (comando di velocita' -> posizione) e modulo minimo di 1 + L_v, con
    1 + L_v = 1 / (1 - T_v H), T_v = s y / u, H filtro in retroazione del drive (osservatore compreso)."""
    pv = dataclasses.replace(p, Kpp=0.0, Kip_hz=0.0, kvff=1.0)
    an = aac.AnelloChiuso(a.m, g.DriveVelocita(Drive(pv)), MOT, a.att, g.PASSO_SIM)
    s = 2j * np.pi * F
    Pv = g.frf(an.A, an.B[:, 0], an.uscite["p"], F)
    H = Drive(pv).F_fb(s)
    return {"Pv": Pv, "mod_v": float(np.min(1 / np.abs(1 - s * Pv * H))), "stabile": stabile(an.A, integratore=True) and stabile(stuck(an, a), integratore=True)}


def lineare(a: g.Asse, p, Pv: np.ndarray) -> dict:
    an = cascata(a, p)
    L = p.Kpp * Pv
    T = g.frf(an.A, an.B[:, 0], an.uscite["p"], F)
    return {"mod": float(np.abs(1 + L).min()), "stabile": stabile(an.A) and stabile(stuck(an, a)), "banda_Hz": g.banda(F, T),
            "picco_T_dB": float(20 * np.log10(np.abs(T).max())), "scarto_T": float(np.abs(T - L / (1 + L)).max())}


def corrente(a: g.Asse, p, W) -> dict:
    b = dataclasses.replace(a, cascata=cascata(a, p), drive=Drive(p))
    ms, picco = [], 0.0
    for w in W:
        r = g.simula(b, None, w)
        ms.append(g.metriche(w, r["p"], r["i"]))
        picco = max(picco, r["max_i_cmd"])
    return {"i_max_A": picco, "err_rms_median_mdeg": float(np.median([m["err_rms_mdeg"] for m in ms])),
            "t_set_median_s": float(np.median([m["t_set_s"] for m in ms])), "over_max_mdeg": float(max(m["over_mdeg"] for m in ms))}


def kpp_massimo(a: g.Asse, p, Pv) -> tuple[float, dict] | None:
    """K_pp piu' alto della griglia che rispetta i vincoli lineari, fermandosi alla prima violazione."""
    ok = None
    for k in KPP:
        an = cascata(a, dataclasses.replace(p, Kpp=float(k)))
        if np.abs(1 + k * Pv).min() < MODULO or not (stabile(an.A) and stabile(stuck(an, a))):
            break
        ok = float(k)
    return None if ok is None else (ok, lineare(a, dataclasses.replace(p, Kpp=ok), Pv))


def asse_ritarato(a: g.Asse) -> dict:
    p0 = a.drive.p
    W = p24.comandi_tipici(a.nome)
    v0 = anello_velocita(a, p0)
    soglia_v = min(MODULO, v0["mod_v"])
    base = {"Kpp": p0.Kpp, "Kpv": p0.Kpv, "Kiv_hz": p0.Kiv_hz, "mod_v": v0["mod_v"]} | lineare(a, p0, v0["Pv"]) | corrente(a, p0, W)
    print(a.nome, "attuale", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in base.items()}, flush=True)
    griglia = []
    for nf, fw in filtri(a).items():
        for sv in S_KPV:
            for si in S_KIV:
                p = dataclasses.replace(p0, Kpv=p0.Kpv * sv, Kiv_hz=p0.Kiv_hz * si, forward=fw)
                v = anello_velocita(a, p)
                if v["mod_v"] < soglia_v or not v["stabile"]:
                    continue
                r = kpp_massimo(a, p, v["Pv"])
                if r is None:
                    continue
                griglia.append({"filters": nf, "s_Kpv": float(sv), "s_Kiv": si, "Kpv": p.Kpv, "Kiv_hz": p.Kiv_hz, "Kpp": r[0], "mod_v": v["mod_v"]} | r[1])
        print(a.nome, nf, "punti ammessi", sum(x["filters"] == nf for x in griglia), flush=True)
    griglia.sort(key=lambda x: -x["banda_Hz"])
    candidati = []
    for c in griglia[:N_CANDIDATI]:
        p = dataclasses.replace(p0, Kpv=c["Kpv"], Kiv_hz=c["Kiv_hz"], forward=filtri(a)[c["filters"]])
        lo, hi = KPP[0], c["Kpp"]
        x = corrente(a, dataclasses.replace(p, Kpp=hi), W)
        if x["i_max_A"] > I_MAX:  # la corrente cresce con K_pp: bisezione in scala logaritmica
            for _ in range(10):
                k = float(np.sqrt(lo * hi))
                if corrente(a, dataclasses.replace(p, Kpp=k), W)["i_max_A"] <= I_MAX:
                    lo = k
                else:
                    hi = k
            hi = lo
            x = corrente(a, dataclasses.replace(p, Kpp=hi), W)
        q = c | {"Kpp": hi} | lineare(a, dataclasses.replace(p, Kpp=hi), anello_velocita(a, p)["Pv"]) | x
        q["current_limited"] = hi < c["Kpp"]
        candidati.append(q)
        print(a.nome, "candidato", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in q.items()}, flush=True)
    migliore = max(candidati, key=lambda x: x["banda_Hz"])
    return {"current": base, "velocity_loop_threshold": soglia_v, "grid": griglia, "candidates": candidati, "best": migliore}


def main() -> None:
    global MOT
    par, Kt, S, T, mot, assi = g.carica()
    MOT = mot
    CARTELLA.mkdir(parents=True, exist_ok=True)
    lqg_banda = json.loads((RADICE / "results_lqg_tuning" / "no_bearings" / "thesis_bandwidth.json").read_text())
    out = {}
    for asse, ch, cat, mod in assi:
        a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        out[asse] = asse_ritarato(a)
        b = out[asse]["best"]
        print(f"{asse}: banda cascata attuale {out[asse]['current']['banda_Hz']:.2f} Hz, ritarata {b['banda_Hz']:.2f} Hz "
              f"({b['filters']}, Kpv x{b['s_Kpv']:.2f}, Kiv x{b['s_Kiv']:g}, Kpp {b['Kpp']:.1f}, {b['i_max_A']:.1f} A), "
              f"LQG-PI {lqg_banda[asse]['lqg']:.2f} Hz", flush=True)
    (CARTELLA / "summary.json").write_text(json.dumps(out, indent=1))


MOT = None

if __name__ == "__main__":
    main()
