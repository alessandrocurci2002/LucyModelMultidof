"""Doppietti V2 nel tempo con i cuscinetti lineari di Maxwell al posto dell'attrito di Dahl.

Richiesta dell'utente (2026-11-01): verificare se la modellazione di Maxwell del modello di sintesi dell'LQG perde
correlazione nel tempo, in anello aperto (corrente misurata) e in anello chiuso (riferimento di posizione), con i guadagni
dei file .akd e con quelli stimati sul doppietto da 1 grado come negli altri casi. Riferimento: stessi calcoli con Dahl
(griglia uniforme in ln f). Il modello lineare non ha saturazioni: i picchi si confrontano dopo con i limiti del drive.

Il modello e' lineare, quindi la corrente di tenuta iniziale I0 si tratta per sovrapposizione: posizione relativa all'inizio,
corrente del modello = I0 + risposta. Con Maxwell si prova anche ``dahl_gains``: i guadagni stimati per Dahl, senza ristima.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_maxwell_tempo.py
"""

from __future__ import annotations

import csv
import dataclasses
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE / "src"), str(RADICE), str(RADICE / "closedLoopModel")]

import matplotlib.pyplot as plt  # noqa: E402
from scipy.linalg import expm  # noqa: E402
from scipy.signal import savgol_filter  # noqa: E402

import analisi_anello_chiuso_ibrido as aac  # noqa: E402
from analisi_lqg_gawronski import ANOMALA, carica  # noqa: E402
from gimbal_mdof import doppietti_tempo as dt  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402
from lucy_cl.driver import Drive, DriveParams  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning" / "time_domain_maxwell"
PUNTO = RADICE / "results_lqg_tuning" / "operating_point.json"
ETI = {"dahl": "model, Dahl bearings", "maxwell": "model, Maxwell bearings"}


class AnelloMaxwell(aac.AnelloChiuso):
    """``AnelloChiuso`` lineare: cuscinetti k_0 + k_1 s/(s + a) al posto dell'attrito; ingresso il comando di posizione."""

    cuscinetti: dict = {}

    def __init__(self, m, drive, motore, att, passo=aac.PASSO):
        super().__init__(m, drive, motore, att, passo)
        n, N = self.n, self.N
        vx = np.r_[m.v, np.zeros(N - n)]
        vv = np.r_[np.zeros(n), m.v, np.zeros(N - 2 * n)]
        c = self.cuscinetti
        A2, B2, _ = lqg.con_cuscinetti(self.A, self.B[:, 0], self.B[:, 1], np.zeros(N), vx, vv, c["k_0"], c["k_1"], c["a"])
        E = expm(np.block([[A2, B2[:, None]], [np.zeros((1, N + 2))]]) * passo)
        self.Ad2, self.Bd2 = E[:N + 1, :N + 1], E[:N + 1, N + 1]
        est = lambda r: np.r_[r, 0.0]  # noqa: E731
        self.r_p, self.r_i, self.r_ic = est(self.uscite["p"]), est(self.uscite["i"]), est(self.ic_row)
        self.r_vp, self.r_vi, self.i_ic, self.i_vp = est(self.vp_row), est(self.vint_row), self.ic_in[0], self.vp_in[0]

    def corri(self, p_cmd: np.ndarray, I0: float) -> dict[str, np.ndarray]:
        s = np.zeros(self.N + 1)
        P, I, IC, VP, VI = (np.empty(len(p_cmd)) for _ in range(5))
        for k, pc in enumerate(p_cmd):
            P[k], I[k], IC[k] = self.r_p @ s, I0 + self.r_i @ s, I0 + self.r_ic @ s + self.i_ic * pc
            VP[k], VI[k] = self.r_vp @ s + self.i_vp * pc, I0 + self.r_vi @ s
            s = self.Ad2 @ s + self.Bd2 * pc
        return {"p": P, "i": I, "i_cmd": IC, "v_phase": VP, "i_int": VI, "Tb0": 0.0}


class AnelloApertoMaxwell:
    """Struttura e Maxwell con la corrente misurata in ingresso, stesso passo e stessa struttura di ``Simulatore``."""

    def __init__(self, m, c, passo=dt.DT):
        n = len(m.M)
        S = dt.Simulatore(m, passo)
        C = m.M @ m.Phi @ np.diag(2 * m.zeta * m.omega) @ m.Phi.T @ m.M
        Mi = np.linalg.inv(m.M)
        A = np.block([[np.zeros((n, n)), np.eye(n)], [-Mi @ S.K0, -Mi @ C]])
        Bt, Bw = np.r_[np.zeros(n), Mi @ m.e], np.r_[np.zeros(n), -Mi @ m.v]
        A2, B2, _ = lqg.con_cuscinetti(A, Bt, Bw, np.zeros(2 * n), np.r_[m.v, np.zeros(n)], np.r_[np.zeros(n), m.v], c["k_0"], c["k_1"], c["a"])
        E = expm(np.block([[A2, B2[:, None]], [np.zeros((1, 2 * n + 2))]]) * passo)
        self.Ad, self.Bd, self.c, self.Kt, self.dt = E[:2 * n + 1, :2 * n + 1], E[:2 * n + 1, 2 * n + 1], np.r_[m.c, np.zeros(n + 1)], m.Kt, passo

    def corri(self, I: np.ndarray) -> np.ndarray:
        tau = self.Kt * (I - I[:int(0.4 / self.dt)].mean())
        s, out = np.zeros(len(self.Ad)), np.empty(len(I))
        for k, t in enumerate(tau):
            out[k] = self.c @ s
            s = self.Ad @ s + self.Bd * t
        return out


def maschera(w):
    vc = np.abs(np.gradient(w.cmd, dt.DT)) > 1e-6
    return tuple(np.flatnonzero(vc)[[0, -1]])


def fig_aperto(nome, curve, seg):
    fig, ax = plt.subplots(3, 3, figsize=(15, 9), sharex="col", gridspec_kw={"height_ratios": [2, 1.2, 1]})
    for j, key in enumerate(aac.FIGURA):
        w, c = curve[key]
        i0 = maschera(w)[0]
        k = w.t <= w.t[i0] + 2.5
        t = w.t[k]
        md = lambda q: np.rad2deg(q[k]) * 1e3  # noqa: E731
        ax[0, j].plot(t, md(w.cmd), color=pt.INCHIOSTRO, lw=0.8, ls=":", label="command")
        ax[0, j].plot(t, md(w.pos), color=pt.INCHIOSTRO, lw=1.6, label="measured encoder")
        ax[1, j].plot(t, savgol_filter(w.pos, 41, 3, deriv=1, delta=dt.DT)[k], color=pt.INCHIOSTRO, lw=1.2)
        for i, (n, y) in enumerate(c.items(), 1):
            ax[0, j].plot(t, md(y), lw=1.3, label=ETI[n], **pt._stile(i))
            ax[1, j].plot(t, np.gradient(y, dt.DT)[k], lw=1.0, **pt._stile(i))
        ax[2, j].plot(t, w.I[k], color=pt.INCHIOSTRO, lw=1.0)
        ax[0, j].set_title(f"{key[0].split('_')[0].replace('01deg', '0.1deg').replace('05deg', '0.5deg')} doublet, move {key[1]}", fontsize=9)
        ax[2, j].set_xlabel("time from window start [s]")
    ax[0, 0].set_ylabel("position [mdeg]"); ax[1, 0].set_ylabel("velocity [rad/s]"); ax[2, 0].set_ylabel("input current [A rms]")
    ax[0, 1].legend(loc="best", fontsize=7)
    pt._titolo(fig, f"{nome}: open loop, measured current input, Dahl and Maxwell bearings", seg)
    return fig


def fig_chiuso(nome, curve, tit, seg):
    fig, ax = plt.subplots(3, 3, figsize=(15, 9.5), sharex="col", gridspec_kw={"height_ratios": [2, 1.3, 1.6]})
    for j, key in enumerate(aac.FIGURA):
        w, c = curve[key]
        i0 = maschera(w)[0]
        k = (w.t >= w.t[i0] - 0.2) & (w.t <= w.t[i0] + 2.0)
        t = w.t[k]
        md = lambda q: np.rad2deg(q[k]) * 1e3  # noqa: E731
        rif = w.cmd - w.cmd[0]
        ax[0, j].plot(t, md(w.cmd), color=pt.INCHIOSTRO, lw=0.8, ls=":", label="position command (model input)")
        ax[0, j].plot(t, md(w.pos), color=pt.INCHIOSTRO, lw=1.6, label="measured encoder")
        ax[1, j].plot(t, md(w.pos - w.cmd), color=pt.INCHIOSTRO, lw=1.3)
        ax[2, j].plot(t, w.I[k], color=pt.INCHIOSTRO, lw=1.1, label="measured (AMPS / 2)")
        for i, (n, (y, im)) in enumerate(c.items(), 1):
            ax[0, j].plot(t, md(y), lw=1.3, label=ETI[n], **pt._stile(i))
            ax[1, j].plot(t, md(y - rif), lw=1.1, **pt._stile(i))
            ax[2, j].plot(t, im[k], lw=1.0, **pt._stile(i))
        ax[0, j].set_title(f"{key[0].split('_')[0].replace('01deg', '0.1deg').replace('05deg', '0.5deg')} doublet, move {key[1]}", fontsize=9)
        ax[2, j].set_xlabel("time from window start [s]")
    ax[0, 0].set_ylabel("position [mdeg]"); ax[1, 0].set_ylabel("following error [mdeg]"); ax[2, 0].set_ylabel("motor current [A rms]")
    ax[0, 1].legend(loc="lower right", fontsize=7); ax[2, 0].legend(fontsize=7)
    pt._titolo(fig, f"{nome}: closed loop, measured position command input, {tit}", seg)
    return fig


def main() -> None:
    par, Kt, S, T, mot, assi = carica()
    seg = par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    O = json.loads(PUNTO.read_text())
    r_ol, r_cl, sommario = [], [], {}
    for asse, chiave, cat, mod in assi:
        nome = "Azimuth" if asse == "azimut" else "Elevation"
        c = O[chiave]["maxwell"]
        v = S[chiave]["uniform_lnf"]["hybrid"]["parameters"]
        m = mod(cat, v, Kt, {}).struttura
        att = dt.Attrito(v["sigma0"], v["Tc"])
        akd, wb = aac.DRIVE[asse]
        p_file = DriveParams.from_files(chiave, aac.VCL / akd, aac.VCL / wb, pl_kp_unit="hz")
        fin = {a: dt.finestre(aac.VCL / "V2" / (aac.PROVE[asse] % a), asse) for a in aac.AMPIEZZE}
        classe = type("Anello", (AnelloMaxwell,), {"cuscinetti": c})
        g_max, costo = aac.stima_guadagni(m, dataclasses.replace(p_file, **aac.SETTEMBRE[asse]), mot, att, fin[aac.STIMA], anello=classe)
        g_dahl = T["gains"][chiave]["uniform_lnf"]
        print(asse, "guadagni Maxwell", g_max, "Dahl", g_dahl, flush=True)
        anelli = {("dahl", "estimated"): aac.AnelloChiuso(m, Drive(dataclasses.replace(p_file, **g_dahl)), mot, att),
                  ("dahl", "akd"): aac.AnelloChiuso(m, Drive(p_file), mot, att),
                  ("maxwell", "estimated"): classe(m, Drive(dataclasses.replace(p_file, **g_max)), mot, att),
                  ("maxwell", "akd"): classe(m, Drive(p_file), mot, att),
                  ("maxwell", "dahl_gains"): classe(m, Drive(dataclasses.replace(p_file, **g_dahl)), mot, att)}
        stab = {}
        for (n, g), an in anelli.items():
            Ad = an.Ad if n == "dahl" else an.Ad2
            stab[f"{n}_{g}"] = float(np.abs(np.linalg.eigvals(Ad)).max())
        sim_d, sim_m = dt.Simulatore(m), AnelloApertoMaxwell(m, c)
        c_ol, c_cl = {}, {}
        picchi = {f"{n}_{g}": 0.0 for n, g in anelli}
        for a, ww in fin.items():
            for w in ww:
                key = (a, w.nome[-2:])
                i0, i1 = maschera(w)
                yd, ym = sim_d.corri(w.I, att), sim_m.corri(w.I)
                mv = slice(i0, min(i1 + int(0.5 / dt.DT), len(yd)))
                for n, y in (("dahl", yd), ("maxwell", ym)):
                    r_ol.append([asse, a, w.nome[-2:], n, dt.fit(y, w.pos), dt.fit(y[mv], w.pos[mv])])
                c_ol[key] = (w, {"dahl": yd, "maxwell": ym})
                t = np.arange(0.0, w.t[-1], aac.PASSO)
                pc, I0 = np.interp(t, w.t, w.cmd - w.cmd[0]), float(w.I[:int(0.4 / dt.DT)].mean())
                mv = slice(i0, min(i1 + int(1.0 / dt.DT), len(w.t)))
                for n, g in anelli:
                    c_cl.setdefault(g, {}).setdefault(key, (w, {}))
                    r = anelli[(n, g)].corri(pc, I0)
                    y, im = np.interp(w.t, t, r["p"]), np.interp(w.t, t, r["i"])
                    picchi[f"{n}_{g}"] = max(picchi[f"{n}_{g}"], float(np.abs(r["i_cmd"]).max()))
                    r_cl.append([asse, a, w.nome[-2:], n, g, dt.fit(y, w.pos), dt.fit(y[mv], w.pos[mv]), dt.fit(im[mv], w.I[mv]),
                                 float(np.rad2deg(np.max(np.abs(w.pos - w.cmd)[mv])) * 1e3), float(np.rad2deg(np.max(np.abs(y - (w.cmd - w.cmd[0]))[mv])) * 1e3),
                                 float(np.abs(w.I[mv]).max()), float(np.abs(im[mv]).max())])
                    c_cl[g][key][1][n] = (y, im)
                print(asse, key, flush=True)
        sommario[chiave] = {"maxwell": c, "gains": {"maxwell_estimated": g_max, "dahl_estimated": g_dahl, "files": {x: getattr(p_file, x) for x in ("Kpp", "Kpv", "Kiv_hz")},
                                                    "start": aac.SETTEMBRE[asse], "cost_maxwell": costo},
                            "max_abs_pole": stab, "peak_i_cmd_A": picchi, "i_limit_A": p_file.i_limit}
        pt.salva(fig_aperto(nome, c_ol, seg), CARTELLA / "figures", f"{chiave}_open_loop")
        pt.salva(fig_chiuso(nome, c_cl["estimated"], "gains estimated on the 1 deg doublet for each bearing model", seg), CARTELLA / "figures", f"{chiave}_closed_loop_estimated")
        pt.salva(fig_chiuso(nome, c_cl["akd"], "gains of the .akd and WorkBench files", seg), CARTELLA / "figures", f"{chiave}_closed_loop_akd")
        plt.close("all")
    with open(CARTELLA / "open_loop_fit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["axis", "test", "move", "bearings", "fit_window_pct", "fit_motion_pct"])
        wr.writerows([r[:4] + [f"{x:.2f}" for x in r[4:]] for r in r_ol])
    with open(CARTELLA / "closed_loop_fit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["axis", "test", "move", "bearings", "gains", "fit_position_window_pct", "fit_position_motion_pct", "fit_current_motion_pct",
                     "max_following_error_meas_mdeg", "max_following_error_model_mdeg", "max_current_meas_A", "max_current_model_A"])
        wr.writerows([r[:5] + [f"{x:.2f}" for x in r[5:]] for r in r_cl])

    def med(rr, i, **f):
        x = [q[i] for q in rr if tuple(q[:3]) != ANOMALA and all(q[{"asse": 0, "n": 3, "g": 4}[k]] == v for k, v in f.items())]
        return float(np.median(x))
    tab = {}
    for asse in ("azimut", "polare"):
        tab[asse] = {"open_loop": {n: {"fit_motion": med(r_ol, 5, asse=asse, n=n)} for n in ("dahl", "maxwell")},
                     "closed_loop": {g: {n: {"fit_position_window": med(r_cl, 5, asse=asse, n=n, g=g), "fit_position_motion": med(r_cl, 6, asse=asse, n=n, g=g),
                                             "fit_current_motion": med(r_cl, 7, asse=asse, n=n, g=g), "err_max_model": med(r_cl, 9, asse=asse, n=n, g=g),
                                             "i_max_model": med(r_cl, 11, asse=asse, n=n, g=g)} for n in ns} for g, ns in (("estimated", ("dahl", "maxwell")), ("akd", ("dahl", "maxwell")), ("dahl_gains", ("maxwell",)))},
                     "measured": {"err_max": med(r_cl, 8, asse=asse, n="dahl", g="estimated"), "i_max": med(r_cl, 10, asse=asse, n="dahl", g="estimated")}}
    sommario["medians_without_anomalous_window"] = tab
    (CARTELLA / "summary.json").write_text(json.dumps(sommario, indent=1))
    print(json.dumps(sommario, indent=0))


if __name__ == "__main__":
    main()
