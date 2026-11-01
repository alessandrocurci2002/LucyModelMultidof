"""Doppietti V2 in anello chiuso: comando di posizione misurato -> drive AKD2G -> motore -> modello ibrido.

Richiesta dell'utente (2026-10-31): dare in ingresso il riferimento di posizione del controllore e lasciare che la
corrente venga dall'anello di controllo, per capire se il comportamento del sistema analitico e di quello
sperimentale e' simile.

Anello di un asse (l'altro e' fermo nei doppietti V2):

* drive: modello AKD2G di ``closedLoopModel/lucy_cl`` (anelli di posizione, velocita' e corrente, biquad, osservatore)
  con i parametri ``Parameters_*.akd`` e filtri dell'export WorkBench della sessione V2
  (``VClosedLoop/V2/BodeMeasurement_*.csv``); senza saturazioni il drive e' lineare e se ne ricavano le matrici da
  ``Drive.step``; dopo la simulazione si controlla che i limiti (corrente, integratore, tensione) non siano raggiunti;
* guadagni Kpp, Kpv, Kiv: quelli dei file (``--guadagni akd``) non riproducono le sequenze del PLC (errore di
  inseguimento ~1.3 gradi a 10 gradi/s su entrambi gli assi: anello di posizione efficace ~7.7 1/s, come trovato in
  closedLoopModel il 2026-10-19, dove per decisione dell'utente i guadagni si stimano dai dati); con ``--guadagni
  stima`` (default) si stimano Kpp, Kpv, Kiv per asse sull'errore di inseguimento delle prove da 1 grado, con il plant
  fisso, partendo dai valori di closedLoopModel; la corrente non entra nella stima ed e' il confronto indipendente;
* motore: TBM(S)-7646-A (``data/motor_TBM7646A.yaml``, a freddo), corrente q per fase rms, f.c.e.m.; coppia
  g_Kt K_t i con il g_Kt del modello;
* plant: struttura modale della stima ibrida (``results_fork_modes_dahl/time_domain_bearings``) senza cuscinetti, piu'
  l'attrito di Dahl puro sui cuscinetti (sigma0 dalle FRF, T_c dai doppietti).

La parte lineare (struttura, motore, drive) e' discretizzata esattamente (ZOH, passo 0.1 ms) con ingressi comando di
posizione e coppia d'attrito; lo stato di Dahl e' aggiornato esattamente a velocita' costante sul passo. Stato
iniziale: equilibrio dell'anello con l'attrito che regge la corrente di tenuta misurata (limitata a T_c). Si
confrontano posizione dell'encoder e corrente (misura: AMPS / 2, la convenzione di closedLoopModel; scala non
documentata, vedi il report).

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_anello_chiuso_ibrido.py
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))
sys.path.insert(0, str(RADICE / "closedLoopModel"))

import matplotlib.pyplot as plt  # noqa: E402
from scipy.linalg import expm, solve  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import doppietti_tempo as dt  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl.driver import Drive, DriveParams  # noqa: E402
from lucy_cl.motor import MotorParams  # noqa: E402

IBRIDO = RADICE / "results_fork_modes_dahl" / "time_domain_bearings" / "summary.json"
VCL = RADICE.parent / "LucyIdentification" / "VClosedLoop"
CARTELLA = RADICE / "results_fork_modes_dahl" / "closed_loop_check"
PASSO = 1e-4
PROVE = {"azimut": "MPCLucy_doublets_azimut_%s.csv", "polare": "MPCLucy_doublets_elevation_%s.csv"}
DRIVE = {"azimut": ("Parameters_azimuth.akd", "V2/BodeMeasurement_azimuth.csv"),
         "polare": ("Parameters_elevation.akd", "V2/BodeMeasurement_elevation.csv")}
AMPIEZZE = ("01deg_minus01deg", "05deg_minus05deg", "1deg_minus1deg", "3deg_minus3deg", "5deg_minus5deg",
            "8deg_minus8deg")
FIGURA = (("01deg_minus01deg", "#1"), ("1deg_minus1deg", "#1"), ("5deg_minus5deg", "#2"))
SETTEMBRE = {"azimut": {"Kpp": 7.43, "Kpv": 198.0, "Kiv_hz": 0.148}, "polare": {"Kpp": 8.58, "Kpv": 61.4, "Kiv_hz": 0.73}}
"""Guadagni efficaci stimati da closedLoopModel (validate.py, run 3, sequenze di settembre): partenza della stima."""
STIMA = "1deg_minus1deg"
PASSO_STIMA = 2e-4


def stima_guadagni(m, p0: DriveParams, mot, att, finestre_stima) -> tuple[dict[str, float], float]:
    """Kpp, Kpv, Kiv_hz sull'errore di inseguimento (residuo normalizzato come il fit NRMSE di ogni finestra)."""
    nomi = ("Kpp", "Kpv", "Kiv_hz")
    dati = []
    for w in finestre_stima:
        vc = np.abs(np.gradient(w.cmd, dt.DT)) > 1e-6
        i0 = np.flatnonzero(vc)[0]
        k = slice(max(i0 - int(0.3 / dt.DT), 0), min(i0 + int(3.0 / dt.DT), len(w.t)))
        tt = w.t[k] - w.t[k][0]
        t = np.arange(0.0, tt[-1], PASSO_STIMA)
        err = (w.cmd - w.pos)[k] - (w.cmd - w.pos)[k][0]
        dati.append((t, np.interp(t, tt, (w.cmd - w.cmd[0])[k] - (w.cmd - w.cmd[0])[k][0]),
                     float(w.I[k][:int(0.25 / dt.DT)].mean()), tt, err))

    def res(x):
        p = dataclasses.replace(p0, **dict(zip(nomi, np.exp(x))))
        an = AnelloChiuso(m, Drive(p), mot, att, PASSO_STIMA)
        out = []
        for t, pc, I0, tt, err in dati:
            r = an.corri(pc, I0)
            e_mod = np.interp(tt, t, pc - r["p"])
            out.append((e_mod - err)[::4] / np.linalg.norm(err - err.mean()) * np.sqrt(4))
        return np.concatenate(out)
    x0 = np.log([getattr(p0, k) for k in nomi])
    r = least_squares(res, x0, bounds=(np.log([0.3, 1.0, 0.02]), np.log([200.0, 2000.0, 60.0])), x_scale="jac",
                      diff_step=1e-3, max_nfev=60)
    return dict(zip(nomi, map(float, np.exp(r.x)))), float(r.cost)


class AnelloChiuso:
    """Asse in anello chiuso: stato [x, xdot, i, drive]; ingressi [p_cmd, T_b]; attrito di Dahl esterno."""

    def __init__(self, m, drive: Drive, motore: MotorParams, att: dt.Attrito, passo: float = PASSO):
        n, nd = len(m.M), drive.n_states
        self.m, self.drive, self.mot, self.att = m, drive, motore, att
        # drive lineare: d xd/dt = A xd + B u, v_phase = Cv xd + Dv u, i_cmd = Ci xd + Di u; u = [p_cmd, p_fb, v_fb, i]
        x0, u0 = np.zeros(nd), np.zeros(4)

        def f(xd, u):
            v, ic, dx = drive.step(xd, u[0], 0.0, u[1], u[2], u[3])
            return np.r_[dx, v, ic]
        y0 = f(x0, u0)
        h = 1e-6
        Jx = np.column_stack([(f(x0 + h * np.eye(nd)[k], u0) - y0) / h for k in range(nd)])
        Ju = np.column_stack([(f(x0, u0 + h * np.eye(4)[k]) - y0) / h for k in range(4)])
        Ad, Bd, Cv, Dv, Ci, Di = Jx[:nd], Ju[:nd], Jx[nd], Ju[nd], Jx[nd + 1], Ju[nd + 1]
        # struttura
        C = m.M @ m.Phi @ np.diag(2 * m.zeta * m.omega) @ m.Phi.T @ m.M
        Mi = np.linalg.inv(m.M)
        N = 2 * n + 1 + nd
        ix, iv, ii, idr = slice(0, n), slice(n, 2 * n), 2 * n, slice(2 * n + 1, N)
        A = np.zeros((N, N)); B = np.zeros((N, 2))
        A[ix, iv] = np.eye(n)
        A[iv, ix] = -Mi @ m.K
        A[iv, iv] = -Mi @ C - motore.b_viscous * np.outer(Mi @ m.e, m.e)
        A[iv, ii] = Mi @ m.e * m.Kt
        B[iv, 1] = -Mi @ m.v
        # ingressi del drive in funzione dello stato: p_fb = c x, v_fb = c xdot, i
        G = np.zeros((4, N)); G[1, ix] = m.c; G[2, iv] = m.c; G[3, ii] = 1.0
        g = np.zeros((4, 2)); g[0, 0] = 1.0
        A[idr, :] += Bd @ G; A[idr, idr] += Ad; B[idr, :] += Bd @ g
        # motore: L di/dt = v_phase - R i - Ke omega_rel, omega_rel = e^T xdot (rotore rispetto allo statore)
        rv = np.zeros(N); rv[idr] = Cv; rv += Dv @ G; rv[ii] -= motore.R_ph; rv[iv] -= motore.Ke_ph * m.e
        A[ii, :] = rv / motore.L_ph
        B[ii, :] = (Dv @ g) / motore.L_ph
        self.A, self.B = A, B
        self.uscite = {"p": np.r_[m.c, np.zeros(N - n)], "i": np.eye(N)[ii]}
        self.ic_row = np.zeros(N); self.ic_row[idr] = Ci; self.ic_row += Di @ G
        self.ic_in = Di @ g
        self.vp_row = rv + motore.R_ph * np.eye(N)[ii] + motore.Ke_ph * np.r_[np.zeros(n), m.e, np.zeros(N - 2 * n)]
        self.vp_in = Dv @ g
        self.vint_row = np.zeros(N)
        wi = 2 * np.pi * drive.p.Kiv_hz
        self.vint_row[2 * n + 1 + drive.i_zv] = drive.p.Kpv * wi
        E = expm(np.block([[A, B], [np.zeros((2, N + 2))]]) * passo)
        self.Ad, self.Bd, self.N, self.n, self.passo = E[:N, :N], E[:N, N:], N, n, passo

    def equilibrio(self, Tb: float) -> np.ndarray:
        """Stato a riposo con p_cmd = 0 e coppia d'attrito Tb (il velocity integrator regge Tb / (g K_t))."""
        return solve(self.A, -self.B @ np.array([0.0, Tb]))

    def corri(self, p_cmd: np.ndarray, I0: float) -> dict[str, np.ndarray]:
        att, s0, n, v = self.att, self.att.sigma0, self.n, self.m.v
        Tb0 = float(np.clip(self.m.Kt * I0, -att.Tc, att.Tc))
        z = Tb0 / s0
        s = self.equilibrio(Tb0)
        p0 = self.uscite["p"] @ s
        P, I, IC, VP, VI = (np.empty(len(p_cmd)) for _ in range(5))
        vv = np.r_[np.zeros(n), v, np.zeros(self.N - 2 * n)]
        for k, pc in enumerate(p_cmd):
            u = np.array([pc, s0 * z])
            P[k], I[k] = self.uscite["p"] @ s - p0, self.uscite["i"] @ s
            IC[k], VP[k], VI[k] = self.ic_row @ s + self.ic_in @ u, self.vp_row @ s + self.vp_in @ u, self.vint_row @ s
            s_new = self.Ad @ s + self.Bd @ u
            vd = 0.5 * (vv @ s + vv @ s_new)
            a = s0 * abs(vd) / att.g(vd)
            if a > 0:
                e = np.exp(-a * self.passo)
                z = z * e + vd / a * (1 - e)
            s = s_new
        return {"p": P, "i": I, "i_cmd": IC, "v_phase": VP, "i_int": VI, "Tb0": Tb0}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--guadagni", choices=("stima", "akd"), default="stima")
    guadagni = ap.parse_args().guadagni
    CARTELLA_G = CARTELLA / guadagni
    pt.applica_stile()
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    Kt, seg = par.sperimentale.Kt, par.segnaposto
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    ib = json.loads(IBRIDO.read_text())
    mot = MotorParams.from_yaml(RADICE / "closedLoopModel" / "data" / "motor_TBM7646A.yaml", "cold")
    (CARTELLA_G / "figures").mkdir(parents=True, exist_ok=True)
    righe, limiti, stimati = [], {}, {}
    for asse, chiave, cat, mod in (("azimut", "azimuth", cat_az, fl.azimut), ("polare", "elevation", cat_pol, fl.polare)):
        h = ib[chiave]["hybrid"]
        v = h["parameters"] | {"Tc": h["Tc_from_time"]}
        m = mod(cat, v, Kt, {}).struttura
        akd, wb = DRIVE[asse]
        att = dt.Attrito(v["sigma0"], v["Tc"])
        p_file = DriveParams.from_files(chiave, VCL / akd, VCL / wb, pl_kp_unit="hz")
        if guadagni == "stima":
            g, costo = stima_guadagni(m, dataclasses.replace(p_file, **SETTEMBRE[asse]), mot, att,
                                      dt.finestre(VCL / "V2" / (PROVE[asse] % STIMA), asse))
            stimati[asse] = {"estimated": g, "cost": costo, "start_closedLoopModel": SETTEMBRE[asse],
                             "files": {"Kpp": p_file.Kpp, "Kpv": p_file.Kpv, "Kiv_hz": p_file.Kiv_hz}}
            print(f"{asse}: guadagni stimati {g} (file: Kpp {p_file.Kpp:.2f}, Kpv {p_file.Kpv}, Kiv {p_file.Kiv_hz})",
                  flush=True)
            p_file = dataclasses.replace(p_file, **g)
        drive = Drive(p_file)
        an = AnelloChiuso(m, drive, mot, att)
        print(f"{asse}: poli del'anello lineare (attrito bloccato) stabili: "
              f"{bool(np.all(np.abs(np.linalg.eigvals(an.Ad)) < 1 + 1e-9))}", flush=True)
        curve, picchi = {}, {"i_cmd": 0.0, "v_phase": 0.0, "i_int": 0.0}
        for a in AMPIEZZE:
            for w in dt.finestre(VCL / "V2" / (PROVE[asse] % a), asse):
                t = np.arange(0.0, w.t[-1], PASSO)
                pc = np.interp(t, w.t, w.cmd - w.cmd[0])
                I0 = float(w.I[:int(0.4 / dt.DT)].mean())
                r = an.corri(pc, I0)
                for k in picchi:
                    picchi[k] = max(picchi[k], float(np.abs(r[k]).max()))
                y = np.interp(w.t, t, r["p"])
                i_mod = np.interp(w.t, t, r["i"])
                vc = np.abs(np.gradient(w.cmd, dt.DT)) > 1e-6
                i0, i1 = np.flatnonzero(vc)[[0, -1]]
                mv = slice(i0, min(i1 + int(1.0 / dt.DT), len(y)))
                righe.append([asse, a, w.nome[-2:], dt.fit(y, w.pos), dt.fit(y[mv], w.pos[mv]),
                              dt.fit(i_mod[mv], w.I[mv]),
                              np.rad2deg(np.max(np.abs(w.pos - w.cmd)[mv])) * 1e3,
                              np.rad2deg(np.max(np.abs(y - (w.cmd - w.cmd[0]))[mv])) * 1e3,
                              float(np.abs(w.I[mv]).max()), float(np.abs(i_mod[mv]).max()),
                              I0, r["Tb0"] / m.Kt])
                curve[(a, w.nome[-2:])] = (w, y, i_mod, i0)
        limiti[asse] = picchi | {"i_limit": drive.p.i_limit, "vel_int_limit": drive.p.vel_int_limit,
                                 "V_phase_max": drive.V_phase_max}

        fig, ax = plt.subplots(3, 3, figsize=(15, 9.5), sharex="col", gridspec_kw={"height_ratios": [2, 1.3, 1.6]})
        for j, key in enumerate(FIGURA):
            w, y, im, i0 = curve[key]
            k = (w.t >= w.t[i0] - 0.2) & (w.t <= w.t[i0] + 2.0)
            t = w.t[k]
            md = lambda q: np.rad2deg(q[k]) * 1e3  # noqa: E731
            ax[0, j].plot(t, md(w.cmd), color=pt.INCHIOSTRO, lw=0.8, ls=":", label="position command (model input)")
            ax[0, j].plot(t, md(w.pos), color=pt.INCHIOSTRO, lw=1.6, label="measured encoder")
            ax[0, j].plot(t, md(y), lw=1.4, label="model in closed loop: AKD2G + TBM7646 + hybrid Dahl plant",
                          **pt._stile(1))
            ax[1, j].plot(t, md(w.pos - w.cmd), color=pt.INCHIOSTRO, lw=1.3)
            ax[1, j].plot(t, md(y - (w.cmd - w.cmd[0])), lw=1.2, **pt._stile(1))
            ax[2, j].plot(t, w.I[k], color=pt.INCHIOSTRO, lw=1.1, label="measured (AMPS / 2)")
            ax[2, j].plot(t, im[k], lw=1.1, label="model", **pt._stile(1))
            ax[0, j].set_title(f"{key[0].split('_')[0].replace('01deg', '0.1deg').replace('05deg', '0.5deg')} doublet, "
                               f"move {key[1]}", fontsize=9)
            ax[2, j].set_xlabel("time from window start [s]")
        ax[0, 0].set_ylabel("position [mdeg]"); ax[1, 0].set_ylabel("following error [mdeg]")
        ax[2, 0].set_ylabel("motor current [A rms]")
        ax[0, 0].legend(loc="upper left", fontsize=7); ax[2, 0].legend(fontsize=7)
        nome = "Azimuth" if asse == "azimut" else "Elevation"
        gtxt = ("drive gains of the .akd/WorkBench files" if guadagni == "akd" else
                f"effective gains Kpp {drive.p.Kpp:.2f} 1/s, Kpv {drive.p.Kpv:.0f} A s/rad, Kiv {drive.p.Kiv_hz:.2f} Hz "
                f"estimated on the 1 deg doublet")
        pt._titolo(fig, f"{nome}: closed loop driven by the measured position command ({gtxt})", seg)
        pt.salva(fig, CARTELLA_G / "figures", f"{chiave}_closed_loop_hybrid")
        plt.close(fig)

    with open(CARTELLA_G / "closed_loop_fit.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["axis", "test", "move", "fit_position_window_pct", "fit_position_motion_pct",
                     "fit_current_motion_pct", "max_following_error_meas_mdeg", "max_following_error_model_mdeg",
                     "max_current_meas_A", "max_current_model_A", "holding_current_meas_A", "holding_current_model_A"])
        wr.writerows([r[:3] + [f"{x:.2f}" for x in r[3:]] for r in righe])
    for r in righe:
        print(f"{r[0]:7s} {r[1]:18s} {r[2]}: pos finestra {r[3]:6.1f} %, pos movimento {r[4]:6.1f} %, corrente "
              f"{r[5]:6.1f} % | err. inseg. max mis {r[6]:7.1f} mod {r[7]:7.1f} mdeg | I max mis {r[8]:5.2f} "
              f"mod {r[9]:5.2f} A | I tenuta mis {r[10]:5.2f} mod {r[11]:5.2f} A")
    print("limiti:", json.dumps(limiti, indent=1))
    (CARTELLA_G / "limits.json").write_text(json.dumps(limiti, indent=1))
    if stimati:
        (CARTELLA_G / "estimated_gains.json").write_text(json.dumps(stimati, indent=1))


if __name__ == "__main__":
    main()
