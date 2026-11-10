"""Confronto dell'LQG-PI tarato con i doppietti misurati e con la cascata, indagine sull'azimut e verifiche sull'elevazione.

La misura e' fatta con la cascata reale: sull'LQG non esiste una misura, quindi il confronto passa dalla cascata simulata, che dice quanto
il banco riproduce la misura. Azimut: catena di controllori dal P sull'encoder all'LQG finale, poi modello interno, V e q_i con q_p
rialzato fino alla corrente. Elevazione: rumore e quantizzazione dell'encoder, ritardo, K_t, plant dalla FRF misurata.

    .venv/bin/python -I analisi_lqg_confronto_misure.py [--rumore]

``--rumore`` ripete solo la prova con rumore assunto, partendo dai risultati gia' salvati: la posizione misurata a riposo e' costante, quindi il
rumore non si puo' stimare dai dati.

Risultati in results_lqg_tuning/comparison_measured/.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE)]

import matplotlib.pyplot as plt  # noqa: E402

import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_passi_2_4 as p24  # noqa: E402
import analisi_lqg_passo1_flessibili as p1  # noqa: E402
from gimbal_mdof import doppietti_tempo as dt  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402

CARTELLA = RADICE / "results_lqg_tuning" / "comparison_measured"
LSB = 2 * np.pi / 2 ** 26
"""FB1.BITS = 26 nei file .akd."""
RHO_AZ = (1.0, 1e2, 1e4, 1e6)
QI_AZ = (0.1, 41.88)
FATTORI = (1.0, 1.3, 1.7, 2.2, 3.0, 4.0, 5.5, 7.5, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0)
K_P_EL = (1e2, 1e3, 1e4, 1e5, 1e6)
SEME = 1
SIGMA_ASSUNTE = (1e-6, 1e-5)
"""Rumore bianco assunto sulla misura, rad."""


class PEncoder:
    """Solo proporzionale sull'encoder, con il ritardo del master: l'equivalente del position loop della cascata."""

    def __init__(self, k):
        self.k = k

    def stato_iniziale(self):
        return {"y": 0.0}

    def uscita(self, st, r):
        return float(self.k * (r - st["y"]))

    def stima(self, st, y, u):
        st["y"] = y


class ConRumore:
    """Rumore bianco e quantizzazione sulla misura che arriva al controllore."""

    def __init__(self, c, sigma, rng):
        self.c, self.sigma, self.rng = c, sigma, rng

    def stato_iniziale(self):
        return self.c.stato_iniziale()

    def uscita(self, st, r):
        return self.c.uscita(st, r)

    def __getattr__(self, nome):
        # ``misura`` e gli altri metodi del controllore avvolto (senza rumore: la misura del PI e' quella dello stimatore solo se sigma = 0)
        return getattr(self.c, nome)

    def stima(self, st, y, u):
        self.c.stima(st, float(np.round((y + self.sigma * self.rng.standard_normal()) / LSB) * LSB), u)


def simula(a, c, w, ritardo=g.RITARDO):
    if c is None:
        return g.simula(a, None, w)
    t = np.arange(0.0, w.t[-1], g.PASSO_SIM)
    pc, I0 = np.interp(t, w.t, w.cmd - w.cmd[0]), float(w.I[:int(0.4 / dt.DT)].mean())
    r = g.corri_lqg(a, c, pc, I0, ritardo)
    return {n: np.interp(w.t, t, r[n]) for n in ("p", "i")} | {f"max_{n}": float(np.abs(r[n]).max()) for n in ("i_cmd", "v_phase", "i_int")}


def riassunto(ms: list[dict], picchi=None) -> dict:
    return {"i_max_A": float(max(picchi or [m["i_max_A"] for m in ms])), "err_rms_median_mdeg": float(np.median([m["err_rms_mdeg"] for m in ms])), "err_max_median_mdeg": float(np.median([m["err_max_mdeg"] for m in ms])),
            "t_set_median_s": float(np.median([m["t_set_s"] for m in ms])), "t_set_max_s": float(max(m["t_set_s"] for m in ms)), "over_max_mdeg": float(max(m["over_mdeg"] for m in ms)), "i_rms_median_A": float(np.median([m["i_rms_A"] for m in ms]))}


def prova(a, c, W, ritardo=g.RITARDO) -> dict:
    ms, pk = [], []
    for w in W:
        r = simula(a, c, w, ritardo)
        ms.append(g.metriche(w, r["p"], r["i"]))
        pk.append(r["max_i_cmd"])
    return riassunto(ms, pk) | {"per_amplitude": [{"name": w.nome, "err_rms_mdeg": m["err_rms_mdeg"], "err_max_mdeg": m["err_max_mdeg"], "t_set_s": m["t_set_s"], "over_mdeg": m["over_mdeg"], "i_max_A": m["i_max_A"]} for w, m in zip(W, ms)]}


def misurato(W) -> dict:
    ms = [g.metriche(w, w.pos - w.pos[0], w.I) for w in W]
    return riassunto(ms) | {"per_amplitude": [{"name": w.nome, "err_rms_mdeg": m["err_rms_mdeg"], "err_max_mdeg": m["err_max_mdeg"], "t_set_s": m["t_set_s"], "over_mdeg": m["over_mdeg"], "i_max_A": m["i_max_A"]} for w, m in zip(W, ms)]}


def rumore_a_riposo(W) -> dict:
    """Scarto della posizione e della corrente misurate prima del moto, sui 0.4 s di ogni finestra."""
    n = int(0.4 / dt.DT)
    return {"position_std_urad": float(np.median([1e6 * np.std(np.diff(w.pos[:n])) / np.sqrt(2) for w in W])), "current_std_A": float(np.median([np.std(w.I[:n]) for w in W])), "sample_s": float(dt.DT)}


def sintesi(A, q_p, q_i, f_modo=None, q_f=None, rho=1.0, base=p1.Q_F):
    q = np.full(len(A.Ap) - 1, float(base))
    if q_f:
        m = next(m for m in A.modi if abs(m.f_hz - f_modo) < 3.0)
        q[[i - 1 for i in m.indici]] = q_f
    return lqg.sintesi(A.Ap, A.Bp, A.Cp, g.PASSO_MASTER, q_i, q_p, q, rho_e=rho, v_continua=True)


def poli_lenti(a, c) -> list:
    an = a.anello_v
    z = p1.poli_chiusi(c, g.varianti_cuscinetto(a)["sliding"], an.B[:, 0], an.uscite["p"])
    l = np.log(z.astype(complex)) / g.PASSO_MASTER
    l = l[(l.imag > 0) & (np.abs(l) / (2 * np.pi) < 20.0)]
    return [{"f_Hz": float(abs(q) / (2 * np.pi)), "zeta": float(-q.real / abs(q))} for q in sorted(l, key=lambda q: -q.real / abs(q))]


def banda_piante(a, c, Pm_f, Pm) -> dict:
    an = a.anello_v
    C, B = an.uscite["p"], an.B[:, 0]
    out = {}
    P = g.frf(g.varianti_cuscinetto(a)["sliding"], B, C, p1.F)
    L, N = g.anello_lqg(c, P, p1.F)
    out["bandwidth_sliding_Hz"] = g.banda(p1.F, N / (1 + L))
    out["modulus_sliding"] = float(np.abs(1 + L).min())
    Lm, Nm = g.anello_lqg(c, Pm, Pm_f)
    out["T_ry_measured_max_dB"] = float(pt._db(Nm / (1 + Lm)).max())
    out["modulus_measured"] = float(np.abs(1 + Lm).min())
    out["modulus_measured_Hz"] = float(Pm_f[np.argmin(np.abs(1 + Lm))])
    out["L_measured_max_dB"] = float(pt._db(Lm).max())
    return out


def catena(a, a0, W, f_modo, fi, Kpp) -> dict:
    """Dal P sull'encoder all'LQG finale, una cosa alla volta."""
    qp0 = Kpp ** 2
    c = sintesi(a, qp0, p1.Q_I, base=0.0)
    c0 = sintesi(a0, qp0, p1.Q_I, base=0.0)
    qf = float(fi["q_f"]["0"]) if fi["q_f"] else None
    cF = sintesi(a, fi["q_p"], fi["q_i"], f_modo, qf)
    nulla = lambda x, kp=None: dataclasses.replace(x, k_p=kp if kp else x.k_p, k_i=0.0, K_f=0 * x.K_f)
    voci = (("Cascade", None), ("P on the encoder, $k=K_{pp}$", PEncoder(Kpp)), ("P on the estimate, no-bearing model", nulla(c0, Kpp)), ("P on the estimate, Maxwell model", nulla(c, Kpp)),
            ("PI on the estimate, Maxwell model", dataclasses.replace(c, K_f=0 * c.K_f)), ("LQG, Maxwell model, $k_p\\approx K_{pp}$", c), ("LQG, no-bearing model, $k_p\\approx K_{pp}$", c0), ("Final LQG (Maxwell)", cF), ("Final weights, no-bearing model", sintesi(a0, fi["q_p"], fi["q_i"], f_modo, qf)))
    out = {}
    for n, x in voci:
        out[n] = prova(a, x, W)
        if isinstance(x, lqg.ControlloreLQG):
            out[n] |= {"slow_poles_sliding": poli_lenti(a, x), "K_e_position_per_s": float(x.K_e[0, 0] / g.PASSO_MASTER), "k_p": x.k_p, "k_i": x.k_i}
        print(a.nome, n, round(out[n]["i_max_A"], 1), round(out[n]["err_rms_median_mdeg"]), round(out[n]["over_max_mdeg"]), round(out[n]["t_set_median_s"], 2), flush=True)
    return out


def configurazioni(a, a0, W, Pm_f, Pm, Kpp) -> list[dict]:
    """Modello interno, V e q_i; per ciascuna q_p sale finche' la corrente sui comandi tipici resta sotto 18 A."""
    righe = []
    for et, A in (("Maxwell", a), ("no bearings", a0)):
        for rho in RHO_AZ:
            for qi in QI_AZ:
                ultimo = None
                for fat in FATTORI:
                    qp = max(Kpp ** 2 / 2, 10.0) * fat
                    c = sintesi(A, qp, qi, rho=rho)
                    x = prova(a, c, W)
                    if x["i_max_A"] > p24.I_LIMITE or abs(p1.poli_chiusi(c, *(g.varianti_cuscinetto(a)["sliding"], a.anello_v.B[:, 0], a.anello_v.uscite["p"]))).max() > 1 + 1e-9:
                        break
                    ultimo = {"internal_model": et, "rho": rho, "q_i": qi, "q_p": qp, "k_p": c.k_p, "k_i": c.k_i, "K_e_position_per_s": float(c.K_e[0, 0] / g.PASSO_MASTER), "doublets": x, "slow_poles_sliding": poli_lenti(a, c)} | banda_piante(a, c, Pm_f, Pm)
                if ultimo:
                    righe.append(ultimo)
                    print(a.nome, et, "rho", rho, "qi", qi, "qp %.3g" % ultimo["q_p"], "kp %.2f" % ultimo["k_p"], "i %.1f err %.0f over %.0f tset %.2f bw %.2f" % (x["i_max_A"] if False else ultimo["doublets"]["i_max_A"], ultimo["doublets"]["err_rms_median_mdeg"], ultimo["doublets"]["over_max_mdeg"], ultimo["doublets"]["t_set_median_s"], ultimo["bandwidth_sliding_Hz"]), flush=True)
    return righe


def verifiche(asse, ch, cat, mod, Kt, S, T, mot, mx, a, W, fi, f_modo, rum, Pm_f, Pm) -> list[dict]:
    """Stesso controllore: nominale, rumore e quantizzazione, un ciclo di ritardo in piu', K_t +-10 %. Anche a guadagni piu' bassi."""
    righe = []
    qf = float(fi["q_f"]["0"]) if fi["q_f"] else None
    punti = [(float(q), fi["q_i"], qf if q >= 1e5 else None) for q in K_P_EL if q < fi["q_p"]] + [(fi["q_p"], fi["q_i"], qf)]
    rng = np.random.default_rng(SEME)
    for qp, qi, q_f in punti:
        c = sintesi(a, qp, qi, f_modo, q_f)
        r = {"q_p": qp, "k_p": c.k_p, "k_i": c.k_i, "K_e_position_per_s": float(c.K_e[0, 0] / g.PASSO_MASTER), "nominal": prova(a, c, W)} | banda_piante(a, c, Pm_f, Pm)
        r["quantisation_only"] = prova(a, ConRumore(c, 0.0, rng), W)
        r["noise_rest"] = prova(a, ConRumore(c, rum["position_std_urad"] * 1e-6, rng), W)
        r["one_more_cycle"] = prova(a, c, W, g.RITARDO + 1)
        r["two_more_cycles"] = prova(a, c, W, g.RITARDO + 2)
        for gs in (0.9, 1.1):
            b = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot, gs)
            r[f"g_Kt_{gs}"] = prova(b, c, W)
            r[f"g_Kt_{gs}_stable"] = bool(all(abs(p1.poli_chiusi(c, A, b.anello_v.B[:, 0], b.anello_v.uscite["p"])).max() < 1 + 1e-9 for A in g.varianti_cuscinetto(b).values()))
        righe.append(r)
        print(ch, "qp %.3g kp %.1f nom err %.1f i %.1f; noise err %.1f i %.1f irms %.2f; +1 err %.1f i %.1f; g0.9 %s %.1f; g1.1 %s %.1f; meas mod %.2f Tmax %.1f" % (qp, c.k_p, r["nominal"]["err_rms_median_mdeg"], r["nominal"]["i_max_A"], r["noise_rest"]["err_rms_median_mdeg"], r["noise_rest"]["i_max_A"], r["noise_rest"]["i_rms_median_A"], r["one_more_cycle"]["err_rms_median_mdeg"], r["one_more_cycle"]["i_max_A"], r["g_Kt_0.9_stable"], r["g_Kt_0.9"]["err_rms_median_mdeg"], r["g_Kt_1.1_stable"], r["g_Kt_1.1"]["err_rms_median_mdeg"], r["modulus_measured"], r["T_ry_measured_max_dB"]), flush=True)
    return righe


def fig_misure(D, seg):
    fig, ax = plt.subplots(2, 4, figsize=(17, 7.5), sharex="col")
    for j, d in enumerate(D):
        a, W, c = d["a"], d["W"], d["c"]
        for kk, iw in enumerate((2, 4)):
            w = W[iw]
            q = j * 2 + kk
            e0 = np.rad2deg((w.cmd - w.cmd[0]) - (w.pos - w.pos[0])) * 1e3
            ax[0, q].plot(w.t, e0, color=pt.INCHIOSTRO, lw=1.3, label="measured (real cascade)")
            ax[1, q].plot(w.t, w.I, color=pt.INCHIOSTRO, lw=0.8)
            for n, (et, x) in enumerate((("cascade, model", None), ("final LQG-PI, model", c)), 1):
                r = simula(a, x, w)
                ax[0, q].plot(w.t, np.rad2deg((w.cmd - w.cmd[0]) - r["p"]) * 1e3, lw=1.1, label=et, **pt._stile(n))
                ax[1, q].plot(w.t, r["i"], lw=0.9, **pt._stile(n))
            ax[0, q].set_title(("Azimuth" if a.nome == "azimut" else "Elevation") + ", " + w.nome.split("_", 1)[1], fontsize=9)
            ax[1, q].set_xlabel("time [s]"); ax[0, q].set_xlim(0, 4)
    ax[0, 0].set_ylabel("following error [mdeg]"); ax[1, 0].set_ylabel("current [A]"); ax[0, 0].legend(fontsize=7)
    pt._titolo(fig, "Measured doublets against the model with the cascade and with the final LQG-PI (no measurement exists with LQG)", seg)
    return fig


def fig_azimut(a, W, voci, seg):
    fig, ax = plt.subplots(2, 2, figsize=(13, 7), sharex="col")
    for j, iw in enumerate((2, 3)):
        w = W[iw]
        ax[0, j].plot(w.t, np.rad2deg((w.cmd - w.cmd[0]) - (w.pos - w.pos[0])) * 1e3, color=pt.INCHIOSTRO, lw=1.2, label="measured (real cascade)")
        ax[1, j].plot(w.t, w.I, color=pt.INCHIOSTRO, lw=0.8)
        for n, (et, x) in enumerate(voci, 1):
            r = simula(a, x, w)
            ax[0, j].plot(w.t, np.rad2deg((w.cmd - w.cmd[0]) - r["p"]) * 1e3, lw=1.1, label=et, **pt._stile(n))
            ax[1, j].plot(w.t, r["i"], lw=0.9, **pt._stile(n))
        ax[0, j].set_title(w.nome, fontsize=9); ax[1, j].set_xlabel("time [s]"); ax[0, j].set_xlim(0, 5)
    ax[0, 0].set_ylabel("following error [mdeg]"); ax[1, 0].set_ylabel("current [A]"); ax[0, 0].legend(fontsize=7)
    pt._titolo(fig, "Azimuth on the Dahl plant: the same P gain through different estimators", seg)
    return fig


def fig_elevazione(V, R, cas, mis, seg):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    k = [r["k_p"] for r in V]
    for n, (et, y) in enumerate((("nominal", [r["nominal"] for r in V]), ("noise 1 urad (assumed)", [q["noise"]["1"] for q in R]), ("noise 10 urad (assumed)", [q["noise"]["10"] for q in R]), ("one more cycle of delay", [r["one_more_cycle"] for r in V])), 1):
        ax[0].plot(k, [x["err_rms_median_mdeg"] for x in y], marker="o", ms=4, lw=1.2, label=et, **pt._stile(n))
        ax[1].plot(k, [x["i_max_A"] for x in y], marker="o", ms=4, lw=1.2, **pt._stile(n))
    ax[2].plot(k, [r["modulus_measured"] for r in V], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO)
    ax[0].axhline(cas["err_rms_median_mdeg"], color=pt.INCHIOSTRO_SECONDARIO, lw=0.8, ls=":", label="cascade, model"); ax[0].axhline(mis["err_rms_median_mdeg"], color=pt.INCHIOSTRO, lw=0.8, ls="--", label="measured")
    ax[1].axhline(p24.I_LIMITE, color=pt.INCHIOSTRO, lw=0.7, ls="--"); ax[2].axhline(0.5, color=pt.INCHIOSTRO, lw=0.7, ls="--")
    for q in ax:
        q.set_xscale("log"); q.set_xlabel("k_p [1/s]")
    ax[0].set_yscale("log"); ax[1].set_yscale("log")
    ax[0].set_ylabel("median rms error, doublets [mdeg]"); ax[1].set_ylabel("largest current [A]; dashed: 18 A"); ax[2].set_ylabel("min |1 + L|, measured plant")
    ax[0].legend(fontsize=7)
    pt._titolo(fig, "Elevation: what the result depends on, at growing gains", seg)
    return fig


def con_rumore(a, x, W) -> dict:
    """Errore e corrente con il rumore assunto; tensione di fase e integratore di velocita' senza rumore, con i loro limiti."""
    out = {}
    for sg in SIGMA_ASSUNTE:
        rng = np.random.default_rng(SEME)
        out[f"{1e6 * sg:.0f}"] = prova(a, ConRumore(x, sg, rng), W)
    pk = [simula(a, x, w) for w in W]
    return {"noise": out, "v_phase_max": float(max(r["max_v_phase"] for r in pk)), "v_bus": a.drive.p.v_bus, "vel_integrator_max": float(max(r["max_i_int"] for r in pk)), "vel_integrator_limit": a.drive.p.vel_int_limit}


def solo_rumore() -> None:
    par, Kt, S, T, mot, assi = g.carica()
    pt.applica_stile()
    O = json.loads((RADICE / "results_lqg_tuning" / "operating_point.json").read_text())
    Fi = json.loads((p24.CARTELLA / "summary.json").read_text())
    Q = json.loads((CARTELLA / "summary.json").read_text())
    out = {}
    for asse, ch, cat, mod in assi:
        b = p1.Banco(asse, ch, cat, mod, Kt, S, T, mot, O[ch]["maxwell"])
        a, W, fi, f_modo = b.a, p24.comandi_tipici(asse), Fi[ch]["final"], b.modi[0]["f_Hz"]
        a0 = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        qf = float(fi["q_f"]["0"]) if fi["q_f"] else None
        r = {"P_encoder_K_pp": con_rumore(a, PEncoder(a.drive.p.Kpp), W)}
        r["final"] = con_rumore(a, sintesi(a, fi["q_p"], fi["q_i"], f_modo, qf), W)
        if asse == "azimut":
            x = next(q for q in Q[ch]["configurations"] if q["internal_model"] == "no bearings" and q["rho"] == 1.0 and q["q_i"] == 0.1)
            r["no_bearing_internal"] = con_rumore(a, sintesi(a0, x["q_p"], x["q_i"]), W)
        else:
            r["ladder"] = []
            for v in Q[ch]["checks"]:
                q_f = qf if v["q_p"] >= 1e5 else None
                r["ladder"].append({"q_p": v["q_p"], "k_p": v["k_p"]} | con_rumore(a, sintesi(a, v["q_p"], fi["q_i"], f_modo, q_f), W))
                print(ch, "qp %.3g" % v["q_p"], {n: (round(y["err_rms_median_mdeg"], 1), round(y["i_max_A"], 1), round(y["i_rms_median_A"], 2)) for n, y in r["ladder"][-1]["noise"].items()}, round(r["ladder"][-1]["v_phase_max"], 1), flush=True)
            fe = fig_elevazione(Q[ch]["checks"], r["ladder"], Q[ch]["chain"]["Cascade"], Q[ch]["measured"], par.segnaposto)
            pt.salva(fe, CARTELLA / "figures", "elevation_checks"); plt.close(fe)
        out[ch] = r
        print(ch, {k: ({n: (round(y["err_rms_median_mdeg"], 1), round(y["i_max_A"], 1), round(y["i_rms_median_A"], 2)) for n, y in v["noise"].items()}, round(v["v_phase_max"], 1), round(v["vel_integrator_max"], 1)) for k, v in r.items() if k != "ladder"}, flush=True)
    (CARTELLA / "noise_check.json").write_text(json.dumps(p24.senza_array(out), indent=1))


def main() -> None:
    if "--rumore" in sys.argv:
        solo_rumore()
        return
    par, Kt, S, T, mot, assi = g.carica()
    seg = par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    O = json.loads((RADICE / "results_lqg_tuning" / "operating_point.json").read_text())
    Fi = json.loads((p24.CARTELLA / "summary.json").read_text())
    out, D = {}, []
    for asse, ch, cat, mod in assi:
        mx = O[ch]["maxwell"]
        b = p1.Banco(asse, ch, cat, mod, Kt, S, T, mot, mx)
        a, W = b.a, p24.comandi_tipici(asse)
        a0 = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot)
        fi, Kpp, f_modo = Fi[ch]["final"], a.drive.p.Kpp, b.modi[0]["f_Hz"]
        mis = g.misura_asse(par.sperimentale, asse)
        k = (mis.f >= 1.0) & (mis.f <= 150.0)
        Pm_f, Pm = mis.f[k], g.plant_misurato(a, mis.f[k], mis.G[k])
        rum = rumore_a_riposo(W)
        cF = sintesi(a, fi["q_p"], fi["q_i"], f_modo, float(fi["q_f"]["0"]) if fi["q_f"] else None)
        r = {"measured": misurato(W), "rest_noise": rum, "chain": catena(a, a0, W, f_modo, fi, Kpp), "final_on_plants": banda_piante(a, cF, Pm_f, Pm)}
        r["P_encoder_ladder"] = [{"k": x} | prova(a, PEncoder(x), W) for x in (Kpp, 10.0, 12.0, 15.0, 20.0, 30.0)]
        D.append({"a": a, "W": W, "c": cF})
        if asse == "azimut":
            r["configurations"] = configurazioni(a, a0, W, Pm_f, Pm, Kpp)
            c0 = sintesi(a0, Kpp ** 2, p1.Q_I, base=0.0)
            c1 = sintesi(a, Kpp ** 2, p1.Q_I, base=0.0)
            fa = fig_azimut(a, W, (("cascade, model", None), ("P on the estimate, no-bearing model", dataclasses.replace(c0, k_p=Kpp, k_i=0.0, K_f=0 * c0.K_f)), ("P on the estimate, Maxwell model", dataclasses.replace(c1, k_p=Kpp, k_i=0.0, K_f=0 * c1.K_f)), ("final LQG-PI", cF)), seg)
        else:
            r["checks"] = verifiche(asse, ch, cat, mod, Kt, S, T, mot, mx, a, W, fi, f_modo, rum, Pm_f, Pm)
        out[ch] = r
        (CARTELLA / "summary.json").write_text(json.dumps(p24.senza_array(out), indent=1))
    salva = lambda f, n: (pt.salva(f, CARTELLA / "figures", n), plt.close(f))
    salva(fig_misure(D, seg), "measured_doublets")
    salva(fa, "azimuth_chain")
    solo_rumore()


if __name__ == "__main__":
    main()
