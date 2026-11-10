"""Passo 1 della taratura dell'LQG-PI (Gawronski 2008, par. 8.2), seguito alla lettera: pesi dei modi flessibili.

Funzioni di trasferimento in anello chiuso e gradino sul solo modello lineare ridotto con Maxwell, analogo del modello ridotto del
libro, con Q = V (8.26) e V intensita' della (8.8). Pesi piccoli su integrale, posizione e su tutti gli stati flessibili; si guarda T_ry
e si alza il peso di un modo solo se il suo picco e' troppo alto. Il riferimento per "troppo alto" e' la Fig. 8.5 del libro. Il plant
Dahl e i plant con i cuscinetti non entrano nella taratura: sono una verifica finale.

    .venv/bin/python -I analisi_lqg_passo1_flessibili.py [--senza-cuscinetti]

Risultati in results_lqg_tuning/step1_flexible/; con ``--senza-cuscinetti`` il modello interno non ha i cuscinetti, come i modelli del libro, e
i risultati vanno in results_lqg_tuning/no_bearings/step1_flexible/.
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
from scipy.linalg import expm  # noqa: E402
from scipy.signal import butter, sosfiltfilt  # noqa: E402

import analisi_lqg_gawronski as g  # noqa: E402
import analisi_lqg_punto_lavoro as apl  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from lucy_cl import lqg  # noqa: E402

SENZA = "--senza-cuscinetti" in sys.argv
CARTELLA = RADICE / "results_lqg_tuning" / ("no_bearings" if SENZA else "") / "step1_flexible"
MODELLO = "reduced no-bearing model" if SENZA else "reduced Maxwell model"
Q_I, Q_P = 0.1, 0.3
"""Primo tentativo del libro, par. 8.4.2 punto 2."""
Q_F = 0.1
"""Peso piccolo su ogni stato flessibile: minore di q_p e q_i, come chiede la (8.22)."""
UNIFORMI = (0.0, 0.01, 0.1, 0.3, 1.0, 5.0, 10.0, 100.0)
"""Stesso peso su tutti gli stati flessibili, come Q_f = 5 I in Fig. 8.5."""
UNIFORMI_FIG = (0.1, 1.0, 10.0, 100.0)
SCALA_MODO = (0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 3e4, 1e5)
"""Peso del solo modo piu' evidente, gli altri a Q_F: mostra cosa accadrebbe alzandolo."""
SCALA_FIG = (0.1, 1e3, 1e4, 1e5)
MOD_MIN = g.MOD_MIN
F = np.geomspace(0.005, 150.0, 3000)
T_GRADINO, T_AVVIO = 20.0, 0.5
GRADINO = np.deg2rad(1.0)
ANTEPRIMA = ((25.0, 4.0), (100.0, 25.0), (400.0, 100.0))
"""q_p, q_i per la relazione pesi-guadagni e per vedere come crescono i picchi; non e' una taratura."""


def poli_chiusi(c: lqg.ControlloreLQG, A, B, C, ritardo: int = g.RITARDO) -> np.ndarray:
    """Poli discreti dell'anello chiuso per un plant qualunque: stati [x, x_hat, e_i, u in coda]."""
    N, n = len(A), len(c.Ad)
    E = expm(np.block([[A, B[:, None]], [np.zeros((1, N + 1))]]) * g.PASSO_MASTER)
    tot = N + n + 1 + ritardo
    M = np.zeros((tot, tot))
    ie, iq = N + n, N + n + 1
    M[:N, :N], M[:N, iq] = E[:N, :N], E[:N, N]
    M[N:N + n, :N] = c.K_e[:, [0]] @ C[None, :]
    M[N:N + n, N:N + n] = c.Ad - c.K_e[:, [0]] @ np.eye(n)[:1]
    M[N:N + n, iq] = c.Bd[:, 0]
    M[ie, ie], M[ie, N] = 1.0, -g.PASSO_MASTER
    M[iq, ie], M[iq, N] = c.k_i, -c.k_p
    if isinstance(c, lqg.ControlloreLQGPIMisurato):  # Fig. 8.4b: il PI legge y = C x invece di y_hat
        M[ie, N], M[iq, N] = 0.0, 0.0
        M[ie, :N], M[iq, :N] = -g.PASSO_MASTER * C, -c.k_p * C
    M[iq, N + 1:N + n] = -c.K_f
    return np.linalg.eigvals(M)


def peggiore(z: np.ndarray) -> dict[str, float]:
    i = int(np.argmax(np.abs(z)))
    return {"abs_z": float(abs(z[i])), "f_Hz": float(abs(np.log(complex(z[i])) / g.PASSO_MASTER) / (2 * np.pi))}


def massimo_locale(f, db):
    """Massimo locale piu' alto in banda, con la sua altezza sul minimo piu' alto dei due lati; ``None`` se la curva e' monotona."""
    k = [i for i in range(1, len(db) - 1) if db[i] > db[i - 1] and db[i] >= db[i + 1]]
    if not k:
        return None
    i = max(k, key=lambda j: db[j])
    return {"f_Hz": float(f[i]), "dB": float(db[i]), "height_dB": float(db[i] - max(db[:i].min(), db[i:].min()))}


class Banco:
    """Un asse, modello ridotto con Maxwell."""

    def __init__(self, asse, chiave, cat, mod, Kt, S, T, mot, mx):
        self.asse, self.chiave, self.mx = asse, chiave, mx
        self.a = a = g.costruisci(asse, chiave, cat, mod, Kt, S, T, mot, cuscinetti=mx)
        self.n = len(a.Ap)
        self.P = np.array([(a.Cp @ np.linalg.solve(2j * np.pi * x * np.eye(self.n) - a.Ap, a.Bp))[0, 0] for x in F])
        self.modi = []
        for m in a.modi:
            if len(m.indici) == 2 and m.zeta < 0.1:
                lo, hi = next((x, y) for x, y in apl.BANDE_RISONANZA[asse] if x <= m.f_hz <= y)
                sg = m.zeta * 2 * np.pi * m.f_hz
                self.modi.append({"f_Hz": m.f_hz, "zeta": m.zeta, "stati": list(m.indici), "banda": (lo, hi), "q_rif": float(4 * sg ** 2 / np.sum(a.Bp[list(m.indici), 0] ** 2))})

    def sintesi(self, q_f: dict[int, float] | None = None, q_i=Q_I, q_p=Q_P, base=Q_F, v_continua=True):
        """Peso ``base`` su ogni stato flessibile; ``q_f``: indice del modo -> peso della sua coppia di stati."""
        q = np.full(self.n - 1, float(base))
        for i, x in (q_f or {}).items():
            q[[j - 1 for j in self.modi[i]["stati"]]] = x
        return lqg.sintesi(self.a.Ap, self.a.Bp, self.a.Cp, g.PASSO_MASTER, q_i, q_p, q, v_continua=v_continua)

    def vettore(self, c_base=Q_F) -> list[dict]:
        """Il vettore q, stato per stato."""
        righe = [{"state": "q_i", "weight": Q_I}, {"state": "q_p", "weight": Q_P}]
        for m in self.a.modi:
            if m.f_hz > 0:
                righe += [{"state": f"{m.f_hz:.2f} Hz" + (f", {k}" if len(m.indici) == 2 else ""), "weight": c_base} for k in (("a", "b") if len(m.indici) == 2 else ("",))]
        return righe

    def tf(self, c):
        L, N = g.anello_lqg(c, self.P, F)
        return N / (1 + L), L

    def valuta(self, c) -> dict:
        T, L = self.tf(c)
        db = pt._db(T)
        lam = np.log(c.poli_regolatore.astype(complex)) / g.PASSO_MASTER
        lam = lam[lam.imag > 0]
        modi = []
        for m in self.modi:
            b = (F >= m["banda"][0]) & (F <= m["banda"][1])
            d = lam[np.argmin(np.abs(np.abs(lam) / (2 * np.pi) - m["f_Hz"]))]
            modi.append({"f_Hz": m["f_Hz"], "band_max_dB": float(db[b].max()), "peak": massimo_locale(F[b], db[b]), "zeta_regulator": float(-d.real / abs(d)),
                         "K_f_norm": float(np.linalg.norm(c.K_f[[j - 1 for j in m["stati"]]])), "K_e_norm": float(np.linalg.norm(c.K_e[m["stati"], 0]))})
        km = int(np.argmin(np.abs(1 + L)))
        return {"k_p": c.k_p, "k_i": c.k_i, "bandwidth_Hz": g.banda(F, T), "modulus": float(np.abs(1 + L[km])), "modulus_Hz": float(F[km]), "T_ry_max_dB": float(db.max()), "modes": modi,
                "K_e_position_per_s": float(c.K_e[0, 0] / g.PASSO_MASTER), "largest_pole": peggiore(poli_chiusi(c, self.a.Ap, self.a.Bp[:, 0], self.a.Cp[0]))}

    def gradino(self, c, durata=T_GRADINO) -> dict:
        """Gradino di 1 grado sul modello ridotto, a passi del master, con lo stesso ritardo di ``g.corri_lqg``."""
        a = self.a
        N = self.n
        E = expm(np.block([[a.Ap, a.Bp], [np.zeros((1, N + 1))]]) * g.PASSO_MASTER)
        Phi, Gam, C = E[:N, :N], E[:N, N], a.Cp[0]
        t = np.arange(0.0, durata, g.PASSO_MASTER)
        x, st, u_app = np.zeros(N), c.stato_iniziale(), 0.0
        y, u = np.empty(len(t)), np.empty(len(t))
        coda = [0.0] * g.RITARDO
        for k, tk in enumerate(t):
            y[k] = C @ x
            if hasattr(c, "misura"):  # aggancio: un controllore con ``misura`` (PI sull'errore misurato, Gawronski Fig. 8.4b) riceve y(k) prima di ``uscita``
                c.misura(st, y[k])
            coda.append(c.uscita(st, GRADINO if tk >= T_AVVIO else 0.0))
            u_app = coda.pop(0)
            u[k] = u_app
            x = Phi @ x + Gam * u_app
            c.stima(st, y[k], u_app)
        return {"t": t, "y": y, "u": u} | metriche_gradino(t, y, 1 / g.PASSO_MASTER)


def metriche_gradino(t, y, fs) -> dict:
    fuori = np.flatnonzero(np.abs(y - GRADINO)[t >= T_AVVIO] > 0.02 * GRADINO)
    vib = sosfiltfilt(butter(4, (30.0, min(150.0, 0.45 * fs)), "bandpass", fs=fs, output="sos"), y) if np.abs(y).max() < 1e3 else np.full_like(y, np.nan)
    return {"vib": vib, "metriche": {"overshoot_pct": float(100 * (y.max() / GRADINO - 1)), "t_set_2pct_s": float(fuori[-1] * (t[1] - t[0])) if fuori.size else 0.0,
                                     "vib_peak_urad": float(1e6 * np.nanmax(np.abs(vib))), "vib_after_100ms_urad": float(1e6 * np.nanmax(np.abs(vib[t >= T_AVVIO + 0.1])))}}


def pesi_guadagni(Kt, S, T, mot, assi, O) -> list[dict]:
    """k_p, k_i da q_p, q_i: formula del plant rigido, dlqr su y' = u, modelli senza cuscinetti e con Maxwell."""
    uno = np.ones((1, 1))
    righe = []
    for qp, qi in ((Q_P, Q_I),) + ANTEPRIMA:
        ki = float(np.sqrt(qi))
        rig = lqg.sintesi(0 * uno, uno, uno, g.PASSO_MASTER, qi, qp, np.zeros(0))
        r = {"q_p": qp, "q_i": qi, "formula": {"k_p": float(np.sqrt(qp + 2 * ki)), "k_i": ki}, "rigid_dlqr": {"k_p": rig.k_p, "k_i": rig.k_i}}
        for asse, ch, cat, mod in assi:
            for et, mx in (("no_bearings", None), ("maxwell", O[ch]["maxwell"])):
                a = g.costruisci(asse, ch, cat, mod, Kt, S, T, mot, cuscinetti=mx)
                c = lqg.sintesi(a.Ap, a.Bp, a.Cp, g.PASSO_MASTER, qi, qp, np.zeros(len(a.Ap) - 1))
                r[f"{ch}_{et}"] = {"k_p": c.k_p, "k_i": c.k_i}
        righe.append(r)
    return righe


class Registra:
    """Controllore che conserva la stima di posizione."""

    def __init__(self, c: lqg.ControlloreLQG):
        self.c, self.y_hat = c, 0.0

    def stato_iniziale(self):
        return self.c.stato_iniziale()

    def uscita(self, st, r):
        self.y_hat = float(st["x"][0])
        return self.c.uscita(st, r)

    def __getattr__(self, nome):
        # ``misura`` e gli altri metodi del controllore avvolto
        return getattr(self.c, nome)

    def stima(self, st, y, u):
        self.c.stima(st, y, u)


def verifica_dahl(b: Banco, c) -> dict:
    """Fuori dalla procedura del libro: poli sul modello completo e sui plant con i cuscinetti, gradino di 1 grado su Dahl (senza limite di corrente)."""
    a, an = b.a, b.a.anello_v
    pl = {"sliding": (g.varianti_cuscinetto(a)["sliding"], an.B[:, 0], an.uscite["p"]), "stuck": (g.varianti_cuscinetto(a)["stuck"], an.B[:, 0], an.uscite["p"])}
    if b.mx:
        pl["full_maxwell"] = lqg.con_cuscinetti(an.A, an.B[:, 0], an.B[:, 1], an.uscite["p"], np.r_[a.m.v, np.zeros(an.N - an.n)], np.r_[np.zeros(an.n), a.m.v, np.zeros(an.N - 2 * an.n)], 0.0, b.mx["k_1"], b.mx["a"])
    t = np.arange(0.0, 120.0, g.PASSO_SIM)
    reg = Registra(c)
    r = g.corri_lqg(a, reg, np.where(t >= T_AVVIO, GRADINO, 0.0), 0.0)
    y = r["p"]
    k = t <= T_GRADINO
    m = metriche_gradino(t[k], y[k], 1 / g.PASSO_SIM)["metriche"]
    return {"largest_pole": {n: peggiore(poli_chiusi(c, *v)) for n, v in pl.items()}, "y_100_120s_mdeg": [float(np.rad2deg(y[t >= 100].min()) * 1e3), float(np.rad2deg(y[t >= 100].max()) * 1e3)],
            "y_minus_y_hat_end_mdeg": float(np.rad2deg(y[-1] - reg.y_hat) * 1e3), "step_20s": m, "i_max_A": float(np.abs(r["i"]).max()), "t": t[::10], "y": y[::10], "i": r["i"][::10]}


def correzione_v(b: Banco) -> list[dict]:
    """Stimatore con V per campione (come usato prima) e con V della (8.8), alla stessa q."""
    out = []
    for x in (0.0, 1e4):
        for et, vc in (("V per sample", False), ("V of (8.8)", True)):
            c = b.sintesi({0: x} if x else None, base=0.0, v_continua=vc)
            an = b.a.anello_v
            ps = {n: peggiore(poli_chiusi(c, A, an.B[:, 0], an.uscite["p"])) for n, A in g.varianti_cuscinetto(b.a).items()}
            out.append({"q_f": x, "estimator": et, "K_e_position_per_s": float(c.K_e[0, 0] / g.PASSO_MASTER), "largest_pole": ps})
    return out


def nome(b):
    return "Azimuth" if b.asse == "azimut" else "Elevation"


def fig_primo(B, C, G, seg):
    fig, ax = plt.subplots(3, 2, figsize=(13, 9.5))
    for j, (b, c, r) in enumerate(zip(B, C, G)):
        T, _ = b.tf(c)
        k = F >= 30.0
        ax[0, j].plot(F, pt._db(T), lw=1.4, color=pt.INCHIOSTRO); ax[0, j].set_xscale("log"); ax[0, j].set_xlabel("frequency [Hz]")
        ax[0, j].axhline(0.0, color=pt.INCHIOSTRO_SECONDARIO, lw=0.8, ls=":"); ax[0, j].axhline(-6.0, color=pt.INCHIOSTRO_SECONDARIO, lw=0.8, ls="--")
        ax[1, j].plot(F[k], pt._db(T[k]), lw=1.4, color=pt.INCHIOSTRO); ax[1, j].set_xscale("log"); ax[1, j].set_xlabel("frequency [Hz]")
        for m in b.modi:
            ax[1, j].axvline(m["f_Hz"], color=pt.INCHIOSTRO_SECONDARIO, lw=0.6, ls=":")
        ax[2, j].plot(r["t"], np.rad2deg(np.where(r["t"] >= T_AVVIO, GRADINO, 0.0)), color=pt.INCHIOSTRO_SECONDARIO, lw=0.9, ls=":", label="command")
        ax[2, j].plot(r["t"], np.rad2deg(r["y"]), lw=1.4, color=pt.INCHIOSTRO, label="encoder y"); ax[2, j].set_xlabel("time [s]")
        ax[0, j].set_title(f"{nome(b)}: k_p = {c.k_p:.2f}, k_i = {c.k_i:.2f}", fontsize=9)
    ax[0, 0].set_ylabel("T_ry [dB]; lines: 0 and -6 dB, peaks of Fig. 8.5"); ax[1, 0].set_ylabel("T_ry, zoom on the modes [dB]"); ax[2, 0].set_ylabel("step response [deg]"); ax[2, 0].legend(fontsize=7)
    pt._titolo(fig, f"First guess on the {MODELLO}: q_i = 0.1, q_p = 0.3, q_f = 0.1 on all flexible states, Q = V", seg)
    return fig


def fig_uniformi(B, seg):
    fig, ax = plt.subplots(3, 2, figsize=(13, 9.5))
    for j, b in enumerate(B):
        k = F >= 30.0
        for n, x in enumerate(UNIFORMI_FIG, 1):
            c = b.sintesi(base=x)
            T = b.tf(c)[0]
            r = b.gradino(c)
            ax[0, j].plot(F, pt._db(T), lw=1.2, label=f"q_f = {x:g}", **pt._stile(n))
            ax[1, j].plot(F[k], pt._db(T[k]), lw=1.1, **pt._stile(n))
            ax[2, j].plot(r["t"], np.rad2deg(r["y"]), lw=1.2, **pt._stile(n))
        ax[0, j].set_xscale("log"); ax[1, j].set_xscale("log"); ax[0, j].set_xlabel("frequency [Hz]"); ax[1, j].set_xlabel("frequency [Hz]"); ax[2, j].set_xlabel("time [s]")
        ax[0, j].set_title(nome(b), fontsize=9); ax[0, j].legend(fontsize=7)
    ax[0, 0].set_ylabel("T_ry [dB]"); ax[1, 0].set_ylabel("T_ry, zoom on the modes [dB]"); ax[2, 0].set_ylabel("step response [deg]")
    pt._titolo(fig, f"Process on the {MODELLO}: the same weight on all flexible states", seg)
    return fig


def fig_scala(B, R, seg):
    fig, ax = plt.subplots(2, 2, figsize=(13, 7.5), sharex="col")
    for j, (b, r) in enumerate(zip(B, R)):
        im = r["ordine"][0]
        k = F >= 30.0
        for n, x in enumerate(SCALA_FIG, 1):
            c = b.sintesi({im: x})
            T = b.tf(c)[0]
            ax[0, j].plot(F[k], pt._db(T[k]), lw=1.2, label=f"q_f = {x:g}", **pt._stile(n))
            ax[1, j].plot(b.gradino(c)["t"], np.rad2deg(b.gradino(c)["y"]), lw=1.1, **pt._stile(n))
        ax[0, j].axvline(b.modi[im]["f_Hz"], color=pt.INCHIOSTRO_SECONDARIO, lw=0.6, ls=":")
        ax[0, j].set_xscale("log"); ax[1, j].set_xscale("linear"); ax[0, j].set_xlabel("frequency [Hz]"); ax[1, j].set_xlabel("time [s]"); ax[0, j].legend(fontsize=7)
        ax[1, j].tick_params(labelbottom=True)
        ax[0, j].set_title(f"{nome(b)}: weight of the {b.modi[im]['f_Hz']:.1f} Hz mode, others at {Q_F:g}", fontsize=9)
    fig.subplots_adjust(hspace=0.3)
    ax[0, 0].set_ylabel("T_ry [dB]"); ax[1, 0].set_ylabel("step response [deg]")
    pt._titolo(fig, f"What raising one flexible weight would do, on the {MODELLO}", seg)
    return fig


def fig_traccia(B, R, seg):
    fig, ax = plt.subplots(4, 2, figsize=(13, 10), sharex="col")
    x = np.array(SCALA_MODO)
    for j, (b, r) in enumerate(zip(B, R)):
        im = r["ordine"][0]
        s = r["scala"]
        m = [v["modes"][im] for v in s]
        ax[0, j].plot(x, [q["band_max_dB"] for q in m], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO)
        ax[1, j].plot(x, [q["zeta_regulator"] for q in m], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[1, j].set_yscale("log")
        ax[2, j].plot(x, [v["modulus"] for v in s], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[2, j].axhline(MOD_MIN, color=pt.INCHIOSTRO, lw=0.7, ls="--")
        ax[3, j].plot(x, [v["bandwidth_Hz"] for v in s], marker="o", ms=4, lw=1.2, color=pt.INCHIOSTRO); ax[3, j].set_ylim(0, 1.3 * r["primo"]["bandwidth_Hz"])
        ax[3, j].set_xscale("log"); ax[3, j].set_xlabel("q_f of the most evident mode, others at 0.1")
        ax[0, j].set_title(f"{nome(b)}, {b.modi[im]['f_Hz']:.1f} Hz mode", fontsize=9)
    ax[0, 0].set_ylabel("largest T_ry in the band [dB]"); ax[1, 0].set_ylabel("damping of the mode"); ax[2, 0].set_ylabel("min |1 + L|"); ax[3, 0].set_ylabel("bandwidth [Hz]")
    pt._titolo(fig, f"Effect of the weight of one mode on the {MODELLO}", seg)
    return fig


def fig_dahl(B, Vd, seg):
    fig, ax = plt.subplots(2, 2, figsize=(13, 7), sharex="col")
    for j, (b, vd) in enumerate(zip(B, Vd)):
        for n, (et, r) in enumerate(vd.items(), 1):
            ax[0, j].plot(r["t"], np.rad2deg(r["y"]), lw=1.1, label=et, **pt._stile(n))
            ax[1, j].plot(r["t"], r["i"], lw=0.8, **pt._stile(n))
        ax[0, j].set_title(nome(b), fontsize=9); ax[0, j].legend(fontsize=7); ax[1, j].set_xlabel("time [s]")
    ax[0, 0].set_ylabel("encoder y [deg]"); ax[1, 0].set_ylabel("current [A]")
    pt._titolo(fig, "Check outside the procedure of the book: 1 degree step on the Dahl plant, no current limit", seg)
    return fig


def main() -> None:
    par, Kt, S, T, mot, assi = g.carica()
    seg = par.segnaposto
    pt.applica_stile()
    (CARTELLA / "figures").mkdir(parents=True, exist_ok=True)
    O = json.loads((RADICE / "results_lqg_tuning" / "operating_point.json").read_text())
    (CARTELLA / "weights_to_gains.json").write_text(json.dumps(pesi_guadagni(Kt, S, T, mot, assi, O), indent=1))
    B, C0, G0, R, Vd, out, segni = [], [], [], [], [], {}, []
    for asse, ch, cat, mod in assi:
        b = Banco(asse, ch, cat, mod, Kt, S, T, mot, None if SENZA else O[ch]["maxwell"])
        c0 = b.sintesi()
        primo = b.valuta(c0)
        ordine = sorted(range(len(b.modi)), key=lambda i: -primo["modes"][i]["band_max_dB"])
        im = ordine[0]
        g0 = b.gradino(c0)
        uniformi = []
        for x in UNIFORMI:
            c = b.sintesi(base=x)
            uniformi.append({"q_f": x} | b.valuta(c) | {"step": b.gradino(c)["metriche"]})
        scala = {str(i): [b.valuta(b.sintesi({i: x})) for x in SCALA_MODO] for i in range(len(b.modi))}
        vd = {"q_f = 0.1 on all flexible states": verifica_dahl(b, c0), "q_f = 0": verifica_dahl(b, b.sintesi(base=0.0))}
        for qp, qi in ((Q_P, Q_I), (100.0, 25.0)):
            cq = b.sintesi(q_i=qi, q_p=qp)
            d = dataclasses.replace(cq, k_i=-cq.k_i, k_p=-cq.k_p)
            segni.append({"axis": ch, "q_p": qp, "q_i": qi, "law_8_17": peggiore(poli_chiusi(cq, b.a.Ap, b.a.Bp[:, 0], b.a.Cp[0])), "law_8_15_as_printed": peggiore(poli_chiusi(d, b.a.Ap, b.a.Bp[:, 0], b.a.Cp[0]))})
        anteprima = [{"q_p": qp, "q_i": qi} | b.valuta(b.sintesi(q_i=qi, q_p=qp)) for qp, qi in ANTEPRIMA]
        out[ch] = {"modes": [{k: m[k] for k in ("f_Hz", "zeta", "banda", "q_rif")} for m in b.modi], "order_most_evident_first": ordine, "q_vector": b.vettore(), "first_guess": primo, "first_guess_step": g0["metriche"],
                   "uniform_ladder": uniformi, "mode_ladder_levels": list(SCALA_MODO), "mode_ladder": scala, "preview_higher_gains": anteprima, "V_discretisation": correzione_v(b),
                   "dahl_check": {k: {q: w[q] for q in ("largest_pole", "y_100_120s_mdeg", "y_minus_y_hat_end_mdeg", "step_20s", "i_max_A")} for k, w in vd.items()}}
        print(ch, "fatto", flush=True)
        B.append(b); C0.append(c0); G0.append(g0); Vd.append(vd); R.append({"ordine": ordine, "scala": scala[str(im)], "primo": primo})
    (CARTELLA / "sign_check.json").write_text(json.dumps(segni, indent=1))
    salva = lambda f, n: (pt.salva(f, CARTELLA / "figures", n), plt.close(f))
    salva(fig_primo(B, C0, G0, seg), "first_guess")
    salva(fig_uniformi(B, seg), "process_uniform")
    salva(fig_scala(B, R, seg), "process_mode_ladder")
    salva(fig_traccia(B, R, seg), "process_trace")
    salva(fig_dahl(B, Vd, seg), "dahl_check")
    (CARTELLA / "summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
