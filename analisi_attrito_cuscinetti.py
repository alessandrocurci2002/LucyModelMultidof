"""Cuscinetti polari con attrito di Dahl/LuGre nel modello a 2 masse sullo statore: compare una risonanza a 6-9 Hz?

Richiesta dell'utente (2026-10-31): nel modello a masse sullo statore (2 masse, ``results_fork_modes_2masses``) il
vincolo lineare isteretico k_g(1+j eta_g) dei cuscinetti polari, collegati alla struttura dello statore, e' sostituito
dalla rigidezza equivalente (funzione descrittiva, prima armonica) dell'attrito di Dahl identificato sui doppietti V2
(``closedLoopModel/results_v2_identification``), frequenza per frequenza all'ampiezza di moto del cuscinetto. Nessuna
massa aggiunta. Ampiezza del cuscinetto:

* ``misurata``: |G_mis| I |S_mis| / omega, con I la corrente iniettata e S_mis = 1/(1+L_mis) la sensibilita' misurata
  (la prova e' in anello chiuso: la corrente nel motore e' I S); il moto dello statore e' trascurabile sotto 20 Hz;
* ``coerente``: bilancio armonico in anello chiuso, l'ampiezza che il modello stesso da' con la rigidezza equivalente a
  quell'ampiezza e con l'anello chiuso dal modello (L = C B G, controllore e biquad di retroazione del file).

Varianti dell'attrito: Dahl (sigma_1 = 0) e LuGre con sigma_1 da zeta = 0.15 sulla risonanza di pre-sliding. Gli altri
parametri restano quelli del fit a 2 masse (nessuna nuova stima).

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_attrito_cuscinetti.py
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))

import matplotlib.pyplot as plt  # noqa: E402

import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402

CARTELLA = RADICE / "results_fork_modes_2masses" / "bearing_friction"
DAHL = RADICE / "closedLoopModel" / "results_v2_identification" / "summary.json"
BANDA = (1.0, 30.0)
ZETA_LUGRE = 0.15
"""Smorzamento della risonanza di pre-sliding che fissa sigma_1 (ordine di grandezza di Piasek et al. 2019)."""
GRIGLIA_A = np.geomspace(5e-6, 2e-2, 48)
"""Ampiezze della tabella delle funzioni descrittive, rad."""


def funzione_descrittiva(p: dict, A: float, f: float, sigma1: float, n: int = 1500, cicli: int = 3) -> complex:
    """Prima armonica coppia/spostamento dell'attrito sotto x = A sin(omega t), N m/rad (regime: ultimo ciclo).

    Stato z (deformazione media delle setole): dz/dt = v - sigma0 |v| z / g(v), g(v) = T_c + (T_s - T_c) exp(-(v/v_s)^2);
    coppia sigma0 z + sigma1 dz/dt. Integrazione esatta a velocita' costante sul passo."""
    w = 2 * np.pi * f
    t = np.arange(cicli * n) / (n * f)
    dt = t[1]
    x, v = A * np.sin(w * t), A * w * np.cos(w * t)
    s0, Tc, Ts, vs = p["sigma0"], p["Tc"], p["Tc"] * p["Ts_ratio"], p["vs"]
    z, coppia = 0.0, np.empty_like(t)
    for i, vi in enumerate(v):
        g = Tc + (Ts - Tc) * np.exp(-(vi / vs) ** 2)
        a = s0 * abs(vi) / g
        e = np.exp(-a * dt)
        z = z * e + vi / a * (1 - e) if a > 0 else z + vi * dt
        coppia[i] = s0 * z + sigma1 * (vi - a * z)
    k = t >= (cicli - 1) / f
    ph = np.exp(-1j * w * t[k])
    return complex(np.mean(coppia[k] * ph) / np.mean(x[k] * ph))


def tabella(p: dict, f: np.ndarray, sigma1: float) -> np.ndarray:
    return np.array([[funzione_descrittiva(p, A, fi, sigma1) for A in GRIGLIA_A] for fi in f])


def interpola(tab: np.ndarray, i: int, A: float) -> complex:
    la = np.log(np.clip(A, GRIGLIA_A[0], GRIGLIA_A[-1]))
    lg = np.log(GRIGLIA_A)
    return complex(np.interp(la, lg, tab[i].real) + 1j * np.interp(la, lg, tab[i].imag))


class Cuscinetti:
    """Modello a 2 masse con la rigidezza dei cuscinetti assegnata per frequenza."""

    def __init__(self, mod: mf.ModelloPolare):
        senza = mf.ModelloPolare(**{**mod.__dict__, "k_g": 0.0, "eta_g": 0.0})
        self.M, self.C, self.K, self.e, self.c = senza.matrici()
        n0 = mod.catena.M.shape[0]
        self.v = np.zeros(self.K.shape[0])
        self.v[mod.catena.j("h")] = 1.0
        self.v[n0:] = -1.0
        self.mod = mod

    def risposta(self, f: float, Kg: complex) -> tuple[complex, complex]:
        """(pianta WorkBench, rotazione relativa del cuscinetto per unita' di corrente) a ``f``."""
        s = 2j * np.pi * f
        P = self.M * s ** 2 + self.C * s + self.K + Kg * np.outer(self.v, self.v)
        x = np.linalg.solve(P, self.e.astype(complex))
        m = self.mod
        fatt = m.Kt * (2 * np.pi * m.f_corrente) / (s + 2 * np.pi * m.f_corrente)
        G = s * fatt * (self.c @ x) * np.exp(-s * m.ritardo)
        return G, fatt * (self.v @ x)

    def coerente(self, f: np.ndarray, tab: np.ndarray, I: float, A0: np.ndarray,
                 CB: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Bilancio armonico in anello chiuso: A = |rotazione del cuscinetto| I |1/(1 + C B G(A))| con K_eq(A); punto
        fisso in log A dalla stima A0. ``CB``: controllore per biquad di retroazione alle frequenze ``f``."""
        G, A = np.empty(f.size, complex), np.empty(f.size)
        for i, fi in enumerate(f):
            a = A0[i]
            for _ in range(200):
                g, r = self.risposta(fi, interpola(tab, i, a))
                nuovo = abs(r) * I / abs(1 + CB[i] * g)
                if abs(np.log(nuovo / a)) < 1e-4:
                    break
                a = float(np.exp(0.7 * np.log(a) + 0.3 * np.log(nuovo)))
            A[i] = a
            G[i] = self.risposta(fi, interpola(tab, i, a))[0]
        return G, A


def main() -> None:
    pt.applica_stile()
    s2 = json.loads((RADICE / "results_fork_modes_2masses" / "summary.json").read_text())
    par = pr.parametri_da_config(RADICE / s2["config"])
    cfg = par.sperimentale
    Kt, seg = cfg.Kt, par.segnaposto
    _, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    mod = mf.modello_polare(cat_pol, s2["elevation"]["parameters"], len(s2["stator_mass_starts_Hz"]), Kt)
    dahl = json.loads(DAHL.read_text())["final"]["elevation"]
    fn, et = afm.MISURE_POLARI[0]
    m = afm.misura(cfg, fn, "polare", et)
    I = m.eccitazione_A
    sel = (m.f >= BANDA[0]) & (m.f <= BANDA[1])
    f, Gm = m.f[sel], m.G[sel]
    _, tracce, _ = spm.carica_bode(cfg.cartella / fn)
    traccia = lambda nome: (10 ** (tracce[nome][0] / 20) * np.exp(1j * np.deg2rad(tracce[nome][1])))[sel]  # noqa: E731
    CB = traccia("Controller") * traccia("Feedback Biquad")
    S_mis = 1 / (1 + traccia("Open Loop"))
    A_mis = np.abs(Gm) * I * np.abs(S_mis) / (2 * np.pi * f)
    J = cat_pol.J
    sigma1 = 2 * ZETA_LUGRE * np.sqrt(dahl["sigma0"] * J)
    print(f"Dahl elevazione: sigma0 {dahl['sigma0']:.0f} N m/rad, Tc {dahl['Tc']:.3f} N m, pre-sliding "
          f"{1e6 * dahl['Tc'] / dahl['sigma0']:.0f} urad; sqrt(sigma0/J)/2pi = "
          f"{np.sqrt(dahl['sigma0'] / J) / (2 * np.pi):.2f} Hz; sigma1 LuGre {sigma1:.1f} N m s/rad", flush=True)
    tab_d = tabella(dahl, f, 0.0)
    tab_l = tabella(dahl, f, sigma1)
    print("tabelle delle funzioni descrittive pronte", flush=True)

    cu = Cuscinetti(mod)
    curve = {"linear hysteretic k_g(1+j eta) (fitted)": mod.frf(f)}
    curve["Dahl, measured amplitude"] = np.array([cu.risposta(fi, interpola(tab_d, i, A_mis[i]))[0]
                                                  for i, fi in enumerate(f)])
    G_dc, A_dc = cu.coerente(f, tab_d, I, A_mis, CB)
    curve["Dahl, self-consistent amplitude"] = G_dc
    G_lc, A_lc = cu.coerente(f, tab_l, I, A_mis, CB)
    curve[f"LuGre (sigma1 = {sigma1:.0f}), self-consistent amplitude"] = G_lc

    # --- figura -------------------------------------------------------------------------------------------------
    CARTELLA.mkdir(parents=True, exist_ok=True)
    fig, (ax_m, ax_f, ax_a) = plt.subplots(3, 1, figsize=(10.5, 9.2), sharex=True,
                                           gridspec_kw={"height_ratios": [1.3, 1.0, 0.8]})
    ax_m.plot(f, pt._db(Gm), color=pt.INCHIOSTRO, marker="o", ms=2.0, lw=0.8, label=f"measured: {et}")
    ax_f.plot(f, np.rad2deg(np.angle(Gm)), color=pt.INCHIOSTRO, marker="o", ms=2.0, lw=0.8)
    for k, (nome, G) in enumerate(curve.items()):
        st = pt._stile(k)
        ax_m.plot(f, pt._db(G), lw=1.6, label=nome, **st)
        ax_f.plot(f, np.rad2deg(np.angle(G)), lw=1.3, **st)
    p0 = dahl["Tc"] / dahl["sigma0"]
    for A, nome, k in ((A_mis, "measured (|G| I |S| / omega)", None), (A_dc, "Dahl, self-consistent", 2),
                       (A_lc, "LuGre, self-consistent", 3)):
        st = {"color": pt.INCHIOSTRO, "marker": "o", "ms": 2.0, "lw": 0.8} if k is None else pt._stile(k)
        ax_a.plot(f, 1e6 * A, label=nome, **st)
    ax_a.axhline(1e6 * p0, color=pt.INCHIOSTRO_SECONDARIO, lw=0.9, ls=":")
    ax_a.text(1.05, 1e6 * p0 * 1.15, "pre-sliding scale T_c / sigma_0", fontsize=7, color=pt.INCHIOSTRO_SECONDARIO)
    for ax in (ax_m, ax_f, ax_a):
        ax.set_xscale("log")
        ax.set_xlim(*BANDA)
    ax_a.set_yscale("log")
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_a.set_ylabel("bearing amplitude [urad]")
    ax_a.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7)
    ax_a.legend(loc="lower left", fontsize=7)
    pt._titolo(fig, "Elevation, 2 stator masses: bearings on the stator structure with Dahl/LuGre friction "
                    "(describing function at the bearing amplitude)", seg)
    pt.salva(fig, CARTELLA, "bearing_friction_bode")

    # --- numeri -------------------------------------------------------------------------------------------------
    out = {"dahl": dahl, "sigma1_lugre": sigma1, "J_axis": J, "excitation_A": I,
           "presliding_scale_rad": p0, "f_presliding_Hz": float(np.sqrt(dahl["sigma0"] / J) / (2 * np.pi)),
           "fit_1_30Hz": {}, "features": {}}
    for nome, G in [("measured", Gm)] + list(curve.items()):
        if nome != "measured":
            out["fit_1_30Hz"][nome] = round(spm.fit_nrmse(Gm, G), 1)
        out["features"][nome] = [[tipo, round(float(fx), 2), round(float(g), 1)] for tipo, fx, g in
                                 spm.estremi(f, pt._db(G), prominenza_db=1.0, finestra_dec=0.03)]
    (CARTELLA / "summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("fit_1_30Hz", "features")}, indent=1))


if __name__ == "__main__":
    main()
