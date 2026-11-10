"""Figure della tesi (cap. "Pointer Model System Identification"): modello della forcella con cuscinetti di Dahl.

Richiesta dell'utente (2026-11-05): figure del capitolo di identificazione nello stile della tesi (``style.py``), con
l'asse delle frequenze senza numeri (policy aziendale). Nessuna stima qui: parametri della stima su griglia uniforme in
ln f (``results_fork_modes_dahl_resampled/summary.json``, variante ``uniform_lnf``) e guadagni efficaci salvati in
``time_domain_summary.json``. Simulatori e dati sono quelli di ``analisi_fit_bode_ricampionato.py`` e
``analisi_tempo_ricampionato.py``.

* ``ident_frf_pol``, ``ident_frf_az``: risposte misurate e modello al canale dell'encoder;
* ``ident_cl``: mossa da 1 grado dei doppietti in anello chiuso (riferimento misurato in ingresso), entrambi gli assi.

Il fit in anello aperto (corrente misurata in ingresso) e' stampato a terminale per le stesse mosse.

    .venv/bin/python -I analisi_tesi_identificazione.py

Figure in Thesis-Writing/plots/identification/.
"""

from __future__ import annotations

import dataclasses
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

RADICE = Path(__file__).resolve().parent
sys.path[:0] = [str(RADICE / "src"), str(RADICE), str(RADICE / "closedLoopModel")]

import analisi_anello_chiuso_ibrido as aac  # noqa: E402
import analisi_forcella_lugre as afl  # noqa: E402
import analisi_modi_forcella as afm  # noqa: E402
from gimbal_mdof import doppietti_tempo as dt  # noqa: E402
from gimbal_mdof import forcella_lugre as fl  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from lucy_cl.driver import Drive, DriveParams  # noqa: E402
from lucy_cl.motor import MotorParams  # noqa: E402

STILE = Path("/home/alessandro/Desktop/Thesis/OpenVins/OpenVincenzo/colcon_ws/src/open_vins/ov_data/sc_cantina01_thesis/output/plots/scripts")
FIGURA = Path("/home/alessandro/Desktop/Thesis/Thesis-Writing/plots/identification")
CARTELLA = RADICE / "results_fork_modes_dahl_resampled"
VAR = "uniform_lnf"
MODELLO = dict(color="#0072B2", ls="-", lw=1.3)
MISURA = dict(color="#000000", ls="-", lw=1.1)
COMANDO = dict(color="#000000", ls=":", lw=0.8)
GRIGI = ("#000000", "#777777", "#BBBBBB")
ETICHETTE_POL = ("Measurement, full band", "Measurement, resonance band", "Measurement, resonance band, re-assembly")
MOSSE = (("01deg_minus01deg", "#1", r"$0.1^\circ$"), ("1deg_minus1deg", "#1", r"$1^\circ$"),
         ("5deg_minus5deg", "#2", r"$5^\circ$"))


def asse_frequenze(ax, f0, f1):
    for x in ax:
        x.set_xlim(f0, f1)
        # frequenze non divulgabili: griglia logaritmica senza numeri
        x.tick_params(axis="x", which="both", labelbottom=False)
        x.grid(True, which="minor", axis="x", color="#EEEEEE", lw=0.4)
    ax[0].set_ylabel("Magnitude [dB (rad/s)/A]")
    ax[1].set_ylabel("Phase [deg]")
    ax[1].set_ylim(-180, 180)
    ax[1].set_yticks(np.arange(-180, 181, 90))
    ax[1].set_xlabel("Frequency (logarithmic scale)")


def legenda(fig, ax, ncol, alto=0.93):
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=ncol, frameon=True, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, alto))


def figura_frf(style, dati, modello, etichette, nome, ms=2.0):
    fig, ax = style.figure(aspect=0.5, nrows=2, ncols=1, sharex=True)
    for d, col, et in zip(dati, GRIGI, etichette):
        ax[0].semilogx(d.f, 20 * np.log10(np.abs(d.G)), color=col, marker="o", ms=ms, lw=0.6, label=et)
        ax[1].semilogx(d.f, np.rad2deg(np.angle(d.G)), color=col, marker="o", ms=ms, lw=0.6)
    for k, d in enumerate(dati):
        G = modello.frf_dato(d)
        ax[0].semilogx(d.f, 20 * np.log10(np.abs(G)), **MODELLO, label="Model" if k == 0 else None)
        ax[1].semilogx(d.f, np.rad2deg(np.angle(G)), **MODELLO)
    asse_frequenze(ax, min(d.f[0] for d in dati), max(d.f[-1] for d in dati))
    legenda(fig, ax[0], ncol=2, alto=0.88 if len(dati) > 1 else 0.92)
    style.save(fig, FIGURA / nome)


def finestra(finestre, prova, mossa):
    return next(w for w in finestre[prova] if w.nome.endswith(mossa))


def inizio(w):
    vc = np.abs(np.gradient(w.cmd, dt.DT)) > 1e-6
    return int(np.flatnonzero(vc)[0])


def main() -> None:
    sys.path.insert(0, str(STILE))
    import style  # noqa: E402
    style.setup()
    FIGURA.mkdir(parents=True, exist_ok=True)
    par = pr.parametri_da_config(RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    cfg, Kt = par.sperimentale, par.sperimentale.Kt
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    s = json.loads((CARTELLA / "summary.json").read_text())
    g_salvati = json.loads((CARTELLA / "time_domain_summary.json").read_text())["gains"]
    mot = MotorParams.from_yaml(RADICE / "closedLoopModel" / "data" / "motor_TBM7646A.yaml", "cold")

    # risposte in frequenza
    pol = [afm.misura(cfg, fn, "polare", et) for fn, et in afm.MISURE_POLARI]
    az = [afm.misura(cfg, afm.MISURA_AZIMUT[0], "azimut", afm.MISURA_AZIMUT[1])]
    dp, ap = afl.dati_con_ampiezza(cfg, pol)
    da, aa = afl.dati_con_ampiezza(cfg, az)
    figura_frf(style, dp, fl.polare(cat_pol, s["elevation"][VAR]["hybrid"]["parameters"], Kt, ap), ETICHETTE_POL,
               "ident_frf_pol", ms=0.9)
    figura_frf(style, da, fl.azimut(cat_az, s["azimuth"][VAR]["hybrid"]["parameters"], Kt, aa), ("Measurement",),
               "ident_frf_az")

    # doppietti
    curve = {}
    for asse, chiave, cat, mod in (("azimut", "azimuth", cat_az, fl.azimut), ("polare", "elevation", cat_pol, fl.polare)):
        v = s[chiave][VAR]["hybrid"]["parameters"]
        m = mod(cat, v, Kt, {}).struttura
        att = dt.Attrito(v["sigma0"], v["Tc"])
        akd, wb = aac.DRIVE[asse]
        p_file = DriveParams.from_files(chiave, aac.VCL / akd, aac.VCL / wb, pl_kp_unit="hz")
        an = aac.AnelloChiuso(m, Drive(dataclasses.replace(p_file, **g_salvati[chiave][VAR])), mot, att)
        sim = dt.Simulatore(m)
        finestre = {p: dt.finestre(aac.VCL / "V2" / (aac.PROVE[asse] % p), asse) for p, _, _ in MOSSE}
        for prova, mossa, _ in MOSSE:
            w = finestra(finestre, prova, mossa)
            y_ol = sim.corri(w.I, att)
            t = np.arange(0.0, w.t[-1], aac.PASSO)
            r = an.corri(np.interp(t, w.t, w.cmd - w.cmd[0]), float(w.I[:int(0.4 / dt.DT)].mean()))
            curve[(asse, prova)] = (w, y_ol, np.interp(w.t, t, r["p"]), np.interp(w.t, t, r["i"]))
            print(f"{asse} {prova} {mossa}: fit anello aperto {dt.fit(y_ol, w.pos):.1f} %, "
                  f"anello chiuso {dt.fit(curve[(asse, prova)][2], w.pos):.1f} %", flush=True)

    md = lambda q: np.rad2deg(q) * 1e3  # noqa: E731
    # anello chiuso, mossa da 1 grado: posizione, errore di inseguimento, corrente
    fig, ax = style.figure(aspect=0.78, nrows=3, ncols=2, sharex="col")
    for j, (asse, nome) in enumerate((("azimut", "Azimuth"), ("polare", "Polar axis"))):
        w, _, y, im = curve[(asse, MOSSE[1][0])]
        i0 = inizio(w)
        k = (w.t >= w.t[i0] - 0.2) & (w.t <= w.t[i0] + 2.0)
        t, c = w.t[k], w.cmd[k] - w.cmd[0]
        ax[0, j].plot(t, md(c), **COMANDO, label="Command")
        ax[0, j].plot(t, md(w.pos[k]), **MISURA, label="Measurement")
        ax[0, j].plot(t, md(y[k]), **MODELLO, label="Model")
        ax[1, j].plot(t, md(w.pos[k] - c), **MISURA)
        ax[1, j].plot(t, md(y[k] - c), **MODELLO)
        ax[2, j].plot(t, w.I[k], **MISURA)
        ax[2, j].plot(t, im[k], **MODELLO)
        ax[0, j].set_title(f"{nome}, {MOSSE[1][2]}", fontsize=9)
        ax[2, j].set_xlabel("Time [s]")
    ax[0, 0].set_ylabel("Position [mdeg]")
    ax[1, 0].set_ylabel("Error [mdeg]")
    ax[2, 0].set_ylabel(r"Current [A$_\mathrm{rms}$]")
    legenda(fig, ax[0, 0], ncol=3, alto=0.94)
    style.save(fig, FIGURA / "ident_cl")

if __name__ == "__main__":
    main()
