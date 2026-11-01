"""Verifica: il giunto anello-telescopio ha la stessa rigidezza in polare (k_b) e in yaw (k di azimut)?

Richiesta dell'utente (2026-10-26), sul modello H7 (``azimuth.telescope_attachment``). Per ogni coppia di
misure prese nella stessa sessione (azimut + polare) si cerca separatamente la rigidezza che meglio
riproduce ciascun asse:

* polare: ``polar.hub.telescope_attachment.k_b_Nm_rad``. Encoder e motore stanno ai due lati del mozzo, quindi
  l'antirisonanza all'encoder e' quella del telescopio sul giunto con il mozzo fermo, sqrt(k_b / I_JJ): non
  dipende dal resto del modello;
* azimut: ``azimuth.telescope_attachment.k_Nm_rad`` (coppia polo-zero del telescopio in yaw).

L'indice e' l'errore RMS del modulo in dB fra modello (canale encoder, con il vincolo a terra) e misura nella
banda :data:`BANDA`. Se i minimi dei due assi cadono sulla stessa k, la rigidezza puo' essere uguale.

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_rigidezza_assi.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))
sys.path.insert(0, str(RADICE))

import matplotlib.pyplot as plt  # noqa: E402

from analisi_nodo_azimut import sistema  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402
from gimbal_mdof import vincolo_terra as vt  # noqa: E402

BANDA = (20.0, 150.0)
"""Banda del confronto, Hz: sopra il vincolo a terra e lo zero a 12 Hz, sotto i modi strutturali alti."""
SESSIONI = (
    ("23 Oct 2026, 14 screws", "VClosedLoop/VTesi/BodeMeasurement_Azimuth_14viti.csv",
     "VClosedLoop/VTesi/BodeMeasurement_Polar_14viti.csv"),
    ("20 Oct 2026", "VClosedLoop/BodeMeasurement_azimuth.csv", "VClosedLoop/BodeMeasurement_elevation.csv"),
)
"""``(etichetta, misura di azimut, misura polare)`` prese nella stessa sessione."""
CHIAVE = {"azimut": ("azimuth", "telescope_attachment", "k_Nm_rad"),
          "polare": ("polar", "hub", "telescope_attachment", "k_b_Nm_rad")}


def argomenti() -> argparse.Namespace:
    a = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    a.add_argument("--config", type=Path, default=RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    a.add_argument("--k", type=float, nargs=3, default=[1.0e4, 3.0e6, 61], metavar=("K_MIN", "K_MAX", "N"),
                   help="griglia logaritmica di k [N m/rad], uguale per i due assi")
    return a.parse_args()


def imposta(dati: dict, asse: str, k: float) -> None:
    nodo = dati
    *strada, ultima = CHIAVE[asse]
    for s in strada:
        nodo = nodo[s]
    nodo[ultima] = float(k)


def frf(percorso: Path, dati: dict, k: float, f: np.ndarray, Kt: float) -> dict[str, np.ndarray]:
    """Pianta velocita'/corrente all'encoder dei due assi con la stessa ``k`` su entrambi i giunti."""
    for asse in CHIAVE:
        imposta(dati, asse, k)
    par = pr.parametri_da_config(percorso, dati=dati)
    ml = vt.applica(lz.linearizza(sistema(par).accoppiato, 0.0, 0.0), par)
    return {asse: spm.frf_velocita(ml, asse, f, Kt)[spm.canali(asse)[0]] for asse in CHIAVE}


def errore_db(m: spm.Misura, G: np.ndarray) -> float:
    ok = (m.f >= BANDA[0]) & (m.f <= BANDA[1])
    return float(np.sqrt(np.mean((m.modulo_db[ok] - 20 * np.log10(np.abs(G[ok]))) ** 2)))


def main() -> None:
    args = argomenti()
    pt.applica_stile()
    par = pr.parametri_da_config(args.config)
    cfg = par.sperimentale
    assert cfg is not None
    voci = {v.file: v for v in cfg.misure}
    misure = {s[0]: {"azimut": spm.carica_misura(cfg, voci[s[1]]), "polare": spm.carica_misura(cfg, voci[s[2]])}
              for s in SESSIONI}
    out = RADICE / par.cartella_risultati
    (out / "tables").mkdir(parents=True, exist_ok=True)

    # le misure di una sessione hanno la stessa griglia di frequenze solo a volte: si valuta il modello
    # sull'unione dei punti e si interpola sul singolo file
    f_tutte = np.unique(np.concatenate([m.f for d in misure.values() for m in d.values()]))
    f_tutte = f_tutte[(f_tutte >= BANDA[0]) & (f_tutte <= BANDA[1])]
    ks = np.geomspace(*args.k[:2], int(args.k[2]))
    dati = pr.carica_config(args.config)
    errori = {s: {a: [] for a in CHIAVE} for s in misure}
    for k in ks:
        G = frf(args.config, dati, k, f_tutte, cfg.Kt)
        for s, d in misure.items():
            for asse, m in d.items():
                ok = (m.f >= BANDA[0]) & (m.f <= BANDA[1])
                Gm = np.full(m.f.shape, np.nan, dtype=complex)
                Gm[ok] = np.interp(m.f[ok], f_tutte, G[asse].real) + 1j * np.interp(m.f[ok], f_tutte, G[asse].imag)
                errori[s][asse].append(errore_db(m, Gm))

    righe, sommario = [], {}
    for s in misure:
        sommario[s] = {}
        for asse in CHIAVE:
            e = np.array(errori[s][asse])
            i = int(np.argmin(e))
            # intervallo di k entro 1 dB dal minimo: quanto la misura vincola k
            entro = ks[e <= e[i] + 1.0]
            sommario[s][asse] = {"k_best": float(ks[i]), "err_db": float(e[i]), "k_within_1dB": [float(entro.min()),
                                 float(entro.max())], "err_db_k_116e3": float(np.interp(np.log(1.16e5), np.log(ks), e))}
        e_tot = np.array(errori[s]["azimut"]) + np.array(errori[s]["polare"])
        i = int(np.argmin(e_tot))
        sommario[s]["k_uguale"] = {"k": float(ks[i]), "err_db_azimut": float(errori[s]["azimut"][i]),
                                   "err_db_polare": float(errori[s]["polare"][i])}
        righe += [[s, f"{k:.5g}", f"{a:.3f}", f"{p:.3f}"] for k, a, p in zip(ks, errori[s]["azimut"], errori[s]["polare"])]
    with open(out / "tables" / "equal_stiffness_scan.csv", "w", newline="") as fh:
        csv.writer(fh).writerows([["session", "k_Nm_rad", "rms_db_azimuth", "rms_db_elevation"]] + righe)
    (out / "equal_stiffness.json").write_text(json.dumps(sommario, indent=2))

    # --- figure: errore contro k e Bode con le k migliori di ciascun asse -------------------------------
    n = len(misure)
    fig, assi = plt.subplots(n, 3, figsize=(15.0, 4.2 * n), gridspec_kw={"width_ratios": [1.0, 1.3, 1.3]})
    assi = np.atleast_2d(assi)
    f = np.geomspace(*BANDA, 800)
    for r, (s, d) in enumerate(misure.items()):
        ax = assi[r, 0]
        for c, (asse, nome) in enumerate((("azimut", "azimuth (yaw k)"), ("polare", "elevation (k_b)"))):
            ax.plot(ks, errori[s][asse], color=pt.SERIE[c], ls=pt.TRATTI[c], lw=1.6, label=nome)
            ax.axvline(sommario[s][asse]["k_best"], color=pt.SERIE[c], lw=0.8, ls=":")
        ax.axvline(1.16e5, color=pt.INCHIOSTRO_TENUE, lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("k [N m/rad] (same on both joints)")
        ax.set_ylabel(f"RMS magnitude error {BANDA[0]:g}-{BANDA[1]:g} Hz [dB]")
        ax.set_title(s, fontsize=9)
        ax.legend(fontsize=7)
        casi = {f"best yaw k = {sommario[s]['azimut']['k_best']:.3g}": sommario[s]["azimut"]["k_best"],
                f"best k_b = {sommario[s]['polare']['k_best']:.3g}": sommario[s]["polare"]["k_best"]}
        G = {nome: frf(args.config, dati, k, f, cfg.Kt) for nome, k in casi.items()}
        for c, asse in enumerate(("azimut", "polare"), start=1):
            ax = assi[r, c]
            m = d[asse]
            ax.plot(m.f, m.modulo_db, color=pt.INCHIOSTRO, marker="o", ms=2.0, lw=0.8, label=f"measured: {m.etichetta}")
            for j, (nome, Gk) in enumerate(G.items()):
                ax.plot(f, pt._db(Gk[asse]), lw=1.5, label=f"model, {nome}", **pt._stile(j))
            ax.set_xscale("log")
            ax.set_xlim(*BANDA)
            ax.set_ylim(-95, -25)
            ax.set_xlabel("frequency [Hz]")
            ax.set_ylabel("magnitude [dB (rad/s)/A]")
            ax.set_title(f"{'azimuth' if asse == 'azimut' else 'elevation'} encoder", fontsize=9)
            ax.legend(fontsize=6.5, loc="lower left")
    pt._titolo(fig, "Same ring-telescope stiffness on both axes? Best k of each axis per test session (theta0 = 0)",
               par.segnaposto)
    figure = pt.salva(fig, out / "figures", "equal_stiffness")
    print(json.dumps(sommario, indent=2))
    print("figure:", *[str(p.relative_to(RADICE)) for p in figure])


if __name__ == "__main__":
    main()
