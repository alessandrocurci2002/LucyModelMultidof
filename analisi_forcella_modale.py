"""Stima modale ai minimi quadrati (parti reale e immaginaria) dei due assi, telescopio rigido sull'anello.

Richiesta dell'utente (2026-10-31): stessi dati e banda del modello a masse sullo statore (1-150 Hz; polare: le tre
misure nominali; azimut: ``BodeMeasurement_Azimuth_14vit_600_20261028.csv``), con
:mod:`gimbal_mdof.forcella_modale`:

* formulazione modale: modi non smorzati, smorzamento viscoso modale zeta_k <= 0.05 per ogni modo, cuscinetti lineari
  isteretici k_g (1 + j eta_g) proiettati sulle coordinate modali;
* inerzie e rigidezze del CAD entro +-20 %; K_t entro la tolleranza di catalogo +-10 % (fattore g_Kt per asse);
  rigidezze a terra, base - sommita' della forcella e parametri delle due masse sullo statore liberi; telescopio
  rigido sull'anello;
* minimi quadrati sulle differenze delle parti reali e immaginarie (nessun logaritmo, perdita quadratica).

Passi precedenti: zeta_k <= 0.01 (``results_fork_modes_modal_zeta001``), poi zeta_k <= 0.05 con CAD entro +-15 %
(``results_fork_modes_modal_zeta005_J15``), poi inerzie entro +-20 % e rigidezze entro +-15 %
(``results_fork_modes_modal_zeta005_J20``), poi CAD entro +-20 % (``results_fork_modes_modal_cad20``); la prima
partenza e' l'ultimo passo (``--partenza``).

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_forcella_modale.py
    .venv/bin/python -I analisi_forcella_modale.py --cad 0.5 --cad-J 0.5 --kt 0 --cartella results_fork_modes_modal_diag_cad50
"""

from __future__ import annotations

import argparse
import csv
import itertools
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
from gimbal_mdof import forcella_modale as fm  # noqa: E402
from gimbal_mdof import linearize as lz  # noqa: E402
from gimbal_mdof import modi_forcella as mf  # noqa: E402
from gimbal_mdof import params as pr  # noqa: E402
from gimbal_mdof import plotting as pt  # noqa: E402
from gimbal_mdof import sperimentale as spm  # noqa: E402

CARTELLA = "results_fork_modes_modal"
PRECEDENTE = RADICE / "results_fork_modes_2masses" / "summary.json"
"""Stima precedente (errore logaritmico, 2 masse): partenze e confronto."""
PASSO_PRECEDENTE = RADICE / "results_fork_modes_modal_cad20" / "summary.json"
"""Stima modale con zeta_k <= 0.05 e CAD entro +-20 %, K_t nominale: prima partenza."""
ZETA_0 = 0.005
MODULI_RESIDUO = (0.3, 1.5)
ZOOM = (40.0, 150.0)


def argomenti() -> argparse.Namespace:
    a = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    a.add_argument("--config", type=Path, default=RADICE / "config_hub_telescopeonring_alhub_fit_ground_aztel.yaml")
    a.add_argument("--f-min", type=float, default=1.0)
    a.add_argument("--f-max", type=float, default=150.0)
    a.add_argument("--zeta-max", type=float, default=fm.ZETA_MAX, help="limite dello smorzamento modale")
    a.add_argument("--cad", type=float, default=fm.VARIAZIONE_CAD, help="variazione relativa massima delle rigidezze CAD")
    a.add_argument("--cad-J", type=float, default=fm.VARIAZIONE_CAD_J, help="variazione relativa massima delle inerzie CAD")
    a.add_argument("--kt", type=float, default=fm.VARIAZIONE_KT, help="tolleranza di K_t (0: K_t nominale)")
    a.add_argument("--cartella", default=CARTELLA)
    a.add_argument("--partenza", type=Path, default=PASSO_PRECEDENTE,
                   help="summary.json di una stima modale da usare come prima partenza (se esiste)")
    return a.parse_args()


def partenze_polari(prec: dict) -> list[dict[str, float]]:
    """La stima precedente e tutte le combinazioni di segno dei residui con due moduli."""
    base = {p.nome: 1.0 for p in fm.parametri_polari() if p.nome.startswith(("s_", "g_"))}
    base |= {k: prec[k] for k in ("k_g", "eta_g", "tau_d", "f_1", "f_2")}
    base |= {f"zeta_{k}": ZETA_0 for k in range(1, len(fm.NODI_POLARE) + 1)}
    out = [base | {"a_1": prec["a_1"], "a_2": prec["a_2"]}]
    out += [base | {"a_1": s1 * a, "a_2": s2 * a} for s1, s2 in itertools.product((1.0, -1.0), repeat=2)
            for a in MODULI_RESIDUO]
    return out


def partenze_azimut(prec: dict, cat: mf.CatenaAzimut) -> list[dict[str, float]]:
    """La stima precedente e tre ripartizioni del load per due frequenze del modo base - sommita'."""
    base = {p.nome: 1.0 for p in fm.parametri_azimut() if p.nome.startswith(("s_", "g_"))}
    base |= {k: prec[k] for k in ("k_g", "eta_g", "tau_d")}
    base |= {f"zeta_{k}": ZETA_0 for k in range(1, len(fm.NODI_AZIMUT) + 1)}
    out = [base | {"alpha": prec["alpha"], "k_LA": prec["k_LA"]}]
    for alfa, f0 in itertools.product((0.2, 0.4, 0.7), (49.0, 56.0)):
        JL, JA = (1 - alfa) * cat.J_L, alfa * cat.J_L + cat.J_p
        out.append(base | {"alpha": alfa, "k_LA": JL * JA / (JL + JA) * (2 * np.pi * f0) ** 2})
    return out


# ---------------------------------------------------------------------------
# Figure (banda del fit)
# ---------------------------------------------------------------------------


def figura_bode(f, misure, G, nome, banda, lim, lim_zoom, titolo, seg):
    fig, assi = plt.subplots(2, 2, figsize=(12.5, 6.8), sharex="col", gridspec_kw={"width_ratios": [1.6, 1.0]})
    for col, xl in ((0, banda), (1, ZOOM)):
        ax_m, ax_f = assi[0, col], assi[1, col]
        for k, m in enumerate(misure):
            b = (m.f >= banda[0]) & (m.f <= banda[1])
            st = {"color": afm.COLORI_MISURE[k % len(afm.COLORI_MISURE)], "marker": "o", "ms": 1.6, "lw": 0.7,
                  "alpha": 0.85}
            ax_m.plot(m.f[b], m.modulo_db[b], label=f"measured: {m.etichetta}", **st)
            ax_f.plot(m.f[b], m.fase_deg[b], **st)
        st = pt._stile(1)
        ax_m.plot(f, pt._db(G), lw=1.8, label=nome, **st)
        ax_f.plot(f, pt._fase_avvolta(G), lw=1.4, **st)
        for ax in (ax_m, ax_f):
            ax.set_xscale("log")
            ax.set_xlim(*xl)
        ax_m.set_ylim(*(lim if col == 0 else lim_zoom))
        ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
        ax_f.set_ylabel("phase [deg]")
        ax_f.set_ylim(-180, 180)
        ax_f.set_yticks(np.arange(-180, 181, 90))
        ax_f.set_xlabel("frequency [Hz]")
    assi[0, 0].legend(loc="lower left", fontsize=7)
    assi[0, 1].set_title(f"zoom {ZOOM[0]:g}-{ZOOM[1]:g} Hz", fontsize=9)
    pt._titolo(fig, titolo, seg)
    return fig


def figura_contributi(f, m, mod: fm.ModelloModale, banda, lim, titolo, seg):
    """Termini modali della pianta (con i fattori WorkBench) e correzione isteretica dei cuscinetti."""
    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(10.5, 6.8), sharex=True)
    b = (m.f >= banda[0]) & (m.f <= banda[1])
    ax_m.plot(m.f[b], m.modulo_db[b], color=pt.INCHIOSTRO_TENUE, marker="o", ms=1.4, lw=0.6,
              label=f"measured: {m.etichetta}")
    ax_f.plot(m.f[b], m.fase_deg[b], color=pt.INCHIOSTRO_TENUE, marker="o", ms=1.4, lw=0.6)
    G = mod.frf(f)
    ax_m.plot(f, pt._db(G), color=pt.INCHIOSTRO, lw=2.0, label="model: sum")
    ax_f.plot(f, pt._fase_avvolta(G), color=pt.INCHIOSTRO, lw=1.6)
    for k, (nome, Gk) in enumerate(mod.contributi(f).items()):
        st = {"color": pt.INCHIOSTRO_SECONDARIO, "linestyle": ":"} if nome.startswith("bearing") else pt._stile(k)
        ax_m.plot(f, pt._db(Gk), lw=1.2, label=nome, **st)
        ax_f.plot(f, pt._fase_avvolta(Gk), lw=1.0, **st)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(*banda)
    ax_m.set_ylim(*lim)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7, ncol=2)
    pt._titolo(fig, titolo, seg)
    return fig


def figura_errori(curve, banda, seg):
    """Differenza complessa modello - misura nel piano, normalizzata sul modulo della misura, e in dB/gradi."""
    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(10.5, 6.2), sharex=True)
    for k, (nome, f, G, Gm) in enumerate(curve):
        e = np.log(G / Gm)
        st = pt._stile(k)
        ax_m.plot(f, 20 / np.log(10) * e.real, lw=1.0, marker="o", ms=1.4, label=nome, **st)
        ax_f.plot(f, np.rad2deg(e.imag), lw=1.0, marker="o", ms=1.4, **st)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(*banda)
        ax.axhline(0.0, color=pt.INCHIOSTRO_SECONDARIO, lw=0.8)
    ax_m.set_ylim(-15, 15)
    ax_f.set_ylim(-120, 120)
    ax_m.set_ylabel("model - measured [dB]")
    ax_f.set_ylabel("model - measured [deg]")
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower left", fontsize=7)
    pt._titolo(fig, f"Error of the modal models on the measurements, fit band {banda[0]:g}-{banda[1]:g} Hz", seg)
    return fig


# ---------------------------------------------------------------------------


def righe_parametri(e: fm.Esito, iniz: dict[str, float]) -> list[list]:
    out = [["parameter", "unit", "estimate", "start", "lower", "upper", "std_percent", "at_bound", "note"]]
    for p in e.parametri:
        out.append([p.nome, p.unita, f"{e.valori[p.nome]:.5g}", f"{iniz.get(p.nome, float('nan')):.5g}",
                    f"{p.lo:.4g}", f"{p.hi:.4g}", f"{100 * e.std_rel[p.nome]:.1f}",
                    "yes" if e.al_limite(p.nome) else "", p.nota])
    return out


def righe_modi(mod: fm.ModelloModale) -> list[list]:
    out = [["mode", "f_Hz", "zeta", "residue_1_kgm2", "eta_bearing"] + [f"shape_{n}" for n in mod.nodi]]
    out += [[k + 1, f"{m.f:.3f}", f"{m.zeta:.5f}", f"{m.residuo:+.5g}", f"{m.eta_cuscinetti:.4g}"]
            + [f"{m.forma[n]:+.3f}" for n in mod.nodi] for k, m in enumerate(mod.modi())]
    return out


def scrivi_csv(percorso: Path, righe: list[list]) -> None:
    with open(percorso, "w", newline="") as fh:
        csv.writer(fh).writerows(righe)


def modi_json(mod: fm.ModelloModale) -> list[dict]:
    return [{"f_Hz": m.f, "zeta": m.zeta, "residue": m.residuo, "eta_bearing": m.eta_cuscinetti, "shape": m.forma}
            for m in mod.modi()]


def main() -> None:
    args = argomenti()
    pt.applica_stile()
    par = pr.parametri_da_config(args.config)
    cfg = par.sperimentale
    Kt, seg = cfg.Kt, par.segnaposto
    banda = (args.f_min, args.f_max)
    out = RADICE / args.cartella
    p_pol, p_az = (fm.parametri_polari(args.zeta_max, args.cad, args.cad_J, args.kt),
                    fm.parametri_azimut(args.zeta_max, args.cad, args.cad_J, args.kt))
    for d in ("figures", "tables"):
        (out / d).mkdir(parents=True, exist_ok=True)
    cat_az, cat_pol = mf.catene_h7(lz.linearizza(afm.sistema(par).accoppiato, 0.0, 0.0))
    prec = json.loads(PRECEDENTE.read_text())
    pol = [afm.misura(cfg, fn, "polare", et) for fn, et in afm.MISURE_POLARI]
    az = afm.misura(cfg, afm.MISURA_AZIMUT[0], "azimut", afm.MISURA_AZIMUT[1])
    dati_pol = [mf.dato(m.etichetta, m.f, m.G, banda) for m in pol]
    dati_az = [mf.dato(az.etichetta, az.f, az.G, banda)]

    # --- stime ------------------------------------------------------------------------------------------------
    passo = json.loads(args.partenza.read_text()) if args.partenza and args.partenza.exists() else None
    pp = partenze_polari(prec["elevation"]["parameters"])
    if passo:
        pp.insert(0, {"g_Kt": 1.0} | passo["elevation"]["parameters"])
    ep = fm.stima(lambda v: fm.polare(cat_pol, v, Kt), p_pol, dati_pol, pp)
    print(f"polare: costo {ep.costo:.6g} ({len(pp)} partenze)", flush=True)
    pa = partenze_azimut(prec["azimuth"]["parameters"], cat_az)
    if passo:
        pa.insert(0, {"g_Kt": 1.0} | passo["azimuth"]["parameters"])
    ea = fm.stima(lambda v: fm.azimut(cat_az, v, Kt), p_az, dati_az, pa)
    print(f"azimut: costo {ea.costo:.6g} ({len(pa)} partenze)", flush=True)
    mod_p, mod_a = fm.polare(cat_pol, ep.valori, Kt), fm.azimut(cat_az, ea.valori, Kt)

    # --- figure -----------------------------------------------------------------------------------------------
    f = np.geomspace(*banda, 4000)
    fig = out / "figures"
    pt.salva(figura_bode(f, pol, mod_p.frf(f), "modal model: H7 + 2 stator masses, rigid telescope", banda,
                         (-75, -30), (-70, -35), "Elevation plant, modal least-squares fit (real and imaginary parts), "
                         f"{banda[0]:g}-{banda[1]:g} Hz", seg), fig, "polar_bode")
    pt.salva(figura_bode(f, [az], mod_a.frf(f), "modal model: flexible load, rigid telescope", banda, (-95, -20),
                         (-70, -25), "Azimuth plant, modal least-squares fit (real and imaginary parts), "
                         f"{banda[0]:g}-{banda[1]:g} Hz", seg), fig, "azimuth_bode")
    pt.salva(figura_contributi(f, pol[0], mod_p, banda, (-110, -30),
                               "Elevation: modal terms R_k/(s^2 + 2 zeta_k omega_k s + omega_k^2) and bearing "
                               "hysteresis term", seg), fig, "polar_modal_terms")
    pt.salva(figura_contributi(f, az, mod_a, banda, (-120, -10),
                               "Azimuth: modal terms R_k/(s^2 + 2 zeta_k omega_k s + omega_k^2) and bearing "
                               "hysteresis term", seg), fig, "azimuth_modal_terms")
    curve = [(d.etichetta, d.f, mod_p.frf(d.f), d.G) for d in dati_pol] + [
        (d.etichetta, d.f, mod_a.frf(d.f), d.G) for d in dati_az]
    pt.salva(figura_errori(curve, banda, seg), fig, "fit_errors")
    plt.close("all")

    # --- tabelle e riepilogo ----------------------------------------------------------------------------------
    scrivi_csv(out / "tables" / "polar_parameters.csv", righe_parametri(ep, pp[0]))
    scrivi_csv(out / "tables" / "azimuth_parameters.csv", righe_parametri(ea, pa[0]))
    scrivi_csv(out / "tables" / "polar_modes.csv", righe_modi(mod_p))
    scrivi_csv(out / "tables" / "azimuth_modes.csv", righe_modi(mod_a))
    q = [afm.qualita("modal LS: H7 + 2 stator masses", d, mod_p.frf(d.f)) for d in dati_pol]
    q += [afm.qualita("modal LS: flexible load", d, mod_a.frf(d.f)) for d in dati_az]
    scrivi_csv(out / "tables" / "fit_quality.csv", [list(q[0])] + [list(r.values()) for r in q])
    scrivi_csv(out / "tables" / "starts.csv", [["axis", "start", "cost"]] + [
        [asse, json.dumps({k: round(v, 6) for k, v in p0.items()}), f"{c:.6g}"]
        for asse, e in (("elevation", ep), ("azimuth", ea)) for p0, c in e.partenze])
    sommario = {
        "config": args.config.name, "fit_band_Hz": list(banda), "Kt_Nm_A": Kt, "current_loop_Hz": mf.F_CORRENTE,
        "zeta_max": args.zeta_max, "cad_variation": args.cad, "cad_inertia_variation": args.cad_J,
        "kt_variation": args.kt,
        "first_start": str(args.partenza.relative_to(RADICE)) if passo else None,
        "residual": "w_i [Re(G - G_mis), Im(G - G_mis)], w_i = sqrt(dln f / 0.01), quadratic loss",
        "elevation": {"measurements": [fn for fn, _ in afm.MISURE_POLARI], "cost": ep.costo, "nfev": ep.nfev,
                      "parameters": ep.valori, "std_rel": ep.std_rel,
                      "at_bound": [p.nome for p in ep.parametri if ep.al_limite(p.nome)],
                      "starts": [c for _, c in ep.partenze], "modes": modi_json(mod_p),
                      "M_diag": np.diag(mod_p.M).tolist(), "nodes": list(mod_p.nodi),
                      "stator_masses": [{"f_Hz": ep.valori[f"f_{r}"], "a": ep.valori[f"a_{r}"],
                                         "J_kgm2": 1 / abs(ep.valori[f"a_{r}"]),
                                         "kappa": float(np.sign(ep.valori[f"a_{r}"])),
                                         "k_Nm_rad": (2 * np.pi * ep.valori[f"f_{r}"]) ** 2 / abs(ep.valori[f"a_{r}"])}
                                        for r in (1, 2)]},
        "azimuth": {"measurement": afm.MISURA_AZIMUT[0], "cost": ea.costo, "nfev": ea.nfev,
                    "parameters": ea.valori, "std_rel": ea.std_rel,
                    "at_bound": [p.nome for p in ea.parametri if ea.al_limite(p.nome)],
                    "starts": [c for _, c in ea.partenze], "modes": modi_json(mod_a),
                    "M_diag": np.diag(mod_a.M).tolist(), "nodes": list(mod_a.nodi)},
        "fit_quality": q,
        "cad": {"elevation_J": np.diag(cat_pol.M).tolist(), "elevation_roles": list(cat_pol.ruoli),
                "elevation_springs": [[i, j, k] for i, j, k in fm.molle(cat_pol.K)[0]],
                "azimuth": {"J_e": cat_az.J_e, "J_m": cat_az.J_m, "J_L": cat_az.J_L, "J_p": cat_az.J_p,
                            "k": list(cat_az.k)}},
    }
    (out / "summary.json").write_text(json.dumps(sommario, indent=1))
    print(json.dumps({k: {kk: sommario[k][kk] for kk in ("cost", "parameters", "at_bound")}
                      for k in ("elevation", "azimuth")}, indent=1))
    for r in q:
        print(r)


if __name__ == "__main__":
    main()
