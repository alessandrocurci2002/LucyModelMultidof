"""Esegue l'intera catena del modello multi-GdL e scrive i risultati in ``results/``.

Esempi::

    .venv/bin/python run_all.py
    .venv/bin/python run_all.py --quick --no-latex
    .venv/bin/python run_all.py --theta0 -15 15 --omega0 0 10 --mass-model consistent
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from gimbal_mdof import analysis as an
from gimbal_mdof import coupled as cp
from gimbal_mdof import fem
from gimbal_mdof import geometry as geo
from gimbal_mdof import latex_report as lr
from gimbal_mdof import linearize as lz
from gimbal_mdof import params as pr
from gimbal_mdof import payload as plm
from gimbal_mdof import plotting as pt
from gimbal_mdof import report as rp
from gimbal_mdof import simulate as sm
from gimbal_mdof import tf

RADICE = Path(__file__).resolve().parent
RISULTATI = RADICE / "results"
FIG = RISULTATI / "figures"
TAB = RISULTATI / "tables"


def argomenti() -> argparse.Namespace:
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--config", type=Path, default=RADICE / "config.yaml")
    a.add_argument("--theta0", type=float, nargs="+", help="punti di lavoro theta0 [deg]")
    a.add_argument("--omega0", type=float, nargs="+", help="velocita' di inseguimento [deg/s]")
    a.add_argument("--mass-model", choices=["lumped", "consistent", "all"], default="all")
    a.add_argument("--mass-form", choices=["guyan", "guyan_lumped", "all"], default="all")
    a.add_argument("--quick", action="store_true", help="griglie ridotte e simulazioni brevi")
    a.add_argument("--no-latex", action="store_true", help="non generare ne' compilare il LaTeX")
    return a.parse_args()


def log(msg: str, t0: float) -> None:
    print(f"[{time.time() - t0:7.1f} s] {msg}", flush=True)


def main() -> None:
    args = argomenti()
    t0 = time.time()
    par = pr.parametri_da_config(args.config)
    theta0 = tuple(np.deg2rad(args.theta0)) if args.theta0 else par.theta0
    omega0 = tuple(np.deg2rad(args.omega0)) if args.omega0 else par.omega0
    mm_list = ["lumped", "consistent"] if args.mass_model == "all" else [args.mass_model]
    mf_list = ["guyan", "guyan_lumped"] if args.mass_form == "all" else [args.mass_form]
    # il modello di riferimento per figure e LaTeX e' sempre consistent/guyan
    combinazioni = sorted({(m, f) for m in mm_list for f in mf_list} | {("consistent", "guyan"),
                                                                        ("consistent", "guyan_lumped"),
                                                                        ("lumped", "guyan")})
    n_punti = 800 if args.quick else par.n_punti
    seg = par.segnaposto
    pt.applica_stile()
    figure: list[Path] = []
    tabelle: list[Path] = []

    # --- 1-2. elementi e FEM ------------------------------------------------------
    assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
            fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    for mod in assi:
        tabelle.append(geo.scrivi_scomposizione(mod.elementi_originali,
                                                TAB / f"element_breakdown_{mod.nome}.csv"))
    log("assi FEM costruiti", t0)

    # --- payload e sistemi ----------------------------------------------------------
    payload = plm.costruisci_payload(par, assi[1])
    avvisi = list(par.avvisi) + [a for m in assi for a in m.avvisi] + list(payload.avvisi)
    pr.stampa_avvisi(avvisi)
    sistemi = {c: cp.costruisci_sistema(par, *c, payload=payload, assi=assi) for c in combinazioni}
    s = sistemi[("consistent", "guyan")]
    log(f"sistemi ridotti e accoppiati: {len(sistemi)} combinazioni", t0)

    # --- validazione della riduzione -------------------------------------------------
    f_full_max = max(float(np.max(fem.frequenze_proprie(*sistemi[c].isolato_completo(a, 0.0)[:2])))
                     for c in sistemi for a in ("azimut", "polare"))
    w_full = 2 * np.pi * np.logspace(np.log10(par.f_min), np.log10(max(par.f_max, 2 * f_full_max)),
                                     n_punti)
    confronti = [an.confronta_riduzione(sistemi[c], a, 0.0, w_full)
                 for c in sorted(sistemi) for a in ("azimut", "polare")]
    tabelle.append(an.salva_csv(TAB / "reduction_validation.csv", an.intestazione_validazione(confronti),
                                an.righe_validazione(confronti), segnaposto=seg))
    righe = []
    for mm in ("lumped", "consistent"):
        for a in ("azimut", "polare"):
            for th in theta0:
                f = fem.frequenze_proprie(*sistemi[(mm, "guyan")].isolato_completo(a, th)[:2])
                righe.append([a, mm, f"{np.rad2deg(th):g}", *[f"{v:.6g}" for v in f]])
    tabelle.append(an.salva_csv(TAB / "full_model_frequencies.csv",
                                ["asse", "modello_massa", "theta0_deg", "f_1..f_n [Hz]"], righe,
                                segnaposto=seg))
    log("validazione della riduzione", t0)

    # --- linearizzazione, poli, zeri, esportazione --------------------------------------
    lineari: dict[tuple[str, str, float, float], lz.ModelloLineare] = {}
    righe_pz: list[list[object]] = []
    for c in sorted(sistemi):
        for th in theta0:
            for Om in omega0:
                ml = lz.linearizza(sistemi[c].accoppiato, th, Om)
                lineari[(c[0], c[1], th, Om)] = ml
                righe_pz += an.righe_poli_zeri(ml, *c)
                if c == ("consistent", "guyan"):
                    tf.esporta(ml, RISULTATI / "tf",
                               f"ss_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}")
    tabelle.append(an.salva_csv(TAB / "poles_zeros.csv", an.INTESTAZIONE_RADICI, righe_pz,
                                segnaposto=seg))
    righe = [[c[0], c[1], f"{np.rad2deg(th):g}", f"{np.rad2deg(Om):g}", f"{ml.punto.tau0[0]:.6g}",
              f"{ml.punto.tau0[1]:.6g}", *[f"{v:.6e}" for v in ml.punto.q0]]
             for (c0, c1, th, Om), ml in lineari.items() for c in [(c0, c1)]]
    tabelle.append(an.salva_csv(TAB / "operating_points.csv",
                                ["massa", "forma", "theta0_deg", "Omega0_deg_s", "tau_psi0_Nm",
                                 "tau_theta0_Nm", *s.accoppiato.struttura.nomi], righe, segnaposto=seg))
    log(f"linearizzazioni: {len(lineari)}", t0)

    # --- figure --------------------------------------------------------------------
    f_max_acc = max(float(np.max(an.frequenze_modi(ml))) for ml in lineari.values())
    w = 2 * np.pi * np.logspace(np.log10(par.f_min), np.log10(max(par.f_max, 2 * f_max_acc)), n_punti)
    for a in ("azimut", "polare"):
        figure += pt.salva(pt.disegno_geometria(par, a), FIG, f"geometria_{a}")
        figure += pt.salva(pt.forme_modali(s, a, 0.0, segnaposto=seg), FIG, f"modi_{a}")
        figure += pt.salva(pt.bode_confronto_asse({"lumped": sistemi[("lumped", "guyan")],
                                                   "consistent": s}, a, 0.0, w_full, segnaposto=seg),
                           FIG, f"confronto_completo_ridotto_{a}")
    nomi_fig: dict[str, str] = {"confronto_azimut": "confronto_completo_ridotto_azimut",
                                "confronto_polare": "confronto_completo_ridotto_polare",
                                "modi_azimut": "modi_azimut", "modi_polare": "modi_polare"}
    for th in theta0:
        modelli = [(f"{mm}, {mf}", lineari[(mm, mf, th, 0.0)])
                   for mm, mf in (("consistent", "guyan"), ("lumped", "guyan"), ("consistent", "guyan_lumped"))
                   if (mm, mf, th, 0.0) in lineari]
        nome = f"bode_encoder_th{np.rad2deg(th):+.0f}_om0"
        figure += pt.salva(pt.bode_2x2(modelli, w, segnaposto=seg,
                                       titolo=f"Coppia -> encoder, theta0 = {np.rad2deg(th):g} deg, "
                                              "Omega0 = 0"), FIG, nome)
        if th == theta0[0]:
            nomi_fig["bode"] = nome
        Om = max(omega0)
        if Om > 0:
            nome = f"bode_encoder_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}"
            figure += pt.salva(pt.bode_2x2([("consistent, guyan", lineari[("consistent", "guyan", th, Om)])],
                                           w, segnaposto=seg,
                                           titolo=f"Coppia -> encoder, theta0 = {np.rad2deg(th):g} deg, "
                                                  f"Omega0 = {np.rad2deg(Om):g} deg/s"), FIG, nome)
            nome = f"errore_misura_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}"
            figure += pt.salva(pt.bode_errore_misura(lineari[("consistent", "guyan", th, Om)], w,
                                                     segnaposto=seg), FIG, nome)
            nomi_fig.setdefault("errore", nome)
    sweep_grid = par.theta0_sweep if not args.quick else par.theta0_sweep[::3]
    sw = an.sweep_theta0(s, sweep_grid, omega0)
    figure += pt.salva(pt.frequenze_vs_theta0(sw, segnaposto=seg), FIG, "frequenze_vs_theta0")
    nomi_fig["sweep"] = "frequenze_vs_theta0"
    th_x = theta0[len(theta0) // 2]
    G = an.frf_incrociate(s, th_x, sorted(set(omega0) | {0.0, max(omega0) / 2}), w)
    figure += pt.salva(pt.frf_incrociate_omega(G, w, th_x, segnaposto=seg), FIG, "frf_incrociate_vs_omega")
    log(f"figure: {len(figure) // 2}", t0)

    # --- simulazioni ------------------------------------------------------------------
    pl = s.payload
    stabili = [th for th in theta0 if pl.k_g(th) > 0] or [np.deg2rad(15.0)]
    th_sim = min(stabili, key=lambda t: abs(t - np.deg2rad(15.0)))
    t_fine = 0.1 if args.quick else 0.3
    simulazioni = [sm.S1(s.accoppiato, th_sim, t_fine=t_fine, n_punti=int(t_fine * 5000) + 1),
                   sm.S2(s.accoppiato, th_sim, par.omega_max, t_fine=t_fine,
                         n_punti=int(t_fine * 5000) + 1),
                   sm.S3(s.accoppiato, th_sim, par.omega_max, t_fine=0.1 if args.quick else 0.2,
                         n_punti=401)]
    for r in simulazioni[:2]:
        figure += pt.salva(pt.storie_temporali(r, segnaposto=seg), FIG, f"simulazione_{r.nome}")
    figure += pt.salva(pt.conservazione(simulazioni[2], segnaposto=seg), FIG, "simulazione_S3")
    log("simulazioni S1, S2, S3", t0)

    # --- LaTeX ------------------------------------------------------------------------
    testo_latex = ""
    if not args.no_latex:
        cont = lr.ContenutoLatex(
            par=par, sistemi=sistemi, confronti=confronti,
            modelli_lineari=[lineari[k] for k in lineari if k[:2] == ("consistent", "guyan")],
            figure={k: f"../figures/{v}.pdf" for k, v in nomi_fig.items()},
            avvisi=avvisi, ml_tabella=lineari[("consistent", "guyan", theta0[0], 0.0)])
        tex = lr.genera(cont, RISULTATI / "latex")
        if par.latex_compila:
            esito = lr.compila(tex)
            testo_latex = (f"`{tex.relative_to(RADICE)}`: {esito.messaggio}"
                           + (f", PDF `{esito.pdf.relative_to(RADICE)}`" if esito.pdf else "")
                           + (". ATTENZIONE: riferimenti non risolti." if esito.riferimenti_irrisolti else "")
                           + (f" Errori: {list(esito.errori)}" if esito.errori else ""))
        else:
            testo_latex = f"`{tex.relative_to(RADICE)}` generato (compilazione disattivata)."
        log(f"LaTeX: {testo_latex}", t0)

    # --- report -----------------------------------------------------------------------
    percorso = rp.genera(rp.ContenutoReport(
        par=par, sistemi=sistemi, confronti=confronti, lineari=lineari, simulazioni=simulazioni,
        figure=[f.relative_to(RADICE) for f in figure if f.suffix == ".pdf"],
        tabelle=[t.relative_to(RADICE) for t in tabelle], avvisi=avvisi, latex=testo_latex),
        RISULTATI / "report.md")
    log(f"report: {percorso.relative_to(RADICE)}", t0)


if __name__ == "__main__":
    main()
