"""Esegue l'intera catena del modello multi-GdL e scrive i risultati in ``results/``.

Esempi::

    .venv/bin/python run_all.py
    .venv/bin/python run_all.py --quick --no-latex
    .venv/bin/python run_all.py --theta0 -15 15 --omega0 0 10 --mass-model consistent
    .venv/bin/python run_all.py --reduction guyan --f-band 1500
    .venv/bin/python run_all.py --hub    # modello hub (config_hub.yaml -> results_hub/)

La riduzione del modello principale (linearizzazioni, esportazione, simulazioni) e'
``reduction.method`` di config.yaml (predefinita: troncamento modale), sovrascrivibile
con ``--reduction``; Guyan resta sempre calcolato come confronto.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import replace
from types import SimpleNamespace
from pathlib import Path

import numpy as np

from gimbal_mdof import analysis as an
from gimbal_mdof import concentrato as cc
from gimbal_mdof import concentrato_analisi as ca
from gimbal_mdof import concentrato_output as co
from gimbal_mdof import simulatore_concentrato as scm
from gimbal_mdof import sperimentale as spm
from gimbal_mdof import sperimentale_output as spo
from gimbal_mdof import coupled as cp
from gimbal_mdof import fem
from gimbal_mdof import geometry as geo
from gimbal_mdof import latex_report as lr
from gimbal_mdof import linearize as lz
from gimbal_mdof import modal as md
from gimbal_mdof import modal_output as mo
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
    a.add_argument("--hub", action="store_true",
                   help="modello hub (mozzo polare e telescopio appeso tramite k_b): usa config_hub.yaml, "
                        "risultati in results_hub/; senza il flag, modello originale (config.yaml, results/)")
    a.add_argument("--theta0", type=float, nargs="+", help="punti di lavoro theta0 [deg]")
    a.add_argument("--omega0", type=float, nargs="+", help="velocita' di inseguimento [deg/s]")
    a.add_argument("--mass-model", choices=["lumped", "consistent", "all"], default="all")
    a.add_argument("--mass-form", choices=["guyan", "guyan_lumped", "all"], default="all")
    a.add_argument("--quick", action="store_true", help="griglie ridotte e simulazioni brevi")
    a.add_argument("--no-latex", action="store_true", help="non generare ne' compilare il LaTeX")
    a.add_argument("--reduction", choices=["modal", "guyan", "craig_bampton"],
                   help="riduzione del modello principale (default: reduction.method di config.yaml)")
    a.add_argument("--f-band", type=float, help="banda di interesse [Hz] (reduction.modal.f_band_Hz e "
                                                 "limite delle griglie di frequenza)")
    return a.parse_args()


def log(msg: str, t0: float) -> None:
    print(f"[{time.time() - t0:7.1f} s] {msg}", flush=True)


def _dati_modali(ml: lz.ModelloLineare) -> dict[str, object]:
    """Grandezze modali da esportare con lo stato: basi, frequenze, smorzamenti, vettori
    d'aggancio e stato a blocchi di Gawronski (eq. 2.25/2.32)."""
    acc = ml.accoppiato
    sg = md.stato_gawronski(ml)
    out: dict[str, object] = {
        "forma": "modale (Gawronski)", "A_gawronski": sg.A, "B_gawronski": sg.B, "C_gawronski": sg.C,
        "S_gawronski": sg.S, "stati_gawronski": np.array(sg.nomi), "E_eta": acc.E_eta,
        "M_tilde": acc.M_tilde, "K_tilde": acc.K_tilde, "D_tilde": acc.D_tilde,
        "theta_ref": float(ml.theta0)}
    for b in (acc.base_psi, acc.base_theta):
        s = "psi" if b.asse == "azimut" else "theta"
        out |= {f"Phi_{s}": b.T, f"omega_{s}": b.omega, f"zeta_{s}": b.zeta,
                f"modi_fem_{s}": np.array(b.modi_fem)}
    return out


def main() -> None:
    args = argomenti()
    if args.hub:
        args.config = RADICE / "config_hub.yaml"
    t0 = time.time()
    par = pr.parametri_da_config(args.config)
    # cartella dei risultati della config (results/ per il modello originale, results_hub/ per il
    # modello hub): i due modelli non si sovrascrivono
    global RISULTATI, FIG, TAB
    RISULTATI = RADICE / par.cartella_risultati
    FIG, TAB = RISULTATI / "figures", RISULTATI / "tables"
    if args.f_band:
        par = replace(par, f_max=args.f_band, riduzione=replace(
            par.riduzione, modale=replace(par.riduzione.modale, f_banda=args.f_band)))
    riduzione = args.reduction or par.riduzione.metodo
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

    # --- troncamento modale (sempre: e' il confronto di Guyan e il modello predefinito) ---
    s_mod = md.costruisci_sistema_modale(par, "consistent", payload=payload, assi=assi,
                                         punti=[(th, Om) for th in theta0 for Om in omega0])
    avvisi += s_mod.avvisi()
    ris_mod = mo.calcola(s_mod)
    tabelle += mo.tabelle(ris_mod, TAB, seg)
    s_conc = cc.costruisci_sistema_concentrato(par, payload, assi)
    variante = "mozzo" if par.mozzo is not None else "catena"
    cc.verifica_contro_lagrange(variante)  # Newton-Eulero = Lagrange (identita' simbolica)
    log(f"modello a parametri concentrati ({'3 + 5 GdL, telescopio appeso al mozzo' if par.mozzo else '3 + 4 GdL'}): "
        "Newton-Eulero verificato contro Lagrange", t0)
    log(f"troncamento modale: {ris_mod.stati} stati, {s_mod.n_aggiunti} modi aggiunti dalla validazione, "
        f"T5 {'rispettato' if all(e.passa for e in s_mod.validazione) else 'NON rispettato'}", t0)
    if riduzione == "modal":
        lin_pr, acc_pr, sistema_pr = s_mod.lineare, s_mod.modello(), s_mod
    elif riduzione == "craig_bampton":
        cb = s_mod.craig_bampton(par.riduzione.n_modi_interfaccia)
        lin_pr, acc_pr, sistema_pr = cb.linearizza, cb, SimpleNamespace(lineare=cb.linearizza)
    else:
        lin_pr = lambda th, Om: lz.linearizza(s.accoppiato, th, Om)  # noqa: E731
        acc_pr, sistema_pr = s.accoppiato, s

    # --- validazione della riduzione (griglia fino a f_max, senza estensioni) -----------
    w_full = 2 * np.pi * np.logspace(np.log10(par.f_min), np.log10(par.f_max), n_punti)
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
                if c == ("consistent", "guyan") and riduzione == "guyan":
                    tf.esporta(ml, RISULTATI / "tf",
                               f"ss_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}")
    lineari_pr: dict[tuple[float, float], lz.ModelloLineare] = {}
    for th in theta0:
        for Om in omega0:
            ml = lin_pr(th, Om)
            lineari_pr[(th, Om)] = ml
            nome = f"ss_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}"
            if riduzione != "guyan":
                righe_pz += an.righe_poli_zeri(ml, "consistent", riduzione)
                tf.esporta(ml, RISULTATI / "tf", nome, extra=_dati_modali(ml) if riduzione == "modal" else None)
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
    w = w_full
    for a in ("azimut", "polare"):
        figure += pt.salva(pt.disegno_geometria(par, a), FIG, f"geometria_{a}")
        figure += pt.salva(pt.forme_modali(s, a, 0.0, segnaposto=seg), FIG, f"modi_{a}")
        figure += pt.salva(pt.bode_confronto_asse({"lumped": sistemi[("lumped", "guyan")],
                                                   "consistent": s}, a, 0.0, w_full, segnaposto=seg),
                           FIG, f"confronto_completo_ridotto_{a}")
    nomi_fig: dict[str, str] = {"confronto_azimut": "confronto_completo_ridotto_azimut",
                                "confronto_polare": "confronto_completo_ridotto_polare",
                                "modi_azimut": "modi_azimut", "modi_polare": "modi_polare"}
    if par.mozzo is not None:  # modello hub: geometria del mozzo usata per le rigidezze
        from gimbal_mdof import mozzo as mz

        figure += pt.salva(pt.disegno_mozzo(par.mozzo, mz.rigidezze(par.mozzo), segnaposto=seg), FIG,
                           "geometria_mozzo")
        nomi_fig["mozzo"] = "geometria_mozzo"
    for th in theta0:
        modelli = [(f"{mm}, {mf}", lineari[(mm, mf, th, 0.0)])
                   for mm, mf in (("consistent", "guyan"), ("lumped", "guyan"), ("consistent", "guyan_lumped"))
                   if (mm, mf, th, 0.0) in lineari]
        nome = f"bode_encoder_th{np.rad2deg(th):+.0f}_om0"
        figure += pt.salva(pt.bode_2x2(modelli, w, segnaposto=seg,
                                       titolo=f"Torque -> encoder, theta0 = {np.rad2deg(th):g} deg, "
                                              "Omega0 = 0"), FIG, nome)
        if th == theta0[0]:
            nomi_fig["bode"] = nome
        Om = max(omega0)
        if Om > 0:
            nome = f"bode_encoder_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}"
            figure += pt.salva(pt.bode_2x2([("consistent, guyan", lineari[("consistent", "guyan", th, Om)])],
                                           w, segnaposto=seg,
                                           titolo=f"Torque -> encoder, theta0 = {np.rad2deg(th):g} deg, "
                                                  f"Omega0 = {np.rad2deg(Om):g} deg/s"), FIG, nome)
            nome = f"errore_misura_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}"
            figure += pt.salva(pt.bode_errore_misura(lineari_pr[(th, Om)], w,  # modello principale
                                                     segnaposto=seg), FIG, nome)
            nomi_fig.setdefault("errore", nome)
    for th in theta0:
        nome = f"bode_modale_th{np.rad2deg(th):+.0f}_om0"
        modelli = [(f"modal ({2 * s_mod.lineare(th, 0.0).n} states)", s_mod.lineare(th, 0.0)),
                   (f"full ({s_mod.completo().n} DOF)", s_mod.lineare_completo(th, 0.0))]
        if par.latex_guyan:  # latex.include_guyan: la curva di Guyan solo se richiesta
            modelli.append(("Guyan (comparison)", lineari[("consistent", "guyan", th, 0.0)]))
        figure += pt.salva(pt.bode_2x2(modelli, w, segnaposto=seg,
                                       titolo=f"Torque -> encoder, theta0 = {np.rad2deg(th):g} deg, Omega0 = 0: "
                                              "modal truncation and full model"
                                              + (" (with Guyan)" if par.latex_guyan else "")), FIG, nome)
        nomi_fig.setdefault("bode_modale", nome)
    for a in ("azimut", "polare"):
        figure += pt.salva(pt.forme_modali_modale(s_mod, a, segnaposto=seg), FIG, f"modi_modale_{a}")
        nomi_fig[f"modi_modale_{a}"] = f"modi_modale_{a}"
    sweep_grid = par.theta0_sweep if not args.quick else par.theta0_sweep[::3]
    sw = an.sweep_theta0(sistema_pr, sweep_grid, omega0)
    figure += pt.salva(pt.frequenze_vs_theta0(sw, segnaposto=seg), FIG, "frequenze_vs_theta0")
    nomi_fig["sweep"] = "frequenze_vs_theta0"
    th_x = theta0[len(theta0) // 2]
    G = an.frf_incrociate(sistema_pr, th_x, sorted(set(omega0) | {0.0, max(omega0) / 2}), w)
    figure += pt.salva(pt.frf_incrociate_omega(G, w, th_x, segnaposto=seg), FIG, "frf_incrociate_vs_omega")
    log(f"figure: {len(figure) // 2}", t0)

    # --- simulazioni ------------------------------------------------------------------
    pl = s.payload
    stabili = [th for th in theta0 if pl.k_g(th) > 0] or [np.deg2rad(15.0)]
    th_sim = min(stabili, key=lambda t: abs(t - np.deg2rad(15.0)))
    t_fine = 0.1 if args.quick else 0.3
    simulazioni = [sm.S1(acc_pr, th_sim, t_fine=t_fine, n_punti=int(t_fine * 20000) + 1),
                   sm.S2(acc_pr, th_sim, par.omega_max, t_fine=t_fine,
                         n_punti=int(t_fine * 20000) + 1),
                   sm.S3(acc_pr, th_sim, par.omega_max, t_fine=0.1 if args.quick else 0.2,
                         n_punti=401)]
    # riferimento a tutti i GdL del FEM: modi fino a centinaia di kHz, DOP853 limitato dalla
    # stabilita' -> tolleranze 1e-8/1e-10 e finestra piu' corta (stesso passo di campionamento)
    t_rif = 0.05 if args.quick else 0.1
    s1_rif = sm.S1(s_mod.completo(), th_sim, t_fine=t_rif, n_punti=int(t_rif * 20000) + 1,
                   rtol=1e-8, atol=1e-10)
    for r in simulazioni[:2]:
        figure += pt.salva(pt.storie_temporali(r, segnaposto=seg), FIG, f"simulazione_{r.nome}")
    figure += pt.salva(pt.conservazione(simulazioni[2], segnaposto=seg), FIG, "simulazione_S3")
    nomi_fig.update({"S1": "simulazione_S1", "S2": "simulazione_S2", "S3": "simulazione_S3"})
    log("simulazioni S1, S2, S3", t0)

    # --- modello a parametri concentrati (3 + 4 GdL) ---------------------------------
    punti = [(th, Om) for th in theta0 for Om in omega0]
    an_conc = ca.analizza(s_conc, s_mod, punti)
    for (th, Om), ml in an_conc.lineari.items():
        righe_pz += an.righe_poli_zeri(ml, "consistent", cc.FORMA)
        tf.esporta(ml, RISULTATI / "tf_concentrato", f"ss_th{np.rad2deg(th):+.0f}_om{np.rad2deg(Om):.0f}")
    an.salva_csv(TAB / "poles_zeros.csv", an.INTESTAZIONE_RADICI, righe_pz, segnaposto=seg)  # con il concentrato
    righe_el = [[e.asse, e.nome, e.nodo_a, e.nodo_b, " ".join(e.elementi), f"{e.J:.6e}", f"{e.k:.6e}"]
                for m in assi for e in cc.elementi_equivalenti(m)]
    righe_el += [["nodo", r, rid.etichette[i], "", "", f"{rid.M_r[i, i]:.6e}", ""]
                 for rid, ruoli in ((s_conc.rid_psi, cc.RUOLI_PSI), (s_conc.rid_theta, cc.RUOLI_THETA))
                 for i, r in enumerate(ruoli)]
    tabelle.append(an.salva_csv(TAB / "lumped_parameters.csv", ["asse_o_nodo", "elemento_o_ruolo", "nodo_a",
                                                                "nodo_b", "elementi_FEM", "J_kgm2", "k_Nm_rad"],
                                righe_el, segnaposto=seg))
    tabelle.append(an.salva_csv(TAB / "lumped_vs_full.csv", ["voce", "theta0_deg", "Omega0_deg_s", "grandezza",
                                                             "valore"], ca.righe_csv(an_conc), segnaposto=seg))
    fig_conc: dict[str, str] = {}
    for th in theta0:
        nome = f"bode_concentrato_th{np.rad2deg(th):+.0f}_om0"
        modelli = [("lumped (3 + 4 DOF)", an_conc.lineari[(th, 0.0)] if (th, 0.0) in an_conc.lineari
                    else lz.linearizza(s_conc.accoppiato, th, 0.0)),
                   (f"full ({s_mod.completo().n} DOF)", s_mod.lineare_completo(th, 0.0)),
                   (f"modal ({2 * s_mod.lineare(th, 0.0).n} states)", s_mod.lineare(th, 0.0))]
        figure += pt.salva(pt.bode_2x2(modelli, w, segnaposto=seg,
                                       titolo=f"Torque -> encoder, theta0 = {np.rad2deg(th):g} deg, Omega0 = 0: "
                                              "lumped-parameter model"), FIG, nome)
        fig_conc.setdefault("bode_concentrato", f"../figures/{nome}.pdf")
    sw_c = an.sweep_theta0(s_conc, sweep_grid, omega0)
    figure += pt.salva(pt.frequenze_vs_theta0(sw_c, segnaposto=seg), FIG, "frequenze_vs_theta0_concentrato")
    fig_conc["sweep_concentrato"] = "../figures/frequenze_vs_theta0_concentrato.pdf"
    # catena polare disegnata per i canali coppia -> encoder e coppia -> motore (PNG nel LaTeX)
    if par.mozzo is None:
        fig_catene = pt.catene_poli_zeri(ca.catena_polare(s_conc, 0.0), segnaposto=seg)
    else:  # modello hub: telescopio appeso al centro del mozzo
        fig_catene = pt.catena_ramo_poli_zeri(ca.catena_ramo_modello(s_conc, 0.0), segnaposto=seg, stimata=False)
    figure += pt.salva(fig_catene, FIG, "catene_polari_poli_zeri")
    fig_conc["catene_polari"] = "../figures/catene_polari_poli_zeri.png"
    sim_c = scm.costruisci_simulatore(s_conc)
    file_codice = [p.relative_to(RADICE).as_posix()
                   for p in scm.esporta_codice(sim_c, RISULTATI / "simulatore_concentrato")]
    t4 = 1.5 if args.quick else 3.0
    sim_conc = [sm.S1(sim_c, th_sim, t_fine=t_fine, n_punti=int(t_fine * 20000) + 1),
                sm.S3(sim_c, th_sim, par.omega_max, t_fine=0.1 if args.quick else 0.2, n_punti=401)]
    for chiave, (th4, dth) in {"S4a": (-45.0, 2.0), "S4b": (-45.0, 10.0), "S4c": (0.0, 10.0)}.items():
        r4 = scm.S4(sim_c, np.deg2rad(th4), par.omega_max, delta_theta=np.deg2rad(dth), t_fine=t4,
                    n_punti=int(t4 * 1000) + 1)
        r4 = replace(r4, nome=chiave)
        sim_conc.append(r4)
        figure += pt.salva(pt.storie_grande_ampiezza(r4, segnaposto=seg), FIG, f"simulazione_concentrato_{chiave}")
        fig_conc[chiave] = f"../figures/simulazione_concentrato_{chiave}.pdf"
    figure += pt.salva(pt.storie_temporali(sim_conc[0], segnaposto=seg), FIG, "simulazione_concentrato_S1")
    log(f"modello concentrato: FdT, {len(sim_conc)} simulazioni, codice in {par.cartella_risultati}/simulatore_concentrato", t0)

    # --- confronto con le misure sperimentali (solo il modello a 7 GdL concentrato) -----
    tex_sper, md_sper = "", []
    if par.sperimentale is not None:
        cfg_s = par.sperimentale
        sis_s = {True: s_conc, False: cc.costruisci_sistema_concentrato(
            par, spm.payload_senza_telescopio(payload), spm.assi_senza_telescopio(par, assi))}
        misure = [spm.carica_misura(cfg_s, v) for v in cfg_s.misure]
        f_s = np.logspace(np.log10(cfg_s.f_min), np.log10(cfg_s.f_max), n_punti)
        confronti_s = [spm.confronta(par, m, sis_s, f_s) for m in misure]
        fig_s: dict[str, str] = {}
        for c in confronti_s:
            nome = f"exp_{Path(c.misura.file).stem}"
            figure += pt.salva(pt.bode_misura(c, segnaposto=seg, f_lim=(cfg_s.f_min, cfg_s.f_max)), FIG, nome)
            fig_s[nome] = f"../figures/{nome}.pdf"
        az = [m for m in misure if m.asse == "azimut"]
        if len(az) > 1:
            figure += pt.salva(pt.bode_misure_sovrapposte(az, titolo="Azimuth-axis measurements compared",
                                                          segnaposto=False, f_lim=(cfg_s.f_min, cfg_s.f_max)),
                               FIG, "exp_azimut_sovrapposte")
            fig_s["exp_azimut_sovrapposte"] = "../figures/exp_azimut_sovrapposte.pdf"
        tabelle.append(an.salva_csv(TAB / "experimental_comparison.csv", spo.INTESTAZIONE_CSV,
                                    spo.righe_csv(confronti_s), segnaposto=seg))
        disc = spo.discussione(confronti_s, sis_s, cfg_s.Kt)
        if par.mozzo is None and disc.el is not None and disc.el.catena is not None:
            # catena polare con il telescopio appeso al mozzo tramite k_b (stima dalla misura)
            figure += pt.salva(pt.catena_ramo_poli_zeri(disc.el.catena, segnaposto=seg), FIG,
                               "catena_polare_ramo_laterale")
            fig_s["catena_ramo"] = "../figures/catena_polare_ramo_laterale.png"
        tex_sper = spo.latex(confronti_s, disc, fig_s, cfg_s.Kt, seg,
                             str(cfg_s.cartella.relative_to(RADICE.parent)) if cfg_s.cartella.is_relative_to(
                                 RADICE.parent) else str(cfg_s.cartella))
        md_sper = spo.markdown(confronti_s, disc)
        log(f"confronto sperimentale: {len(misure)} misure", t0)

    # --- LaTeX ------------------------------------------------------------------------
    testo_latex = ""
    if not args.no_latex:
        cont = lr.ContenutoLatex(
            par=par, sistemi=sistemi, confronti=confronti,
            modelli_lineari=[lineari[k] for k in lineari if k[:2] == ("consistent", "guyan")],
            figure={k: f"../figures/{v}.pdf" for k, v in nomi_fig.items()},
            avvisi=avvisi, ml_tabella=lineari[("consistent", "guyan", theta0[0], 0.0)],
            simulazioni=simulazioni, riduzione=riduzione, modale=ris_mod, s1_riferimento=s1_rif,
            includi_guyan=par.latex_guyan, ml_principale=lineari_pr[(theta0[0], omega0[0])],
            concentrato=s_conc,
            concentrato_extra=[ca.latex(s_conc, an_conc, fig_conc, seg),
                               co.latex_simulatore(sim_conc, fig_conc, seg, file_codice)],
            sperimentale_tex=tex_sper)
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
        tabelle=[t.relative_to(RADICE) for t in tabelle], avvisi=avvisi, latex=testo_latex,
        riduzione=riduzione, modale=ris_mod, lineari_principali=lineari_pr, s1_riferimento=s1_rif,
        extra_md=co.markdown(s_conc, an_conc, sim_conc, seg, file_codice) + md_sper),
        RISULTATI / "report.md")
    log(f"report: {percorso.relative_to(RADICE)}", t0)


if __name__ == "__main__":
    main()
