"""Uscite del troncamento modale: tabelle CSV, sezione di ``report.md``, sezione LaTeX.

Ogni tabella e ogni sezione riportano, se i dati sono segnaposto, l'avviso
"risultato da parametri segnaposto". Le grandezze sono quelle di
:mod:`gimbal_mdof.modal`; le equazioni citate sono del documento
``MODELLO_unificato_gawronski_dedonatis.md`` (§4-§8) e di Gawronski (2008).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from . import modal as md
from .analysis import salva_csv
from .linearize import ModelloLineare
from .params import NOTA_REALISMO, in_inglese

MOTIVO_TENUTO: dict[str, str] = {"forzato": "forced", "taglio": "cut-off", "hankel": "Hankel",
                                 "validazione": "validation"}
"""Motivo per cui un modo e' tenuto, in inglese per il report."""

NOMI_AGGANCIO = {"azimut": "c_L", "polare": "c_P"}


@dataclass(frozen=True)
class RisultatiModali:
    """Tutto cio' che serve a tabelle, report e LaTeX, calcolato una volta sola."""

    sistema: md.SistemaModale
    analisi: dict[str, md.AnalisiModaleAsse]
    analisi_lumped: dict[str, np.ndarray]
    """Frequenze (Hz) di tutti i modi con massa concentrata, per il controllo di maglia."""
    residui: tuple[md.ResiduoAsse, md.ResiduoAsse]
    maglia: tuple[md.ControlloMaglia, ...]
    basi_correzione: dict[str, md.BaseAsse]
    """Basi con vettori residui (per riportare lo spostamento dei modi tenuti)."""
    stati: int
    stati_guyan: int
    validazione: tuple[md.EsitoConfronto, ...]
    confronto_guyan: tuple[md.EsitoConfronto, ...]
    """Stessi criteri T5 applicati al modello di Guyan (confronto)."""


def calcola(s: md.SistemaModale, *, guyan: bool = True) -> RisultatiModali:
    from dataclasses import replace

    an = {a: s.analisi(a) for a in ("azimut", "polare")}
    lumped = {}
    for k, a in enumerate(("azimut", "polare")):
        ast_l = md.assi_strutturali(s.par, s.payload, s.assi_fem, "lumped")[k]
        lumped[a] = md.analizza_asse(ast_l, s.payload, an[a].theta_ref, s.pm).f
    s_rv = replace(s, pm=replace(s.pm, correzione_statica="residual_vectors"))
    basi = {a: s_rv.base(a) for a in ("azimut", "polare")}
    conf_g: list[md.EsitoConfronto] = []
    if guyan:
        mg = s.guyan()
        for e in s.validazione:
            conf_g.append(md.confronta(mg.linearizza(e.theta0, e.Omega0),
                                       s.lineare_completo(e.theta0, e.Omega0), s.par.f_min, s.pm.f_banda))
    return RisultatiModali(sistema=s, analisi=an, analisi_lumped=lumped, residui=s.residui(),
                           maglia=s.controllo_maglia(), basi_correzione=basi, stati=2 * s.modello().n,
                           stati_guyan=2 * s.guyan().n, validazione=s.validazione,
                           confronto_guyan=tuple(conf_g))


# ---------------------------------------------------------------------------
# Tabelle CSV
# ---------------------------------------------------------------------------


def tabelle(r: RisultatiModali, cartella: Path, segnaposto: bool) -> list[Path]:
    s = r.sistema
    out: list[Path] = []
    for asse, an in r.analisi.items():
        f_l = r.analisi_lumped[asse]
        mm = an.masse_modali
        righe = []
        for k in range(an.Phi.shape[1]):
            righe.append([k, "rigido" if k == 0 else "flessibile", f"{an.f[k]:.6g}", f"{an.zeta[k]:.4g}",
                          f"{mm[k]:.6g}", "kg m^2" if k == 0 else "- (normalizzato)",
                          f"{an.c_e[k]:.6g}", f"{an.c_m[k]:.6g}", f"{an.c_agg[k]:.6g}",
                          "inf" if not np.isfinite(an.hankel[k]) else f"{an.hankel[k]:.4e}",
                          "si" if k == 0 or k in an.tenuti else "no",
                          "sempre (rigido)" if k == 0 else an.motivo.get(k, ""),
                          f"{f_l[k]:.6g}", "" if k == 0 else f"{f_l[k] / an.f[k] - 1:+.3e}"])
        out.append(salva_csv(cartella / f"modal_modes_{asse}.csv",
                             ["modo", "tipo", "f_Hz", "zeta", "massa_modale", "unita_massa", "c_e", "c_m",
                              NOMI_AGGANCIO[asse], "hankel_rad_Nm", "tenuto", "motivo", "f_lumped_Hz",
                              "scarto_lumped"], righe, segnaposto=segnaposto))
    righe = [[i + 1, f"{f:.6g}", a, k, "rigido" if k == 0 else "flessibile", "si" if t else "no"]
             for i, (f, a, k, t) in enumerate(s.spettro())]
    out.append(salva_csv(cartella / "modal_spectrum.csv",
                         ["ordine", "f_Hz", "asse", "modo", "tipo", "tenuto"], righe, segnaposto=segnaposto))
    righe = []
    for nome, esiti in (("modale", r.validazione), ("guyan", r.confronto_guyan)):
        for e in esiti:
            base = [nome, f"{np.rad2deg(e.theta0):g}", f"{np.rad2deg(e.Omega0):g}"]
            righe += [base + ["err_dB tau_psi->psi_e", f"{e.err_db[0]:.4e}"],
                      base + ["err_dB tau_theta->theta_e", f"{e.err_db[1]:.4e}"],
                      base + ["err_deg tau_psi->psi_e", f"{e.err_deg[0]:.4e}"],
                      base + ["err_deg tau_theta->theta_e", f"{e.err_deg[1]:.4e}"],
                      base + ["err_rel_2x2", f"{e.err_incrociato:.4e}"],
                      base + ["passa_T5", "si" if e.passa else "no"]]
            righe += [base + [f"zero {c} {fr:.6g} Hz: f ridotto", f"{fm:.6g} ({err:.2e})"]
                      for c, fr, fm, err in e.zeri]
    for res in r.residui:
        righe += [["flessibilita", res.asse, "", "F_collocata_scartati_rad_Nm", f"{res.F_collocata:.6e}"],
                  ["flessibilita", res.asse, "", "F_non_collocata_scartati_rad_Nm", f"{res.F_non_collocata:.6e}"],
                  ["flessibilita", res.asse, "", "F_collocata_totale_rad_Nm", f"{res.F_collocata_tot:.6e}"],
                  ["flessibilita", res.asse, "", "F_non_collocata_totale_rad_Nm",
                   f"{res.F_non_collocata_tot:.6e}"]]
        righe += [["zero_asse_isolato", res.asse, "", f"{fdt} completo {fz:.6g} Hz",
                   f"puro {fm:.6g} ({fm / fz - 1:+.2e}); con correzione {fc:.6g} ({fc / fz - 1:+.2e})"]
                  for fdt, fz, fm, fc in res.zeri]
    for a, b in r.basi_correzione.items():
        if b.spostamenti is not None:
            righe.append(["vettori_residui", a, "", "spostamento_max_modi_tenuti",
                          f"{np.max(np.abs(b.spostamenti)) if b.spostamenti.size else 0.0:.3e}"])
    righe.append(["ciclo_validazione", "", "", "modi_aggiunti", str(s.n_aggiunti)])
    righe += [["ciclo_validazione", a, "", f"aggiunto modo {k}", f"{f:.6g} Hz"] for a, k, f in s.storia_aggiunte]
    righe += [["stati", "", "", "modale", str(r.stati)], ["stati", "", "", "guyan", str(r.stati_guyan)]]
    out.append(salva_csv(cartella / "modal_vs_full.csv", ["voce", "asse_o_modello", "theta0_deg_o_vuoto",
                                                          "grandezza", "valore"], righe, segnaposto=segnaposto))
    return out


# ---------------------------------------------------------------------------
# report.md
# ---------------------------------------------------------------------------


def _g(x: float, c: int = 4) -> str:
    return "--" if x is None or not np.isfinite(x) else f"{x:.{c}g}"


def _tab(intest: Sequence[str], righe: Sequence[Sequence[object]]) -> list[str]:
    out = ["| " + " | ".join(intest) + " |", "|" + "|".join("---" for _ in intest) + "|"]
    return out + ["| " + " | ".join(str(c) for c in rr) + " |" for rr in righe] + [""]


def markdown(r: RisultatiModali, segnaposto: bool) -> list[str]:
    s = r.sistema
    pm = s.pm
    si = lambda v: "yes" if v else "no"  # noqa: E731
    L = ["## Gawronski modal truncation (default reduction)", ""]
    if segnaposto:
        L += [f"> Every number in this section is a **{NOTA_REALISMO}**.", ""]
    L += [f"Modal basis from the FULL FEM of each axis with the payload frozen at theta_ref = "
          f"{np.rad2deg(pm.theta_ref):g} deg (k_g = 0); rigid mode 1 not normalised (eta_0 = angle, modal "
          f"mass J_tot), flexible modes mass-normalised. Band {pm.f_banda:g} Hz, cut-off "
          f"{pm.f_taglio:g} Hz, hankel_min_ratio = {pm.hankel_min}, n_flex = { {in_inglese(k): v for k, v in dict(pm.n_flessibili).items()} }, "
          f"static correction: {pm.correzione_statica}.", ""]
    for asse, an in r.analisi.items():
        L += [f"### {in_inglese(asse).capitalize()} axis (J_tot = {an.J_tot:.5g} kg m^2, beta = {an.beta:.4e} s)", ""]
        righe = []
        for k in range(min(an.Phi.shape[1], 8)):
            righe.append([k, _g(an.f[k], 6), _g(an.zeta[k], 3), _g(an.c_e[k]), _g(an.c_m[k]), _g(an.c_agg[k]),
                          "inf" if not np.isfinite(an.hankel[k]) else f"{an.hankel[k]:.3e}",
                          "**yes**" if k == 0 or k in an.tenuti else "no",
                          "rigid" if k == 0 else MOTIVO_TENUTO.get(an.motivo.get(k, ""), an.motivo.get(k, ""))])
        L += _tab(["mode", "f [Hz]", "zeta", "c_e", "c_m", NOMI_AGGANCIO[asse], "Hankel [rad/(N m)]",
                   "kept", "reason"], righe)
    L += ["### Spectrum of the two axes (rigid modes first)", ""]
    L += _tab(["#", "f [Hz]", "axis", "mode", "kept"],
              [[i + 1, _g(f, 6), in_inglese(a), "rigid" if k == 0 else f"flexible {k}", si(t)]
               for i, (f, a, k, t) in enumerate(s.spettro()[:12])])
    L += [f"States of the modal model: **{r.stati}** (Guyan: {r.stati_guyan}). Modes added by the "
          f"validation loop: **{s.n_aggiunti}**"
          + (": " + ", ".join(f"{in_inglese(a)} mode {k} ({f:.1f} Hz)" for a, k, f in s.storia_aggiunte)
             if s.storia_aggiunte else "") + ".", ""]
    L += ["### Comparison with the all-DOF model (T5 criteria, up to "
          f"{pm.f_banda:g} Hz)", "",
          "Direct channels: largest deviation outside +/-2 % around the reference zeros; zeros: largest "
          "relative frequency deviation; 2x2: max ||G - G_ref||_2 / ||G_ref||_2.", ""]
    righe = []
    for e, eg in zip(r.validazione, r.confronto_guyan or [None] * len(r.validazione)):
        z = max((err for *_, err in e.zeri), default=0.0)
        riga = [f"{np.rad2deg(e.theta0):g}", f"{np.rad2deg(e.Omega0):g}",
                f"{max(e.err_db):.3g} dB / {max(e.err_deg):.3g} deg", f"{z:.2e}", f"{e.err_incrociato:.2e}",
                "yes" if e.passa else "**no**"]
        if eg is not None:
            zg = max((err for *_, err in eg.zeri), default=0.0)
            riga += [f"{max(eg.err_db):.3g} dB / {max(eg.err_deg):.3g} deg", f"{zg:.2e}", si(eg.passa)]
        righe.append(riga)
    intest = ["theta0 [deg]", "Omega0 [deg/s]", "modal: direct", "modal: zeros", "modal: 2x2", "modal T5"]
    if r.confronto_guyan:
        intest += ["Guyan: direct", "Guyan: zeros", "Guyan T5"]
    L += _tab(intest, righe)
    L += ["### Residual flexibility of the discarded modes and zeros of the isolated axis", ""]
    righe = []
    for res in r.residui:
        righe.append([in_inglese(res.asse),
                      f"{res.F_collocata:.3e} ({res.F_collocata / res.F_collocata_tot:.1e} of the total)",
                      f"{res.F_non_collocata:.3e} ({res.F_non_collocata / res.F_non_collocata_tot:+.1e})"])
    L += _tab(["axis", "collocated F [rad/(N m)]", "non-collocated F [rad/(N m)]"], righe)
    righe = []
    for res in r.residui:
        for fdt, fz, fm, fc in res.zeri:
            righe.append([in_inglese(res.asse), fdt.replace("motore", "motor"), _g(fz, 6),
                          f"{_g(fm, 6)} ({fm / fz - 1:+.2e})", f"{_g(fc, 6)} ({fc / fz - 1:+.2e})"])
    L += _tab(["axis", "TF", "f full [Hz]", "pure truncation", "truncation + static correction"], righe)
    sp = [f"{in_inglese(a)}: {np.max(np.abs(b.spostamenti)):.1e}" for a, b in r.basi_correzione.items()
          if b.spostamenti is not None and b.spostamenti.size]
    L += ["With `static_correction: residual_vectors` (residual static vectors at the motor and attachment "
          "nodes) the largest shift of the retained frequencies is " + (", ".join(sp) or "zero") + ".", ""]
    L += ["### Mesh check (lumped against consistent mass, retained modes)", ""]
    L += _tab(["axis", "mode", "f consistent [Hz]", "f lumped [Hz]", "deviation", "warning (> 2 %)"],
              [[in_inglese(c.asse), c.modo, _g(c.f_consistente, 6), _g(c.f_concentrata, 6),
                f"{c.scarto * 100:+.2f} %", "**yes**" if c.avviso else "no"] for c in r.maglia])
    return L


# ---------------------------------------------------------------------------
# LaTeX
# ---------------------------------------------------------------------------


def latex(r: RisultatiModali, figure: dict[str, str], segnaposto: bool, *, guyan: bool = True) -> str:
    """Sezione LaTeX; ``guyan = False`` toglie ogni riferimento e confronto con Guyan."""
    from .latex_report import num

    s = r.sistema
    conf_g = r.confronto_guyan if guyan else ()
    pm = s.pm
    nota = rf" \emph{{{NOTA_REALISMO.capitalize()}.}}" if segnaposto else ""
    out = [
        r"\section{Modal truncation (Gawronski): the default reduced model}\label{sec:modal}",
        (r"The default reduction (\texttt{reduction.method: modal}, module \texttt{modal.py}) replaces the "
         r"Guyan chains of Section~\ref{sec:guyan} with a modal truncation in the sense of Gawronski "
         r"(\emph{Modeling and Control of Antennas and Telescopes}, 2008, ch.~2 and~4), built from the "
         r"\emph{full} FEM of each axis (9 azimuth nodes, 10 polar nodes with P). The reason is measured: the "
         r"Guyan masters carry the large inertias, but the flexible modes are formed by light slave parts (the "
         r"discs of azimuth element~2 and polar element~9), so Guyan overestimates their frequencies "
         r"(Table~\ref{tab:validation}); the truncation keeps frequencies, shapes, modal masses and damping of "
         r"the retained modes exactly." if guyan else
         r"Each axis is reduced by a modal truncation in the sense of Gawronski (\emph{Modeling and Control of "
         r"Antennas and Telescopes}, 2008, ch.~2 and~4; \texttt{reduction.method: modal}, module "
         r"\texttt{modal.py}), built from the \emph{full} FEM of each axis (9 azimuth nodes, 10 polar nodes "
         r"with P): the retained modes keep their frequencies, shapes, modal masses and damping exactly, and the "
         r"modes above the band of interest are dropped."),
        r"\paragraph{Modal basis.} For each axis $a$, with the telescope inertia frozen at the attachment node "
        r"in a reference configuration $\theta_{ref}$ and $k_g = 0$,",
        r"\begin{gather}",
        r"\bar{\mathbf M}_\psi = \mathbf M_\psi + \bigl[J^{ax} + J_{\psi,p}(\theta_{ref})\bigr]\mathbf e_L"
        r"\mathbf e_L^{\mathsf T},\qquad \bar{\mathbf M}_\theta = \mathbf M_\theta + I_{JJ}\,\mathbf e_P"
        r"\mathbf e_P^{\mathsf T},\qquad (\mathbf K_a - \omega^2\bar{\mathbf M}_a)\boldsymbol\phi = \mathbf 0,"
        r"\label{eq:modal_eig}\\",
        r"\boldsymbol\Phi_a = [\mathbf 1,\ \boldsymbol\phi_1,\ \dots,\ \boldsymbol\phi_k],\qquad "
        r"\boldsymbol\Phi_a^{\mathsf T}\bar{\mathbf M}_a\boldsymbol\Phi_a = \operatorname{diag}(J_{a,tot}, 1, "
        r"\dots, 1),\qquad \boldsymbol\Phi_a^{\mathsf T}\mathbf K_a\boldsymbol\Phi_a = \operatorname{diag}(0, "
        r"\omega_1^2, \dots, \omega_k^2).",
        r"\end{gather}",
        r"The rigid mode $\mathbf 1$ ($\mathbf K\mathbf 1 = \mathbf 0$, Gawronski eq.~2.30) is not normalised, "
        r"so $\eta_0$ is an angle (the inertia-weighted mean rotation of the chain) with modal mass "
        r"$J_{a,tot} = \mathbf 1^{\mathsf T}\bar{\mathbf M}_a\mathbf 1$; the flexible modes are mass-normalised. "
        r"The rows of $\boldsymbol\Phi$ at the encoder, the motor and the attachment node are Gawronski's "
        r"output and input vectors (eq.~2.21--2.22), $\mathbf c_e$, $\mathbf c_m$, and the new "
        r"$\mathbf c_L$ (azimuth node 9) and $\mathbf c_P$, all with first component 1. With $\mathbf C = "
        r"\beta\mathbf K + \mathbf C_b$: $\zeta_k = \beta\omega_k/2$ exactly, and the bearing friction is kept "
        r"as a full matrix, $\tilde{\mathbf D} = \operatorname{diag}(0, 2\zeta_k\omega_k) + "
        r"\boldsymbol\Phi^{\mathsf T}\mathbf C_b\boldsymbol\Phi$.",
        r"\paragraph{Nonlinear equations in modal coordinates.} The telescope is De Donatis' rigid body attached "
        r"at L and P through $[\psi_L,\ \theta_p]^{\mathsf T} = \mathbf E_\eta\boldsymbol\eta$, "
        r"$\mathbf E_\eta = \operatorname{blkdiag}(\mathbf c_L^{\mathsf T}, \mathbf c_P^{\mathsf T})$:",
        r"\begin{gather}",
        r"\tilde{\mathbf M}\ddot{\boldsymbol\eta} + \tilde{\mathbf D}\dot{\boldsymbol\eta} + \tilde{\mathbf K}"
        r"\boldsymbol\eta + \mathbf E_\eta^{\mathsf T}\Bigl\{\bigl[\mathbf M_p(\theta_p) - \bar{\mathbf M}_p"
        r"\bigr]\ddot{\mathbf q}_p + \mathbf C_p\dot{\mathbf q}_p + \mathbf g_p(\theta_p)\Bigr\} = "
        r"\boldsymbol\Phi^{\mathsf T}\mathbf B\boldsymbol\tau,\label{eq:modal_eom}",
        r"\end{gather}",
        r"with $\bar{\mathbf M}_p = \operatorname{diag}(J_{\psi,p}(\theta_{ref}), I_{JJ})$ the inertia already "
        r"in the modes. $\mathbf M_p$, $\mathbf C_p$, $\mathbf g_p$ come from the Lagrangian of the telescope "
        r"alone, derived by SymPy with the same kinematics as Section~\ref{sec:coupled}, and are compiled with "
        r"\texttt{lambdify}. Keeping all modes, \eqref{eq:modal_eom} is $\boldsymbol\Phi^{\mathsf T}$ times the "
        r"equations of the coupled model with \emph{all} FEM degrees of freedom (the exact reference, 9 + 10 "
        r"coordinates, verified to $10^{-10}$ on random states, test T2)."
        + (r" The same code with the Guyan basis $\mathbf T$ reproduces the SymPy model of "
           r"Section~\ref{sec:coupled} to $10^{-12}$." if guyan else ""),
        r"\paragraph{Linear model.} At each operating point the basis is rebuilt with $\theta_{ref} = \theta_0$, "
        r"and the linearisation is the modal Gawronski model perturbed by three low-rank corrections built from "
        r"the attachment vectors (document \texttt{MODELLO\_unificato}, eq.~7.2):",
        r"\begin{gather}",
        r"\bigl[\tilde{\mathbf M} + M_{12}(\hat{\mathbf c}_L\hat{\mathbf c}_P^{\mathsf T} + \hat{\mathbf c}_P"
        r"\hat{\mathbf c}_L^{\mathsf T})\bigr]\delta\ddot{\boldsymbol\eta} + \bigl[\tilde{\mathbf D} + \gamma"
        r"(\hat{\mathbf c}_L\hat{\mathbf c}_P^{\mathsf T} - \hat{\mathbf c}_P\hat{\mathbf c}_L^{\mathsf T})"
        r"\bigr]\delta\dot{\boldsymbol\eta} + \bigl[\tilde{\mathbf K} + (k_c + k_g)\hat{\mathbf c}_P"
        r"\hat{\mathbf c}_P^{\mathsf T}\bigr]\delta\boldsymbol\eta = \boldsymbol\Phi^{\mathsf T}\mathbf B\,"
        r"\delta\boldsymbol\tau.",
        r"\end{gather}",
        r"$k_c + k_g$ is kept as a rank-1 coupling on $\theta_p$, not inside the modes. In state space, "
        r"$\mathbf x_i = [\omega_i\eta_i,\ \dot\eta_i]$ for the flexible modes and $[\eta_0,\ \dot\eta_0]$ for "
        r"the rigid ones give Gawronski's diagonal blocks $\bigl[\begin{smallmatrix}0 & \omega_i\\ -\omega_i & "
        r"-2\zeta_i\omega_i\end{smallmatrix}\bigr]$ (eq.~2.25) and $\bigl[\begin{smallmatrix}0 & 1\\ 0 & 0"
        r"\end{smallmatrix}\bigr]$ (eq.~2.32), plus off-diagonal terms from the corrections; with all modes the "
        r"FRF equals the full model to $10^{-8}$ (test T3), and with the rigid modes alone it is the "
        r"linearised two-DOF (De Donatis) model (test T4).",
        r"\paragraph{Choice of the modes.} Dynamics above "
        rf"$f_{{band}} = {num(pm.f_banda, 4)}$\,Hz are not of interest. The rigid modes are always kept "
        r"(infinite norm, Gawronski eq.~4.11--4.12); flexible modes are kept if $f_i \le "
        rf"{num(pm.fattore_taglio, 3)}\,f_{{band}}$ and, optionally, if their Hankel norm for the position "
        r"output exceeds a fraction of the largest,",
        r"\begin{gather}",
        r"\|G_i\|_h \simeq \frac{|\phi_i(m)\,\phi_i(e)|}{4\zeta_i\omega_i^2}\qquad\text{(Gawronski "
        r"eq.~4.4--4.6, $\boldsymbol\Phi$ mass-normalised).}",
        r"\end{gather}",
        r"A validation loop then compares the truncated linear model with the full one at every operating point "
        r"(torque$\to$encoder, up to $f_{band}$): direct channels within 1\,dB and $5^\circ$ outside $\pm2\,\%$ "
        r"bands around the reference zeros, each excluded zero within 1\,\% in frequency (near a narrow "
        r"antiresonance the FRF ratio is hypersensitive and is not used alone), and the $2\times2$ matrix "
        r"within 12\,\% in the 2-norm. If a point fails, the discarded mode with the largest Hankel norm is "
        r"added and the loop repeats. "
        + (rf"Here {s.n_aggiunti} mode(s) were added: "
           + ", ".join(rf"{in_inglese(a)} mode {k} ({num(f, 5)}\,Hz)" for a, k, f in s.storia_aggiunte) + ". "
           if s.n_aggiunti else r"No mode had to be added. ")
        + rf"The model has {r.stati} states"
        + (rf", against {r.stati_guyan} of the Guyan chains." if guyan else "."),
    ]
    # tabella dei modi
    out += [r"\begin{table}[H]\centering\small",
            r"\caption{Modes of the full FEM of each axis at $\theta_{ref}$ (consistent mass, telescope frozen, "
            r"$k_g = 0$): frequency, damping, participation at encoder, motor and attachment node, Hankel norm "
            r"(rad/(N\,m)), and whether the mode is kept. Modal mass: $J_{tot}$ for the rigid mode, 1 for the "
            r"flexible ones." + nota + r"}\label{tab:modal_modes}",
            r"\begin{tabular}{llrrrrrrl}\toprule",
            r"axis & mode & $f$ [Hz] & $\zeta$ & $c_e$ & $c_m$ & $c_{L}$/$c_P$ & $\|G_i\|_h$ & kept \\ \midrule"]
    for asse, an in r.analisi.items():
        for k in range(min(an.Phi.shape[1], 6)):
            h = r"$\infty$" if not np.isfinite(an.hankel[k]) else f"${num(an.hankel[k], 3)}$"
            stato = "rigid" if k == 0 else (MOTIVO_TENUTO.get(an.motivo[k], an.motivo[k]) if k in an.tenuti else "no")
            out.append(rf"{in_inglese(asse) if k == 0 else ''} & {k} & ${num(an.f[k], 6)}$ & ${num(an.zeta[k], 3)}$ & "
                       rf"${num(an.c_e[k], 4)}$ & ${num(an.c_m[k], 4)}$ & ${num(an.c_agg[k], 4)}$ & {h} & "
                       rf"{stato} \\")
        out.append(r"\midrule" if asse == "azimut" else "")
    out += [r"\bottomrule\end{tabular}\end{table}"]
    # confronto con il completo e con Guyan
    righe = []
    for e, eg in zip(r.validazione, conf_g or [None] * len(r.validazione)):
        z = max((err for *_, err in e.zeri), default=0.0)
        riga = (rf"${np.rad2deg(e.theta0):g}$ & ${np.rad2deg(e.Omega0):g}$ & ${num(max(e.err_db), 2)}$ & "
                rf"${num(max(e.err_deg), 2)}$ & ${num(z, 2)}$ & ${num(e.err_incrociato, 2)}$")
        if eg is not None:
            zg = max((err for *_, err in eg.zeri), default=0.0)
            riga += (rf" & ${num(max(eg.err_db), 2)}$ & ${num(max(eg.err_deg), 2)}$ & ${num(zg, 2)}$ & "
                     rf"${num(eg.err_incrociato, 2)}$")
        righe.append(riga + r" \\")
    n_col = 10 if conf_g else 6
    out += [r"\begin{table}[H]\centering\small\setlength{\tabcolsep}{4pt}",
            r"\caption{Truncated modal model" + (r" and Guyan chains" if conf_g else "")
            + r" against the full-DOF model, torque$\to$encoder up "
            rf"to {num(pm.f_banda, 4)}\,Hz: worst direct-channel error (dB, degrees) outside the zero bands, "
            r"worst relative error of the excluded zeros, worst relative 2-norm error of the $2\times2$ "
            r"matrix." + nota + r"}\label{tab:modal_vs_full}",
            r"\begin{tabular}{" + "r" * n_col + r"}\toprule",
            r"& & \multicolumn{4}{c}{modal}" + (r" & \multicolumn{4}{c}{Guyan}" if conf_g else "")
            + r" \\ \cmidrule(lr){3-6}" + (r"\cmidrule(lr){7-10}" if conf_g else ""),
            r"$\theta_0$ & $\Omega_0$ & dB & deg & zeros & $2\times2$"
            + (r" & dB & deg & zeros & $2\times2$" if conf_g else "") + r" \\ \midrule",
            *righe,
            r"\bottomrule\end{tabular}\end{table}"]
    # flessibilita' residua e zeri
    out += [r"\paragraph{Residual flexibility of the discarded modes.} A modal truncation keeps the poles of "
            r"the retained modes but moves the zeros, which depend on the static flexibility of the discarded "
            r"ones. For each axis the report gives the collocated and non-collocated residual flexibility "
            r"$\sum_{i\ \text{disc.}}\phi_i(m)^2/\omega_i^2$ and $\sum_{i\ \text{disc.}}\phi_i(e)\phi_i(m)/"
            r"\omega_i^2$ and the zero shift they cause on the isolated axis. The option "
            r"\texttt{static\_correction: residual\_vectors} adds, for every loaded node (motor and attachment), "
            r"the residual static vector $\mathbf r_j = \sum_{i\ \text{disc.}}\boldsymbol\phi_i\phi_i(j)/"
            r"\omega_i^2$ (residual flexibility with inertia relief, the only one defined for a free chain), "
            r"orthogonalises it and solves the eigenproblem in the enlarged space: the model stays a valid ODE "
            r"also in the nonlinear case, and the retained frequencies do not move.",
            r"\begin{table}[H]\centering\small",
            r"\caption{Zeros of the isolated axes (Hz): full FEM, pure truncation, truncation plus static "
            r"correction; residual flexibilities in rad/(N\,m)." + nota + r"}\label{tab:modal_zeros}",
            r"\begin{tabular}{llrrr}\toprule",
            r"axis & transfer function & full & truncated & + static corr. \\ \midrule"]
    for res in r.residui:
        for fdt, fz, fm, fc in res.zeri:
            fdt_tex = fdt.replace("tau -> ", r"$\tau\to$ ").replace("motore", "motor")
            out.append(rf"{in_inglese(res.asse)} & {fdt_tex} & ${num(fz, 6)}$ & ${num(fm, 6)}$ & ${num(fc, 6)}$ \\")
        out.append(rf"\multicolumn{{5}}{{r}}{{\footnotesize $F_{{coll}} = {num(res.F_collocata, 3)}$, "
                   rf"$F_{{non\text{{-}}coll}} = {num(res.F_non_collocata, 3)}$}} \\ \midrule")
    out += [r"\bottomrule\end{tabular}\end{table}"]
    maglia = [c for c in r.maglia]
    out += [r"\paragraph{Mesh check.} For the retained modes the lumped and consistent mass give "
            + ", ".join(rf"{in_inglese(c.asse)} {c.modo}: ${c.scarto * 100:+.2f}\,\%$" for c in maglia)
            + r" (a warning is raised above 2\,\%)."]
    if "bode_modale" in figure:
        out += [r"\begin{figure}[H]\centering",
                rf"\includegraphics[width=\linewidth,height=0.8\textheight,keepaspectratio]{{{figure['bode_modale']}}}",
                r"\caption{Torque$\to$encoder $2\times2$ matrix up to the band of interest: truncated modal model "
                r"(with its poles, dotted, and invariant zeros, circles) and full-DOF reference"
                + (r"; Guyan chains for comparison." if guyan else ".")
                + nota + r"}\label{fig:bode_modal}\end{figure}"]
    return "\n".join(out)
