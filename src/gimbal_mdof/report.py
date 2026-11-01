"""Generazione di ``results/report.md`` (sezione 14).

Contenuti: nota di realismo, parametri con i segnaposto evidenziati,
scomposizione per elemento, validazione della riduzione, frequenze e zeri,
simulazioni, elenco delle figure, stato delle domande aperte della sezione 3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from . import tf
from .analysis import ConfrontoRiduzione, frequenze_modi
from .coupled import Sistema
from .linearize import USCITE, ModelloLineare
from .params import MM, NOTA_REALISMO, Avviso, Parametri, in_inglese
from .simulate import RisultatoScenario

NOTA_REALISMO_MD = f"""\
> **Realism note (section 0 of the specification).** Shaft geometries,
> materials, joint stiffnesses, rotor and fork inertias and damping values
> are currently **placeholders**. Every frequency, transfer function and figure of
> this report is a **{NOTA_REALISMO}** and must not be presented as a
> characterisation of the real system.
"""


@dataclass
class ContenutoReport:
    par: Parametri
    sistemi: Mapping[tuple[str, str], Sistema]
    confronti: Sequence[ConfrontoRiduzione]
    lineari: Mapping[tuple[str, str, float, float], ModelloLineare]
    simulazioni: Sequence[RisultatoScenario]
    figure: Sequence[Path]
    tabelle: Sequence[Path]
    avvisi: Sequence[Avviso]
    latex: str = ""
    note: list[str] = field(default_factory=list)
    riduzione: str = "guyan"
    """Riduzione del modello principale (lineari e simulazioni): modal | guyan | craig_bampton."""
    modale: object | None = None
    """:class:`gimbal_mdof.modal_output.RisultatiModali`, se calcolati."""
    lineari_principali: Mapping[tuple[float, float], ModelloLineare] = field(default_factory=dict)
    """Modelli lineari della riduzione principale, per (theta0, Omega0)."""
    s1_riferimento: RisultatoScenario | None = None
    """S1 sul modello a tutti i GdL del FEM, per il confronto con il modello principale."""
    extra_md: Sequence[str] = ()
    """Sezioni aggiuntive gia' in Markdown (modello a parametri concentrati)."""


def _tab(intest: Sequence[str], righe: Sequence[Sequence[object]]) -> list[str]:
    out = ["| " + " | ".join(intest) + " |", "|" + "|".join("---" for _ in intest) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in righe]
    return out + [""]


def _g(x: float, c: int = 4) -> str:
    return "--" if x is None or not np.isfinite(x) else f"{x:.{c}g}"


STATO_DOMANDE: dict[str, str] = {
    "D1": "Confirmed: telescope on polar element 3, node P at the midpoint.",
    "D2": "Confirmed: polar elements 2 and 4 on the bearings.",
    "D3": "Confirmed: polar masters [1, P, 7, 11].",
    "Q4": "User decision (2026-10-13): Gawronski modal truncation of the full FEM "
          "(reduction.method: modal); Guyan (decision of 2026-10-12) and Craig-Bampton remain as comparisons.",
    "Q6": "User answer (2026-10-12): values of the .m file (steel 8000 kg/m^3, 77 GPa; "
          "aluminium 7075 2810 kg/m^3, 26.9 GPa). Still placeholders per part.",
    "Q10": "Solved: z_O = 54.9 cm measured in the CAD and already confirmed in the 2-DOF model.",
    "Q12": "Assumption adopted and declared: rigid fork.",
}
DEFAULT_DOMANDE: dict[str, str] = {
    "Q1": "Default: xi_p = 0.5 (telescope over the whole element 3).",
    "Q2": "Default: no stiffening by the clamp (payload_clamp_rigid_mm = 0).",
    "Q5": "Default: stiffnesses from the geometry, no joint compliance.",
    "Q7": "Default: null, J^ax without Steiner terms (irrelevant if the azimuth load includes the polar shafts).",
    "Q8": "Default: C = beta K with zeta = 0.01 on the first flexible mode; bearings b = 0.",
    "Q11": "Default: Omega_max = 10 deg/s.",
}


def _si(v: bool) -> str:
    return "**yes**" if v else "no"


def genera(c: ContenutoReport, percorso: str | Path) -> Path:
    par = c.par
    s = c.sistemi[("consistent", "guyan")]
    L: list[str] = ["# Multi-DOF model of the alt-azimuth gimbal: report", ""]
    if par.segnaposto:
        L += [NOTA_REALISMO_MD, ""]

    # --- parametri ------------------------------------------------------------
    tr = par.terne
    L += ["## Reference frames (De Donatis convention)", "",
          "- N inertial: origin O_N (origin of the Onshape assembly, base of the azimuth shaft), "
          "unit vectors x_hat, y_hat, z_hat;",
          f"- body 1 (azimuth): origin o_1 on the azimuth axis, dist(O_N, o_1) = {tr.z_o1 * 1e3:g} mm "
          "along z_hat; unit vectors i_hat, j_hat, k_hat = z_hat;",
          f"- body 2 (polar): origin o_2 = centre of rotation of the polar shaft, dist(O_N, o_2) = "
          f"{tr.z_o2 * 1e3:g} mm along z_hat; unit vectors K_hat (pointing, towards the optics), "
          "J_hat = j_hat (polar axis), I_hat = J_hat x K_hat;",
          "- the origins of the bodies lie on the rotation axes, not at the centres of mass.", ""]
    L += ["## Parameters used", "", "### Materials", ""]
    L += _tab(["material", "rho [kg/m^3]", "G [GPa]", "placeholder"],
              [[in_inglese(m.nome), _g(m.rho), _g(m.G / 1e9), _si(m.segnaposto)]
               for m in par.materiali.values()])
    for nome, mod, a in (("azimut", s.asse_psi, par.azimut), ("polare", s.asse_theta, par.polare)):
        L += [f"### {in_inglese(nome).capitalize()} axis", ""]
        L += _tab(["element", "mode", "L [mm]", "J [kg m^2]", "k [N m/rad]", "m [kg]", "I_t [kg m^2]",
                   "placeholder"],
                  [[p.etichetta, p.modo, _g(p.lunghezza / MM), _g(p.J), _g(p.k), _g(p.m), _g(p.I_t),
                    _si(p.segnaposto)] for p in mod.elementi_originali])
        L += _tab(["node", "inertia", "J [kg m^2]", "source", "placeholder"],
                  [[n.nodo, n.nome, _g(n.J), n.fonte, _si(n.segnaposto)]
                   for n in a.inerzie_nodali])
        L += [f"Breakdown by layer: `results/tables/element_breakdown_{nome}.csv`.", ""]
    pl = s.payload
    pp = par.payload
    origine = ("body B of the CAD of the 2-DOF model" if pp.sorgente == "gimbal_2dof" else
               f"config.yaml, frame `{pp.terna}`" + (f" ({pp.nota})" if pp.nota else "")
               + (" **[placeholder]**" if pp.segnaposto else ""))
    L += ["### Telescope (payload, Q3)", "",
          f"- source: {origine}; subtraction of the polar parts: "
          f"{'yes' if pp.sottrai_alberi else 'no'}",
          "- resulting values, about o_2 in the frame of body 2 (I_hat, J_hat, K_hat):",
          f"  - mass m_p = {pl.m:.4f} kg, eccentricity d_2 = (d_I, d_J, d_K) = ({pl.d_2[0] * 1e3:+.2f}, "
          f"{pl.d_2[1] * 1e3:+.2f}, {pl.d_2[2] * 1e3:+.2f}) mm (K_hat = pointing)",
          f"  - I_II = {pl.I_2[0, 0]:.6f}, I_JJ = {pl.I_2[1, 1]:.6f}, I_KK = {pl.I_2[2, 2]:.6f} kg m^2; "
          f"I_IJ = {pl.I_2[0, 1]:.3e}, I_IK = {pl.I_2[0, 2]:.3e}, I_JK = {pl.I_2[1, 2]:.3e} kg m^2",
          f"- subtracted: mass {pl.sottratto.get('m', 0):.4f} kg, sum J = "
          f"{pl.sottratto.get('J', 0):.3e} kg m^2",
          f"- J^ax = {pl.J_ax:.4e} kg m^2"
          + (f" (J^ax = {pl.J_ax_nel_carico:.4e} kg m^2 from the polar elements is already in the azimuth "
             "load: not added)" if pl.J_ax_nel_carico else ""), ""]

    # --- riduzione -------------------------------------------------------------
    L += ["## Reduction", "",
          f"Stiffnesses of the reduced chain: k_psi = {np.array2string(s.rid_psi.k_catena, precision=5)} N m/rad, "
          f"k_theta = {np.array2string(s.rid_theta.k_catena, precision=5)} N m/rad.",
          f"Damping: beta_psi = {s.beta_psi:.4e} s, beta_theta = {s.beta_theta:.4e} s.", "",
          "### Validation (isolated axes, theta0 = 0)", ""]
    righe = []
    for cf in c.confronti:
        for k in range(len(cf.f_red)):
            righe.append([in_inglese(cf.asse), cf.modello_massa, cf.forma_massa, k + 1, _g(cf.f_full[k], 6),
                          _g(cf.f_red[k], 6), f"{cf.errore_relativo[k] * 100:+.3f} %",
                          _g(cf.f_limite, 4) if k == 0 else ""])
    L += _tab(["axis", "mass", "reduced mass", "mode", "f full [Hz]", "f reduced [Hz]", "error",
               "3 dB limit f [Hz]"], righe)

    # --- punti di lavoro ---------------------------------------------------------
    L += ["## Operating points, frequencies and zeros", ""]
    righe = []
    for (mm, mf, th, Om), ml in c.lineari.items():
        f = frequenze_modi(ml)
        # il doppio polo nullo dell'azimut (difettivo) e' sparpagliato da eig a ~1e-9 |p|_max:
        # non e' un'instabilita'
        soglia = 1e-7 * float(np.max(np.abs(ml.poli)))
        instabili = [p for p in ml.poli if np.real(p) > soglia and abs(np.imag(p)) < soglia]
        righe.append([mm, mf, f"{np.rad2deg(th):g}", f"{np.rad2deg(Om):g}",
                      _g(ml.punto.tau0[1], 5), ", ".join(_g(v, 5) for v in f[:5]),
                      f"{np.real(instabili[0]):.4g} rad/s" if instabili else "no"])
    L += _tab(["mass", "form", "theta0 [deg]", "Omega0 [deg/s]", "tau_theta0 [N m]",
               "first frequencies [Hz]", "unstable pole"], righe)
    ml0 = c.lineari.get(("consistent", "guyan", par.theta0[0], 0.0))
    if ml0 is not None:
        L += [f"Invariant zeros of the torque -> encoder matrix (consistent/guyan, theta0 = "
              f"{np.rad2deg(ml0.theta0):g} deg, Omega0 = 0), frequencies in Hz:", ""]
        righe = []
        for u in ("psi_e", "theta_e"):
            for j, ing in enumerate(("tau_psi", "tau_theta")):
                z = tf.zeri_invarianti(ml0.A, ml0.Bu[:, j], ml0.Cy[ml0.riga(u)])
                fz = sorted({round(abs(v) / (2 * np.pi), 4) for v in z if np.imag(v) >= 0})
                righe.append([f"{u}/{ing}", ", ".join(_g(v, 5) for v in fz[:6]) or "--"])
        L += _tab(["TF", "zeros"], righe)
    L += ["Full tables: `results/tables/poles_zeros.csv`; state matrices in `results/tf/`.", ""]

    # --- troncamento modale -----------------------------------------------------------
    if c.modale is not None:
        from . import modal_output as mo

        L += mo.markdown(c.modale, par.segnaposto)  # type: ignore[arg-type]
    if c.lineari_principali:
        L += [f"### Main model ({c.riduzione}): operating points", ""]
        righe = []
        for (th, Om), ml in c.lineari_principali.items():
            f = frequenze_modi(ml)
            righe.append([f"{np.rad2deg(th):g}", f"{np.rad2deg(Om):g}", _g(ml.punto.tau0[1], 5),
                          ", ".join(_g(v, 5) for v in f[:6])])
        L += _tab(["theta0 [deg]", "Omega0 [deg/s]", "tau_theta0 [N m]", "first frequencies [Hz]"], righe)

    # --- simulazioni ---------------------------------------------------------------
    L += ["## Nonlinear simulations", "",
          f"Model: {c.riduzione}" + (" (modal truncation, theta_ref from config.yaml)"
                                     if c.riduzione == "modal" else "") + ".", ""]
    for r in c.simulazioni:
        if r.nome == "S3":
            L.append(f"- **S3** ({r.descrizione}): relative drift of E = {r.derive['energia']:.2e}, "
                     f"of p_psi = {r.derive['momento_psi']:.2e} (threshold 1e-8).")
        else:
            err = ", ".join(f"{k}: {v:.2e}" for k, v in r.errori.items() if np.isfinite(v))
            picco_th = np.max(np.abs(r.y_nl[:, 1])) * 206264.8
            L.append(f"- **{r.nome}** ({r.descrizione}, theta0 = {np.rad2deg(r.theta0):g} deg, "
                     f"Omega0 = {np.rad2deg(r.Omega0):g} deg/s): relative nonlinear/linear error "
                     f"{err}. Peak of theta_e: {picco_th:.3g} arcsec.")
    if c.s1_riferimento is not None and c.simulazioni:
        rf_, rm_ = c.s1_riferimento, c.simulazioni[0]
        n_t = min(rf_.t.size, rm_.t.size)
        diff = []
        for u in ("psi_e", "psi_m", "psi_L"):
            i = USCITE.index(u)
            picco = float(np.max(np.abs(rf_.y_nl[:n_t, i])))
            diff.append(f"{u}: {np.max(np.abs(rm_.y_nl[:n_t, i] - rf_.y_nl[:n_t, i])) / picco:.2e}")
        L.append(f"- **S1, comparison with the all-DOF model** (DOP853 with rtol 1e-8, atol 1e-10: the modes "
                 f"of the full FEM up to hundreds of kHz make the explicit integrator stability-limited; "
                 f"t_end = {rf_.t[-1]:g} s): maximum deviation relative to the peak " + ", ".join(diff) + ".")
    L += ["", "In the polar channels of S1 the linear response exists only through the products of inertia: "
          "the quadratic term (1/2) J' dpsi^2 is of the same order there, so the relative error is not "
          "a measure of the quality of the linearisation.", ""]

    L += list(c.extra_md)

    # --- figure -----------------------------------------------------------------
    L += ["## Figures", ""]
    L += [f"- `{p.as_posix()}`" for p in c.figure] + [""]
    L += ["## Tables", ""] + [f"- `{p.as_posix()}`" for p in c.tabelle] + [""]
    if c.latex:
        L += ["## LaTeX document", "", c.latex, ""]

    # --- domande aperte -------------------------------------------------------
    L += ["## Status of the decisions and open questions (section 3)", ""]
    righe = [[k, v] for k, v in STATO_DOMANDE.items()] + [[k, f"**{v}** (assumption)"]
                                                          for k, v in DEFAULT_DOMANDE.items()]
    if par.payload.sorgente == "manual":
        q3 = (f"User decision (2026-10-12): telescope tensor entered in config.yaml "
              f"(frame `{par.payload.terna}`), subtraction of the polar parts: "
              f"{'yes' if par.payload.sottrai_alberi else 'no'}"
              + (". Current values are placeholders." if par.payload.segnaposto else "."))
    else:
        q3 = "**Default: tensor of body B of the 2-DOF model minus the polar elements of the FEM.** (assumption)"
    righe.append(["Q3", q3])
    carico = [n for n in par.azimut.inerzie_nodali if n.nodo == par.azimut.nodo_estremo]
    if carico and carico[0].J_G is not None:
        n9 = carico[0]
        q9 = (f"User data (2026-10-13): load '{n9.nome}' from Onshape ({n9.fonte}), m = {n9.m:.3f} kg, "
              f"J = I_zz,G + m d^2 = {n9.J_G:.4f} + {n9.m:.3f} x {n9.distanza_asse * 1e3:.3f}^2 mm^2 = "
              f"{n9.J:.5f} kg m^2"
              + ("; it includes the polar shafts: J^ax not added (already in the CAD value)."
                 if n9.contiene_alberi_polari else "."))
    else:
        q9 = ("User answer (2026-10-12): manual input of the fork alone from Onshape"
              + (". Current value is a placeholder." if carico and carico[0].segnaposto else "."))
    righe.append(["Q9", q9])
    righe.sort(key=lambda r: (r[0][0], int(r[0][1:])))
    L += _tab(["ID", "status"], righe)

    # --- avvisi -----------------------------------------------------------------
    L += ["## Warnings", ""]
    etichette = {"segnaposto": "PLACEHOLDER", "incoerenza": "INCONSISTENCY", "informazione": "NOTE"}
    L += [f"- **[{etichette[a.severita]}] {a.codice}**: {a.messaggio}" for a in c.avvisi]
    L += [""] + [f"- {n}" for n in c.note]

    percorso = Path(percorso)
    percorso.parent.mkdir(parents=True, exist_ok=True)
    percorso.write_text("\n".join(L) + "\n", encoding="utf-8")
    return percorso
