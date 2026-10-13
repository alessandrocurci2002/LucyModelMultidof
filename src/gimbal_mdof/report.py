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
from .linearize import ModelloLineare
from .params import MM, NOTA_REALISMO, Avviso, Parametri
from .simulate import RisultatoScenario

NOTA_REALISMO_MD = f"""\
> **Nota di realismo (sezione 0 della specifica).** Geometrie degli alberi,
> materiali, rigidezze di giunzione, inerzie di rotori e forcella e smorzamenti
> sono oggi **segnaposto**. Ogni frequenza, FdT e figura di questo report e' un
> **{NOTA_REALISMO}** e non va presentata come caratterizzazione del sistema reale.
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


def _tab(intest: Sequence[str], righe: Sequence[Sequence[object]]) -> list[str]:
    out = ["| " + " | ".join(intest) + " |", "|" + "|".join("---" for _ in intest) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in righe]
    return out + [""]


def _g(x: float, c: int = 4) -> str:
    return "--" if x is None or not np.isfinite(x) else f"{x:.{c}g}"


STATO_DOMANDE: dict[str, str] = {
    "D1": "Confermata: telescopio sull'elemento polare 3, nodo P nella mezzeria.",
    "D2": "Confermata: elementi polari 2 e 4 sui cuscinetti.",
    "D3": "Confermata: master polari [1, P, 7, 11].",
    "Q4": "Risposta dell'utente (2026-10-12): Guyan sui master, Craig-Bampton solo come verifica.",
    "Q6": "Risposta dell'utente (2026-10-12): valori del file .m (acciaio 8000 kg/m^3, 77 GPa; "
          "alluminio 7075 2810 kg/m^3, 26.9 GPa). Restano segnaposto per pezzo.",
    "Q9": "Risposta dell'utente (2026-10-12): input manuale della sola forcella da Onshape. "
          "Valore attuale segnaposto.",
    "Q10": "Risolta: z_O = 54.9 cm misurato nel CAD e gia' confermato nel modello a 2 GdL.",
    "Q12": "Ipotesi adottata e dichiarata: forcella rigida.",
}
DEFAULT_DOMANDE: dict[str, str] = {
    "Q1": "Default: xi_p = 0.5 (telescopio su tutto l'elemento 3).",
    "Q2": "Default: nessun irrigidimento del morsetto (payload_clamp_rigid_mm = 0).",
    "Q5": "Default: rigidezze dalla geometria, nessuna cedevolezza di giunzione.",
    "Q7": "Default: null, J^ax senza Steiner (ininfluente sulla somma J^ax + J_psi,p).",
    "Q8": "Default: C = beta K con zeta = 0.01 sul primo modo flessibile; cuscinetti b = 0.",
    "Q11": "Default: Omega_max = 10 deg/s.",
}


def genera(c: ContenutoReport, percorso: str | Path) -> Path:
    par = c.par
    s = c.sistemi[("consistent", "guyan")]
    L: list[str] = ["# Modello multi-GdL del gimbal alt-azimutale: report", ""]
    if par.segnaposto:
        L += [NOTA_REALISMO_MD, ""]

    # --- parametri ------------------------------------------------------------
    tr = par.terne
    L += ["## Sistemi di riferimento (convenzione De Donatis)", "",
          "- N inerziale: origine O_N (origine dell'assieme Onshape, base dell'albero di azimut), "
          "versori x_hat, y_hat, z_hat;",
          f"- corpo 1 (azimut): origine o_1 sull'asse di azimut, dist(O_N, o_1) = {tr.z_o1 * 1e3:g} mm "
          "lungo z_hat; versori i_hat, j_hat, k_hat = z_hat;",
          f"- corpo 2 (polare): origine o_2 = centro di rotazione dell'albero polare, dist(O_N, o_2) = "
          f"{tr.z_o2 * 1e3:g} mm lungo z_hat; versori K_hat (puntamento, verso le ottiche), "
          "J_hat = j_hat (asse polare), I_hat = J_hat x K_hat;",
          "- le origini dei corpi stanno sugli assi di rotazione, non nei baricentri.", ""]
    L += ["## Parametri usati", "", "### Materiali", ""]
    L += _tab(["materiale", "rho [kg/m^3]", "G [GPa]", "segnaposto"],
              [[m.nome, _g(m.rho), _g(m.G / 1e9), "**si**" if m.segnaposto else "no"]
               for m in par.materiali.values()])
    for nome, mod, a in (("azimut", s.asse_psi, par.azimut), ("polare", s.asse_theta, par.polare)):
        L += [f"### Asse {nome}", ""]
        L += _tab(["elemento", "modo", "L [mm]", "J [kg m^2]", "k [N m/rad]", "m [kg]", "I_t [kg m^2]",
                   "segnaposto"],
                  [[p.etichetta, p.modo, _g(p.lunghezza / MM), _g(p.J), _g(p.k), _g(p.m), _g(p.I_t),
                    "**si**" if p.segnaposto else "no"] for p in mod.elementi_originali])
        L += _tab(["nodo", "inerzia", "J [kg m^2]", "fonte", "segnaposto"],
                  [[n.nodo, n.nome, _g(n.J), n.fonte, "**si**" if n.segnaposto else "no"]
                   for n in a.inerzie_nodali])
        L += [f"Scomposizione per strato: `results/tables/element_breakdown_{nome}.csv`.", ""]
    pl = s.payload
    pp = par.payload
    origine = ("corpo B del CAD del modello a 2 GdL" if pp.sorgente == "gimbal_2dof" else
               f"config.yaml, frame `{pp.terna}`" + (f" ({pp.nota})" if pp.nota else "")
               + (" **[segnaposto]**" if pp.segnaposto else ""))
    L += ["### Telescopio (payload, Q3)", "",
          f"- origine: {origine}; sottrazione delle parti polari: "
          f"{'si' if pp.sottrai_alberi else 'no'}",
          "- valori risultanti, rispetto a o_2 nella terna del corpo 2 (I_hat, J_hat, K_hat):",
          f"  - massa m_p = {pl.m:.4f} kg, eccentricita' d_2 = (d_I, d_J, d_K) = ({pl.d_2[0] * 1e3:+.2f}, "
          f"{pl.d_2[1] * 1e3:+.2f}, {pl.d_2[2] * 1e3:+.2f}) mm (K_hat = puntamento)",
          f"  - I_II = {pl.I_2[0, 0]:.6f}, I_JJ = {pl.I_2[1, 1]:.6f}, I_KK = {pl.I_2[2, 2]:.6f} kg m^2; "
          f"I_IJ = {pl.I_2[0, 1]:.3e}, I_IK = {pl.I_2[0, 2]:.3e}, I_JK = {pl.I_2[1, 2]:.3e} kg m^2",
          f"- sottratti: massa {pl.sottratto.get('m', 0):.4f} kg, sum J = "
          f"{pl.sottratto.get('J', 0):.3e} kg m^2",
          f"- J^ax = {pl.J_ax:.4e} kg m^2", ""]

    # --- riduzione -------------------------------------------------------------
    L += ["## Riduzione", "",
          f"Rigidezze della catena ridotta: k_psi = {np.array2string(s.rid_psi.k_catena, precision=5)} N m/rad, "
          f"k_theta = {np.array2string(s.rid_theta.k_catena, precision=5)} N m/rad.",
          f"Smorzamento: beta_psi = {s.beta_psi:.4e} s, beta_theta = {s.beta_theta:.4e} s.", "",
          "### Validazione (assi isolati, theta0 = 0)", ""]
    righe = []
    for cf in c.confronti:
        for k in range(len(cf.f_red)):
            righe.append([cf.asse, cf.modello_massa, cf.forma_massa, k + 1, _g(cf.f_full[k], 6),
                          _g(cf.f_red[k], 6), f"{cf.errore_relativo[k] * 100:+.3f} %",
                          _g(cf.f_limite, 4) if k == 0 else ""])
    L += _tab(["asse", "massa", "massa ridotta", "modo", "f completo [Hz]", "f ridotto [Hz]", "errore",
               "f limite 3 dB [Hz]"], righe)

    # --- punti di lavoro ---------------------------------------------------------
    L += ["## Punti di lavoro, frequenze e zeri", ""]
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
    L += _tab(["massa", "forma", "theta0 [deg]", "Omega0 [deg/s]", "tau_theta0 [N m]",
               "prime frequenze [Hz]", "polo instabile"], righe)
    ml0 = c.lineari.get(("consistent", "guyan", par.theta0[0], 0.0))
    if ml0 is not None:
        L += [f"Zeri invarianti della matrice coppia -> encoder (consistent/guyan, theta0 = "
              f"{np.rad2deg(ml0.theta0):g} deg, Omega0 = 0), frequenze in Hz:", ""]
        righe = []
        for u in ("psi_e", "theta_e"):
            for j, ing in enumerate(("tau_psi", "tau_theta")):
                z = tf.zeri_invarianti(ml0.A, ml0.Bu[:, j], ml0.Cy[ml0.riga(u)])
                fz = sorted({round(abs(v) / (2 * np.pi), 4) for v in z if np.imag(v) >= 0})
                righe.append([f"{u}/{ing}", ", ".join(_g(v, 5) for v in fz[:6]) or "--"])
        L += _tab(["FdT", "zeri"], righe)
    L += ["Tabelle complete: `results/tables/poles_zeros.csv`; matrici di stato in `results/tf/`.", ""]

    # --- simulazioni ---------------------------------------------------------------
    L += ["## Simulazioni non lineari", ""]
    for r in c.simulazioni:
        if r.nome == "S3":
            L.append(f"- **S3** ({r.descrizione}): deriva relativa di E = {r.derive['energia']:.2e}, "
                     f"di p_psi = {r.derive['momento_psi']:.2e} (soglia 1e-8).")
        else:
            err = ", ".join(f"{k}: {v:.2e}" for k, v in r.errori.items() if np.isfinite(v))
            picco_th = np.max(np.abs(r.y_nl[:, 1])) * 206264.8
            L.append(f"- **{r.nome}** ({r.descrizione}, theta0 = {np.rad2deg(r.theta0):g} deg, "
                     f"Omega0 = {np.rad2deg(r.Omega0):g} deg/s): errore relativo non lineare/lineare "
                     f"{err}. Picco di theta_e: {picco_th:.3g} arcsec.")
    L += ["", "Nei canali polari di S1 la risposta lineare esiste solo per i prodotti d'inerzia: il "
          "termine quadratico (1/2) J' dpsi^2 vi e' dello stesso ordine, quindi l'errore relativo non e' "
          "un indice della qualita' della linearizzazione.", ""]

    # --- figure -----------------------------------------------------------------
    L += ["## Figure", ""]
    L += [f"- `{p.as_posix()}`" for p in c.figure] + [""]
    L += ["## Tabelle", ""] + [f"- `{p.as_posix()}`" for p in c.tabelle] + [""]
    if c.latex:
        L += ["## Documento LaTeX", "", c.latex, ""]

    # --- domande aperte -------------------------------------------------------
    L += ["## Stato delle decisioni e delle domande aperte (sezione 3)", ""]
    righe = [[k, v] for k, v in STATO_DOMANDE.items()] + [[k, f"**{v}** (ipotesi)"]
                                                          for k, v in DEFAULT_DOMANDE.items()]
    if par.payload.sorgente == "manual":
        q3 = (f"Decisione dell'utente (2026-10-12): tensore del telescopio inserito in config.yaml "
              f"(frame `{par.payload.terna}`), sottrazione delle parti polari: "
              f"{'si' if par.payload.sottrai_alberi else 'no'}"
              + (". Valori attuali segnaposto." if par.payload.segnaposto else "."))
    else:
        q3 = "**Default: tensore del corpo B del 2 GdL meno gli elementi polari del FEM.** (ipotesi)"
    righe.append(["Q3", q3])
    righe.sort(key=lambda r: (r[0][0], int(r[0][1:])))
    L += _tab(["ID", "stato"], righe)

    # --- avvisi -----------------------------------------------------------------
    L += ["## Avvisi", ""]
    etichette = {"segnaposto": "SEGNAPOSTO", "incoerenza": "INCOERENZA", "informazione": "NOTA"}
    L += [f"- **[{etichette[a.severita]}] {a.codice}**: {a.messaggio}" for a in c.avvisi]
    L += [""] + [f"- {n}" for n in c.note]

    percorso = Path(percorso)
    percorso.parent.mkdir(parents=True, exist_ok=True)
    percorso.write_text("\n".join(L) + "\n", encoding="utf-8")
    return percorso
