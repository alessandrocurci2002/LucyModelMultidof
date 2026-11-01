"""Documento LaTeX dell'analisi FEM, generato e compilato (sezione 12).

Genera ``results/latex/multidof_fem_model.tex`` (in inglese, coerente con
``gimbal_2dof/modello_matematico.tex``) e lo compila con ``latexmk -pdf`` se
disponibile, altrimenti con ``pdflatex`` tre volte. La geometria inserita in
``config.yaml`` e' disegnata in TikZ, una figura per asse.

Note tecniche obbligatorie della specifica, rispettate qui:
``\\setcounter{MaxMatrixCols}{20}``; forzanti scritte come
:math:`\\tau\\mathbf e_p`; blocchi di equazioni in ``gather``;
``\\texorpdfstring`` per la matematica nei titoli; log controllato.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import sympy as sp

from . import coupled as cp
from . import fem
from . import tf
from .analysis import ConfrontoRiduzione
from .coupled import Sistema
from .geometry import ProprietaElemento
from .linearize import ModelloLineare
from .params import MM, NOTA_REALISMO, R_P2, Avviso, ParametriAsse, Parametri, in_inglese

NOME_FILE = "multidof_fem_model"


# Notazione del documento (terne della tesi De Donatis). Le equazioni sono derivate nella
# terna interna di puntamento P (x_P = K_hat, y_P = -J_hat, z_P = I_hat); per la stampa le
# componenti si riscrivono nella terna 2 con i segni dati da R_P2.
_I_II, _I_JJ, _I_KK = sp.symbols("I_{II} I_{JJ} I_{KK}", positive=True)
_I_IJ, _I_IK, _I_JK = sp.symbols("I_{IJ} I_{IK} I_{JK}", real=True)
_d_I, _d_J, _d_K = sp.symbols("d_{I} d_{J} d_{K}", real=True)
SOSTITUZIONE_TERNA_2 = {
    cp.I_xx: _I_KK, cp.I_yy: _I_JJ, cp.I_zz: _I_II,
    cp.I_xy: -_I_JK, cp.I_xz: _I_IK, cp.I_yz: -_I_IJ,
    cp.d_x: _d_K, cp.d_y: -_d_J, cp.d_z: _d_I,
}


def L(espr) -> str:
    """``sympy.latex`` con le componenti nella terna 2 del documento."""
    if isinstance(espr, sp.MatrixBase):  # le matrici mutabili non sono sp.Basic
        espr = espr.applyfunc(lambda e: e.xreplace(SOSTITUZIONE_TERNA_2))
    elif isinstance(espr, sp.Basic):
        espr = espr.xreplace(SOSTITUZIONE_TERNA_2)
    return sp.latex(espr)


# ---------------------------------------------------------------------------
# Utilita'
# ---------------------------------------------------------------------------


def esc(testo: str) -> str:
    """Escape dei caratteri speciali di LaTeX nel testo semplice."""
    sost = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
            "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\^{}"}
    return "".join(sost.get(c, c) for c in str(testo))


def num(x: float, cifre: int = 4) -> str:
    """Numero in notazione LaTeX (mantissa per 10 alla potenza)."""
    if x is None or not np.isfinite(x):
        return "--"
    if x == 0:
        return "0"
    s = f"{x:.{cifre}g}"
    if "e" in s:
        m, e = s.split("e")
        return rf"{m}\times10^{{{int(e)}}}"
    return s


def matrice_numerica(A: np.ndarray, cifre: int = 3) -> str:
    righe = [" & ".join(num(v, cifre) if abs(v) > 1e-14 * np.max(np.abs(A)) else "0" for v in r)
             for r in A]
    return r"\begin{bmatrix}" + r" \\ ".join(righe) + r"\end{bmatrix}"


def riduci_se_serve(eq: str) -> str:
    """Equazione centrata, ridotta alla larghezza del testo solo se la supera."""
    return (r"\begin{center}\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{$\displaystyle "
            + eq + r"$}\end{center}")


def ridimensiona(eq: str) -> str:
    """Equazione centrata ridimensionata alla larghezza del testo."""
    return (r"\begin{center}\resizebox{\linewidth}{!}{$\displaystyle " + eq
            + r"$}\end{center}")


# ---------------------------------------------------------------------------
# Geometria in TikZ
# ---------------------------------------------------------------------------

LARGHEZZA_CM = 15.5
LARGHEZZA_MIN_CM = 1.05
ALTEZZA_MAX_CM = 3.0


def _scala_orizzontale(lunghezze: Sequence[float]) -> float:
    """Scala s (cm/mm) tale che sum(max(L s, w_min)) = larghezza totale."""
    lo, hi = 1e-6, 10.0
    for _ in range(100):
        s = 0.5 * (lo + hi)
        tot = sum(max(L * s, LARGHEZZA_MIN_CM) for L in lunghezze)
        lo, hi = (s, hi) if tot < LARGHEZZA_CM else (lo, s)
    return 0.5 * (lo + hi)


def _colore(strato) -> str:
    if not strato.portante:
        return "orange!45"
    return "blue!30" if strato.materiale.nome == "acciaio" else "cyan!20"


def tikz_geometria(a: ParametriAsse, mod: fem.ModelloAsse) -> tuple[str, bool]:
    """Vista laterale stilizzata di un asse. Restituisce ``(tikz, scala_per_elemento)``.

    Ogni segmento e' un rettangolo largo in proporzione a L e alto in
    proporzione a :math:`D_o`; i fori dei cilindri cavi sono tratteggiati, gli
    strati non portanti in arancio. Gli elementi troppo corti sono allargati a
    una larghezza minima (scala per elemento, dichiarata in didascalia).
    """
    L_mm = [el.lunghezza / MM for el in a.elementi]
    s = _scala_orizzontale(L_mm)
    larghezze = [max(L * s, LARGHEZZA_MIN_CM) for L in L_mm]
    allargati = any(L * s < LARGHEZZA_MIN_CM - 1e-9 for L in L_mm)
    D_max = max([st.d_out / MM for el in a.elementi for sg in el.segmenti for st in sg.strati]
                + [20.0])
    sy = ALTEZZA_MAX_CM / D_max  # cm per mm di diametro
    r_man = 0.5 * 16.0 * sy

    t: list[str] = [r"\begin{tikzpicture}[font=\scriptsize]"]
    x = 0.0
    x_nodi = [0.0]
    for el, w in zip(a.elementi, larghezze):
        # sopra gli elementi su cuscinetto c'e' il simbolo del cuscinetto: la quota parte piu' in alto
        alzo = 0.45 if el.id in a.elementi_cuscinetto else 0.08
        if el.modo == "telescope":
            h = 0.5 * ALTEZZA_MAX_CM
            t.append(rf"\draw[rounded corners=3pt, fill=purple!10, draw=purple!70!black] "
                     rf"({x + 0.04:.3f},{-h:.3f}) rectangle ({x + w - 0.04:.3f},{h:.3f});")
            t.append(rf"\node[purple!70!black, align=center, font=\scriptsize] at "
                     rf"({x + w / 2:.3f},{0.45 * h:.3f}) {{telescope\\(no shaft)}};")
            t.append(rf"\node[rotate=90, anchor=west] at ({x + w / 2:.3f},{h + 0.1:.3f}) "
                     rf"{{$L{el.lunghezza / MM:g}$}};")
        elif el.modo == "manual":
            t.append(rf"\draw[pattern=crosshatch, pattern color=gray!60] ({x:.3f},{-r_man:.3f}) "
                     rf"rectangle ({x + w:.3f},{r_man:.3f});")
            t.append(rf"\node[rotate=90, anchor=west] at ({x + w / 2:.3f},{r_man + 0.1:.3f}) "
                     rf"{{CAD: $J={num(el.J_manuale, 3)}$, $k={num(el.k_manuale, 3)}$, "
                     rf"$L={el.lunghezza / MM:g}$}};")
        else:
            sc = w / el.lunghezza
            xs = x
            for sg in el.segmenti:
                ws = sg.lunghezza * sc
                strati = sorted(sg.strati, key=lambda st: -st.d_out)
                for st in strati:
                    R, r = 0.5 * st.d_out / MM * sy, 0.5 * st.d_in / MM * sy
                    bande = ((-R, R),) if st.forma == "solid" else ((r, R), (-R, -r))
                    for y0, y1 in bande:
                        t.append(rf"\filldraw[fill={_colore(st)}, draw=black, line width=0.3pt] "
                                 rf"({xs:.3f},{y0:.3f}) rectangle ({xs + ws:.3f},{y1:.3f});")
                r_foro = min((st.d_in for st in sg.strati), default=0.0)
                if all(st.forma == "hollow" for st in sg.strati) and r_foro > 0:
                    rf_ = 0.5 * r_foro / MM * sy
                    t.append(rf"\fill[pattern=north east lines, pattern color=gray] "
                             rf"({xs:.3f},{-rf_:.3f}) rectangle ({xs + ws:.3f},{rf_:.3f});")
                R_top = 0.5 * max(st.d_out for st in sg.strati) / MM * sy
                quote = " + ".join(
                    (rf"$\varnothing{st.d_out / MM:g}$" if st.forma == "solid"
                     else rf"$\varnothing{st.d_out / MM:g}/{st.d_in / MM:g}$")
                    for st in strati)
                t.append(rf"\node[rotate=90, anchor=west, inner sep=1pt] at "
                         rf"({xs + ws / 2:.3f},{R_top + alzo:.3f}) {{{quote}, $L{sg.lunghezza / MM:g}$}};")
                xs += ws
        x += w
        x_nodi.append(x)

    y_bot = -0.5 * ALTEZZA_MAX_CM - (0.6 if a.elementi_cuscinetto else 0.35)
    for i, xn in enumerate(x_nodi):
        t.append(rf"\draw[dashed, gray] ({xn:.3f},{y_bot + 0.2:.3f}) -- ({xn:.3f},{-y_bot - 0.6:.3f});")
        t.append(rf"\node[circle, draw, inner sep=1pt] at ({xn:.3f},{y_bot:.3f}) {{{i + 1}}};")
    for el, x0, x1 in zip(a.elementi, x_nodi[:-1], x_nodi[1:]):
        t.append(rf"\node[gray] at ({0.5 * (x0 + x1):.3f},{y_bot:.3f}) {{e{el.id}}};")

    simbolo = r"\psi" if a.nome == "azimut" else r"\theta"
    xm = x_nodi[a.nodo_coppia - 1]
    t.append(rf"\draw[->, thick, red!70!black] ({xm - 0.25:.3f},{y_bot - 0.45:.3f}) arc[start angle=200, "
             rf"end angle=-20, x radius=0.27, y radius=0.12] node[right] {{$\tau_{simbolo}$}};")
    xe = x_nodi[a.nodo_encoder - 1]
    t.append(rf"\node[draw, fill=green!15, anchor=east] at ({xe - 0.1:.3f},0) {{encoder}};")
    xl = x_nodi[a.nodo_estremo - 1]
    etichetta = "load (fork)" if a.nome == "azimut" else "free end"
    t.append(rf"\node[draw, fill=gray!15, anchor=west] at ({xl + 0.1:.3f},0) {{{etichetta}}};")

    if a.nome == "polare":
        for e in a.elementi_cuscinetto:
            x0, x1 = x_nodi[e - 1], x_nodi[e]
            xc = 0.5 * (x0 + x1)
            R = 0.5 * max((st.d_out for sg in a.elemento(e).segmenti for st in sg.strati),
                          default=16 * MM) / MM * sy
            for segno in (1, -1):
                yc = segno * (R + 0.22)
                t.append(rf"\draw[fill=white] ({xc - 0.25:.3f},{yc - 0.15:.3f}) rectangle "
                         rf"({xc + 0.25:.3f},{yc + 0.15:.3f});")
                t.append(rf"\draw ({xc:.3f},{yc:.3f}) circle (0.1);")
            t.append(rf"\node[font=\scriptsize] at ({xc:.3f},{y_bot - 0.42:.3f}) {{bearing}};")
        ep = a.elemento_payload
        assert ep is not None and a.xi_payload is not None
        x0, x1 = x_nodi[ep - 1], x_nodi[ep]
        xP = x0 + a.xi_payload * (x1 - x0)
        yT = 0.5 * ALTEZZA_MAX_CM + 1.9
        if a.elemento(ep).modo == "telescope":
            t.append(rf"\draw[thick, purple] ({x0:.3f},{y_bot + 0.25:.3f}) -- ({x1:.3f},{y_bot + 0.25:.3f});")
            t.append(rf"\node[purple!70!black, font=\scriptsize] at ({xP:.3f},{y_bot - 0.42:.3f}) "
                     rf"{{{ep} = P = {ep + 1}}};")
        else:
            t.append(rf"\draw[decorate, decoration={{brace, amplitude=4pt}}] ({x0:.3f},{yT - 0.35:.3f}) -- "
                     rf"({x1:.3f},{yT - 0.35:.3f}) node[midway, above=4pt] {{telescope clamp}};")
            t.append(rf"\draw[thick, purple] ({xP:.3f},{-0.5 * ALTEZZA_MAX_CM:.3f}) -- "
                     rf"({xP:.3f},{yT - 0.45:.3f});")
        t.append(rf"\node[circle, fill=purple, text=white, inner sep=1pt] at ({xP:.3f},{y_bot:.3f}) {{P}};")
    t.append(r"\end{tikzpicture}")
    return "\n".join(t), allargati


# ---------------------------------------------------------------------------
# Tabelle
# ---------------------------------------------------------------------------


def tabella_elementi(asse: str, props: Sequence[ProprietaElemento]) -> str:
    """Contributi per strato e totali di ogni elemento (longtable).

    13 colonne: corpo ``\\footnotesize`` e ``\\tabcolsep`` ridotto, altrimenti la
    tabella esce dal margine destro (una longtable non si puo' scalare).
    """
    r: list[str] = [
        r"\begingroup\footnotesize\setlength{\tabcolsep}{3.2pt}",
        r"\begin{longtable}{llllrrrlcrrrr}",
        rf"\caption{{Element properties, {asse} axis: layer contributions and element totals "
        r"($m$ in g, inertias in kg\,m$^2$, $k$ in N\,m/rad).}\\",
        r"\toprule",
        r"el. & seg. & lay. & shape & $D_o$ & $D_i$ & $L$ & mat. & carr. & $m$ & $J$ & $I_t$ & $k$ \\",
        r"\midrule\endfirsthead",
        r"\toprule el. & seg. & lay. & shape & $D_o$ & $D_i$ & $L$ & mat. & carr. & $m$ & $J$ & $I_t$ & $k$ \\ "
        r"\midrule\endhead",
    ]
    for p in props:
        for st in p.strati:
            k_s = st.GJ_p / st.lunghezza if st.GJ_p > 0 else None
            r.append(rf"{p.etichetta} & {st.segmento} & {st.strato} & {st.forma} & {st.d_out / MM:g} & "
                     rf"{st.d_in / MM:g} & {st.lunghezza / MM:g} & {esc(in_inglese(st.materiale))} & "
                     rf"{'yes' if st.portante else 'no'} & ${num(st.m * 1e3, 3)}$ & ${num(st.J, 3)}$ & "
                     rf"${num(st.I_t, 3)}$ & ${num(k_s, 3) if k_s else '--'}$ \\")
        for pz in p.pezzi:
            if "manuale" in pz.nome or "extra" in pz.nome:
                r.append(rf"{p.etichetta} & \multicolumn{{8}}{{l}}{{{esc(pz.nome)}}} & ${num(pz.m * 1e3, 3)}$ & "
                         rf"${num(pz.J, 3)}$ & ${num(pz.I_t, 3)}$ & -- \\")
        r.append(rf"\midrule \textbf{{{p.etichetta}}} & \multicolumn{{5}}{{l}}{{total ({p.modo})}} & "
                 rf"{p.lunghezza / MM:g} & & & ${num(p.m * 1e3, 4)}$ & ${num(p.J, 4)}$ & ${num(p.I_t, 4)}$ & "
                 rf"${num(p.k, 4)}$ \\ \midrule")
    r.append(r"\bottomrule\end{longtable}")
    r.append(r"\endgroup")
    return "\n".join(r)


def tabella_validazione(confronti: Sequence[ConfrontoRiduzione]) -> str:
    r = [r"\begin{longtable}{llllrrrr}",
         r"\caption{Natural frequencies of the isolated axes: full FEM model against the Guyan-reduced "
         r"chain (Hz), relative error, Craig--Bampton with 4 fixed-interface modes, and the frequency "
         r"above which the reduced torque$\to$encoder FRF departs by more than 3\,dB.}\label{tab:validation}\\",
         r"\toprule axis & mass & reduced mass & mode & full & reduced & error & CB4 \\ \midrule\endfirsthead",
         r"\toprule axis & mass & reduced mass & mode & full & reduced & error & CB4 \\ \midrule\endhead"]
    for c in confronti:
        n_cb = max(c.f_cb) if c.f_cb else None
        for k in range(len(c.f_red)):
            cb = num(c.f_cb[n_cb][k], 5) if n_cb is not None else "--"
            r.append(rf"{in_inglese(c.asse)} & {c.modello_massa} & {esc(c.forma_massa)} & {k + 1} & "
                     rf"${num(c.f_full[k], 5)}$ & ${num(c.f_red[k], 5)}$ & "
                     rf"${num(c.errore_relativo[k] * 100, 3)}\,\%$ & ${cb}$ \\")
        r.append(rf"\multicolumn{{8}}{{r}}{{\footnotesize reduced FRF valid up to "
                 rf"${num(c.f_limite, 4)}$\,Hz}} \\ \midrule")
    r.append(r"\bottomrule\end{longtable}")
    return "\n".join(r)


def tabella_poli_zeri(ml: ModelloLineare) -> str:
    r = [r"\begin{longtable}{lrrrl}",
         rf"\caption{{Poles and invariant zeros of the torque$\to$encoder matrix at "
         rf"$\theta_0 = {np.rad2deg(ml.theta0):g}^\circ$, $\Omega_0 = {np.rad2deg(ml.Omega0):g}^\circ$/s "
         r"(only roots with $\operatorname{Im} \ge 0$).}\label{tab:polezero}\\",
         r"\toprule quantity & Re [rad/s] & Im [rad/s] & $|\lambda|/2\pi$ [Hz] & type \\ \midrule\endfirsthead",
         r"\toprule quantity & Re [rad/s] & Im [rad/s] & $|\lambda|/2\pi$ [Hz] & type \\ \midrule\endhead"]
    for rad in tf.classifica(ml.poli):
        if np.imag(rad.valore) < 0:
            continue
        r.append(rf"pole & ${num(np.real(rad.valore), 4)}$ & ${num(np.imag(rad.valore), 5)}$ & "
                 rf"${num(rad.f_hz, 5)}$ & {esc(rad.tipo)} \\")
    nomi = {("psi_e", 0): r"$G_{\psi_e\tau_\psi}$", ("psi_e", 1): r"$G_{\psi_e\tau_\theta}$",
            ("theta_e", 0): r"$G_{\theta_e\tau_\psi}$", ("theta_e", 1): r"$G_{\theta_e\tau_\theta}$"}
    for (u, j), nome in nomi.items():
        r.append(r"\midrule")
        z = tf.zeri_invarianti(ml.A, ml.Bu[:, j], ml.Cy[ml.riga(u)])
        for rad in tf.classifica(z):
            if np.imag(rad.valore) < 0:
                continue
            r.append(rf"zero {nome} & ${num(np.real(rad.valore), 4)}$ & ${num(np.imag(rad.valore), 5)}$ & "
                     rf"${num(rad.f_hz, 5)}$ & {esc(rad.tipo)} \\")
    r.append(r"\bottomrule\end{longtable}")
    return "\n".join(r)


# ---------------------------------------------------------------------------
# Matrici simboliche
# ---------------------------------------------------------------------------


def _etichetta_k(e: str) -> str:
    return f"k_{{{e}}}"


def matrici_simboliche_fem(mod: fem.ModelloAsse) -> tuple[str, str, str]:
    """K, M concentrata e M consistente globali in forma simbolica.

    Per la consistente si usano le voci di elemento :math:`a_e = \\int N_1^2\\bar J`,
    :math:`b_e = \\int N_1 N_2\\bar J`, :math:`c_e = \\int N_2^2\\bar J` (piu'
    :math:`J^{x}_e/2` sulla diagonale per le inerzie extra).
    """
    n = mod.n
    K = sp.zeros(n, n)
    Ml = sp.zeros(n, n)
    Mc = sp.zeros(n, n)
    for p, (i, j) in zip(mod.elementi, mod.collegamenti):
        e = p.etichetta
        k = sp.Symbol(_etichetta_k(e), positive=True)
        J = sp.Symbol(f"J_{{{e}}}", positive=True)
        a, b, c = (sp.Symbol(f"{s}_{{{e}}}", positive=True) for s in "abc")
        for (r, s_), v in (((i, i), k), ((j, j), k), ((i, j), -k), ((j, i), -k)):
            K[r, s_] += v
        Ml[i, i] += J / 2
        Ml[j, j] += J / 2
        Mc[i, i] += a
        Mc[j, j] += c
        Mc[i, j] += b
        Mc[j, i] += b
    for i, Jn in enumerate(mod.J_nodali):
        if Jn > 0:
            s = sp.Symbol(f"J^{{n}}_{{{mod.etichette[i]}}}", positive=True)
            Ml[i, i] += s
            Mc[i, i] += s
    return L(K), L(Ml), L(Mc)


# ---------------------------------------------------------------------------
# Documento
# ---------------------------------------------------------------------------


@dataclass
class ContenutoLatex:
    par: Parametri
    sistemi: Mapping[tuple[str, str], Sistema]
    confronti: Sequence[ConfrontoRiduzione]
    modelli_lineari: Sequence[ModelloLineare]
    figure: Mapping[str, str]
    """Nome logico -> percorso relativo del PDF della figura (da results/latex)."""
    avvisi: Sequence[Avviso]
    ml_tabella: ModelloLineare | None = None
    note: list[str] = field(default_factory=list)
    simulazioni: Sequence = ()
    riduzione: str = "guyan"
    """Riduzione del modello principale (simulazioni): modal | guyan | craig_bampton."""
    modale: object | None = None
    """:class:`gimbal_mdof.modal_output.RisultatiModali` (sezione del troncamento modale)."""
    s1_riferimento: object | None = None
    """S1 sul modello a tutti i GdL del FEM (confronto)."""
    includi_guyan: bool = True
    """``latex.include_guyan``: con ``False`` (e il modello modale disponibile) il documento non
    contiene la riduzione di Guyan ne' i confronti con essa."""
    ml_principale: ModelloLineare | None = None
    """Modello lineare principale (modale) per la tabella di poli e zeri senza Guyan."""
    concentrato: Sistema | None = None
    """Modello a parametri concentrati (3 + 4 GdL), sezione ``sec:lumped``."""
    concentrato_extra: Sequence[str] = ()
    """Sottosezioni aggiuntive della sezione del modello concentrato (FdT, simulatore)."""
    sperimentale_tex: str = ""
    """Sezione del confronto con le misure sperimentali (``sec:experimental``), gia' in LaTeX."""


def _con_guyan(c: ContenutoLatex) -> bool:
    """Guyan nel documento: se richiesto, oppure se manca il modello modale (unica riduzione)."""
    return c.includi_guyan or c.modale is None


PREAMBOLO = r"""\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[english]{babel}
\usepackage[margin=2.2cm]{geometry}
\usepackage{amsmath,amssymb,mathtools,bm}
\usepackage{booktabs,longtable,array,tabularx}
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}
\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
\usepackage{graphicx}
\usepackage{float}
\usepackage[font=small]{caption}
\usepackage{xcolor}
\usepackage{tikz}
\usetikzlibrary{patterns,decorations.pathreplacing}
\usepackage[hidelinks]{hyperref}
\setcounter{MaxMatrixCols}{20}
\numberwithin{equation}{section}
\newcommand{\avviso}[1]{\begin{center}\fcolorbox{red!70!black}{red!5}{\parbox{0.92\linewidth}{\small\textbf{Warning.} #1}}\end{center}}
"""


def _sezione_sistema(par: Parametri, modale: bool = False, guyan: bool = True) -> str:
    if modale and not guyan:
        stadi = r"""\item reduction of each axis by modal truncation in the sense of Gawronski of the full FEM (rigid
      mode plus the flexible modes in the band of interest; \texttt{modal.py}, Section~\ref{sec:modal});
\item coupling of the two reduced axes through the kinetic energy of the whole system
      (\texttt{coupled.py}, \texttt{modal.py}, Sections~\ref{sec:coupled} and~\ref{sec:modal})."""
        coord = r"""\paragraph{Reduced coordinates.} $\boldsymbol\eta = [\eta_{\psi0}, \eta_{\psi k}, \dots,\
\eta_{\theta0}, \eta_{\theta k}, \dots]^{\mathsf T}$, with $\eta_{a0}$ the rigid rotation of axis $a$ and $\eta_{ak}$
the amplitudes of the retained flexible modes (Section~\ref{sec:modal}). The full FEM coordinates are the
rotations of all nodes, $\mathbf q_\psi = [\psi_1, \dots, \psi_9]^{\mathsf T}$ and $\mathbf q_\theta$ (polar
nodes, with P)."""
    elif modale:
        stadi = r"""\item reduction of each axis: by default a modal truncation in the sense of Gawronski of the full FEM
      (rigid mode plus the flexible modes in the band of interest; \texttt{modal.py},
      Section~\ref{sec:modal}); the Guyan reduction to a short chain on master nodes (\texttt{reduction.py},
      Section~\ref{sec:guyan}) and Craig--Bampton are kept for comparison;
\item coupling of the two reduced axes through the kinetic energy of the whole system
      (\texttt{coupled.py}, \texttt{modal.py}, Sections~\ref{sec:coupled} and~\ref{sec:modal})."""
        coord = r"""\paragraph{Reduced coordinates.} Modal model (default): $\boldsymbol\eta = [\eta_{\psi0},
\eta_{\psi k}, \dots,\ \eta_{\theta0}, \eta_{\theta k}, \dots]^{\mathsf T}$, with $\eta_{a0}$ the rigid rotation of
axis $a$ and $\eta_{ak}$ the amplitudes of the retained flexible modes (Section~\ref{sec:modal}). Guyan chains
(comparison):
\begin{gather}
\mathbf q = [\psi_e,\ \psi_m,\ \psi_L,\ \theta_e,\ \theta_p,\ \theta_m,\ \theta_r]^{\mathsf T},
\end{gather}
with masters $\{1, 5, 9\}$ (azimuth) and $\{1, \mathrm P, 7, 11\}$ (polar)."""
    else:
        stadi = r"""\item Guyan reduction of each axis to a short chain on master nodes (\texttt{reduction.py},
      Section~\ref{sec:guyan});
\item coupling of the two reduced chains through the kinetic energy of the whole system
      (\texttt{coupled.py}, Section~\ref{sec:coupled})."""
        coord = r"""\paragraph{Reduced coordinates.}
\begin{gather}
\mathbf q = [\psi_e,\ \psi_m,\ \psi_L,\ \theta_e,\ \theta_p,\ \theta_m,\ \theta_r]^{\mathsf T},
\end{gather}
with masters $\{1, 5, 9\}$ (azimuth) and $\{1, \mathrm P, 7, 11\}$ (polar)."""
    return r"""
\section{System, conventions and modelling chain}
The model is the second link of the chain of increasingly complex models of the alt-azimuth mount. The rigid
two-degree-of-freedom model (\texttt{gimbal\_2dof}, document \texttt{modello\_matematico.tex}) is extended with
the torsional flexibility of the azimuth and polar (elevation) shafts. The computation has four stages, each
implemented in a module of \texttt{src/gimbal\_mdof/}:
\begin{enumerate}
\item element properties from the composition of solid and hollow cylinders, or from CAD data
      (\texttt{geometry.py}, Section~\ref{sec:elements});
\item a torsional FEM of each axis in isolation, with lumped and consistent mass (\texttt{fem.py},
      Section~\ref{sec:fem});
@@STADI@@
\end{enumerate}

\paragraph{Frames and angles.} The frames and sign conventions of the rigid model are kept unchanged; they are
collected in Section~\ref{sec:frames}, together with the frame in which every inertia in \texttt{config.yaml}
must be entered. SI units are used internally.

\paragraph{Numbering and coordinates.} Nodes are numbered from 1; element $i$ joins node $i$ to node $i+1$.
Azimuth: nodes $1\dots9$, node 1 at the bottom (encoder), torque $\tau_\psi$ on node 5, node 9 carries the fork
(body 1). Polar: nodes $1\dots11$, encoder on node 1, torque $\tau_\theta$ on node 7 (motor rotor; the stator is
on the fork), free end on node 11; elements 2 and 4 run in the bearings housed in the fork; the telescope
(body 2) sits on element 3 at the node P. $\psi_i$ is the absolute rotation of azimuth node $i$ about $\hat k$;
$\theta_i$ is the rotation of polar node $i$ relative to the fork about $\hat J$, which coincides with its
absolute rotation about $\hat J$ because the fork only rotates about $\hat k$. The azimuth encoder reads
$\psi_1$, the polar encoder (stator on the fork) reads $\theta_1$.

@@COORD@@
""".replace("@@STADI@@", stadi).replace("@@COORD@@", coord)


def _direzione_N(v: np.ndarray) -> str:
    """Nome di un versore di N (es. ``-\\hat y``) o le sue componenti."""
    for i, nome in enumerate("xyz"):
        for s, segno in ((1.0, "+"), (-1.0, "-")):
            e = np.zeros(3)
            e[i] = s
            if np.allclose(v, e, atol=1e-9):
                return rf"{segno}\hat{{{nome}}}"
    return "(" + ",\\ ".join(num(c, 3) for c in v) + ")"


def _simbolico(x: float) -> sp.Expr:
    """Valore numerico di una matrice di rotazione in forma pulita (interi se esatti)."""
    return sp.Integer(int(round(x))) if abs(x - round(x)) < 1e-12 else sp.Float(x, 6)


SIMBOLI_I_N = sp.symbols(r"I^{N}_{xx} I^{N}_{yy} I^{N}_{zz} I^{N}_{xy} I^{N}_{xz} I^{N}_{yz}", real=True)
SIMBOLI_R = sp.symbols(r"\Delta{x} \Delta{y} \Delta{z}", real=True)


def mappatura_N_2(R_N2_0: np.ndarray) -> tuple[sp.Matrix, sp.Matrix]:
    """Componenti di :math:`\\mathbf I_{p,o_2}|_2` e di :math:`\\mathbf d_2` in funzione di quelle in N.

    E' la stessa trasformazione di :func:`gimbal.inertia.ruota`, in forma simbolica per
    il documento: :math:`\\mathbf I|_2 = \\mathbf R_{N2}^{(0)\\mathsf T}\\mathbf I|_N\\mathbf R_{N2}^{(0)}`,
    :math:`\\mathbf d_2 = \\mathbf R_{N2}^{(0)\\mathsf T}\\mathbf r`.
    """
    Ixx, Iyy, Izz, Ixy, Ixz, Iyz = SIMBOLI_I_N
    IN = sp.Matrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    R = sp.Matrix(3, 3, lambda i, j: _simbolico(float(R_N2_0[i, j])))
    I2 = (R.T * IN * R).applyfunc(sp.expand)
    d2 = (R.T * sp.Matrix(SIMBOLI_R)).applyfunc(sp.expand)
    return I2, d2


def tikz_terne(tr, s: float) -> str:
    """Viste schematiche della posa CAD (dall'alto e laterale), dai dati di ``frames``."""
    ia, ja = tr.i_hat_0, tr.J_hat_0
    ang_x = float(np.degrees(np.arctan2(ia[1], ia[0])))
    ang_y = float(np.degrees(np.arctan2(ja[1], ja[0])))
    th_cad = float(np.degrees(tr.theta_cad))
    z1, z2 = tr.z_o1 * 1e3, tr.z_o2 * 1e3
    h2 = 3.2
    h1 = max(0.85, h2 * z1 / z2) if z2 > 0 else 0.85
    lab_x = _direzione_N(ia)
    freccia = ang_y if s > 0 else ang_y + 180.0
    cx, sx_ = np.cos(np.radians(ang_x)), np.sin(np.radians(ang_x))
    cy, sy_ = np.cos(np.radians(ang_y)), np.sin(np.radians(ang_y))

    def p(a: float, b: float, dx: float = 0.0, dy: float = 0.0) -> str:
        """Punto a*i_hat_0 + b*J_hat_0 nella vista dall'alto, piu' uno spostamento sullo schermo."""
        return f"({a * cx + b * cy + dx:.3f},{a * sx_ + b * sy_ + dy:.3f})"

    return "\n".join([
        r"\begin{tikzpicture}[font=\small]",
        r"\begin{scope}",
        r"\draw[->, gray] (-2.4,0) -- (2.6,0) node[right] {$\hat x$};",
        r"\draw[->, gray] (0,-2.6) -- (0,2.4) node[above] {$\hat y$};",
        rf"\draw[fill=gray!15] {p(-0.4, -0.25)} -- {p(1.8, -0.25)} -- {p(1.8, 0.25)} -- {p(-0.4, 0.25)} -- cycle;",
        rf"\node[font=\scriptsize] at {p(1.25, -0.75)} {{telescope}};",
        rf"\draw[line width=2.2pt, blue!60!black] {p(0, -1.7)} -- {p(0, -0.25)};",
        rf"\draw[line width=2.2pt, blue!60!black] {p(0, 0.25)} -- {p(0, 1.7)};",
        rf"\draw[very thick, ->] (0,0) -- {p(2.2, 0)} node[anchor={'north' if ia[1] < -0.5 else 'west'}] "
        r"{$\hat\imath_0 = \hat K_0$};",
        rf"\draw[very thick, ->, blue!60!black] (0,0) -- {p(0, 2.2)};",
        rf"\node[blue!60!black, anchor=south] at {p(0, 1.8, 0.0, 0.3)} {{$\hat J_0 = -\hat\jmath_0$}};",
        rf"\draw[->, blue!60!black, dashed] {p(0, -1.6 if s > 0 else -0.7, 0.0, -0.4)} -- ++({freccia:.3f}:0.9);",
        rf"\node[blue!60!black, font=\scriptsize, anchor=north] at {p(0, -1.15, 0.0, -0.5)} {{polar node index}};",
        r"\fill (0,0) circle (1.6pt);",
        rf"\node[anchor=south east, font=\scriptsize] at {p(0, 0, -0.3, 0.08)} {{$O_N, o_1, o_2$}};",
        r"\node at (0,-3.1) {(a) from above, CAD pose ($\psi = 0$)};",
        r"\end{scope}",
        r"\begin{scope}[xshift=6.9cm, yshift=-2.5cm]",
        rf"\draw[->, gray] (-1.2,0) -- (4.1,0) node[below] {{$\hat\imath_0$ $({lab_x})$}};",
        r"\draw[->, gray] (0,-0.2) -- (0,4.4) node[right] {$\hat z = \hat k$};",
        r"\fill (0,0) circle (1.6pt) node[below left] {$O_N$};",
        rf"\draw[line width=2.2pt, blue!60!black] (0,0.1) -- (0,{h2 - 0.6:.3f});",
        rf"\fill ({0:.3f},{h1:.3f}) circle (1.6pt) node[above right] {{$o_1$}};",
        rf"\node[right, font=\scriptsize] at (0.08,{h2 - 0.75:.3f}) {{az. node 9 (fork)}};",
        r"\node[right, font=\scriptsize] at (0.08,0.12) {az. node 1};",
        rf"\draw[fill=gray!15, rotate around={{{-th_cad:.3f}:(0,{h2})}}] (-0.35,{h2 - 0.18}) rectangle (3.0,{h2 + 0.18});",
        rf"\draw[very thick, ->, rotate around={{{-th_cad:.3f}:(0,{h2})}}] (0,{h2}) -- (3.7,{h2}) node[right] {{$\hat K_0$}};",
        rf"\draw[dashed, ->] (0,{h2}) -- ++(15:3.4) node[right] {{$\hat K,\ \theta>0$}};",
        rf"\draw[thick, ->] (0,{h2}) ++(0.0,0.15) -- ++(0,0.75) node[left, font=\scriptsize] {{$\hat I_0$}};",
        rf"\draw[->] (1.5,{h2}) arc (0:15:1.5) node[midway, right] {{$\theta$}};",
        rf"\draw[fill=white] (0,{h2}) circle (0.12); \fill (0,{h2}) circle (0.035);",
        rf"\node at (-0.35,{h2 + 0.35}) {{$o_2$}};",
        rf"\node[blue!60!black, font=\scriptsize, anchor=east] at (-0.45,{h2}) {{$\hat J$ (out of page)}};",
        rf"\draw[<->] (-2.0,0) -- (-2.0,{h2 - 0.2}) node[midway, left] {{${z2:.6g}$\,mm}};",
        rf"\draw[<->] (-0.35,0) -- (-0.35,{h1:.3f}) node[midway, left, font=\scriptsize] {{${z1:.6g}$\,mm}};",
        r"\node at (1.4,-0.7) {(b) side view, plane $(\hat\imath_0, \hat z)$};",
        r"\end{scope}",
        r"\end{tikzpicture}",
    ])


def _sezione_terne(par: Parametri, sistema: Sistema) -> str:
    """Sistemi di riferimento (convenzione della tesi De Donatis), coerenti con il codice."""
    tr = par.terne
    pol = par.polare
    s = pol.verso_numerazione
    I2, d2 = mappatura_N_2(tr.R_N2_0)
    nomi = [("II", 0, 0), ("JJ", 1, 1), ("KK", 2, 2), ("IJ", 0, 1), ("IK", 0, 2), ("JK", 1, 2)]
    mappa = r",\quad ".join(rf"I_{{{n}}} = {sp.latex(I2[i, j])}" for n, i, j in nomi[:3])
    mappa2 = r",\quad ".join(rf"I_{{{n}}} = {sp.latex(I2[i, j])}" for n, i, j in nomi[3:])
    mappa_d = r",\quad ".join(rf"d_{c} = {sp.latex(d2[k])}" for k, c in enumerate("IJK"))
    r1, r2 = tr.r_o1 * 1e3, tr.r_o2 * 1e3
    i_lab, j_lab = _direzione_N(tr.i_hat_0), _direzione_N(tr.J_hat_0)
    verso_lab = r"+\hat J" if s > 0 else r"-\hat J"
    y1 = pol.y_nodo1_da_O
    y1_testo = (rf"$y_1 = {num(y1 * 1e3, 4)}$\,mm" if y1 is not None else
                r"$y_1$ is not given (\texttt{null}, Q7): the code then sets $y_c = 0$ for every piece, "
                r"i.e.\ no Steiner term in $J^{ax}$")
    pp = par.payload
    corrente = (r"\texttt{source: gimbal\_2dof} (body B of the rigid model)" if pp.sorgente == "gimbal_2dof"
                else rf"\texttt{{source: manual}}, \texttt{{frame: {esc(pp.terna)}}}, "
                     rf"\texttt{{subtract\_shaft\_parts: {str(pp.sottrai_alberi).lower()}}}"
                     + (r" (\textbf{placeholder values})" if pp.segnaposto else ""))
    mod_pol = sistema.asse_theta
    telescopio = getattr(mod_pol, "telescopio_rigido", False)
    ep = pol.elemento_payload
    return "\n".join([
        r"\section{Reference frames}\label{sec:frames}",
        r"The frames follow the convention of the De Donatis thesis: an inertial frame and one body frame "
        r"for each rotating body. The origins of the body frames lie on the rotation axes and \emph{not} at "
        r"the centres of mass, so that the eccentricity of each centre of mass with respect to its centre of "
        r"rotation is a separate, explicit quantity. The data are in the section \texttt{frames} of "
        r"\texttt{config.yaml}; the frames are built by \texttt{gimbal.params.costruisci\_terne} and the "
        r"inertia tensors are transported and rotated by \texttt{gimbal.inertia.trasporta} and "
        r"\texttt{gimbal.inertia.ruota}, called by \texttt{gimbal\_mdof.payload.tensore\_da\_config}.",
        r"\begin{table}[H]\centering\small",
        r"\caption{Frames of the model.}\label{tab:frames}",
        r"\begin{tabularx}{\linewidth}{L{0.16\linewidth}L{0.2\linewidth}YL{0.2\linewidth}}\toprule",
        r"frame & origin & unit vectors & role \\ \midrule",
        r"$N$, inertial & $O_N$, centre of the base of the azimuth shaft & $\hat x, \hat y, \hat z$: axes of "
        r"the Onshape assembly, $\hat z$ vertical upwards along the azimuth axis & reference for all CAD "
        r"exports \\",
        r"1, body 1 (azimuth, fork) & $o_1$, on the azimuth axis & $\hat\imath$ horizontal projection of the "
        r"pointing direction, $\hat k = \hat z$, $\hat\jmath = \hat k\times\hat\imath$ (along the polar axis) "
        r"& rotates with $\psi_L$ about "
        r"$\hat k$ \\",
        r"2, body 2 (polar, telescope) & $o_2$, centre of rotation of the polar shaft & $\hat K$ pointing "
        r"direction (positive towards the optics), $\hat I$ upwards at $\theta = 0$, "
        r"$\hat J = \hat K\times\hat I = -\hat\jmath$ polar axis & rotates with "
        r"$\theta_p$ about $\hat J$ relative to body 1 \\ \bottomrule",
        r"\end{tabularx}\end{table}",
        r"\paragraph{Why $\hat J = -\hat\jmath$.} Body 2 is fixed by three requirements: $\hat I$ coincides "
        r"with $\hat k = \hat z$ when the telescope is horizontal, $\hat K$ points out of the optics, and "
        r"$\theta > 0$ raises the pointing. A right-handed frame then leaves no choice for the third axis: "
        r"$\hat J = \hat K\times\hat I$, which at $\theta = 0$ ($\hat K = \hat\imath$, $\hat I = \hat k$) is "
        r"$\hat\imath\times\hat k = -\hat\jmath$. A positive rotation $\theta$ about this $\hat J$ gives "
        r"$\hat K = \cos\theta\,\hat\imath + \sin\theta\,\hat k$, i.e.\ it raises the pointing, as required. "
        r"$\hat\jmath$ is fixed independently by body 1 ($\hat\jmath = \hat k\times\hat\imath$, with $\hat\imath$ "
        r"the horizontal projection of the pointing), so it points the opposite way. $\hat J = +\hat\jmath$ "
        r"would require either $\theta > 0$ to lower the pointing (the convention of the two-DOF model) or a "
        r"left-handed frame. The sign only affects how $\theta$ and the tensor components are written; the "
        r"physics is unchanged.",
        r"$O_N$, $o_1$ and $o_2$ all lie on $\hat z$ and differ only by a translation along $\hat z$; $o_2$ is "
        r"the intersection of the azimuth and polar axes:",
        r"\begin{gather}",
        rf"\mathbf r_{{o_1}}\big|_N = [{num(r1[0], 5)},\ {num(r1[1], 5)},\ {num(r1[2], 6)}]^{{\mathsf T}}\ "
        r"\text{mm}, \qquad",
        rf"\mathbf r_{{o_2}}\big|_N = [{num(r2[0], 5)},\ {num(r2[1], 5)},\ {num(r2[2], 6)}]^{{\mathsf T}}\ "
        r"\text{mm}, \\",
        r"\mathbf R_{N1} = \mathbf R_z(\psi_L),\qquad \mathbf R_{12} = \mathbf R_y(-\theta_p)\,\boldsymbol\Pi"
        r" = \begin{bmatrix}-\sin\theta_p & 0 & \cos\theta_p\\ 0&-1&0\\ \cos\theta_p & 0 & \sin\theta_p"
        r"\end{bmatrix},\qquad \mathbf R_{N2} = \mathbf R_{N1}\mathbf R_{12}, \\",
        r"\mathbf R_z(\psi) = \begin{bmatrix}\cos\psi & -\sin\psi & 0\\ \sin\psi & \cos\psi & 0\\ 0&0&1"
        r"\end{bmatrix}, \qquad",
        r"\mathbf R_y(\theta) = \begin{bmatrix}\cos\theta & 0 & \sin\theta\\ 0&1&0\\ -\sin\theta & 0 & "
        r"\cos\theta\end{bmatrix}.",
        r"\end{gather}",
        r"with $\boldsymbol\Pi = \left[\begin{smallmatrix}0&0&1\\0&-1&0\\1&0&0\end{smallmatrix}\right]$ the "
        r"constant rotation that takes body 1 to body 2 at $\theta = 0$ ($\hat K = \hat\imath$, "
        r"$\hat I = \hat k$, $\hat J = -\hat\jmath$). $\mathbf R_{XY}$ maps components in frame $Y$ to "
        r"components in frame $X$; a tensor changes frame as "
        r"$\mathbf I|_Y = \mathbf R_{XY}^{\mathsf T}\,\mathbf I|_X\,\mathbf R_{XY}$. $\psi$ is positive "
        r"counter-clockwise about $\hat k$ seen from above. $\theta$ is positive by the right-hand rule about "
        r"$\hat J$, so $\theta>0$ \emph{raises} the pointing: $\theta$ is the elevation ($\theta = 0$ horizon, "
        r"$\theta = 90^\circ$ zenith). This is the opposite of the rigid two-degree-of-freedom model, where "
        r"$\theta>0$ lowered the pointing; internally the code keeps the pointing frame of that model and uses "
        r"$\theta_{2dof} = -\theta$.",
        r"\paragraph{CAD pose.} The CAD assembly is a snapshot of the mount at $\psi = 0$ and "
        r"$\theta = \theta_{cad}$. Given the pointing direction $\hat\imath_0$ in $N$ "
        r"(\texttt{pointing\_cad}),",
        r"\begin{gather}",
        r"\hat\jmath_0 = \hat z\times\hat\imath_0, \qquad \mathbf R_{N1}^{(0)} = \bigl[\,\hat\imath_0\ \big|\ "
        r"\hat\jmath_0\ \big|\ \hat z\,\bigr], \qquad \mathbf R_{N2}^{(0)} = \bigl[\,\hat I_0\ \big|\ "
        r"\hat J_0\ \big|\ \hat K_0\,\bigr] = \mathbf R_{N1}^{(0)}\,\mathbf R_{12}(\theta_{cad}),",
        r"\end{gather}",
        r"so that $\hat K_0$ is the pointing direction, $\hat I_0$ points upwards when the telescope is "
        r"horizontal and $\hat J_0 = \hat K_0\times\hat I_0 = -\hat\jmath_0$. Internally the code works in the "
        r"pointing frame of the rigid model (first axis along $\hat K$, second along $\hat\jmath$) and "
        r"$\boldsymbol\Pi$ converts every component printed in this document.",
        rf"The code checks that $\hat J_0 = -\hat\jmath_0$ agrees with \texttt{{polar\_axis\_cad}}. Currently "
        rf"$\hat\imath_0 = {i_lab}$, $\hat J_0 = {j_lab}$, $\theta_{{cad}} = "
        rf"{num(np.degrees(tr.theta_cad), 4)}^\circ$:",
        r"\begin{gather}",
        r"\mathbf R_{N1}^{(0)} = " + matrice_numerica(tr.R_N1_0, 4) + r", \qquad "
        r"\mathbf R_{N2}^{(0)} = " + matrice_numerica(tr.R_N2_0, 4) + r".",
        r"\end{gather}",
        rf"In words: in the CAD pose the telescope points along $\hat K_0 = {_direzione_N(tr.R_N2_0[:, 2])}$ "
        rf"(towards the optics), the polar axis is $\hat J_0 = {_direzione_N(tr.R_N2_0[:, 1])}$ and "
        rf"$\hat I_0 = {_direzione_N(tr.R_N2_0[:, 0])}$.",
        r"\begin{figure}[H]\centering",
        r"\resizebox{0.95\linewidth}{!}{" + tikz_terne(tr, s) + "}",
        r"\caption{CAD pose of the mount, drawn from the data of \texttt{frames}. (a) View from above: "
        r"telescope along $\hat\imath_0 = \hat K_0$, polar shaft along $\hat J_0 = -\hat\jmath_0$ (blue) with "
        r"the direction of increasing polar node index; $O_N$, $o_1$, $o_2$ project onto the same point. (b) "
        r"Side view in the vertical plane of the pointing direction: origins $O_N$, $o_1$, $o_2$ on "
        r"$\hat z$, azimuth shaft from node 1 to node 9 (fork), $\hat I_0$ upwards, $\hat J$ out of the page, "
        r"$\theta>0$ raising the pointing.}"
        r"\label{fig:frames}",
        r"\end{figure}",
        r"\subsection{Eccentricities}\label{subsec:eccentricity}",
        r"Because the body origins are on the rotation axes, the centre of mass of each body has an "
        r"eccentricity with respect to its centre of rotation:",
        r"\begin{gather}",
        r"\mathbf e_1 = \mathbf r_{G_1} - \mathbf r_{o_1}\ \text{(fork, in frame 1)}, \qquad "
        r"\mathbf d_2 = \mathbf R_{N2}^{(0)\mathsf T}\bigl(\mathbf r_{G_2} - \mathbf r_{o_2}\bigr)\ "
        r"\text{(telescope, in frame 2)}.",
        r"\end{gather}",
        r"$\mathbf d_2$ enters the dynamics through the tensor $\mathbf I_{p,o_2}|_2$ (Huygens--Steiner) and the "
        r"gravitational energy $m_p g\,(d_K\sin\theta_p + d_I\cos\theta_p)$ (height of $G_2$ above $o_2$); "
        r"$d_J$ lies along the polar axis "
        r"and never changes the height of the centre of mass. The fork turns only about $\hat k$, so in the "
        r"current model its eccentricity enters only through its moment of inertia about the azimuth axis, "
        r"$J_1 = I_{kk,G_1} + m_1\,(e_{1,x}^2 + e_{1,y}^2)$; keeping $o_1$ on the axis prepares the extension "
        r"in which $\mathbf e_1$ acts explicitly.",
        r"\subsection{Coordinates along the shafts}",
        r"\begin{itemize}",
        r"\item \textbf{Azimuth shaft:} along $\hat k = \hat z$, node index increasing upwards (node 1 encoder "
        r"at the bottom, node 9 fork at the top). $\psi_i$ is the absolute rotation of node $i$ about $\hat k$. "
        r"The azimuth shaft turns only about its own axis, so only the axial inertias of its parts enter.",
        r"\item \textbf{Polar shaft:} along $\hat J = -\hat\jmath$. $\theta_i$ is the rotation of node $i$ about "
        r"$\hat J$ relative to body 1, which is also its absolute rotation about $\hat J$. The polar node index "
        rf"increases along ${verso_lab}$ (\texttt{{nodes\_increase\_along}}), so the centre of mass of a polar "
        r"piece at abscissa $x$ from node 1 sits at",
        r"\begin{gather}",
        rf"y_c = y_1 + s\,x, \qquad s = {'+1' if s > 0 else '-1'},",
        r"\end{gather}",
        r"from $o_2$ along $\hat J$, where $y_1$ is \texttt{y\_node1\_from\_o2\_mm}. Currently " + y1_testo
        + r". $y_c$ enters only $J^{ax}$ (Section~\ref{sec:coupled}).",
        (rf"\item \textbf{{Telescope:}} there is no shaft through the telescope (element {ep} in "
         r"\texttt{mode: telescope}). The rigid telescope joins the polar shaft ends at nodes "
         rf"{ep} and {ep + 1}, which therefore coincide with the telescope node P; only the telescope inertia "
         r"(\texttt{payload}) is attached to P." if telescopio and par.mozzo is None else
         rf"\item \textbf{{Telescope (hub model):}} element {ep} has no shaft; its nodes coincide with the hub "
         r"centre H, and the telescope is the node P joined to H by the attachment spring $k_b$ "
         r"(Section~\ref{sec:hub}); only the telescope inertia (\texttt{payload}) is attached to P." if telescopio else
         rf"\item \textbf{{Telescope:}} clamped on the shaft element {ep}, split at P."),
        r"\item \textbf{Element data} (\texttt{mode: manual}): $J$ about the shaft's own axis; $I_t$ about a "
        r"transverse axis through the element's own centre of mass; \texttt{y\_c\_mm} measured from node $i$ in "
        r"the direction of increasing node index.",
        r"\end{itemize}",
        r"\subsection{In which frame each inertia is entered}\label{subsec:input_frames}",
        r"\paragraph{Telescope (\texttt{payload}).} The dynamics needs the tensor of the telescope about "
        r"$o_2$ in frame 2, $\mathbf I_{p,o_2}|_2$, with components $I_{II}, I_{JJ}, I_{KK}, I_{IJ}, I_{IK}, "
        r"I_{JK}$ along $\hat I, \hat J, \hat K$, the mass $m_p$ and the eccentricity $\mathbf d_2 = "
        r"[d_I, d_J, d_K]$. All tensors are in tensor convention, $I_{ab} = -\int a\,b\,dm$ for $a\neq b$. Two "
        r"input frames are accepted (key \texttt{frame}):",
        r"\begin{enumerate}",
        r"\item \texttt{frame: cad\_cg} (recommended). Enter the tensor \emph{exactly as Onshape exports it with "
        r"no mate connector selected}: about the centre of mass $G_2$ of the telescope, with axes parallel to "
        r"$\hat x, \hat y, \hat z$ (\texttt{I\_kgm2}, kg\,m$^2$), and the centre of mass $\mathbf r_{G_2}|_N$ "
        r"measured from $O_N$ (\texttt{r\_G\_mm}, mm). The code applies the sign convention "
        r"(\texttt{products\_sign}), then",
        r"\begin{gather}",
        r"\mathbf r = \mathbf r_{G_2} - \mathbf r_{o_2}, \qquad"
        r"\mathbf I_{p,o_2}|_N = \mathbf I_{p,G_2}|_N + m_p\left(\|\mathbf r\|^2\mathbf I - \mathbf r\,"
        r"\mathbf r^{\mathsf T}\right), \\"
        r"\mathbf I_{p,o_2}|_2 = \mathbf R_{N2}^{(0)\mathsf T}\,\mathbf I_{p,o_2}|_N\,\mathbf R_{N2}^{(0)}, "
        r"\qquad \mathbf d_2 = \mathbf R_{N2}^{(0)\mathsf T}\,\mathbf r .",
        r"\end{gather}",
        r"With the current CAD pose this amounts to (components of $\mathbf I_{p,o_2}|_N$ on the right, "
        r"$\Delta x, \Delta y, \Delta z$ the components of $\mathbf r$ in $N$):",
        r"\begin{gather}",
        mappa + r", \\",
        mappa2 + r", \\",
        mappa_d + ".",
        r"\end{gather}",
        r"\item \texttt{frame: body2\_o2}. Enter $\mathbf I_{p,o_2}|_2$ directly (about $o_2$, axes "
        r"$\hat I, \hat J, \hat K$) and $\mathbf d_2$ in mm (\texttt{d\_mm}); no transformation is applied.",
        r"\end{enumerate}",
        r"Checklist for an Onshape export with \texttt{frame: cad\_cg}:",
        r"\begin{itemize}",
        r"\item select only the telescope (there is no shaft through it) and keep "
        r"\texttt{subtract\_shaft\_parts: false};",
        r"\item select no mate connector, so that the tensor is about the centre of mass with axes parallel to "
        r"$\hat x, \hat y, \hat z$, and the centre of mass is measured from $O_N$;",
        rf"\item export in the CAD pose of \texttt{{frames}}: telescope pointing along ${i_lab}$ at "
        rf"$\theta_{{cad}} = {num(np.degrees(tr.theta_cad), 4)}^\circ$; if the assembly is in another pose, "
        r"update \texttt{pointing\_cad}, \texttt{polar\_axis\_cad} and \texttt{theta\_cad\_deg};",
        r"\item convert to kg\,m$^2$ and mm, and set \texttt{products\_sign} (\texttt{tensor} if the "
        r"off-diagonal terms are $-\int a\,b\,dm$, \texttt{positive} if they are $+\int a\,b\,dm$).",
        r"\end{itemize}",
        r"The code then checks that the tensor is symmetric, positive definite and satisfies the triangle "
        r"inequality, and warns (\texttt{PAYLOAD\_ASSE\_SOSPETTO}) if the smallest centroidal moment is not "
        r"about $\hat K$: for an optical tube this usually reveals a wrong pose or a wrong frame. Current "
        r"input: " + corrente + ".",
        _testo_payload(sistema),
        r"\paragraph{Fork (body 1, azimuth node 9).} A scalar: the moment of inertia of the fork alone about "
        r"the azimuth axis $\hat k$. From an export about the centre of mass with axes parallel to $N$, "
        r"$J_1 = I_{zz,G_1} + m_1\,(x_{G_1}^2 + y_{G_1}^2)$ (the azimuth axis passes through $O_N$, so the "
        r"horizontal components of $\mathbf r_{G_1}$ are those of the eccentricity $\mathbf e_1$).",
        r"\paragraph{Rotors and other node inertias.} $J$ about the axis of the shaft they sit on; for the polar "
        r"axis also $m$ and the transverse $I_t$ about their own centre of mass, located at the node.",
    ])


def _sezione_geometria(par: Parametri, sistema: Sistema) -> str:
    out = [r"\section{Input geometry}\label{sec:geometry}",
           r"The figures are generated from \texttt{config.yaml}. Each segment is drawn as a rectangle "
           r"whose width is proportional to its length and whose height is proportional to its outer "
           r"diameter; the bore of hollow cylinders is hatched, torque-carrying steel layers are blue, "
           r"aluminium light blue and non-carrying layers (spacers) orange; elements given as CAD data "
           r"(\texttt{mode: manual}) are cross-hatched. Dimensions are in mm."]
    for a, mod in ((par.azimut, sistema.asse_psi), (par.polare, sistema.asse_theta)):
        tikz, allargati = tikz_geometria(a, mod)
        nota = (" Elements shorter than the readable minimum are stretched: lengths are to scale within "
                "each element, not across elements." if allargati else
                " Lengths are to scale.")
        nome = "Azimuth" if a.nome == "azimut" else "Polar"
        ep = a.elemento_payload
        if a.nome == "azimut":
            extra = " (drawn horizontally; node 1 is at the bottom)"
        elif ep is not None and a.elemento(ep).modo == "telescope" and par.mozzo is not None:
            pm = par.mozzo
            extra = (rf"; elements {pm.elementi_lato[0]} and {pm.elementi_lato[1]} are the side hubs (tube and "
                     rf"flange) and element {ep} has no shaft: its nodes coincide with the hub centre H, from "
                     r"which the telescope hangs through $k_b$ (the support ring is in Figure~\ref{fig:hub_geometry})")
        elif ep is not None and a.elemento(ep).modo == "telescope":
            extra = (rf"; element {ep} has no shaft: the rigid telescope joins nodes {ep} and {ep + 1}, "
                     r"which coincide with the telescope node P")
        else:
            extra = (r"; P is the telescope node at $\xi_p = " + num(a.xi_payload or 0, 3)
                     + r"$ of element " + str(ep))
        out.append(r"\begin{figure}[H]\centering")
        out.append(r"\resizebox{\linewidth}{!}{" + tikz + "}")
        out.append(rf"\caption{{{nome} shaft as entered in \texttt{{config.yaml}}{extra}.{nota} "
                   r"Vertical scale: proportional to the diameter.}")
        out.append(rf"\label{{fig:geom_{a.nome}}}\end{{figure}}")
    return "\n".join(out)


def _testo_tratto_telescopio(pm, par: Parametri) -> str:
    """Paragrafo sul tratto di tubo del telescopio solidale all'anello (``hub.telescope_on_ring``)."""
    from . import mozzo as mz

    t = pm.tratto_telescopio
    m_t, J_ax, J_t = mz.inerzia_tratto(pm)
    return (
        r"\paragraph{Telescope section clamped in the ring.} The part of the telescope tube that sits inside "
        rf"the support ring ({esc(in_inglese(t.materiale.nome))}, $D_o/D_i = {num(t.d_out * 1e3, 4)}/"
        rf"{num(t.d_in * 1e3, 4)}$\,mm, length ${num(t.lunghezza * 1e3, 4)}$\,mm, centred at $o_2$, axis along "
        r"$\hat K$) is taken as rigid with the ring: it rotates with the hub and does not oscillate on $k_b$. "
        r"Its inertia is moved from the telescope to the hub centre H:"
        r"\begin{gather}"
        rf"m_t = \rho\,\tfrac{{\pi}}{{4}}(D_o^2 - D_i^2)L = {num(m_t, 3)}\ \text{{kg}},\qquad "
        rf"J_t = m_t\Bigl[\frac{{D_o^2 + D_i^2}}{{16}} + \frac{{L^2}}{{12}}\Bigr] = {num(J_t, 3)}\ \text{{kg\,m}}^2,"
        r"\end{gather}"
        rf"about the polar axis (a diameter of the tube) and about $\hat z$ at $\theta = 0$; about the pointing "
        rf"axis it carries $m_t(D_o^2 + D_i^2)/8 = {num(J_ax, 3)}$\,kg\,m$^2$. The telescope tensor loses these "
        r"terms and the mass $m_t$ (its centre of mass moves accordingly), the hub gains them, and the total "
        r"inertia of the polar axis is unchanged. $k_b$ is not modified."
    )


def _sezione_mozzo(par: Parametri, sistema: Sistema, figure: Mapping[str, str]) -> str:
    """Modello ``hub``: geometria, rigidezze in serie e inerzia del mozzo, attacco k_b del telescopio."""
    from . import mozzo as mz

    pm = par.mozzo
    assert pm is not None
    rg = mz.rigidezze(pm)
    im = mz.inerzia(pm, par.terne)
    rt = sistema.rid_theta
    ks = dict(zip(rt.archi, rt.k_catena))
    a = rt.archi
    seg = rf" \emph{{{NOTA_REALISMO.capitalize()}.}}" if par.segnaposto else ""
    righe = [rf"{esc(p)} & {esc(g)} & ${num(k, 4)}$ \\" for p, g, k in mz.righe_tabella(pm)]
    righe = [r.replace("acciaio", "steel").replace("alluminio", "aluminium") for r in righe]
    I = pm.I_G_N
    tensore = (r"\begin{bmatrix}" + r" \\ ".join(" & ".join(num(I[i, j], 4) for j in range(3)) for i in range(3))
               + r"\end{bmatrix}")
    fig = []
    if "mozzo" in figure:
        fig = [r"\begin{figure}[H]\centering",
               rf"\includegraphics[width=\linewidth,height=0.75\textheight,keepaspectratio]{{{figure['mozzo']}}}",
               r"\caption{Hub of the polar axis as used for the stiffness: (a) view from above, with the polar "
               r"axis horizontal and the telescope (ring) axis vertical; (b) view along the telescope axis; "
               r"(c) springs in series from bearing to bearing, with the telescope attached to the hub centre H "
               r"through $k_b$. Dimensions in mm, stiffnesses in N\,m/rad." + seg + r"}\label{fig:hub_geometry}"
               r"\end{figure}"]
    return "\n".join([
        r"\section{Hub and compliant telescope attachment (hub model)}\label{sec:hub}",
        r"\paragraph{Purpose and switch.} This document is generated with \texttt{polar.hub.enabled: true} "
        r"(file \texttt{config\_hub.yaml}, results in \texttt{results\_hub/}). The original model, with the "
        r"telescope rigidly attached to the shaft node P, is kept in \texttt{config.yaml} and "
        r"\texttt{results/}; the flag switches between the two. The hub model follows the antiresonance "
        r"measured at the elevation encoder (Section~\ref{sec:experimental}): a zero at an encoder on the load "
        r"side requires a large inertia outside the encoder--motor path, so the telescope is no longer rigid on "
        r"the shaft but hangs from the hub through a compliance $k_b$.",
        r"\paragraph{Topology.} The two side hubs and the support ring form the hub, which carries the torque "
        r"between the two polar bearings; its centre H (at $o_2$) carries the inertia of the hub. The two halves "
        rf"of the hub replace polar elements {pm.elementi_lato[0]} and {pm.elementi_lato[1]} as massless springs, "
        rf"element {pm.elemento_centro} (no shaft) collapses into H, and the telescope is a node P off the chain:",
        r"\begin{gather}",
        r"1 \overset{e_1}{\longleftrightarrow} 2 \overset{\text{half hub}}{\longleftrightarrow} \mathrm H "
        r"\overset{\text{half hub}}{\longleftrightarrow} 5 \overset{e_5}{\longleftrightarrow} 6 \cdots 11, \qquad "
        r"\mathrm H \overset{k_b}{\longleftrightarrow} \mathrm P\ (\text{telescope}).",
        r"\end{gather}",
        rf"\paragraph{{Geometry.}} Each side hub is a single {in_inglese(pm.tubo.materiale.nome)} body: a solid "
        r"flange fixed to the ring and a tube that seats the bearings. "
        rf"The support ring ({in_inglese(pm.anello.materiale.nome)}) holds the telescope, which is slid into "
        r"it, so its axis is the telescope axis $\hat K$, perpendicular to the polar axis; the side hubs are "
        r"fixed to it at two diametrically opposite points. With these materials the nominal geometry gives "
        rf"${num(im.V_geometria * 1e9, 4)}$\,mm$^3$ and ${num(im.m_geometria, 3)}$\,kg, against "
        rf"${num(pm.m, 4)}$\,kg of the inertia data"
        + (r" (the Onshape export, made with steel side hubs, corrected by hand for their actual material: "
           r"placeholder until a new export is available)." if pm.inerzia_segnaposto else r" (Onshape)."),
        *fig,
        r"\paragraph{Stiffness.} The parts are springs in series. The tube and the flange are shafts in torsion, "
        r"$k = G J_p/L$. The ring, loaded by a couple $M$ about $\hat J$ at each attachment point, works as two "
        r"half rings in parallel; along each half the moment is the constant vector $\tfrac12 M\hat J$, which "
        r"splits into an out-of-plane bending moment $\tfrac12 M\cos\varphi$ about the radius and a torque "
        r"$\tfrac12 M\sin\varphi$ about the tangent. The complementary energy (Castigliano) gives",
        r"\begin{gather}",
        r"\theta = \frac{\pi R}{4}\left(\frac{1}{E I_r} + \frac{1}{G J_t}\right) M,\qquad "
        r"k_{ring} = \frac{4}{\pi R}\left(\frac{1}{E I_r} + \frac{1}{G J_t}\right)^{-1},\\",
        r"I_r = \frac{t\,b^3}{12},\qquad J_t = \beta\,b\,t^3,\quad "
        r"\beta = \tfrac13\bigl(1 - 0.63\,\tfrac tb + 0.052\,(\tfrac tb)^5\bigr),\qquad E = 2G(1 + \nu),",
        r"\end{gather}",
        rf"with the mean radius $R = {num(rg.R_medio * 1e3, 4)}$\,mm, the wall thickness $t = {num(rg.spessore * 1e3, 3)}$\,mm, "
        rf"the width $b = {num(pm.anello.larghezza * 1e3, 4)}$\,mm and $\nu = {num(pm.anello.poisson, 2)}$: "
        rf"$E I_r = {num(rg.EI_r, 3)}$\,N\,m$^2$ and $G J_t = {num(rg.GJ_t, 3)}$\,N\,m$^2$. The thin wall twists "
        r"easily, so the ring dominates the series. The ring is taken alone, without the stiffening of the "
        r"telescope tube slid into it (user's choice): its value is a lower bound."
        + (rf" \emph{{In this version the ring stiffness is multiplied by "
           rf"$f_{{ring}} = {num(rg.fattore, 3)}$}} (\texttt{{ring.stiffness\_factor}}): $k_{{ring}} = "
           rf"{num(rg.fattore, 3)} \times {num(rg.k_semianelli, 3)} = {num(rg.k_anello, 3)}$\,N\,m/rad. "
           r"With the ring alone the model has, without the telescope, a structural mode near 175\,Hz that the "
           r"measurement does not show (no peak below about 400\,Hz); with a much stiffer ring the hub is "
           r"practically rigid and the first structural mode without the telescope is set by the shaft between "
           r"hub and motor." if rg.fattore != 1.0 else ""),
        r"\begin{table}[H]\centering\small",
        r"\caption{Hub parts, geometry (mm) and torsional stiffness about the polar axis (N\,m/rad)." + seg
        + r"}\label{tab:hub_parts}",
        r"\begin{tabular}{llr}\toprule",
        r"part & geometry & $k$ \\ \midrule",
        *righe,
        r"\midrule",
        rf"half hub (tube, flange, half ring: $2k_{{ring}}$) & in series & ${num(rg.k_lato, 4)}$ \\",
        rf"whole hub, bearing to bearing & in series & ${num(rg.k_totale, 4)}$ \\",
        rf"telescope attachment $k_b$ & {'placeholder (see text)' if pm.k_b_segnaposto else 'given'} & "
        rf"${num(pm.k_b, 4)}$ \\",
        r"\bottomrule\end{tabular}\end{table}",
        r"\paragraph{Inertia.} From the Onshape export of the three parts (" + esc(pm.fonte_inerzia) + r"), "
        rf"mass ${num(pm.m, 4)}$\,kg, centre of mass at ${num(im.scarto_baricentro * 1e3, 3)}$\,mm from $o_2$, "
        r"tensor about the centre of mass in the axes of N (kg\,m$^2$)",
        r"\begin{gather}",
        r"\mathbf I_{hub,G}|_N = " + tensore + r".",
        r"\end{gather}",
        rf"Transported to $o_2$ and rotated into the frame of body 2, it gives $J_{{hub}} = {num(pm.J_polare, 4)}$"
        r"\,kg\,m$^2$ about the polar axis, at node H, and "
        rf"${num(pm.J_azimut, 4)}$\,kg\,m$^2$ about $\hat z$ at $\theta = 0$"
        + (r", added to $J^{ax}$ because the azimuth load exported from Onshape does not contain the hub "
           r"(\texttt{in\_azimuth\_load: false})." if not pm.nel_load_azimut else
           r", already contained in the azimuth load."),
        *([_testo_tratto_telescopio(pm, par)] if pm.tratto_telescopio is not None else []),
        r"\paragraph{Telescope attachment.} The spring $k_b$ joins H to the telescope node P, which carries "
        r"only the telescope (payload section). Its value is not computed from a geometry: "
        + (r"it is a placeholder, the two-mass estimate from the antiresonance and resonance measured at "
           r"the elevation encoder." if pm.k_b_segnaposto else esc(pm.fonte_k_b) + ".")
        + rf" With the shafts, the reduced polar tree has $k_{{\theta,1}} = {num(ks[a[0]], 4)}$ (encoder--H), "
        rf"$k_{{\theta,2}} = {num(ks[a[1]], 4)}$ (H--motor), $k_{{\theta,3}} = {num(ks[a[2]], 4)}$ (motor--rear end) "
        rf"and $k_{{\theta,4}} = k_b = {num(ks[a[3]], 4)}$\,N\,m/rad (H--P).",
    ])


def _sezione_elementi(par: Parametri, sistema: Sistema) -> str:
    return "\n".join([
        r"\section{Element properties}\label{sec:elements}",
        r"For an annular layer of length $L$, outer and inner diameters $D_o$, $D_i$ ($D_i = 0$ for a solid "
        r"cylinder), density $\rho$ and shear modulus $G$:",
        r"\begin{gather}",
        r"m = \rho\,\frac{\pi}{4}\left(D_o^2 - D_i^2\right)L, \qquad",
        r"J = \frac{\pi\rho L}{32}\left(D_o^4 - D_i^4\right) = \frac{m}{8}\left(D_o^2 + D_i^2\right), \\",
        r"I_t = m\left[\frac{D_o^2 + D_i^2}{16} + \frac{L^2}{12}\right], \qquad",
        r"J_p = \frac{\pi}{32}\left(D_o^4 - D_i^4\right).",
        r"\end{gather}",
        r"$J$ is the axial mass moment of inertia, $I_t$ the transverse one about the layer's own centre of mass "
        r"and $J_p$ the polar second moment of area. Coaxial layers of a segment add their masses and inertias, "
        r"and the torque-carrying layers act in parallel; the segments of an element add their inertias and act "
        r"in series, together with an optional joint compliance:",
        r"\begin{gather}",
        r"k_s = \frac{1}{L_s}\sum_{\ell\ \text{carrying}} G_\ell J_{p,\ell}, \qquad",
        r"\frac{1}{k_e} = \sum_s \frac{1}{k_s} + \frac{1}{k_{joint}}, \qquad",
        r"J_e = \sum_s\sum_\ell J_{s\ell} + \sum J^{x}_e.",
        r"\end{gather}",
        r"$J^{x}_e$ are concentrated extra inertias attached to the element. An element given as CAD data "
        r"enters with its $J_e$ and $k_e$ directly. The transverse inertia of the element about its own centre "
        r"of mass follows from the layers by the parallel-axis theorem along the shaft.",
        tabella_elementi("azimuth", sistema.asse_psi.elementi_originali),
        tabella_elementi("polar", sistema.asse_theta.elementi_originali),
    ])


def _sezione_fem(par: Parametri, sistemi: Mapping[tuple[str, str], Sistema]) -> str:
    s = sistemi[("consistent", "guyan")]
    out = [r"\section{\texorpdfstring{Finite element model of each axis}{Finite element model of each axis}}\label{sec:fem}",
           r"\paragraph{Element matrices.} The two-node torsional element has",
           r"\begin{gather}",
           r"\mathbf K_e = k_e\begin{bmatrix}1&-1\\-1&1\end{bmatrix}, \qquad",
           r"\mathbf M_e^{\ell} = \frac{J_e}{2}\begin{bmatrix}1&0\\0&1\end{bmatrix} \quad\text{(lumped)}.",
           r"\end{gather}",
           r"For the consistent mass, the inertia per unit length $\bar J(x)$ is piecewise constant (one value "
           r"$\bar J_s = J_s/L_s$ per segment) and the shape functions are linear, $N_1 = 1 - x/L_e$, "
           r"$N_2 = x/L_e$:",
           r"\begin{gather}",
           r"(\mathbf M_e^{c})_{ab} = \sum_s \bar J_s\int_{x_s}^{x_s+L_s} N_a N_b\,dx",
           r"= \sum_s \bar J_s\,\frac{L_s}{6}\Bigl(2N_a^0N_b^0 + N_a^0N_b^1 + N_a^1N_b^0 + 2N_a^1N_b^1\Bigr)",
           r"+ \frac{J^{x}_e}{2}\,\mathbf I,",
           r"\end{gather}",
           r"where $N_a^0$, $N_a^1$ are the values of $N_a$ at the ends of segment $s$: the formula is exact because "
           r"the integrand is the product of two linear functions on each segment. For one uniform segment it gives "
           r"$\mathbf M_e^{c} = \frac{J_e}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$. "
           r"We write $a_e$, $b_e$, $c_e$ for the entries $(1,1)$, $(1,2)$, $(2,2)$ of $\mathbf M_e^{c}$.",
           _testo_nodo_P(par, s.asse_theta),
           r"\paragraph{Assembly.} The global matrices are assembled by summing on the shared nodes; the "
           r"concentrated inertias $J^{n}$ (rotors, fork) are added on the diagonal and the generalised force is "
           r"$\mathbf Q = \tau\,\mathbf e_p$ with $p$ the motor node. In symbolic form:"]
    for nome, mod in (("azimuth", s.asse_psi), ("polar", s.asse_theta)):
        K, Ml, Mc = matrici_simboliche_fem(mod)
        out.append(rf"\paragraph{{{nome.capitalize()} axis}} (node order {', '.join(mod.etichette)}):")
        out.append(ridimensiona(r"\mathbf K = " + K))
        out.append(ridimensiona(r"\mathbf M^{\ell} = " + Ml))
        out.append(ridimensiona(r"\mathbf M^{c} = " + Mc))
        out.append(r"Numerically (SI units, consistent mass, shafts and node inertias only):")
        out.append(ridimensiona(r"\mathbf K = " + matrice_numerica(mod.K)))
        out.append(ridimensiona(r"\mathbf M^{c} = " + matrice_numerica(mod.M["consistent"])))
    out += [
        r"\paragraph{Checks.} Both matrices are symmetric, $\mathbf M$ is positive definite, "
        r"$\mathbf K\mathbf 1 = \mathbf 0$ (one rigid mode per axis), and "
        r"$\mathbf 1^{\mathsf T}\mathbf M\mathbf 1 = \sum_e J_e + \sum J^{n}$ for both mass models.",
        r"\paragraph{Lumped versus consistent mass.} For one element,",
        r"\begin{gather}",
        r"\mathbf M_e^{\ell} - \mathbf M_e^{c} = \begin{bmatrix} b_e + \delta_e & -b_e \\ -b_e & b_e - \delta_e"
        r"\end{bmatrix}, \qquad \delta_e = \frac{c_e - a_e}{2}, \qquad",
        r"\det\left(\mathbf M_e^{\ell} - \mathbf M_e^{c}\right) = -\delta_e^2 ,",
        r"\end{gather}",
        r"using $a_e + 2b_e + c_e = J_e$. Hence the difference is positive semidefinite \emph{if and only if} "
        r"$\delta_e = 0$, i.e.\ when the inertia of the element is symmetric about its mid-span (in particular "
        r"for uniform elements, where $v^{\mathsf T}(\mathbf M^{\ell}-\mathbf M^{c})v = \sum_e \frac{J_e}{6}"
        r"(v_i - v_{i+1})^2$). In that case $\mathbf M^{\ell}\succeq\mathbf M^{c}$ and, by the Courant--Fischer "
        r"min-max characterisation,",
        r"\begin{gather}",
        r"\omega_k^2 = \min_{\dim S = k}\ \max_{v\in S}\frac{v^{\mathsf T}\mathbf K v}{v^{\mathsf T}\mathbf M v}"
        r"\quad\Longrightarrow\quad \omega_k^{\ell} \le \omega_k^{c}\quad\forall k,",
        r"\end{gather}",
        r"because a larger mass lowers every Rayleigh quotient. With non-symmetric elements the ordering is "
        r"\emph{not} guaranteed and may fail for individual modes; the code reports it as a warning "
        r"(Section~\ref{sec:assumptions}).",
        r"\paragraph{Damping.} $\mathbf C = \beta\mathbf K + \mathbf C_b$, with $\beta = 2\zeta/\omega_1$ tuned "
        rf"on the first flexible mode of each isolated axis ($\zeta = {num(par.zeta)}$, placeholder; "
        rf"$\beta_\psi = {num(s.beta_psi)}$\,s, $\beta_\theta = {num(s.beta_theta)}$\,s). The mass-proportional "
        r"part is deliberately omitted: it would damp the rigid rotation of a free shaft. $\mathbf C_b$ is the "
        r"viscous friction of the bearings, $\frac{b}{2}\operatorname{diag}(1,1)$ on the two nodes of a "
        r"bearing element.",
    ]
    return "\n".join(out)


def _testo_nodo_P(par: Parametri, mod: fem.ModelloAsse) -> str:
    ep = par.polare.elemento_payload
    if getattr(mod, "ramo", False):
        pm = par.mozzo
        return (rf"\paragraph{{Nodes H and P (hub model).}} Element {ep} has no shaft: its nodes {ep} and "
                rf"{ep + 1} coincide with the hub centre H, which carries the hub inertia $J^n_H$. Elements "
                rf"{pm.elementi_lato[0]} and {pm.elementi_lato[1]} are the two halves of the hub, massless "
                r"springs (their inertia is in $J^n_H$). The telescope node P is an extra node, last in the "
                r"ordering, joined to H by the massless spring $k_b$: it carries no inertia of its own in "
                r"$\mathbf M$, the telescope enters through the payload (Section~\ref{sec:hub}). The matrices are "
                r"therefore those of a tree, not of a chain.")
    if getattr(mod, "telescopio_rigido", False):
        return (rf"\paragraph{{Node P.}} There is no shaft through the telescope: element {ep} is in "
                r"\texttt{mode: telescope}, the rigid telescope joins the shaft ends at nodes "
                rf"{ep} and {ep + 1}, and the two nodes coincide with the telescope node P. The polar model has "
                r"therefore one node less, no element between the two bearings, and the telescope inertia is "
                r"attached to P only through the payload (Section~\ref{subsec:input_frames}).")
    return (r"\paragraph{Node P.} The telescope element is split at $x_P = \xi_p L_3$ into 3a and 3b. In "
            r"geometry mode the segments are cut at $x_P$ and the properties recomputed, so that "
            r"$1/k_{3a} + 1/k_{3b} = 1/k_3$ and $J_{3a} + J_{3b} = J_3$; for CAD elements the element is "
            r"assumed uniform, $k_{3a} = k_3/\xi_p$, $k_{3b} = k_3/(1-\xi_p)$.")


def _gruppi_catena(mod: fem.ModelloAsse) -> list[list[str]]:
    """Etichette degli elementi fra master consecutivi (molle in serie della catena ridotta)."""
    return [[p.etichetta for p, (i, j) in zip(mod.elementi, mod.collegamenti) if a <= i and j <= b]
            for a, b in zip(mod.masters[:-1], mod.masters[1:])]


def _sezione_guyan(par: Parametri, sistemi: Mapping[tuple[str, str], Sistema],
                   confronti: Sequence[ConfrontoRiduzione]) -> str:
    s = sistemi[("consistent", "guyan")]
    sl = sistemi[("consistent", "guyan_lumped")]
    kp, kt = s.rid_psi.k_catena, s.rid_theta.k_catena
    out = [r"\section{Guyan reduction}\label{sec:guyan}",
           r"Partitioning the degrees of freedom into masters (m) and slaves (s),",
           r"\begin{gather}",
           r"\mathbf T = \begin{bmatrix}\mathbf I\\ -\mathbf K_{ss}^{-1}\mathbf K_{sm}\end{bmatrix}, \qquad",
           r"\mathbf M_r = \mathbf T^{\mathsf T}\mathbf M\,\mathbf T, \quad",
           r"\mathbf K_r = \mathbf T^{\mathsf T}\mathbf K\,\mathbf T, \quad",
           r"\mathbf C_r = \mathbf T^{\mathsf T}\mathbf C\,\mathbf T.",
           r"\end{gather}",
           r"For a chain the condensation has a closed form: between consecutive masters $a<b$ the springs are in "
           r"series, and a slave $j$ is interpolated with compliance weights,",
           r"\begin{gather}",
           r"\frac{1}{k_{ab}} = \sum_{i=a}^{b-1}\frac{1}{k_i}, \qquad",
           r"\varphi_j = (1-w_j)\varphi_a + w_j\varphi_b, \qquad",
           r"w_j = \frac{\sum_{i=a}^{j-1}1/k_i}{\sum_{i=a}^{b-1}1/k_i},",
           r"\end{gather}",
           r"slaves beyond the last master follow it rigidly. Thus $\mathbf K_r$ and $\mathbf M_r$ are tridiagonal. "
           r"$\mathbf T$ depends only on the shaft stiffness: the telescope, the fork and the gravitational stiffness "
           r"sit on masters and pass through unchanged. The form \texttt{guyan\_lumped} sums the rows of "
           r"$\mathbf M_r$ on the diagonal: it preserves the total inertia and the rigid mode and yields the pure "
           r"chain \emph{inertia--spring--inertia}.",
           r"\paragraph{Reduced chains.} Azimuth: $\psi_e \overset{k_{\psi,1}}{\longleftrightarrow} \psi_m "
           r"\overset{k_{\psi,2}}{\longleftrightarrow} \psi_L$; polar: $\theta_e \overset{k_{\theta,1}}"
           r"{\longleftrightarrow}\theta_p\overset{k_{\theta,2}}{\longleftrightarrow}\theta_m"
           r"\overset{k_{\theta,3}}{\longleftrightarrow}\theta_r$, with "
           + "; ".join(rf"$k_{{\theta,{n + 1}}}$ collecting elements {', '.join(g)}"
                       for n, g in enumerate(_gruppi_catena(s.asse_theta))) + ".",
           r"\begin{table}[H]\centering\small",
           r"\caption{Reduced parameters (consistent mass). Stiffnesses in N\,m/rad, inertias in kg\,m$^2$; the "
           r"fork inertia is included in $M_{LL}$, the telescope is not.}\label{tab:reduced}",
           r"\begin{tabular}{lrlr}\toprule",
           rf"$k_{{\psi,1}}$ & ${num(kp[0], 5)}$ & $k_{{\theta,1}}$ & ${num(kt[0], 5)}$ \\",
           rf"$k_{{\psi,2}}$ & ${num(kp[1], 5)}$ & $k_{{\theta,2}}$ & ${num(kt[1], 5)}$ \\",
           rf" & & $k_{{\theta,3}}$ & ${num(kt[2], 5)}$ \\ \midrule"]
    for i, r in enumerate(s.rid_psi.ruoli):
        out.append(rf"$M^{{\psi}}_{{{r}{r}}}$ (guyan) & ${num(s.rid_psi.M_r[i, i], 5)}$ & "
                   rf"$M^{{\psi}}_{{{r}{r}}}$ (lumped) & ${num(sl.rid_psi.M_r[i, i], 5)}$ \\")
    for i, r in enumerate(s.rid_theta.ruoli):
        out.append(rf"$M^{{\theta}}_{{{r}{r}}}$ (guyan) & ${num(s.rid_theta.M_r[i, i], 5)}$ & "
                   rf"$M^{{\theta}}_{{{r}{r}}}$ (lumped) & ${num(sl.rid_theta.M_r[i, i], 5)}$ \\")
    out += [r"\bottomrule\end{tabular}\end{table}",
            r"\paragraph{Validation.} Guyan is a Rayleigh--Ritz approximation of the full model, so with the "
            r"\texttt{guyan} mass $\omega_k^{r}\ge\omega_k^{full}$ for every $k$ (verified by the test suite). "
            r"Table~\ref{tab:validation} compares the natural frequencies of the isolated axes "
            rf"(telescope and fork included, $\theta_0 = 0$).",
            tabella_validazione(confronti)]
    return "\n".join(out)


def _testo_payload(sistema: Sistema) -> str:
    """Dati del telescopio come inseriti e valori risultanti rispetto a o_2 in terna 2."""
    pp = sistema.par.payload  # type: ignore[attr-defined]
    pl = sistema.payload
    if not pl.abilitato:
        return r"\paragraph{Current telescope data.} The telescope is disabled (\texttt{enabled: false})."
    out = [r"\paragraph{Current telescope data.}"]
    if pp.sorgente == "manual":
        assert pp.I is not None and pp.m is not None
        if pp.terna == "cad_cg":
            assert pp.r_G_N is not None
            pos = (r"\mathbf r_{G_2}\big|_N = [" + r",\ ".join(num(v * 1e3, 5) for v in pp.r_G_N)
                   + r"]\ \text{mm}")
            nome_I = r"\mathbf I_{p,G_2}\big|_N"
        else:
            assert pp.d is not None
            pos = r"\mathbf d_2 = [" + r",\ ".join(num(v * 1e3, 5) for v in pp.d) + r"]\ \text{mm}"
            nome_I = r"\mathbf I_{p,o_2}\big|_2"
        out.append(r"As entered in \texttt{config.yaml}"
                   + (rf" ({esc(pp.nota)})" if pp.nota else "") + ":")
        out.append(r"\begin{gather}" + nome_I + " = " + matrice_numerica(pp.I, 6)
                   + rf"\ \text{{kg\,m}}^2, \\ m = {num(pp.m, 6)}\ \text{{kg}}, \qquad " + pos
                   + r".\end{gather}")
    else:
        out.append(r"From body B of the CAD data of the rigid model.")
    if pp.sottrai_alberi:
        out.append(rf"The polar pieces of the FEM are subtracted: mass ${num(pl.sottratto['m'], 4)}$\,kg, "
                   rf"$\sum J = {num(pl.sottratto['J'], 4)}$\,kg\,m$^2$, "
                   rf"$\sum(I_t + m y_c^2) = {num(pl.sottratto['trasversale'], 4)}$\,kg\,m$^2$.")
    out.append(rf"Values used by the dynamics: $m_p = {num(pl.m, 5)}$\,kg, eccentricity "
               rf"$\mathbf d_2 = [{num(pl.d_2[0] * 1e3, 4)},\ {num(pl.d_2[1] * 1e3, 4)},\ {num(pl.d_2[2] * 1e3, 4)}]$"
               r"\,mm (components along $\hat I, \hat J, \hat K$),")
    out.append(r"\begin{gather}\mathbf I_{p,o_2}\big|_2 = " + matrice_numerica(pl.I_2, 5)
               + r"\ \text{kg\,m}^2 .\end{gather}")
    return "\n".join(out)


def _sezione_accoppiato(sistema: Sistema) -> str:
    mod = sistema.accoppiato.simbolico
    st = mod.struttura
    sost, _ = cp.sostituzione_catena_lumped(mod)
    sost.update(cp.PRODOTTI_NULLI)
    eq_L = sp.collect(sp.expand(mod.eom[st.j_psi_L].subs(sost)), [mod.ddq[st.j_psi_L]])
    eq_p = sp.collect(sp.expand(mod.eom[st.j_theta_p].subs(sost)), [mod.ddq[st.j_theta_p]])
    Cnz = [(i, j) for i in range(st.n) for j in range(st.n) if mod.C[i, j] != 0]
    C_righe = r", \\ ".join(rf"C_{{{i + 1}{j + 1}}} = {L(mod.C[i, j])}" for i, j in Cnz)
    return "\n".join([
        r"\section{Coupled model}\label{sec:coupled}",
        r"\paragraph{Kinematics.} The fork rotates with $\psi_L$ (rigid fork). Every body of the polar axis has "
        r"angular velocity $\boldsymbol\omega_j = \dot\psi_L\hat k + \dot\theta_j\hat J$. For the "
        r"telescope, Poisson's formula $[\boldsymbol\omega]_\times = \mathbf R^{\mathsf T}\dot{\mathbf R}$ with "
        r"$\mathbf R = \mathbf R_{N2}$ gives, with components along "
        r"$\hat I, \hat J, \hat K$ (derived with SymPy and checked against the closed form),",
        r"\begin{gather}",
        r"\boldsymbol\omega_2\big|_2 = " + L(sp.Matrix(R_P2.T.astype(int)) * mod.omega_p) + ".",
        r"\end{gather}",
        r"\paragraph{Kinetic energy.}",
        r"\begin{gather}",
        r"T = \tfrac12\dot{\mathbf q}_\psi^{\mathsf T}\mathbf M_{r,\psi}\dot{\mathbf q}_\psi + \tfrac12 J^{ax}"
        r"\dot\psi_L^2 + \tfrac12\dot{\mathbf q}_\theta^{\mathsf T}\mathbf M_{r,\theta}\dot{\mathbf q}_\theta"
        r" + \tfrac12\boldsymbol\omega_2|_2^{\mathsf T}\,\mathbf I_{p,o_2}|_2\,\boldsymbol\omega_2|_2 .",
        r"\end{gather}",
        r"\emph{Justification of $J^{ax}$.} Each axisymmetric part of the polar shaft has tensor "
        r"$\operatorname{diag}(I_t, J, I_t)$ in any frame rotated about $\hat J$, and its centre of mass lies on "
        r"$\hat J$ at distance $y_c$ from $o_2$. Since $\hat k\perp\hat J$, its kinetic energy is",
        r"\begin{gather}",
        r"\tfrac12\left(I_t + m\,y_c^2\right)\dot\psi_L^2 + \tfrac12 J\,\dot\theta_j^2,",
        r"\end{gather}",
        r"without mixed terms and without dependence on $\theta$. The axial part is already in "
        r"$\mathbf M_{r,\theta}$; the transverse part adds up to the constant "
        rf"$J^{{ax}} = \sum(I_t + m y_c^2) = {num(sistema.payload.J_ax, 4)}$\,kg\,m$^2$. "
        r"\textbf{The polar shafts therefore produce no gyroscopic coupling}: the dynamic coupling comes only "
        r"from the telescope (anisotropic tensor, off-axis centre of mass), which is why the telescope must sit on "
        r"a master node. The fork inertia is the master inertia $M_{r,\psi,LL}$. The telescope tensor "
        r"$\mathbf I_{p,o_2}|_2$ is obtained from the input as described in Section~\ref{subsec:input_frames}.",
        r"\paragraph{Potential energy and dissipation.}",
        r"\begin{gather}",
        r"V = \tfrac12\mathbf q_\psi^{\mathsf T}\mathbf K_{r,\psi}\mathbf q_\psi + \tfrac12\mathbf q_\theta^{\mathsf T}"
        r"\mathbf K_{r,\theta}\mathbf q_\theta + m_p g\left(d_K\sin\theta_p + d_I\cos\theta_p\right), \\",
        r"\mathcal D = \tfrac12\dot{\mathbf q}_\psi^{\mathsf T}\mathbf C_{r,\psi}\dot{\mathbf q}_\psi + \tfrac12"
        r"\dot{\mathbf q}_\theta^{\mathsf T}\mathbf C_{r,\theta}\dot{\mathbf q}_\theta .",
        r"\end{gather}",
        r"\paragraph{Generalised forces and the reaction torque.} The polar motor applies $+\tau_\theta\hat J$ "
        r"to the rotor (node 7) and $-\tau_\theta\hat J$ to the stator (fork). The virtual rotations are "
        r"$\delta\boldsymbol\phi_{rot} = \delta\psi_L\hat k + \delta\theta_m\hat J$ and "
        r"$\delta\boldsymbol\phi_{stat} = \delta\psi_L\hat k$, so",
        r"\begin{gather}",
        r"\delta W = \tau_\psi\,\hat k\cdot\delta\psi_m\hat k + \tau_\theta\,\hat J\cdot\left("
        r"\delta\boldsymbol\phi_{rot}-\delta\boldsymbol\phi_{stat}\right) = " + L(mod.lavoro_virtuale)
        + r", \qquad \mathbf Q = \tau_\psi\,\mathbf e_{\psi_m} + \tau_\theta\,\mathbf e_{\theta_m}.",
        r"\end{gather}",
        r"The reaction torque on azimuth node 9 exists physically, but it is directed along "
        r"$\hat J\perp\hat k$: on the only degree of freedom of the fork ($\psi_L$, rotation about $\hat k$) "
        r"it does no work and is carried by the azimuth bearings as a bending moment. With polar coordinates "
        r"relative to the fork the exchange is already included and must \emph{not} be added as a force on "
        r"$\psi_L$; it would enter only with a fork bending degree of freedom (out of scope, Q12).",
        r"\paragraph{Equations of motion.} Lagrange's equations, "
        r"$\frac{d}{dt}\frac{\partial\mathcal L}{\partial\dot{\mathbf q}} - \frac{\partial\mathcal L}{\partial"
        r"\mathbf q} + \frac{\partial\mathcal D}{\partial\dot{\mathbf q}} = \mathbf Q$, are derived with SymPy "
        r"and put in the form $\mathbf M(\mathbf q)\ddot{\mathbf q} + \mathbf C(\mathbf q,\dot{\mathbf q})"
        r"\dot{\mathbf q} + \mathbf C_d\dot{\mathbf q} + \partial V/\partial\mathbf q = \mathbf Q$ with "
        r"$\mathbf C$ from the Christoffel symbols ($\dot{\mathbf M} - 2\mathbf C$ is verified to be "
        r"skew-symmetric). The mass matrix is",
        ridimensiona(r"\mathbf M(\mathbf q) = " + L(mod.M)),
        r"and the non-zero entries of $\mathbf C$ are",
        r"\begin{gather*}" + C_righe + r".\end{gather*}",
        r"With the pure-chain reduced mass (\texttt{guyan\_lumped}) and zero products of inertia, the two "
        r"nonlinear equations read",
        ridimensiona(L(eq_L) + " = 0"),
        ridimensiona(L(eq_p) + " = 0"),
        r"with $J_{\psi,p}(\theta) = I_{KK}\sin^2\theta + I_{II}\cos^2\theta$; the equations of "
        r"$\psi_e,\psi_m,\theta_e,\theta_m,\theta_r$ are linear and identical to those of the isolated chains.",
    ])


def _sezione_linearizzazione(contenuto: ContenutoLatex) -> str:
    s = contenuto.sistemi[("consistent", "guyan")]
    ml = contenuto.ml_tabella or contenuto.modelli_lineari[0]
    st = ml.accoppiato.struttura
    torsione = ml.punto.q0[st.j_theta_m] - ml.punto.q0[st.j_theta_p]
    rp, rt = s.rid_psi, s.rid_theta
    fa = tf.fdt_simboliche_catena(rp.n, rp.i_motore, rp.i_encoder, "psi", rp.i_estremo, rp.ruoli)
    fp = tf.fdt_simboliche_albero(rt.n, rt.archi, rt.i_motore, rt.i_encoder, "theta", rt.i_P, rt.ruoli,
                                  con_k_g=True)
    out = [r"\section{Linearisation and transfer functions}\label{sec:linear}",
           r"\paragraph{Operating point.} $\theta_p = \theta_0$, all $\dot\psi = \Omega_0$, all $\dot\theta = 0$, "
           r"$\ddot{\mathbf q} = 0$, $\psi_e = 0$ (cyclic). With $\theta_p$ fixed the residual is affine in the "
           r"other coordinates, so the static torsions and the nominal torques $\boldsymbol\tau_0$ follow from a "
           r"linear system. The holding torque reaches the telescope only through $k_{\theta,2}$: "
           r"$\theta_e = \theta_p$, $\theta_r = \theta_m$ and $\theta_m - \theta_p = \tau_{\theta0}/k_{\theta,2}$; "
           r"the polar encoder, on the opposite side of the telescope with respect to the motor, reads the "
           r"telescope angle exactly in steady state and departs from it only dynamically. "
           rf"At $\theta_0 = {np.rad2deg(ml.theta0):g}^\circ$, $\Omega_0 = {np.rad2deg(ml.Omega0):g}^\circ$/s: "
           rf"$\tau_{{\theta 0}} = {num(ml.punto.tau0[1], 5)}$\,N\,m, "
           rf"$\theta_m - \theta_p = {num(torsione, 4)}$\,rad.",
           r"\paragraph{Linearised matrices.} The Jacobians of the residual at the operating point give "
           r"$\mathbf M_0\delta\ddot{\mathbf q} + \mathbf C_0\delta\dot{\mathbf q} + \mathbf K_0\delta\mathbf q = "
           r"\mathbf E\,\delta\boldsymbol\tau$. With zero products of inertia, the only changes with respect to the "
           r"isolated chains are",
           r"\begin{gather}",
           r"\mathbf C_0[\psi_L,\theta_p] = -\mathbf C_0[\theta_p,\psi_L] = \gamma = (I_{KK}-I_{II})\,\Omega_0"
           r"\sin2\theta_0, \\",
           r"\mathbf K_0[\theta_p,\theta_p] \mathrel{+}= k_c + k_g, \qquad k_c = -(I_{KK}-I_{II})\,\Omega_0^2"
           r"\cos2\theta_0, \qquad k_g = -m_p g\,(d_K\sin\theta_0 + d_I\cos\theta_0).",
           r"\end{gather}",
           r"Numerically, at the operating point above:",
           ridimensiona(r"\mathbf M_0 = " + matrice_numerica(ml.M0, 4)),
           ridimensiona(r"\mathbf C_0 = " + matrice_numerica(ml.C0, 4)),
           ridimensiona(r"\mathbf K_0 = " + matrice_numerica(ml.K0, 4)),
           r"\paragraph{Symbolic transfer functions of the isolated chains} ($\Omega_0 = 0$, pure chain, no "
           r"damping). With $\mathbf P(s) = \operatorname{diag}(J_i)\,s^2 + \mathbf K$ and "
           r"$\Delta(s) = \det\mathbf P(s)$, every transfer function from the motor torque is "
           r"$G_{im}(s) = [\operatorname{adj}\mathbf P(s)]_{im}/\Delta(s)$. Azimuth (torque on $\psi_m$; "
           r"$J_{\psi,L}$ includes the fork, $J^{ax}$ and $J_{\psi,p}(\theta_0)$):",
           r"\begin{gather}",
           r"\mathbf P_\psi(s) = " + L(fa["P"]) + r", \\",
           r"\frac{\Psi_e}{T_\psi} = \frac{" + L(fa["N_encoder"]) + r"}{\Delta_\psi}, \qquad",
           r"\frac{\Psi_m}{T_\psi} = \frac{" + L(fa["N_motore"]) + r"}{\Delta_\psi}, \qquad",
           r"\frac{\Psi_L}{T_\psi} = \frac{" + L(fa["N_carico"]) + r"}{\Delta_\psi},",
           r"\end{gather}",
           ridimensiona(r"\Delta_\psi(s) = " + L(fa["det"])),
           r"The azimuth encoder transfer function has the antiresonance "
           r"$\omega_a = \sqrt{k_{\psi,2}/J_{\psi,L}}$, where motor and encoder stand still while the load "
           r"oscillates on $k_{\psi,2}$; with a light motor and a heavy fork it lies far below the first "
           r"resonance. Polar (torque on $\theta_m$; $J_{\theta,p}$ includes $I_{yy,p}$, and $k_g$ acts on the "
           r"telescope node):",
           r"\begin{gather}",
           r"\mathbf P_\theta(s) = " + L(fp["P"]) + r", \\",
           r"\frac{\Theta_e}{T_\theta} = \frac{" + L(fp["N_encoder"]) + r"}{\Delta_\theta}, \qquad",
           r"\frac{\Theta_p}{T_\theta} = \frac{" + L(fp["N_carico"]) + r"}{\Delta_\theta},",
           r"\end{gather}",
           ridimensiona(r"\frac{\Theta_m}{T_\theta} = \frac{" + L(fp["N_motore"]) + r"}{\Delta_\theta},"),
           r"with $\Delta_\theta(s) = \det\mathbf P_\theta(s)$, of degree 8 in $s$ (not expanded).",
           r"Both encoder transfer functions are \emph{non-collocated}: motor and encoder are on different "
           r"nodes, so the alternation of poles and zeros is not guaranteed. The torque$\to$motor transfer "
           r"functions are collocated and do alternate (test suite). The polar encoder numerator "
           r"$k_{\theta,1}k_{\theta,2}(J_{\theta,r}s^2 + k_{\theta,3})$ shows that its only zeros are those of the "
           r"rear branch beyond the motor.",
           ]
    for chiave, didascalia in (
        ("bode", r"Torque$\to$encoder $2\times2$ matrix of the coupled Guyan-reduced model (comparison), "
                 r"with poles (dotted) and invariant zeros (circles)."),
        ("confronto_azimut", r"Isolated azimuth axis: full against reduced, lumped against consistent mass."),
        ("confronto_polare", r"Isolated polar axis: full against reduced, lumped against consistent mass."),
        ("errore", r"Torque to true pointing ($\psi_L$, $\theta_p$), to the encoders and to the measurement "
                   r"error (encoder minus true angle)."),
        ("modi_azimut", r"Mode shapes of the isolated azimuth axis, full and reduced (expanded with $\mathbf T$)."),
        ("modi_polare", r"Mode shapes of the isolated polar axis, full and reduced (expanded with $\mathbf T$)."),
        ("sweep", r"Natural frequencies of the coupled model as a function of $\theta_0$."),
    ):
        if chiave in contenuto.figure:
            out.append(r"\begin{figure}[H]\centering")
            out.append(rf"\includegraphics[width=\linewidth,height=0.8\textheight,keepaspectratio]"
                       rf"{{{contenuto.figure[chiave]}}}")
            out.append(rf"\caption{{{didascalia} \emph{{{NOTA_REALISMO.capitalize()}.}}}}\end{{figure}}")
    out.append(tabella_poli_zeri(ml))
    return "\n".join(out)


def _sezione_accoppiato_completo(contenuto: ContenutoLatex) -> str:
    """Modello accoppiato con tutti i GdL del FEM (riferimento del troncamento modale), senza Guyan."""
    from . import modal as md

    sm = contenuto.modale.sistema  # type: ignore[union-attr]
    pl = sm.payload
    fp = md.funzioni_payload()
    _, th, dpsi, dth = fp.simboli
    omega = cp.velocita_angolare_forma_chiusa(th, dpsi, dth)
    Cq = (fp.C_simbolica * sp.Matrix([dpsi, dth])).applyfunc(sp.simplify)
    n_psi, n_th = sm.ast_psi.n, sm.ast_theta.n
    et_th = ", ".join(sm.ast_theta.etichetta(i) for i in range(n_th))
    carico = [n for n in contenuto.par.azimut.inerzie_nodali if n.contiene_alberi_polari]
    nota_jax = (rf" In the current data the azimuth load (node {carico[0].nodo}) is a CAD export that already "
                r"contains the polar shafts, so their contribution about $\hat k$, including the true "
                r"$m y_c^2$ terms, is already in the node inertia and $J^{ax}$ is \emph{not} added (the value "
                rf"computed from the polar elements, ${num(pl.J_ax_nel_carico, 4)}$\,kg\,m$^2$, would count "
                r"the shafts twice). The axial inertia of the shafts about $\hat J$ stays in the polar FEM."
                if carico else "")
    return "\n".join([
        r"\section{Coupled model}\label{sec:coupled}",
        r"\paragraph{Coordinates.} The coupled model is written first with \emph{all} FEM degrees of freedom: "
        rf"$\mathbf q_\psi$ ({n_psi} azimuth nodes, absolute rotations about $\hat k$) and $\mathbf q_\theta$ "
        rf"({n_th} polar nodes {et_th}, rotations about $\hat J$ relative to the fork). The fork rotates with "
        r"$\psi_L = \psi_9$ (rigid fork) and the telescope with $\theta_p$, the rotation of node P: "
        r"$\mathbf q_p = [\psi_L,\ \theta_p]^{\mathsf T} = \mathbf E\,\mathbf q$, with $\mathbf E$ selecting the "
        r"two attachment nodes. This is the reference model (exact within the torsional FEM) of which the "
        r"reduced model of Section~\ref{sec:modal} is the projection $\mathbf q = \boldsymbol\Phi\boldsymbol\eta$.",
        r"\paragraph{Kinematics.} Every body of the polar axis has angular velocity $\boldsymbol\omega_j = "
        r"\dot\psi_L\hat k + \dot\theta_j\hat J$. For the telescope, Poisson's formula $[\boldsymbol\omega]_\times "
        r"= \mathbf R^{\mathsf T}\dot{\mathbf R}$ with $\mathbf R = \mathbf R_{N2}$ gives, with components along "
        r"$\hat I, \hat J, \hat K$ (derived with SymPy and checked against the closed form),",
        r"\begin{gather}",
        r"\boldsymbol\omega_2\big|_2 = " + L(sp.Matrix(R_P2.T.astype(int)) * omega) + ".",
        r"\end{gather}",
        r"\paragraph{Energies and dissipation.}",
        r"\begin{gather}",
        r"T = \tfrac12\dot{\mathbf q}_\psi^{\mathsf T}\mathbf M_{\psi}\dot{\mathbf q}_\psi + \tfrac12 J^{ax}"
        r"\dot\psi_L^2 + \tfrac12\dot{\mathbf q}_\theta^{\mathsf T}\mathbf M_{\theta}\dot{\mathbf q}_\theta"
        r" + \tfrac12\boldsymbol\omega_2|_2^{\mathsf T}\,\mathbf I_{p,o_2}|_2\,\boldsymbol\omega_2|_2, \\",
        r"V = \tfrac12\mathbf q_\psi^{\mathsf T}\mathbf K_{\psi}\mathbf q_\psi + \tfrac12\mathbf q_\theta^{\mathsf T}"
        r"\mathbf K_{\theta}\mathbf q_\theta + m_p g\left(d_K\sin\theta_p + d_I\cos\theta_p\right),\qquad "
        r"\mathcal D = \tfrac12\dot{\mathbf q}^{\mathsf T}\mathbf C\,\dot{\mathbf q},",
        r"\end{gather}",
        r"with $\mathbf M$, $\mathbf K$, $\mathbf C = \beta\mathbf K + \mathbf C_b$ the FEM matrices of "
        r"Section~\ref{sec:fem} (the fork inertia is the concentrated inertia of azimuth node 9).",
        r"\emph{Justification of $J^{ax}$.} Each axisymmetric part of the polar shaft has tensor "
        r"$\operatorname{diag}(I_t, J, I_t)$ in any frame rotated about $\hat J$, and its centre of mass lies on "
        r"$\hat J$ at distance $y_c$ from $o_2$. Since $\hat k\perp\hat J$, its kinetic energy is "
        r"$\tfrac12(I_t + m y_c^2)\dot\psi_L^2 + \tfrac12 J\dot\theta_j^2$, without mixed terms and without "
        r"dependence on $\theta$: the axial part is in $\mathbf M_\theta$, the transverse part is the constant "
        rf"$J^{{ax}} = \sum(I_t + m y_c^2)$ (here ${num(pl.J_ax, 4)}$\,kg\,m$^2$). \textbf{{The polar shafts "
        r"therefore produce no gyroscopic coupling}: the dynamic coupling comes only from the telescope." + nota_jax,
        r"\paragraph{Generalised forces and the reaction torque.} The polar motor applies $+\tau_\theta\hat J$ "
        r"to the rotor and $-\tau_\theta\hat J$ to the stator (fork); with $\delta\boldsymbol\phi_{rot} = "
        r"\delta\psi_L\hat k + \delta\theta_m\hat J$ and $\delta\boldsymbol\phi_{stat} = \delta\psi_L\hat k$,",
        r"\begin{gather}",
        r"\delta W = \tau_\psi\,\delta\psi_m + \tau_\theta\,\hat J\cdot\left(\delta\boldsymbol\phi_{rot} - "
        r"\delta\boldsymbol\phi_{stat}\right) = \tau_\psi\,\delta\psi_m + \tau_\theta\,\delta\theta_m,\qquad "
        r"\mathbf Q = \mathbf B\boldsymbol\tau = \tau_\psi\,\mathbf e_{\psi_m} + \tau_\theta\,\mathbf e_{\theta_m}.",
        r"\end{gather}",
        r"The reaction torque on the fork is along $\hat J\perp\hat k$: it does no work on $\psi_L$ and is carried "
        r"by the azimuth bearings as a bending moment (rigid fork, Q12).",
        r"\paragraph{Telescope terms.} Only the telescope depends nonlinearly on the state. Its Lagrangian "
        r"$T_p - V_g$, a function of $\mathbf q_p$ alone, is differentiated by SymPy, which also checks the "
        r"matrix form against Lagrange's equations:",
        riduci_se_serve(r"\mathbf M_p(\theta_p) = " + L(fp.M_simbolica) + ","),
        riduci_se_serve(r"\mathbf C_p\dot{\mathbf q}_p = " + L(Cq) + ","),
        r"\begin{gather}\mathbf g_p = " + L(fp.g_simbolica) + r".\end{gather}",
        r"$J_{\psi,p}(\theta) = M_{p,11}$ is the telescope inertia about the azimuth axis. Lagrange's equations "
        r"of the whole system are then",
        r"\begin{gather}",
        r"\bigl[\mathbf M_s + \mathbf E^{\mathsf T}\mathbf M_p(\theta_p)\mathbf E\bigr]\ddot{\mathbf q} + "
        r"\mathbf C\dot{\mathbf q} + \mathbf K\mathbf q + \mathbf E^{\mathsf T}\bigl(\mathbf C_p\dot{\mathbf q}_p"
        r" + \mathbf g_p\bigr) = \mathbf B\boldsymbol\tau,\qquad \mathbf M_s = \operatorname{blkdiag}(\mathbf M_\psi"
        r" + J^{ax}\mathbf e_L\mathbf e_L^{\mathsf T},\ \mathbf M_\theta),\label{eq:full_eom}",
        r"\end{gather}",
        r"and the state dependence is only through $\theta_p$ and the velocities $\dot{\mathbf q}_p$.",
    ])


def _sezione_linearizzazione_modale(contenuto: ContenutoLatex) -> str:
    """Punto di lavoro e linearizzazione del payload (senza le catene di Guyan)."""
    sm = contenuto.modale.sistema  # type: ignore[union-attr]
    pl = sm.payload
    par = contenuto.par
    righe = []
    for th in par.theta0:
        for Om in par.omega0:
            ml = sm.lineare(th, Om)
            righe.append(rf"${np.rad2deg(th):g}$ & ${np.rad2deg(Om):g}$ & ${num(ml.punto.tau0[1], 5)}$ & "
                         rf"${num(float(pl.k_g(th)), 4)}$ & ${num(pl.k_c(th, Om), 3)}$ & "
                         rf"${num(pl.gamma(th, Om), 3)}$ & ${num(float(pl.M_12(th)), 3)}$ \\")
    return "\n".join([
        r"\section{Operating point and linearisation}\label{sec:linear}",
        r"\paragraph{Operating point.} A steady tracking motion: $\theta_p = \theta_0$, rigid azimuth rotation at "
        r"$\dot\psi = \Omega_0$, all $\dot\theta = 0$, $\ddot{\mathbf q} = \mathbf 0$, $\psi_e = 0$ (cyclic). With "
        r"$\theta_p$ fixed the residual of \eqref{eq:full_eom} is affine in the remaining coordinates, so the "
        r"static torsions and the nominal torques $\boldsymbol\tau_0$ follow from a linear system "
        r"(Section~\ref{sec:workflow}).",
        r"\paragraph{Linearised telescope terms.} At the operating point the Jacobians of the telescope terms "
        r"with respect to $\ddot{\mathbf q}_p$, $\dot{\mathbf q}_p$ and $\mathbf q_p$ are",
        r"\begin{gather}",
        r"\mathbf M_{p0} = \mathbf M_p(\theta_0),\qquad \mathbf G_{p0} = \begin{bmatrix}0 & \gamma\\ -\gamma & 0"
        r"\end{bmatrix},\qquad \mathbf K_{p0} = \begin{bmatrix}0 & 0\\ 0 & k_c + k_g\end{bmatrix},\\",
        r"\gamma = \Omega_0\,J_{\psi,p}'(\theta_0),\qquad k_c = -\tfrac12\Omega_0^2\,J_{\psi,p}''(\theta_0),\qquad "
        r"k_g = -m_p g\,(d_K\sin\theta_0 + d_I\cos\theta_0),",
        r"\end{gather}",
        r"($\gamma$ gyroscopic, $k_c$ centrifugal, $k_g$ gravitational stiffness; $\mathbf G_{p0}$ stays purely "
        r"skew-symmetric also with products of inertia, which enter only through $M_{p,12}$). With zero products, "
        r"$\gamma = (I_{KK}-I_{II})\,\Omega_0\sin2\theta_0$ and $k_c = -(I_{KK}-I_{II})\,\Omega_0^2\cos2\theta_0$. "
        r"The linear model of the full system is",
        r"\begin{gather}",
        r"\bigl[\mathbf M_s + \mathbf E^{\mathsf T}\mathbf M_{p0}\mathbf E\bigr]\delta\ddot{\mathbf q} + "
        r"\bigl[\mathbf C + \mathbf E^{\mathsf T}\mathbf G_{p0}\mathbf E\bigr]\delta\dot{\mathbf q} + "
        r"\bigl[\mathbf K + \mathbf E^{\mathsf T}\mathbf K_{p0}\mathbf E\bigr]\delta\mathbf q = \mathbf B\,"
        r"\delta\boldsymbol\tau;",
        r"\end{gather}",
        r"its projection on the modal basis is the reduced linear model of Section~\ref{sec:modal}. "
        r"Table~\ref{tab:op} lists the operating points ($k_g < 0$: the unbalanced telescope is statically "
        r"unstable without control).",
        r"\begin{table}[H]\centering\small",
        r"\caption{Operating points: nominal polar torque and linearised telescope terms (N\,m, N\,m/rad, "
        rf"N\,m\,s/rad, kg\,m$^2$). \emph{{{NOTA_REALISMO.capitalize()}.}}}}\label{{tab:op}}",
        r"\begin{tabular}{rrrrrrr}\toprule",
        r"$\theta_0$ [deg] & $\Omega_0$ [deg/s] & $\tau_{\theta0}$ & $k_g$ & $k_c$ & $\gamma$ & $M_{p,12}$ \\ "
        r"\midrule",
        *righe,
        r"\bottomrule\end{tabular}\end{table}",
    ])


def _sezione_concentrato(contenuto: ContenutoLatex) -> str:
    from . import concentrato_output as co

    s_mod = contenuto.modale.sistema if contenuto.modale is not None else None  # type: ignore[union-attr]
    return co.latex(contenuto.concentrato, contenuto.par.segnaposto, s_mod=s_mod,  # type: ignore[arg-type]
                    extra=contenuto.concentrato_extra)


def _sezione_modale(contenuto: ContenutoLatex) -> str:
    from . import modal_output as mo

    guyan = _con_guyan(contenuto)
    out = [mo.latex(contenuto.modale, dict(contenuto.figure), contenuto.par.segnaposto,  # type: ignore[arg-type]
                    guyan=guyan)]
    if guyan:
        return "\n".join(out)
    # senza Guyan i risultati del modello lineare stanno qui (figure e poli/zeri del modale)
    out.append(r"\subsection{Results of the reduced linear model}")
    for chiave, didascalia in (
        ("errore", r"Torque to true pointing ($\psi_L$, $\theta_p$), to the encoders and to the measurement "
                   r"error (encoder minus true angle), modal model."),
        ("sweep", r"Natural frequencies of the coupled modal model as a function of $\theta_0$."),
        ("modi_modale_azimut", r"First modes of the full FEM of the azimuth axis (basis of the truncation), "
                               r"retained or discarded."),
        ("modi_modale_polare", r"First modes of the full FEM of the polar axis (basis of the truncation), "
                               r"retained or discarded."),
    ):
        if chiave in contenuto.figure:
            out += [r"\begin{figure}[H]\centering",
                    rf"\includegraphics[width=\linewidth,height=0.8\textheight,keepaspectratio]"
                    rf"{{{contenuto.figure[chiave]}}}",
                    rf"\caption{{{didascalia} \emph{{{NOTA_REALISMO.capitalize()}.}}}}\end{{figure}}"]
    if contenuto.ml_principale is not None:
        out.append(tabella_poli_zeri(contenuto.ml_principale))
    return "\n".join(out)


def _sezione_ipotesi(avvisi: Sequence[Avviso]) -> str:
    out = [r"\section{Assumptions and placeholders}\label{sec:assumptions}",
           r"\begin{longtable}{p{0.31\linewidth}p{0.12\linewidth}p{0.49\linewidth}}",
           r"\toprule code & status & description \\ \midrule\endhead"]
    stato = {"segnaposto": "placeholder", "incoerenza": "inconsistency", "informazione": "note"}
    for a in avvisi:
        out.append(rf"\texttt{{{esc(a.codice)}}} & {stato[a.severita]} & {esc(a.messaggio)} \\")
    out.append(r"\bottomrule\end{longtable}")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Autovalori e autovettori: soluzione analitica e numerica
# ---------------------------------------------------------------------------


def _tabella_barra(par: Parametri) -> tuple[str, float]:
    """Barra uniforme libera-libera: soluzione analitica contro FEM concentrato e consistente."""
    mat = par.materiali.get("acciaio") or next(iter(par.materiali.values()))
    rho, G, Lb, D = mat.rho, mat.G, 0.3, 0.02
    Jp = np.pi * D**4 / 32.0
    c = float(np.sqrt(G / rho))
    righe = [r"\begin{table}[H]\centering\small",
             rf"\caption{{Uniform free--free steel bar ($L = {Lb:g}$\,m, $D = {D * 1e3:g}$\,mm, "
             rf"$\rho = {mat.rho:g}$\,kg/m$^3$, $G = {mat.G / 1e9:g}$\,GPa, $c = {c:.5g}$\,m/s): exact "
             r"frequencies $f_n = nc/2L$ against the torsional FEM with $N$ equal elements. The errors follow the "
             r"closed forms $\mp(n\pi/N)^2/24$.}\label{tab:bar}",
             r"\begin{tabular}{rrrrrrr}\toprule",
             r"$N$ & $n$ & $f_n$ [Hz] & lumped [Hz] & error & consistent [Hz] & error \\ \midrule"]
    pic = 0.0
    for N in (4, 16, 64):
        h = Lb / N
        k, J = G * Jp / h, rho * Jp * h
        f = {m: fem.frequenze_proprie(*fem.catena_uniforme([k] * N, [J] * N, m))  # type: ignore[arg-type]
             for m in ("lumped", "consistent")}
        for n in (1, 2, 3):
            fe = n * c / (2 * Lb)
            fl, fc = f["lumped"][n], f["consistent"][n]
            pic = max(pic, abs(fl / fe - 1), abs(fc / fe - 1))
            righe.append(rf"{N if n == 1 else ''} & {n} & ${num(fe, 6)}$ & ${num(fl, 6)}$ & "
                         rf"${num((fl / fe - 1) * 100, 3)}\,\%$ & ${num(fc, 6)}$ & "
                         rf"${num((fc / fe - 1) * 100, 3)}\,\%$ \\")
        righe.append(r"\midrule" if N != 64 else "")
    righe.append(r"\bottomrule\end{tabular}\end{table}")
    return "\n".join(righe), pic


def _catena_azimut(sistemi: Mapping[tuple[str, str], Sistema]) -> str:
    """Catena ridotta di azimut a tre inerzie: polinomio caratteristico e confronto analitico/numerico."""
    import scipy.linalg as sla

    s = sistemi.get(("consistent", "guyan_lumped"))
    if s is None:
        return ""
    M, K, _ = s.isolato_ridotto("azimut", 0.0)
    if np.max(np.abs(M - np.diag(np.diag(M)))) > 1e-12 * np.max(np.abs(M)):
        return ""
    Je, Jm, JL, k1, k2, lam = sp.symbols("J_e J_m J_L k_1 k_2 lambda", positive=True)
    Ms = sp.diag(Je, Jm, JL)
    Ks = sp.Matrix([[k1, -k1, 0], [-k1, k1 + k2, -k2], [0, -k2, k2]])
    p = sp.factor(sp.det(Ks - lam * Ms))
    a, b, c = sp.Poly(sp.cancel(p / lam), lam).all_coeffs()
    valori = {Je: M[0, 0], Jm: M[1, 1], JL: M[2, 2], k1: -K[0, 1], k2: -K[1, 2]}
    an, bn, cn = (float(v.subs(valori)) for v in (a, b, c))
    disc = bn**2 - 4 * an * cn
    lam_an = np.sort(np.concatenate([[0.0], np.array([(-bn - np.sqrt(disc)) / (2 * an),
                                                      (-bn + np.sqrt(disc)) / (2 * an)])]))
    w2, Phi = sla.eigh(K, M)
    # modi analitici: dalle prime due righe di (K - lambda M) phi = 0
    righe = []
    for i, l_a in enumerate(lam_an):
        ph = np.array([1.0, 0.0, 0.0])
        ph[1] = 1.0 - l_a * M[0, 0] / (-K[0, 1])
        ph[2] = ((K[1, 1] - l_a * M[1, 1]) * ph[1] + K[0, 1] * ph[0]) / (-K[1, 2])
        pn = Phi[:, i]
        mac = float((ph @ M @ pn) ** 2 / ((ph @ M @ ph) * (pn @ M @ pn)))
        fa = np.sqrt(max(l_a, 0.0)) / (2 * np.pi)
        fn = np.sqrt(max(w2[i], 0.0)) / (2 * np.pi)
        righe.append(rf"{i + 1} & ${num(fa, 8)}$ & ${num(fn, 8)}$ & "
                     rf"${num(abs(fa - fn) / max(fa, 1.0), 2)}$ & ${mac:.12f}$ \\")
    return "\n".join([
        r"\paragraph{A three-inertia chain in closed form.} The reduced azimuth chain (encoder, motor, load, "
        r"\texttt{guyan\_lumped}) has $\mathbf M = \operatorname{diag}(J_e, J_m, J_L)$ and the chain stiffness "
        r"matrix with springs $k_1$, $k_2$. Its characteristic polynomial, expanded with SymPy, is",
        r"\begin{gather}",
        r"\det(\mathbf K - \lambda\mathbf M) = " + sp.latex(p) + r",",
        r"\end{gather}",
        r"so one root is the rigid mode, $\lambda = 0$, and the other two solve a quadratic, "
        r"$a\lambda^2 + b\lambda + c = 0$, with",
        r"\begin{gather}",
        r"a = " + sp.latex(a) + r", \qquad b = " + sp.latex(b) + r", \qquad c = " + sp.latex(c) + r",",
        r"\end{gather}",
        r"hence $\lambda_{\pm} = \bigl(-b \pm \sqrt{b^2 - 4ac}\bigr)/(2a)$ and $f = \sqrt\lambda/2\pi$. The mode "
        r"shape follows from the first two rows of $(\mathbf K - \lambda\mathbf M)\boldsymbol\varphi = \mathbf 0$ "
        r"by forward substitution,",
        r"\begin{gather}",
        r"\varphi_1 = 1,\qquad \varphi_2 = 1 - \frac{\lambda J_e}{k_1}\,,\qquad "
        r"\varphi_3 = \frac{(k_1 + k_2 - \lambda J_m)\,\varphi_2 - k_1\varphi_1}{k_2}\,.",
        r"\end{gather}",
        r"Table~\ref{tab:chain3} compares this closed-form solution, evaluated with the current numbers, with "
        r"the numerical solver of the next subsection; the last column is the modal assurance criterion "
        r"$\mathrm{MAC} = (\boldsymbol\varphi_a^{\mathsf T}\mathbf M\boldsymbol\varphi_n)^2/"
        r"[(\boldsymbol\varphi_a^{\mathsf T}\mathbf M\boldsymbol\varphi_a)(\boldsymbol\varphi_n^{\mathsf T}"
        r"\mathbf M\boldsymbol\varphi_n)]$.",
        r"\begin{table}[H]\centering\small",
        r"\caption{Reduced azimuth chain at $\theta_0 = 0$: closed form against \texttt{scipy.linalg.eigh}.}"
        r"\label{tab:chain3}",
        r"\begin{tabular}{rrrrr}\toprule",
        r"mode & $f$ closed form [Hz] & $f$ \texttt{eigh} [Hz] & relative difference & MAC \\ \midrule",
        *righe,
        r"\bottomrule\end{tabular}\end{table}",
    ])


def _sezione_autovalori(par: Parametri, sistemi: Mapping[tuple[str, str], Sistema],
                        guyan: bool = True) -> str:
    tab_barra, _ = _tabella_barra(par)
    ridotti = (
        r"\item \textbf{Reduced axes.} The reduced chain is solved with the same \texttt{eigh}, on "
        r"$\mathbf K_r$ and $\mathbf M_r$. The Guyan transformation is computed by "
        r"\texttt{numpy.linalg.solve} on $\mathbf K_{ss}\mathbf X = -\mathbf K_{sm}$ (LU with partial pivoting, "
        r"\texttt{dgesv}), never by forming $\mathbf K_{ss}^{-1}$. Because $\mathbf T$ is a Ritz basis, the "
        r"reduced eigenvalues are Rayleigh--Ritz approximations and bound the exact ones from above; the "
        r"Craig--Bampton basis adds the fixed-interface modes of the slave partition, obtained with "
        r"\texttt{eigh} on $(\mathbf K_{ss}, \mathbf M_{ss})$." if guyan else
        r"\item \textbf{Modal basis} (\texttt{modal.modi\_con\_rigido}). The same \texttt{eigh} on the full FEM "
        r"with the telescope frozen, $(\mathbf K, \bar{\mathbf M})$; the column of the null eigenvalue is checked to "
        r"be proportional to $\mathbf 1$ and replaced by $\mathbf 1$, the others are kept mass-normalised. The "
        r"retained columns are a Ritz basis that contains the exact modes, so their frequencies are exact "
        r"(Section~\ref{sec:modal}).")
    verifica = (r"the Rayleigh--Ritz bound of Guyan, " if guyan else
                r"the exactness of the retained modal frequencies, ")
    return "\n".join([
        r"\section{Natural frequencies and mode shapes}\label{sec:eigen}",
        r"\subsection{The eigenvalue problem}",
        r"Free, undamped torsional vibration of an axis obeys $\mathbf M\ddot{\mathbf q} + \mathbf K\mathbf q = "
        r"\mathbf 0$. Looking for harmonic solutions $\mathbf q = \boldsymbol\varphi\,e^{j\omega t}$ gives the "
        r"generalised symmetric eigenvalue problem",
        r"\begin{gather}",
        r"\left(\mathbf K - \omega^2\mathbf M\right)\boldsymbol\varphi = \mathbf 0 \quad\Longleftrightarrow\quad "
        r"\mathbf K\boldsymbol\varphi = \lambda\,\mathbf M\boldsymbol\varphi,\qquad \lambda = \omega^2,",
        r"\end{gather}",
        r"whose non-trivial solutions are the roots of the characteristic polynomial "
        r"$\det(\mathbf K - \lambda\mathbf M) = 0$, of degree $n$ in $\lambda$. Because $\mathbf K$ is symmetric "
        r"positive semidefinite and $\mathbf M$ symmetric positive definite, all eigenvalues are real and "
        r"non-negative, and the eigenvectors can be chosen mass-orthonormal:",
        r"\begin{gather}",
        r"\boldsymbol\varphi_i^{\mathsf T}\mathbf M\boldsymbol\varphi_j = \delta_{ij},\qquad "
        r"\boldsymbol\varphi_i^{\mathsf T}\mathbf K\boldsymbol\varphi_j = \omega_i^2\delta_{ij},\qquad "
        r"\lambda_i = \frac{\boldsymbol\varphi_i^{\mathsf T}\mathbf K\boldsymbol\varphi_i}"
        r"{\boldsymbol\varphi_i^{\mathsf T}\mathbf M\boldsymbol\varphi_i}\ (\text{Rayleigh quotient}).",
        r"\end{gather}",
        r"Since $\mathbf K\mathbf 1 = \mathbf 0$, every free axis has exactly one rigid mode, "
        r"$\lambda_1 = 0$, $\boldsymbol\varphi_1\propto\mathbf 1$; the others are the torsional modes, with "
        r"$f_i = \sqrt{\lambda_i}/2\pi$. When the gravitational stiffness $k_g$ is added at the telescope node "
        r"(isolated polar axis) the matrix $\mathbf K + k_g\mathbf e_P\mathbf e_P^{\mathsf T}$ may have a "
        r"\emph{negative} eigenvalue, reported by the code with a negative sign and interpreted as the static "
        r"instability of the unbalanced telescope rather than as a mode. With the stiffness-proportional damping "
        r"$\mathbf C = \beta\mathbf K$ the modal equations decouple, "
        r"$\ddot\eta_i + 2\zeta_i\omega_i\dot\eta_i + \omega_i^2\eta_i = \boldsymbol\varphi_i^{\mathsf T}"
        r"\mathbf Q$ with $\zeta_i = \beta\omega_i/2$, which is why $\beta$ is tuned on the first flexible mode.",
        r"\subsection{Analytical solutions}",
        r"\paragraph{Uniform shaft.} A uniform free--free shaft obeys the wave equation "
        r"$\partial_t^2\vartheta = c^2\partial_x^2\vartheta$, $c = \sqrt{G/\rho}$, with free ends "
        r"$\partial_x\vartheta = 0$ (no torque), so",
        r"\begin{gather}",
        r"\varphi_n(x) = \cos\frac{n\pi x}{L},\qquad \omega_n = \frac{n\pi}{L}\sqrt{\frac{G}{\rho}},\qquad "
        r"n = 0, 1, 2,\dots",
        r"\end{gather}",
        r"independent of the diameter. The FEM of $N$ equal elements has the exact discrete spectrum "
        r"$\theta_n = n\pi/N$ and",
        r"\begin{gather}",
        r"\omega_n^{\ell} = \frac{2c}{h}\sin\frac{\theta_n}{2},\qquad "
        r"\omega_n^{c} = \frac{c}{h}\sqrt{\frac{6\,(1 - \cos\theta_n)}{2 + \cos\theta_n}},\qquad h = L/N,",
        r"\end{gather}",
        r"for lumped and consistent mass respectively; expanding, $\omega_n^{\ell,c}\approx\omega_n\,[1\mp"
        r"(n\pi/N)^2/24]$, so the two models bound the exact frequency from below and above and converge with "
        r"second order.",
        tab_barra,
        _catena_azimut(sistemi) if guyan else "",
        r"\subsection{Numerical solution}",
        r"The code never solves the characteristic polynomial: finding the roots of a polynomial whose "
        r"coefficients span many decades (here $\lambda$ ranges from about $10^{0}$ to $10^{10}$\,s$^{-2}$) "
        r"is ill-conditioned. It uses the following linear-algebra tools instead.",
        r"\begin{enumerate}",
        r"\item \textbf{Torsional modes of an axis} (\texttt{fem.modi}). \texttt{scipy.linalg.eigh(K, M)} calls "
        r"the LAPACK driver \texttt{dsygvd} for the generalised symmetric-definite problem: Cholesky "
        r"factorisation $\mathbf M = \mathbf L\mathbf L^{\mathsf T}$, reduction to the standard symmetric "
        r"problem $\mathbf L^{-1}\mathbf K\mathbf L^{-\mathsf T}\mathbf y = \lambda\mathbf y$, Householder "
        r"tridiagonalisation, a divide-and-conquer eigensolver, and back-transformation "
        r"$\boldsymbol\varphi = \mathbf L^{-\mathsf T}\mathbf y$, which returns mass-orthonormal modes. The "
        r"cost is $O(n^3)$ on matrices of size at most 12. The rigid mode comes out as $|\lambda| \sim "
        r"10^{-16}\lambda_{max}$ and is set to zero below the threshold $10^{-12}\lambda_{max}$; "
        r"frequencies are returned in Hz.",
        ridotti,
        r"\item \textbf{Coupled linearised system} (\texttt{ModelloLineare.poli}). With gyroscopic and "
        r"centrifugal terms the problem is no longer self-adjoint, so the poles are the eigenvalues of the "
        r"state matrix $\mathbf A = \begin{bmatrix}\mathbf 0 & \mathbf I\\ -\mathbf M_0^{-1}\mathbf K_0 & "
        r"-\mathbf M_0^{-1}\mathbf C_0\end{bmatrix}$ computed with \texttt{numpy.linalg.eigvals} (LAPACK "
        r"\texttt{dgeev}: balancing, Hessenberg reduction, shifted QR). A pair $\lambda = -\zeta\omega_n \pm "
        r"j\omega_n\sqrt{1-\zeta^2}$ gives $f_n = |\lambda|/2\pi$ and $\zeta = -\mathrm{Re}\,\lambda/|\lambda|$; "
        r"\texttt{tf.classifica} labels each root as null, real stable, real unstable or complex. The "
        r"azimuth integrator gives a double null pole, a defective eigenvalue that \texttt{dgeev} returns "
        r"scattered at $\sim10^{-9}|\lambda|_{max}$; the code treats anything below that threshold as zero.",
        r"\item \textbf{Zeros.} The invariant zeros of each SISO channel are the finite generalised "
        r"eigenvalues of the Rosenbrock pencil "
        r"$\left(\begin{bmatrix}\mathbf A & \mathbf b_j\\ \mathbf c_i & 0\end{bmatrix}, "
        r"\begin{bmatrix}\mathbf I & \mathbf 0\\ \mathbf 0 & 0\end{bmatrix}\right)$, computed with "
        r"\texttt{scipy.linalg.eig} (QZ algorithm, \texttt{dggev}) in homogeneous form, discarding the infinite "
        r"eigenvalues ($\beta \approx 0$).",
        r"\item \textbf{Frequency response.} $\mathbf G(j\omega) = \mathbf C_q\,\mathbf P(j\omega)^{-1}"
        r"\mathbf E$, with $\mathbf P(s) = \mathbf M_0 s^2 + \mathbf C_0 s + \mathbf K_0$, is evaluated "
        r"directly by a batched \texttt{numpy.linalg.solve} over the whole frequency grid, not through the "
        r"polynomial coefficients, and cross-checked against $\mathbf C(j\omega\mathbf I - \mathbf A)^{-1}"
        r"\mathbf B$.",
        r"\end{enumerate}",
        r"\paragraph{Verification.} The test suite checks the numerical solution against the analytical one: "
        r"the shaft convergence above (lumped from below, consistent from above), the "
        r"identity $\mathbf v^{\mathsf T}(\mathbf M^{\ell} - \mathbf M^{c})\mathbf v = \sum_e J_e(v_i - "
        r"v_{i+1})^2/6$, " + verifica + r"the equality of the poles with the roots of "
        r"$\det\mathbf P(s)$, the alternation of poles and zeros on the collocated channels, and the "
        r"equality of the zeros of the collocated channel with the natural frequencies of the structure with "
        r"the motor node clamped.",
    ])


# ---------------------------------------------------------------------------
# Funzionamento del modello e simulazione
# ---------------------------------------------------------------------------


def _tikz_flusso(modale: bool = False, guyan: bool = True) -> str:
    box = r"draw, rounded corners=3pt, fill=#1, text width=3.5cm, align=center, minimum height=1.55cm, font=\footnotesize"
    return "\n".join([
        r"\begin{tikzpicture}[box/.style={" + box.replace("#1", "blue!8") + r"},"
        r"box2/.style={" + box.replace("#1", "green!10") + r"},"
        r"box3/.style={" + box.replace("#1", "orange!12") + r"}]",
        r"\node[box] (cfg) at (0,0) {\textbf{Input}\\ \texttt{config.yaml}\\ \scriptsize PyYAML, validation};",
        r"\node[box] (s1) at (4.4,0) {\textbf{1. Elements}\\ $m, J, I_t, k$ per element\\ "
        r"\scriptsize closed formulas, NumPy};",
        r"\node[box] (s2) at (8.8,0) {\textbf{2. Axis FEM}\\ $\mathbf K,\mathbf M$, $\mathbf K\boldsymbol\varphi = "
        r"\lambda\mathbf M\boldsymbol\varphi$\\ \scriptsize SciPy \texttt{eigh}};",
        (r"\node[box] (s3) at (13.2,0) {\textbf{3. Modal truncation}\\ $\boldsymbol\Phi = [\mathbf 1,"
         r"\boldsymbol\phi_1,\dots]$, Hankel, validation\\ \scriptsize SciPy \texttt{eigh}"
         + (r"; Guyan for comparison};" if guyan else r"};")
         if modale else
         r"\node[box] (s3) at (13.2,0) {\textbf{3. Guyan reduction}\\ $\mathbf T = [\mathbf I;\,-\mathbf K_{ss}^{-1}"
         r"\mathbf K_{sm}]$\\ \scriptsize NumPy \texttt{solve}};"),
        r"\node[box2] (s4) at (13.2,-2.7) {\textbf{4. Coupled model}\\ Lagrange: $T, V, \mathcal D\to"
        r"\mathbf M(\mathbf q)\ddot{\mathbf q} + \dots$\\ \scriptsize SymPy, \texttt{lambdify}};"
        if guyan else
        r"\node[box2] (s4) at (13.2,-2.7) {\textbf{4. Coupled model}\\ telescope Lagrangian on $\mathbf E_\eta"
        r"\boldsymbol\eta$\\ \scriptsize SymPy, \texttt{lambdify}};",
        r"\node[box2] (s5) at (8.8,-2.7) {\textbf{5. Linear analysis}\\ equilibrium, $\mathbf M_0,\mathbf C_0,"
        r"\mathbf K_0$, FRF, poles, zeros\\ \scriptsize "
        + (r"SymPy Jacobians" if guyan else r"closed form (7.2)") + r", \texttt{eigvals}, QZ};",
        r"\node[box3] (s6) at (4.4,-2.7) {\textbf{6. Time simulation}\\ $\dot{\mathbf x} = \mathbf f(\mathbf x,"
        r"\boldsymbol\tau)$, S1--S3\\ \scriptsize SciPy \texttt{solve\_ivp}, DOP853};",
        r"\node[box3] (out) at (0,-2.7) {\textbf{Output}\\ tables, figures, report, PDF\\ \scriptsize "
        r"Matplotlib, \texttt{pdflatex}};",
        r"\draw[->, thick] (cfg) -- (s1); \draw[->, thick] (s1) -- (s2); \draw[->, thick] (s2) -- (s3);",
        r"\draw[->, thick] (s3) -- (s4); \draw[->, thick] (s4) -- (s5); \draw[->, thick] (s5) -- (s6);",
        r"\draw[->, thick] (s6) -- (out);",
        r"\end{tikzpicture}",
    ])


def _scarto_riferimento(contenuto: "ContenutoLatex") -> float:
    """Scarto massimo (relativo al picco) di S1 fra il modello principale e quello a tutti i GdL."""
    from .linearize import USCITE

    rf_, rm_ = contenuto.s1_riferimento, contenuto.simulazioni[0]
    n = min(rf_.t.size, rm_.t.size)  # type: ignore[union-attr]
    out = 0.0
    for u in ("psi_e", "psi_m", "psi_L"):
        i = USCITE.index(u)
        picco = float(np.max(np.abs(rf_.y_nl[:n, i])))  # type: ignore[union-attr]
        out = max(out, float(np.max(np.abs(rm_.y_nl[:n, i] - rf_.y_nl[:n, i]))) / picco)  # type: ignore[union-attr]
    return out


def _blocco_equazioni(contenuto: "ContenutoLatex", guyan: bool) -> list[str]:
    """Sottosezioni "Equations" e "Operating point" del flusso di calcolo.

    ``guyan``: modello accoppiato SymPy a 7 GdL delle catene di Guyan; altrimenti il
    modello modale (eq. ``eq:modal_eom``), integrato numericamente.
    """
    if not guyan:
        return _blocco_equazioni_modale(contenuto)
    st = contenuto.sistemi[("consistent", "guyan")].accoppiato.struttura
    return [
            r"\subsection{Equations of the coupled model}",
            r"The two reduced chains are coupled in a single Lagrangian system with "
            rf"$n = {st.n}$ coordinates, $\mathbf q = [\psi_e, \psi_m, \psi_L, \theta_e, \theta_p, \theta_m, "
            r"\theta_r]^{\mathsf T}$ for the Guyan chains. The kinetic energy $T$, the potential energy $V$ and the dissipation function "
            r"$\mathcal D$ of Section~\ref{sec:coupled} are written symbolically, the reduced matrices entering "
            r"with the tridiagonal structure of a chain; SymPy differentiates "
            r"$\mathcal L = T - V$ and returns",
            r"\begin{gather}",
            r"\mathbf M(\mathbf q)\,\ddot{\mathbf q} + \bigl[\mathbf C(\mathbf q,\dot{\mathbf q}) + \mathbf C_d\bigr]"
            r"\dot{\mathbf q} + \frac{\partial V}{\partial\mathbf q} = \mathbf E\,\boldsymbol\tau,\qquad "
            r"\boldsymbol\tau = [\tau_\psi,\ \tau_\theta]^{\mathsf T},",
            r"\end{gather}",
            r"where $\mathbf C(\mathbf q,\dot{\mathbf q})$ comes from the Christoffel symbols of $\mathbf M$ "
            r"(centrifugal and gyroscopic terms), $\mathbf C_d = \partial^2\mathcal D/\partial\dot{\mathbf q}^2$ is "
            r"the structural damping and $\mathbf E$ places $\tau_\psi$ on $\psi_m$ and $\tau_\theta$ on $\theta_m$. "
            r"All state dependence is through $\theta_p$ (the telescope node) and the velocities: "
            r"$\mathbf M$ differs from $\operatorname{diag}(\mathbf M_{r,\psi},\mathbf M_{r,\theta})$ only in the "
            r"entries $(\psi_L,\psi_L)$, $(\psi_L,\theta_p)$ and $(\theta_p,\theta_p)$. SymPy also checks that "
            r"$\dot{\mathbf M} - 2\mathbf C$ is skew-symmetric and that the matrix form equals the Lagrange "
            r"residual. The expressions are compiled with \texttt{lambdify} (NumPy, common subexpression "
            r"elimination), so that evaluating $\mathbf M(\mathbf q)$ and the force vector costs tens of "
            r"microseconds. In first-order form, with $\mathbf x = [\mathbf q,\ \dot{\mathbf q}]$,",
            r"\begin{gather}",
            r"\dot{\mathbf x} = \mathbf f(\mathbf x,\boldsymbol\tau) = \begin{bmatrix}\dot{\mathbf q}\\ "
            r"\mathbf M(\mathbf q)^{-1}\Bigl[\mathbf E\boldsymbol\tau - \bigl(\mathbf C + \mathbf C_d\bigr)"
            r"\dot{\mathbf q} - \partial V/\partial\mathbf q\Bigr]\end{bmatrix},",
            r"\end{gather}",
            r"where the inverse is never formed: the code solves $\mathbf M\ddot{\mathbf q} = \dots$ with "
            r"\texttt{numpy.linalg.solve} at every right-hand-side call.",
            r"\subsection{Operating point and linearisation}",
            r"An operating point is a steady tracking motion: $\theta_p = \theta_0$, every $\dot\psi_i = \Omega_0$, "
            r"every $\dot\theta_i = 0$, $\ddot{\mathbf q} = \mathbf 0$ and $\psi_e = 0$ (cyclic). The residual "
            r"$\mathbf r(\mathbf q) = \mathbf M\ddot{\mathbf q} + (\mathbf C + \mathbf C_d)\dot{\mathbf q} + "
            r"\partial V/\partial\mathbf q$ is affine in the free coordinates once $\theta_p$ is fixed, so the "
            r"static torsions $\mathbf q_0$ and nominal torques $\boldsymbol\tau_0$ solve the linear system "
            r"$\mathbf K_0[:,\text{free}]\,\mathbf q_{free} - \mathbf E\boldsymbol\tau_0 = -\mathbf r_0$, "
            r"by \texttt{numpy.linalg.solve}; the code then verifies that the residual is below $10^{-8}$ of its "
            r"scale. Jacobians of $\mathbf r$ with respect to $\ddot{\mathbf q}$, $\dot{\mathbf q}$, $\mathbf q$ "
            r"(symbolic, evaluated at the operating point) give $\mathbf M_0$, $\mathbf C_0$, $\mathbf K_0$ and "
            r"the linear model",
            r"\begin{gather}",
            r"\mathbf M_0\,\delta\ddot{\mathbf q} + \mathbf C_0\,\delta\dot{\mathbf q} + \mathbf K_0\,"
            r"\delta\mathbf q = \mathbf E\,\delta\boldsymbol\tau,\qquad \delta\dot{\mathbf x} = \mathbf A\,"
            r"\delta\mathbf x + \mathbf B\,\delta\boldsymbol\tau,\qquad \mathbf y = \mathbf C_q\,\delta\mathbf q,",
            r"\end{gather}",
            r"with outputs $\psi_e, \theta_e$ (encoders), $\psi_m, \theta_m$ (motors), $\psi_L, \theta_p$ (true "
            r"pointing), $\theta_r$ and the measurement errors $\psi_e - \psi_L$, $\theta_e - \theta_p$. "
            r"Frequency responses, poles and zeros of this model are computed as in Section~\ref{sec:eigen}.",
    ]


def _blocco_equazioni_modale(contenuto: "ContenutoLatex") -> list[str]:
    r = contenuto.modale
    n = r.stati // 2 if r is not None else 0  # type: ignore[union-attr]
    return [
        r"\subsection{Equations of the reduced model}",
        rf"The default model has $n = {n}$ modal coordinates $\boldsymbol\eta$ (Section~\ref{{sec:modal}}) and "
        r"obeys \eqref{eq:modal_eom}. Its mass matrix",
        r"\begin{gather}",
        r"\mathbf M(\boldsymbol\eta) = \tilde{\mathbf M} + \mathbf E_\eta^{\mathsf T}\bigl[\mathbf M_p(\theta_p) - "
        r"\bar{\mathbf M}_p\bigr]\mathbf E_\eta,\qquad \theta_p = \mathbf c_P^{\mathsf T}\boldsymbol\eta_\theta,",
        r"\end{gather}",
        r"depends on the state only through $\theta_p$, and the force vector is $\tilde{\mathbf D}\dot{\boldsymbol"
        r"\eta} + \tilde{\mathbf K}\boldsymbol\eta + \mathbf E_\eta^{\mathsf T}(\mathbf C_p\dot{\mathbf q}_p + "
        r"\mathbf g_p)$. The $2\times2$ payload blocks $\mathbf M_p$, $\mathbf C_p\dot{\mathbf q}_p$ and "
        r"$\mathbf g_p$ of Section~\ref{sec:coupled} are derived once by SymPy (which also checks them against "
        r"Lagrange's equations of $T_p - V_g$) and compiled with \texttt{lambdify}; the structural matrices are "
        r"constant projections $\mathbf T^{\mathsf T}\mathbf M\mathbf T$, $\mathbf T^{\mathsf T}\mathbf K\mathbf T$, "
        r"$\mathbf T^{\mathsf T}\mathbf C\mathbf T$ of the full FEM. In first-order form, with $\mathbf x = "
        r"[\boldsymbol\eta,\ \dot{\boldsymbol\eta}]$,",
        r"\begin{gather}",
        r"\dot{\mathbf x} = \mathbf f(\mathbf x,\boldsymbol\tau) = \begin{bmatrix}\dot{\boldsymbol\eta}\\ "
        r"\mathbf M(\boldsymbol\eta)^{-1}\Bigl[\boldsymbol\Phi^{\mathsf T}\mathbf B\boldsymbol\tau - \tilde{\mathbf D}"
        r"\dot{\boldsymbol\eta} - \tilde{\mathbf K}\boldsymbol\eta - \mathbf E_\eta^{\mathsf T}(\mathbf C_p\dot{\mathbf q}_p"
        r" + \mathbf g_p)\Bigr]\end{bmatrix},",
        r"\end{gather}",
        r"where the inverse is never formed: the code solves $\mathbf M\ddot{\boldsymbol\eta} = \dots$ with "
        r"\texttt{numpy.linalg.solve} at every right-hand-side call. The same code, with the identity as basis, "
        r"is the full-DOF reference model used for validation.",
        r"\subsection{Operating point and linearisation}",
        r"An operating point is a steady tracking motion: $\theta_p = \theta_0$, rigid azimuth rotation at "
        r"$\Omega_0$ ($\dot\eta_{\psi0} = \Omega_0$, all other velocities zero), $\ddot{\boldsymbol\eta} = \mathbf 0$ "
        r"and $\psi_e = 0$ (cyclic). With $\theta_p$ fixed the residual is affine in $\boldsymbol\eta$, so the "
        r"static deformation and the nominal torques solve the square linear system",
        r"\begin{gather}",
        r"\begin{bmatrix}\tilde{\mathbf K} & -\boldsymbol\Phi^{\mathsf T}\mathbf B\\ \hat{\mathbf c}_P^{\mathsf T} & "
        r"\mathbf 0\\ \hat{\mathbf c}_e^{\mathsf T} & \mathbf 0\end{bmatrix}\begin{bmatrix}\boldsymbol\eta_0\\ "
        r"\boldsymbol\tau_0\end{bmatrix} = \begin{bmatrix}-\mathbf r_0\\ \theta_0\\ 0\end{bmatrix},",
        r"\end{gather}",
        r"with $\mathbf r_0$ the payload and bearing forces of the steady motion, solved by "
        r"\texttt{numpy.linalg.solve}; the code verifies that the residual is below $10^{-8}$ of its scale. "
        r"The linear matrices are the closed form of Section~\ref{sec:modal}, with the payload Jacobians of "
        r"Section~\ref{sec:linear}; the test suite checks them against finite differences of the nonlinear "
        r"residual. The linear model is",
        r"\begin{gather}",
        r"\mathbf M_0\,\delta\ddot{\boldsymbol\eta} + \mathbf C_0\,\delta\dot{\boldsymbol\eta} + \mathbf K_0\,"
        r"\delta\boldsymbol\eta = \boldsymbol\Phi^{\mathsf T}\mathbf B\,\delta\boldsymbol\tau,\qquad \delta\dot{\mathbf x}"
        r" = \mathbf A\,\delta\mathbf x + \mathbf B\,\delta\boldsymbol\tau,\qquad \mathbf y = \mathbf C_q\,"
        r"\delta\boldsymbol\eta,",
        r"\end{gather}",
        r"with outputs $\psi_e, \theta_e$ (encoders), $\psi_m, \theta_m$ (motors), $\psi_L, \theta_p$ (true "
        r"pointing), $\theta_r$ and the measurement errors $\psi_e - \psi_L$, $\theta_e - \theta_p$, all "
        r"obtained from the rows of the modal basis at the corresponding nodes. Frequency responses, poles and "
        r"zeros of this model are computed as in Section~\ref{sec:eigen}.",
    ]


def _sezione_funzionamento(contenuto: "ContenutoLatex") -> str:
    import inspect
    import platform

    import scipy

    from . import simulate as sm

    modale = contenuto.modale is not None
    guyan = _con_guyan(contenuto)
    firma = inspect.signature(sm.S1).parameters
    ampiezza, durata = float(firma["ampiezza"].default), float(firma["durata"].default)
    righe_sim = []
    for r in contenuto.simulazioni:
        if r.nome == "S3":
            righe_sim.append(rf"S3 & {num(np.degrees(r.theta0), 3)} & {num(np.degrees(r.Omega0), 3)} & "
                             rf"conservation, $\mathbf C = \mathbf 0$, $\boldsymbol\tau = \mathbf 0$ & "
                             rf"$|\Delta E|/E = {num(r.derive['energia'], 2)}$, "
                             rf"$|\Delta p_\psi|/p_\psi = {num(r.derive['momento_psi'], 2)}$ \\")
        else:
            e = r.errori
            righe_sim.append(rf"{r.nome} & {num(np.degrees(r.theta0), 3)} & {num(np.degrees(r.Omega0), 3)} & "
                             rf"doublet on $\tau_\psi$ & $e(\psi_e) = {num(e['psi_e'], 2)}$, "
                             rf"$e(\psi_L) = {num(e['psi_L'], 2)}$ \\")
    fig = []
    for chiave, didascalia in (
        ("S1", r"S1: nonlinear (solid) and linearised (dashed) response to the $\tau_\psi$ doublet at "
               r"$\Omega_0 = 0$."),
        ("S2", r"S2: the same doublet while tracking at $\Omega_{max}$; the cross response on "
               r"$\theta_e$ and $\theta_p$ is the gyroscopic coupling."),
        ("S3", r"S3: drift of $E = T + V$ and $p_\psi$ in the conservative case."),
    ):
        if chiave in contenuto.figure:
            fig += [r"\begin{figure}[H]\centering",
                    rf"\includegraphics[width=0.9\linewidth,height=0.45\textheight,keepaspectratio]"
                    rf"{{{contenuto.figure[chiave]}}}",
                    rf"\caption{{{didascalia} \emph{{{NOTA_REALISMO.capitalize()}.}}}}\end{{figure}}"]
    return "\n".join([
        r"\section{How the model works and how it is simulated}\label{sec:workflow}",
        r"\subsection{Computational workflow}",
        r"Everything is computed from \texttt{config.yaml} by \texttt{run\_all.py}; Figure~\ref{fig:flow} shows "
        r"the chain of stages, each with its module and numerical tool, and Table~\ref{tab:stages} lists the "
        r"equations solved at each stage.",
        r"\begin{figure}[H]\centering",
        r"\resizebox{\linewidth}{!}{" + _tikz_flusso(modale, guyan) + "}",
        r"\caption{Computational workflow. Blue: structural stages (Sections~\ref{sec:elements}--"
        + (r"\ref{sec:guyan}" if guyan else r"\ref{sec:eigen} and~\ref{sec:modal}")
        + r"); green: coupled model and its analysis; orange: simulation and output.}"
        r"\label{fig:flow}\end{figure}",
        r"\begin{table}[H]\centering\small",
        r"\caption{Equations solved at each stage and numerical tools.}\label{tab:stages}",
        r"\begin{tabular}{p{0.04\linewidth}p{0.4\linewidth}p{0.46\linewidth}}\toprule",
        r"& equations & tool \\ \midrule",
        r"1 & $m, J, I_t, J_p$ of every layer; $k_s$ in parallel, $1/k_e$ in series (Section~\ref{sec:elements}) "
        r"& closed formulas in NumPy \\",
        r"2 & assembly of $\mathbf K,\mathbf M^{\ell},\mathbf M^{c}$; $\mathbf K\boldsymbol\varphi = "
        r"\lambda\mathbf M\boldsymbol\varphi$ (Section~\ref{sec:eigen}) & \texttt{scipy.linalg.eigh} "
        r"(LAPACK \texttt{dsygvd}) \\",
        (r"3 & $(\mathbf K - \omega^2\bar{\mathbf M})\boldsymbol\phi = \mathbf 0$ of the full FEM; Hankel "
         r"norms; validation against the full model (Section~\ref{sec:modal})"
         + (r"; Guyan $\mathbf K_{ss}\mathbf X = -\mathbf K_{sm}$ for comparison" if guyan else "")
         + r" & \texttt{scipy.linalg.eigh}, \texttt{numpy.linalg.solve} \\" if modale else
         r"3 & $\mathbf K_{ss}\mathbf X = -\mathbf K_{sm}$; $\mathbf M_r = \mathbf T^{\mathsf T}\mathbf M\mathbf T$ "
         r"& \texttt{numpy.linalg.solve} (LU) \\"),
        (r"4 & Lagrange's equations of the 7-DOF system (this section) & SymPy for the derivation; "
         r"\texttt{lambdify} to NumPy functions \\" if guyan else
         r"4 & Lagrange's equations of the telescope, \eqref{eq:modal_eom} (this section) & SymPy for the "
         r"derivation; \texttt{lambdify} to NumPy functions \\"),
        r"5 & steady state $\mathbf r(\mathbf q_0) = \mathbf E\boldsymbol\tau_0$; Jacobians; "
        r"$\mathbf G(j\omega)$; poles; zeros & NumPy \texttt{solve}, \texttt{eigvals}; SymPy "
        r"\texttt{jacobian}; \texttt{scipy.linalg.eig} (QZ) \\",
        r"6 & $\dot{\mathbf x} = \mathbf f(\mathbf x,\boldsymbol\tau)$ and its linearisation & "
        r"\texttt{scipy.integrate.solve\_ivp}, method DOP853 \\",
        r"7 & tables (CSV), figures, \texttt{.npz}/\texttt{.mat}/JSON exports, this document & Matplotlib, "
        r"\texttt{scipy.io.savemat}, \texttt{pdflatex} \\ \bottomrule",
        r"\end{tabular}\end{table}",
        *_blocco_equazioni(contenuto, guyan),
        r"\subsection{Time-domain simulation}",
        (r"Three scenarios are run on the default reduced model, the modal truncation of Section~\ref{sec:modal} "
         r"with a single $\theta_{ref}$ (\texttt{reduction.modal.theta\_ref\_deg}), for which \eqref{eq:modal_eom} "
         r"is exact for any $\theta_p$ "
         if contenuto.riduzione == "modal" and modale else r"Three scenarios are run on the coupled reduced model ")
        + rf"(Table~\ref{{tab:sim}}). In S1 and S2 the system starts from the equilibrium of the operating point "
        r"with the nominal torque applied as a feedforward, so that the perturbation is only the doublet",
        r"\begin{gather}",
        rf"\tau_\psi(t) = \tau_{{\psi0}} + a\,\bigl[\,\chi_{{[t_i,\,t_i+T)}}(t) - "
        rf"\chi_{{[t_i+T,\,t_i+2T)}}(t)\,\bigr],\qquad a = {num(ampiezza, 3)}\ \text{{N\,m}},\ "
        rf"T = {num(durata * 1e3, 3)}\ \text{{ms}},\ t_i = 0.1\,T,",
        r"\end{gather}",
        r"of zero mean so that it excites the dynamics without shifting the operating point. The nonlinear "
        r"system $\dot{\mathbf x} = \mathbf f(\mathbf x,\boldsymbol\tau)$ and the linear one "
        r"$\delta\dot{\mathbf x} = \mathbf A\delta\mathbf x + \mathbf B\delta\boldsymbol\tau$ are integrated "
        r"with the same solver and compared as $\mathbf q - \mathbf q_0$ against $\delta\mathbf q$ (the "
        r"rigid drift $\Omega_0 t$ of the $\psi$'s is removed). In S3 all damping and torques are removed, "
        r"the state starts from a rigid motion ($\dot\psi_i = \Omega_{max}$, $\dot\theta_i$ small, unstrained "
        r"chains) and the conserved quantities $E = T + V$ and $p_\psi = \sum_i\partial\mathcal L/\partial"
        r"\dot\psi_i$ are monitored.",
        r"\paragraph{Integrator.} \texttt{scipy.integrate.solve\_ivp} with method \texttt{DOP853}: an explicit "
        r"Runge--Kutta scheme of order 8 (Dormand--Prince 8(5,3)) with adaptive step and embedded error control, "
        r"relative tolerance $10^{-10}$ and absolute tolerance $10^{-12}$. Since the doublet is "
        r"discontinuous, the integration is split at $t_i$, $t_i+T$, $t_i+2T$, so that the adaptive solver never "
        r"steps across a jump (which would degrade its order). The step is limited by the highest natural "
        r"frequency of the reduced chains, a few kHz, which sets the cost of a run (tens of seconds); the "
        r"explicit method is adequate because the damping is light and the horizon short. The agreement is "
        r"measured by $e(y) = \max_t|y_{nl}(t) - y_{lin}(t)|/\max_t|y_{lin}(t)|$ for each output.",
        *([r"\paragraph{Full-DOF reference.} S1 is also run on the coupled model with all FEM degrees of freedom. "
           r"Its modes reach hundreds of kHz, so the explicit DOP853 step is limited by stability and "
           r"the reference is integrated with relative tolerance $10^{-8}$ and absolute tolerance $10^{-10}$; the "
           rf"largest difference from the modal model on $\psi_e$, $\psi_m$, $\psi_L$ is "
           rf"${num(_scarto_riferimento(contenuto), 2)}$ of the peak response."]
          if contenuto.s1_riferimento is not None and contenuto.simulazioni else []),
        r"\begin{table}[H]\centering\small",
        r"\caption{Simulation scenarios and results (relative errors of the nonlinear against the linearised "
        r"response, and conservation drifts).}\label{tab:sim}",
        r"\begin{tabular}{llllp{0.4\linewidth}}\toprule",
        r"id & $\theta_0$ [deg] & $\Omega_0$ [deg/s] & forcing & result \\ \midrule",
        *righe_sim,
        r"\bottomrule\end{tabular}\end{table}",
        r"In S1 the polar channels respond only through the products of inertia, so their linear response is "
        r"of the same order as the quadratic term $\tfrac12 J_{\psi,p}'\dot\psi^2$ that the linearisation "
        r"drops: their relative error is therefore not an index of the quality of the linearisation, and "
        r"only the azimuth channels are used for it. Their amplitude (about $10^{-12}$\,rad) is also "
        r"comparable with the absolute tolerance of the integrator, so the isolated spikes visible on the "
        r"nonlinear polar curves are integration noise, not dynamics.",
        *fig,
        r"\subsection{Software}",
        r"\begin{table}[H]\centering\small",
        r"\caption{Software used to produce this document.}\label{tab:software}",
        r"\begin{tabular}{ll}\toprule",
        rf"Python & {platform.python_version()} \\",
        rf"NumPy & {np.__version__} \\",
        rf"SciPy & {scipy.__version__} (LAPACK through \texttt{{scipy.linalg}}) \\",
        rf"SymPy & {sp.__version__} \\",
        r"Matplotlib, PyYAML, pytest & figures, input parsing, test suite \\",
        r"\texttt{pdflatex} & this document (TikZ for the drawings) \\ \bottomrule",
        r"\end{tabular}\end{table}",
        r"All results are reproducible with \texttt{python run\_all.py}; the symbolic derivation is "
        r"deterministic and the numerical tolerances above are fixed in \texttt{simulate.py}.",
    ])


def genera(contenuto: ContenutoLatex, cartella: str | Path) -> Path:
    """Scrive ``multidof_fem_model.tex`` nella cartella indicata."""
    cartella = Path(cartella)
    cartella.mkdir(parents=True, exist_ok=True)
    par = contenuto.par
    s = contenuto.sistemi[("consistent", "guyan")]
    corpo = [
        PREAMBOLO,
        (r"\title{Torsional multi-DOF model of the alt-azimuth mount:\\ FEM, modal truncation and coupling"
         if contenuto.modale is not None else
         r"\title{Torsional multi-DOF model of the alt-azimuth mount:\\ FEM, Guyan reduction and coupling")
        + ((r"\\ \large hub model: telescope hanging from the hub of the polar axis"
            + (r",\\ tube section clamped in the ring rigid with the hub}" if par.mozzo.tratto_telescopio is not None
               else "}")) if par.mozzo is not None else "}"),
        r"\author{Generated by \texttt{gimbal\_mdof.latex\_report}}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
    ]
    if par.segnaposto:
        corpo.append(r"\avviso{All shaft geometries, materials, joint stiffnesses, rotor and fork inertias and "
                     r"damping values are \textbf{placeholders}. Every frequency, transfer function and figure "
                     r"in this document is a \emph{" + NOTA_REALISMO + r"} and must not be read as a "
                     r"characterisation of the real mount (Section~\ref{sec:assumptions}).}")
    corpo += [r"\tableofcontents",
              _sezione_sistema(par, contenuto.modale is not None, _con_guyan(contenuto)),
              _sezione_terne(par, s),
              _sezione_geometria(par, s),
              _sezione_elementi(par, s),
              *([_sezione_mozzo(par, s, contenuto.figure)] if par.mozzo is not None else []),
              _sezione_fem(par, contenuto.sistemi),
              *([_sezione_guyan(par, contenuto.sistemi, contenuto.confronti)] if _con_guyan(contenuto) else []),
              _sezione_autovalori(par, contenuto.sistemi, _con_guyan(contenuto)),
              *([_sezione_accoppiato(s), _sezione_linearizzazione(contenuto)] if _con_guyan(contenuto) else
                [_sezione_accoppiato_completo(contenuto), _sezione_linearizzazione_modale(contenuto)]),
              *([_sezione_modale(contenuto)] if contenuto.modale is not None else []),
              *([_sezione_concentrato(contenuto)] if contenuto.concentrato is not None else []),
              *([contenuto.sperimentale_tex] if contenuto.sperimentale_tex else []),
              _sezione_funzionamento(contenuto),
              _sezione_ipotesi(contenuto.avvisi),
              r"\end{document}"]
    percorso = cartella / f"{NOME_FILE}.tex"
    percorso.write_text("\n".join(corpo), encoding="utf-8")
    return percorso


@dataclass(frozen=True)
class EsitoCompilazione:
    compilato: bool
    pdf: Path | None
    errori: tuple[str, ...]
    riferimenti_irrisolti: bool
    messaggio: str


def compila(tex: Path) -> EsitoCompilazione:
    """``latexmk -pdf`` se disponibile, altrimenti ``pdflatex`` tre volte; controlla il log."""
    tex = Path(tex)
    cartella = tex.parent
    if shutil.which("latexmk"):
        comandi = [["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", tex.name]]
    elif shutil.which("pdflatex"):
        comandi = [["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name]] * 3
    else:
        return EsitoCompilazione(False, None, (), False,
                                 "LaTeX non installato: lasciato il .tex senza compilarlo")
    for c in comandi:
        subprocess.run(c, cwd=cartella, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=600, check=False)
    log = (cartella / f"{tex.stem}.log").read_text(encoding="utf-8", errors="replace")
    errori = tuple(m.group(0) for m in re.finditer(r"^! .*$", log, flags=re.M))
    irrisolti = bool(re.search(r"There were undefined references|Reference .* undefined|"
                               r"Label\(s\) may have changed", log))
    pdf = cartella / f"{tex.stem}.pdf"
    ok = pdf.exists() and not errori
    return EsitoCompilazione(ok, pdf if pdf.exists() else None, errori, irrisolti,
                             "compilato" if ok else "errori di compilazione")
