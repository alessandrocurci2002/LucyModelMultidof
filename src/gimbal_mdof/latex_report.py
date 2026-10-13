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
from .params import MM, NOTA_REALISMO, R_P2, Avviso, ParametriAsse, Parametri

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
    if isinstance(espr, sp.Basic):
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
    """Contributi per strato e totali di ogni elemento (longtable)."""
    r: list[str] = [
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
                     rf"{st.d_in / MM:g} & {st.lunghezza / MM:g} & {esc(st.materiale)} & "
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
            r.append(rf"{c.asse} & {c.modello_massa} & {esc(c.forma_massa)} & {k + 1} & "
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


PREAMBOLO = r"""\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[english]{babel}
\usepackage[margin=2.2cm]{geometry}
\usepackage{amsmath,amssymb,mathtools,bm}
\usepackage{booktabs,longtable,array}
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


def _sezione_sistema(par: Parametri) -> str:
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
\item Guyan reduction of each axis to a short chain on master nodes (\texttt{reduction.py},
      Section~\ref{sec:guyan});
\item coupling of the two reduced chains through the kinetic energy of the whole system
      (\texttt{coupled.py}, Section~\ref{sec:coupled}).
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

\paragraph{Reduced coordinates.}
\begin{gather}
\mathbf q = [\psi_e,\ \psi_m,\ \psi_L,\ \theta_e,\ \theta_p,\ \theta_m,\ \theta_r]^{\mathsf T},
\end{gather}
with masters $\{1, 5, 9\}$ (azimuth) and $\{1, \mathrm P, 7, 11\}$ (polar).
"""


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
        r"\begin{tabular}{p{0.15\linewidth}p{0.17\linewidth}p{0.34\linewidth}p{0.22\linewidth}}\toprule",
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
        r"\end{tabular}\end{table}",
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
         r"(\texttt{payload}) is attached to P." if telescopio else
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
    fp = tf.fdt_simboliche_catena(rt.n, rt.i_motore, rt.i_encoder, "theta", rt.i_P, rt.ruoli,
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
        ("bode", r"Torque$\to$encoder $2\times2$ matrix of the coupled reduced model, with poles (dotted) and "
                 r"invariant zeros (circles)."),
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


def _sezione_ipotesi(avvisi: Sequence[Avviso]) -> str:
    out = [r"\section{Assumptions and placeholders}\label{sec:assumptions}",
           r"\begin{longtable}{p{0.27\linewidth}p{0.12\linewidth}p{0.53\linewidth}}",
           r"\toprule code & status & description \\ \midrule\endhead"]
    stato = {"segnaposto": "placeholder", "incoerenza": "inconsistency", "informazione": "note"}
    for a in avvisi:
        out.append(rf"\texttt{{{esc(a.codice)}}} & {stato[a.severita]} & {esc(a.messaggio)} \\")
    out.append(r"\bottomrule\end{longtable}")
    return "\n".join(out)


def genera(contenuto: ContenutoLatex, cartella: str | Path) -> Path:
    """Scrive ``multidof_fem_model.tex`` nella cartella indicata."""
    cartella = Path(cartella)
    cartella.mkdir(parents=True, exist_ok=True)
    par = contenuto.par
    s = contenuto.sistemi[("consistent", "guyan")]
    corpo = [
        PREAMBOLO,
        r"\title{Torsional multi-DOF model of the alt-azimuth mount:\\ FEM, Guyan reduction and coupling}",
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
              _sezione_sistema(par),
              _sezione_terne(par, s),
              _sezione_geometria(par, s),
              _sezione_elementi(par, s),
              _sezione_fem(par, contenuto.sistemi),
              _sezione_guyan(par, contenuto.sistemi, contenuto.confronti),
              _sezione_accoppiato(s),
              _sezione_linearizzazione(contenuto),
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
