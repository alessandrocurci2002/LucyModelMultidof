"""Sezione LaTeX del modello a parametri concentrati (3 + 4 GdL).

Contenuto: definizione degli elementi equivalenti (inerzia somma, rigidezza in serie),
parametri numerici, equazioni di Newton-Eulero non lineari per ogni GdL (forma compatta
e forma esplicita con le funzioni trigonometriche), verifica contro Lagrange. Le
sottosezioni sulle funzioni di trasferimento e sul simulatore non lineare si aggiungono
con gli argomenti opzionali di :func:`latex`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import sympy as sp

from . import concentrato as cc
from . import coupled as cp
from . import fem
from .params import NOTA_REALISMO, in_inglese


@dataclass(frozen=True)
class FormaCompatta:
    """Equazioni di Newton-Eulero in forma compatta, verificate contro la forma espansa."""

    J_psi_p: sp.Expr
    M_12: sp.Expr
    M_gJ: sp.Expr
    residuo: sp.Matrix


def forma_compatta(variante: str = "catena") -> FormaCompatta:
    """Riscrive le equazioni dei nodi L e P con :math:`J_{\\psi,p}(\\theta_p)`, :math:`M_{12}(\\theta_p)`
    e le loro derivate, e verifica che coincidano con quelle di :func:`concentrato.newton_eulero`."""
    ne = cc.newton_eulero(variante)
    q, dq, ddq = ne.q, ne.dq, ne.ddq
    iL, ip = ne.nomi.index("psi_L"), ne.i_theta_p
    th = q[ip]
    Jp = cp.I_xx * sp.sin(th) ** 2 + cp.I_zz * sp.cos(th) ** 2 + cp.I_xz * sp.sin(2 * th)
    M12 = -(cp.I_xy * sp.sin(th) + cp.I_yz * sp.cos(th))
    MgJ = -cp.m_p * cp.g * (cp.d_x * sp.cos(th) - cp.d_z * sp.sin(th))
    J, k, c, b = ne.J, ne.k, ne.c, ne.b

    # momento delle molle su ciascun GdL (+ sul secondo nodo, - sul primo), come in newton_eulero
    mom = [sp.S.Zero] * ne.n
    off = {"psi": 0, "theta": len(cc.RUOLI_PSI)}
    for asse, archi in (("psi", cc.ARCHI_PSI), ("theta", cc.ARCHI_THETA[variante])):
        for i, (a, bb) in enumerate(archi, 1):
            ia, ib = off[asse] + a, off[asse] + bb
            T = k[f"{asse}_{i}"] * (q[ia] - q[ib]) + c[f"{asse}_{i}"] * (dq[ia] - dq[ib])
            mom[ia] -= T
            mom[ib] += T
    righe = []
    for i, nome in enumerate(ne.nomi):
        if i == iL:
            r = ((J["psi_L"] + cp.J_ax + Jp) * ddq[iL] + M12 * ddq[ip] + sp.diff(Jp, th) * dq[iL] * dq[ip]
                 + sp.diff(M12, th) * dq[ip] ** 2)
        elif i == ip:
            r = (J["theta_p"] + cp.I_yy) * ddq[ip] + M12 * ddq[iL] - sp.Rational(1, 2) * sp.diff(Jp, th) * dq[iL] ** 2 - MgJ
        else:
            r = J[nome] * ddq[i]
        r += -mom[i] + b[nome] * dq[i]
        if nome == "psi_m":
            r -= cp.tau_psi
        if nome == "theta_m":
            r -= cp.tau_theta
        righe.append(r)
    res = sp.Matrix(righe)
    diff = (res - ne.residuo).applyfunc(lambda e: sp.simplify(sp.expand(sp.expand_trig(e))))
    if diff != sp.zeros(ne.n, 1):
        raise AssertionError(f"forma compatta diversa da Newton-Eulero: {diff.T}")
    return FormaCompatta(J_psi_p=Jp, M_12=M12, M_gJ=MgJ, residuo=res)


def _L(e) -> str:
    from .latex_report import L

    return L(e)


def _tabella_parametri(s: cp.Sistema, segnaposto: bool) -> str:
    from .latex_report import num

    nota = rf" \emph{{{NOTA_REALISMO.capitalize()}.}}" if segnaposto else ""
    righe = []
    for mod, simb, beta in ((s.asse_psi, r"\psi", s.beta_psi), (s.asse_theta, r"\theta", s.beta_theta)):
        for e in cc.elementi_equivalenti(mod):
            righe.append(rf"${simb}$ & {e.nome} & {e.nodo_a}--{e.nodo_b} & {', '.join(e.elementi)} & "
                         rf"${num(e.J, 4)}$ & ${num(e.k, 5)}$ & ${num(beta * e.k, 4)}$ \\")
    nodi = []
    for rid, simb, ruoli in ((s.rid_psi, r"\psi", cc.RUOLI_PSI),
                             (s.rid_theta, r"\theta", cc.ruoli_theta(cc.variante_di(s)))):
        for i, r in enumerate(ruoli):
            nodi.append(rf"$J_{{{simb},{r}}}$ & {rid.etichette[i]} & ${num(rid.M_r[i, i], 5)}$ & "
                        rf"${num(rid.completo.J_nodali[rid.completo.masters[i]], 4)}$ & "
                        rf"${num(rid.C_b_r[i, i], 3)}$ \\")
    return "\n".join([
        r"\begin{table}[H]\centering\small",
        r"\caption{Equivalent elements: condensed FEM elements, inertia $J_{eq} = \sum J_i$ (kg\,m$^2$), "
        r"series stiffness $1/k_{eq} = \sum 1/k_i$ (N\,m/rad) and structural damping $c = \beta k$ "
        r"(N\,m\,s/rad)." + nota + r"}\label{tab:lump_elements}",
        r"\begin{tabular}{llllrrr}\toprule",
        r"axis & element & nodes & FEM elements & $J_{eq}$ & $k_{eq}$ & $c$ \\ \midrule",
        *righe,
        r"\bottomrule\end{tabular}\end{table}",
        r"\begin{table}[H]\centering\small",
        r"\caption{Lumped node inertias $J = J^n + \tfrac12\sum_{e\ni\text{node}}J_{eq,e}$ (kg\,m$^2$), with "
        r"$J^n$ the concentrated inertia of the node (fork at azimuth node 9), and bearing friction $b$ "
        r"(N\,m\,s/rad). The telescope is not included: it enters through its own terms." + nota
        + r"}\label{tab:lump_nodes}",
        r"\begin{tabular}{llrrr}\toprule",
        r"symbol & node & $J$ & $J^n$ & $b$ \\ \midrule",
        *nodi,
        r"\bottomrule\end{tabular}\end{table}",
    ])


def _frequenze(s: cp.Sistema, s_mod, segnaposto: bool) -> str:
    """Frequenze degli assi isolati (theta0 = 0): concentrato, FEM completo."""
    from .latex_report import num

    nota = rf" \emph{{{NOTA_REALISMO.capitalize()}.}}" if segnaposto else ""
    righe = []
    for asse in ("azimut", "polare"):
        fc = fem.frequenze_proprie(*s.isolato_ridotto(asse, 0.0)[:2])
        an = s_mod.analisi(asse, 0.0) if s_mod is not None else None
        for i in range(1, len(fc)):
            ff = float(an.f[i]) if an is not None else float("nan")
            righe.append(rf"{in_inglese(asse) if i == 1 else ''} & {i} & ${num(fc[i], 6)}$ & ${num(ff, 6)}$ & "
                         rf"${num((fc[i] / ff - 1) * 100, 3)}\,\%$ \\")
    return "\n".join([
        r"\begin{table}[H]\centering\small",
        r"\caption{Flexible modes of the isolated lumped axes against the full FEM (telescope and fork "
        r"included, $\theta_0 = 0$, $k_g = 0$ for the polar axis)." + nota + r"}\label{tab:lump_freq}",
        r"\begin{tabular}{llrrr}\toprule",
        r"axis & mode & lumped [Hz] & full FEM [Hz] & difference \\ \midrule",
        *righe,
        r"\bottomrule\end{tabular}\end{table}",
    ])


def latex(s: cp.Sistema, segnaposto: bool, *, s_mod=None, extra: Sequence[str] = ()) -> str:
    """Sezione "Lumped-parameter model"; ``extra``: sottosezioni aggiunte in coda (FdT, simulatore)."""
    from .latex_report import num

    variante = cc.variante_di(s)
    mozzo = variante == "mozzo"
    ne = cc.newton_eulero(variante)
    fc = forma_compatta(variante)
    ip = ne.i_theta_p
    th = ne.q[ip]
    sim_q = {ne.q[ip]: sp.Symbol(r"\theta_p", real=True)}
    Jp_lat = _L(fc.J_psi_p.subs(sim_q))
    M12_lat = _L(fc.M_12.subs(sim_q))
    dJp_lat = _L(sp.simplify(sp.diff(fc.J_psi_p, th)).subs(sim_q))
    dM12_lat = _L(sp.diff(fc.M_12, th).subs(sim_q))
    MgJ_lat = _L(fc.M_gJ.subs(sim_q))
    Hz_psi = _L(sp.simplify(ne.H_z.coeff(ne.dq[2])).subs(sim_q))
    Hz_th = _L(sp.simplify(ne.H_z.coeff(ne.dq[ip])).subs(sim_q))
    # forma esplicita: coefficienti trigonometrici scritti per esteso
    Jp, M12 = fc.J_psi_p.subs(sim_q), fc.M_12.subs(sim_q)
    thp = sp.Symbol(r"\theta_p", real=True)
    c_LL = _L(Jp)
    c_Lp = _L(M12)
    c_Lpt = _L(sp.simplify(sp.diff(Jp, thp)))
    c_Ltt = _L(sp.diff(M12, thp))
    c_pLL = _L(sp.expand(-sp.Rational(1, 2) * sp.simplify(sp.diff(Jp, thp))))
    c_g = _L(fc.M_gJ.subs(sim_q))
    testo = "\n".join([
        r"\section{Lumped-parameter model (3 + 4 DOF)}\label{sec:lumped}",
        r"\paragraph{Topology.} A reduced model with the nodes of the original topology only:",
        r"\begin{gather}",
        r"\text{azimuth: } 1 \overset{e_{\psi,1}}{\longleftrightarrow} 5 \overset{e_{\psi,2}}"
        r"{\longleftrightarrow} 9, \qquad \text{polar: } 1 \overset{e_{\theta,1}}{\longleftrightarrow} \mathrm P "
        r"\overset{e_{\theta,2}}{\longleftrightarrow} 7 \overset{e_{\theta,3}}{\longleftrightarrow} 11,",
        r"\end{gather}",
        r"with coordinates $\mathbf q = [\psi_e,\ \psi_m,\ \psi_L,\ \theta_e,\ \theta_p,\ \theta_m,\ "
        r"\theta_r]^{\mathsf T}$ (encoder, motor, load / fork; encoder, telescope, motor, rear end). Each "
        r"equivalent element condenses the FEM elements between two consecutive nodes (module "
        r"\texttt{concentrato.py}):",
        r"\begin{gather}",
        r"J_{eq} = \sum_i J_i, \qquad \frac{1}{k_{eq}} = \sum_i\frac{1}{k_i}\quad\text{(springs in series)},",
        r"\end{gather}",
        r"and any concentrated inertia of an intermediate node is added to $J_{eq}$. The model is a lumped one: "
        r"each node carries its own concentrated inertia plus half of the inertia of the adjacent equivalent "
        r"elements, and the elements are massless springs with the structural damper $c = \beta k$ in "
        r"parallel ($\beta$ tuned on the first flexible mode of each isolated lumped axis, as for the other "
        r"models). The dynamic coupling between the axes is the one of the coupled model "
        r"(Section~\ref{sec:coupled}): the fork turns with $\psi_L$, the telescope is the rigid body attached "
        r"at P, and the polar shafts add the constant $J^{ax}$ about $\hat k$ "
        rf"(here ${num(s.payload.J_ax, 4)}$\,kg\,m$^2$).",
        _tabella_parametri(s, segnaposto),
        _frequenze(s, s_mod, segnaposto),
        r"Halving the inertia of an equivalent element puts part of a heavy disc on the node at the other end "
        r"(polar element 2 on the encoder node, element 9 on the rear end): this, and the missing intermediate "
        r"masses, is why the lumped frequencies differ from the full FEM. The lumped stiffness chain is the "
        r"exact static condensation (the same $\mathbf K_r$ as a static reduction on these nodes).",
        r"\subsection{Newton--Euler equations}\label{subsec:lump_ne}",
        r"Each degree of freedom gets a moment balance on its node: the rate of change of the angular "
        r"momentum equals the external moments about the axis of the node. The torque transmitted by an "
        r"equivalent element, positive on its second node, is",
        r"\begin{gather}",
        r"T_{\psi,1} = k_{\psi,1}(\psi_e - \psi_m) + c_{\psi,1}(\dot\psi_e - \dot\psi_m),\qquad "
        r"T_{\psi,2} = k_{\psi,2}(\psi_m - \psi_L) + c_{\psi,2}(\dot\psi_m - \dot\psi_L),\\",
        r"T_{\theta,1} = k_{\theta,1}(\theta_e - \theta_p) + c_{\theta,1}(\dot\theta_e - \dot\theta_p),\qquad "
        r"T_{\theta,2} = k_{\theta,2}(\theta_p - \theta_m) + c_{\theta,2}(\dot\theta_p - \dot\theta_m),\\ "
        r"T_{\theta,3} = k_{\theta,3}(\theta_m - \theta_r) + c_{\theta,3}(\dot\theta_m - \dot\theta_r),",
        r"\end{gather}",
        r"and $b\,\dot q$ is the viscous friction of a bearing (azimuth: to the ground; polar: to the fork).",
        r"\paragraph{Chain nodes.} For the azimuth nodes the axis $\hat k$ is fixed. The polar nodes are "
        r"axisymmetric about $\hat J$, which turns with the fork: $\dot{\hat J} = \dot\psi_L\hat k\times\hat J"
        r"\perp\hat J$, so the axial component of Euler's equation of each polar node is simply "
        r"$J\ddot\theta$ (its gyroscopic terms are orthogonal to $\hat J$ and are carried by the bearings).",
        r"\paragraph{Telescope.} Both rotation axes pass through $o_2$, a fixed point of the telescope, so its "
        r"Euler equation can be written about $o_2$ with the tensor $\mathbf I_{p,o_2}$: "
        r"$\dot{\mathbf H} = \mathbf M_{o_2}$, $\mathbf H = \mathbf I_{p,o_2}\boldsymbol\omega_2$, with "
        r"$\boldsymbol\omega_2$ of Section~\ref{sec:coupled} and, in body components, "
        r"$\dot{\mathbf H}|_N = \dot{\mathbf H}|_2 + \boldsymbol\omega_2\times\mathbf H$. "
        r"\emph{Azimuth node L} (fork, polar shafts and telescope, rotating about the fixed $\hat z$): only the "
        r"component along $\hat z$ is needed, and the reaction of the polar motor and of the polar bearings is "
        r"orthogonal to $\hat z$, so",
        r"\begin{gather}",
        r"\frac{d}{dt}\Bigl[(J_{\psi,L} + J^{ax})\dot\psi_L + \hat z\cdot\mathbf H\Bigr] = T_{\psi,2} - "
        r"b_{\psi,L}\dot\psi_L,\\",
        r"\hat z\cdot\mathbf H = \bigl[" + Hz_psi + r"\bigr]\,\dot\psi_L + \bigl[" + Hz_th
        + r"\bigr]\,\dot\theta_p.",
        r"\end{gather}",
        r"\emph{Polar node P} (shaft node and telescope): the component of Euler's equation along $\hat J$ "
        r"contains the shaft torques at P and the moment of gravity about $o_2$, "
        r"$\mathbf M_g = \mathbf r_{G} \times (-m_p g\,\hat z)$; the moments orthogonal to $\hat J$ are "
        r"carried by the bearings in the fork.",
        r"\paragraph{Equations of motion, one per degree of freedom.} With",
        r"\begin{gather}",
        r"J_{\psi,p}(\theta_p) = " + Jp_lat + r",\qquad M_{12}(\theta_p) = " + M12_lat + r",\\",
        r"J_{\psi,p}'(\theta_p) = " + dJp_lat + r",\qquad M_{12}'(\theta_p) = " + dM12_lat + r",\\",
        r"M_{g,J}(\theta_p) = " + MgJ_lat + r",",
        r"\end{gather}",
        r"the nonlinear equations are",
        r"\begin{gather}",
        r"(\psi_e)\quad J_{\psi,e}\,\ddot\psi_e = -T_{\psi,1} - b_{\psi,e}\dot\psi_e, \label{eq:ne_psie}\\",
        r"(\psi_m)\quad J_{\psi,m}\,\ddot\psi_m = \tau_\psi + T_{\psi,1} - T_{\psi,2} - b_{\psi,m}\dot\psi_m, "
        r"\label{eq:ne_psim}\\",
        r"(\psi_L)\quad \bigl[J_{\psi,L} + J^{ax} + J_{\psi,p}(\theta_p)\bigr]\ddot\psi_L + M_{12}(\theta_p)\,"
        r"\ddot\theta_p + J_{\psi,p}'(\theta_p)\,\dot\psi_L\dot\theta_p + M_{12}'(\theta_p)\,\dot\theta_p^2 = "
        r"T_{\psi,2} - b_{\psi,L}\dot\psi_L, \label{eq:ne_psiL}\\",
        r"(\theta_e)\quad J_{\theta,e}\,\ddot\theta_e = -T_{\theta,1} - b_{\theta,e}\dot\theta_e, "
        r"\label{eq:ne_the}\\",
        r"(\theta_p)\quad \bigl[J_{\theta,p} + I_{JJ}\bigr]\ddot\theta_p + M_{12}(\theta_p)\,\ddot\psi_L - "
        r"\tfrac12 J_{\psi,p}'(\theta_p)\,\dot\psi_L^2 = T_{\theta,1} - T_{\theta,2} + M_{g,J}(\theta_p) - "
        r"b_{\theta,p}\dot\theta_p, \label{eq:ne_thp}\\",
        r"(\theta_m)\quad J_{\theta,m}\,\ddot\theta_m = \tau_\theta + T_{\theta,2} - T_{\theta,3} - "
        r"b_{\theta,m}\dot\theta_m, \label{eq:ne_thm}\\",
        r"(\theta_r)\quad J_{\theta,r}\,\ddot\theta_r = T_{\theta,3} - b_{\theta,r}\dot\theta_r. "
        r"\label{eq:ne_thr}",
        r"\end{gather}",
        r"The only nonlinearities are the trigonometric functions of $\theta_p$ in $J_{\psi,p}$, $M_{12}$, "
        r"$M_{g,J}$ and the products of velocities: $J_{\psi,p}'\dot\psi_L\dot\theta_p$ (Coriolis, on the "
        r"azimuth), $-\tfrac12 J_{\psi,p}'\dot\psi_L^2$ (centrifugal, on the polar axis) and "
        r"$M_{12}'\dot\theta_p^2$ (from the products of inertia). Written out with the trigonometric "
        r"functions, \eqref{eq:ne_psiL} and \eqref{eq:ne_thp} read",
        r"\begin{gather}",
        r"\begin{split}",
        r"\Bigl[J_{\psi,L} + J^{ax} + " + c_LL + r"\Bigr]\,\ddot\psi_L + \Bigl[" + c_Lp + r"\Bigr]\,"
        r"\ddot\theta_p \\",
        r"{}+ \Bigl[" + c_Lpt + r"\Bigr]\,\dot\psi_L\dot\theta_p + \Bigl[" + c_Ltt + r"\Bigr]\,"
        r"\dot\theta_p^2 = T_{\psi,2} - b_{\psi,L}\dot\psi_L,",
        r"\end{split}\label{eq:ne_psiL_esteso}\\",
        r"\begin{split}",
        r"\bigl[J_{\theta,p} + I_{JJ}\bigr]\,\ddot\theta_p + \Bigl[" + c_Lp + r"\Bigr]\,\ddot\psi_L + "
        r"\Bigl[" + c_pLL + r"\Bigr]\,\dot\psi_L^2 \\",
        r"{}= T_{\theta,1} - T_{\theta,2} " + ("" if c_g.lstrip().startswith("-") else "+ ") + c_g
        + r" - b_{\theta,p}\dot\theta_p.",
        r"\end{split}\label{eq:ne_thp_esteso}",
        r"\end{gather}",
        r"\paragraph{Verification.} The equations above are derived by SymPy from the moment balances "
        r"(function \texttt{newton\_eulero} of \texttt{concentrato.py}), rewritten in the compact form "
        r"\eqref{eq:ne_psie}--\eqref{eq:ne_thr} and checked term by term against both the expanded form and "
        r"Lagrange's equations of the coupled model (Section~\ref{sec:coupled}) with the lumped chain: the "
        r"differences are identically zero. In matrix form they are",
        r"\begin{gather}",
        r"\mathbf M(\theta_p)\ddot{\mathbf q} + \mathbf C(\theta_p, \dot{\mathbf q})\dot{\mathbf q} + "
        r"\mathbf C_d\dot{\mathbf q} + \mathbf K\mathbf q - M_{g,J}(\theta_p)\,\mathbf e_{\theta_p} = "
        r"\mathbf E\boldsymbol\tau,\label{eq:ne_matrix}",
        r"\end{gather}",
        r"where $\mathbf M$ differs from "
        r"$\operatorname{diag}(J_{\psi,e}, \dots, J_{\theta,r})$ only in the entries $(\psi_L,\psi_L)$, "
        r"$(\psi_L,\theta_p)$, $(\theta_p,\theta_p)$.",
        *extra,
    ])
    return _testo_mozzo(testo, s) if mozzo else testo


def _testo_mozzo(testo: str, s: cp.Sistema) -> str:
    """Adatta il testo della sezione al modello ``hub`` (telescopio appeso al centro del mozzo)."""
    from .latex_report import num

    k_b = float(s.rid_theta.k_catena[-1])
    sost = [
        (r"\section{Lumped-parameter model (3 + 4 DOF)}", r"\section{Lumped-parameter model (3 + 5 DOF)}"),
        (r"\text{polar: } 1 \overset{e_{\theta,1}}{\longleftrightarrow} \mathrm P "
         r"\overset{e_{\theta,2}}{\longleftrightarrow} 7 \overset{e_{\theta,3}}{\longleftrightarrow} 11,",
         r"\text{polar: } 1 \overset{e_{\theta,1}}{\longleftrightarrow} \mathrm H "
         r"\overset{e_{\theta,2}}{\longleftrightarrow} 7 \overset{e_{\theta,3}}{\longleftrightarrow} 11,"
         r"\quad \mathrm H \overset{e_{\theta,4}}{\longleftrightarrow} \mathrm P,"),
        (r"with coordinates $\mathbf q = [\psi_e,\ \psi_m,\ \psi_L,\ \theta_e,\ \theta_p,\ \theta_m,\ "
         r"\theta_r]^{\mathsf T}$ (encoder, motor, load / fork; encoder, telescope, motor, rear end).",
         r"with coordinates $\mathbf q = [\psi_e,\ \psi_m,\ \psi_L,\ \theta_e,\ \theta_h,\ \theta_m,\ "
         r"\theta_r,\ \theta_p]^{\mathsf T}$ (encoder, motor, load / fork; encoder, hub centre, motor, rear end, "
         r"telescope). In the hub model the polar chain is a tree: the telescope node P hangs from the hub "
         rf"centre H through the attachment spring $e_{{\theta,4}}$ ($k_{{\theta,4}} = k_b = {num(k_b, 3)}$\,N\,m/rad), "
         r"outside the path encoder--motor."),
        (r"the telescope is the rigid body attached at P,",
         r"the telescope is the rigid body attached at P, which hangs from the hub centre H through $k_b$,"),
        (r"T_{\theta,1} = k_{\theta,1}(\theta_e - \theta_p) + c_{\theta,1}(\dot\theta_e - \dot\theta_p),\qquad "
         r"T_{\theta,2} = k_{\theta,2}(\theta_p - \theta_m) + c_{\theta,2}(\dot\theta_p - \dot\theta_m),\\ "
         r"T_{\theta,3} = k_{\theta,3}(\theta_m - \theta_r) + c_{\theta,3}(\dot\theta_m - \dot\theta_r),",
         r"T_{\theta,1} = k_{\theta,1}(\theta_e - \theta_h) + c_{\theta,1}(\dot\theta_e - \dot\theta_h),\qquad "
         r"T_{\theta,2} = k_{\theta,2}(\theta_h - \theta_m) + c_{\theta,2}(\dot\theta_h - \dot\theta_m),\\ "
         r"T_{\theta,3} = k_{\theta,3}(\theta_m - \theta_r) + c_{\theta,3}(\dot\theta_m - \dot\theta_r),\qquad "
         r"T_{\theta,4} = k_{\theta,4}(\theta_h - \theta_p) + c_{\theta,4}(\dot\theta_h - \dot\theta_p),"),
        (r"\emph{Polar node P} (shaft node and telescope): the component of Euler's equation along $\hat J$ "
         r"contains the shaft torques at P",
         r"\emph{Polar node P} (telescope, off the shaft chain): the component of Euler's equation along "
         r"$\hat J$ contains the torque $T_{\theta,4}$ of the attachment spring"),
        (r"(\theta_e)\quad J_{\theta,e}\,\ddot\theta_e = -T_{\theta,1} - b_{\theta,e}\dot\theta_e, "
         r"\label{eq:ne_the}\\",
         r"(\theta_e)\quad J_{\theta,e}\,\ddot\theta_e = -T_{\theta,1} - b_{\theta,e}\dot\theta_e, "
         r"\label{eq:ne_the}\\" "\n"
         r"(\theta_h)\quad J_{\theta,h}\,\ddot\theta_h = T_{\theta,1} - T_{\theta,2} - T_{\theta,4} - "
         r"b_{\theta,h}\dot\theta_h, \label{eq:ne_thh}\\"),
        (r"\tfrac12 J_{\psi,p}'(\theta_p)\,\dot\psi_L^2 = T_{\theta,1} - T_{\theta,2} + M_{g,J}(\theta_p) - ",
         r"\tfrac12 J_{\psi,p}'(\theta_p)\,\dot\psi_L^2 = T_{\theta,4} + M_{g,J}(\theta_p) - "),
        (r"{}= T_{\theta,1} - T_{\theta,2} ", r"{}= T_{\theta,4} "),
        (r"$\operatorname{diag}(J_{\psi,e}, \dots, J_{\theta,r})$", r"$\operatorname{diag}(J_{\psi,e}, \dots, J_{\theta,p})$"),
    ]
    for vecchio, nuovo in sost:
        if vecchio not in testo:
            raise AssertionError(f"testo da adattare al modello hub non trovato: {vecchio[:60]}")
        testo = testo.replace(vecchio, nuovo)
    return testo


def latex_simulatore(simulazioni: Sequence, figure: dict[str, str], segnaposto: bool,
                     file_esportati: Sequence[str] = ()) -> str:
    """Sottosezione "Nonlinear simulator" della sezione del modello concentrato."""
    from .latex_report import esc, num
    from .simulate import ATOL, RTOL

    nota = rf" \emph{{{NOTA_REALISMO.capitalize()}.}}" if segnaposto else ""
    righe = []
    for r in simulazioni:
        if r.nome == "S3":
            esito = (rf"$|\Delta E|/E = {num(r.derive['energia'], 2)}$, "
                     rf"$|\Delta p_\psi|/p_\psi = {num(r.derive['momento_psi'], 2)}$")
            forz = r"$\mathbf C = \mathbf 0$, $\boldsymbol\tau = \mathbf 0$"
        else:
            e = r.errori
            esito = (rf"$e(\psi_e) = {num(e['psi_e'], 2)}$, $e(\psi_L) = {num(e['psi_L'], 2)}$, "
                     rf"$e(\theta_p) = {num(e['theta_p'], 2)}$")
            forz = (r"doublet on $\tau_\psi$" if r.nome in ("S1", "S2") else
                    rf"polar offset ${np.rad2deg(r.ampiezza):g}^\circ$, $\boldsymbol\tau = \boldsymbol\tau_0$")
        righe.append(rf"{esc(r.nome)} & ${np.rad2deg(r.theta0):g}$ & ${np.rad2deg(r.Omega0):g}$ & {forz} & "
                     rf"{esito} \\")
    fig = []
    for chiave, didascalia in (
        ("S4a", r"S4 at small amplitude: nonlinear (solid) and linearised (dashed) responses almost coincide."),
        ("S4b", r"S4 at larger amplitude: the trigonometric terms detune the polar pendulum and modulate the "
                r"azimuth inertia $J_{\psi,p}(\theta_p)$; the linear model drifts away."),
        ("S4c", r"S4 near the stability limit of the unbalanced telescope: the nonlinear response leaves the "
                r"potential well that the linear model assumes."),
    ):
        if chiave in figure:
            fig += [r"\begin{figure}[H]\centering",
                    rf"\includegraphics[width=0.9\linewidth,height=0.42\textheight,keepaspectratio]{{{figure[chiave]}}}",
                    rf"\caption{{{didascalia}{nota}}}\label{{fig:lump_{chiave}}}\end{{figure}}"]
    cartelle = sorted({f.rsplit("/", 1)[0] for f in file_esportati})
    file_txt = (r" The exported files (" + ", ".join(rf"\texttt{{{esc(f.rsplit('/', 1)[-1])}}}"
                                                     for f in file_esportati)
                + r") are in \texttt{" + esc(", ".join(cartelle)) + "}." if file_esportati else "")
    return "\n".join([
        r"\subsection{Nonlinear simulator}\label{subsec:lump_sim}",
        r"The simulator (module \texttt{simulatore\_concentrato.py}) integrates the Newton--Euler equations "
        r"\eqref{eq:ne_psie}--\eqref{eq:ne_thr} as they are, with all the trigonometric functions of "
        r"$\theta_p$. Writing their residual as $\mathbf r(\mathbf q, \dot{\mathbf q}, \ddot{\mathbf q}, "
        r"\boldsymbol\tau)$,",
        r"\begin{gather}",
        r"\mathbf M(\theta_p) = \frac{\partial\mathbf r}{\partial\ddot{\mathbf q}},\qquad "
        r"\mathbf h(\mathbf q, \dot{\mathbf q}) = \mathbf r\big|_{\ddot{\mathbf q} = \mathbf 0,\ "
        r"\boldsymbol\tau = \mathbf 0},\qquad \dot{\mathbf x} = \begin{bmatrix}\dot{\mathbf q}\\ "
        r"\mathbf M(\theta_p)^{-1}\bigl[\mathbf E\boldsymbol\tau(t) - \mathbf h(\mathbf q, \dot{\mathbf q})\bigr]"
        r"\end{bmatrix},\qquad \mathbf x = [\mathbf q,\ \dot{\mathbf q}],",
        r"\end{gather}",
        r"where $\mathbf h$ collects the gyroscopic and centrifugal terms ($J_{\psi,p}'\dot\psi_L\dot\theta_p$, "
        r"$-\tfrac12 J_{\psi,p}'\dot\psi_L^2$, $M_{12}'\dot\theta_p^2$), the spring and damper torques, the "
        r"bearing friction and the gravity moment $-M_{g,J}(\theta_p)$. Both are compiled from SymPy with "
        r"\texttt{lambdify}; the linear system for $\ddot{\mathbf q}$ is solved at every call. Integration: "
        rf"\texttt{{solve\_ivp}}, DOP853, relative tolerance ${num(RTOL, 1)}$, absolute ${num(ATOL, 1)}$, split at "
        r"the discontinuities of the input. Energy, $p_\psi$ and the linearisation are taken from the "
        r"equivalent Lagrangian model, whose accelerations the simulator reproduces to $10^{-16}$ (test suite).",
        r"\paragraph{Exported code.} \texttt{simulatore\_concentrato.esporta\_codice} writes the same model as "
        r"stand-alone code: a Python module (NumPy only: \texttt{massa(q)}, \texttt{termini(q, dq)}, "
        r"\texttt{rhs(t, x, tau)} for \texttt{solve\_ivp}) and a MATLAB/Octave function for \texttt{ode45}, "
        r"with the parameters of Tables~\ref{tab:lump_elements} and~\ref{tab:lump_nodes} and the telescope in "
        r"the components of frame 2 ($I_{II}, \dots, d_K$). The Python version is checked against the "
        r"simulator to $10^{-16}$; the MATLAB file is generated by the same expressions but has not been run "
        r"(no MATLAB or Octave available here)." + file_txt,
        r"\paragraph{Scenarios.} S1--S3 are those of Section~\ref{sec:workflow}. S4 starts from the steady "
        r"tracking motion at $(\theta_0, \Omega_0)$ with the nominal torques applied and the whole polar axis "
        r"rotated by $\Delta\theta$: the telescope swings as a pendulum about $\theta_0$, with the restoring "
        r"moment $M_{g,J}(\theta_p)$, and the azimuth responds through $J_{\psi,p}(\theta_p)$ and the Coriolis "
        r"term. The same initial condition is integrated with the linearised model; the error "
        r"$e(y) = \max_t|y_{nl} - y_{lin}|/\max_t|y_{lin}|$ grows with $\Delta\theta$, which measures the "
        r"effect of the trigonometric nonlinearities.",
        r"\begin{table}[H]\centering\small",
        r"\caption{Simulations of the lumped model (relative error of the nonlinear against the linearised "
        r"response; conservation drifts for S3)." + nota + r"}\label{tab:lump_sim}",
        r"\begin{tabular}{lrrp{0.27\linewidth}p{0.36\linewidth}}\toprule",
        r"id & $\theta_0$ [deg] & $\Omega_0$ [deg/s] & forcing & result \\ \midrule",
        *righe,
        r"\bottomrule\end{tabular}\end{table}",
        *fig,
    ])


def markdown(s: cp.Sistema, analisi, simulazioni: Sequence, segnaposto: bool,
             file_esportati: Sequence[str] = ()) -> list[str]:
    """Riepilogo del modello a parametri concentrati per ``report.md``."""
    L = ["## Lumped-parameter model (3 + 4 DOF)", ""]
    if segnaposto:
        L += [f"> Every number in this section is a **{NOTA_REALISMO}**.", ""]
    L += ["Equivalent elements: J_eq = sum of the inertias, 1/k_eq = sum of the compliances (springs in series); "
          "masses lumped at the nodes (half of each adjacent element + node inertia). Coupling between the "
          "axes as in the coupled model (telescope at P, fork at L). Newton-Euler equations "
          "verified identical to the Lagrange ones.", "",
          "| axis | element | nodes | FEM elements | J_eq [kg m^2] | k_eq [N m/rad] |", "|---|---|---|---|---|---|"]
    for mod in (s.asse_psi, s.asse_theta):
        for e in cc.elementi_equivalenti(mod):
            L.append(f"| {in_inglese(e.asse)} | {e.nome} | {e.nodo_a}-{e.nodo_b} | {', '.join(e.elementi)} | {e.J:.5g} | {e.k:.6g} |")
    L += ["", "| axis | mode | lumped [Hz] | full FEM [Hz] |", "|---|---|---|---|"]
    for asse in ("azimut", "polare"):
        fc = fem.frequenze_proprie(*s.isolato_ridotto(asse, 0.0)[:2])
        ff = fem.frequenze_proprie(*s.isolato_completo(asse, 0.0)[:2])
        L += [f"| {in_inglese(asse)} | {i} | {fc[i]:.6g} | {ff[i]:.6g} |" for i in range(1, len(fc))]
    L += ["", "Closed-form antiresonances against the zeros of the linearised coupled model (Hz):", ""]
    L += [f"- {nome}: {fc:.6g} (linear {fl:.6g})" for nome, (fc, fl) in analisi.zeri_chiusa.items()]
    e0 = analisi.confronto[0]
    L += ["", f"Against the all-DOF model (T5 criteria, theta0 = {np.rad2deg(e0.theta0):g} deg): direct channels "
          f"{max(e0.err_db):.3g} dB / {max(e0.err_deg):.3g} deg, 2x2 matrix {e0.err_incrociato:.3g}: the lumped "
          "model reproduces the rigid-body motion and the load-spring antiresonance, not the structural modes.", "",
          "Simulations (nonlinear Newton-Euler simulator):", ""]
    for r in simulazioni:
        if r.nome == "S3":
            L.append(f"- S3: drift of E = {r.derive['energia']:.2e}, of p_psi = {r.derive['momento_psi']:.2e}")
        else:
            L.append(f"- {r.nome} (theta0 = {np.rad2deg(r.theta0):g} deg, Omega0 = {np.rad2deg(r.Omega0):g} deg/s"
                     + (f", offset {np.rad2deg(r.ampiezza):g} deg" if r.nome.startswith("S4") else "")
                     + f"): e(psi_e) = {r.errori['psi_e']:.2e}, e(theta_p) = {r.errori['theta_p']:.2e}")
    if file_esportati:
        L += ["", "Simulator code: " + ", ".join(f"`{f}`" for f in file_esportati) + "."]
    return L + [""]
