"""Funzioni di trasferimento del modello a parametri concentrati (3 + 4 GdL).

Stesso schema dell'analisi degli altri modelli:

* FdT **simboliche** delle catene isolate (:math:`\\Omega_0 = 0`, senza smorzamento e senza
  prodotti d'inerzia), :math:`G_{im}(s) = [\\operatorname{adj}\\mathbf P(s)]_{im}/\\det\\mathbf P(s)` con
  :func:`gimbal_mdof.tf.fdt_simboliche_catena`; per il modello concentrato sono esatte;
  zeri e poli in forma chiusa (quadratiche in :math:`s^2`);
* modello **linearizzato** in ogni punto di lavoro (:func:`gimbal_mdof.linearize.linearizza`),
  poli, zeri invarianti (pencil di Rosenbrock), Bode 2x2, sensibilita' a :math:`\\theta_0`;
* **confronto** con il modello a tutti i GdL del FEM con gli stessi criteri del troncamento
  modale (:func:`gimbal_mdof.modal.confronta`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import sympy as sp

from functools import lru_cache

from . import concentrato as cc
from . import coupled as cp
from . import linearize as lz
from . import modal as md
from . import tf


@dataclass(frozen=True)
class FdtSimboliche:
    """FdT simboliche delle catene concentrate isolate e loro forme chiuse."""

    azimut: dict[str, sp.Expr]
    polare: dict[str, sp.Expr]
    zeri: dict[str, sp.Expr]
    """Antirisonanze in forma chiusa, :math:`\\omega^2`."""


@lru_cache(maxsize=2)
def fdt_simboliche(variante: str = "catena") -> FdtSimboliche:
    fa = tf.fdt_simboliche_catena(3, 1, 0, "psi", 2, ("e", "m", "L"))
    rt = cc.ruoli_theta(variante)
    fp = tf.fdt_simboliche_albero(len(rt), cc.ARCHI_THETA[variante], rt.index("m"), rt.index("e"), "theta",
                                  rt.index("p"), rt, con_k_g=True)
    J = {r: sp.Symbol(f"J_{{\\psi,{r}}}", positive=True) for r in ("e", "m", "L")}
    k1, k2 = sp.Symbol("k_{\\psi,1}", positive=True), sp.Symbol("k_{\\psi,2}", positive=True)
    Jt = {r: sp.Symbol(f"J_{{\\theta,{r}}}", positive=True) for r in rt}
    kt3 = sp.Symbol("k_{\\theta,3}", positive=True)
    zeri = {
        "psi_e/tau_psi": k2 / J["L"],        # carico che oscilla su k_psi,2, motore ed encoder fermi
        "psi_m/tau_psi (1)": k1 / J["e"],    # collocata: catena con il nodo motore bloccato,
        "psi_m/tau_psi (2)": k2 / J["L"],    # due rami indipendenti
        "theta_e/tau_theta": kt3 / Jt["r"],  # ramo posteriore oltre il motore
    }
    if variante == "mozzo":
        # telescopio appeso al mozzo: ramo fuori dal percorso encoder-motore, mozzo bloccato
        zeri["theta_e/tau_theta (tel)"] = (sp.Symbol("k_{\\theta,4}", positive=True)
                                           + sp.Symbol("k_g", real=True)) / Jt["p"]
    return FdtSimboliche(azimut=fa, polare=fp, zeri=zeri)


@dataclass(frozen=True)
class AnalisiFdt:
    """Risultati numerici delle FdT del modello concentrato."""

    lineari: dict[tuple[float, float], lz.ModelloLineare]
    confronto: tuple[md.EsitoConfronto, ...]
    """Criteri T5 contro il modello a tutti i GdL (fino a f_band)."""
    zeri_chiusa: dict[str, tuple[float, float]]
    """Per ogni zero in forma chiusa: (valore numerico in Hz, zero invariante del lineare in Hz)."""


def zeri_numerici(s: cp.Sistema, theta0: float) -> dict[str, float]:
    """Valori (Hz) delle antirisonanze in forma chiusa con i parametri concentrati in ``theta0``."""
    rp, rt = s.rid_psi, s.rid_theta
    pl = s.payload
    J_L = rp.M_r[2, 2] + pl.J_ax + float(pl.J_psi(theta0))
    k1, k2 = rp.k_catena
    kt3 = rt.k_catena[2]
    i_r = rt.ruoli.index("r")
    f = lambda w2: float(np.sqrt(w2) / (2 * np.pi))  # noqa: E731
    out = {"psi_e/tau_psi": f(k2 / J_L), "psi_m/tau_psi (1)": f(k1 / rp.M_r[0, 0]),
           "psi_m/tau_psi (2)": f(k2 / J_L), "theta_e/tau_theta": f(kt3 / rt.M_r[i_r, i_r])}
    if rt.completo.ramo:
        i_p = rt.ruoli.index("p")
        out["theta_e/tau_theta (tel)"] = f((rt.k_catena[3] + float(pl.k_g(theta0)))
                                           / (rt.M_r[i_p, i_p] + pl.I_yy))
    return out


@dataclass(frozen=True)
class CatenaPolare:
    """Catena polare concentrata isolata (e, p, m, r) in theta0, con poli e zeri dei canali
    coppia -> encoder e coppia -> motore letti dalla catena (rami bloccati nel nodo d'attacco)."""

    theta0: float
    J: dict[str, float]
    """Inerzie dei nodi e, p, m, r (kg m^2); quella di p comprende il telescopio."""
    k: tuple[float, float, float]
    """k_theta,1 (e-p), k_theta,2 (p-m), k_theta,3 (m-r), N m/rad."""
    k_g: float
    poli: tuple[float, ...]
    """Frequenze proprie della catena libera (Hz), comuni a tutti i canali."""
    zeri_e: tuple[tuple[float, str], ...]
    """Zeri di theta_e/tau_theta: (Hz, ramo che li genera)."""
    zeri_m: tuple[tuple[float, str], ...]
    """Zeri di theta_m/tau_theta: (Hz, ramo che li genera)."""


def catena_polare(s: cp.Sistema, theta0: float) -> CatenaPolare:
    """Poli della catena libera e zeri per blocco dei rami fuori dal percorso motore -> sensore."""
    M, K = (np.asarray(X, dtype=float) for X in s.isolato_ridotto("polare", theta0)[:2])
    ruoli = ("e", "p", "m", "r")
    i = {r: n for n, r in enumerate(ruoli)}
    k1, k2, k3 = -K[i["e"], i["p"]], -K[i["p"], i["m"]], -K[i["m"], i["r"]]
    k_g = float(K[i["p"], i["p"]] - k1 - k2)

    def freq(nodi: Sequence[str]) -> list[float]:
        idx = [i[n] for n in nodi]
        w2 = np.linalg.eigvals(np.linalg.solve(M[np.ix_(idx, idx)], K[np.ix_(idx, idx)])).real
        return sorted(float(np.sqrt(abs(v)) / (2 * np.pi)) for v in w2)

    dietro = tuple((f, "rear end on $k_{\\theta,3}$, motor clamped") for f in freq(["r"]))
    sinistra = freq(["e", "p"])
    nomi_sx = ("telescope on $k_{\\theta,2}$, motor clamped", "encoder on $k_{\\theta,1}$ against the telescope")
    return CatenaPolare(
        theta0=theta0, J={r: float(M[i[r], i[r]]) for r in ruoli}, k=(float(k1), float(k2), float(k3)),
        k_g=k_g, poli=tuple(freq(ruoli)), zeri_e=dietro,
        zeri_m=tuple(sorted(dietro + tuple(zip(sinistra, nomi_sx)))))


@dataclass(frozen=True)
class CatenaRamo:
    """Catena polare con il telescopio appeso al mozzo P tramite una cedevolezza k_b (ramo laterale):
    nodi e, P (mozzo), m, r sul percorso e T (telescopio) fuori dal percorso."""

    theta0: float
    J: dict[str, float]
    """Inerzie dei nodi e, P, m, r, T (kg m^2)."""
    k: tuple[float, float, float, float]
    """k_theta,1 (e-P), k_theta,2 (P-m), k_theta,3 (m-r), k_b (P-T), N m/rad."""
    k_g: float
    poli: tuple[float, ...]
    zeri_e: tuple[tuple[float, str], ...]
    zeri_m: tuple[tuple[float, str], ...]


def _catena_ramo(M: np.ndarray, K: np.ndarray, theta0: float) -> CatenaRamo:
    """Poli e zeri della catena ad albero (e, h, m, r, p): matrici 5 x 5 nell'ordine dei ruoli del
    modello ``hub``; nel disegno il mozzo e' "P" e il telescopio "T"."""
    ruoli = ("e", "h", "m", "r", "p")
    i = {r: n for n, r in enumerate(ruoli)}
    k1, k2, k3, k_b = -K[0, 1], -K[1, 2], -K[2, 3], -K[1, 4]
    k_g = float(K[4, 4] - k_b)

    def freq(nodi) -> list[float]:
        idx = [i[n] for n in nodi]
        w2 = np.linalg.eigvals(np.linalg.solve(M[np.ix_(idx, idx)], K[np.ix_(idx, idx)])).real
        return sorted(float(np.sqrt(abs(v)) / (2 * np.pi)) for v in w2)

    dietro = tuple((f, "rear end on $k_{\\theta,3}$, motor clamped") for f in freq(["r"]))
    tel = tuple((f, "telescope on $k_b$, hub clamped") for f in freq(["p"]))
    sinistra = tuple((f, "encoder-hub-telescope piece, motor clamped") for f in freq(["e", "h", "p"]))
    J = {"e": M[0, 0], "P": M[1, 1], "m": M[2, 2], "r": M[3, 3], "T": M[4, 4]}
    return CatenaRamo(theta0=theta0, J={k_: float(v) for k_, v in J.items()},
                      k=(float(k1), float(k2), float(k3), float(k_b)), k_g=k_g, poli=tuple(freq(ruoli)),
                      zeri_e=tuple(sorted(dietro + tel)), zeri_m=tuple(sorted(dietro + sinistra)))


def catena_ramo_modello(s: cp.Sistema, theta0: float) -> CatenaRamo:
    """Catena concentrata del modello ``hub`` (telescopio appeso al mozzo), isolata in ``theta0``."""
    M, K = (np.asarray(X, dtype=float) for X in s.isolato_ridotto("polare", theta0)[:2])
    return _catena_ramo(M, K, theta0)


def catena_polare_ramo(s: cp.Sistema, theta0: float, J_c: float, k_b: float) -> CatenaRamo:
    """Divide l'inerzia del nodo P della catena concentrata in un mozzo (J_p - J_c) e nel telescopio J_c,
    collegato al mozzo dalla molla k_b; la rigidezza gravitazionale k_g passa al telescopio."""
    c4 = catena_polare(s, theta0)
    k1, k2, k3 = c4.k
    J = [c4.J["e"], c4.J["p"] - J_c, c4.J["m"], c4.J["r"], J_c]
    M = np.diag(J)
    K = np.zeros((5, 5))
    for (u, v), kk in (((0, 1), k1), ((1, 2), k2), ((2, 3), k3), ((1, 4), k_b)):
        K[np.ix_([u, v], [u, v])] += kk * np.array([[1.0, -1.0], [-1.0, 1.0]])
    K[4, 4] += c4.k_g
    return _catena_ramo(M, K, theta0)


def analizza(s: cp.Sistema, s_mod: md.SistemaModale, punti: Sequence[tuple[float, float]]) -> AnalisiFdt:
    lineari = {(th, Om): lz.linearizza(s.accoppiato, th, Om) for th, Om in punti}
    conf = tuple(md.confronta(ml, s_mod.lineare_completo(th, Om), s.par.f_min, s_mod.pm.f_banda)
                 for (th, Om), ml in lineari.items())
    th0 = punti[0][0]
    ml0 = lineari.get((th0, 0.0)) or lz.linearizza(s.accoppiato, th0, 0.0)
    zc = zeri_numerici(s, th0)
    canali = {"psi_e/tau_psi": ("psi_e", 0), "psi_m/tau_psi (1)": ("psi_m", 0), "psi_m/tau_psi (2)": ("psi_m", 0),
              "theta_e/tau_theta": ("theta_e", 1), "theta_e/tau_theta (tel)": ("theta_e", 1)}
    out = {}
    for nome, f in zc.items():
        u, j = canali[nome]
        z = tf.zeri_invarianti(ml0.A, ml0.Bu[:, j], ml0.Cy[ml0.riga(u)])
        fz = np.abs(z[np.imag(z) > 0]) / (2 * np.pi)
        out[nome] = (f, float(fz[np.argmin(np.abs(fz - f))]) if fz.size else float("nan"))
    return AnalisiFdt(lineari=lineari, confronto=conf, zeri_chiusa=out)


def righe_csv(a: AnalisiFdt) -> list[list[object]]:
    righe = []
    for e in a.confronto:
        base = ["confronto_completo", f"{np.rad2deg(e.theta0):g}", f"{np.rad2deg(e.Omega0):g}"]
        righe += [base + ["err_dB_max", f"{max(e.err_db):.4e}"], base + ["err_deg_max", f"{max(e.err_deg):.4e}"],
                  base + ["err_rel_2x2", f"{e.err_incrociato:.4e}"], base + ["passa_T5", "si" if e.passa else "no"]]
        righe += [base + [f"zero {c} {fr:.6g} Hz", f"{fm:.6g} ({err:.2e})"] for c, fr, fm, err in e.zeri]
    for nome, (fc, fl) in a.zeri_chiusa.items():
        righe.append(["zero_forma_chiusa", "", "", nome, f"{fc:.6g} Hz (lineare {fl:.6g} Hz)"])
    return righe


def _lettura_catena_mozzo(s: cp.Sistema, figure: Mapping[str, str], nota: str) -> str:
    """Paragrafo e figura della catena ad albero del modello ``hub`` (telescopio appeso al mozzo)."""
    from .latex_report import num

    c = catena_ramo_modello(s, 0.0)
    hz = lambda v: rf"${num(v, 4)}$\,Hz"  # noqa: E731
    z_T = next(f for f, d in c.zeri_e if d.startswith("telescope"))
    z_r = next(f for f, d in c.zeri_e if d.startswith("rear"))
    z_sx = [f for f, d in c.zeri_m if not d.startswith("rear")]
    return "\n".join([
        r"\paragraph{Reading poles and zeros from the chain.} Figure~\ref{fig:catene_polari} draws the lumped "
        r"polar chain of the hub model at $\theta_0 = 0^\circ$ for the torque$\to$encoder and torque$\to$motor "
        r"channels. The poles belong to the whole free chain and are the same in every channel, "
        + ", ".join(hz(f) for f in c.poli) + r". The zeros are the natural frequencies of the branches "
        r"outside the path from the motor to the sensor, each with its attachment node clamped; the springs "
        r"on the path give only a constant factor, so a chain of $n$ inertias has "
        r"$n - 1 - (\text{springs on the path})$ zeros: two for the encoder, four for the motor. "
        rf"At the encoder the telescope, hanging from the hub through $k_b = {num(c.k[3], 3)}$\,N\,m/rad, gives "
        rf"the zero at {hz(z_T)} with the hub clamped, and the rear end gives the zero at {hz(z_r)} with the "
        r"motor clamped. At the motor, clamping the motor leaves the encoder--hub--telescope piece and the rear "
        r"end, with zeros " + ", ".join(hz(f) for f in sorted(z_sx + [z_r])) + ".",
        r"\begin{figure}[H]\centering",
        rf"\includegraphics[width=\linewidth,height=0.75\textheight,keepaspectratio]{{{figure['catene_polari']}}}",
        r"\caption{Lumped polar chain of the hub model drawn for the two channels: (a) torque$\to$encoder, "
        r"(b) torque$\to$motor. Blue: path from the actuator (torque $\tau_\theta$ on the motor) to the "
        r"sensor, whose springs create no zero. Orange, dashed: branches outside the path, each with the "
        r"zeros given by its natural frequencies with the attachment node clamped (hatched). The poles, "
        r"listed under each panel, are those of the whole free chain." + nota
        + r"}\label{fig:catene_polari}\end{figure}",
    ])


def _lettura_catena(s: cp.Sistema, figure: Mapping[str, str], nota: str) -> str:
    """Paragrafo e figura: come si leggono poli e zeri della catena polare (coppia -> encoder e motore)."""
    from .latex_report import num

    if "catene_polari" not in figure:
        return ""
    if cc.variante_di(s) == "mozzo":
        return _lettura_catena_mozzo(s, figure, nota)
    c = catena_polare(s, 0.0)
    hz = lambda v: rf"${num(v, 4)}$\,Hz"  # noqa: E731
    poli = ", ".join(hz(f) for f in c.poli)
    z_r = c.zeri_e[0][0]
    z_m = [f for f, _ in c.zeri_m]
    f_tel = next(f for f, d in c.zeri_m if d.startswith("telescope"))
    f_enc = next(f for f, d in c.zeri_m if d.startswith("encoder"))
    gravita = (rf" The lowest one, {hz(c.poli[0])}, is the rigid-body mode lifted by the gravity stiffness "
               rf"$k_g = {num(c.k_g, 2)}$\,N\,m/rad of the unbalanced telescope." if c.poli[0] < 1.0 else "")
    sequenza = sorted([(f, "pole") for f in c.poli] + [(f, "zero") for f in z_m])
    # coppie di poli consecutivi senza zeri del canale encoder fra loro
    coppie_e = [(p, q) for p, q in zip(c.poli[:-1], c.poli[1:])
                if not any(p < z < q for z, _ in c.zeri_e)]
    seq_txt = r" $<$ ".join(rf"{n} {num(f, 4)}" for f, n in sequenza)
    return "\n".join([
        r"\paragraph{Reading poles and zeros from the chain.} Figure~\ref{fig:catene_polari} draws the lumped "
        r"polar chain at $\theta_0 = 0^\circ$ twice, once for the torque$\to$encoder channel and once for the "
        r"torque$\to$motor channel; the torque is always applied at the motor. "
        r"\emph{Poles} belong to the whole free chain: they are the roots of $\det(\mathbf K - \omega^2\mathbf M) = 0$, "
        r"involve every inertia and every spring, and are the same in every channel, here " + poli + "." + gravita
        + r" \emph{Zeros} depend on where the torque is applied and where the angle is read. For a chain the "
        r"numerator of $G_{sa}(s) = [\operatorname{adj}\mathbf P(s)]_{sa}/\Delta(s)$ is the product of the springs on "
        r"the path from the actuator $a$ to the sensor $s$ times the determinants of the pieces \emph{outside} the "
        r"path, each with its attachment node clamped:",
        r"\begin{gather}",
        r"\frac{\Theta_e}{T_\theta} = \frac{k_{\theta,1}\,k_{\theta,2}\,\bigl(J_{\theta,r}s^2 + k_{\theta,3}\bigr)}"
        r"{\Delta_\theta(s)},\qquad"
        r"\frac{\Theta_m}{T_\theta} = \frac{\bigl(J_{\theta,r}s^2 + k_{\theta,3}\bigr)\,\Delta_{ep|m}(s)}"
        r"{\Delta_\theta(s)},",
        r"\end{gather}",
        r"where $\Delta_{ep|m}(s)$ is the determinant of the encoder--telescope piece with the motor clamped. "
        r"The springs on the path ($k_{\theta,1}k_{\theta,2}$ for the encoder, none for the collocated motor) "
        r"enter only as a constant gain and create no zero; each of them removes one zero, so a chain of $n$ "
        r"inertias has $n - 1 - (\text{springs on the path})$ zeros: one for the encoder, three for the motor. "
        r"Each off-path branch gives the zeros of its own natural frequencies with the attachment node "
        r"clamped, $\omega_z \approx \sqrt{k_{branch}/J_{branch}}$: at that frequency the branch oscillates in "
        r"antiphase and its spring returns to the attachment node a torque equal and opposite to the motor "
        r"torque (a dynamic absorber), so the attachment node and everything beyond it, sensor included, "
        r"stand still.",
        r"\begin{itemize}",
        r"\item \textbf{Torque$\to$encoder (non-collocated).} The telescope is \emph{on} the path between the "
        r"encoder and the motor: its inertia appears only in $\Delta_\theta$, i.e.\ in the poles. The only "
        rf"off-path branch is the rear end on $k_{{\theta,3}}$, which gives the zero at {hz(z_r)}: there, "
        r"in steady state, the rear end oscillates with $\theta_r = -T_\theta/k_{\theta,3}$, the torque of its "
        r"spring cancels the motor torque, and motor, telescope and encoder stand still. The channel has "
        rf"{len(c.poli)} poles but only {len(c.zeri_e)} zero" + ("s" if len(c.zeri_e) != 1 else "")
        + r": the pole pairs " + " and ".join(rf"({num(p, 4)}, {num(q, 4)})\,Hz" for p, q in coppie_e)
        + r" have no zero between them, so after each of them the phase falls by a further $180^\circ$ "
        r"instead of recovering: this is the signature of a non-collocated sensor.",
        r"\item \textbf{Torque$\to$motor (collocated).} Clamping the motor splits the chain into the "
        r"encoder--telescope piece and the rear end, and all their frequencies are zeros: "
        rf"{hz(f_tel)} (telescope on $k_{{\theta,2}}$), {hz(z_r)} (rear end on $k_{{\theta,3}}$) and {hz(f_enc)} "
        r"(the light encoder on $k_{\theta,1}$ against the heavy telescope, which barely moves the motor: this "
        r"zero almost cancels the pole at the same frequency). Poles and zeros interlace, " + seq_txt
        + r"\,Hz, and the phase of the collocated channel stays bounded.",
        r"\end{itemize}",
        r"Heavy branches on soft springs give zeros in the control band; light branches (the encoder side, "
        r"the rear end) give zeros far above it. The antiresonance measured at the encoder of the elevation "
        r"(Section~\ref{subsec:exp_discussion}) therefore requires a large inertia \emph{outside} the "
        r"encoder--motor path, which this chain does not contain: here the telescope can produce a zero only "
        rf"in the motor channel ({hz(f_tel)}).",
        r"\begin{figure}[H]\centering",
        rf"\includegraphics[width=\linewidth,height=0.75\textheight,keepaspectratio]{{{figure['catene_polari']}}}",
        r"\caption{Lumped polar chain drawn for the two channels: (a) torque$\to$encoder, (b) torque$\to$motor. "
        r"Blue: path from the actuator (torque $\tau_\theta$ on the motor) to the sensor, whose springs create no "
        r"zero. Orange, dashed: branches outside the path, each with the zeros given by its natural frequencies "
        r"with the motor clamped (hatched). The poles, listed under each panel, are those of the whole free "
        r"chain." + nota + r"}\label{fig:catene_polari}\end{figure}",
    ])


def latex(s: cp.Sistema, a: AnalisiFdt, figure: Mapping[str, str], segnaposto: bool) -> str:
    """Sottosezione "Transfer functions" della sezione del modello concentrato."""
    from .latex_report import L, num, ridimensiona, tabella_poli_zeri
    from .params import NOTA_REALISMO

    mozzo = cc.variante_di(s) == "mozzo"
    fs = fdt_simboliche(cc.variante_di(s))
    fa, fp = fs.azimut, fs.polare
    s_ = sp.Symbol("s")
    poli = sp.Poly(sp.expand(sp.cancel(fa["det"] / s_**2)), s_)
    coeff = {k: sp.factor(poli.coeff_monomial(s_**k)) for k in (0, 2, 4)}
    nota = rf" \emph{{{NOTA_REALISMO.capitalize()}.}}" if segnaposto else ""
    th0 = next(iter(a.lineari))[0]
    righe_conf = []
    for e in a.confronto:
        z = max((err for *_, err in e.zeri), default=0.0)
        righe_conf.append(rf"${np.rad2deg(e.theta0):g}$ & ${np.rad2deg(e.Omega0):g}$ & ${num(max(e.err_db), 3)}$ & "
                          rf"${num(max(e.err_deg), 3)}$ & ${num(z, 2)}$ & ${num(e.err_incrociato, 2)}$ \\")
    nomi_z = {"psi_e/tau_psi": r"$\Psi_e/T_\psi$: $\sqrt{k_{\psi,2}/J_{\psi,L}^*}$",
              "psi_m/tau_psi (1)": r"$\Psi_m/T_\psi$: $\sqrt{k_{\psi,1}/J_{\psi,e}}$",
              "psi_m/tau_psi (2)": r"$\Psi_m/T_\psi$: $\sqrt{k_{\psi,2}/J_{\psi,L}^*}$",
              "theta_e/tau_theta": r"$\Theta_e/T_\theta$: $\sqrt{k_{\theta,3}/J_{\theta,r}}$",
              "theta_e/tau_theta (tel)": r"$\Theta_e/T_\theta$: $\sqrt{(k_{\theta,4} + k_g)/J_{\theta,p}}$"}
    righe_z = [rf"{nomi_z[c]} & ${num(fc, 6)}$ & ${num(fl, 6)}$ \\" for c, (fc, fl) in a.zeri_chiusa.items()]
    fig = []
    for chiave, didascalia in (
        ("bode_concentrato", r"Torque$\to$encoder $2\times2$ matrix: lumped model (with its poles, dotted, and "
                             r"invariant zeros, circles), full-DOF reference and modal model."),
        ("sweep_concentrato", r"Natural frequencies of the coupled lumped model as a function of $\theta_0$."),
    ):
        if chiave in figure:
            fig += [r"\begin{figure}[H]\centering",
                    rf"\includegraphics[width=\linewidth,height=0.8\textheight,keepaspectratio]{{{figure[chiave]}}}",
                    rf"\caption{{{didascalia}{nota}}}\label{{fig:{chiave}}}\end{{figure}}"]
    ml0 = a.lineari.get((th0, 0.0)) or next(iter(a.lineari.values()))
    tab_pz = tabella_poli_zeri(ml0).replace(r"\label{tab:polezero}", r"\label{tab:polezero_lump}")
    return "\n".join([
        r"\subsection{Transfer functions}\label{subsec:lump_tf}",
        r"\paragraph{Symbolic transfer functions of the isolated chains} ($\Omega_0 = 0$, no damping, zero "
        r"products of inertia). For the lumped model the chain is exactly the inertia--spring--inertia "
        r"system, so $\mathbf P(s) = \operatorname{diag}(J_i)s^2 + \mathbf K$ and $G_{im}(s) = "
        r"[\operatorname{adj}\mathbf P(s)]_{im}/\det\mathbf P(s)$ hold without approximation. Azimuth (torque "
        r"on $\psi_m$; $J_{\psi,L}^* = J_{\psi,L} + J^{ax} + J_{\psi,p}(\theta_0)$ written $J_{\psi,L}$):",
        r"\begin{gather}",
        r"\mathbf P_\psi(s) = " + L(fa["P"]) + r",\qquad \Delta_\psi(s) = \det\mathbf P_\psi(s), \\",
        r"\frac{\Psi_e}{T_\psi} = \frac{" + L(fa["N_encoder"]) + r"}{\Delta_\psi},\qquad "
        r"\frac{\Psi_m}{T_\psi} = \frac{" + L(fa["N_motore"]) + r"}{\Delta_\psi},\qquad "
        r"\frac{\Psi_L}{T_\psi} = \frac{" + L(fa["N_carico"]) + r"}{\Delta_\psi},",
        r"\end{gather}",
        r"\begin{gather}",
        r"\Delta_\psi(s) = s^2\bigl(a_4 s^4 + a_2 s^2 + a_0\bigr),\qquad a_4 = " + L(coeff[4]) + r", \\",
        r"a_2 = " + L(coeff[2]) + r",\qquad a_0 = " + L(coeff[0]) + r",",
        r"\end{gather}",
        r"so the azimuth has the rigid double pole and two flexible poles, $\omega^2$ being the roots of "
        r"the quadratic $\Delta_\psi/s^2$ in $-s^2$. Polar (torque on $\theta_m$; $J_{\theta,p}$ includes "
        r"$I_{JJ}$; $k_g$ acts on the telescope node):",
        *([ridimensiona(r"\mathbf P_\theta(s) = " + L(fp["P"]) + ",")] if mozzo else []),
        r"\begin{gather}",
        *([] if mozzo else [r"\mathbf P_\theta(s) = " + L(fp["P"]) + r", \\"]),
        r"\frac{\Theta_e}{T_\theta} = \frac{" + L(fp["N_encoder"]) + r"}{\Delta_\theta},\qquad "
        r"\frac{\Theta_p}{T_\theta} = \frac{" + L(fp["N_carico"]) + r"}{\Delta_\theta},",
        r"\end{gather}",
        *([r"\begin{gather}",
            r"\frac{\Theta_m}{T_\theta} = \frac{\bigl(J_{\theta,r}s^2 + k_{\theta,3}\bigr)\,\Delta_{ehp|m}(s)}"
            r"{\Delta_\theta}, \\ \Delta_{ehp|m}(s) = \det" + L(fp["P"].extract([0, 1, 4], [0, 1, 4])) + ",",
            r"\end{gather}",
            r"where $\Delta_{ehp|m}$ is the determinant of the encoder--hub--telescope piece with the motor "
            r"clamped."] if mozzo else
          [ridimensiona(r"\frac{\Theta_m}{T_\theta} = \frac{" + L(fp["N_motore"]) + r"}{\Delta_\theta}.")]),
        r"The zeros follow by inspection (Figure~\ref{fig:catene_polari}). In the azimuth the load is beyond "
        r"the motor with respect to the encoder, so $\Psi_e/T_\psi$ vanishes when the load resonates on "
        r"$k_{\psi,2}$ with the motor clamped. "
        + (r"In the hub model the telescope hangs from the hub centre through $k_{\theta,4} = k_b$, outside "
           r"the path encoder--motor, so $\Theta_e/T_\theta$ vanishes both when the rear end resonates with "
           r"the motor clamped and when the telescope resonates on $k_b$ with the hub clamped. " if mozzo else
           r"In the polar axis the telescope lies between the encoder and "
           r"the motor, so $\Theta_e/T_\theta$ keeps only the zero of the rear end beyond the motor. ")
        + 
        r"Table~\ref{tab:lump_zeros} compares these closed forms, evaluated with the lumped "
        rf"parameters at $\theta_0 = {np.rad2deg(th0):g}^\circ$, with the invariant zeros of the coupled "
        r"linearised model.",
        r"\begin{table}[H]\centering\small",
        r"\caption{Antiresonances of the lumped model: closed form against the invariant zeros of the "
        r"coupled linearised model (Hz)." + nota + r"}\label{tab:lump_zeros}",
        r"\begin{tabular}{lrr}\toprule",
        r"zero & closed form & coupled model \\ \midrule",
        *righe_z,
        r"\bottomrule\end{tabular}\end{table}",
        _lettura_catena(s, figure, nota),
        r"\paragraph{Linearised coupled model.} The lumped model is linearised at every operating point exactly "
        r"as the other models (Section~\ref{sec:linear}): $\mathbf M_0$, $\mathbf C_0$, $\mathbf K_0$ from the "
        r"Jacobians of the residual, FRF from $\mathbf P(j\omega)^{-1}$, poles from $\operatorname{eig}(\mathbf A)$, "
        r"invariant zeros from the Rosenbrock pencil. Table~\ref{tab:lump_vs_full} applies the criteria of "
        r"Section~\ref{sec:modal} against the full-DOF model: the lumped model reproduces the rigid and "
        r"low-frequency behaviour, including the load antiresonance, but not the structural modes, whose "
        r"frequencies it shifts (Table~\ref{tab:lump_freq}).",
        r"\begin{table}[H]\centering\small",
        r"\caption{Lumped model against the full-DOF model, torque$\to$encoder up to the band of interest: "
        r"worst direct-channel error outside the zero bands, worst relative error of the zeros, worst relative "
        r"2-norm error of the $2\times2$ matrix." + nota + r"}\label{tab:lump_vs_full}",
        r"\begin{tabular}{rrrrrr}\toprule",
        r"$\theta_0$ & $\Omega_0$ & dB & deg & zeros & $2\times2$ \\ \midrule",
        *righe_conf,
        r"\bottomrule\end{tabular}\end{table}",
        *fig,
        tab_pz,
    ])
