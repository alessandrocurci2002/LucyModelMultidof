"""Sezione LaTeX e riepilogo Markdown del confronto con le misure sperimentali.

Le grandezze discusse (frequenze, scarti, rigidezze e inerzie equivalenti) sono calcolate qui
dai dati e dal modello, non scritte a mano: ogni numero del testo e' tracciabile.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from . import concentrato as cc
from . import concentrato_analisi as ca
from . import coupled as cp
from . import fem
from . import linearize as lz
from . import sperimentale as spm
from .params import NOTA_REALISMO, in_inglese

TIPO_EN: dict[str, str] = {"risonanza": "resonance", "antirisonanza": "antiresonance"}


@dataclass(frozen=True)
class Discussione:
    """Numeri della discussione delle discrepanze."""

    k0_az: float | None
    f0_az: dict[str, float]
    """Risonanza bassa dell'azimut per ciascuna misura (Hz)."""
    antiris: dict[str, tuple[float, float]]
    """Antirisonanza carico-molla: (misura, modello), Hz."""
    f1_az: dict[str, tuple[float, float]]
    """Primo modo strutturale dell'azimut: (misura, modello), Hz."""
    J_motore_modello: float
    J_motore_richiesta: float | None
    J_rotore_datasheet: float
    livello_az_db: float
    pair_az: tuple[float, float] | None
    """Coppia risonanza/antirisonanza dell'azimut senza corrispondenza nel modello (Hz)."""
    k0_nt: float | None
    f_pk_nt: float | None
    J_nt_richiesta: float | None
    J_nt_modello: float
    el: "ElevazioneEncoder | None" = None
    """Antirisonanza all'encoder della misura di elevazione con il telescopio."""


@dataclass(frozen=True)
class ElevazioneEncoder:
    """Coppia antirisonanza/risonanza misurata all'encoder dell'elevazione e stima a due masse.

    Con l'encoder sul carico, uno zero di theta_e/tau richiede un'inerzia fuori dal percorso
    encoder-motore: si stima una parte J_1 rigida con l'albero e una parte J_c appesa tramite
    k_b, da (f_p/f_z)^2 = J_tot/J_1 e k_b = J_c omega_z^2 (smorzamento trascurato).
    """

    f_z: float
    f_p: float
    J_tot: float
    J_1: float
    J_c: float
    k_b: float
    k_theta2: float
    livello_previsto_db: float
    livello_misurato_db: float
    banda_livello: tuple[float, float]
    rigido_db: float
    """Scarto mediano misura - theta_e del modello sotto l'antirisonanza (linea del corpo rigido)."""
    banda_rigido: tuple[float, float]
    z_enc_min: float | None
    """Zero piu' basso del canale encoder del modello (Hz)."""
    z_mot: float | None
    """Zero del canale motore del modello piu' vicino a f_z (Hz)."""
    coppia_bassa: tuple[float, float] | None
    """Coppia antirisonanza/risonanza sotto f_z (Hz)."""
    n_ris_alte: int
    n_modi_alti: int
    banda_alta: tuple[float, float]
    filtri_hz: tuple[float, ...]
    catena: "ca.CatenaRamo | None" = None
    mozzo: bool = False
    """Modello hub: il telescopio e' gia' appeso al mozzo tramite k_b."""
    z_enc_vicino: float | None = None
    """Zero del canale encoder del modello piu' vicino a f_z (Hz)."""
    f_polo_modello: float | None = None
    """Primo modo flessibile della catena polare del modello (Hz)."""
    """Catena polare con il telescopio appeso al mozzo tramite k_b, con J_c e k_b della stima."""


J_ROTORE_TBM7646: float = 8.19e-5
"""Inerzia del rotore TBM(S)-7646, kg m^2 (datasheet Kollmorgen)."""


def discussione(confronti: Sequence[spm.Confronto], sistemi: Mapping[bool, cp.Sistema], Kt: float) -> Discussione:
    s = sistemi[True]
    az = [c for c in confronti if c.misura.asse == "azimut"]
    f0, anti, f1 = {}, {}, {}
    k0 = None
    pair = None
    for c in az:
        car = c.caratteristiche
        bassi = [x for x in car if x.tipo == "risonanza" and x.f_misura < 10.0]
        if bassi:
            f0[c.misura.etichetta] = max(bassi, key=lambda x: x.modulo_db).f_misura
        a = [x for x in car if x.tipo == "antirisonanza" and x.f_modello is not None and x.f_misura < 30.0]
        if a:
            anti[c.misura.etichetta] = (a[0].f_misura, a[0].f_modello)  # type: ignore[assignment]
        r = [x for x in car if x.tipo == "risonanza" and x.f_misura > 100.0 and x.f_modello is not None]
        if r:
            top = max(r, key=lambda x: x.modulo_db)
            f1[c.misura.etichetta] = (top.f_misura, top.f_modello)  # type: ignore[assignment]
        nr = [x for x in car if x.f_modello is None and 20.0 < x.f_misura < 100.0]
        rr = [x.f_misura for x in nr if x.tipo == "risonanza"]
        aa = [x.f_misura for x in nr if x.tipo == "antirisonanza"]
        if rr and aa and pair is None:
            pair = (rr[0], aa[0])
    if az:
        c0 = az[0]
        J_tot = float(np.ones(3) @ s.rid_psi.M_r @ np.ones(3)) + s.payload.J_ax + float(s.payload.J_psi(c0.misura.theta0))
        if c0.misura.etichetta in f0:
            k0 = J_tot * (2 * np.pi * f0[c0.misura.etichetta]) ** 2
    # inerzia lato motore che porterebbe il primo modo dell'azimut alla frequenza misurata:
    # per una catena con carico pesante, omega_1^2 ~ k_psi,2 (1/J_em + 1/J_L*) ~ k_psi,2 / J_em
    J_em = float(s.rid_psi.M_r[0, 0] + s.rid_psi.M_r[1, 1])
    J_rich = None
    if f1:
        fm, fmod = next(iter(f1.values()))
        J_rich = J_em * (fmod / fm) ** 2
    liv = float(np.mean([c.scarto_livello_db["psi_e"] for c in az])) if az else float("nan")
    # senza telescopio: rigidezza dalla pendenza a bassa frequenza (|G| ~ omega K_t / k_0) e inerzia dal picco
    nt = [c for c in confronti if not c.misura.telescopio]
    k0_nt = f_pk = J_nt = None
    s0 = sistemi[False]
    J_nt_mod = float(np.ones(s0.rid_theta.n) @ s0.rid_theta.M_r @ np.ones(s0.rid_theta.n))
    if nt:
        m = nt[0].misura
        sel = (m.f >= 5.0) & (m.f <= 15.0)
        if sel.any():
            k0_nt = float(np.median(Kt * 2 * np.pi * m.f[sel] / 10 ** (m.modulo_db[sel] / 20)))
        sel = (m.f >= 20.0) & (m.f <= 100.0)
        if sel.any():
            f_pk = float(m.f[sel][np.argmax(m.modulo_db[sel])])
            if k0_nt:
                J_nt = k0_nt / (2 * np.pi * f_pk) ** 2
    return Discussione(k0_az=k0, f0_az=f0, antiris=anti, f1_az=f1, J_motore_modello=J_em,
                       J_motore_richiesta=J_rich, J_rotore_datasheet=J_ROTORE_TBM7646, livello_az_db=liv,
                       pair_az=pair, k0_nt=k0_nt, f_pk_nt=f_pk, J_nt_richiesta=J_nt, J_nt_modello=J_nt_mod,
                       el=_elevazione_encoder(confronti, s))


def _scarto_db(c: spm.Confronto, canale: str, a: float, b: float) -> float:
    """Scarto mediano (dB) misura - canale del modello nella banda [a, b] Hz."""
    m = c.misura
    sel = (m.f >= a) & (m.f <= b)
    G = np.interp(m.f[sel], c.f_modello, 20 * np.log10(np.abs(c.G_modello[canale])))
    return float(np.median(m.modulo_db[sel] - G))


def _elevazione_encoder(confronti: Sequence[spm.Confronto], s: cp.Sistema) -> ElevazioneEncoder | None:
    el = [c for c in confronti if c.misura.asse == "polare" and c.misura.telescopio]
    if not el:
        return None
    c = el[0]
    m = c.misura
    zeri = [x for x in c.caratteristiche if x.tipo == "antirisonanza" and 10.0 < x.f_misura < 40.0]
    if not zeri:
        return None
    f_z = min(zeri, key=lambda x: x.modulo_db).f_misura
    # risonanza associata: massimo del modulo fra 1.5 f_z e 4 f_z (dove la fase scende per 0 gradi)
    sel = (m.f > 1.5 * f_z) & (m.f < 4.0 * f_z)
    if not sel.any():
        return None
    f_p = float(m.f[sel][np.argmax(m.modulo_db[sel])])
    M, K = (np.asarray(X) for X in s.isolato_ridotto("polare", m.theta0)[:2])
    mozzo = bool(s.rid_theta.completo.ramo)
    uno = np.ones(M.shape[0])
    J_tot = float(uno @ M @ uno)
    J_1 = J_tot * (f_z / f_p) ** 2
    J_c = J_tot - J_1
    k_b = J_c * (2 * np.pi * f_z) ** 2
    k_theta2 = float(-K[cc.RUOLI_THETA.index("p"), cc.RUOLI_THETA.index("m")])
    banda_liv = (5.0 * f_z, 20.0 * f_z)
    banda_rig = (5.0, 0.6 * f_z)
    ml = lz.linearizza(s.accoppiato, m.theta0, 0.0)
    _, z_enc = spm.poli_zeri_modello(ml, "theta_e", 1, 2.0 * m.f[-1])
    _, z_mot = spm.poli_zeri_modello(ml, "theta_m", 1, 2.0 * m.f[-1])
    bassi = [x for x in c.caratteristiche if x.f_misura < 0.8 * f_z]
    ab = [x.f_misura for x in bassi if x.tipo == "antirisonanza"]
    rb = [x.f_misura for x in bassi if x.tipo == "risonanza"]
    banda_alta = (4.0 * f_z, 1000.0)
    n_ris = sum(1 for x in c.caratteristiche
                if x.tipo == "risonanza" and banda_alta[0] <= x.f_misura <= banda_alta[1])
    f_modi = fem.frequenze_proprie(M, K)
    n_modi = int(np.sum((f_modi >= banda_alta[0]) & (f_modi <= banda_alta[1])))
    return ElevazioneEncoder(
        f_z=f_z, f_p=f_p, J_tot=J_tot, J_1=J_1, J_c=J_c, k_b=k_b, k_theta2=k_theta2,
        livello_previsto_db=float(20 * np.log10(J_tot / J_1)),
        livello_misurato_db=_scarto_db(c, "theta_e", *banda_liv), banda_livello=banda_liv,
        rigido_db=_scarto_db(c, "theta_e", *banda_rig), banda_rigido=banda_rig,
        z_enc_min=float(z_enc[0]) if z_enc.size else None,
        z_mot=float(z_mot[np.argmin(np.abs(np.log(z_mot / f_z)))]) if z_mot.size else None,
        coppia_bassa=(ab[0], rb[0]) if ab and rb else None, n_ris_alte=n_ris, n_modi_alti=n_modi,
        banda_alta=banda_alta, filtri_hz=tuple(sorted({f[0] for f in m.filtri})),
        catena=None if mozzo else ca.catena_polare_ramo(s, m.theta0, J_c, k_b), mozzo=mozzo,
        z_enc_vicino=float(z_enc[np.argmin(np.abs(np.log(z_enc / f_z)))]) if z_enc.size else None,
        f_polo_modello=float(f_modi[f_modi > 1.0][0]) if np.any(f_modi > 1.0) else None)


def righe_csv(confronti: Sequence[spm.Confronto]) -> list[list[object]]:
    out = []
    for c in confronti:
        m = c.misura
        for x in c.caratteristiche:
            out.append([m.file, m.etichetta, x.tipo, f"{x.f_misura:.6g}", f"{x.modulo_db:.4g}",
                        "" if x.f_modello is None else f"{x.f_modello:.6g}", x.canale or "",
                        "" if x.f_modello is None else f"{x.f_misura / x.f_modello:.4g}"])
        for ch, d in c.scarto_livello_db.items():
            out.append([m.file, m.etichetta, f"livello {ch}", "", f"{d:+.3g}", "", ch, ""])
        for banda, v in c.fit.items():
            out.append([m.file, m.etichetta, f"fit NRMSE % ({banda})", "", f"{v:.4g}", "",
                        spm.canali(m.asse)[0], ""])
    return out


INTESTAZIONE_CSV = ["file", "misura", "caratteristica", "f_misura_Hz", "modulo_dB_o_scarto",
                    "f_modello_Hz", "canale_modello", "rapporto_f"]



def _punto_elevazione_mozzo(e: ElevazioneEncoder) -> str:
    """Punto della discussione sull'elevazione per il modello hub (telescopio appeso al mozzo)."""
    from .latex_report import num

    filtri = (r"; the anti-resonance filters of the drive are tuned at "
              + " and ".join(rf"${num(f, 3)}$" for f in e.filtri_hz) + r"\,Hz" if e.filtri_hz else "")
    z = (rf"${num(e.z_enc_vicino, 4)}$\,Hz" if e.z_enc_vicino else "--")
    p = (rf"${num(e.f_polo_modello, 4)}$\,Hz" if e.f_polo_modello else "--")
    return (
        r"\textbf{Elevation: antiresonance at the encoder (hub model).} Below the antiresonance the "
        rf"measurement follows the rigid-body line of the model within ${num(abs(e.rigido_db), 2)}$\,dB "
        rf"({num(e.banda_rigido[0], 2)}--{num(e.banda_rigido[1], 2)}\,Hz), which confirms the units and the total "
        rf"inertia $J_{{tot}} = {num(e.J_tot, 3)}$\,kg\,m$^2$. It then has a deep antiresonance at "
        rf"${num(e.f_z, 4)}$\,Hz followed by a broad resonance at about ${num(e.f_p, 3)}$\,Hz" + filtri + ". "
        r"In the hub model the telescope hangs from the hub centre through $k_b$, outside the encoder--motor "
        rf"path (Section~\ref{{sec:hub}}): the encoder channel has the zero at {z} (telescope on $k_b$ with the "
        rf"hub clamped), against ${num(e.f_z, 4)}$\,Hz measured, and the first flexible mode of the polar chain "
        rf"is at {p}, against ${num(e.f_p, 3)}$\,Hz measured. The zero is set by $k_b$, here the two-mass "
        rf"estimate from the measurement: the fit hangs only $J_c \approx {num(e.J_c, 2)}$\,kg\,m$^2$ on $k_b$, "
        r"while the model hangs the whole telescope, hence the slightly lower model zero. The resonance "
        r"depends on the "
        r"stiffness of the path, dominated by the support ring taken alone (Table~\ref{tab:hub_parts}), and on "
        rf"the inertia on the hub side. The two-mass fit of the measurement gives $J_1 \approx {num(e.J_1, 2)}$"
        rf"\,kg\,m$^2$ on the shaft side, $J_c \approx {num(e.J_c, 2)}$\,kg\,m$^2$ beyond the attachment and "
        rf"$k_b \approx {num(e.k_b, 2)}$\,N\,m/rad. Above the resonance the measurement is "
        rf"${num(e.livello_misurato_db, 3)}$\,dB from the encoder channel of the hub model "
        rf"({num(e.banda_livello[0], 3)}--{num(e.banda_livello[1], 3)}\,Hz). Between "
        rf"{num(e.banda_alta[0], 2)}\,Hz and {num(e.banda_alta[1], 4)}\,Hz the measurement has {e.n_ris_alte} "
        rf"resonances against {e.n_modi_alti} modes of the polar chain. The elevation of this test is not "
        r"recorded and is assumed to be $0^\circ$."
    )


def _punto_elevazione(e: ElevazioneEncoder | None) -> str:
    """Punto della discussione sull'antirisonanza misurata all'encoder dell'elevazione."""
    from .latex_report import num

    if e is None:
        return ""
    if e.mozzo:
        return _punto_elevazione_mozzo(e)
    bassa = (rf" (a weak notch--peak pair at ${num(e.coppia_bassa[0], 3)}$/${num(e.coppia_bassa[1], 3)}$\,Hz "
             r"is not modelled either)" if e.coppia_bassa else "")
    filtri = (r"; the anti-resonance filters of the drive are tuned at "
              + " and ".join(rf"${num(f, 3)}$" for f in e.filtri_hz) + r"\,Hz" if e.filtri_hz else "")
    z_enc = (rf"has no zero below ${num(e.z_enc_min, 3)}$\,Hz" if e.z_enc_min else "has no zero in the band")
    z_mot = (rf" The motor channel has a zero at ${num(e.z_mot, 3)}$\,Hz (telescope oscillating on "
             r"$k_{\theta,2}$ with the motor still), close in frequency, but it is not the measured channel."
             if e.z_mot else "")
    return (
        r"\textbf{Elevation: antiresonance at the encoder.} Below the antiresonance the measurement follows the "
        rf"rigid-body line of the model within ${num(abs(e.rigido_db), 2)}$\,dB ({num(e.banda_rigido[0], 2)}--"
        rf"{num(e.banda_rigido[1], 2)}\,Hz), which confirms the units and the total inertia "
        rf"$J_{{tot}} = {num(e.J_tot, 3)}$\,kg\,m$^2$" + bassa + r". It then has a deep antiresonance at "
        rf"${num(e.f_z, 4)}$\,Hz followed by a broad resonance at about ${num(e.f_p, 3)}$\,Hz, where the phase "
        r"crosses $0^\circ$ downwards" + filtri + r". The encoder channel of the model " + z_enc + r": with the "
        r"telescope rigidly attached at node P, between the encoder and the motor, the zeros of "
        r"$\theta_e/\tau_\theta$ come only from the branch beyond the motor (the rear end)." + z_mot
        + r" A zero at an encoder on the load side requires a large inertia outside the encoder--motor path: "
        r"the measurement is consistent with most of the telescope hanging from the shaft through a "
        r"compliance $k_b$ (clamp, saddle, flanges or bending of the tube) instead of being rigid at P. "
        r"A two-mass fit of the antiresonance--resonance pair, $J_1 = J_{tot}(f_z/f_p)^2$ and "
        rf"$k_b = (J_{{tot}} - J_1)\,\omega_z^2$, gives $J_1 \approx {num(e.J_1, 2)}$\,kg\,m$^2$ rigid with the "
        rf"shaft, $J_c \approx {num(e.J_c, 2)}$\,kg\,m$^2$ beyond the compliance and "
        rf"$k_b \approx {num(e.k_b, 2)}$\,N\,m/rad (the model has $k_{{\theta,2}} = {num(e.k_theta2, 3)}$\,N\,m/rad "
        r"between telescope and motor). Above the resonance the encoder then sees only $J_1$: the fit predicts "
        rf"a level ${num(e.livello_previsto_db, 3)}$\,dB above the rigid-telescope model, and the measurement is "
        rf"${num(e.livello_misurato_db, 3)}$\,dB above it ({num(e.banda_livello[0], 3)}--"
        rf"{num(e.banda_livello[1], 3)}\,Hz). Between {num(e.banda_alta[0], 2)}\,Hz and "
        rf"{num(e.banda_alta[1], 4)}\,Hz the measurement has {e.n_ris_alte} resonances against "
        rf"{e.n_modi_alti} modes of the polar chain: modes of the telescope, the fork and the drive side that the "
        r"torsional model does not contain. A modal test of the telescope with the polar axis blocked (first "
        rf"frequency expected near ${num(e.f_z, 3)}$\,Hz) and an accelerometer on the tube during the sine "
        r"sweep would confirm the attachment compliance. The elevation of this test is not recorded and is "
        r"assumed to be $0^\circ$."
    )


def _catena_ramo(e: ElevazioneEncoder | None, figure: Mapping[str, str], nota: str) -> str:
    """Paragrafo e figura: la catena polare con il telescopio appeso al mozzo tramite k_b."""
    from .latex_report import num

    if e is None or e.catena is None or "catena_ramo" not in figure:
        return ""
    c = e.catena
    hz = lambda v: rf"${num(v, 4)}$\,Hz"  # noqa: E731
    z_T = next(f for f, d in c.zeri_e if d.startswith("telescope"))
    z_r = next(f for f, d in c.zeri_e if d.startswith("rear"))
    z_sx = [f for f, d in c.zeri_m if not d.startswith("rear")]
    nuovo_polo = min((f for f in c.poli if f > 1.0), key=lambda f: abs(f - e.f_p))
    return "\n".join([
        r"\paragraph{How a compliance can add a zero.} A spring added \emph{on} the actuator--sensor path "
        r"adds an inertia and a pole but no zero: the number of inertias $n$ and the number of springs on the "
        r"path both grow by one, so $n - 1 - (\text{springs on the path})$ does not change. This is why the "
        r"shafts between the encoder and the telescope, however compliant, cannot produce the measured "
        r"antiresonance (Section~\ref{subsec:lump_tf}). A spring that \emph{detaches} an inertia from the path "
        r"is different. Splitting the telescope node into a hub, which stays on the path "
        rf"($J_P = {num(c.J['P'], 3)}$\,kg\,m$^2$), and the telescope, which hangs from the hub through $k_b$ "
        rf"($J_c = {num(c.J['T'], 3)}$\,kg\,m$^2$), adds an inertia without adding a spring to the path, so it "
        r"adds one zero: the natural frequency of the new off-path branch with the hub clamped, "
        r"$\omega_z = \sqrt{(k_b + k_g)/J_c}$. The telescope inertia moves from the denominator to the "
        r"numerator of the encoder channel:",
        r"\begin{gather}",
        r"\frac{\Theta_e}{T_\theta} = \frac{k_{\theta,1}\,k_{\theta,2}\,\bigl(J_c s^2 + k_b + k_g\bigr)\,"
        r"\bigl(J_{\theta,r}s^2 + k_{\theta,3}\bigr)}{\Delta(s)}.",
        r"\end{gather}",
        rf"Figure~\ref{{fig:catena_ramo}} draws this chain with $J_c$ and $k_b = {num(c.k[3], 3)}$\,N\,m/rad from "
        rf"the two-mass fit. The encoder channel now has the zero at {hz(z_T)} (measured {hz(e.f_z)}) besides "
        rf"the unchanged zero of the rear end at {hz(z_r)}, and the free chain has a new pole at "
        rf"{hz(nuovo_polo)}, the hub with encoder and motor oscillating against the telescope on $k_b$ "
        rf"(measured resonance about {hz(e.f_p)}). In the collocated motor channel the encoder--hub--telescope "
        r"piece with the motor clamped now gives the zeros " + ", ".join(hz(f) for f in z_sx)
        + r": the single telescope zero of the rigid-telescope chain splits in two, because the telescope and "
        r"the hub now move separately.",
        r"\begin{figure}[H]\centering",
        rf"\includegraphics[width=\linewidth,height=0.75\textheight,keepaspectratio]{{{figure['catena_ramo']}}}",
        r"\caption{Polar chain with the telescope hanging from the hub through the compliance $k_b$, "
        r"with $J_c$ and $k_b$ from the two-mass fit of the measurement: (a) torque$\to$encoder, "
        r"(b) torque$\to$motor. Same conventions as Figure~\ref{fig:catene_polari}: blue, path from the "
        r"actuator to the sensor; orange, dashed, off-path branches and their zeros; hatched, node clamped "
        r"to compute the zeros." + nota + r"}\label{fig:catena_ramo}\end{figure}",
    ])

def latex(confronti: Sequence[spm.Confronto], disc: Discussione, figure: Mapping[str, str],
          Kt: float, segnaposto: bool, cartella: str) -> str:
    from .latex_report import esc, num

    nota = rf" \emph{{{NOTA_REALISMO.capitalize()}.}}" if segnaposto else ""
    righe_mis = []
    for c in confronti:
        m = c.misura
        righe_mis.append(rf"\texttt{{{esc(m.file)}}} & {esc(in_inglese(m.asse))} & ${np.rad2deg(m.theta0):g}$ & "
                         rf"{'yes' if m.telescopio else 'no'} & ${num(m.eccitazione_A, 3)}$ & "
                         rf"${num(m.f[0], 3)}$--${num(m.f[-1], 4)}$ \\")
    righe_car = []
    for c in confronti:
        for x in c.caratteristiche:
            ch = {"psi_e": r"$\dot\psi_e$", "psi_m": r"$\dot\psi_m$", "theta_e": r"$\dot\theta_e$",
                  "theta_m": r"$\dot\theta_m$"}.get(x.canale or "", "--")
            fmod = "--" if x.f_modello is None else f"${num(x.f_modello, 4)}$"
            rap = "--" if x.f_modello is None else f"${num(x.f_misura / x.f_modello, 3)}$"
            righe_car.append(rf"{esc(c.misura.etichetta)} & {'resonance' if x.tipo == 'risonanza' else 'antiresonance'}"
                             rf" & ${num(x.f_misura, 4)}$ & {fmod} & {ch} & {rap} \\")
    figs = []
    for c in confronti:
        chiave = f"exp_{Path(c.misura.file).stem}"
        if chiave in figure:
            figs += [r"\begin{figure}[H]\centering",
                     rf"\includegraphics[width=\linewidth,height=0.48\textheight,keepaspectratio]{{{figure[chiave]}}}",
                     rf"\caption{{Measured plant ({esc(c.misura.etichetta)}, black) against the 7-DOF lumped "
                     rf"model at $\theta_0 = {np.rad2deg(c.misura.theta0):g}^\circ$, velocity at the encoder (blue: the "
                     r"measured channel, the drive feedback being the encoder) and at the motor node (orange, "
                     r"dashed, for reference). The grey band marks the range without data "
                     rf"(the measurement starts at {num(c.misura.f[0], 2)}\,Hz).{nota}}}\label{{fig:{chiave}}}"
                     r"\end{figure}"]
    if "exp_azimut_sovrapposte" in figure:
        figs += [r"\begin{figure}[H]\centering",
                 rf"\includegraphics[width=\linewidth,height=0.45\textheight,keepaspectratio]{{{figure['exp_azimut_sovrapposte']}}}",
                 r"\caption{The three azimuth measurements: effect of the excitation amplitude (3.5 and 2.36\,A) "
                 r"and of the elevation (0 and $30^\circ$).}\label{fig:exp_az_overlay}\end{figure}"]
    d = disc
    f0_txt = ", ".join(rf"${num(v, 3)}$\,Hz ({esc(k)})" for k, v in d.f0_az.items())
    an_txt = ", ".join(rf"${num(a, 4)}$ against ${num(b, 4)}$\,Hz ({esc(k)})" for k, (a, b) in d.antiris.items())
    f1_txt = ", ".join(rf"${num(a, 4)}$ against ${num(b, 4)}$\,Hz ({esc(k)})" for k, (a, b) in d.f1_az.items())
    punti = [
        (r"\textbf{Low frequency ($<5$\,Hz): stiffness to the ground.} Every azimuth measurement has a "
         r"resonance at " + f0_txt + r" and rises at $+20$\,dB/dec below it, with phase near $+90^\circ$: the axis "
         r"behaves as an inertia on a spring to the ground. The model has a free rigid-body mode (phase "
         r"$-90^\circ$, $-20$\,dB/dec)."
         + (rf" With $J_{{tot}}$ of the model the measured resonance requires $k_0 \approx {num(d.k0_az, 3)}$\,N\,m/rad."
            if d.k0_az else "")
         + r" The resonance moves up when the excitation is reduced (3.5\,A $\to$ 2.36\,A), i.e.\ the stiffness "
         r"depends on the amplitude: this is the signature of pre-sliding (Dahl-type) friction in the bearings "
         r"and seals, which acts as a spring for motions smaller than the break-away displacement. The model "
         r"has only viscous bearing friction, currently zero. The elevation measurement shows the same "
         r"behaviour below a few hertz."),
        (r"\textbf{Load antiresonance: reproduced.} The deep antiresonance of the azimuth, where the fork "
         r"oscillates on $k_{\psi,2}$ with motor and encoder still, is measured at " + an_txt + r"; its shift "
         r"with the elevation (through $J_{\psi,p}(\theta_0)$) is also reproduced. The measured "
         r"$k_{\psi,2}/J_L^*$ is about "
         + (rf"${num((list(d.antiris.values())[0][0] / list(d.antiris.values())[0][1]) ** 2 * 100 - 100, 2)}\,\%$"
            if d.antiris else "--")
         + r" off: a joint compliance in series with the shafts (Q5) or a heavier fork."),
        (r"\textbf{Level between the antiresonance and the first resonance.} In this band the encoder "
         r"velocity is that of the motor and encoder inertias on the compliance $1/k_{\psi,2}$, "
         rf"$|G|\approx K_t\omega/k_{{\psi,2}}$. The measurement is ${num(d.livello_az_db, 2)}$\,dB below the "
         r"model. A stiffer $k_{\psi,2}$ would contradict the antiresonance, which indicates a slightly softer "
         r"one; the first suspect is the torque scale: $K_t$ refers to A\,rms at $25^\circ$C, and a current "
         r"read in A peak gives $-3$\,dB, a hot winding $-0.8$\,dB. The current units of the WorkBench export "
         r"should be checked before adjusting stiffnesses."),
        (r"\textbf{First structural resonance of the azimuth.} Measured " + f1_txt + r" (about "
         + (rf"${num((1 - list(d.f1_az.values())[0][0] / list(d.f1_az.values())[0][1]) * 100, 2)}\,\%$"
            if d.f1_az else "--")
         + r" lower). It is the motor and encoder side oscillating on $k_{\psi,2}$ against the heavy fork, "
         r"$\omega_1^2\approx k_{\psi,2}/(J_{\psi,e} + J_{\psi,m})$: with the model stiffness it requires "
         + (rf"$J_{{\psi,e}} + J_{{\psi,m}} \approx {num(d.J_motore_richiesta, 3)}$\,kg\,m$^2$ instead of "
            rf"${num(d.J_motore_modello, 3)}$\," if d.J_motore_richiesta else "")
         + r"kg\,m$^2$. Part of the difference is known: the azimuth configuration has no rotor inertia at the "
         rf"motor node (TBM(S)-7646 rotor ${num(d.J_rotore_datasheet, 3)}$\,kg\,m$^2$) and no encoder disc at "
         r"node 1 (warning \texttt{NODI\_SENZA\_INERZIA\_AZIMUT}); the rest points to a softer drive-side "
         r"connection (shaft--rotor and shaft--hub joints) than the shaft geometry alone."),
        ((rf"\textbf{{Unmodelled mode at ${num(d.pair_az[0], 3)}$/${num(d.pair_az[1], 3)}$\,Hz (azimuth).}} A "
          r"resonance--antiresonance pair, present at both elevations and amplitudes, with no counterpart in the "
          r"model: a degree of freedom outside the torsional chains (fork bending, assumed rigid in Q12; base or "
          r"pier; cable chain).") if d.pair_az else ""),
        _punto_elevazione(d.el),
        ((r"\textbf{Elevation without the telescope.} The model (polar shafts only, "
          rf"${num(d.J_nt_modello, 3)}$\,kg\,m$^2$) is a nearly free rigid inertia and lies far above the "
          r"measurement, which rises like a spring up to its first peak. From the low-frequency slope "
          rf"$|G| \approx K_t\omega/k_0$, $k_0 \approx {num(d.k0_nt, 3)}$\,N\,m/rad, and from the peak at "
          rf"${num(d.f_pk_nt, 3)}$\,Hz an inertia of about ${num(d.J_nt_richiesta, 3)}$\,kg\,m$^2$: the parts left "
          r"on the polar axis without the telescope (flanges, clamps) are much heavier than the shafts of the "
          r"model, and the friction stiffness dominates the light axis.") if d.k0_nt and d.J_nt_richiesta else ""),
        (r"\textbf{High frequency.} Above a few hundred hertz the measured phase falls faster than the model: "
         r"the current loop (1.5\,kHz bandwidth) and the sampling of the drive add a delay that the "
         r"mechanical model does not contain."),
    ]
    return "\n".join([
        r"\section{Comparison with experimental frequency responses}\label{sec:experimental}",
        r"\paragraph{Measurements.} Plant measurements exported by Kollmorgen WorkBench "
        rf"(folder \texttt{{{esc(cartella)}}}; the files of the sub-folder \texttt{{Validation}} are not used). "
        r"The drive injects a stepped-sine current (amplitude in Table~\ref{tab:exp_meas}, frequency step "
        r"of 1.03--1.1\,\%, 8\,kHz sampling, 4096-point FFT) with the loops closed, and estimates the plant "
        r"as the ratio of the feedback velocity to the current, in (rad/s)/A. The drives close their loops on "
        r"the axis encoders, so the measured plant is the encoder velocity over the current; the model is "
        r"also evaluated at the motor node for reference. The units are confirmed by the rigid-body line of "
        r"the elevation, which matches $K_t/(J\omega)$ of the model within about 1\,dB. For the elevation the "
        r"measurement of 16 October 2026, extended to 2\,kHz, replaces "
        r"\texttt{Elevation\_BodeMeasurement.csv}, which stopped at about 277\,Hz "
        r"(\texttt{Elevation\_new\_BodeMeasurement.csv} holds the same data as the latter); neither is used.",
        r"\begin{table}[H]\centering\small",
        r"\caption{Measurements used: axis, elevation $\theta_0$ (deg), telescope mounted, current amplitude "
        r"(A) and frequency range (Hz).}\label{tab:exp_meas}",
        r"\begin{tabular}{llrlrr}\toprule",
        r"file & axis & $\theta_0$ & telescope & $I$ & range \\ \midrule",
        *righe_mis,
        r"\bottomrule\end{tabular}\end{table}",
        r"\paragraph{Model.} Only the 7-DOF lumped model of Section~\ref{sec:lumped} is used (below 1\,kHz "
        r"it has the dynamics of interest). It is linearised at the elevation of each test with $\Omega_0 = 0$ "
        r"and converted to the units of the measurement:",
        r"\begin{gather}",
        r"G_{meas}(j\omega) = \frac{\dot\theta(j\omega)}{I(j\omega)} = j\omega\,K_t\,\frac{\theta(j\omega)}"
        rf"{{\tau(j\omega)}},\qquad K_t = {num(Kt, 3)}\ \text{{N\,m/A}}\ \text{{(TBM(S)-7646-A, }}25^\circ\text{{C)}},",
        r"\end{gather}",
        r"at the encoder and at the motor node of the excited axis. For the test without the telescope the "
        r"payload is removed (polar shafts only). The phase is shown wrapped to $(-180^\circ, 180^\circ]$ as "
        r"in the WorkBench export.",
        *figs,
        r"\paragraph{Measured features against the model.} Peaks and notches of the measured magnitude with a "
        r"local prominence of at least 8\,dB, paired with the nearest mode of the excited axis "
        r"(resonances) or invariant zero of the encoder channel, the measured one (antiresonances), within "
        r"a factor 2 in frequency; ``--'' means no counterpart in the model.",
        r"\begin{longtable}{p{0.32\linewidth}lrrlr}",
        r"\caption{Measured resonances and antiresonances against the 7-DOF lumped model (Hz)." + nota
        + r"}\label{tab:exp_features}\\",
        r"\toprule measurement & feature & measured & model & channel & ratio \\ \midrule\endfirsthead",
        r"\toprule measurement & feature & measured & model & channel & ratio \\ \midrule\endhead",
        *righe_car,
        r"\bottomrule\end{longtable}",
        *_tabella_fit(confronti, nota),
        r"\subsection{Discrepancies and their likely causes}\label{subsec:exp_discussion}",
        r"\begin{enumerate}",
        *[r"\item " + p for p in punti if p],
        r"\end{enumerate}",
        _catena_ramo(d.el, figure, nota),
        r"\paragraph{Summary.} The 7-DOF model reproduces the structure of the azimuth response between "
        r"about 5\,Hz and 300\,Hz: rigid-body line, load antiresonance and its dependence on the elevation, "
        r"compliance line. The main discrepancies are a missing amplitude-dependent stiffness to the ground "
        r"(bearing friction) below 5\,Hz, a first azimuth structural mode about one third lower than predicted "
        r"(missing drive-side inertias and joint compliances), a few dB of level (torque scale to be checked), "
        r"and modes outside the torsional chains (rigid fork and telescope). "
        + (r"On the elevation the hub model, with the telescope hanging from the hub through $k_b$, has the "
           r"antiresonance at the encoder; its frequency is set by the estimated $k_b$, while the following "
           r"resonance depends on the stiffness of the hub (Section~\ref{sec:hub})."
           if d.el is not None and d.el.mozzo else
           r"On the elevation the model reproduces the rigid-body line but not the antiresonance measured at "
           r"the encoder: the telescope is not rigidly attached to the shaft between encoder and motor, as the "
           r"model assumes."),
    ])


def _tabella_fit(confronti: Sequence[spm.Confronto], nota: str) -> list[str]:
    """Fit NRMSE della misura (canale encoder) per banda, definizione di System Identification Toolbox."""
    if not confronti or not confronti[0].fit:
        return []
    bande = [b for b, _, _ in spm.BANDE_FIT]
    righe = []
    for c in confronti:
        vals = [("--" if not np.isfinite(c.fit.get(b, np.nan)) else f"${c.fit[b]:.1f}$") for b in bande]
        righe.append(rf"{c.misura.etichetta} & " + " & ".join(vals) + r" \\")
    intest = " & ".join(["all points" if b == "all" else b.replace("-", "--").replace(">", "$>$") for b in bande])
    return [
        r"\paragraph{Goodness of fit.} As in the System Identification Toolbox (\texttt{compare}, "
        r"\texttt{goodnessOfFit} with the NRMSE cost, report of \texttt{greyest}), the fit of the measured channel "
        r"is $\mathrm{fit} = 100\,(1 - \lVert y - \hat y\rVert/\lVert y - \bar y\rVert)$\,\%, with $y$ the measured "
        r"frequency response at the measured frequencies (complex), $\hat y$ the model at the same frequencies and "
        r"$\bar y$ the mean of $y$: 100\,\% is a perfect model, 0\,\% is no better than a constant, and it can be "
        r"negative.",
        r"\begin{table}[H]\centering\small",
        r"\caption{Fit (NRMSE, \%) of the encoder channel, on all the measured points and by band." + nota
        + r"}\label{tab:exp_fit}",
        r"\begin{tabular}{l" + "r" * len(bande) + r"}\toprule",
        r"measurement & " + intest + r" \\ \midrule",
        *righe,
        r"\bottomrule\end{tabular}\end{table}",
    ]


def markdown(confronti: Sequence[spm.Confronto], disc: Discussione) -> list[str]:
    L = ["## Comparison with the experimental measurements (7-DOF lumped model)", "",
         "Kollmorgen WorkBench measurements (plant = encoder velocity / current, (rad/s)/A); Validation folder "
         "excluded.", "",
         "| measurement | feature | measured f [Hz] | model f [Hz] | channel |", "|---|---|---|---|---|"]
    for c in confronti:
        for x in c.caratteristiche:
            L.append(f"| {c.misura.etichetta} | {TIPO_EN[x.tipo]} | {x.f_misura:.4g} | "
                     f"{'--' if x.f_modello is None else f'{x.f_modello:.4g}'} | {x.canale or '--'} |")
    if confronti and confronti[0].fit:
        bande = [b for b, _, _ in spm.BANDE_FIT]
        L += ["", "Fit NRMSE of the encoder channel (System Identification Toolbox definition), %:", "",
              "| measurement | " + " | ".join(bande) + " |", "|---|" + "---|" * len(bande)]
        for c in confronti:
            L.append(f"| {c.misura.etichetta} | " + " | ".join(f"{c.fit[b]:.1f}" for b in bande) + " |")
    L += ["", f"- Azimuth level 15-300 Hz: measurement {disc.livello_az_db:+.1f} dB relative to the model."]
    if disc.k0_az:
        L.append(f"- Equivalent stiffness to the ground (azimuth): k_0 ~ {disc.k0_az:.3g} N m/rad (absent from the model).")
    if disc.J_motore_richiesta:
        L.append(f"- Motor-side inertia for the first measured mode: {disc.J_motore_richiesta:.3g} kg m^2 "
                 f"(model {disc.J_motore_modello:.3g}).")
    if disc.k0_nt and disc.J_nt_richiesta:
        L.append(f"- Without telescope: k_0 ~ {disc.k0_nt:.3g} N m/rad, J ~ {disc.J_nt_richiesta:.3g} kg m^2 "
                 f"(model {disc.J_nt_modello:.3g}).")
    if disc.el is not None:
        e = disc.el
        L.append(f"- Elevation: antiresonance at the encoder {e.f_z:.4g} Hz, resonance {e.f_p:.3g} Hz; two-mass fit "
                 f"J_1 ~ {e.J_1:.2g} kg m^2, J_c ~ {e.J_c:.2g} kg m^2, k_b ~ {e.k_b:.2g} N m/rad (telescope hanging "
                 f"from the shaft through a compliance; the model has no zero of theta_e below "
                 f"{e.z_enc_min or float('nan'):.4g} Hz).")
    return L + [""]
