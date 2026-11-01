"""Figure del modello multi-GdL (sezione 11).

Lo stile e' quello del modello a 2 GdL (:mod:`gimbal.plotting`, riusato): stessa
palette categorica in ordine fisso, colore mai unica codifica (tratto diverso
per ogni serie), PDF vettoriale e PNG a 300 dpi. Testi delle figure in inglese, Hz e
gradi. Ogni figura che dipende da dati segnaposto porta l'avviso nel titolo.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from gimbal.plotting import (GRIGLIA, INCHIOSTRO, INCHIOSTRO_SECONDARIO, INCHIOSTRO_TENUE, SERIE,
                             TRATTI, applica_stile)
from gimbal.plotting import salva as _salva_2dof
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, NullLocator

from . import fem
from . import linearize as lz
from . import tf
from .analysis import SweepTheta
from .coupled import Sistema
from .params import MM, NOME_ASSE, NOTA_REALISMO, Parametri
from .simulate import RisultatoScenario

__all__ = ["applica_stile", "salva"]

GRIGLIA_MINORE = "#efeeeb"
"""Griglia secondaria, piu' tenue della principale :data:`GRIGLIA`."""


def rifinisci_assi(fig: Figure) -> None:
    """Griglia principale e secondaria, valori ed etichetta con unita' su ogni asse.

    Con ``sharex`` matplotlib nasconde i valori e l'etichetta dell'asse x dei
    pannelli interni: qui si rimettono, copiando l'etichetta del pannello della
    stessa colonna che la porta. Gli assi con ``gid = "disegno"`` (disegni in
    scala, non grafici) restano senza griglia.
    """
    etichette_x: dict[tuple[int, int], str] = {}
    for ax in fig.axes:
        spec = ax.get_subplotspec()
        if spec is not None and ax.get_xlabel():
            etichette_x[(spec.colspan.start, spec.colspan.stop)] = ax.get_xlabel()
    for ax in fig.axes:
        if ax.get_gid() == "disegno":
            continue
        ax.minorticks_on()
        if isinstance(ax.xaxis.get_major_locator(), FixedLocator):
            ax.xaxis.set_minor_locator(NullLocator())  # tacche con nome (nodi): niente secondarie
        ax.grid(True, which="major", color=GRIGLIA, lw=0.7)
        ax.grid(True, which="minor", color=GRIGLIA_MINORE, lw=0.45, ls="-")
        ax.tick_params(which="both", labelbottom=True, labelleft=True)
        spec = ax.get_subplotspec()
        if spec is not None and not ax.get_xlabel():
            ax.set_xlabel(etichette_x.get((spec.colspan.start, spec.colspan.stop), ""))
    fig.tight_layout()


def salva(fig: Figure, cartella: str | Path, nome: str) -> list[Path]:
    """Rifinisce gli assi (:func:`rifinisci_assi`) e salva in PDF e PNG (stile del 2 GdL)."""
    rifinisci_assi(fig)
    return _salva_2dof(fig, cartella, nome)


def _titolo(fig: Figure, testo: str, segnaposto: bool) -> None:
    if segnaposto:
        testo = f"{testo}\n[{NOTA_REALISMO.upper()}]"
    fig.suptitle(testo, fontsize=10, color=INCHIOSTRO)


def _stile(i: int) -> dict[str, object]:
    return {"color": SERIE[i % len(SERIE)], "linestyle": TRATTI[i % len(TRATTI)]}


def _db(G: np.ndarray) -> np.ndarray:
    return 20 * np.log10(np.maximum(np.abs(G), 1e-300))


def _fase(G: np.ndarray) -> np.ndarray:
    return np.rad2deg(np.unwrap(np.angle(G)))


def _assi_bode(ax_m: plt.Axes, ax_f: plt.Axes, titolo: str) -> None:
    ax_m.set_title(titolo)
    ax_m.set_ylabel("magnitude [dB rad/(N m)]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_xlabel("frequency [Hz]")
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")


# ---------------------------------------------------------------------------
# Bode 2x2 coppia -> encoder
# ---------------------------------------------------------------------------


def bode_2x2(modelli: Sequence[tuple[str, lz.ModelloLineare]], w: np.ndarray, *, titolo: str,
             segnaposto: bool, uscite: tuple[str, str] = ("psi_e", "theta_e"),
             segna_radici: bool = True) -> Figure:
    """Matrice 2x2 coppia -> uscite: per ogni canale modulo e fase.

    Sul primo modello si segnano i poli (linee verticali tratteggiate) e gli
    zeri invarianti (cerchi sul modulo) delle singole FdT SISO.
    """
    f = w / (2 * np.pi)
    fig, assi = plt.subplots(4, 2, figsize=(10.0, 10.0), sharex=True)
    ingressi = (r"$\tau_\psi$", r"$\tau_\theta$")
    fasi_rif: dict[tuple[int, int], np.ndarray] = {}
    for k, (nome, ml) in enumerate(modelli):
        G = ml.frf(w, list(uscite))
        for i in range(2):
            for j in range(2):
                ax_m, ax_f = assi[2 * i, j], assi[2 * i + 1, j]
                ax_m.plot(f, _db(G[:, i, j]), label=nome, **_stile(k))
                fase = _fase(G[:, i, j])
                # a un polo non smorzato il salto di +/-180 gradi e' ambiguo e unwrap puo'
                # scegliere rami diversi per modelli diversi: si allinea al primo modello
                if (i, j) in fasi_rif:
                    fase = fase - 360.0 * np.round(np.median(fase - fasi_rif[(i, j)]) / 360.0)
                else:
                    fasi_rif[(i, j)] = fase
                ax_f.plot(f, fase, **_stile(k))
    ml0 = modelli[0][1]
    for i, u in enumerate(uscite):
        for j in range(2):
            ax_m, ax_f = assi[2 * i, j], assi[2 * i + 1, j]
            _assi_bode(ax_m, ax_f, f"{lz.ETICHETTE_USCITE[u]} / {ingressi[j]}")
            if not segna_radici:
                continue
            fp = np.abs(ml0.poli[np.imag(ml0.poli) > 0]) / (2 * np.pi)
            for v in fp[(fp > f[0]) & (fp < f[-1])]:
                ax_m.axvline(v, color=INCHIOSTRO_TENUE, lw=0.6, ls=":")
            z = tf.zeri_invarianti(ml0.A, ml0.Bu[:, j], ml0.Cy[ml0.riga(u)])
            fz = np.abs(z[np.imag(z) > 0]) / (2 * np.pi)
            fz = fz[(fz > f[0]) & (fz < f[-1])]
            if fz.size:
                Gz = ml0.frf(2 * np.pi * fz, [u])[:, 0, j]
                ax_m.plot(fz, _db(Gz), "o", mfc="none", color=INCHIOSTRO, ms=4,
                          label="zeros" if (i, j) == (0, 0) else None)
    assi[0, 0].legend(loc="lower left")
    _titolo(fig, titolo, segnaposto)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Confronto completo/ridotto e concentrata/consistente (asse isolato)
# ---------------------------------------------------------------------------


def bode_confronto_asse(sistemi: Mapping[str, Sistema], asse: str, theta0: float, w: np.ndarray,
                        *, segnaposto: bool) -> Figure:
    """FRF coppia -> encoder e coppia -> motore dell'asse isolato:
    completo e ridotto, massa concentrata e consistente."""
    from .analysis import frf_encoder_isolato

    f = w / (2 * np.pi)
    fig, assi = plt.subplots(2, 2, figsize=(10.0, 6.5), sharex=True)
    k = 0
    for mm in ("lumped", "consistent"):
        s = sistemi[mm]
        mod = s.asse_psi if asse == "azimut" else s.asse_theta
        rid = s.rid_psi if asse == "azimut" else s.rid_theta
        M, K, C = s.isolato_completo(asse, theta0)
        Mr, Kr, Cr = s.isolato_ridotto(asse, theta0)
        for nome, (MM_, KK, CC, im, ie) in (
            (f"full, {mm}", (M, K, C, mod.i_motore, mod.i_encoder)),
            (f"reduced ({s.forma_massa}), {mm}", (Mr, Kr, Cr, rid.i_motore, rid.i_encoder)),
        ):
            for col, uscita in enumerate((ie, im)):
                G = frf_encoder_isolato(MM_, CC, KK, im, uscita, w)
                assi[0, col].plot(f, _db(G), label=nome, **_stile(k))
                assi[1, col].plot(f, _fase(G), **_stile(k))
            k += 1
    simbolo = r"\psi" if asse == "azimut" else r"\theta"
    _assi_bode(assi[0, 0], assi[1, 0], f"torque -> encoder ${simbolo}_e$ (non-collocated)")
    _assi_bode(assi[0, 1], assi[1, 1], f"torque -> motor ${simbolo}_m$ (collocated)")
    assi[0, 0].legend(loc="lower left", fontsize=7)
    _titolo(fig, f"Isolated {NOME_ASSE[asse]} axis, theta0 = {np.rad2deg(theta0):g} deg: full vs reduced",
            segnaposto)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Forme modali
# ---------------------------------------------------------------------------


def forme_modali(s: Sistema, asse: str, theta0: float, *, n_modi: int = 4, segnaposto: bool) -> Figure:
    """Forme modali del modello completo isolato, con i master evidenziati, e
    forme del ridotto riespanse con :math:`\\mathbf T` e sovrapposte."""
    mod = s.asse_psi if asse == "azimut" else s.asse_theta
    rid = s.rid_psi if asse == "azimut" else s.rid_theta
    M, K, _ = s.isolato_completo(asse, theta0)
    Mr, Kr, _ = s.isolato_ridotto(asse, theta0)
    w2, Phi = fem.modi(M, K)
    w2r, Phir = fem.modi(Mr, Kr)
    x = mod.x / MM
    # il moto rigido (autovalore nullo a meno dell'arrotondamento) non e' una forma modale
    scala = float(np.max(np.abs(w2)))
    indici = [k for k in range(len(w2)) if abs(w2[k]) > 1e-12 * scala][:n_modi]
    fig, assi = plt.subplots(len(indici), 1, figsize=(8.0, 2.0 * len(indici)), sharex=True)
    assi = np.atleast_1d(assi)
    for ax, k in zip(assi, indici):
        phi = Phi[:, k] / np.max(np.abs(Phi[:, k]))
        f_k = np.sign(w2[k]) * np.sqrt(abs(w2[k])) / (2 * np.pi)
        ax.plot(x, phi, marker="o", ms=3, label=f"full, f = {f_k:.4g} Hz", **_stile(0))
        if k < len(w2r):
            phir = rid.T @ Phir[:, k]
            phir = phir / np.max(np.abs(phir)) * (np.sign(phir @ phi) or 1.0)
            f_r = np.sign(w2r[k]) * np.sqrt(abs(w2r[k])) / (2 * np.pi)
            ax.plot(x, phir, label=f"reduced, re-expanded, f = {f_r:.4g} Hz", **_stile(1))
        ax.plot(x[list(mod.masters)], phi[list(mod.masters)], "s", ms=7, mfc="none",
                color=INCHIOSTRO, label="master nodes" if k == indici[0] else None)
        ax.set_ylabel(f"mode {k + 1}\nnormalised shape [-]")
        ax.legend(loc="best", fontsize=7)
        ax.axhline(0.0, color=GRIGLIA, lw=0.8)
    assi[-1].set_xticks(x)
    assi[-1].set_xticklabels(mod.etichette, fontsize=7)
    assi[-1].set_xlabel("node (position along the shaft, to scale)")
    _titolo(fig, f"Mode shapes, isolated {NOME_ASSE[asse]} axis, theta0 = {np.rad2deg(theta0):g} deg "
                 f"({s.modello_massa}, {s.forma_massa})", segnaposto)
    fig.tight_layout()
    return fig


def forme_modali_modale(s_mod, asse: str, *, n_modi: int = 4, segnaposto: bool) -> Figure:
    """Modi flessibili del FEM completo (base del troncamento modale, payload congelato in
    theta_ref), con l'indicazione tenuto/scartato e i nodi encoder, motore e d'aggancio."""
    an = s_mod.analisi(asse)
    ast = s_mod.asse(asse)
    mod = ast.modello
    x = mod.x / MM
    indici = list(range(1, min(n_modi, an.Phi.shape[1] - 1) + 1))
    fig, assi = plt.subplots(len(indici), 1, figsize=(8.0, 2.0 * len(indici)), sharex=True)
    assi = np.atleast_1d(assi)
    nodi = {"encoder": ast.i_e, "motor": ast.i_m, "attachment " + ("L" if asse == "azimut" else "P"): ast.i_agg}
    for ax, k in zip(assi, indici):
        phi = an.Phi[:, k] / np.max(np.abs(an.Phi[:, k]))
        stato = "retained" if k in an.tenuti else "discarded"
        ax.plot(x, phi, marker="o", ms=3, label=f"mode {k}, f = {an.f[k]:.4g} Hz ({stato})", **_stile(0))
        for j, (nome, i) in enumerate(nodi.items()):
            ax.plot(x[i], phi[i], ["s", "D", "^"][j], ms=7, mfc="none", color=INCHIOSTRO,
                    label=nome if k == indici[0] else None)
        ax.set_ylabel(f"mode {k}\nnormalised shape [-]")
        ax.legend(loc="best", fontsize=7)
        ax.axhline(0.0, color=GRIGLIA, lw=0.8)
    assi[-1].set_xticks(x)
    assi[-1].set_xticklabels(mod.etichette, fontsize=7)
    assi[-1].set_xlabel("node (position along the shaft, to scale)")
    _titolo(fig, f"Modes of the full FEM, {NOME_ASSE[asse]} axis (payload frozen at theta_ref = "
                 f"{np.rad2deg(an.theta_ref):g} deg)", segnaposto)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Sensibilita'
# ---------------------------------------------------------------------------


def frequenze_vs_theta0(sw: SweepTheta, *, segnaposto: bool, n_modi: int = 6) -> Figure:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    th = np.rad2deg(sw.theta0)
    for k, Om in enumerate(sw.omega0):
        F = sw.frequenze[Om]
        for m in range(min(n_modi, F.shape[1])):
            ax.plot(th, F[:, m], label=f"Omega0 = {np.rad2deg(Om):g} deg/s" if m == 0 else None,
                    **_stile(k))
    ax.set_yscale("log")
    ax.set_xlabel(r"$\theta_0$ = elevation [deg] ($\theta > 0$ raises the pointing)")
    ax.set_ylabel("frequency of the complex modes [Hz]")
    ax.legend()
    _titolo(fig, "Natural frequencies of the coupled system against theta0", segnaposto)
    fig.tight_layout()
    return fig


def frf_incrociate_omega(G: Mapping[float, np.ndarray], w: np.ndarray, theta0: float, *,
                         segnaposto: bool) -> Figure:
    f = w / (2 * np.pi)
    fig, assi = plt.subplots(1, 2, figsize=(10.0, 4.0), sharex=True)
    for k, (Om, Gm) in enumerate(G.items()):
        assi[0].plot(f, _db(Gm[:, 0, 1]), label=f"Omega0 = {np.rad2deg(Om):g} deg/s", **_stile(k))
        assi[1].plot(f, _db(Gm[:, 1, 0]), **_stile(k))
    assi[0].set_title(r"$\psi_e / \tau_\theta$")
    assi[1].set_title(r"$\theta_e / \tau_\psi$")
    for ax in assi:
        ax.set_xscale("log")
        ax.set_xlabel("frequency [Hz]")
        ax.set_ylabel("magnitude [dB rad/(N m)]")
    assi[0].legend()
    _titolo(fig, f"Cross FRFs torque -> encoder, theta0 = {np.rad2deg(theta0):g} deg",
            segnaposto)
    fig.tight_layout()
    return fig


def bode_errore_misura(ml: lz.ModelloLineare, w: np.ndarray, *, segnaposto: bool) -> Figure:
    """Coppia -> puntamento reale, encoder ed errore di misura, per ciascun asse."""
    f = w / (2 * np.pi)
    fig, assi = plt.subplots(2, 2, figsize=(10.0, 6.5), sharex=True)
    gruppi = ((0, ("psi_L", "psi_e", "err_psi")), (1, ("theta_p", "theta_e", "err_theta")))
    for col, (j, uscite) in enumerate(gruppi):
        G = ml.frf(w, list(uscite))
        for k, u in enumerate(uscite):
            assi[0, col].plot(f, _db(G[:, k, j]), label=lz.ETICHETTE_USCITE[u], **_stile(k))
            assi[1, col].plot(f, _fase(G[:, k, j]), **_stile(k))
        ingresso = r"$\tau_\psi$" if j == 0 else r"$\tau_\theta$"
        _assi_bode(assi[0, col], assi[1, col], f"response to {ingresso}")
        assi[0, col].legend(loc="lower left", fontsize=7)
    _titolo(fig, f"Measurement error encoder - true angle, theta0 = {np.rad2deg(ml.theta0):g} deg, "
                 f"Omega0 = {np.rad2deg(ml.Omega0):g} deg/s", segnaposto)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Simulazioni
# ---------------------------------------------------------------------------


def storie_temporali(r: RisultatoScenario, *, segnaposto: bool) -> Figure:
    canali = (("psi_e", "psi_L"), ("theta_e", "theta_p"))
    fig, assi = plt.subplots(2, 1, figsize=(8.0, 6.0), sharex=True)
    for ax, gruppo in zip(assi, canali):
        for k, u in enumerate(gruppo):
            i = lz.USCITE.index(u)
            ax.plot(r.t * 1e3, np.rad2deg(r.y_nl[:, i]) * 3600, label=f"{lz.ETICHETTE_USCITE[u]} nonlinear",
                    color=SERIE[k], linestyle="-")
            if r.y_lin is not None:
                ax.plot(r.t * 1e3, np.rad2deg(r.y_lin[:, i]) * 3600, label="linear",
                        color=SERIE[k], linestyle="--", lw=1.0)
        ax.set_ylabel("deviation from equilibrium [arcsec]")
        ax.legend(fontsize=7, loc="best")
    assi[-1].set_xlabel("t [ms]")
    _titolo(fig, f"{r.nome}: {r.descrizione} (theta0 = {np.rad2deg(r.theta0):g} deg, "
                 f"Omega0 = {np.rad2deg(r.Omega0):g} deg/s)", segnaposto)
    fig.tight_layout()
    return fig


def storie_grande_ampiezza(r: RisultatoScenario, *, segnaposto: bool) -> Figure:
    """Scenario a grande ampiezza (S4): deviazioni in gradi, non lineare contro lineare."""
    fig, assi = plt.subplots(2, 1, figsize=(8.0, 6.0), sharex=True)
    for ax, (u, nome) in zip(assi, (("theta_p", r"$\theta_p - \theta_0$ (telescope)"),
                                    ("psi_L", r"$\psi_L - \psi_0 - \Omega_0 t$ (fork)"))):
        i = lz.USCITE.index(u)
        ax.plot(r.t, np.rad2deg(r.y_nl[:, i]), label="nonlinear (Newton-Euler)", color=SERIE[0], ls="-")
        if r.y_lin is not None:
            ax.plot(r.t, np.rad2deg(r.y_lin[:, i]), label="linearised", color=SERIE[1], ls="--", lw=1.2)
        ax.set_ylabel(f"{nome} [deg]")
        ax.legend(fontsize=7, loc="best")
    assi[-1].set_xlabel("t [s]")
    _titolo(fig, f"{r.nome}: {r.descrizione} (theta0 = {np.rad2deg(r.theta0):g} deg, "
                 f"Omega0 = {np.rad2deg(r.Omega0):g} deg/s)", segnaposto)
    fig.tight_layout()
    return fig


ETICHETTE_CANALI: dict[str, str] = {
    "psi_e": r"7-DOF model, encoder $\dot\psi_e$", "psi_m": r"7-DOF model, motor $\dot\psi_m$",
    "theta_e": r"7-DOF model, encoder $\dot\theta_e$", "theta_m": r"7-DOF model, motor $\dot\theta_m$",
}


def _fase_avvolta(G: np.ndarray) -> np.ndarray:
    """Fase in (-180, 180] gradi, come nelle esportazioni di WorkBench."""
    return np.rad2deg(np.angle(G))


def bode_misura(conf, *, segnaposto: bool, f_lim: tuple[float, float] = (0.1, 2000.0)) -> Figure:
    """Pianta misurata (velocita'/corrente) contro il modello a 7 GdL ai nodi encoder e motore."""
    mis = conf.misura
    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(9.0, 6.8), sharex=True)
    for k, (c, G) in enumerate(conf.G_modello.items()):
        st = _stile(k)
        ax_m.plot(conf.f_modello, _db(G), label=ETICHETTE_CANALI.get(c, c), lw=1.5, **st)
        ax_f.plot(conf.f_modello, _fase_avvolta(G), lw=1.2, **st)
    ax_m.plot(mis.f, mis.modulo_db, color=INCHIOSTRO, marker="o", ms=2.5, lw=1.0,
              label=f"measurement ({mis.eccitazione_A:g} A)")
    ax_f.plot(mis.f, mis.fase_deg, color=INCHIOSTRO, marker="o", ms=2.5, lw=1.0)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(*f_lim)
        ax.axvspan(f_lim[0], mis.f[0], color=GRIGLIA_MINORE, zorder=0)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower right", fontsize=7)
    _titolo(fig, f"Measured plant against model: {mis.etichetta} (theta0 = {np.rad2deg(mis.theta0):g} deg)",
            segnaposto)
    fig.tight_layout()
    return fig


def bode_misure_sovrapposte(misure: Sequence, *, titolo: str, segnaposto: bool,
                            f_lim: tuple[float, float] = (0.1, 2000.0)) -> Figure:
    """Piu' misure dello stesso asse sovrapposte (effetto di ampiezza e di elevazione)."""
    fig, (ax_m, ax_f) = plt.subplots(2, 1, figsize=(9.0, 6.8), sharex=True)
    colori = [INCHIOSTRO] + list(SERIE[2:])
    for k, m in enumerate(misure):
        st = {"color": colori[k % len(colori)], "linestyle": TRATTI[k % len(TRATTI)]}
        ax_m.plot(m.f, m.modulo_db, marker="o", ms=2.0, lw=1.0, label=m.etichetta, **st)
        ax_f.plot(m.f, m.fase_deg, marker="o", ms=2.0, lw=1.0, **st)
    for ax in (ax_m, ax_f):
        ax.set_xscale("log")
        ax.set_xlim(*f_lim)
    ax_m.set_ylabel("magnitude [dB (rad/s)/A]")
    ax_f.set_ylabel("phase [deg]")
    ax_f.set_ylim(-180, 180)
    ax_f.set_yticks(np.arange(-180, 181, 90))
    ax_f.set_xlabel("frequency [Hz]")
    ax_m.legend(loc="lower right", fontsize=7)
    _titolo(fig, titolo, segnaposto)
    fig.tight_layout()
    return fig


def conservazione(r: RisultatoScenario, *, segnaposto: bool) -> Figure:
    assert r.energia is not None and r.momento is not None
    fig, assi = plt.subplots(2, 1, figsize=(8.0, 5.0), sharex=True)
    assi[0].plot(r.t * 1e3, (r.energia - r.energia[0]) / abs(r.energia[0]), color=SERIE[0])
    assi[0].set_ylabel(r"$(E - E_0)/|E_0|$ [-]")
    assi[1].plot(r.t * 1e3, (r.momento - r.momento[0]) / abs(r.momento[0]), color=SERIE[1], ls="--")
    assi[1].set_ylabel(r"$(p_\psi - p_{\psi,0})/|p_{\psi,0}|$ [-]")
    assi[1].set_xlabel("t [ms]")
    _titolo(fig, f"S3: conservation with C = 0, tau = 0 (drifts {r.derive['energia']:.1e}, "
                 f"{r.derive['momento_psi']:.1e})", segnaposto)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Catene polari: da dove vengono poli e zeri (coppia -> encoder, coppia -> motore)
# ---------------------------------------------------------------------------


def _molla(ax: plt.Axes, x0: float, x1: float, y: float, testo: str, colore: str) -> None:
    """Molla a zig-zag fra x0 e x1 con l'etichetta sopra."""
    n, h = 6, 0.14
    xs = np.linspace(x0 + 0.2, x1 - 0.2, 2 * n + 1)
    ys = y + h * np.array([0.0] + [(-1) ** k for k in range(2 * n - 1)] + [0.0])
    ax.plot([x0, xs[0]], [y, y], color=colore, lw=1.3)
    ax.plot(xs, ys, color=colore, lw=1.3)
    ax.plot([xs[-1], x1], [y, y], color=colore, lw=1.3)
    ax.text((x0 + x1) / 2, y + 0.25, testo, ha="center", va="bottom", fontsize=7.5, color=colore)


def _pannello_catena(ax: plt.Axes, c, sensore: str, titolo: str) -> None:
    """Un canale della catena polare: percorso attuatore -> sensore (blu), rami fuori dal percorso
    (arancio, tratteggiati) con gli zeri che generano, motore bloccato per il calcolo degli zeri."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    x = {"e": 0.0, "p": 3.4, "m": 6.8, "r": 9.8}
    lato = {"e": 0.9, "p": 1.5, "m": 0.9, "r": 0.9}
    nome = {"e": "encoder", "p": "telescope", "m": "motor", "r": "rear end"}
    simbolo = {"e": r"$\theta_e$", "p": r"$\theta_p$", "m": r"$\theta_m$", "r": r"$\theta_r$"}
    blu, arancio, verde = SERIE[0], SERIE[1], SERIE[2]
    percorso = ["e", "p", "m"] if sensore == "e" else ["m"]
    rami = [["r"]] if sensore == "e" else [["e", "p"], ["r"]]
    zeri = c.zeri_e if sensore == "e" else c.zeri_m

    def estremi(nodi: Sequence[str], margine: float) -> tuple[float, float]:
        return x[nodi[0]] - lato[nodi[0]] / 2 - margine, x[nodi[-1]] + lato[nodi[-1]] / 2 + margine

    # percorso attuatore -> sensore
    a, b = estremi(percorso, 0.3 if sensore == "e" else 0.12)
    ax.add_patch(FancyBboxPatch((a, -1.15), b - a, 2.3, boxstyle="round,pad=0.02,rounding_size=0.3",
                                facecolor=blu, alpha=0.12, edgecolor=blu, lw=1.0, zorder=1))
    fattore = r"$k_{\theta,1}\,k_{\theta,2}$" if sensore == "e" else "1 (collocated)"
    ax.text((a + b) / 2, -1.25, "path actuator → sensor\nno zeros, only the constant factor " + fattore,
            ha="center", va="top", fontsize=7.5, color=blu)
    # rami fuori dal percorso
    for ramo in rami:
        # il ramo comprende la molla che lo collega al motore bloccato
        a, b = estremi(ramo, 0.2)
        if ramo == ["r"]:
            a = x["m"] + lato["m"] / 2 + 0.2
        else:
            b = x["m"] - lato["m"] / 2 - 0.2
        ax.add_patch(FancyBboxPatch((a, -1.0), b - a, 2.0, boxstyle="round,pad=0.02,rounding_size=0.3",
                                    facecolor="none", edgecolor=arancio, lw=1.4, ls="--", zorder=2))
        testi = [f"{f:.4g} Hz: {d}" for f, d in zeri if ("rear" in d) == (ramo == ["r"])]
        ax.text((a + b) / 2, 1.55, "off-path branch → zero" + ("s" if len(testi) > 1 else "") + "\n"
                + "\n".join(testi), ha="center", va="bottom", fontsize=7.5, color=arancio)
    # molle
    for (u, v), kk, nm in ((("e", "p"), c.k[0], r"$k_{\theta,1}$"), (("p", "m"), c.k[1], r"$k_{\theta,2}$"),
                           (("m", "r"), c.k[2], r"$k_{\theta,3}$")):
        _molla(ax, x[u] + lato[u] / 2, x[v] - lato[v] / 2, 0.0, f"{nm} = {kk:.3g}", INCHIOSTRO_SECONDARIO)
    # inerzie
    for n in ("e", "p", "m", "r"):
        L = lato[n]
        ax.add_patch(Rectangle((x[n] - L / 2, -L / 2), L, L, facecolor="white", edgecolor=INCHIOSTRO, lw=1.3,
                               zorder=3))
        ax.text(x[n], 0.08, simbolo[n], ha="center", va="bottom", fontsize=9, zorder=4)
        ax.text(x[n], -0.02, nome[n], ha="center", va="top", fontsize=6.5, zorder=4)
        ax.text(x[n], -L / 2 + 0.05, f"J = {c.J[n]:.3g}", ha="center", va="bottom", fontsize=6,
                color=INCHIOSTRO_SECONDARIO, zorder=4)
    # attuatore: coppia sul motore
    ax.add_patch(FancyArrowPatch((x["m"] - 0.5, 0.62), (x["m"] + 0.5, 0.62), connectionstyle="arc3,rad=-0.7",
                                 arrowstyle="-|>", mutation_scale=12, color=INCHIOSTRO, lw=1.4, zorder=5))
    ax.text(x["m"], 1.0, r"$\tau_\theta$", ha="center", va="bottom", fontsize=10, zorder=5)
    # motore bloccato nel calcolo degli zeri
    ax.add_patch(Rectangle((x["m"] - 0.45, -0.62), 0.9, 0.1, hatch="////", facecolor="none",
                           edgecolor=arancio, lw=0.8, zorder=4))
    # sensore
    ax.annotate("sensor", xy=(x[sensore] - lato[sensore] / 2, lato[sensore] / 2),
                xytext=(x[sensore] - lato[sensore] / 2 - 0.9, 1.0), fontsize=8, color=verde, zorder=6,
                arrowprops={"arrowstyle": "-|>", "color": verde, "lw": 1.2})
    n_molle = len(percorso) - 1
    ax.text(-1.4, -2.15,
            f"poles (free chain, the same in every channel): {', '.join(f'{f:.4g}' for f in c.poli)} Hz"
            + (f"   ({c.poli[0]:.2g} Hz: telescope on the gravity stiffness k_g = {c.k_g:.2g} N m/rad)"
               if c.poli[0] < 1.0 and c.k_g > 0 else "") + "\n"
            f"zeros = n − 1 − (springs on the path) = 4 − 1 − {n_molle} = {3 - n_molle}:  "
            + ", ".join(f"{f:.4g}" for f, _ in zeri) + " Hz   (hatched: motor clamped to compute the zeros)",
            ha="left", va="top", fontsize=8, color=INCHIOSTRO)
    ax.set_title(titolo, fontsize=10, loc="left")
    ax.set_xlim(-1.5, 11.0)
    ax.set_ylim(-2.9, 2.9)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_gid("disegno")


def catene_poli_zeri(c, *, segnaposto: bool) -> Figure:
    """Catena polare concentrata disegnata per i canali coppia -> encoder e coppia -> motore."""
    fig, assi = plt.subplots(2, 1, figsize=(11.0, 9.5))
    _pannello_catena(assi[0], c, "e", r"(a) $\tau_\theta \to \theta_e$: sensor at the encoder (non-collocated)")
    _pannello_catena(assi[1], c, "m", r"(b) $\tau_\theta \to \theta_m$: sensor at the motor (collocated)")
    _titolo(fig, f"Lumped polar chain at theta0 = {np.rad2deg(c.theta0):g} deg: where poles and zeros come from "
                 "(J in kg m^2, k in N m/rad)", segnaposto)
    fig.tight_layout()
    return fig


def _molla_verticale(ax: plt.Axes, x: float, y0: float, y1: float, testo: str, colore: str) -> None:
    """Molla a zig-zag verticale fra y0 e y1 con l'etichetta a destra."""
    n, h = 5, 0.14
    ys = np.linspace(y0 + 0.15, y1 - 0.15, 2 * n + 1)
    xs = x + h * np.array([0.0] + [(-1) ** k for k in range(2 * n - 1)] + [0.0])
    ax.plot([x, x], [y0, ys[0]], color=colore, lw=1.3)
    ax.plot(xs, ys, color=colore, lw=1.3)
    ax.plot([x, x], [ys[-1], y1], color=colore, lw=1.3)
    ax.text(x + 0.3, (y0 + y1) / 2, testo, ha="left", va="center", fontsize=7.5, color=colore)


def _pannello_ramo(ax: plt.Axes, c, sensore: str, titolo: str) -> None:
    """Catena polare con il telescopio T appeso al mozzo P tramite k_b: percorso attuatore -> sensore (blu),
    rami fuori dal percorso (arancio) con i loro zeri, motore bloccato per il calcolo degli zeri."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    pos = {"e": (0.0, 0.0), "P": (3.4, 0.0), "m": (6.8, 0.0), "r": (9.8, 0.0), "T": (3.4, 2.7)}
    lato = {"e": 0.9, "P": 1.0, "m": 0.9, "r": 0.9, "T": 1.5}
    nome = {"e": "encoder", "P": "hub", "m": "motor", "r": "rear end", "T": "telescope"}
    simbolo = {"e": r"$\theta_e$", "P": r"$\theta_P$", "m": r"$\theta_m$", "r": r"$\theta_r$", "T": r"$\theta_c$"}
    blu, arancio, verde = SERIE[0], SERIE[1], SERIE[2]
    zeri = c.zeri_e if sensore == "e" else c.zeri_m

    def box(x0: float, y0: float, x1: float, y1: float, **kw) -> None:
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0.02,rounding_size=0.3", **kw))

    xe, xP, xm, xr = (pos[n][0] for n in ("e", "P", "m", "r"))
    yT = pos["T"][1]
    if sensore == "e":
        box(xe - 0.75, -1.05, xm + 0.75, 1.05, facecolor=blu, alpha=0.12, edgecolor=blu, lw=1.0, zorder=1)
        ax.text((xe + xm) / 2, -1.15, "path actuator → sensor: no zeros, only the constant factor "
                r"$k_{\theta,1}\,k_{\theta,2}$", ha="center", va="top", fontsize=7.5, color=blu)
        # rami fuori dal percorso: telescopio su k_b (attaccato al mozzo) e estremita' posteriore
        box(xP - 1.0, 0.75, xP + 2.0, yT + 1.0, facecolor="none", edgecolor=arancio, lw=1.4, ls="--", zorder=2)
        box(xm + 0.65, -0.95, xr + 0.65, 0.95, facecolor="none", edgecolor=arancio, lw=1.4, ls="--", zorder=2)
        f_T = [f for f, d in zeri if d.startswith("telescope")]
        f_r = [f for f, d in zeri if d.startswith("rear")]
        ax.text(xP + 2.2, yT + 0.3, "off-path branch → zero\n" + "\n".join(
            f"{f:.4g} Hz: telescope on $k_b$, hub clamped" for f in f_T), ha="left", va="center", fontsize=7.5,
            color=arancio)
        ax.text((xm + xr) / 2 + 0.5, 1.15, "off-path branch → zero\n" + "\n".join(
            f"{f:.4g} Hz: rear end on $k_{{\\theta,3}}$" for f in f_r), ha="center", va="bottom", fontsize=7.5,
            color=arancio)
        # mozzo bloccato per lo zero del telescopio
        ax.add_patch(Rectangle((xP - 0.5, 0.52), 1.0, 0.1, hatch="////", facecolor="none", edgecolor=arancio,
                               lw=0.8, zorder=4))
    else:
        box(xm - 0.6, -1.05, xm + 0.6, 1.05, facecolor=blu, alpha=0.12, edgecolor=blu, lw=1.0, zorder=1)
        ax.text(xm, -1.15, "path actuator → sensor\nfactor 1 (collocated)", ha="center", va="top",
                fontsize=7.5, color=blu)
        box(xe - 0.65, -0.95, xm - 0.65, yT + 1.0, facecolor="none", edgecolor=arancio, lw=1.4, ls="--", zorder=2)
        box(xm + 0.65, -0.95, xr + 0.65, 0.95, facecolor="none", edgecolor=arancio, lw=1.4, ls="--", zorder=2)
        f_l = [f for f, d in zeri if not d.startswith("rear")]
        f_r = [f for f, d in zeri if d.startswith("rear")]
        ax.text(xe - 0.5, yT + 0.3, "off-path piece (encoder, hub, telescope)\n→ zeros "
                + ", ".join(f"{f:.4g}" for f in f_l) + " Hz\n(motor clamped)", ha="left", va="center",
                fontsize=7.5, color=arancio)
        ax.text((xm + xr) / 2 + 0.5, 1.15, "off-path branch → zero\n" + "\n".join(
            f"{f:.4g} Hz: rear end on $k_{{\\theta,3}}$" for f in f_r), ha="center", va="bottom", fontsize=7.5,
            color=arancio)
    # molle
    for (u, v), kk, nm in ((("e", "P"), c.k[0], r"$k_{\theta,1}$"), (("P", "m"), c.k[1], r"$k_{\theta,2}$"),
                           (("m", "r"), c.k[2], r"$k_{\theta,3}$")):
        _molla(ax, pos[u][0] + lato[u] / 2, pos[v][0] - lato[v] / 2, 0.0, f"{nm} = {kk:.3g}",
               INCHIOSTRO_SECONDARIO)
    _molla_verticale(ax, xP, lato["P"] / 2, yT - lato["T"] / 2, f"$k_b$ = {c.k[3]:.3g}", arancio)
    # inerzie
    for n, (xn, yn) in pos.items():
        L = lato[n]
        ax.add_patch(Rectangle((xn - L / 2, yn - L / 2), L, L, facecolor="white", edgecolor=INCHIOSTRO, lw=1.3,
                               zorder=3))
        ax.text(xn, yn + 0.08, simbolo[n], ha="center", va="bottom", fontsize=9, zorder=4)
        ax.text(xn, yn - 0.02, nome[n], ha="center", va="top", fontsize=6.5, zorder=4)
        ax.text(xn, yn - L / 2 + 0.05, f"J = {c.J[n]:.3g}", ha="center", va="bottom", fontsize=6,
                color=INCHIOSTRO_SECONDARIO, zorder=4)
    # attuatore e motore bloccato
    ax.add_patch(FancyArrowPatch((xm - 0.5, 0.62), (xm + 0.5, 0.62), connectionstyle="arc3,rad=-0.7",
                                 arrowstyle="-|>", mutation_scale=12, color=INCHIOSTRO, lw=1.4, zorder=5))
    ax.text(xm, 1.0, r"$\tau_\theta$", ha="center", va="bottom", fontsize=10, zorder=5)
    ax.add_patch(Rectangle((xm - 0.45, -0.62), 0.9, 0.1, hatch="////", facecolor="none", edgecolor=arancio,
                           lw=0.8, zorder=4))
    # sensore
    xs, Ls = pos[sensore][0], lato[sensore]
    if sensore == "e":
        xy, xyt = (xs - Ls / 2, -Ls / 2), (xs - Ls / 2 - 0.9, -1.0)
    else:
        xy, xyt = (xs + Ls / 2, -Ls / 2), (xs + Ls / 2 + 0.5, -1.45)
    ax.annotate("sensor", xy=xy, xytext=xyt, fontsize=8, color=verde, zorder=6,
                arrowprops={"arrowstyle": "-|>", "color": verde, "lw": 1.2})
    n_molle = 2 if sensore == "e" else 0
    ax.text(-1.4, -1.85,
            f"poles (free chain, the same in every channel): {', '.join(f'{f:.4g}' for f in c.poli)} Hz\n"
            f"zeros = n − 1 − (springs on the path) = 5 − 1 − {n_molle} = {4 - n_molle}:  "
            + ", ".join(f"{f:.4g}" for f, _ in zeri) + " Hz   (hatched: node clamped to compute the zeros)",
            ha="left", va="top", fontsize=8, color=INCHIOSTRO)
    ax.set_title(titolo, fontsize=10, loc="left")
    ax.set_xlim(-1.5, 11.0)
    ax.set_ylim(-2.6, yT + 1.3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_gid("disegno")


def catena_ramo_poli_zeri(c, *, segnaposto: bool, stimata: bool = True) -> Figure:
    """Catena polare con il telescopio come ramo laterale (k_b), canali coppia -> encoder e coppia -> motore."""
    fig, assi = plt.subplots(2, 1, figsize=(11.0, 11.0))
    _pannello_ramo(assi[0], c, "e", r"(a) $\tau_\theta \to \theta_e$: sensor at the encoder (non-collocated)")
    _pannello_ramo(assi[1], c, "m", r"(b) $\tau_\theta \to \theta_m$: sensor at the motor (collocated)")
    _titolo(fig, "Polar chain with the telescope hanging from the hub on k_b "
                 + ("(fitted to the measurement), " if stimata else "(hub model), ") +
                 f"theta0 = {np.rad2deg(c.theta0):g} deg (J in kg m^2, k in N m/rad)", segnaposto)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Mozzo dell'asse polare (modello hub)
# ---------------------------------------------------------------------------


def disegno_mozzo(pm, rg, *, segnaposto: bool) -> Figure:
    """Geometria del mozzo usata per le rigidezze: (a) vista dall'alto (asse polare orizzontale, asse
    del telescopio verticale), (b) vista lungo l'asse polare, (c) molle in serie. Quote in mm."""
    from matplotlib.patches import Circle, FancyArrowPatch

    tu, fl, an = pm.tubo, pm.flangia, pm.anello
    Ro, Ri, b = an.d_out / 2 / MM, an.d_in / 2 / MM, an.larghezza / MM
    Lf, Lt = fl.lunghezza / MM, tu.lunghezza / MM
    rf, rto, rti = fl.d_out / 2 / MM, tu.d_out / 2 / MM, tu.d_in / 2 / MM
    colore = {"acciaio": "#5b7fb5", "alluminio": "#cfd9e6"}
    acciaio, alluminio = colore.get(tu.materiale.nome, "#5b7fb5"), colore.get(an.materiale.nome, "#cfd9e6")
    nome_mat = {"acciaio": "steel", "alluminio": "aluminium"}
    fig = plt.figure(figsize=(12.0, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 0.65], width_ratios=[1.6, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    # --- (a) vista dall'alto: J orizzontale, K verticale -----------------------------
    for segno in (-1, 1):
        # parete dell'anello (sezione a J = +-R)
        x0 = segno * Ri if segno > 0 else -Ro
        ax.add_patch(Rectangle((x0, -b / 2), Ro - Ri, b, facecolor=alluminio, edgecolor=INCHIOSTRO, lw=0.8))
        # flangia piena e tubo sede dei cuscinetti
        xf = Ro if segno > 0 else -Ro - Lf
        ax.add_patch(Rectangle((xf, -rf), Lf, 2 * rf, facecolor=acciaio, edgecolor=INCHIOSTRO, lw=0.8, alpha=0.85))
        xt = Ro + Lf if segno > 0 else -Ro - Lf - Lt
        for s2 in (-1, 1):
            y0 = rti if s2 > 0 else -rto
            ax.add_patch(Rectangle((xt, y0), Lt, rto - rti, facecolor=acciaio, edgecolor=INCHIOSTRO, lw=0.8,
                                   alpha=0.85))
        ax.add_patch(Rectangle((xt, -rti), Lt, 2 * rti, facecolor="none", edgecolor=INCHIOSTRO_TENUE,
                               hatch="///", lw=0.4))
    # anello visto di taglio: bordi superiore e inferiore tratteggiati
    for y in (-b / 2, b / 2):
        ax.plot([-Ri, Ri], [y, y], color=INCHIOSTRO_SECONDARIO, lw=0.8, ls="--")
    # telescopio infilato (contorno)
    ax.add_patch(Rectangle((-Ri, -2.3 * b), 2 * Ri, 4.6 * b, facecolor="none", edgecolor=SERIE[2], lw=1.0,
                           ls=":"))
    ax.text(0, 2.05 * b, f"telescope tube (Ø{2 * Ri:g}) slid into the ring", ha="center", fontsize=7,
            color=SERIE[2])
    xmax = Ro + Lf + Lt
    ax.axhline(0, color=INCHIOSTRO_TENUE, lw=0.6, ls="-.")
    ax.text(xmax + 8, 0, r"polar axis $\hat J$", va="center", fontsize=8)
    ax.annotate("", xy=(0, 2.4 * b), xytext=(0, -2.4 * b),
                arrowprops={"arrowstyle": "-", "color": INCHIOSTRO_TENUE, "ls": "-."})
    ax.text(4, -2.45 * b, r"telescope axis $\hat K$ (ring axis)", fontsize=7, va="top")
    # quote
    def quota(x0, x1, y, testo):
        ax.annotate("", xy=(x0, y), xytext=(x1, y), arrowprops={"arrowstyle": "<->", "lw": 0.7})
        ax.text((x0 + x1) / 2, y + 4, testo, ha="center", va="bottom", fontsize=7)
    quota(-Ro, Ro, -b / 2 - 30, f"ring Ø{2 * Ro:g}/{2 * Ri:g}")
    quota(Ro, Ro + Lf, rf + 45, f"flange {Lf:g}")
    quota(Ro + Lf, xmax, rf + 18, f"tube {Lt:g}")
    ax.annotate("", xy=(-Ro - 25, -b / 2), xytext=(-Ro - 25, b / 2), arrowprops={"arrowstyle": "<->", "lw": 0.7})
    ax.text(-Ro - 30, 0, f"width {b:g}", rotation=90, ha="right", va="center", fontsize=7)
    ax.text(Ro + Lf + Lt / 2, -rto - 8, f"Ø{2 * rto:g}/{2 * rti:g}\nbearing seat", ha="center", va="top", fontsize=7)
    ax.text(Ro + Lf / 2, rf + 3, f"Ø{2 * rf:g}", ha="center", va="bottom", fontsize=6)
    ax.plot(0, 0, marker="+", color=INCHIOSTRO, ms=10)
    ax.text(5, 6, "H ($o_2$)", fontsize=8)
    ax.set_xlim(-xmax - 45, xmax + 70)
    ax.set_ylim(-2.6 * b, 2.6 * b)
    ax.set_aspect("equal")
    ax.set_title(f"(a) view from above: side hubs ({nome_mat.get(tu.materiale.nome, tu.materiale.nome)}) and "
                 f"support ring ({nome_mat.get(an.materiale.nome, an.materiale.nome)})", fontsize=9, loc="left")
    ax.axis("off")
    ax.set_gid("disegno")

    # --- (b) vista lungo l'asse polare -------------------------------------------------
    ax2.add_patch(Circle((0, 0), Ro, facecolor=alluminio, edgecolor=INCHIOSTRO, lw=0.8))
    ax2.add_patch(Circle((0, 0), Ri, facecolor="white", edgecolor=INCHIOSTRO, lw=0.8))
    ax2.add_patch(Circle((0, 0), Ri - 2, facecolor="none", edgecolor=SERIE[2], lw=1.0, ls=":"))
    ax2.text(0, Ri * 0.55, "telescope", ha="center", fontsize=7, color=SERIE[2])
    ax2.text(0, -Ro - 12, f"ring seen along $\\hat K$: Ø{2 * Ro:g}/{2 * Ri:g},\nwall {Ro - Ri:g} mm",
             ha="center", va="top", fontsize=7)
    for segno in (-1, 1):
        ax2.add_patch(Rectangle((segno * Ro if segno > 0 else -Ro - Lf - Lt, -rf), Lf + Lt, 2 * rf,
                                facecolor=acciaio, edgecolor=INCHIOSTRO, lw=0.8, alpha=0.6))
    ax2.annotate("", xy=(-Ro - Lf - Lt, 0), xytext=(Ro + Lf + Lt, 0),
                 arrowprops={"arrowstyle": "-", "color": INCHIOSTRO_TENUE, "ls": "-."})
    ax2.text(0, -6, r"$\hat J$", ha="center", va="top", fontsize=8)
    ax2.add_patch(FancyArrowPatch((-Ro - Lf - Lt / 2, -rf - 18), (-Ro - Lf - Lt / 2, rf + 18),
                                  connectionstyle="arc3,rad=0.9", arrowstyle="-|>", mutation_scale=10,
                                  color=SERIE[1]))
    ax2.text(0, Ro + 10, "torque about $\\hat J$ at the two attachment points:\n"
             "each half ring twists (wall torsion) and bends out of plane",
             ha="center", va="bottom", fontsize=7)
    ax2.set_xlim(-Ro - Lf - Lt - 15, Ro + Lf + Lt + 15)
    ax2.set_ylim(-Ro - 60, Ro + 60)
    ax2.set_aspect("equal")
    ax2.set_title(r"(b) view along the telescope axis $\hat K$", fontsize=9, loc="left")
    ax2.axis("off")
    ax2.set_gid("disegno")

    # --- (c) molle in serie ------------------------------------------------------------
    pezzi = [("tube", rg.k_tubo), ("flange", rg.k_flangia), ("1/2 ring", 2 * rg.k_anello),
             ("1/2 ring", 2 * rg.k_anello), ("flange", rg.k_flangia), ("tube", rg.k_tubo)]
    xs = np.linspace(0.0, 12.0, len(pezzi) + 1)
    nodi = ["bearing\n(node 2)", "", "", "H\n(hub centre)", "", "", "bearing\n(node 5)"]
    for x, n in zip(xs, nodi):
        ax3.plot(x, 0, "o", color=INCHIOSTRO, ms=5 if n else 3)
        if n:
            ax3.text(x, -0.55, n, ha="center", va="top", fontsize=7)
    for (nome, k), x0, x1 in zip(pezzi, xs[:-1], xs[1:]):
        _molla(ax3, x0, x1, 0.0, f"{nome}\n{k:.3g}", INCHIOSTRO_SECONDARIO)
    ax3.text(6.0, 1.15, f"half hub k = {rg.k_lato:.4g} N m/rad   |   whole hub (bearing to bearing) "
             f"k = {rg.k_totale:.4g} N m/rad   |   telescope on H through k_b = {pm.k_b:.3g} N m/rad",
             ha="center", fontsize=8)
    _molla_verticale(ax3, xs[3], 0.15, 0.85, "", SERIE[1])
    ax3.set_xlim(-0.8, 12.8)
    ax3.set_ylim(-1.4, 1.5)
    ax3.axis("off")
    ax3.set_gid("disegno")
    ax3.set_title("(c) springs in series along the polar axis (N m/rad); the hub inertia sits at H",
                  fontsize=9, loc="left")
    _titolo(fig, "Hub of the polar axis: geometry used for the stiffness (mm)", segnaposto)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Disegno della geometria (controllo rapido, PNG)
# ---------------------------------------------------------------------------


def disegno_geometria(par: Parametri, asse: str) -> Figure:
    """Vista laterale stilizzata in scala reale: segmenti, strati, nodi."""
    a = par.asse(asse)  # type: ignore[arg-type]
    fig, ax = plt.subplots(figsize=(12.0, 3.2))
    x = 0.0
    for el in a.elementi:
        if el.modo == "telescope":
            h = 0.6 * max(_d_max(a), 20)
            ax.add_patch(Rectangle((x, -h), el.lunghezza / MM, 2 * h, fill=True, alpha=0.15,
                                   color="purple"))
            ax.text(x + el.lunghezza / MM / 2, 0, "telescope\n(no shaft)", ha="center",
                    va="center", fontsize=6)
        if el.modo == "manual":
            ax.add_patch(Rectangle((x, -6), el.lunghezza / MM, 12, fill=False, hatch="//",
                                   edgecolor=INCHIOSTRO_SECONDARIO))
            ax.text(x + el.lunghezza / MM / 2, 8, "manual", ha="center", fontsize=6)
        xs = x
        for seg in el.segmenti:
            L = seg.lunghezza / MM
            for st in sorted(seg.strati, key=lambda s: -s.d_out):
                colore = SERIE[0] if st.portante else SERIE[3]
                R, r = st.d_out / MM / 2, st.d_in / MM / 2
                for segno in (1, -1):
                    ax.add_patch(Rectangle((xs, segno * r if segno > 0 else -R), L, R - r,
                                           facecolor=colore, edgecolor=INCHIOSTRO, lw=0.4,
                                           alpha=0.85 if st.materiale.nome == "acciaio" else 0.5))
            xs += L
        ax.axvline(x, color=INCHIOSTRO_TENUE, lw=0.5, ls=":")
        ax.text(x, -0.62 * max(_d_max(a), 20), str(el.id), ha="center", fontsize=7)
        ax.text(x + el.lunghezza / MM / 2, 0.55 * max(_d_max(a), 20), f"e{el.id}", ha="center",
                fontsize=7, color=INCHIOSTRO_SECONDARIO)
        x += el.lunghezza / MM
    ax.axvline(x, color=INCHIOSTRO_TENUE, lw=0.5, ls=":")
    ax.text(x, -0.62 * max(_d_max(a), 20), str(a.n_nodi), ha="center", fontsize=7)
    ax.set_xlim(-5, x + 5)
    lim = 0.7 * max(_d_max(a), 20)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("axial position [mm]")
    ax.set_ylabel("radius [mm]")
    ax.grid(False)
    ax.set_gid("disegno")
    _titolo(fig, f"Geometry as entered, {NOME_ASSE[asse]} axis (blue: torque-carrying, yellow: not torque-carrying, "
                 "transparent: aluminium)", a.segnaposto)
    fig.tight_layout()
    return fig


def _d_max(a) -> float:
    return max((st.d_out / MM for el in a.elementi for s in el.segmenti for st in s.strati),
               default=20.0)
