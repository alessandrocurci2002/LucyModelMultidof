"""Figure del modello multi-GdL (sezione 11).

Lo stile e' quello del modello a 2 GdL (:mod:`gimbal.plotting`, riusato): stessa
palette categorica in ordine fisso, colore mai unica codifica (tratto diverso
per ogni serie), PDF vettoriale e PNG a 300 dpi. Etichette in italiano, Hz e
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
                             TRATTI, applica_stile, salva)
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from . import fem
from . import linearize as lz
from . import tf
from .analysis import SweepTheta
from .coupled import Sistema
from .params import MM, NOTA_REALISMO, Parametri
from .simulate import RisultatoScenario

__all__ = ["applica_stile", "salva"]


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
    ax_m.set_ylabel("modulo [dB rad/(N m)]")
    ax_f.set_ylabel("fase [deg]")
    ax_f.set_xlabel("frequenza [Hz]")
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
    for k, (nome, ml) in enumerate(modelli):
        G = ml.frf(w, list(uscite))
        for i in range(2):
            for j in range(2):
                ax_m, ax_f = assi[2 * i, j], assi[2 * i + 1, j]
                ax_m.plot(f, _db(G[:, i, j]), label=nome, **_stile(k))
                ax_f.plot(f, _fase(G[:, i, j]), **_stile(k))
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
                          label="zeri" if (i, j) == (0, 0) else None)
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
            (f"completo, {mm}", (M, K, C, mod.i_motore, mod.i_encoder)),
            (f"ridotto ({s.forma_massa}), {mm}", (Mr, Kr, Cr, rid.i_motore, rid.i_encoder)),
        ):
            for col, uscita in enumerate((ie, im)):
                G = frf_encoder_isolato(MM_, CC, KK, im, uscita, w)
                assi[0, col].plot(f, _db(G), label=nome, **_stile(k))
                assi[1, col].plot(f, _fase(G), **_stile(k))
            k += 1
    simbolo = r"\psi" if asse == "azimut" else r"\theta"
    _assi_bode(assi[0, 0], assi[1, 0], f"coppia -> encoder ${simbolo}_e$ (non collocata)")
    _assi_bode(assi[0, 1], assi[1, 1], f"coppia -> motore ${simbolo}_m$ (collocata)")
    assi[0, 0].legend(loc="lower left", fontsize=7)
    _titolo(fig, f"Asse {asse} isolato, theta0 = {np.rad2deg(theta0):g} deg: completo vs ridotto",
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
        ax.plot(x, phi, marker="o", ms=3, label=f"completo, f = {f_k:.4g} Hz", **_stile(0))
        if k < len(w2r):
            phir = rid.T @ Phir[:, k]
            phir = phir / np.max(np.abs(phir)) * (np.sign(phir @ phi) or 1.0)
            f_r = np.sign(w2r[k]) * np.sqrt(abs(w2r[k])) / (2 * np.pi)
            ax.plot(x, phir, label=f"ridotto riespanso, f = {f_r:.4g} Hz", **_stile(1))
        ax.plot(x[list(mod.masters)], phi[list(mod.masters)], "s", ms=7, mfc="none",
                color=INCHIOSTRO, label="nodi master" if k == indici[0] else None)
        ax.set_ylabel(f"modo {k + 1}")
        ax.legend(loc="best", fontsize=7)
        ax.axhline(0.0, color=GRIGLIA, lw=0.8)
    assi[-1].set_xticks(x)
    assi[-1].set_xticklabels(mod.etichette, fontsize=7)
    assi[-1].set_xlabel("nodo (posizione lungo l'albero in scala)")
    _titolo(fig, f"Forme modali, asse {asse} isolato, theta0 = {np.rad2deg(theta0):g} deg "
                 f"({s.modello_massa}, {s.forma_massa})", segnaposto)
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
    ax.set_xlabel(r"$\theta_0$ = elevazione [deg] ($\theta > 0$ alza il puntamento)")
    ax.set_ylabel("frequenza dei modi complessi [Hz]")
    ax.legend()
    _titolo(fig, "Frequenze proprie del sistema accoppiato in funzione di theta0", segnaposto)
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
        ax.set_xlabel("frequenza [Hz]")
        ax.set_ylabel("modulo [dB rad/(N m)]")
    assi[0].legend()
    _titolo(fig, f"FRF incrociate coppia -> encoder, theta0 = {np.rad2deg(theta0):g} deg",
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
        _assi_bode(assi[0, col], assi[1, col], f"risposta a {ingresso}")
        assi[0, col].legend(loc="lower left", fontsize=7)
    _titolo(fig, f"Errore di misura encoder - angolo vero, theta0 = {np.rad2deg(ml.theta0):g} deg, "
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
            ax.plot(r.t * 1e3, np.rad2deg(r.y_nl[:, i]) * 3600, label=f"{lz.ETICHETTE_USCITE[u]} non lineare",
                    color=SERIE[k], linestyle="-")
            if r.y_lin is not None:
                ax.plot(r.t * 1e3, np.rad2deg(r.y_lin[:, i]) * 3600, label="lineare",
                        color=SERIE[k], linestyle="--", lw=1.0)
        ax.set_ylabel("scarto dall'equilibrio [arcsec]")
        ax.legend(fontsize=7, loc="best")
    assi[-1].set_xlabel("t [ms]")
    _titolo(fig, f"{r.nome}: {r.descrizione} (theta0 = {np.rad2deg(r.theta0):g} deg, "
                 f"Omega0 = {np.rad2deg(r.Omega0):g} deg/s)", segnaposto)
    fig.tight_layout()
    return fig


def conservazione(r: RisultatoScenario, *, segnaposto: bool) -> Figure:
    assert r.energia is not None and r.momento is not None
    fig, assi = plt.subplots(2, 1, figsize=(8.0, 5.0), sharex=True)
    assi[0].plot(r.t * 1e3, (r.energia - r.energia[0]) / abs(r.energia[0]), color=SERIE[0])
    assi[0].set_ylabel(r"$(E - E_0)/|E_0|$")
    assi[1].plot(r.t * 1e3, (r.momento - r.momento[0]) / abs(r.momento[0]), color=SERIE[1], ls="--")
    assi[1].set_ylabel(r"$(p_\psi - p_{\psi,0})/|p_{\psi,0}|$")
    assi[1].set_xlabel("t [ms]")
    _titolo(fig, f"S3: conservazione con C = 0, tau = 0 (derive {r.derive['energia']:.1e}, "
                 f"{r.derive['momento_psi']:.1e})", segnaposto)
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
            ax.text(x + el.lunghezza / MM / 2, 0, "telescopio\n(nessun albero)", ha="center",
                    va="center", fontsize=6)
        if el.modo == "manual":
            ax.add_patch(Rectangle((x, -6), el.lunghezza / MM, 12, fill=False, hatch="//",
                                   edgecolor=INCHIOSTRO_SECONDARIO))
            ax.text(x + el.lunghezza / MM / 2, 8, "manuale", ha="center", fontsize=6)
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
    ax.set_xlabel("posizione assiale [mm]")
    ax.set_ylabel("raggio [mm]")
    ax.grid(False)
    _titolo(fig, f"Geometria inserita, asse {asse} (blu: portante, giallo: non portante, "
                 "trasparente: alluminio)", a.segnaposto)
    fig.tight_layout()
    return fig


def _d_max(a) -> float:
    return max((st.d_out / MM for el in a.elementi for s in el.segmenti for st in s.strati),
               default=20.0)
