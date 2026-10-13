"""FRF, poli, zeri invarianti, FdT simboliche ed esportazione (sezione 10.3).

* **FRF**: da :math:`\\mathbf P(j\\omega)^{-1}`, non dai polinomi.
* **Poli**: autovalori di :math:`\\mathbf A`, classificati riusando
  :func:`gimbal.analysis.classifica_radice` con l'interpretazione aggiornata (qui
  i poli complessi sono anche modi strutturali).
* **Zeri** di ogni FdT SISO: zeri invarianti dal pencil di Rosenbrock
  :math:`\\begin{bmatrix}\\mathbf A - s\\mathbb 1 & \\mathbf b_j\\\\ \\mathbf c_i & 0\\end{bmatrix}`,
  autovalori generalizzati finiti.

La FdT coppia -> encoder e' **non collocata**: l'alternanza poli-zeri non e'
garantita. La FdT coppia -> nodo motore e' collocata e la deve mostrare.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Sequence

import numpy as np
import scipy.io as sio
import scipy.linalg as sla
import sympy as sp
from gimbal.analysis import Radice, classifica_radice

# ---------------------------------------------------------------------------
# FRF
# ---------------------------------------------------------------------------


def frf_secondo_ordine(M: np.ndarray, C: np.ndarray, K: np.ndarray, B: np.ndarray,
                       Cout: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """:math:`\\mathbf C_{out}(-\\omega^2\\mathbf M + j\\omega\\mathbf C + \\mathbf K)^{-1}\\mathbf B`.

    Returns
    -------
    np.ndarray
        Forma ``(len(omega), n_uscite, n_ingressi)``.
    """
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    jw = 1j * omega
    P = M[None] * (jw**2)[:, None, None] + C[None] * jw[:, None, None] + K[None]
    X = np.linalg.solve(P, np.broadcast_to(B.astype(complex), (omega.size, *B.shape)))
    return np.einsum("ij,wjk->wik", Cout, X)


# ---------------------------------------------------------------------------
# Poli e zeri
# ---------------------------------------------------------------------------


def classifica(radici: Sequence[complex], *, tol: float = 1e-6) -> list[Radice]:
    """Classifica poli o zeri; per le coppie complesse l'origine e' strutturale."""
    out = []
    for r in radici:
        c = classifica_radice(r, tol=tol)
        if c.tipo == "complesso":
            c = replace(c, origine="coppia complessa: modo torsionale (risonanza o antirisonanza)")
        out.append(c)
    return sorted(out, key=lambda r: (abs(r.valore), float(np.imag(r.valore))))


def zeri_invarianti(A: np.ndarray, b: np.ndarray, c: np.ndarray, d: float = 0.0) -> np.ndarray:
    """Zeri invarianti di una FdT SISO dal pencil di Rosenbrock.

    Autovalori generalizzati di
    :math:`\\left(\\begin{bmatrix}\\mathbf A & \\mathbf b\\\\ \\mathbf c & d\\end{bmatrix},
    \\begin{bmatrix}\\mathbb 1 & 0\\\\ 0 & 0\\end{bmatrix}\\right)`; si scartano
    quelli infiniti.
    """
    n = A.shape[0]
    Ap = np.zeros((n + 1, n + 1))
    Ap[:n, :n] = A
    Ap[:n, n] = b.reshape(-1)
    Ap[n, :n] = c.reshape(-1)
    Ap[n, n] = d
    Bp = np.zeros((n + 1, n + 1))
    Bp[:n, :n] = np.eye(n)
    alpha, beta = sla.eig(Ap, Bp, right=False, homogeneous_eigvals=True)
    finiti = np.abs(beta) > 1e-10 * np.maximum(np.abs(alpha), 1.0)
    z = alpha[finiti] / beta[finiti]
    scala = max(float(np.max(np.abs(np.linalg.eigvals(A)))), 1.0)
    return z[np.abs(z) < 1e6 * scala]


def zeri_matrice(A: np.ndarray, B: np.ndarray, C: np.ndarray, uscite: Sequence[str]) -> dict[tuple[str, int], np.ndarray]:
    """Zeri invarianti di tutte le FdT SISO (uscita, ingresso)."""
    return {(u, j): zeri_invarianti(A, B[:, j], C[i])
            for i, u in enumerate(uscite) for j in range(B.shape[1])}


def interlacciati(poli: np.ndarray, zeri: np.ndarray, *, tol: float = 1e-6) -> bool:
    """Poli e zeri immaginari puri si alternano sull'asse :math:`j\\omega`."""
    wp = np.sort(np.abs(np.imag(poli[np.imag(poli) >= -tol])))
    wz = np.sort(np.abs(np.imag(zeri[np.imag(zeri) >= -tol])))
    unione = sorted([(w, "p") for w in wp] + [(w, "z") for w in wz])
    tipi = [t for _, t in unione]
    return all(a != b for a, b in zip(tipi[:-1], tipi[1:]))


# ---------------------------------------------------------------------------
# FdT simboliche delle catene isolate (Omega_0 = 0)
# ---------------------------------------------------------------------------


def fdt_simboliche_catena(n: int, i_motore: int, i_encoder: int, nome: str,
                          i_carico: int | None = None, ruoli: Sequence[str] | None = None,
                          con_k_g: bool = False) -> dict[str, sp.Expr]:
    """FdT di una catena pura (masse diagonali, molle in serie, senza smorzamento).

    :math:`G_{ij}(s) = [\\operatorname{adj}\\mathbf P(s)]_{ij}/\\det\\mathbf P(s)`,
    con :math:`\\mathbf P(s) = \\operatorname{diag}(J_i)s^2 + \\mathbf K`. Con
    ``con_k_g`` si aggiunge :math:`k_g` al nodo di carico (polare).
    """
    s = sp.Symbol("s")
    ruoli = list(ruoli) if ruoli is not None else [str(i + 1) for i in range(n)]
    J = [sp.Symbol(f"J_{{\\{nome},{r}}}", positive=True) for r in ruoli]
    k = [sp.Symbol(f"k_{{\\{nome},{i + 1}}}", positive=True) for i in range(n - 1)]
    K = sp.zeros(n, n)
    for i, ki in enumerate(k):
        K[i, i] += ki
        K[i + 1, i + 1] += ki
        K[i, i + 1] -= ki
        K[i + 1, i] -= ki
    if con_k_g and i_carico is not None:
        K[i_carico, i_carico] += sp.Symbol("k_g", real=True)
    P = sp.diag(*J) * s**2 + K
    det = sp.factor(P.det())
    adj = P.adjugate()
    out = {"P": P, "det": det,
           "N_encoder": sp.factor(adj[i_encoder, i_motore]),
           "N_motore": sp.factor(adj[i_motore, i_motore])}
    if i_carico is not None:
        out["N_carico"] = sp.factor(adj[i_carico, i_motore])
    for chiave in [k for k in out if k.startswith("N_")]:
        out[chiave[2:]] = out[chiave] / det
    return out


# ---------------------------------------------------------------------------
# Esportazione
# ---------------------------------------------------------------------------


def esporta(ml, cartella: str | Path, nome: str) -> list[Path]:
    """Matrici di stato in ``.npz`` e ``.mat``, zeri/poli/guadagno e coefficienti in JSON.

    I coefficienti num/den sono ricavati da ``scipy.signal.ss2tf``: servono solo
    come scambio verso altri strumenti; le FRF di riferimento restano quelle da
    :math:`\\mathbf P(j\\omega)^{-1}`.
    """
    from scipy.signal import ss2tf

    from .linearize import USCITE

    cartella = Path(cartella)
    cartella.mkdir(parents=True, exist_ok=True)
    dati = dict(A=ml.A, B=ml.Bu, C=ml.Cy, D=np.zeros((ml.Cy.shape[0], 2)), M0=ml.M0, C0=ml.C0,
                K0=ml.K0, E=ml.E, theta0=ml.theta0, Omega0=ml.Omega0, q0=ml.punto.q0,
                tau0=ml.punto.tau0, uscite=np.array(USCITE), stati=np.array(
                    list(ml.accoppiato.struttura.nomi) + [f"d{n}" for n in ml.accoppiato.struttura.nomi]))
    p_npz = cartella / f"{nome}.npz"
    np.savez(p_npz, **dati)
    p_mat = cartella / f"{nome}.mat"
    sio.savemat(p_mat, {k: (v if not isinstance(v, np.ndarray) or v.dtype.kind != "U" else list(v))
                        for k, v in dati.items()})
    poli = np.linalg.eigvals(ml.A)
    fdt = {}
    for i, u in enumerate(USCITE):
        for j, ingr in enumerate(("tau_psi", "tau_theta")):
            num, den = ss2tf(ml.A, ml.Bu[:, [j]], ml.Cy[[i]], np.zeros((1, 1)))
            z = zeri_invarianti(ml.A, ml.Bu[:, j], ml.Cy[i])
            fdt[f"{u}/{ingr}"] = {
                "num": np.real_if_close(num[0]).tolist(), "den": np.real(den).tolist(),
                "zeri": [[float(np.real(v)), float(np.imag(v))] for v in z],
            }
    p_json = cartella / f"{nome}.json"
    p_json.write_text(json.dumps({
        "theta0_deg": float(np.rad2deg(ml.theta0)), "Omega0_deg_s": float(np.rad2deg(ml.Omega0)),
        "tau0_Nm": ml.punto.tau0.tolist(),
        "poli": [[float(np.real(v)), float(np.imag(v))] for v in poli],
        "fdt": fdt,
        "avviso": "risultato da parametri segnaposto" if ml.accoppiato.rid_psi.completo.segnaposto else "",
    }, indent=1), encoding="utf-8")
    return [p_npz, p_mat, p_json]
