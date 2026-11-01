"""Guyan contro troncamento modale (Gawronski) sul FEM torsionale di ciascun asse.

Per ogni asse isolato (azimut, polare) confronta, sullo STESSO modello completo
(matrici M, K del FEM con il payload sul nodo di aggancio):

* **Guyan**           : condensazione statica sui master (reduction.riduci).
* **Guyan lumped**    : stessa T, massa diagonalizzata per somma di riga.
* **Gawronski**       : troncamento modale del FEM completo, modo rigido + k modi
                        flessibili (cap. 2 e 4 del libro), con e senza correzione
                        statica dei modi scartati.
* **SEREP**           : stessi modi del troncamento, scritti nei GdL dei master
                        (coordinate fisiche).
* **Craig-Bampton**   : Guyan piu' n modi a interfaccia fissa (reduction.riduci_craig_bampton).

Uscite: tabella delle frequenze, MAC, errore statico, norme di Hankel dei modi,
banda di validita' delle FRF coppia->encoder e coppia->motore, figura Bode.

Uso (dalla radice del repository, con l'ambiente del progetto attivo)::

    python confronto_guyan_modale.py --config config.yaml --theta0 0 --out results/guyan_vs_modale

Il payload e' preso da ``coupled.costruisci_sistema`` (pipeline completa). Se il
pacchetto del modello a 2 GdL non e' installato si puo' usare ``--payload-manuale``,
che stima I_yy,p, J_psi,p e J^ax direttamente da config.yaml (solo per
``frame: cad_cg``, asse polare +X nella posa CAD): e' una stima, non il percorso
ufficiale, e il rapporto lo segnala.

Tutti i dati di config.yaml con ``placeholder: true`` restano segnaposto: i numeri
qui sotto sono risultati da parametri segnaposto.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as sla

from gimbal_mdof import fem
from gimbal_mdof import params as pr
from gimbal_mdof import reduction as rd

SOGLIA_DB = 1.0
SOGLIA_DEG = 5.0
SOGLIA_DB_REPO = 3.0  # come analysis.SOGLIA_DB del repository (solo modulo)


# ---------------------------------------------------------------------------
# Matrici dell'asse isolato
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Asse:
    nome: str
    M: np.ndarray
    K: np.ndarray
    masters: list[int]
    i_enc: int
    i_mot: int
    i_agg: int  # nodo di aggancio del payload (P per il polare, nodo 9 per l'azimut)
    etichette: tuple[str, ...]
    beta: float
    nota: str


def _payload_manuale(par: pr.Parametri, mod_polare: fem.ModelloAsse, theta0: float):
    """Stima (J_add_azimut, I_yy_p) da config.yaml, frame cad_cg, polare +X, puntamento +Y."""
    p = par.payload
    if p.terna != "cad_cg":
        raise SystemExit("--payload-manuale supporta solo frame: cad_cg")
    I, m, rG = np.asarray(p.I), float(p.m), np.asarray(p.r_G_N)
    d = rG - np.array([0.0, 0.0, par.terne.z_o2])
    I_yy_p = I[0, 0] + m * (d[1] ** 2 + d[2] ** 2)
    # momento del payload attorno a z (azimut) in funzione dell'elevazione: rotazione attorno a x
    c, s = np.cos(-theta0), np.sin(-theta0)  # theta > 0 alza il puntamento (convenzione config)
    R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    Ir = R @ I @ R.T
    rr = R @ rG
    J_psi_p = Ir[2, 2] + m * (rr[0] ** 2 + rr[1] ** 2)
    J_ax = float(sum(pz.pezzo.I_t for pz in mod_polare.pezzi))  # Q7 = None: nessuno Steiner
    return J_ax + J_psi_p, I_yy_p


def carica_assi(config: Path, theta0_deg: float, modello_massa: str, payload_manuale: bool):
    par = pr.parametri_da_config(config)
    th = np.deg2rad(theta0_deg)
    out: list[Asse] = []
    if not payload_manuale:
        from gimbal_mdof import coupled as cp

        s = cp.costruisci_sistema(par, modello_massa, "guyan")
        for nome, mod, beta in (("azimut", s.asse_psi, s.beta_psi), ("polare", s.asse_theta, s.beta_theta)):
            M, K, _ = s.isolato_completo(nome, th)
            out.append(Asse(nome, M, K, list(mod.masters), mod.i_encoder, mod.i_motore,
                            mod.i_P if nome == "polare" else mod.i_estremo, mod.etichette, beta,
                            "payload da coupled.costruisci_sistema"))
        return par, out
    mod_psi = fem.costruisci_asse(par.azimut, par.fem.forma_consistente)
    mod_th = fem.costruisci_asse(par.polare, par.fem.forma_consistente)
    J_az, I_yy = _payload_manuale(par, mod_th, th)
    for nome, mod, J_add in (("azimut", mod_psi, J_az), ("polare", mod_th, I_yy)):
        M, K = fem.matrici_isolate(mod, modello_massa, J_aggiunta=J_add, k_aggiunta=0.0)  # type: ignore[arg-type]
        beta = fem.beta_smorzamento(M, mod.K, par.zeta)
        out.append(Asse(nome, M, K, list(mod.masters), mod.i_encoder, mod.i_motore,
                        mod.i_P if nome == "polare" else mod.i_estremo, mod.etichette, beta,
                        f"payload STIMATO a mano da config.yaml (J_add = {J_add:.4f} kg m^2, k_g = 0)"))
    return par, out


# ---------------------------------------------------------------------------
# Modi del modello completo e riduzioni
# ---------------------------------------------------------------------------


def modi_completi(M: np.ndarray, K: np.ndarray):
    w2, Phi = sla.eigh(K, M)  # Phi^T M Phi = 1
    w2[np.abs(w2) < 1e-9 * w2.max()] = 0.0
    return w2, Phi


def frf_fisica(M, C, K, i_in, i_out, w):
    n = M.shape[0]
    b = np.zeros(n)
    b[i_in] = 1.0
    G = np.empty(len(w), dtype=complex)
    for k, wk in enumerate(w):
        G[k] = np.linalg.solve(K - wk**2 * M + 1j * wk * C, b)[i_out]
    return G


def frf_modale(w2, Phi, beta, i_in, i_out, w, modi, statica=False):
    """Somma dei modi ``modi`` (Gawronski 2.20); ``statica`` aggiunge la flessibilita' residua."""
    om = np.sqrt(w2)
    z = 0.5 * beta * om
    G = np.zeros(len(w), dtype=complex)
    for i in modi:
        G += Phi[i_out, i] * Phi[i_in, i] / (w2[i] - w**2 + 2j * z[i] * om[i] * w)
    if statica:
        scartati = [i for i in range(len(w2)) if i not in set(modi) and w2[i] > 0]
        G += sum(Phi[i_out, i] * Phi[i_in, i] / w2[i] for i in scartati)
    return G


def serep(M, K, Phi, masters):
    n_m = len(masters)
    Phi_m = Phi[np.ix_(masters, range(n_m))]
    T = Phi[:, :n_m] @ np.linalg.inv(Phi_m)
    return T, T.T @ M @ T, T.T @ K @ T


def mac_massa(Phi_a, Phi_b, M):
    """MAC pesato sulla massa fra le colonne di a e di b (diagonale)."""
    out = []
    for a, b in zip(Phi_a.T, Phi_b.T):
        out.append((a @ M @ b) ** 2 / ((a @ M @ a) * (b @ M @ b)))
    return np.array(out)


def f_limite(w, G, Gref):
    r = G / Gref
    fuori = (np.abs(20 * np.log10(np.abs(r))) > SOGLIA_DB) | (np.abs(np.degrees(np.angle(r))) > SOGLIA_DEG)
    k = np.nonzero(fuori)[0]
    return float(w[k[0]] / (2 * np.pi)) if k.size else float("inf")


def f_limite_db(w, G, Gref, soglia=SOGLIA_DB_REPO):
    """Primo scarto di solo modulo oltre ``soglia`` dB (criterio di analysis.py del repository)."""
    d = np.abs(20 * np.log10(np.abs(G / Gref)))
    k = np.nonzero(d > soglia)[0]
    return float(w[k[0]] / (2 * np.pi)) if k.size else float("inf")


def zeri_da_modi(x, r, D=0.0, n_scan=4001):
    """Zeri non smorzati (Hz) di g(x) = sum_i r_i/(x_i - x) + D, con x = omega^2.

    ``x[0] = 0`` e' il modo rigido. Si cerca un cambio di segno in ogni intervallo
    fra poli consecutivi (e oltre l'ultimo) e si raffina con brentq. ``D`` e' il
    termine statico della correzione dei modi scartati.
    """
    from scipy.optimize import brentq

    ordine = np.argsort(x)
    x, r = np.asarray(x)[ordine], np.asarray(r)[ordine]

    def g(xx):
        return float(np.sum(r / (x - xx)) + D)

    bordi = list(x) + [10.0 * x[-1]]
    t = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0.0, 1.0, n_scan)[1:-1])  # fitto agli estremi
    out = []
    for a, b in zip(bordi[:-1], bordi[1:]):
        pts = a + (b - a) * t
        val = np.array([g(p) for p in pts])
        for k in np.nonzero(np.sign(val[:-1]) * np.sign(val[1:]) < 0)[0]:
            out.append(np.sqrt(brentq(g, pts[k], pts[k + 1], xtol=1e-14 * max(pts[k], 1e-30), rtol=1e-13)) / (2 * np.pi))
    return np.array(sorted(out))


def flessibilita(w2, Phi, i_out, i_in, modi):
    """Flessibilita' elastica F_{out,in} = sum_{i>=1} phi_out phi_in / w_i^2 (modo rigido escluso).

    E' il termine statico che la FRF coppia->uscita ha oltre al modo rigido: il
    troncamento modale lo sottostima perche' perde il contributo dei modi scartati.
    """
    return sum(Phi[i_out, i] * Phi[i_in, i] / w2[i] for i in modi if w2[i] > 0)


# ---------------------------------------------------------------------------
# Confronto per asse
# ---------------------------------------------------------------------------


def analizza(ax: Asse, w: np.ndarray, n_cb: tuple[int, ...] = (1, 2, 4)) -> dict:
    M, K, mast = ax.M, ax.K, ax.masters
    n, n_m = M.shape[0], len(mast)
    C = ax.beta * K
    w2, Phi = modi_completi(M, K)
    f_full = np.sqrt(np.abs(w2)) / (2 * np.pi)
    r: dict = {"asse": ax, "f_full": f_full, "w2": w2, "Phi": Phi}

    # --- Guyan (tridiagonale e lumped)
    T = rd.trasformazione_guyan(K, mast)
    Mg, Kg = T.T @ M @ T, T.T @ K @ T
    Mg, Kg = 0.5 * (Mg + Mg.T), 0.5 * (Kg + Kg.T)
    Ml = rd.massa_lumped(Mg)
    w2g, Pg = sla.eigh(Kg, Mg)
    w2l, _ = sla.eigh(Kg, Ml)
    w2g[np.abs(w2g) < 1e-9 * w2g.max()] = 0.0
    w2l[np.abs(w2l) < 1e-9 * w2l.max()] = 0.0
    r["f_guyan"] = np.sqrt(w2g) / (2 * np.pi)
    r["f_lumped"] = np.sqrt(w2l) / (2 * np.pi)
    r["mac_guyan"] = mac_massa(T @ Pg, Phi[:, :n_m], M)

    # --- Craig-Bampton
    r["f_cb"] = {}
    for k in n_cb:
        _, Mc, Kc = rd.riduci_craig_bampton(M, K, mast, k)
        wc, _ = sla.eigh(Kc, Mc)
        wc[np.abs(wc) < 1e-9 * wc.max()] = 0.0
        r["f_cb"][k] = np.sqrt(wc) / (2 * np.pi)

    # --- SEREP (= troncamento modale scritto nei GdL dei master)
    Ts, Ms, Ks = serep(M, K, Phi, mast)
    r["serep_K"] = Ks
    r["serep_M"] = Ms
    r["cond_Phi_m"] = float(np.linalg.cond(Phi[np.ix_(mast, range(n_m))]))

    # --- FRF coppia -> encoder e coppia -> motore
    modi_k = list(range(n_m))  # rigido + (n_m - 1) flessibili: stesso ordine di Guyan
    r["frf"] = {}
    for tag, i_out in (("enc", ax.i_enc), ("mot", ax.i_mot)):
        G_full = frf_fisica(M, C, K, ax.i_mot, i_out, w)
        G_chk = frf_modale(w2, Phi, ax.beta, ax.i_mot, i_out, w, range(n))
        i_g = {mast.index(ax.i_mot), mast.index(i_out)}
        G_g = frf_fisica(Mg, ax.beta * Kg, Kg, mast.index(ax.i_mot), mast.index(i_out), w)
        G_l = frf_fisica(Ml, ax.beta * Kg, Kg, mast.index(ax.i_mot), mast.index(i_out), w)
        G_m = frf_modale(w2, Phi, ax.beta, ax.i_mot, i_out, w, modi_k)
        G_ms = frf_modale(w2, Phi, ax.beta, ax.i_mot, i_out, w, modi_k, statica=True)
        G_s = frf_fisica(Ms, ax.beta * Ks, Ks, mast.index(ax.i_mot), mast.index(i_out), w)
        cb = {}
        for k in n_cb:
            Tc, Mc, Kc = rd.riduci_craig_bampton(M, K, mast, k)
            nn = Mc.shape[0]
            bi = np.zeros(nn); bi[mast.index(ax.i_mot)] = 1.0
            Cc = ax.beta * Kc
            Gc = np.array([np.linalg.solve(Kc - wk**2 * Mc + 1j * wk * Cc, bi)[mast.index(i_out)] for wk in w])
            cb[k] = Gc
        r["frf"][tag] = dict(full=G_full, check=G_chk, guyan=G_g, lumped=G_l, modale=G_m,
                             modale_corr=G_ms, serep=G_s, cb=cb)
    # --- limiti di validita'
    r["f_lim"] = {tag: {k: f_limite(w, v, d["full"]) for k, v in d.items() if k not in ("full", "cb", "check")}
                  | {f"CB{k}": f_limite(w, v, d["full"]) for k, v in d["cb"].items()}
                  for tag, d in r["frf"].items()}
    r["f_lim_db"] = {tag: {k: f_limite_db(w, v, d["full"]) for k, v in d.items() if k not in ("full", "cb", "check")}
                     | {f"CB{k}": f_limite_db(w, v, d["full"]) for k, v in d["cb"].items()}
                     for tag, d in r["frf"].items()}
    # --- zeri non smorzati di tau -> motore e tau -> encoder
    r["zeri"] = {}
    kept = list(range(n_m))
    for tag, i_out in (("mot", ax.i_mot), ("enc", ax.i_enc)):
        res_full = Phi[i_out, :] * Phi[ax.i_mot, :]
        Pg_all = T @ Pg
        res_g = Pg_all[i_out, :] * Pg_all[ax.i_mot, :]
        D_corr = flessibilita(w2, Phi, i_out, ax.i_mot, range(n_m, n))
        r["zeri"][tag] = dict(
            full=zeri_da_modi(w2, res_full),
            guyan=zeri_da_modi(w2g, res_g),
            modale=zeri_da_modi(w2[kept], res_full[kept]),
            modale_corr=zeri_da_modi(w2[kept], res_full[kept], D_corr),
        )
    # --- flessibilita' statica elastica: collocata (motore) e non collocata (encoder da motore)
    Pg_full = T @ Pg
    r["flex"] = {}
    for tag, i_out in (("mot", ax.i_mot), ("enc", ax.i_enc)):
        f_full_ = flessibilita(w2, Phi, i_out, ax.i_mot, range(1, n))
        f_mod = flessibilita(w2, Phi, i_out, ax.i_mot, range(1, n_m))
        f_g = sum(Pg_full[i_out, i] * Pg_full[ax.i_mot, i] / w2g[i] for i in range(1, n_m))
        r["flex"][tag] = dict(full=f_full_, guyan=f_g, modale=f_mod)
    # --- norme di Hankel per modo (Gawronski 4.6, uscita di posizione): tau -> encoder
    om = np.sqrt(w2[1:])
    z = 0.5 * ax.beta * om
    h = np.abs(Phi[ax.i_mot, 1:] * Phi[ax.i_enc, 1:]) / (4 * z * om**2)
    r["hankel"] = h
    r["ordine_hankel"] = np.argsort(-h) + 1
    r["n_modi_flessibili_full"] = n - 1
    # struttura di K_r
    r["k_catena_guyan"] = np.array([-Kg[i, i + 1] for i in range(n_m - 1)])
    r["Kg"] = Kg
    return r


# ---------------------------------------------------------------------------
# Rapporto e figura
# ---------------------------------------------------------------------------


def rapporto(res: list[dict], mass_model: str, theta0: float, note: list[str]) -> str:
    L: list[str] = []
    L.append(f"# Guyan contro troncamento modale (Gawronski): modello di massa {mass_model}, theta0 = {theta0:g} deg")
    L.append("\n> Risultato da parametri segnaposto: vale per i dati attuali di config.yaml.\n")
    for n in note:
        L.append(f"- {n}")
    for r in res:
        ax: Asse = r["asse"]
        n_m = len(ax.masters)
        L.append(f"\n## Asse {ax.nome}  ({ax.M.shape[0]} GdL completi, {n_m} master: "
                 f"{[ax.etichette[i] for i in ax.masters]})")
        L.append(f"_{ax.nota}_\n")
        L.append("### Frequenze (Hz)\n")
        hdr = "| modo | completo | Guyan | err % | Guyan lumped | err % | " + " | ".join(
            f"CB n={k}" for k in r["f_cb"]) + " | MAC Guyan |"
        L.append(hdr)
        L.append("|" + "---|" * (hdr.count("|") - 1))
        nrow = max(n_m + 2, 6)
        for i in range(1, nrow):
            ff = r["f_full"][i]
            if i >= n_m:
                L.append(f"| {i} | {ff:.1f} | - | - | - | - | " + " | ".join(
                    (f"{f[i]:.1f} ({(f[i]/ff-1)*100:+.2f}%)" if i < len(f) else "-") for f in r["f_cb"].values()) + " | - |")
                continue
            fg = r["f_guyan"][i]
            fl = r["f_lumped"][i]
            cbs = []
            for k, f in r["f_cb"].items():
                cbs.append(f"{f[i]:.1f} ({(f[i]/ff-1)*100:+.2f}%)" if i < len(f) else "-")
            mac = f"{r['mac_guyan'][i]:.4f}" if i < n_m else "-"
            L.append(f"| {i} | {ff:.1f} | {fg:.1f} | {(fg/ff-1)*100:+.2f} | {fl:.1f} | {(fl/ff-1)*100:+.2f} | "
                     + " | ".join(cbs) + f" | {mac} |")
        L.append(f"\nModi flessibili del modello completo: {r['n_modi_flessibili_full']}; "
                 f"modi flessibili che Guyan puo' rappresentare: {n_m - 1}. "
                 f"Il troncamento modale (e SEREP) riproduce i primi {n_m - 1} esattamente (errore 0 per costruzione).")
        L.append(f"\nCondizionamento di Phi sui master (SEREP): {r['cond_Phi_m']:.3g}.")
        L.append("\nK_r di SEREP (fuori diagonale, N m/rad):\n```\n" + np.array2string(r["serep_K"], precision=3) + "\n```")
        L.append("K_r di Guyan (catena di molle in serie, tridiagonale):\n```\n" + np.array2string(r["Kg"], precision=3) + "\n```")
        L.append("\n### Flessibilita' statica elastica F = sum_{i>=1} phi_out phi_mot / w_i^2 (rad/(N m))\n")
        L.append("| uscita | completo | Guyan (err) | modale puro (err) | modale + correzione statica |\n|---|---|---|---|---|")
        for tag, lab in (("mot", "motore (collocata)"), ("enc", "encoder (non collocata)")):
            c = r["flex"][tag]
            L.append(f"| {lab} | {c['full']:.4e} | {c['guyan']:.4e} ({(c['guyan']/c['full']-1)*100:+.3g}%) | "
                     f"{c['modale']:.4e} ({(c['modale']/c['full']-1)*100:+.3g}%) | {c['full']:.4e} (0%) |")
        L.append("\n### Primi zeri non smorzati (antirisonanze), Hz\n")
        L.append("| TF | completo | Guyan | modale puro | modale + correzione statica |\n|---|---|---|---|---|")
        for tag, lab in (("mot", "tau -> motore"), ("enc", "tau -> encoder")):
            z = r["zeri"][tag]
            for j in range(min(3, len(z["full"]))):
                def cella(k):
                    if j >= len(z[k]):
                        return "-"
                    return f"{z[k][j]:.4f} ({(z[k][j]/z['full'][j]-1)*100:+.3g}%)" if k != "full" else f"{z[k][j]:.4f}"
                L.append(f"| {lab} #{j+1} | {cella('full')} | {cella('guyan')} | {cella('modale')} | {cella('modale_corr')} |")
        L.append("\n### Banda di validita' delle FRF, Hz\n")
        L.append("Severo: primo scarto > 1 dB o > 5 deg. Repo: primo scarto di modulo > 3 dB. "
                 "Vicino a un'antirisonanza stretta il rapporto fra FRF e' ipersensibile a uno spostamento minimo dello zero: "
                 "leggere insieme alla tabella degli zeri.\n")
        L.append("| modello | enc, severo | enc, 3 dB | motore, severo | motore, 3 dB |\n|---|---|---|---|---|")
        nomi = {"guyan": "Guyan", "lumped": "Guyan lumped", "modale": f"modale puro (rigido + {n_m-1})",
                "modale_corr": f"modale + correzione statica (rigido + {n_m-1})", "serep": "SEREP"}
        nomi |= {f"CB{k}": f"Craig-Bampton n={k}" for k in r["f_cb"]}
        for k, lab in nomi.items():
            L.append(f"| {lab} | {r['f_lim']['enc'][k]:.1f} | {r['f_lim_db']['enc'][k]:.1f} | "
                     f"{r['f_lim']['mot'][k]:.1f} | {r['f_lim_db']['mot'][k]:.1f} |")
        L.append("\n### Norme di Hankel dei modi flessibili (tau -> encoder), ordinate\n")
        L.append("| ordine | modo (per frequenza) | f (Hz) | h_i (rad/(N m)) |\n|---|---|---|---|")
        for pos, i in enumerate(r["ordine_hankel"][:6], 1):
            L.append(f"| {pos} | {i} | {r['f_full'][i]:.1f} | {r['hankel'][i-1]:.3e} |")
        L.append("\nControllo interno: FRF modale con tutti i modi contro FRF fisica del completo: "
                 f"max errore relativo {max(np.max(np.abs(d['check']/d['full']-1)) for d in r['frf'].values()):.2e}; "
                 f"FRF SEREP contro modale puro: "
                 f"{max(np.max(np.abs(d['serep']/d['modale']-1)) for d in r['frf'].values()):.2e}.")
    # spettro unico dei due assi (caso disaccoppiato: Omega0 = 0, payload bilanciato)
    L.append("\n## Spettro dei due assi insieme (disaccoppiati): ordine dei modi alla Gawronski\n")
    righe = []
    for r in res:
        for i, f in enumerate(r["f_full"]):
            righe.append((f, r["asse"].nome, i, i < len(r["asse"].masters)))
    righe.sort(key=lambda t: (t[0], t[1]))
    L.append("| # | f (Hz) | asse | modo | tenuto da Guyan |\n|---|---|---|---|---|")
    for k, (f, nome, i, tenuto) in enumerate(righe[:14], 1):
        tipo = "rigido" if i == 0 else f"flessibile {i}"
        L.append(f"| {k} | {f:.1f} | {nome} | {tipo} | {'si' if tenuto else 'no'} |")
    n_full = sum(len(r["f_full"]) - 1 for r in res)
    n_guy = sum(len(r["asse"].masters) - 1 for r in res)
    L.append(f"\nModi flessibili totali del FEM completo: {n_full}; rappresentabili con Guyan: {n_guy}.")
    return "\n".join(L)


def figura(res: list[dict], w: np.ndarray, percorso: Path, titolo: str):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, SEC, GRID = "#0b0b0b", "#52514e", "#e3e2de"
    COL = dict(full=INK, guyan="#eb6834", modale="#2a78d6", modale_corr="#1baf7a")
    LAB = dict(full="FEM completo", guyan="Guyan", modale="Gawronski (modale)", modale_corr="Gawronski + correzione statica")
    f = w / (2 * np.pi)
    fig, axs = plt.subplots(2, 2, figsize=(11, 6.8), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    for c, r in enumerate(res):
        ax: Asse = r["asse"]
        d = r["frf"]["enc"]
        a0, a1 = axs[0, c], axs[1, c]
        for k in ("full", "guyan", "modale", "modale_corr"):
            lw = 2.0 if k == "full" else 1.5
            ls = "-" if k in ("full", "guyan", "modale_corr") else (0, (5, 2))
            a0.loglog(f, np.abs(d[k]), color=COL[k], lw=lw, ls=ls, label=LAB[k])
            if k != "full":
                a1.semilogx(f, 20 * np.log10(np.abs(d[k] / d["full"])), color=COL[k], lw=1.4, ls=ls)
        for fr in r["f_full"][1:len(ax.masters)]:
            a0.axvline(fr, color=GRID, lw=0.8, zorder=0)
            a1.axvline(fr, color=GRID, lw=0.8, zorder=0)
        a0.set_title(f"Asse {ax.nome}: coppia motore → encoder", loc="left", fontsize=11, color=INK)
        a0.set_ylabel("|G| [rad/(N m)]", color=SEC)
        a1.set_ylabel("errore [dB]", color=SEC)
        a1.set_xlabel("frequenza [Hz]", color=SEC)
        a1.set_ylim(-12, 12)
        a1.axhline(0, color=SEC, lw=0.6)
        for a in (a0, a1):
            a.grid(True, color=GRID, lw=0.6)
            for s in ("top", "right"):
                a.spines[s].set_visible(False)
            a.tick_params(colors=SEC, labelsize=9)
        a0.set_xlim(f[0], f[-1])
    axs[0, 0].legend(frameon=False, fontsize=9, loc="lower left")
    fig.suptitle(titolo + "  -  risultato da parametri segnaposto", x=0.01, ha="left", fontsize=10, color=SEC)
    fig.tight_layout()
    percorso.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(percorso, dpi=200)
    return percorso


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=Path("config.yaml"))
    ap.add_argument("--theta0", type=float, default=0.0, help="elevazione di lavoro, gradi")
    ap.add_argument("--mass-model", choices=("lumped", "consistent"), default="consistent")
    ap.add_argument("--payload-manuale", action="store_true")
    ap.add_argument("--out", type=Path, default=Path("results/guyan_vs_modale"))
    a = ap.parse_args()

    par, assi = carica_assi(a.config, a.theta0, a.mass_model, a.payload_manuale)
    w = 2 * np.pi * np.logspace(0, np.log10(4e4), 8000)
    res = [analizza(ax, w) for ax in assi]
    note = [f"configurazione: {a.config}", f"smorzamento: C = beta K, zeta = {par.zeta:g} sul primo modo flessibile (segnaposto)"]
    txt = rapporto(res, a.mass_model, a.theta0, note)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "rapporto.md").write_text(txt, encoding="utf-8")
    figura(res, w, a.out / "bode_guyan_vs_modale.png", f"Modello di massa {a.mass_model}, theta0 = {a.theta0:g} deg")
    print(txt)


if __name__ == "__main__":
    main()
