"""Accoppiamento fra gli assi e bilanciamento del telescopio dai doppietti di posizione (VClosedLoop/V2, 20/10/2026).

Richiesta dell'utente (2026-10-28). Ogni file ``MPCLucy_doublets_<asse>_<A>deg_minus<A>deg.csv`` registra a 2 kHz,
con entrambi gli assi in anello chiuso, un doppietto 0 -> +A -> -A -> 0 su un asse (rampe a ~10 deg/s, plateau di
30-45 s) mentre l'altro e' comandato fermo. Canali: posizione e comando dei due assi (millesimi di grado), corrente
e coppia dei due motori (AXIS1 = azimut, AXIS2 = polare; la coppia registrata e' 0.18475 x corrente). Le colonne
hanno istanti diversi all'inizio e alla fine; le righe complete sono simultanee (coppia e corrente dello stesso asse
coincidono con correlazione 1 sulla stessa riga).

Analisi:

* accoppiamento dinamico: spettri incrociati fra la coppia dell'asse mosso e la coppia dell'altro asse (che il suo
  anello usa per tenerlo fermo) nelle finestre dei transitori; coerenza e rapporto. Un modo comune ai due assi da'
  un picco di coerenza alla sua frequenza;
* accoppiamento a bassa frequenza: deviazione dell'encoder dell'altro asse durante ogni movimento, confrontata con
  il verso del movimento e con la coppia di mantenimento dell'altro asse prima del movimento;
* bilanciamento: coppia di mantenimento sui plateau (ultimi 3 s) contro l'angolo e il verso di arrivo,
  T = c + f s + K theta (s = +1/-1 verso di arrivo): f banda d'attrito, K rigidezza di gravita' (o di cavi), c
  sbilanciamento a theta = 0. Anche la coppia di azimut durante i doppietti polari (accoppiamento statico).

Esecuzione (da sim/gimbal_multidof)::

    .venv/bin/python -I analisi_v2_accoppiamento.py
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.signal import csd, welch

RADICE = Path(__file__).resolve().parent
sys.path.insert(0, str(RADICE / "src"))

import matplotlib.pyplot as plt  # noqa: E402

from gimbal_mdof import plotting as pt  # noqa: E402

DATI = RADICE.parent / "LucyIdentification" / "VClosedLoop" / "V2"
CARTELLA = RADICE / "results_V2_coupling"
FS = 2000.0
"""Frequenza di campionamento del recorder, Hz."""
MDEG = np.pi / 180.0 / 1000.0
"""Millesimo di grado in rad (unita' delle posizioni registrate)."""
LSB_MDEG = 360.0 / 2 ** 20 * 1000.0
"""Passo dell'encoder (20 bit), millesimi di grado."""
FINESTRA = (0.2, 2.0)
"""Finestra dei transitori attorno all'inizio di ogni movimento (prima, dopo), s."""
NPERSEG = 2048
PLATEAU_FINE = 3.0
"""Coppia di mantenimento: media sugli ultimi secondi di ogni plateau."""
A_MIN_FIT = 1.0
"""Ampiezza minima dei plateau nel fit del bilanciamento, deg: sotto, l'asse resta nel pre-scorrimento."""
M_TELESCOPIO = 15.3 + 2.05
"""Massa che ruota attorno all'asse polare (telescopio + mozzo e anello), kg."""
RISONANZE_FRF = {"elevation, 1 % grid, 28 Oct": 56.3, "azimuth, 28 Oct": 55.7, "elevation, 20 Oct": 57.4}
"""Risonanze misurate sulle piante (Hz) per confronto nelle figure."""
ASSI = {
    "azimut": {"cmd": "psi_c", "pos": "psi", "tq": "tq_az", "altro_pos": "th", "altro_tq": "tq_el"},
    "elevation": {"cmd": "th_c", "pos": "th", "tq": "tq_el", "altro_pos": "psi", "altro_tq": "tq_az"},
}


def carica(percorso: Path) -> dict[str, np.ndarray]:
    with open(percorso) as fh:
        r = csv.reader(fh)
        h = [c.strip() for c in next(r)]
        righe = [x for x in r]
    a = np.full((len(righe), len(h)), np.nan)
    for i, x in enumerate(righe):
        for j, v in enumerate(x[: len(h)]):
            if v.strip():
                a[i, j] = float(v)
    a = a[~np.isnan(a).any(axis=1)]
    k = {c: j for j, c in enumerate(h)}
    return {"t": (a[:, 0] - a[0, 0]) / 1000.0,
            "psi": a[:, k["AzimuthAxis/ActualPosition"]] * MDEG, "th": a[:, k["PolarAxis/ActualPosition"]] * MDEG,
            "psi_c": a[:, k["AzimuthAxis/CommandPosition"]] * MDEG, "th_c": a[:, k["PolarAxis/CommandPosition"]] * MDEG,
            "tq_az": a[:, k["VM/VARIABLESTOMONITOR.AXIS1_TORQUE"]], "tq_el": a[:, k["VM/VARIABLESTOMONITOR.AXIS2_TORQUE"]]}


def prove(asse: str) -> list[tuple[float, str, dict[str, np.ndarray]]]:
    """``(ampiezza [deg], nome, dati)`` dei doppietti di ``asse``; ``01deg`` = 0.1 deg."""
    out = []
    for p in sorted(DATI.glob(f"MPCLucy_doublets_{asse}_*.csv")):
        m = re.search(r"_(\d+)deg_", p.name)
        assert m is not None
        a = float(m.group(1)) / (10.0 if m.group(1).startswith("0") else 1.0)
        out.append((a, p.name, carica(p)))
    return sorted(out, key=lambda x: x[0])


def movimenti(d: dict[str, np.ndarray], cmd: str) -> list[tuple[int, int]]:
    """Indici ``(inizio, fine)`` delle rampe del comando (separate da piu' di 1 s di comando fermo)."""
    t = d["t"]
    ch = np.flatnonzero(np.abs(np.diff(d[cmd])) > 0)
    tagli = np.flatnonzero(np.diff(t[ch]) > 1.0)
    inizi = np.r_[ch[0], ch[tagli + 1]]
    fini = np.r_[ch[tagli], ch[-1]] + 1
    return list(zip(inizi.tolist(), fini.tolist()))


# ---------------------------------------------------------------------------
# Accoppiamento dinamico
# ---------------------------------------------------------------------------


def spettri(asse: str, dati) -> dict[str, np.ndarray]:
    """Spettri medi sulle finestre dei transitori: coppia dell'asse mosso (x), coppia dell'altro asse (y), loro
    incrociato, e coppia dell'altro asse a riposo (ultimi 8 s di ogni file)."""
    c = ASSI[asse]
    Sxx = Syy = Sxy = Sq = 0.0
    n = nq = 0
    for _, _, d in dati:
        t = d["t"]
        for i0, _ in movimenti(d, c["cmd"]):
            w = (t >= t[i0] - FINESTRA[0]) & (t <= t[i0] + FINESTRA[1])
            x, y = (d[k][w] - d[k][w].mean() for k in (c["tq"], c["altro_tq"]))
            f, pxx = welch(x, FS, nperseg=NPERSEG, nfft=4 * NPERSEG)
            pyy = welch(y, FS, nperseg=NPERSEG, nfft=4 * NPERSEG)[1]
            pxy = csd(x, y, FS, nperseg=NPERSEG, nfft=4 * NPERSEG)[1]
            Sxx, Syy, Sxy, n = Sxx + pxx, Syy + pyy, Sxy + pxy, n + 1
        q = t > t[-1] - 8.0
        Sq, nq = Sq + welch(d[c["altro_tq"]][q] - d[c["altro_tq"]][q].mean(), FS, nperseg=NPERSEG,
                            nfft=4 * NPERSEG)[1], nq + 1
    Sxx, Syy, Sxy = Sxx / n, Syy / n, Sxy / n
    return {"f": f, "Sxx": Sxx, "Syy": Syy, "Sxy": Sxy, "Sq": Sq / nq, "coh": np.abs(Sxy) ** 2 / (Sxx * Syy),
            "H": Sxy / Sxx, "finestre": n}


def picco_coerenza(s: dict[str, np.ndarray], banda=(40.0, 80.0)) -> dict[str, float]:
    ok = (s["f"] >= banda[0]) & (s["f"] <= banda[1])
    i = np.flatnonzero(ok)[int(np.argmax(s["coh"][ok]))]
    j = np.flatnonzero(ok)[int(np.argmax(np.abs(s["H"][ok])))]
    return {"f_coh_Hz": float(s["f"][i]), "coh": float(s["coh"][i]), "ratio_at_coh": float(abs(s["H"][i])),
            "power_ratio_vs_rest": float(s["Syy"][i] / s["Sq"][i]), "f_ratio_max_Hz": float(s["f"][j]),
            "ratio_max": float(abs(s["H"][j]))}


# ---------------------------------------------------------------------------
# Accoppiamento a bassa frequenza: deviazione dell'altro asse in ogni movimento
# ---------------------------------------------------------------------------


def deviazioni(asse: str, dati) -> list[dict]:
    c = ASSI[asse]
    out = []
    for a, nome, d in dati:
        t = d["t"]
        for i0, i1 in movimenti(d, c["cmd"]):
            o = d[c["altro_pos"]] / MDEG
            w = (t >= t[i0]) & (t <= t[i1] + 0.6)
            dev = o[w] - o[i0 - 10]
            j5 = min(int(np.searchsorted(t, t[i1] + 5.0)), len(t) - 1)
            out.append({"file": nome, "A_deg": a, "verso": int(np.sign(d[c["cmd"]][i1] - d[c["cmd"]][i0])),
                        "coppia_prima": float(d[c["altro_tq"]][max(i0 - 200, 0): i0].mean()) if i0 > 0 else np.nan,
                        "coppia_dopo_5s": float(d[c["altro_tq"]][j5]),
                        "deviazione_mdeg": float(dev[int(np.argmax(np.abs(dev)))])})
    return out


# ---------------------------------------------------------------------------
# Bilanciamento: coppia di mantenimento sui plateau
# ---------------------------------------------------------------------------


def plateau(asse: str, dati) -> list[dict]:
    c = ASSI[asse]
    out = []
    for a, nome, d in dati:
        t = d["t"]
        mv = movimenti(d, c["cmd"])
        for k, (i0, i1) in enumerate(mv):
            j1 = mv[k + 1][0] if k + 1 < len(mv) else len(t)
            seg = slice(i1, j1)
            tt = t[seg]
            if tt.size == 0 or tt[-1] - tt[0] < PLATEAU_FINE + 1.0:
                continue
            fine = tt > tt[-1] - PLATEAU_FINE
            inizio = tt < tt[0] + 1.0
            out.append({"file": nome, "A_deg": a, "angolo_deg": float(d[c["cmd"]][i1] / MDEG / 1000.0),
                        "verso": int(np.sign(d[c["cmd"]][i1] - d[c["cmd"]][i0])), "durata_s": float(tt[-1] - tt[0]),
                        "coppia_inizio": float(d[c["tq"]][seg][inizio].mean()),
                        "coppia": float(d[c["tq"]][seg][fine].mean()),
                        "coppia_altro": float(d[c["altro_tq"]][seg][fine].mean()),
                        "t": tt - tt[0], "storia": d[c["tq"]][seg]})
    return out


def fit_mantenimento(pl: list[dict], chiave: str = "coppia", con_verso: bool = True) -> dict[str, float]:
    """Minimi quadrati di ``T = c + f s + K theta`` (theta in rad) sui plateau con ampiezza >= :data:`A_MIN_FIT`."""
    sel = [p for p in pl if p["A_deg"] >= A_MIN_FIT]
    th = np.deg2rad([p["angolo_deg"] for p in sel])
    s = np.array([p["verso"] for p in sel], dtype=float)
    y = np.array([p[chiave] for p in sel])
    X = np.c_[np.ones_like(th), s, th] if con_verso else np.c_[np.ones_like(th), th]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    cov = np.linalg.inv(X.T @ X) * (r @ r) / max(len(y) - X.shape[1], 1)
    e = np.sqrt(np.diag(cov))
    nomi = ("c", "f", "K") if con_verso else ("c", "K")
    out = {n: float(v) for n, v in zip(nomi, b)}
    out.update({f"sigma_{n}": float(v) for n, v in zip(nomi, e)})
    out["rms_residuo"] = float(np.sqrt(np.mean(r ** 2)))
    out["n"] = len(y)
    return out


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def figura_spettri(sp: dict[str, dict], *, segnaposto: bool):
    fig, assi = plt.subplots(3, 2, figsize=(13.0, 9.0), sharex=True)
    titoli = {"azimut": "azimuth moved, elevation held", "elevation": "elevation moved, azimuth held"}
    altro = {"azimut": "elevation", "elevation": "azimuth"}
    for col, (asse, s) in enumerate(sp.items()):
        f = s["f"]
        a = assi[0, col]
        a.loglog(f, s["Sxx"], color=pt.SERIE[0], lw=1.0, label=f"{'azimuth' if asse == 'azimut' else 'elevation'} "
                 "torque (moved axis)")
        a.loglog(f, s["Syy"], color=pt.SERIE[1], lw=1.0, label=f"{altro[asse]} torque, transients")
        a.loglog(f, s["Sq"], color=pt.INCHIOSTRO_TENUE, lw=0.8, label=f"{altro[asse]} torque, at rest")
        a.set_ylabel("PSD [N² m²/Hz]")
        a.legend(fontsize=7, loc="lower left")
        a.set_title(f"{titoli[asse]} ({s['finestre']} transients)", fontsize=9)
        assi[1, col].loglog(f, np.abs(s["H"]), color=pt.SERIE[2], lw=1.0)
        assi[1, col].set_ylabel("|held / moved| torque")
        assi[1, col].set_ylim(1e-3, 3.0)
        assi[2, col].semilogx(f, s["coh"], color=pt.INCHIOSTRO, lw=1.0)
        assi[2, col].set_ylabel("coherence")
        assi[2, col].set_ylim(0, 1)
        assi[2, col].set_xlabel("frequency [Hz]")
        for a in assi[:, col]:
            a.set_xlim(2.0, 1000.0)
            for k, fr in enumerate(RISONANZE_FRF.values()):
                a.axvline(fr, color=pt.SERIE[3], lw=0.7, ls=(":", "--", "-.")[k], alpha=0.8)
    assi[2, 0].plot([], [], color=pt.SERIE[3], lw=0.7, ls=":", label="plant resonances: elevation 56.3, "
                    "azimuth 55.7, elevation 20 Oct 57.4 Hz")
    assi[2, 0].legend(fontsize=7, loc="upper right")
    pt._titolo(fig, "Cross-axis coupling during the position doublets (VClosedLoop/V2, 20 Oct 2026)", segnaposto)
    return fig


def figura_rilascio(dev: dict[str, list[dict]], esempio: dict[str, np.ndarray], *, segnaposto: bool):
    fig = plt.figure(figsize=(13.0, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.4])
    a = fig.add_subplot(gs[:, 0])
    for v, col, mk in ((1, pt.SERIE[0], "^"), (-1, pt.SERIE[1], "v")):
        sel = [r for r in dev["azimut"] if r["verso"] == v and np.isfinite(r["coppia_prima"])]
        a.scatter([r["coppia_prima"] for r in sel], [r["deviazione_mdeg"] for r in sel], color=col, marker=mk, s=28,
                  label=f"azimuth moving {'+' if v > 0 else '-'}")
    a.axhline(0, color=pt.INCHIOSTRO_TENUE, lw=0.7)
    a.axvline(0, color=pt.INCHIOSTRO_TENUE, lw=0.7)
    a.set_xlabel("elevation holding torque before the azimuth move [N m]")
    a.set_ylabel("largest elevation encoder deviation [mdeg]")
    a.legend(fontsize=7)
    a.set_title("elevation slips the way its own motor was pushing", fontsize=9)
    t = esempio["t"]
    i0 = movimenti(esempio, "psi_c")[0][0]
    w = (t >= t[i0] - 0.5) & (t <= t[i0] + 8.0)
    b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 1], sharex=b)
    b.plot(t[w] - t[i0], esempio["th"][w] / MDEG, color=pt.SERIE[1], lw=0.9)
    b.set_ylabel("elevation encoder [mdeg]")
    b.tick_params(labelbottom=False)
    c.plot(t[w] - t[i0], esempio["tq_el"][w], color=pt.INCHIOSTRO, lw=0.9, label="elevation torque [N m]")
    c.plot(t[w] - t[i0], np.rad2deg(esempio["psi"][w]) / 10.0, color=pt.SERIE[0], lw=0.9, ls="--",
           label="azimuth position / 10 [deg]")
    c.set_ylabel("[N m], [deg / 10]")
    c.set_xlabel("time from the start of the azimuth move [s]")
    c.legend(fontsize=7, loc="lower right")
    b.set_title("example: azimuth 0 -> +5 deg", fontsize=9)
    pt._titolo(fig, "Low-frequency interaction: an azimuth move releases the elevation stiction", segnaposto)
    return fig


def figura_bilanciamento(pl_el: list[dict], fit_el: dict, fit_az_x: dict, *, segnaposto: bool):
    fig, assi = plt.subplots(1, 3, figsize=(15.5, 4.8))
    a = assi[0]
    th = np.linspace(-16, 16, 50)
    for v, col, mk in ((1, pt.SERIE[0], "^"), (-1, pt.SERIE[1], "v")):
        sel = [p for p in pl_el if p["verso"] == v]
        big = [p for p in sel if p["A_deg"] >= A_MIN_FIT]
        small = [p for p in sel if p["A_deg"] < A_MIN_FIT]
        a.scatter([p["angolo_deg"] for p in big], [p["coppia"] for p in big], color=col, marker=mk, s=30,
                  label=f"arrived moving {'+' if v > 0 else '-'}")
        a.scatter([p["angolo_deg"] for p in small], [p["coppia"] for p in small], facecolors="none", edgecolors=col,
                  marker=mk, s=30)
        a.plot(th, fit_el["c"] + v * fit_el["f"] + fit_el["K"] * np.deg2rad(th), color=col, lw=1.0)
    a.set_xlabel("elevation angle [deg]")
    a.set_ylabel("elevation holding torque, last 3 s [N m]")
    a.set_title(f"T = {fit_el['c']:+.3f} ± {fit_el['f']:.3f} + {fit_el['K']:.2f}·θ  (N m, θ rad)", fontsize=9)
    a.legend(fontsize=7, title="open: 0.1 and 0.5 deg (not fitted)", title_fontsize=7)
    b = assi[1]
    for p in pl_el:
        if p["A_deg"] >= A_MIN_FIT:
            col = pt.SERIE[0] if p["verso"] > 0 else pt.SERIE[1]
            b.plot(p["t"], p["storia"], color=col, lw=0.5, alpha=0.6)
    b.set_xlabel("time on the plateau [s]")
    b.set_ylabel("elevation torque [N m]")
    b.set_title("holding torque after each move (blue: arrived +, orange: arrived -)", fontsize=9)
    c = assi[2]
    c.scatter([p["angolo_deg"] for p in pl_el], [p["coppia_altro"] for p in pl_el], color=pt.SERIE[2], s=26)
    c.plot(th, fit_az_x["c"] + fit_az_x["K"] * np.deg2rad(th), color=pt.SERIE[2], lw=1.0)
    c.set_xlabel("elevation angle [deg]")
    c.set_ylabel("azimuth holding torque, last 3 s [N m]")
    c.set_title(f"azimuth torque against elevation angle: {fit_az_x['K']:.2f} N m/rad", fontsize=9)
    pt._titolo(fig, "Holding torque on the plateaus of the elevation doublets: friction band, unbalance, "
                    "static coupling", segnaposto)
    return fig


def main() -> None:
    pt.applica_stile()
    seg = False
    for d in ("figures", "tables"):
        (CARTELLA / d).mkdir(parents=True, exist_ok=True)
    dati = {asse: prove(asse) for asse in ASSI}

    sp = {asse: spettri(asse, dati[asse]) for asse in ASSI}
    picchi = {asse: picco_coerenza(s) for asse, s in sp.items()}
    dev = {asse: deviazioni(asse, dati[asse]) for asse in ASSI}
    pl = {asse: plateau(asse, dati[asse]) for asse in ASSI}
    fit_el = fit_mantenimento(pl["elevation"])
    fit_az = fit_mantenimento(pl["azimut"])
    fit_az_x = fit_mantenimento(pl["elevation"], "coppia_altro", con_verso=False)

    # deviazione dell'elevazione: segno uguale al verso dell'azimut o alla coppia polare di prima?
    rr = [r for r in dev["azimut"] if np.isfinite(r["coppia_prima"]) and abs(r["deviazione_mdeg"]) > LSB_MDEG]
    concordi_coppia = sum(np.sign(r["deviazione_mdeg"]) == np.sign(r["coppia_prima"]) for r in rr)
    concordi_verso = sum(np.sign(r["deviazione_mdeg"]) == r["verso"] for r in rr)
    ribaltate = sum(np.sign(r["coppia_dopo_5s"]) != np.sign(r["coppia_prima"]) for r in rr)
    dev_el = [abs(r["deviazione_mdeg"]) for r in dev["elevation"]]

    with open(CARTELLA / "tables" / "plateaus.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["axis", "file", "A_deg", "angle_deg", "arrival", "duration_s", "torque_first_1s_Nm",
                    "torque_last_3s_Nm", "other_axis_torque_last_3s_Nm"])
        for asse, v in pl.items():
            for p in v:
                w.writerow([asse, p["file"], p["A_deg"], f"{p['angolo_deg']:.3f}", p["verso"], f"{p['durata_s']:.1f}",
                            f"{p['coppia_inizio']:.4f}", f"{p['coppia']:.4f}", f"{p['coppia_altro']:.4f}"])
    with open(CARTELLA / "tables" / "cross_deviations.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["moved_axis", "file", "A_deg", "direction", "held_axis_torque_before_Nm",
                    "held_axis_torque_after_5s_Nm", "held_axis_encoder_max_deviation_mdeg"])
        for asse, v in dev.items():
            for r in v:
                w.writerow([asse, r["file"], r["A_deg"], r["verso"], f"{r['coppia_prima']:.4f}",
                            f"{r['coppia_dopo_5s']:.4f}", f"{r['deviazione_mdeg']:.3f}"])

    esempio = next(d for a, n, d in dati["azimut"] if n.startswith("MPCLucy_doublets_azimut_5deg"))
    figure = pt.salva(figura_spettri(sp, segnaposto=seg), CARTELLA / "figures", "cross_coupling_spectra")
    figure += pt.salva(figura_rilascio(dev, esempio, segnaposto=seg), CARTELLA / "figures", "stiction_release")
    figure += pt.salva(figura_bilanciamento(pl["elevation"], fit_el, fit_az_x, segnaposto=seg), CARTELLA / "figures",
                       "holding_torque_balance")

    g = 9.81
    sommario = {
        "coherence_peak": picchi,
        "low_frequency": {"moves": len(rr), "deviation_sign_equals_held_torque_sign": int(concordi_coppia),
                          "deviation_sign_equals_move_direction": int(concordi_verso),
                          "held_torque_sign_flipped": int(ribaltate),
                          "azimuth_encoder_deviation_in_elevation_moves_mdeg_max": float(max(dev_el))},
        "elevation_holding_fit": fit_el, "azimuth_holding_fit": fit_az, "azimuth_torque_vs_elevation_angle": fit_az_x,
        "elevation_unbalance_if_gravity": {"m_kg": M_TELESCOPIO,
                                           "cg_below_axis_mm": 1000 * fit_el["K"] / (M_TELESCOPIO * g),
                                           "cg_along_optical_axis_mm": 1000 * fit_el["c"] / (M_TELESCOPIO * g)},
    }
    (CARTELLA / "summary.json").write_text(json.dumps(sommario, indent=2))
    print(json.dumps(sommario, indent=2))
    print("figure:", *[str(p.relative_to(RADICE)) for p in figure if p.suffix == ".png"])


if __name__ == "__main__":
    main()
