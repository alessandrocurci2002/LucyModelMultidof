"""Confronto con le misure: lettura dei CSV WorkBench, unita', estremi, validazione della configurazione.

I test usano CSV sintetici nel formato di Kollmorgen WorkBench, non i dati reali (che stanno
fuori dal progetto e possono cambiare).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from gimbal_mdof import concentrato as cc
from gimbal_mdof import fem
from gimbal_mdof import linearize as lz
from gimbal_mdof import params as pr
from gimbal_mdof import payload as plm
from gimbal_mdof import sperimentale as spm

from conftest import config_dati, parametri


def _scrivi_csv(percorso: Path, f: np.ndarray, g_db: np.ndarray, p_deg: np.ndarray, corrente: float = 2.5) -> Path:
    intest = ["Frequency[Hz]", "Plant Gain[dB]", "Plant Phase[deg]", "Measured Gain[dB]", "Measured Phase[deg]",
              "", "#ParentName", "#Name", "#Type", "#Value"]
    meta = [("Root", "DataStructureVersion", "String", "V1.0"),
            ("MeasurementConfiguration", "CurrentExcitation", "Double", str(corrente)),
            ("MeasurementConfiguration", "SineExcitationEndFrequency", "Double", str(f[-1])),
            ("AntiResonanceFilters[0]", "NumeratorF", "Double", "65"),
            ("AntiResonanceFilters[0]", "NumeratorQ", "Double", "5"),
            ("AntiResonanceFilters[0]", "DenominatorF", "Double", "65"),
            ("AntiResonanceFilters[0]", "DenominatorQ", "Double", "17"),
            ("AntiResonanceFilters[0]", "ArfGroup", "String", "ArfGroupForward")]
    righe = [",".join(intest)]
    n = max(len(f), len(meta))
    for i in range(n):
        dati = [f"{f[i]}", f"{g_db[i]}", f"{p_deg[i]}", f"{g_db[i]}", f"{p_deg[i]}"] if i < len(f) else [""] * 5
        m = list(meta[i]) if i < len(meta) else ["", "", "", ""]
        righe.append(",".join(dati + [""] + m))
    percorso.write_text("﻿" + "\n".join(righe), encoding="utf-8")
    return percorso


def test_lettura_csv(tmp_path: Path) -> None:
    f = np.logspace(0, 3, 50)
    g = -20 * np.log10(f)
    p = np.full_like(f, -90.0)
    file = _scrivi_csv(tmp_path / "m.csv", f, g, p, corrente=3.5)
    fr, tracce, meta = spm.carica_bode(file)
    assert np.allclose(fr, f) and np.allclose(tracce["Plant"][0], g) and np.allclose(tracce["Plant"][1], p)
    assert meta["MeasurementConfiguration/CurrentExcitation"] == "3.5"
    cfg = pr.ParametriSperimentali(cartella=tmp_path, Kt=0.258, traccia="Plant", f_min=0.1, f_max=2000.0,
                                   banda_livello=(15.0, 300.0), misure=())
    voce = pr.VoceMisura(file="m.csv", asse="azimut", theta0=0.0, telescopio=True, etichetta="prova")
    m = spm.carica_misura(cfg, voce)
    assert m.eccitazione_A == 3.5 and m.filtri == ((65.0, 5.0, 65.0, 17.0),)
    assert np.allclose(np.abs(m.G), 10 ** (g / 20))


def test_unita_velocita_su_corrente(par: pr.Parametri) -> None:
    """G_meas = j omega K_t theta/tau; nel tratto rigido |G| = K_t / (J omega)."""
    assi = (fem.costruisci_asse(par.azimut), fem.costruisci_asse(par.polare))
    s = cc.costruisci_sistema_concentrato(par, plm.costruisci_payload(par, assi[1]), assi)
    ml = lz.linearizza(s.accoppiato, 0.0, 0.0)
    f = np.array([0.2, 0.5])  # ben sotto l'antirisonanza carico-molla (~11 Hz)
    G = spm.frf_velocita(ml, "azimut", f, 0.258)
    J = float(np.ones(3) @ s.rid_psi.M_r @ np.ones(3)) + s.payload.J_ax + float(s.payload.J_psi(0.0))
    assert np.allclose(np.abs(G["psi_e"]), 0.258 / (J * 2 * np.pi * f), rtol=5e-3)
    assert np.allclose(np.angle(G["psi_e"], deg=True), -90.0, atol=2.0)


def test_estremi() -> None:
    f = np.logspace(0, 3, 400)
    s = 2j * np.pi * f
    w1, w2 = 2 * np.pi * 50, 2 * np.pi * 300
    G = (s**2 + 2 * 0.02 * w1 * s + w1**2) / (s**2 + 2 * 0.02 * w2 * s + w2**2) / s
    e = spm.estremi(f, 20 * np.log10(np.abs(G)))
    tipi = {t: fr for t, fr, _ in e}
    assert abs(tipi["antirisonanza"] / 50 - 1) < 0.03 and abs(tipi["risonanza"] / 300 - 1) < 0.03


def test_configurazione_sperimentale(tmp_path: Path) -> None:
    (tmp_path / "Validation").mkdir()
    f = np.logspace(0, 2, 10)
    _scrivi_csv(tmp_path / "a.csv", f, -f, 0 * f)
    _scrivi_csv(tmp_path / "Validation" / "b.csv", f, -f, 0 * f)
    dati = config_dati()
    dati["experimental"] = {"data_dir": str(tmp_path), "torque_constant_Nm_per_A": 0.258,
                            "measurements": [{"file": "a.csv", "axis": "azimut", "theta0_deg": 30}]}
    p = parametri(dati)
    assert p.sperimentale is not None and p.sperimentale.misure[0].theta0 == pytest.approx(np.deg2rad(30))
    dati["experimental"]["measurements"].append({"file": "Validation/b.csv", "axis": "polare"})
    with pytest.raises(pr.ErroreConfigurazione, match="Validation"):
        parametri(dati)
    dati["experimental"]["measurements"] = [{"file": "manca.csv", "axis": "azimut"}]
    with pytest.raises(pr.ErroreConfigurazione, match="inesistente"):
        parametri(dati)
