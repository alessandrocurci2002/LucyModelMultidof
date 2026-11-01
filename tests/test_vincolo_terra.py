"""Vincolo verso terra (ground_constraint): stima su dati sintetici e applicazione al modello lineare."""

from __future__ import annotations

import numpy as np
import pytest

from conftest import config_dati, parametri
from gimbal_mdof import concentrato as cc
from gimbal_mdof import fem, linearize as lz, modal as md
from gimbal_mdof import payload as plm
from gimbal_mdof import vincolo_terra as vt

KT = 0.258


def test_stima_ritrova_k_ed_eta() -> None:
    f = np.logspace(0, np.log10(5.0), 40)
    G = vt.risposta(400.0, 1.5, 0.0, 0.9, f, KT)
    s = vt.stima(f, 20 * np.log10(np.abs(G)), np.rad2deg(np.angle(G)), 0.9, KT)
    assert s.k == pytest.approx(400.0, rel=1e-4) and s.eta == pytest.approx(1.5, rel=1e-4)
    assert s.rms_db < 1e-6 and s.f0 == pytest.approx(np.sqrt(400.0 / 0.9) / (2 * np.pi))
    sv = vt.stima(f, 20 * np.log10(np.abs(G)), np.rad2deg(np.angle(G)), 0.9, KT, perdita="viscosa")
    assert sv.c > 0 and sv.eta == 0.0


def _con_vincolo() -> dict:
    d = config_dati()
    d["ground_constraint"] = {"azimuth": {"k_Nm_rad": 1500.0, "eta": 0.3},
                              "polar": {"k_Nm_rad": 200.0, "eta": 3.0},
                              "polar_without_telescope": {"k_Nm_rad": 2000.0, "eta": 1.0}}
    return d


def test_senza_sezione_il_modello_non_cambia(par) -> None:
    assert par.vincolo_terra is None
    assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
            fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    s = cc.costruisci_sistema_concentrato(par, plm.costruisci_payload(par, assi[1]), assi)
    ml = lz.linearizza(s.accoppiato, 0.0, 0.0)
    assert vt.applica(ml, par) is ml


def test_applica_sul_nodo_dei_cuscinetti() -> None:
    par = parametri(_con_vincolo())
    assert par.vincolo_terra.per("polare", False).k == 2000.0
    assi = (fem.costruisci_asse(par.azimut, par.fem.forma_consistente),
            fem.costruisci_asse(par.polare, par.fem.forma_consistente))
    pl = plm.costruisci_payload(par, assi[1])
    # modello ridotto: forcella (L) e nodo del telescopio (P)
    ml = lz.linearizza(cc.costruisci_sistema_concentrato(par, pl, assi).accoppiato, 0.0, 0.0)
    mv = vt.applica(ml, par)
    st = ml.accoppiato.struttura
    dK = mv.K0 - ml.K0
    assert dK[st.j_psi_L, st.j_psi_L] == pytest.approx(1500.0 * (1 + 0.3j))
    assert dK[st.j_theta_p, st.j_theta_p] == pytest.approx(200.0 * (1 + 3.0j))
    assert np.count_nonzero(dK) == 2
    # FEM: nodo 9 dell'azimut e nodo P del polare
    mf = md.costruisci_sistema_modale(par, "consistent", payload=pl, assi=assi, valida=False).lineare_completo(0, 0)
    dKf = vt.applica(mf, par).K0 - mf.K0
    nomi = list(mf.accoppiato.base_psi.nomi) + list(mf.accoppiato.base_theta.nomi)
    assert dKf[nomi.index("psi_9"), nomi.index("psi_9")] == pytest.approx(1500.0 * (1 + 0.3j))
    assert dKf[nomi.index("theta_P"), nomi.index("theta_P")] == pytest.approx(200.0 * (1 + 3.0j))
    # a bassa frequenza l'asse vincolato e' una molla: modulo crescente, fase positiva
    w = 2 * np.pi * np.array([0.3, 0.6])
    G = mv.frf(w, ["theta_e"])[:, 0, 1] * 1j * w
    assert abs(G[1]) > abs(G[0]) and np.all(np.angle(G) > 0)
