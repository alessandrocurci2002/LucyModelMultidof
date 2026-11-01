"""Modello hub: mozzo polare sul percorso fra i cuscinetti e telescopio appeso al centro del
mozzo tramite k_b (richiesta 2026-10-17).

Si verificano le rigidezze del mozzo (forme chiuse), l'inerzia dal tensore Onshape, la
topologia ad albero del FEM e della riduzione, lo zero all'encoder che il ramo introduce, la
coerenza di Newton-Eulero con Lagrange e che il modello originale resti invariato.
"""

from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from conftest import config_dati, parametri
from gimbal_mdof import concentrato as cc
from gimbal_mdof import coupled as cp
from gimbal_mdof import fem, mozzo as mz, reduction as rd
from gimbal_mdof import linearize as lz
from gimbal_mdof import sperimentale as spm

K_B = 9.8e3
I_HUB = [[0.018808324, -1.5e-8, 2.0e-9], [-1.5e-8, 0.096041034, 2.438e-6], [2.0e-9, 2.438e-6, 0.082369929]]


def _dati_hub() -> dict:
    d = config_dati()
    p = d["polar"]
    p["elements"][2] = {"id": 3, "mode": "telescope", "length_mm": 1.0, "placeholder": False}
    p["payload_clamp_rigid_mm"] = 0.0
    p["hub"] = {
        "enabled": True, "side_elements": [2, 4], "center_element": 3,
        "inertia": {"source": "test", "frame": "cad_cg", "m_kg": 3.534,
                    "r_G_mm": [0.0, 0.0, d["frames"]["o2_z_mm"]], "I_kgm2": I_HUB,
                    "in_azimuth_load": False, "placeholder": False},
        "side": {"tube": {"d_out_mm": 90.0, "d_in_mm": 84.0, "length_mm": 42.0, "material": "acciaio"},
                 "flange": {"d_out_mm": 90.0, "length_mm": 17.0, "material": "acciaio"}},
        "ring": {"d_out_mm": 310.0, "d_in_mm": 300.0, "width_mm": 100.0, "material": "alluminio",
                 "poisson": 0.33},
        "telescope_attachment": {"k_b_Nm_rad": K_B, "source": "test", "placeholder": True},
    }
    return d


@pytest.fixture(scope="module")
def par_hub():
    return parametri(_dati_hub())


def test_flag_disattivato_lascia_il_modello_originale(par) -> None:
    assert par.mozzo is None
    mod = fem.costruisci_asse(par.polare, par.fem.forma_consistente)
    assert not mod.ramo and mod.i_H is None
    d = _dati_hub()
    d["polar"]["hub"]["enabled"] = False
    assert parametri(d).mozzo is None


def test_rigidezze_mozzo_forme_chiuse(par_hub) -> None:
    pm = par_hub.mozzo
    rg = mz.rigidezze(pm)
    G_s, G_a = pm.tubo.materiale.G, pm.anello.materiale.G
    assert rg.k_tubo == pytest.approx(G_s * np.pi * (0.090**4 - 0.084**4) / 32 / 0.042, rel=1e-12)
    assert rg.k_flangia == pytest.approx(G_s * np.pi * 0.090**4 / 32 / 0.017, rel=1e-12)
    R, t, b = 0.1525, 0.005, 0.100
    E = 2 * G_a * 1.33
    beta = (1 - 0.63 * t / b + 0.052 * (t / b) ** 5) / 3
    k_an = 4 / (np.pi * R * (1 / (E * t * b**3 / 12) + 1 / (G_a * beta * b * t**3)))
    assert rg.k_anello == pytest.approx(k_an, rel=1e-12)
    assert 1 / rg.k_lato == pytest.approx(1 / rg.k_tubo + 1 / rg.k_flangia + 1 / (2 * k_an), rel=1e-12)
    assert rg.k_totale == pytest.approx(rg.k_lato / 2, rel=1e-12)


def test_inerzia_mozzo_dal_tensore_onshape(par_hub) -> None:
    """Baricentro in o_2, asse polare = x della posa CAD: J attorno all'asse polare = Lxx."""
    pm = par_hub.mozzo
    assert pm.J_polare == pytest.approx(I_HUB[0][0], rel=1e-9)
    assert pm.J_azimut == pytest.approx(I_HUB[2][2], rel=1e-9)
    im = mz.inerzia(pm, par_hub.terne)
    assert im.scarto_baricentro < 1e-9


def test_fem_ad_albero(par_hub) -> None:
    mod = fem.costruisci_asse(par_hub.polare, par_hub.fem.forma_consistente)
    assert mod.ramo and mod.etichetta(mod.i_P) == "P" and mod.etichette[mod.i_H] == "H"
    # il telescopio P e' collegato solo al centro del mozzo, tramite k_b, e non ha massa propria
    assert mod.K[mod.i_P, mod.i_H] == pytest.approx(-K_B)
    assert np.count_nonzero(mod.K[mod.i_P]) == 2
    assert mod.M["consistent"][mod.i_P, mod.i_P] == 0.0
    assert mod.J_nodali[mod.i_H] == pytest.approx(par_hub.mozzo.J_polare)
    assert np.allclose(mod.K @ np.ones(mod.n), 0.0, atol=1e-6)
    # senza telescopio il ramo sparisce e il centro del mozzo fa da P
    mod0 = fem.costruisci_asse(par_hub.polare, par_hub.fem.forma_consistente, con_telescopio=False)
    assert not mod0.ramo and mod0.n == mod.n - 1


def test_riduzione_ad_albero(par_hub) -> None:
    mod = fem.costruisci_asse(par_hub.polare, par_hub.fem.forma_consistente)
    rid = rd.riduci(mod, "consistent", "guyan")
    assert rid.ruoli == ("e", "h", "m", "r", "p")
    assert rid.archi == ((0, 1), (1, 2), (2, 3), (1, 4))
    assert rid.k_catena[-1] == pytest.approx(K_B)


def test_zero_del_telescopio_all_encoder(par_hub) -> None:
    """Il ramo laterale introduce lo zero sqrt((k_b + k_g)/I_JJ) nel canale coppia -> encoder."""
    s = cp.costruisci_sistema(par_hub, "consistent", "guyan")
    ml = lz.linearizza(s.accoppiato, 0.0, 0.0)
    _, fz = spm.poli_zeri_modello(ml, "theta_e", 1, 1000.0)
    atteso = np.sqrt((K_B + float(s.payload.k_g(0.0))) / s.payload.I_yy) / (2 * np.pi)
    assert np.min(np.abs(fz - atteso)) / atteso < 1e-3


def test_newton_eulero_mozzo_coincide_con_lagrange() -> None:
    assert cc.verifica_contro_lagrange("mozzo") == sp.zeros(8, 1)


def test_concentrato_mozzo(par_hub) -> None:
    assi = (fem.costruisci_asse(par_hub.azimut, par_hub.fem.forma_consistente),
            fem.costruisci_asse(par_hub.polare, par_hub.fem.forma_consistente))
    from gimbal_mdof import payload as plm

    pl = plm.costruisci_payload(par_hub, assi[1])
    s = cc.costruisci_sistema_concentrato(par_hub, pl, assi)
    assert cc.variante_di(s) == "mozzo" and s.accoppiato.n == 8
    ee = cc.elementi_equivalenti(assi[1])
    assert (ee[-1].nodo_a, ee[-1].nodo_b) == ("H", "P") and ee[-1].k == pytest.approx(K_B)
    # il modello senza telescopio (prova senza payload) torna una catena
    s0 = cc.costruisci_sistema_concentrato(par_hub, spm.payload_senza_telescopio(pl),
                                           spm.assi_senza_telescopio(par_hub, assi))
    assert cc.variante_di(s0) == "catena"


def test_tratto_telescopio_sul_mozzo_conserva_l_inerzia(par_hub) -> None:
    """telescope_on_ring: l'inerzia del tratto di tubo passa dal telescopio al mozzo, il totale non cambia."""
    from gimbal_mdof import payload as plm

    d = _dati_hub()
    d["polar"]["hub"]["telescope_on_ring"] = {"enabled": True, "d_out_mm": 300.0, "d_in_mm": 290.0,
                                             "length_mm": 100.0, "material": "alluminio", "placeholder": True}
    par_t = parametri(d)
    m_t, _, J_t = mz.inerzia_tratto(par_t.mozzo)
    ro, ri, L = 0.150, 0.145, 0.100
    rho = par_t.mozzo.tratto_telescopio.materiale.rho
    assert m_t == pytest.approx(rho * np.pi * (ro**2 - ri**2) * L, rel=1e-12)
    assert J_t == pytest.approx(m_t * ((ro**2 + ri**2) / 4 + L**2 / 12), rel=1e-12)
    assert par_t.mozzo.J_polare == pytest.approx(par_hub.mozzo.J_polare + J_t, rel=1e-12)
    mod_a, mod_b = (fem.costruisci_asse(p.polare, p.fem.forma_consistente) for p in (par_hub, par_t))
    pl_a, pl_b = plm.costruisci_payload(par_hub, mod_a), plm.costruisci_payload(par_t, mod_b)
    assert pl_b.I_yy == pytest.approx(pl_a.I_yy - J_t, rel=1e-12)
    assert pl_b.m == pytest.approx(pl_a.m - m_t, rel=1e-12)
    J_tot = lambda mod, pl: float(mod.M["consistent"].sum()) + pl.I_yy  # noqa: E731
    assert J_tot(mod_b, pl_b) == pytest.approx(J_tot(mod_a, pl_a), rel=1e-12)


def test_senza_telescopio_il_tratto_sparisce() -> None:
    d = _dati_hub()
    d["polar"]["hub"]["telescope_on_ring"] = {"enabled": True, "d_out_mm": 300.0, "d_in_mm": 290.0,
                                             "length_mm": 100.0, "material": "alluminio", "placeholder": True}
    par_t = parametri(d)
    mod0 = fem.costruisci_asse(par_t.polare, par_t.fem.forma_consistente, con_telescopio=False)
    J_t = mz.inerzia_tratto(par_t.mozzo)[2]
    assert mod0.J_nodali[mod0.i_P] == pytest.approx(par_t.mozzo.J_polare - J_t, rel=1e-12)


def test_fattore_rigidezza_anello(par_hub) -> None:
    d = _dati_hub()
    d["polar"]["hub"]["ring"]["stiffness_factor"] = 100.0
    rg1, rg100 = mz.rigidezze(par_hub.mozzo), mz.rigidezze(parametri(d).mozzo)
    assert rg100.k_anello == pytest.approx(100.0 * rg1.k_anello, rel=1e-12)
    assert rg100.k_semianelli == pytest.approx(rg1.k_anello, rel=1e-12)
    assert 1 / rg100.k_lato == pytest.approx(1 / rg100.k_tubo + 1 / rg100.k_flangia + 1 / (2 * rg100.k_anello),
                                             rel=1e-12)
