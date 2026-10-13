"""Telescopio inserito in ``config.yaml`` (``payload.source: manual``)."""

from __future__ import annotations

import numpy as np
import pytest

from gimbal_mdof import fem
from gimbal_mdof import params as pr
from gimbal_mdof import payload as plm

from conftest import config_dati, parametri


def _payload(dati: dict) -> plm.Payload:
    par = parametri(dati)
    return plm.costruisci_payload(par, fem.costruisci_asse(par.polare))


def _uguali(a: plm.Payload, b: plm.Payload, rel: float = 1e-12) -> None:
    scala = np.max(np.abs(b.I))
    assert np.allclose(a.I, b.I, rtol=0, atol=rel * scala)
    assert a.m == pytest.approx(b.m, rel=rel)
    assert np.allclose(a.d, b.d, rtol=0, atol=rel * np.max(np.abs(b.d)))
    assert a.J_ax == pytest.approx(b.J_ax, rel=rel)


def test_cad_cg_riproduce_il_corpo_B_del_2dof() -> None:
    """Il corpo B del CAD inserito a mano (kg m^2, mm) == sorgente gimbal_2dof."""
    dati_2dof = config_dati()
    dati_2dof["payload"] = {"enabled": True, "source": "gimbal_2dof",
                            "config_2dof": dati_2dof["payload"]["config_2dof"],
                            "subtract_shaft_parts": True}
    _uguali(_payload(config_dati()), _payload(dati_2dof))


def test_body2_o2_equivale_a_cad_cg() -> None:
    dati = config_dati()
    dati["payload"]["subtract_shaft_parts"] = False
    ref = _payload(dati)
    dati_BO = config_dati()
    dati_BO["payload"].update(frame="body2_o2", subtract_shaft_parts=False, I_kgm2=ref.I_2.tolist(),
                              d_mm=(ref.d_2 * 1e3).tolist())
    dati_BO["payload"].pop("r_G_mm")
    _uguali(_payload(dati_BO), ref)


def test_convenzione_dei_prodotti() -> None:
    dati = config_dati()
    I = np.array(dati["payload"]["I_kgm2"], dtype=float)
    fuori = ~np.eye(3, dtype=bool)
    I[fuori] = -I[fuori]
    dati["payload"].update(I_kgm2=I.tolist(), products_sign="positive")
    _uguali(_payload(dati), _payload(config_dati()))


def test_senza_sottrazione_massa_inserita() -> None:
    dati = config_dati()
    dati["payload"]["subtract_shaft_parts"] = False
    pl = _payload(dati)
    assert pl.m == pytest.approx(dati["payload"]["m_kg"], rel=1e-15)
    assert pl.sottratto == {"m": 0.0, "J": 0.0, "trasversale": 0.0}
    con = _payload(config_dati())
    # la sottrazione tocca solo la diagonale: i prodotti restano invariati
    fuori = ~np.eye(3, dtype=bool)
    assert np.allclose(pl.I[fuori], con.I[fuori], rtol=0, atol=1e-15)
    assert pl.I[1, 1] - con.I[1, 1] == pytest.approx(con.sottratto["J"], rel=1e-10)


@pytest.mark.parametrize("modifica, messaggio", [
    ({"I_kgm2": [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]}, "simmetrico"),
    ({"I_kgm2": [[1.0, 0.0], [0.0, 1.0]]}, "3x3"),
    ({"m_kg": 0.0}, "m_kg"),
    ({"frame": "W"}, "frame"),
    ({"r_G_mm": [1.0, 2.0]}, "r_G_mm"),
    ({"products_sign": "boh"}, "products_sign"),
    ({"source": "onshape"}, "source"),
])
def test_errori_di_validazione(modifica: dict, messaggio: str) -> None:
    dati = config_dati()
    dati["payload"].update(modifica)
    with pytest.raises(pr.ErroreConfigurazione, match=messaggio):
        parametri(dati)


def test_campi_obbligatori_per_terna() -> None:
    dati = config_dati()
    dati["payload"].pop("r_G_mm")
    with pytest.raises(pr.ErroreConfigurazione, match="r_G_mm"):
        parametri(dati)
    dati = config_dati()
    dati["payload"]["frame"] = "body2_o2"
    with pytest.raises(pr.ErroreConfigurazione, match="d_mm"):
        parametri(dati)


def test_tensore_non_fisico() -> None:
    """I_xx + I_yy < I_zz viola la disuguaglianza triangolare."""
    dati = config_dati()
    dati["payload"].update(I_kgm2=[[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 1.0]],
                           subtract_shaft_parts=False)
    with pytest.raises(ValueError, match="triangolare"):
        _payload(dati)


def test_asse_di_puntamento_sospetto() -> None:
    """Se il momento baricentrico minimo non e' attorno a I_hat, avviso di incoerenza."""
    dati = config_dati()
    dati["payload"].update(frame="body2_o2", d_mm=[0.0, 0.0, 0.0], subtract_shaft_parts=False,
                           I_kgm2=[[1.0, 0.0, 0.0], [0.0, 0.6, 0.0], [0.0, 0.0, 1.2]])
    dati["payload"].pop("r_G_mm")
    pl = _payload(dati)
    assert any(a.codice == "PAYLOAD_ASSE_SOSPETTO" for a in pl.avvisi)
    assert not any(a.codice == "PAYLOAD_ASSE_SOSPETTO" for a in _payload(config_dati()).avvisi)


def test_terna_2_K_lungo_il_puntamento() -> None:
    """K_hat = puntamento, I_hat verso l'alto, J_hat = K_hat x I_hat; d_K positivo verso le ottiche."""
    from gimbal_mdof.params import R_P2

    p = parametri(config_dati())
    tr = p.terne
    I0, J0, K0 = tr.R_N2_0.T
    assert np.allclose(K0, tr.R_NP_0[:, 0]) and np.allclose(J0, tr.J_hat_0)
    assert np.allclose(I0, np.cross(J0, K0)) and np.isclose(np.linalg.det(tr.R_N2_0), 1.0)
    assert np.isclose(np.linalg.det(R_P2), 1.0)
    pl = _payload(config_dati())
    # la componente lungo il puntamento e' d_x della terna interna, I_KK il suo momento
    assert pl.d_2[2] == pytest.approx(pl.d_x) and pl.d_2[0] == pytest.approx(pl.d_z)
    assert pl.d_2[1] == pytest.approx(-pl.d_y)
    assert pl.I_2[2, 2] == pytest.approx(pl.I_xx) and pl.I_2[0, 0] == pytest.approx(pl.I_zz)
    assert pl.I_2[0, 2] == pytest.approx(pl.I_xz) and pl.I_2[1, 2] == pytest.approx(-pl.I_xy)
    # I_hat verso l'alto con il telescopio orizzontale
    assert np.allclose(I0, [0.0, 0.0, 1.0])
