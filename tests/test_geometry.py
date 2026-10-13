"""Test della sezione 13, punti 1-3: proprieta' degli elementi e validazione."""

from __future__ import annotations

import numpy as np
import pytest

from gimbal_mdof import geometry as geo
from gimbal_mdof import params as pr

from conftest import config_dati, parametri

ACCIAIO = pr.Materiale(nome="acciaio", rho=8000.0, G=77e9, segnaposto=False)
ALLUMINIO = pr.Materiale(nome="alluminio", rho=2810.0, G=26.9e9, segnaposto=False)


def _elemento(segmenti: list[pr.Segmento], *, k_joint: float | None = None) -> pr.Elemento:
    return pr.Elemento(
        id=1, modo="geometry", lunghezza=sum(s.lunghezza for s in segmenti),
        segmenti=tuple(segmenti), k_joint=k_joint, extra=(), J_manuale=None,
        k_manuale=None, m_manuale=None, I_t_manuale=None, y_c_manuale=None,
        segnaposto=False,
    )


def _strato(d_out: float, d_in: float = 0.0, *, mat: pr.Materiale = ACCIAIO,
            portante: bool = True) -> pr.Strato:
    return pr.Strato(forma="hollow" if d_in > 0 else "solid", d_out=d_out, d_in=d_in,
                     materiale=mat, portante=portante)


# --- test 1 -----------------------------------------------------------------


def test_cavo_tende_al_pieno() -> None:
    L, D = 0.04, 0.03
    pieno = geo.proprieta_strato(_strato(D), L)
    cavo = geo.proprieta_strato(_strato(D, 1e-9), L)
    for nome in ("m", "J", "I_t", "J_p", "GJ_p"):
        assert getattr(cavo, nome) == pytest.approx(getattr(pieno, nome), rel=1e-12)


def test_additivita_cavo_piu_pieno() -> None:
    """Elemento = cavo di lunghezza a + pieno di lunghezza b: somma dei singoli."""
    a, b = 0.025, 0.040
    cavo = pr.Segmento(lunghezza=a, strati=(_strato(0.030, 0.010),))
    pieno = pr.Segmento(lunghezza=b, strati=(_strato(0.020),))
    composto = geo.proprieta_elemento(_elemento([cavo, pieno]))
    solo_cavo = geo.proprieta_elemento(_elemento([cavo]))
    solo_pieno = geo.proprieta_elemento(_elemento([pieno]))
    assert composto.J == pytest.approx(solo_cavo.J + solo_pieno.J, rel=1e-14)
    assert composto.m == pytest.approx(solo_cavo.m + solo_pieno.m, rel=1e-14)
    assert 1.0 / composto.k == pytest.approx(1.0 / solo_cavo.k + 1.0 / solo_pieno.k, rel=1e-14)
    # inerzia trasversale nel baricentro composto: Steiner lungo l'asse
    x_G = composto.x_G
    I_t_atteso = (solo_cavo.I_t + solo_cavo.m * (a / 2 - x_G) ** 2
                  + solo_pieno.I_t + solo_pieno.m * (a + b / 2 - x_G) ** 2)
    assert composto.I_t == pytest.approx(I_t_atteso, rel=1e-12)


# --- test 2 -----------------------------------------------------------------


def test_valori_analitici_cilindro() -> None:
    rho, G = 8000.0, 77e9
    L, Do, Di = 0.05, 0.04, 0.02
    s = geo.proprieta_strato(_strato(Do, Di), L)
    m = rho * np.pi / 4 * (Do**2 - Di**2) * L
    assert s.m == pytest.approx(m, rel=1e-14)
    assert s.J == pytest.approx(np.pi * rho * L / 32 * (Do**4 - Di**4), rel=1e-14)
    assert s.J == pytest.approx(m / 8 * (Do**2 + Di**2), rel=1e-14)
    assert s.I_t == pytest.approx(m * ((Do**2 + Di**2) / 16 + L**2 / 12), rel=1e-14)
    k = G * np.pi / 32 * (Do**4 - Di**4) / L
    el = geo.proprieta_elemento(_elemento([pr.Segmento(lunghezza=L, strati=(_strato(Do, Di),))]))
    assert el.k == pytest.approx(k, rel=1e-14)


def test_serie_e_parallelo() -> None:
    seg = pr.Segmento(lunghezza=0.03, strati=(_strato(0.02),))
    k1 = geo.proprieta_elemento(_elemento([seg])).k
    assert geo.proprieta_elemento(_elemento([seg, seg])).k == pytest.approx(k1 / 2, rel=1e-14)
    # due strati portanti in parallelo: le rigidezze si sommano
    anima = _strato(0.016)
    mozzo = _strato(0.060, 0.016, mat=ALLUMINIO)
    k_anima = geo.proprieta_elemento(_elemento([pr.Segmento(0.017, (anima,))])).k
    k_mozzo = geo.proprieta_elemento(_elemento([pr.Segmento(0.017, (mozzo,))])).k
    k_par = geo.proprieta_elemento(_elemento([pr.Segmento(0.017, (anima, mozzo))])).k
    assert k_par == pytest.approx(k_anima + k_mozzo, rel=1e-14)
    # uno strato non portante non cambia la rigidezza ma aggiunge inerzia
    distanziale = _strato(0.030, 0.016, portante=False)
    con = geo.proprieta_elemento(_elemento([pr.Segmento(0.017, (anima, distanziale))]))
    assert con.k == pytest.approx(k_anima, rel=1e-14)
    assert con.J > geo.proprieta_elemento(_elemento([pr.Segmento(0.017, (anima,))])).J
    # cedevolezza di giunzione in serie
    kj = 1e4
    giunto = geo.proprieta_elemento(_elemento([seg], k_joint=kj)).k
    assert 1 / giunto == pytest.approx(1 / k1 + 1 / kj, rel=1e-14)


def test_taglio_conserva_massa_e_inerzia() -> None:
    segs = [pr.Segmento(0.025, (_strato(0.025),)), pr.Segmento(0.035, (_strato(0.025, 0.01),))]
    intero = geo.proprieta_elemento(_elemento(segs))
    tagliati = geo.taglia_segmenti(segs, [0.03])
    assert len(tagliati) == 3
    diviso = geo.proprieta_elemento(_elemento(tagliati))
    assert diviso.m == pytest.approx(intero.m, rel=1e-14)
    assert diviso.J == pytest.approx(intero.J, rel=1e-14)
    assert diviso.k == pytest.approx(intero.k, rel=1e-14)
    assert diviso.I_t == pytest.approx(intero.I_t, rel=1e-12)


# --- test 3 -----------------------------------------------------------------


def _con_strato(strato: dict) -> dict:
    dati = config_dati()
    dati["azimuth"]["elements"][0]["segments"][0]["layers"] = [strato]
    return dati


def test_errore_pieno_con_foro() -> None:
    with pytest.raises(pr.ErroreConfigurazione, match="pieno"):
        parametri(_con_strato({"shape": "solid", "d_in_mm": 5.0, "d_out_mm": 12.0,
                               "material": "acciaio", "carries_torque": True}))


@pytest.mark.parametrize("d_in", [12.0, 15.0, 0.0])
def test_errore_cavo_diametri(d_in: float) -> None:
    with pytest.raises(pr.ErroreConfigurazione, match="0 < d_in < d_out"):
        parametri(_con_strato({"shape": "hollow", "d_in_mm": d_in, "d_out_mm": 12.0,
                               "material": "acciaio", "carries_torque": True}))


def test_errore_segmento_senza_strati_portanti() -> None:
    with pytest.raises(pr.ErroreConfigurazione, match="carries_torque"):
        parametri(_con_strato({"shape": "solid", "d_out_mm": 12.0, "material": "acciaio",
                               "carries_torque": False}))


def test_altri_errori_di_validazione() -> None:
    dati = config_dati()
    dati["polar"]["elements"] = dati["polar"]["elements"][:-1]
    with pytest.raises(pr.ErroreConfigurazione, match="n_nodes - 1"):
        parametri(dati)
    dati = config_dati()
    dati["polar"]["payload_xi"] = 1.5
    with pytest.raises(pr.ErroreConfigurazione, match="payload_xi"):
        parametri(dati)
    dati = config_dati()
    dati["polar"]["bearing_elements"] = [2, 3]
    with pytest.raises(pr.ErroreConfigurazione, match="cuscinetto"):
        parametri(dati)
    dati = config_dati()
    dati["azimuth"]["node_inertias"][-1]["J_kgm2"] = None
    with pytest.raises(pr.ErroreConfigurazione, match="Q9"):
        parametri(dati)
    dati = config_dati()
    dati["azimuth"]["masters"] = [1, 9]
    with pytest.raises(AssertionError, match="master"):
        parametri(dati)


def test_config_progetto_valida(par: pr.Parametri) -> None:
    assert par.azimut.n_elementi == 8 and par.polare.n_elementi == 10
    assert par.segnaposto  # la geometria e' segnaposto: le uscite vanno marcate
    assert par.polare.xi_payload == 0.5
    codici = {a.codice for a in par.avvisi}
    assert {"GEOMETRIA_AZIMUT", "GEOMETRIA_POLARE", "Q7_POSIZIONE_POLARE"} <= codici
