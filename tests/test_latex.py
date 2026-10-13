"""Test della sezione 13, punto 22: il documento LaTeX compila e contiene la geometria."""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from gimbal_mdof import analysis as an
from gimbal_mdof import coupled as cp
from gimbal_mdof import fem
from gimbal_mdof import latex_report as lr
from gimbal_mdof import linearize as lz
from gimbal_mdof import params as pr
from gimbal_mdof import payload as plm

from conftest import config_dati, parametri


@pytest.fixture(scope="module")
def tex(par: pr.Parametri, tmp_path_factory: pytest.TempPathFactory) -> Path:
    assi = (fem.costruisci_asse(par.azimut), fem.costruisci_asse(par.polare))
    pl = plm.costruisci_payload(par, assi[1])
    combinazioni = [("consistent", "guyan"), ("consistent", "guyan_lumped"), ("lumped", "guyan")]
    sistemi = {c: cp.costruisci_sistema(par, *c, payload=pl, assi=assi) for c in combinazioni}
    w = 2 * np.pi * np.logspace(-1, 5.5, 200)
    confronti = [an.confronta_riduzione(sistemi[c], a, 0.0, w) for c in combinazioni
                 for a in ("azimut", "polare")]
    ml = lz.linearizza(sistemi[("consistent", "guyan")].accoppiato, par.theta0[0], 0.0)
    cont = lr.ContenutoLatex(par=par, sistemi=sistemi, confronti=confronti, modelli_lineari=[ml],
                             figure={}, avvisi=list(par.avvisi) + list(pl.avvisi), ml_tabella=ml)
    return lr.genera(cont, tmp_path_factory.mktemp("latex"))


def test_contenuto(tex: Path) -> None:
    testo = tex.read_text(encoding="utf-8")
    # una figura di geometria per asse, piu' quella dei sistemi di riferimento
    assert testo.count(r"\begin{tikzpicture}") == 3
    assert r"\label{fig:geom_azimut}" in testo and r"\label{fig:geom_polare}" in testo
    assert r"\label{sec:frames}" in testo and r"\label{fig:frames}" in testo
    assert r"\label{subsec:input_frames}" in testo
    assert r"\setcounter{MaxMatrixCols}{20}" in testo
    assert r"\begin{align}" not in testo
    assert "placeholders" in testo


def test_mappatura_documento_coerente_con_il_codice(par: pr.Parametri) -> None:
    """Le formule I_N -> I_2 e r -> d_2 stampate nel PDF riproducono il tensore del codice."""
    from gimbal import inertia as inz

    dati = config_dati()
    dati["payload"]["subtract_shaft_parts"] = False
    p = parametri(dati)
    pl = plm.costruisci_payload(p, fem.costruisci_asse(p.polare))
    tr = p.terne
    I_G = inz.normalizza_prodotti(p.payload.I, p.payload.segno_prodotti)
    r = p.payload.r_G_N - tr.r_o2
    I_O_W = inz.trasporta(I_G, p.payload.m, r)
    IB, dB = lr.mappatura_N_2(tr.R_N2_0)
    sost = dict(zip(lr.SIMBOLI_I_N, [I_O_W[0, 0], I_O_W[1, 1], I_O_W[2, 2], I_O_W[0, 1],
                                     I_O_W[0, 2], I_O_W[1, 2]]))
    sost.update(dict(zip(lr.SIMBOLI_R, r)))
    I_doc = np.array(IB.subs(sost), dtype=float)
    d_doc = np.array(dB.subs(sost), dtype=float).reshape(3)
    assert np.allclose(I_doc, pl.I_2, rtol=0, atol=1e-13)
    assert np.allclose(d_doc, pl.d_2, rtol=0, atol=1e-15)


def test_verso_numerazione_polare_speculare() -> None:
    """Numerazione lungo -J_hat con y_1 -> -y_1: configurazione speculare, stesso J^ax."""
    base = config_dati()
    base["polar"]["y_node1_from_o2_mm"] = -75.0
    speculare = config_dati()
    speculare["polar"].update(y_node1_from_o2_mm=75.0, nodes_increase_along="-J_hat")
    p1, p2 = parametri(base), parametri(speculare)
    J1 = plm.J_ax(fem.costruisci_asse(p1.polare), p1.polare.y_nodo1_da_O, p1.polare.verso_numerazione)
    J2 = plm.J_ax(fem.costruisci_asse(p2.polare), p2.polare.y_nodo1_da_O, p2.polare.verso_numerazione)
    assert J1 == pytest.approx(J2, rel=1e-14)
    senza = plm.J_ax(fem.costruisci_asse(p1.polare), None)
    assert J1 > senza  # il termine di Steiner c'e'
    errato = config_dati()
    errato["polar"]["nodes_increase_along"] = "+x"
    with pytest.raises(pr.ErroreConfigurazione, match="nodes_increase_along"):
        parametri(errato)


@pytest.mark.lento
@pytest.mark.skipif(not (shutil.which("pdflatex") or shutil.which("latexmk")),
                    reason="LaTeX non installato")
def test_compilazione(tex: Path) -> None:
    esito = lr.compila(tex)
    assert esito.compilato, esito.errori
    assert not esito.riferimenti_irrisolti
    log = tex.with_suffix(".log").read_text(encoding="utf-8", errors="replace")
    assert "There were undefined references" not in log
    assert "Citation" not in log or "undefined" not in log
