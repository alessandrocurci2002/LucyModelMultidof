"""Modello multi-GdL di un telescopio alt-azimutale: FEM torsionale degli alberi,
riduzione di Guyan e accoppiamento dei due assi ridotti.

Il pacchetto implementa, modulo per modulo, la specifica
``SPEC_modello_gimbal_multidof.md``. Il flusso rispecchia i quattro stadi della
sezione 0 della specifica:

======================  =====================================================
modulo                  contenuto
======================  =====================================================
:mod:`params`           lettura e validazione di ``config.yaml`` (sez. 5)
:mod:`geometry`         proprieta' degli elementi da cilindri (sez. 6)
:mod:`fem`              modello FEM torsionale di ciascun asse (sez. 7)
:mod:`reduction`        condensazione di Guyan, Craig-Bampton (sez. 8)
:mod:`payload`          adattatore verso ``gimbal_2dof`` (telescopio)
:mod:`coupled`          modello accoppiato, Lagrange con SymPy (sez. 9)
:mod:`linearize`        punto di lavoro e matrici linearizzate (sez. 10.1-2)
:mod:`tf`               FRF, poli, zeri invarianti, esportazione (sez. 10.3)
:mod:`analysis`         tabelle e analisi parametriche (sez. 11)
:mod:`simulate`         scenari S1..S3 (sez. 11)
:mod:`plotting`         figure (sez. 11)
:mod:`latex_report`     documento LaTeX con la geometria in TikZ (sez. 12)
:mod:`report`           ``results/report.md`` (sez. 14)
======================  =====================================================

Convenzioni vincolanti (sezione 2 della specifica): SI internamente, angoli in
rad, nodi numerati da 1 in tutti gli input e output (base 0 solo nel codice),
tensori d'inerzia in convenzione tensoriale.
"""

from __future__ import annotations

__version__ = "0.1.0"
