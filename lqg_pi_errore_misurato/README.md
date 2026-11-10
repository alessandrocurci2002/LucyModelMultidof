# LQG-PI con il PI sull'errore misurato (Gawronski 2008, Fig. 8.4b)

Richiesta dell'utente (2026-11-05): tenere l'implementazione della Fig. 8.4b, con la parte proporzionale-integrale sull'errore
misurato e = r - y (baseline robusta, che non dipende dal modello interno) e la retroazione degli stati flessibili stimati.
L'LQG-PI della tesi usa invece l'errore stimato, e = r - y_hat (legge 8.15, Fig. 8.3). La tesi non cambia.

## File

| File | Contenuto |
|---|---|
| `controllore.py` | `ControlloreLQGPIMisurato`: stessa sintesi (CARE, FARE) di `lucy_cl.lqg.sintesi`, PI su r - y; funzioni di trasferimento e poli in anello chiuso; `attiva()` sostituisce sintesi e poli solo nel processo che la chiama |
| `taratura.py` | ritaratura con `analisi_lqg_passi_2_4` (stessi passi e vincoli della tesi) piu' il vincolo di bassa autorita' q_p / q_f >= 10 |
| `simulazioni.py` | le simulazioni di `analisi_tesi_hil.py` (circonferenza, dataset7/8 con config_v3 ed EKF2, scenari 0-3, tre semi) con questo controllore |
| `funzioni_trasferimento.py` | T_ry, T_wy, T_ra, T_wa di cascata, LQG-PI della tesi e LQG-PI della Fig. 8.4b |
| `confronto.py` | tabella di confronto con i risultati della tesi |
| `config/lqg_tuning.json` | tuning nuovo (copia di `results/steps_2_4/summary.json`) |
| `config/lqg_tuning_estimated_error.json` | tuning della tesi, per riferimento |
| `results/`, `figures/` | risultati e figure (asse delle frequenze senza numeri) |

## Modifica al codice condiviso

`g.corri_lqg`, `cc.simula_coppia` e `p1.Banco.gradino` passano y(k) al controllore con `misura` prima di `uscita`, solo se il
controllore ha quel metodo; i wrapper `Registra` e `ConRumore` lo inoltrano. Il controllore della tesi non lo ha: i suoi
risultati non cambiano (verificato: la cascata e i casi della tesi si riproducono identici).

## Esecuzione

    .venv/bin/python -I lqg_pi_errore_misurato/taratura.py
    LAVORI=2 .venv/bin/python -I lqg_pi_errore_misurato/simulazioni.py
    .venv/bin/python -I lqg_pi_errore_misurato/funzioni_trasferimento.py
    .venv/bin/python -I lqg_pi_errore_misurato/confronto.py

Con 4 processi le simulazioni hanno esaurito la memoria della macchina (15 GB, con Firefox e VS Code aperti): usare 2.
