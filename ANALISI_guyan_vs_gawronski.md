# Riduzione di Guyan contro troncamento modale alla Gawronski

> **Realismo.** I numeri vengono da `config.yaml` del 13/10/2026, in parte segnaposto (forcella 1 kg·m², smorzamento, materiali). Il payload è stato **stimato a mano** da config.yaml (pacchetto `gimbal` del 2 GdL non disponibile in sandbox): vedi §6. Analisi per asse isolato, θ0 = 0, Ω0 = 0, k_g = 0, massa consistente. Script: `confronto_guyan_modale.py`.

## 1. Il punto di partenza

Nel modello ridotto di Guyan i **due modi rigidi esistono già e sono esatti**: da K·1 = 0 segue T·1 = 1, quindi il vettore rigido sta nello span della base di Guyan (frequenza 0 a precisione di macchina, un modo per asse). Nello spettro dei due assi insieme l'ordine è già: azimut rigido, polare rigido, poi i flessibili. Quello che distingue i due approcci non è l'ordine dei modi ma **che cosa sono gli altri modi** e che cosa viene conservato esattamente.

## 2. Le due basi di Ritz

| | Guyan | Gawronski (troncamento modale) |
|---|---|---|
| Base | spostamenti statici dei GdL slave per spostamento unitario di un master: `T = [1; -Kss⁻¹ Ksm]` | autovettori del FEM completo: `Φ = [1, φ1, ..., φk]` |
| Coordinate | **fisiche**: rotazioni dei nodi master (ψe, ψm, ψL, ...) | **modali**: ampiezze η_i; il modo rigido η0 è l'angolo medio pesato sull'inerzia |
| M, K ridotte | piene (tridiagonali per catena); K_r = molle in serie reali dell'albero | diagonali; accoppiamento solo tramite B_m, C_m e i vettori d'aggancio c_L, c_P |
| Esatto | statica fra i master (cedevolezza, zeri legati alla statica), modo rigido | frequenze, forme, masse modali dei k modi tenuti |
| Approssimato | dinamica: l'inerzia degli slave è ignorata (frequenze **per eccesso**, errore crescente con l'ordine) | statica: la flessibilità dei modi scartati manca (zeri, cedevolezze) |
| N. modi flessibili | = n. master − 1, fissato dalla scelta dei nodi | scelto con le norme di Hankel (Gawronski 4.4–4.6) |
| Correzione | Craig–Bampton (modi a interfaccia fissa) | correzione statica (flessibilità residua) |

Se si vogliono coordinate fisiche **e** modi esatti si può usare SEREP: stesse FRF del troncamento modale (verificato, scarto 1e-8), ma K_r perde la forma di catena e le molle non sono più quelle dell'albero (azimut: k12 = 2.16e4 contro 9.55e3 N·m/rad di Guyan; polare k(7,11) = 1.27e4 contro 9.69e3).

## 3. Risultati sul modello (consistente, θ0 = 0)

**Frequenze (Hz).**

| asse | modo | FEM completo | Guyan | errore |
|---|---|---|---|---|
| azimut | 1 | 752.2 | 834.3 | +10.9 % |
| azimut | 2 | 4595.1 | 6167.6 | +34.2 % |
| polare | 1 | 335.6 | 354.4 | +5.6 % |
| polare | 2 | 3557.2 | 3857.3 | +8.4 % |
| polare | 3 | 14568.1 | 14583.0 | +0.1 % |

Il troncamento modale ha errore nullo sui modi tenuti. Guyan rappresenta 5 modi flessibili sui 17 del FEM (azimut 2 su 8, polare 3 su 9). Il MAC di Guyan è alto (0.997–0.9998; 0.977 sul secondo modo azimutale) ma **il MAC non rivela l'errore in frequenza**: la forma è quasi giusta, la frequenza no.

**Zeri non smorzati (antirisonanze), Hz.**

| TF | completo | Guyan | modale puro | modale + correzione |
|---|---|---|---|---|
| azimut τ→encoder #1 | 13.3014 | 13.3014 | 13.3017 (+0.002 %) | 13.3014 |
| polare τ→motore #1 | 17.1788 | 17.1788 | 17.1847 (+0.035 %) | 17.1788 |
| polare τ→encoder #1 | 504.49 | 577.35 (+14.4 %) | 504.49 | 504.49 |
| azimut τ→motore #2 | 1257.06 | 1886.0 (+50 %) | 1258.0 (+0.08 %) | 1257.06 |

Lo zero a bassa frequenza (13.3 Hz azimut, 17.2 Hz polare) è l'inerzia del carico sulla molla fra motore e carico: dipende solo da statica e modo rigido, quindi Guyan lo ritrova **esatto**, mentre il troncamento modale puro lo sposta di 0.03–0.05 % (flessibilità residua mancante, −0.09 % e −0.07 % sulla F collocata). Gli zeri successivi dipendono dalla dinamica degli slave e Guyan li sbaglia del 14–50 %.

**Banda di validità** (primo scarto > 3 dB sulla FRF τ→encoder, criterio di `analysis.py`): Guyan 659 Hz (azimut) e 311 Hz (polare), cioè poco sotto il primo modo; troncamento con correzione statica oltre 29 kHz. Oltre i primi modi il confronto è con un FEM a 8–9 elementi che ha errore di discretizzazione proprio (modo 3 azimutale 18.7 kHz con massa concentrata contro 21.3 kHz con consistente): non leggere le bande alte come caratterizzazione fisica. Vicino a un'antirisonanza stretta il rapporto fra FRF è ipersensibile: la banda "severa" (1 dB / 5°) del modale puro collassa a 13–17 Hz solo per lo spostamento dello zero, per questo conviene guardare la tabella degli zeri.

**Norme di Hankel (τ→encoder).** Un solo modo flessibile conta per asse: azimut 1.96e-3 contro 9.1e-6 del secondo; polare 3.6e-6 contro 2.7e-9. I modi tenuti da Guyan sono quelli giusti: il problema non è la selezione ma il bias di frequenza proprio sul modo dominante.

## 4. Perché Guyan sbaglia: scelta dei master

L'errore dipende dall'inerzia degli slave **rispetto all'inerzia che partecipa al modo**, non rispetto all'inerzia totale. I master ([1, 5, 9] azimut; [1, P, 7, 11] polare) portano le inerzie grandi (forcella, telescopio), ma i modi flessibili sono formati dalle parti leggere, che sono slave: il disco in alluminio dell'elemento 2 (azimut, J = 3.5e-4) e il disco in acciaio dell'elemento 9 (polare, J = 1.26e-3).

Aggiungere **un solo master** nel nodo del disco:

| asse | master | modo 1 | modo 2 | modo 3 |
|---|---|---|---|---|
| azimut | [1,5,9] | +10.91 % | +34.2 % | – |
| azimut | + nodo 3 | +0.01 % | +0.15 % | +38 % (il 3° è oltre la base) |
| polare | [1,P,7,11] | +5.6 % | +8.4 % | +0.1 % |
| polare | + nodo 9 | 0.0 % | +0.08 % | +0.1 % |

Craig–Bampton già nel codice (`riduci_craig_bampton`) con n = 1–2 dà errori ≤ 0.12 % sui primi modi.

## 5. Raccomandazioni

1. **Coerenza del modello unificato.** `MODELLO_unificato_gawronski_dedonatis.md` §3–4 costruisce Φ dalle matrici **già ridotte con Guyan** (3 e 4 GdL). I modi risultanti sono quelli di Guyan, con le frequenze per eccesso di §3, e non i modi del FEM. Per un'identificazione e un confronto con prove modali serve Φ dal FEM completo (9 e 10 GdL con P), poi: modo rigido + modi scelti con Hankel, correzione statica per gli zeri, righe c_L e c_P di Φ nei nodi 9 e P. Nel documento unificato, le tabelle di §8.2 (bande di validità) valgono quindi per il Guyan-modale e non per il troncamento del FEM.
2. **Se si resta con Guyan** (coordinate fisiche, catena di molle): aggiungere i master nei nodi 3 (azimut) e 9 (polare), oppure Craig–Bampton con 1–2 modi. Si perde la topologia a 3 e 4 GdL della spec §0, si guadagna l'accuratezza.
3. **Correzione statica**: migliora gli zeri a bassa frequenza (qui 0.03–0.05 %) ma lascia un termine costante che sopra l'ultimo modo scartato non decade come il modello vero (verde sopra il nero in figura oltre 20 kHz): usarla per il progetto del controllo, non per bande molto alte.
4. **`mass_form: guyan_lumped`** è accettabile sul primo modo (+7.9 %, +3.2 %) ma sbaglia pesantemente i successivi (azimut modo 2 −60 %, polare modo 2 −64 %, modo 3 −21 %): la somma di riga di M_r sposta troppa massa sui master.
5. **Dati di ingresso da controllare** (mie osservazioni, da confermare): `polar.node_inertias` è tutto commentato (rotore motore, disco encoder) e nell'azimut non c'è inerzia di rotore al nodo 5 né di encoder al nodo 1, quindi i nodi motore ed encoder sono quasi privi di massa. Con rotore 8.19e-5 e encoder 2.0e-5 sul polare (valori delle righe commentate) il modo 2 scende da 3557 a 2163 Hz e il 3 da 14568 a 7005 Hz, mentre l'errore di Guyan sul modo 1 resta ~5.3 %.

## 6. Limiti dell'analisi

- Asse isolato, Ω0 = 0, payload bilanciato (k_g = 0): niente accoppiamento non lineare ψ–θ, niente pendolo. Poiché l'accoppiamento passa dalle righe di Φ nei nodi L e P, ci si aspetta che le conclusioni sul bias restino, ma non è verificato sul modello a 7 GdL.
- Payload stimato a mano da config.yaml (I_yy,p = 0.8815 kg·m² per Steiner dal tensore Onshape; J_ψ,p(θ0) per rotazione del tensore; J^ax = Σ I_t dei pezzi polari). Con il pacchetto `gimbal` installato lo script usa `costruisci_sistema` e non serve `--payload-manuale`.
- La frequenza del primo modo azimutale è indipendente da θ0 nel campo ±45° (752.1–752.2 Hz) perché la forcella domina l'inerzia; cambia la posizione dell'antirisonanza (√(k/J)).
- Smorzamento C = βK con ζ = 0.01 sul primo modo (segnaposto): influenza ampiezza dei picchi, non frequenze né zeri non smorzati.
