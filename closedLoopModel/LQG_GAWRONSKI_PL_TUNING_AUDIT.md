# Audit del capitolo 8 di Gawronski: taratura di un controllore PL

PL = PI nella maglia di velocità, LQG nella maglia di posizione. Integra `LQG_GAWRONSKI_AUDIT.md` e confronta il libro con la taratura di `../analisi_lqg_gawronski.py`.

Fonte: W. Gawronski, *Modeling and Control of Antennas and Telescopes*, Springer, 2008. I numeri di pagina sono quelli stampati; per il capitolo 8 nel PDF vanno aggiunti 12.

Struttura:
- §1 conclusioni;
- §2 PL e DSN nel libro;
- §3 procedura di taratura del libro;
- §4 prestazioni misurate e limiti;
- §5 incoerenze del testo;
- §6 confronto con la taratura fatta su LUCY;
- §7 procedura proposta per LUCY;
- §8 audit puntuale del §8.2 e verifica su LUCY (passo 1).

---

## 1. Conclusioni

1. **La taratura del libro è guidata da inseguimento e disturbi, non da mosse veloci.** La misura primaria è l'errore rms nel vento durante l'inseguimento (§8.3, p. 110). La mia taratura usava solo i doppietti, quindi ottimizzava una grandezza che il libro considera secondaria.
2. **Devo correggere una mia affermazione: PL non «dà poco» nel libro.** Il fattore 2 è una simulazione sull'LMT (§8.5.4). Sull'antenna reale DSS55 da 34 m il rapporto PI/LQG misurato nel vento è **4.9 in azimut e 6.5 in elevazione** (Tab. 8.1, p. 117).
3. **Feedforward, preprocessore del comando (CPP) e limiti non risultano nell'LQG della DSN descritto dal capitolo 8.** Sono nell'analisi LMT (§8.5.4, pp. 129–131). Il feedforward sulla 34 m compare solo con il PI, come simulazione (§7.5). Il mio schema senza di essi è quindi coerente con §8.1–8.4, non con l'analisi LMT.
4. **Sulla PI di velocità della DSN il testo non è univoco** (§2): §8.5 dice «simple proportional controller», §8.5.2 elenca PL come PI + LQG e dice che è la configurazione DSN. I guadagni della maglia di velocità sono dati solo per l'LMT.
5. **La relazione pesi → guadagni che ho usato è confermata** dai valori della GUI del libro (§3.3).

---

## 2. PL e DSN nel libro

| Elemento | DSN (§1.2, §8.1–8.4) | LMT, analisi PL (§8.5.4) |
|---|---|---|
| Stato | implementato; misure a DSS55 | simulazione su modello analitico |
| Velocity loop | «velocity controller» con tachimetro (Fig. 1.6–1.7); «simple proportional» (§8.5, p. 124) | PI: $k_p=300$; $k_i=300$ az, $400$ el; 1 kHz; banda 1.0 Hz |
| Plant della sintesi | velocity loop: ingresso $u_c$ [deg/s], uscita encoder [deg] | velocity loop del modello a elementi finiti |
| Identificazione | n4sid dalla prova (§3.2, vedi audit precedente) | modello analitico |
| Position loop | LQG nel calcolatore | LQG + feedforward + CPP + limiti |
| Feedforward e CPP con LQG | non nominati in Fig. 8.1–8.4, §8.2, §8.3.2, §8.3.4, §8.4 | sì; CPP: $k_v=6$, $k_o=0.6$, $\beta=20$; 1 deg/s, 0.5 deg/s² |

Il tachimetro è una differenza rispetto a LUCY, dove velocità e posizione vengono dallo stesso encoder.

Il capitolo 8 non cita un test con PI di velocità sulla DSN e non dà i guadagni di quel loop. **«PI come implementato per la DSN» non si può quindi dedurre con certezza dal libro:** la struttura PI è verificabile solo sull'analisi LMT.

---

## 3. La procedura di taratura

### 3.1 Principi (§8.1.4, 8.1.6, 8.4.1)

- Coordinate modali: stati debolmente accoppiati, quindi $Q$ e $V$ diagonali.
- $Q=V=\mathrm{diag}(q_{ei},q_e,q_{f1},q_{f1},\dots,q_{fn},q_{fn})$ (8.26), $R=W=1$. Con 4 modi sono 10 pesi, ma ogni coppia ha un solo valore: 6 parametri.
- $Q=V$ è dichiarata «opportunistic»; «the tuning of Q and V separately may improve the performance» (p. 119).
- Bassa autorità: $\|Q_t\|\gg\|Q_f\|$ (8.22).

| Peso (Tab. 8.2) | Agisce su |
|---|---|
| Proporzionale | assestamento, banda, reiezione |
| Integrale | sovraelongazione, errore a velocità costante, reiezione |
| Frequenza $i$ | ampiezza di vibrazione del modo $i$ |

### 3.2 Passi (§8.2, pp. 108–110)

1. Pesi $q_i$, $q_p$ e flessibili piccoli. Se $T_{ry}$ mostra risonanze eccessive, alzare **solo** il peso del modo corrispondente.
2. Alzare $q_p$ finché la banda raggiunge la frequenza fondamentale dell'antenna.
3. Alzare $q_i$ finché compaiono oscillazioni, rispettando il guadagno integrale critico.
4. Se le risonanze riemergono, correggere i pesi flessibili.

**I pesi flessibili si fissano all'inizio e si ritoccano alla fine.** Non vanno lasciati solo in coda.

### 3.3 Esempio numerico: 34 m (§8.4.2, pp. 120–122)

| Passo | $q_i$ | $q_p$ | $q_{f1}$ | $t_s$ [s] | Sovr. [%] | Banda [Hz] | $k_p$ | $k_i$ |
|---|---|---|---|---|---|---|---|---|
| partenza | 0.1 | 0.3 | 0 | 8.52 | 18 | 0.21 | – | – |
| modo 1 (Fig. 8.17a) | 0.1 | 0.3 | 5 | 8.66 | 17 | 0.21 | 0.9485 | 0.304 |
| finale (Fig. 8.17b) | 3 | 10.05 | 5 | 4.06 | 13 | 0.79 | 3.549 | 1.61 |

- **Verifica della formula rigida** $k_i=\sqrt{q_i}$, $k_p^2=q_p+2k_i$: dà $k_i=0.316$ e $k_p=0.966$ contro 0.304 e 0.9485 della GUI, e $k_i=1.73$, $k_p=3.68$ contro 1.61 e 3.549 all'ultimo passo.
- **Rapporto tipico:** $k_i=0.5\,k_p$ «typical for the DSN antennas» (p. 111). Le figure coprono $k_i/k_p$ da 0.2 a 0.8.
- **Prova usata:** gradino di 1°, più $T_{ry}$ e $T_{wy}$ (§8.4.2, punti 6–7).

### 3.4 Taratura fine (§8.4.3)

Ottimizzazione vincolata dei pesi (`constr`), partendo dai pesi già tarati a mano:

$$f=w_1e_r+w_2e_{d\max}+w_3e_w+w_4e_o+w_5t_s+w_6(3-f_o)+w_7m_h+w_8m_{d\max}\quad(8.27)$$

| Variabile | Definizione |
|---|---|
| $e_r$ | errore in regime a 0.1 deg/s |
| $e_{d\max}$ | errore massimo per disturbo a gradino unitario |
| $e_w$ | errore rms nel vento a 24 km/h, 100 s |
| $e_o$, $t_s$ | sovraelongazione e assestamento del gradino unitario |
| $f_o$ | banda; il 3 è il limite della 34 m, va adattato |
| $m_h$ | picco di $T_{ry}$ sopra 3.5 Hz |
| $m_{d\max}$ | picco di $T_{wy}$ |

Ogni variabile è ottimizzata ($w_i>0$) o vincolata. Si parte da pesi nulli per verificare i vincoli, si danno priorità e si procede a blocchi di circa 300 iterazioni. Il capitolo fornisce variabili e metodo, non i valori di $w_i$ usati.

---

## 4. Prestazioni e limiti

**Relazioni della 34 m (§8.3.2, Fig. 8.11–8.13, elevazione).** Con $k_p$ da 0.5 a 3: assestamento da 14 a 3 s, banda da 0.1 a 0.7 Hz, errore nel vento (20 mph) da 5 a 1 mdeg. Più banda e meno assestamento significano meno errore nel vento.

**Reiezione e banda (§8.3.3).** $T_{wy}=G(1-T_{ry})$: il miglioramento viene dall'allargamento della banda, ma $T_{ry}$ deve essere stretta a 1 solo in bassa frequenza. La frequenza fondamentale limita la reiezione.

**Misura reale (§8.3.4, Tab. 8.1).** DSS55, 16.1 km/h:

| Controllore | AZ rms [mdeg] | EL rms [mdeg] | Radiale [mdeg] |
|---|---|---|---|
| PI | 0.32 | 2.69 | 2.71 |
| LQG | 0.066 | 0.42 | 0.42 |
| PI / LQG | 4.9 | 6.5 | 6.5 |

Le misure di PI e LQG sono in sequenza e scalate col quadrato della velocità del vento. Le pagine lette non riportano $k_p$, $k_i$ di questa prova.

**Limite principale (§8.3.5).** «The acceleration limit effectively cancels the benefits of a high-gain controller». Un LQG ad alto guadagno, su un gradino di 20 mdeg, chiederebbe 10 deg/s² contro un limite di 0.4 deg/s²: entra in regime non lineare e oscilla (Fig. 8.15).

**Limiti di §8.1.7.** Pesi flessibili eccessivi riducono banda e reiezione; il peso di posizione è limitato dall'accelerazione; il peso integrale è limitato dal guadagno critico.

**Feedforward e CPP (§7.5, §8.5.3–8.5.4).** Il feedforward deriva il comando e lo manda all'ingresso del velocity loop; migliora l'inseguimento ad alta velocità, non la reiezione dei disturbi. Sulla 34 m con PI porta l'assestamento da 15 a 4 s, ma con una risposta più oscillante, che il CPP elimina. Il CPP serve «to avoid large overshoots during target acquisition and limit cycling during slewing». Il §11.1.1 non fa parte di questo audit.

---

## 5. Incoerenze del testo

| Punto | Discrepanza | Cosa uso |
|---|---|---|
| Passo 3, p. 110 | cita «condition (8.15)», che è la legge di controllo; il guadagno critico è la (7.9) | (7.9) |
| Segni della legge | Fig. 8.3–8.4: $K_f$ col più; 8.15: tutti col meno; 8.17 e Fig. 8.22c: $k_p,k_i$ col più, $k_f$ col meno | 8.17 e Fig. 8.22c. La 8.15 con $e=r-\hat y$ dà un anello instabile su entrambi gli assi (`../results_lqg_tuning/step1_flexible/sign_check.json`): vale solo leggendo $e=\hat y-r$ |
| Errore $e$ | Fig. 8.3 e testo: $r-\hat y$; Fig. 8.4, 8.17 e 8.22c: $r-y$ | $\hat e$, su richiesta |
| Stato $\hat y_i$ | stimato in 8.11, ma non osservabile da $y$ | stimatore senza integrale |
| Velocity loop DSN | §8.5 «proportional»; §8.5.2 PL con PI = DSN | verificare su altre fonti |
| $q_i$ finale 34 m | testo: «integral gain increased to 3»; Fig. 8.17b, ingrandita: 3 | 3 |
| PL, LMT | testo: 35 %, 1.2 s, 1.5 Hz; Tab. 8.3: 20 %, 1.4 s, 1.4 Hz | testo per il gradino da 0.01° |
| PP, LMT | testo: nessuna sovraelongazione, 1.8 Hz; Tab. 8.3: 20 %, 1.2 Hz | citare entrambi |
| LL, LMT | testo: 40 Hz, 0.0058 mdeg, «250 volte»; Tab. 8.3: 20 Hz, 0.004 mdeg | citare entrambi |

---

## 6. Confronto con la taratura fatta su LUCY

Taratura finale: azimut $k_p=14.8$, $k_i=24.1$; polare $k_p=129$, $k_i=123$.

| Aspetto | Libro | LUCY, stato attuale | Giudizio |
|---|---|---|---|
| Struttura | PI velocità + LQG posizione | uguale | conforme |
| Sintesi | `dlqr`, `dlqe`, $R=W=1$ | uguale, con `control` | conforme |
| Plant | identificato da prova | modello con $K_{pv}$, $K_{iv}$ stimati | **debole** |
| Coordinate | modali | modali, scala mia ($\lVert B_i\rVert=\lVert C_i\rVert$) | conforme; i $q_f$ non sono confrontabili con quelli del libro |
| $V$ | $=Q$ come primo tentativo | $\rho_e\,\mathrm{diag}(q_p,q_f)$ | ammesso; in azimut serve $\rho_e=10^{-3}$ |
| Ordine dei passi | $q_f$ all'inizio e alla fine | $q_f=0$, poi solo alla fine | **deviazione** |
| Arresto di $q_p$ | banda = fondamentale, accelerazione | 15 A, modulo ≥ 0.5 | adattamento motivato (§7) |
| Arresto di $q_i$ | oscillazioni, critico | indice, modulo, critico | parziale |
| $k_i/k_p$ | 0.2–0.8, tipico 0.5 (con banda ≤ 0.8 Hz) | 1.6 az, 0.95 el | non confrontabile direttamente: frazione del critico 0.44 az, 0.03 el |
| Prova di taratura | gradino, $T_{ry}$, $T_{wy}$ | doppietti reali | **deviazione** |
| Indice | 8 variabili | 3 (errore rms, $t_s$, sovraelongazione) | **manca il più importante** |
| $e_r$, $e_{d\max}$, $e_w$ | sì | no | **mancano** |
| $m_h$, $m_{d\max}$ | sì | non nell'indice | **mancano** |
| Ottimizzazione | vincolata, sui pesi | griglia | semplificazione |
| Feedforward, CPP, limiti | non nell'LQG DSN del cap. 8; sì nell'analisi LMT | no | conforme a DSN; **manca** rispetto a LMT |
| Attrito | non nella taratura | nel simulatore, non nello stimatore | oltre il libro |
| Ritardo | non trattato qui | 1 ciclo | oltre il libro |

Il ciclo limite in polare e la sovraelongazione sono ciò che il libro assegna al CPP. Il confronto finora fatto **dimostra solo che PL senza CPP, valutato sui doppietti, non conviene rispetto a una cascata con $K_{pp}$ più alto.**

---

## 7. Procedura proposta per LUCY

Le prime due voci sono la procedura del libro, le altre sono mie proposte di adattamento.

**A. Prima della taratura**
1. Misurare il plant $u\to y$ con la configurazione reale del PLC e del drive, oppure chiarire $K_{pv}$ e $K_{iv}$ effettivi.
2. Nessun feedforward né CPP nella taratura in stile DSN (§3.2). Se servono, si aggiungono dopo come variante separata, anche alla cascata di confronto; il CPP richiede di leggere §11.1.1.
3. Fissare le metriche di LUCY (8.27 adattata):
   - $e_r$ alla velocità di inseguimento reale, per esempio 0.06 deg/s;
   - $e_{d\max}$ e $m_{d\max}$ per coppia sui cuscinetti;
   - $e_w$ sostituito dall'errore rms sulla traiettoria reale: il vento può non essere il disturbo principale;
   - $e_o$, $t_s$ su un gradino piccolo e sui doppietti;
   - $f_o$ e $m_h$ con la banda da definire.

**B. Taratura (§8.2)**
1. $q_i$, $q_p$ piccoli; alzare $q_f$ solo dei modi che risuonano in $T_{ry}$. Se non risuonano, $q_f$ resta basso (bassa autorità).
2. Alzare $q_p$ fino al limite. Su LUCY la regola «banda = fondamentale» non vincola, perché i modi sono a 56 Hz: il limite è l'accelerazione, cioè la corrente, e il margine con i cuscinetti bloccati.
3. Alzare $q_i$ fino alle oscillazioni, partendo da $k_i/(0.25k_p^2)$ intorno a 0.5, con il critico come tetto. In polare il limite è già noto: il ciclo limite con l'attrito.
4. Ricontrollare le risonanze.
5. Taratura fine con (8.27), usando la griglia solo come punto di partenza.

**C. Criterio di esito**
Confronto con la cascata attuale e con la **cascata con $K_{pp}$ ritarato**, con le stesse aggiunte eventuali su entrambi gli schemi. Verifiche richieste: rumore e quantizzazione dell'encoder, poi prova sull'asse reale.

**Punti da verificare fuori dal capitolo 8**
- Struttura e guadagni del velocity loop DSN (rif. [7] citato in §8.3; Gawronski & Mellstrom 1994, rif. [5] del capitolo 1).
- $k_p$, $k_i$ usati alla DSS55 (rif. [9] di §8.3).
- Se la DSN usa feedforward con LQG in esercizio: il capitolo 8 non lo dice; il rif. [2] del cap. 7 (Gawronski & Mellstrom 1994, satelliti a bassa orbita) tratta il feedforward e non è stato letto.
- Se la pointing di LUCY è dominata dal vento o da altro: determina il peso di $e_w$.

---

## 8. Audit puntuale del §8.2 (pp. 108–110) e verifica su LUCY

Modello interno con Maxwell, plant simulato Dahl. Numeri e figure in `../results_lqg_tuning/step1_flexible/` e in `../report_lqg_tuning/`, sezioni «Audit» e «Step 1».

### 8.1 Cosa dice il testo

1. **Pesi flessibili.** Pesi piccoli su integrale, posizione e modi. Si guarda la funzione di trasferimento in anello chiuso; se le risonanze sono eccessive si alza il peso del modo corrispondente. Niente pesi inutilmente grandi: i modi sovrasmorzati peggiorano l'inseguimento.
2. **Peso proporzionale.** Si alza finché la banda raggiunge la frequenza fondamentale.
3. **Peso integrale.** Si alza finché compaiono oscillazioni, rispettando il guadagno critico.
4. **Correzione.** Se le risonanze riemergono, si alzano i pesi dei modi.

### 8.2 Cosa manca o non torna

| Punto | Libro | Conseguenza su LUCY |
|---|---|---|
| «Excessive» | Nessuna soglia. §8.4.2: «sharp resonance peaks», «visible flexible oscillations» | Confronto con i picchi di Fig. 8.5: circa 0 dB quello rimosso, circa −6 dB quello lasciato |
| Primo tentativo | §8.2 «small»; §8.4.2 $q_{ei}=0.1$, $q_e=0.3$, pesi di frequenza nulli | Stessi $q_i$, $q_p$; $q_f=0.1$ su tutti gli stati flessibili |
| Cosa osservare | §8.2 solo $T_{ry}$; §8.4.2 anche gradino e disturbo | $T_{ry}$ e gradino sul modello ridotto con Maxwell |
| Stimatore | Assente dal §8.2; $Q=V$ in (8.26), separati «may improve» a p. 119 | $Q=V$, con $V$ intensità della (8.8); $V$ separata si rimanda alla taratura fine |
| Scala di $q_f$ | $q_f=5$ contro $q_p=10$; $\lVert Q_t
Vert\gg\lVert Q_f
Vert$ (8.22) | Numeri non trasferibili: dipendono dalla scala modale e dal decadimento del modo. Qui $q_f=0.1$ non tocca il modo a 56 Hz, un effetto compare da circa $10^3$ |
| Indipendenza con $Q_f=0$ | $S_{cf}=S_{ctf}=0$ (§8.1.4) | Vale solo in parte: $K_f\neq0$ a $q_f=0$ e $k_p=1.35$ az, $1.05$ el contro 0.97 del rigido, per lo stato di Maxwell: senza cuscinetti è 0.97 su entrambi gli assi. La dipendenza dal peso del modo a 56 Hz è debole: meno dell'1 % a $q_f=10^4$ |
| Arresto del passo 2 | Banda = fondamentale | Non vincola con modi a 56 Hz; vale il limite di corrente (§8.1.7) |
| Condizione del passo 3 | Cita (8.15), che è la legge di controllo | È la (7.9): l'esempio di p. 104 ($k_p=16$, critico 64) la conferma |
| Plant | Lo stesso modello per progetto e valutazione | Uguale: modello ridotto con Maxwell (decisione dell'utente). Dahl è solo una verifica finale |

### 8.3 Esito del passo 1

Si segue il §8.2 alla lettera: TF e gradino sul solo modello ridotto con Maxwell, $Q=V$, pesi bassi, si alza il peso corrispondente solo se il picco è troppo alto (indicazioni dell'utente).

- **Correzione.** Le versioni precedenti passavano $V$ a `dlqe` per campione, come nel riquadro di p. 98. La (8.8) è continua: a $W_d=1$ il peso è $VT^2$. Lo stimatore era circa 760 volte più veloce (418 contro 0.55 rad/s). Sono **ritirate** le conclusioni «con $Q=V$ si perde la stabilità a pesi bassi» e «ciclo limite in elevazione alla first guess». Lo stesso vale per $\rho_e=10^{-3}$ di §6.
- **Vettore scelto.** $q_i=0.1$, $q_p=0.3$, $q_f=0.1$ su ogni stato flessibile: 11 elementi in azimut, 13 in elevazione.
- **Risultato.** Azimut $k_p=1.37$, $k_i=0.32$, banda 0.12 Hz, sovraelongazione 28 %, $t_s=12.1$ s; elevazione $k_p=1.05$, $k_i=0.32$, banda 0.20 Hz, 21 %, 9.6 s. Il libro: 8.52 s, 18 %, 0.21 Hz.
- **I picchi non sono troppo alti.** −59.5 dB az, −62.7 dB el; nella Fig. 8.5 il picco rimosso è circa a 0 dB e quello lasciato circa a −6 dB (letti dalla figura). Nessun peso corrispondente viene alzato.
- **Sovrasmorzamento.** Con lo stesso peso su tutti gli stati flessibili, da $q_f\geq1$ in azimut la banda scende e il tempo di assestamento cresce: a $q_f=10$ la banda è la metà. La causa è un solo stato, il polo reale lento aggiunto da Maxwell (0.05 Hz az, 0.4 Hz el), non i modi a 56 Hz: pesando uno stato alla volta, solo questo riproduce l'effetto.
- **Alzare il solo modo a 56 Hz** non cambia banda né gradino; sotto circa $10^3$ non ha effetto. Non è selezionato.
- **Verifica su Dahl, fuori dal libro.** Stabile su tutti i plant. Azimut a 1000 mdeg; elevazione circa 4 mdeg sopra il comando a 120 s, con $y-\hat y=4.15$ mdeg: è il bias della stima.

### 8.4 Passi 2-4, iterazione e $V$

Obiettivo dell'utente: assestamento basso, banda alta. Le scelte che il libro non fissa sono A1-A11 nel report (`../report_lqg_tuning/`).

- **Passo 2.** La regola della fondamentale non si applica: in azimut ferma prima la corrente ($q_p=86$, banda 0.34 Hz); in elevazione servirebbe $k_p\approx7400$ a 1 kHz, quindi lo stop è la stabilità a cuscinetti bloccati (scelta mia).
- **Passo 3.** Oscillazioni e (7.9) non arrivano mai per prime: fermano la corrente in azimut e la stabilità a cuscinetti bloccati in elevazione.
- **Passo 4.** In elevazione il picco a 56 Hz riemerge a −2.7 dB, come previsto dal libro: $q_f=3\cdot10^4$ sul solo modo a 56 Hz.
- **Iterazione.** Utile in azimut (scambio $q_p$ con $q_i$ a pari $k_p$), inutile in elevazione.
- **$V$.** A $Q$ fissa non cambia l'obiettivo sul modello. Non adottata: resta $Q=V$. In azimut $\rho=100$-1000 dimezza quasi l'errore su Dahl ma triplica l'assestamento: è un'opzione da riprendere con i disturbi.
