# Audit di Gawronski: LQG con la struttura attuale del driver (P di posizione, PI di velocità, K_f)

Materiale per sviluppare (non implementato) il controllore LQG di Gawronski su LUCY, mantenendo la struttura dei drive AKD2G:
- maglia di posizione P ($K_{pp}$);
- maglia di velocità PI ($K_{pv}$, $K_{iv}$);
- smorzamento attivo $K_f$ sui modi flessibili stimati.

Fonte: W. Gawronski, *Modeling and Control of Antennas and Telescopes*, Springer, 2008 (DOI 10.1007/978-0-387-78793-0). File: `/home/alessandro/Desktop/Thesis/978-0-387-78793-0.pdf`. I numeri di pagina sono quelli stampati nel libro; nel PDF vanno aggiunti 12.

Struttura del documento:
- §1 dove sta nel libro ciò che serve;
- §2 il metodo di Gawronski, fedele al testo e con i numeri delle equazioni;
- §3 l'adattamento a LUCY con la mappatura esatta LQR → guadagni AKD;
- §4 i numeri attuali di LUCY;
- §5 la procedura di sintesi passo passo per lo sviluppo;
- §6 vincoli e punti aperti.

---

## 1. Mappa del libro

| sezione (pagine) | contenuto | uso per LUCY |
|---|---|---|
| 1.2 (5–7) | maglia di velocità (ingresso: comando di velocità; uscita: encoder) e maglia di posizione | il *plant* del controllore di posizione è la maglia di velocità chiusa |
| 2.1 (11–12) | antenna rigida, maglia di velocità P: $G_{rl}=1/(1+Ts)$, $T=J/k_o$ | modello rigido di riferimento |
| 2.2.2–2.2.4 (14–21) | modello modale, stato modale $x_i=[\omega_iq_{mi},\ \dot q_{mi}]$, modi rigidi | coordinate per la sintesi |
| 2.2.5 (21–23) | modello discreto modale, condizione $\omega_i\Delta t\le\pi$ | master EtherCAT a 1 kHz |
| 3.2 (38–40) | modello della maglia di velocità identificato (n4sid, ordine 40 → ≤14) | identificazione per lo stimatore |
| 4.3–4.4 (47–50) | riduzione modale con norme $H_2$, $H_\infty$, Hankel; poli nell'origine mantenuti | scelta dei modi da tenere |
| 6.1 (73–76) | criteri di prestazione | metriche di taratura |
| **6.2 (77–79)** | **coordinate modali → posizione come primo stato → integrale della posizione** | **forma del modello per la sintesi** |
| 7.1–7.3 (81–87) | PI su antenna rigida: funzioni di trasferimento, guadagno integrale critico, passi di taratura, anello chiuso flessibile | analogia PI ↔ LQG, vincoli |
| 7.5 (91–93) | feedforward di velocità | equivalente a $K_{vff}$ dell'AKD |
| **8.1 (95–108)** | **LQG: CARE/FARE, LQG d'inseguimento, $K_c=[k_i\ k_p\ K_f]$, pesi, somiglianza col PI, limiti** | **nucleo della sintesi** |
| **8.2 (108–110)** | **passi di taratura dell'LQG** | **procedura** |
| 8.3.3–8.3.5 (113–117) | reiezione dei disturbi e banda, limiti di accelerazione | vincoli |
| 8.4 (117–124) | scelta dei pesi ($Q=V$ diagonali), GUI, ottimizzazione vincolata (8.27) | taratura fine |
| **8.5 (124–133)** | **LQG nella maglia di velocità; configurazioni PP, PL, LP, LL (LMT)** | **scelta della collocazione: la richiesta corrisponde a LP/LL** |
| 10.1 (145–148) | maglia singola: la PID equivale alla cascata | equivalenza PID ↔ cascata |
| 11.1.2 (162–165) | anti-windup $G_{aw}=k_{aw}/s$ | saturazioni |
| 11.2 (166–175) | modello d'attrito secco (11.9), dither a bassa velocità | attrito di LUCY (elevazione lentissima) |
| 12.1 (183–) | encoder ≠ fascio RF (posizione reale dell'asse ottico) | torsione encoder–asse ottico in azimut |

Non servono: i capitoli 5 (vento), 9 ($H_\infty$), 12.2–12.3 (monopulse/conscan) e 13 (tabella del binario).

---

## 2. Il metodo di Gawronski

### 2.1 Architettura (1.2, 3.2.1)

- Due maglie per asse: la **maglia di velocità** (drive e struttura; ingresso il comando di velocità $u_c$, uscita la posizione all'encoder $y$) e la **maglia di posizione** chiusa dal calcolatore tra comando $r$ ed encoder.
- Il controllore LQG "standard" del libro (DSN) **sostituisce il PI di posizione**: ingresso $r$ e $y$, uscita il comando di velocità. Il plant della sintesi è la maglia di velocità, identificata in prova (n4sid, ordine 40 ridotto a ≤14).
- Gawronski avverte che nella pratica la teoria "non dovrebbe" funzionare (prefazione, p. ix): poli rigidi nell'origine, rumori non gaussiani, riduzioni rigide non ottime. Eppure l'LQG funziona: "every LQG controller is an optimal controller, but not every one is acceptable".

### 2.2 Preparazione del modello (6.2, 2.2, 4.4)

**Tre trasformazioni** (6.2, p. 77):
1. **coordinate modali** $(A_m,B_m,C_m)$, con $A_m$ diagonale a blocchi 2×2 (2.23–2.25):
$$A_{mi}=\begin{bmatrix}0&\omega_i\\-\omega_i&-2\zeta_i\omega_i\end{bmatrix},\quad B_{mi}=\begin{bmatrix}0\\b_{mi}\end{bmatrix},\quad C_{mi}=\begin{bmatrix}\dfrac{c_{mqi}}{\omega_i}&c_{mvi}\end{bmatrix},\quad x_i=\begin{bmatrix}\omega_iq_{mi}\\\dot q_{mi}\end{bmatrix}$$
con $b_{mi}=\phi_{ik}/m_{mi}$ all'attuatore e $c_{mqi}=\phi_{ik}$ al sensore (2.21–2.22). Un modo rigido ha $A_{mi}=\begin{bmatrix}0&1\\0&0\end{bmatrix}$ (2.32);
2. **posizione come primo stato**: $x_p=[y;\ x_f]=Px_m$ con $P=\begin{bmatrix}C_m\\P_o\end{bmatrix}$, $P_o=[0_{n-1,1}\ I_{n-1}]$ (6.1–6.3);
3. **aumento con l'integrale della posizione**, per avere errore nullo a velocità costante e una reiezione migliore (6.4–6.7):
$$x_o=\begin{bmatrix}y_i\\y\\x_f\end{bmatrix},\qquad A_o=\begin{bmatrix}0&C_p\\0&A_p\end{bmatrix},\quad B_o=\begin{bmatrix}0\\B_p\end{bmatrix},\quad C_o=\begin{bmatrix}0&C_p\end{bmatrix}$$

**Riduzione** (4.3–4.4):
- Norme di un modo: $\|G_i\|_2\simeq\dfrac{\|B_{mi}\|\|C_{mi}\|}{2\sqrt{\zeta_i\omega_i}}$, $\|G_i\|_\infty\simeq\dfrac{\|B_{mi}\|\|C_{mi}\|}{2\zeta_i\omega_i}$, Hankel $=\tfrac12 H_\infty$ (4.4–4.6). L'$H_\infty$ è l'altezza del picco di risonanza.
- I modi con poli nell'origine (rigidi) hanno norma infinita e si **tengono sempre** (4.11–4.12); poi si tronca per norme decrescenti.

**Modello discreto** (2.34–2.37): $A_d=e^{A\Delta t}$, $B_d=\int_0^{\Delta t}e^{A\tau}B\,d\tau$. In coordinate modali i blocchi sono
$$A_{dmi}=e^{-\zeta_i\omega_i\Delta t}\begin{bmatrix}\cos\omega_i\Delta t&-\sin\omega_i\Delta t\\\sin\omega_i\Delta t&\cos\omega_i\Delta t\end{bmatrix},$$
validi se $\omega_i\Delta t\le\pi$ (cioè sotto Nyquist).

### 2.3 LQG (8.1.1, pp. 95–98)

| elemento | equazione |
|---|---|
| plant | $\dot x=Ax+Bu+v,\quad y=Cx+w$; $V=E(vv^T)$, $W=E(ww^T)=I$, $E(vw^T)=0$ (8.1) |
| stimatore | $\dot{\hat x}=A\hat x+Bu+K_e(y-C\hat x)$ (8.2) |
| indice | $J^2=E\!\int_0^\infty(x^TQx+u^TRu)\,dt$, con $R=1$ (scala nota di $Q/R$) (8.3) |
| guadagno di controllo | $u=-K_c\hat x$, $K_c=B^TS_c$, $A^TS_c+S_cA-S_cBB^TS_c+Q=0$ (CARE) (8.4–8.6) |
| guadagno dello stimatore | $K_e=S_eC^T$, $AS_e+S_eA^T-S_eC^TCS_e+V=0$ (FARE) (8.7–8.8) |
| controllore $y\to u$ | $A_{lqg}=A-BK_c-K_eC,\ B_{lqg}=K_e,\ C_{lqg}=-K_c$ (8.9–8.10) |
| Matlab | `Kc=lqr(A,B,Q,R)`, `Kc=dlqr(...)`; `Ke=lqe(A,I,C,Q,R)`, `Ke=dlqe(...)` |

Nota: nella (8.5) il libro scrive $B_o^T S_c$; per $R=1$ è $K_c=R^{-1}B^TS_c$.

### 2.4 LQG d'inseguimento (8.1.2–8.1.3, pp. 98–101; Figg. 8.3–8.4)

- Lo stato stimato è $\hat x=[\hat y_i,\ \hat y,\ \hat x_f]$ (8.11). Per inseguire si sostituisce $\hat y$ con l'errore $e=r-\hat y$ e $\hat y_i$ con $e_i=\int e\,dt$ (8.12).
- Il guadagno si divide come $K_c=[k_i\ \ k_p\ \ K_f]$ (8.14). La legge (8.15), letta con la convenzione di segno dell'anello chiuso (8.17), è
$$u=k_p\,e+k_i\,e_i-K_f\hat x_f,$$
cioè una **PI sulla posizione più una retroazione degli stati flessibili stimati** (Fig. 8.4: "PI controller" + "Flexible Mode controller").
- Anello chiuso (8.16–8.18), con $x_{cl}=[e_i;\ x;\ \hat x]$:
$$A_{cl}=\begin{bmatrix}0&0&-C\\Bk_i&A&-Bk_pC-Bk_fC_f\\Bk_i&k_eC&A-k_eC-Bk_pC-Bk_fC_f\end{bmatrix},\ B_{cl}=\begin{bmatrix}1\\Bk_p\\Bk_p\end{bmatrix},\ B_v=\begin{bmatrix}0\\I\\0\end{bmatrix},\ C_{cl}=[0\ \ C\ \ 0].$$

### 2.5 Pesi (8.1.4, 8.1.6, 8.4.1)

- **Partizione** inseguimento/flessibile (8.19–8.22): $x=[x_t;\ x_f]$ con $x_t=[e_i;\ e]$, $Q=\mathrm{diag}(Q_t,Q_f)$.
- **Bassa autorità**: $\|Q_t\|\gg\|Q_f\|$. Per $Q_f=0$ si ha $S_{cf}=S_{ctf}=0$, e il guadagno d'inseguimento non dipende dalla parte flessibile; con $Q_f$ piccolo la dipendenza resta debole (rif. Collins et al. 1994).
- **Q diagonale** in coordinate modali: $Q=\mathrm{diag}(q_i,\ q_p,\ q_f)$, con $q_f=[q_{f1m},q_{f2m},\dots]$ e $q_{fim}=[q_{fi}\ \dot q_{fi}]$ per modo (8.23–8.25).
- **Semplificazione GUI**: $Q=V=\mathrm{diag}(q)$, $q=[q_{ei}\ q_e\ q_{f1}\ q_{f1}\ \cdots\ q_{fn}\ q_{fn}]$ (8.26). Con $n=4$ modi sono 10 pesi.
- **Effetto dei pesi** (Tab. 8.2, p. 121):

| peso | effetto |
|---|---|
| proporzionale $q_e$ | tempo di assestamento, banda, reiezione dei disturbi |
| integrale $q_{ei}$ | sovraelongazione, errore a velocità costante, reiezione a bassa frequenza |
| modo $i$ ($q_{fi}$) | ampiezza della vibrazione di quel modo (smorzamento attivo) |

- Proprietà (8.1.6):
  - aumentare $q_p$ aumenta banda, $k_p$ e anche $k_i$;
  - aumentare $q_i$ migliora la reiezione senza cambiare la banda;
  - $q_p$ e $q_i$ sono accoppiati (Figg. 8.9–8.10, che danno i guadagni in funzione dei pesi per $k_i/k_p=0.2\dots0.8$).

### 2.6 Somiglianza col PI, limiti, passi di taratura (8.1.5, 8.1.7, 8.2)

- **Somiglianza** (p. 104): con i pesi giusti, l'antenna flessibile con LQG si comporta come un'antenna rigida con PI. $K_f$ sopprime le vibrazioni e permette $k_p$ alti che con il solo PI destabilizzerebbero i modi (Fig. 7.7). L'esempio usa $k_p=9.5$, $k_i=6.3$.
- **Limiti** (8.1.7):
  1. pesi flessibili non eccessivi, altrimenti i modi diventano sovrasmorzati e banda e reiezione peggiorano (Fig. 8.8);
  2. il peso di posizione è limitato dall'accelerazione o dalla coppia disponibile;
  3. il peso integrale non deve superare il **guadagno integrale critico**: $k_i\le0.25\,k_p^2$ (7.9–7.10), con frequenza critica $\omega_o=0.5k_p=\sqrt{k_i}$ (7.11).
- **Passi di taratura** (8.2, pp. 108–110):
  1. pesi flessibili piccoli, $q_i$ e $q_p$ piccoli; alzare i pesi dei modi che risuonano nella funzione d'anello chiuso;
  2. alzare $q_p$ finché la banda raggiunge la frequenza fondamentale della struttura;
  3. alzare $q_i$ finché compaiono oscillazioni, rispettando il limite critico;
  4. correggere i pesi flessibili se le risonanze riemergono.
- **Taratura fine** (8.4.3): ottimizzazione vincolata dei pesi sull'indice (8.27),
$$f=w_1e_r+w_2e_{d\max}+w_3e_w+w_4e_o+w_5t_s+w_6(3-f_o)+w_7m_h+w_8m_{d\max}$$
  - $e_r$: errore a velocità costante;
  - $e_{d\max}$: errore massimo per un disturbo a gradino;
  - $e_w$: errore rms nel vento;
  - $e_o$: sovraelongazione;
  - $t_s$: tempo di assestamento;
  - $f_o$: banda;
  - $m_h$: picco ad alta frequenza;
  - $m_{d\max}$: picco della funzione di disturbo.

### 2.7 PI su antenna rigida (7.1–7.3, 7.5)

- $G=1/s$, $K=k_p+k_i/s$:
$$T_{ry}=\frac{k_ps+k_i}{s^2+k_ps+k_i},\quad T_{wy}=\frac{s}{s^2+k_ps+k_i}\quad(7.4)$$
- $k_p$ allarga la banda e riduce $T_{wy}$, ma fa crescere l'accelerazione (7.6–7.7).
- $k_i$ agisce solo su $T_{wy}$ a bassa frequenza ($T_{wy}\simeq s/k_i$).
- Taratura (7.2): alzare $k_p$ fino al limite di accelerazione e tenerlo al 75 %; alzare $k_i$ fino a oscillazioni o sottoelongazione e tenerlo al 75 %, sotto il valore critico.
- Anello chiuso flessibile con PI (7.12–7.14):
$$A_{cl}=\begin{bmatrix}0&-C\\Bk_i&A-Bk_pC\end{bmatrix},\quad B_{cl}=\begin{bmatrix}1\\Bk_p\end{bmatrix}$$
- **Feedforward** (7.5): la derivata del comando inviata alla maglia di velocità approssima l'inversa della maglia ($G_{ff}=s$). Migliora l'inseguimento (banda da 0.1 a 1.4 Hz sulla DSN) ma **non la reiezione dei disturbi**.

### 2.8 Collocazione dell'LQG: maglia di velocità (8.5, pp. 124–133)

- **Banda in cascata** (8.28–8.29): con la maglia di velocità $G=1/(s(Ts+1))$ e la posizione P di guadagno $k$, l'anello chiuso è $H=\dfrac{k/T}{s^2+s/T+k/T}$. *La banda della maglia di posizione dipende da quella della maglia di velocità*: se la seconda è stretta, alzare $k$ dà solo oscillazioni.
- **Quattro configurazioni** sul telescopio LMT (Tab. 8.3, a 1 ms di campionamento):

| sistema | velocità | posizione | assestamento [s] | banda [Hz] | errore rms nel vento [mdeg] |
|---|---|---|---|---|---|
| PP | PI | PI | 3.0 | 1.2 | 1.48 |
| PL | PI | LQG | 1.4 | 1.4 | 0.76 |
| LP | LQG | PID | 0.6 | 20 | 0.15 |
| LL | LQG | LQG | 0.5 | 20 | 0.004 |

- In LP/LL l'LQG di velocità è tarato sul modello ad anello aperto della maglia di velocità (retroazione di velocità tolta), con ingresso coppia/corrente. La posizione è PID con feedforward, preprocessore di comando e limiti.
- **La richiesta per LUCY (PI di velocità, P di posizione e $K_f$, tutti dalla sintesi) corrisponde alla famiglia LP/LL.** È la collocazione che il libro trova migliore: guadagno da 10 a 250 volte sull'errore nel vento, con l'avvertenza che i numeri vengono da modelli analitici.

### 2.9 Aspetti non lineari (10, 11)

- **Maglia singola** (10.1): con il solo anello di posizione e ingresso coppia serve una PID ($k_d>0$). La cascata P–PI con velocità P è equivalente a una PID (10.3). È la base della mappatura del §3.
- **Anti-windup** (11.8): retroazione di $u_o-u$ attraverso $G_{aw}=k_{aw}/s$ sull'integratore. Il modello AKD di `lucy_cl` ha già una back-calculation (`T_AW`).
- **Preprocessore di comando** (11.1.1): limita velocità e accelerazione del comando nei grandi spostamenti.
- **Attrito** (11.9): Coulomb con stiction e soglia di velocità $\nu_t$; **dither** per l'inseguimento a bassissima velocità (11.2.2). È rilevante per l'elevazione di LUCY sulla circonferenza (≤ 0.06 °/s).
- **Encoder ≠ asse del fascio** (8.5 conclusioni, 12): l'encoder misura la propria posizione, non quella del fascio; la cedevolezza tra i due è una fonte d'errore. Su LUCY è la torsione d'azimut tra motore/encoder e forcella ($1/k_{\psi,2}=4.4$ mdeg per N m).

---

## 3. Adattamento a LUCY mantenendo la struttura AKD

### 3.1 La legge dell'AKD è una PID sull'errore

Legge del drive (`lucy_cl/driver.py`, `Drive.step`), trascurando filtri e saturazioni, con $\omega_i=2\pi K_{iv}$:
$$v_{cmd}=K_{pp}\,e+K_{vff}\,\dot r,\qquad i_{cmd}=K_{pv}(v_{cmd}-v)+K_{pv}\,\omega_i\!\int(v_{cmd}-v)\,dt .$$

Con $K_{vff}=1$ (feedforward attivo) e $\int(\dot r-v)\,dt=e$ si ottiene **esattamente**
$$\boxed{\,i_{cmd}=k_D\,\dot e+k_P\,e+k_I\!\int e\,dt,\qquad k_D=K_{pv},\quad k_P=K_{pv}(K_{pp}+\omega_i),\quad k_I=K_{pv}K_{pp}\,\omega_i\,}$$
cioè una PID sull'errore d'inseguimento: lo stesso stato $[e_i,\ e,\ \dot e]$ dell'LQG d'inseguimento di Gawronski, esteso alla derivata perché qui l'ingresso è la corrente (coppia) e non la velocità.

Con $K_{vff}=0$ il termine derivativo diventa $-K_{pv}v$: la legge differisce per $K_{pv}\dot r$, cioè per il feedforward.

### 3.2 Mappatura inversa: dai guadagni LQR ai guadagni AKD

Dati i guadagni LQR sullo stato rigido $[e_i,\ e,\ \dot e]$, cioè $K_{c,rig}=[k_I\ \ k_P\ \ k_D]$:
$$K_{pv}=k_D,\qquad K_{pp},\ \omega_i=\text{radici di }\ z^2-\frac{k_P}{k_D}z+\frac{k_I}{k_D}=0,\qquad K_{iv}=\frac{\omega_i}{2\pi}.$$

- **Rappresentabilità** come cascata (radici reali): $\left(\dfrac{k_P}{k_D}\right)^2\ge4\dfrac{k_I}{k_D}$. È l'analogo del guadagno integrale critico di Gawronski ($k_i\le0.25k_p^2$, 7.9). Se la sintesi lo viola, si riduce il peso integrale (passo 3 del §2.6).
- **Scelta della radice.** Le due assegnazioni ($K_{pp}$ alla radice maggiore o minore) danno la stessa PID lineare, ma differiscono con saturazioni e anti-windup e con i filtri sulla velocità. In I5 l'azimut ha $K_{pp}>\omega_i$, l'elevazione $K_{pp}<\omega_i$. Si sceglie verificando le due versioni nel modello non lineare.
- **$K_f$** non ha un posto nell'AKD: è una corrente aggiuntiva $i_f=-K_f\hat x_f$, calcolata nel master (stimatore) e inviata come offset di corrente. È lo stesso canale già usato per il feedforward dell'osservatore di Hilkert (`i_exc` del modello del drive, RxPDO nel modello EtherCAT).

### 3.3 Formulazione della sintesi (variante richiesta, tipo LP/LL)

Per ogni asse (gli assi sono disaccoppiati, come in Gawronski 1.2):

1. **Plant**: corrente → posizioni, dal plant H6 linearizzato (`GimbalPlant.linear_matrices`, `rigid_inertia`), con
   - ingresso $\kappa_iK_t\,i$ sul nodo motore;
   - uscite l'encoder ($C_{enc}$) e, se si vuole, l'asse ottico ($C_{LOS}$) o il giroscopio sul telescopio ($C_{gyro}$, velocità);
   - attrito **escluso** dalla linearizzazione: in strisciamento non ha rigidezza ed entra come disturbo (in `linearise` la rigidezza $\sigma_0$ del Dahl creerebbe un falso modo a circa 7 Hz);
   - molla dei cavi trascurabile o tenuta come termine statico ($K_c\approx1.7$ N m/rad).
2. **Coordinate modali** (§2.2): modo rigido più i modi flessibili scelti con le norme di Gawronski (§4.2). Posizione come primo stato e aumento con l'integrale, nella forma d'inseguimento:
$$x=\begin{bmatrix}e_i\\e\\\dot e\\x_f\end{bmatrix},\qquad A=\begin{bmatrix}0&1&0&0\\0&0&1&0\\0&0&0&*\\0&0&0&A_f\end{bmatrix},\quad B=\begin{bmatrix}0\\0\\-b\\-B_f\end{bmatrix},\quad b=\frac{\kappa_iK_t}{J}.$$
   Il segno di $B$ dipende dalla convenzione $e=r-y$. Il blocco rigido è un **triplo integratore**, uno in più del caso di Gawronski, che ha ingresso in velocità.
3. **Pesi** $Q=\mathrm{diag}(q_i,\ q_p,\ q_d,\ q_f)$, $R=1$ in A². Si parte dalla bassa autorità $q_f\ll(q_i,q_p,q_d)$ e si usa il §4.1 come punto di partenza (i pesi che riproducono I5).
4. **CARE** → $K_c=[k_I\ k_P\ k_D\ K_f]$ → §3.2 → $K_{pp},K_{pv},K_{iv}$ per l'AKD, e $K_f$ per il master.
5. **Stimatore** (FARE) sugli stati flessibili e di posizione:
   - misura dell'encoder (rumore $W_{enc}=\sigma_{enc}^2$);
   - opzionale il giroscopio ADIS16507 già modellato ($W_g=N^2/\Delta t$, §4 del report sull'osservatore);
   - $V$ dal disturbo di coppia d'attrito ($\sigma_T\approx0.1$–$0.4$ N m → $V_{\dot e}=(b\,\sigma_T)^2$);
   - Gawronski usa $Q=V$ come primo tentativo (8.26).
6. **Discretizzazione** a 1 kHz (master EtherCAT) per stimatore e $K_f$: $A_d=e^{A\Delta t}$ (2.34), `dlqe`/`dlqr`. I guadagni dell'AKD restano continui, perché le maglie interne del drive sono molto più veloci del ciclo EtherCAT.

### 3.4 Variante del libro (PL), come riferimento

- La maglia di velocità AKD resta quella attuale e fa parte del plant (ingresso $v_{cmd}$, uscita encoder; identificabile come in 3.2).
- L'LQG di posizione dà $[k_i\ k_p\ K_f]$: $K_{pp}=k_p$, `PL.KI` $=k_i$ (l'integratore di posizione esiste nell'AKD e nel modello, oggi a zero), e $K_f\hat x_f$ come offset di velocità.
- **Non** fornisce $K_{pv}$, $K_{iv}$: la richiesta corrisponde alla variante del §3.3.

### 3.5 Opzioni utili oltre il libro

- **Stato di disturbo nello stimatore** (estensione suggerita da Hilkert 2011 e coerente con la FARE): un modello interno della coppia d'attrito (passeggiata aleatoria) rende stimabili il disturbo e la **torsione d'azimut** $T_d/k_{\psi,2}$.
- **Uscita sull'asse ottico**: con quello stato, l'errore d'inseguimento si può calcolare sull'asse ottico stimato ($r-C_{LOS}\hat x$) invece che sull'encoder. Toglierebbe il residuo di 2.7 mdeg visto nello studio con l'osservatore, che una legge sull'encoder non può togliere.
- **Ritardo EtherCAT**: $K_f$ arriva al drive un ciclo dopo. Il modello discreto dello stimatore va aumentato con uno stato di ritardo dell'ingresso, altrimenti la fase persa ai modi flessibili annulla lo smorzamento attivo (§4.2).

---

## 4. Numeri attuali di LUCY (modello in uso)

Ottenuti da:
- plant `config_vcl_fit.yaml`, linearizzato senza attrito né cavi;
- inerzie rigide $J=\mathbf r^T\mathbf M\mathbf r$;
- scala di coppia identificata su V2;
- guadagni I5.

### 4.1 Guadagni attuali come PID equivalente (punto di partenza dei pesi)

| asse | $J$ [kg m²] | $b=\kappa_iK_t/J$ [rad/s²/A] | $K_{pp}$ [1/s] | $K_{pv}$ [A s/rad] | $K_{iv}$ [Hz] | $k_D$ [A s/rad] | $k_P$ [A/rad] | $k_I$ [A/(rad s)] | poli rigidi in anello chiuso [Hz / ζ] |
|---|---|---|---|---|---|---|---|---|---|
| azimut | 2.162 | 0.1151 | 6.98 | 378.7 | 0.26 | 378.7 | 3265 | 4330 | 5.17/1.00, 1.51/1.00, 0.26/1.00 |
| elevazione | 0.900 | 0.3072 | 6.78 | 53.8 | 2.88 | 53.8 | 1338 | 6599 | 2.97/0.29 (coppia), 0.93/1.00 |

- Discriminante $(k_P/k_D)^2-4k_I/k_D$: azimut 28.6, elevazione 127.3. Entrambe le tarature sono rappresentabili come cascata.
- I poli sono quelli di $s^3+b(k_Ds^2+k_Ps+k_I)$, con maglia di corrente ideale e senza filtri.
- L'elevazione attuale ha una coppia di poli rigidi a circa 3 Hz poco smorzata ($\zeta=0.29$). Un LQR la sposterebbe naturalmente.

### 4.2 Modi flessibili per asse (uscita: posizione; ingresso: corrente)

| asse | modo | $f$ [Hz] | $\zeta$ | $H_\infty$ corrente → encoder [µrad/A] | → asse ottico [µrad/A] | $\omega\Delta t$ a 1 kHz [rad] | fase di 1.5 ms di ritardo [°] |
|---|---|---|---|---|---|---|---|
| azimut | rigido | 0 | – | ∞ | ∞ | – | – |
| azimut | 1 | 737 | 0.88 | 13.6 | 0.002 | 4.63 (> π) | 398 |
| azimut | 2 | 1647 | 0.85 | 2.8 | 0.0002 | 10.4 (> π) | 890 |
| elevazione | rigido | 0 | – | ∞ | ∞ | – | – |
| elevazione | 1 | **210** | **0.024** | 31.5 | **4.7** | 1.32 | 113 |
| elevazione | 2 | 404 | 0.16 | 1.2 | 0.04 | 2.54 | 218 |
| elevazione | 3 | 871 | 0.24 | 0.05 | 0.0004 | 5.47 (> π) | 470 |

Nodi del modello: azimut $[\psi_e,\psi_m,\psi_L]$, con encoder $\psi_e$, motore $\psi_m$, asse ottico $\psi_L$; elevazione $[\theta_e,\theta_h,\theta_m,\theta_r,\theta_p]$, con encoder $\theta_e$, motore $\theta_m$, asse ottico $\theta_p$. Lo smorzamento è proporzionale alla rigidezza ($\zeta_1=0.1$ di configurazione).

**Lettura secondo Gawronski:**
- **Elevazione: si tiene il modo 1** (210 Hz, poco smorzato, visibile all'asse ottico). Il modo 2 è al limite, gli altri si troncano.
- **Azimut: nessun modo flessibile utile** al master a 1 kHz. I modi sono oltre Nyquist, molto smorzati e quasi invisibili all'asse ottico: $K_f$ in azimut è vuoto. L'errore all'asse ottico d'azimut è la torsione quasi statica, da trattare come al §3.5.
- **Scala dei tempi diversa dal libro.** Le antenne del libro hanno il primo modo a 1–4 Hz, vicino alla banda, ed è lì che $K_f$ permette $k_p$ più alti (§2.6). Su LUCY la banda rigida è di pochi Hz e il primo modo del modello è 40–70 volte più in alto: **la parte flessibile è a bassissima autorità per costruzione**, e i guadagni rigidi si possono progettare quasi indipendenti ($Q_f\to0$, $S_{ctf}\to0$, 8.22).
- **Il modo misurato d'elevazione è circa 95 Hz**, non 210 (Bode di ottobre). A 95 Hz il ritardo di 1.5 ms vale circa 51°, quindi lo smorzamento attivo dal master è possibile, purché lo stimatore includa il ritardo. Per una sintesi realistica serve un modello H6 con il modo a 95 Hz, oppure il modello identificato della maglia di velocità come fa Gawronski (3.2: n4sid sui dati del Bode).

---

## 5. Procedura di sintesi (per lo sviluppo)

Da eseguire per asse. La simulazione e le metriche sono quelle già pronte per lo studio dell'osservatore di Hilkert: circonferenza di Attorrese, EtherCAT a 1 ms, ADIS16507.

1. **Modello di sintesi**:
   - linearizzazione del plant H6 senza attrito né cavi;
   - scala di coppia $\kappa_i$ sull'ingresso;
   - analisi modale per asse;
   - norme dei modi (4.4–4.6) e scelta dei modi (§4.2);
   - forma modale con la posizione come primo stato e l'integrale (6.1–6.7) in forma d'inseguimento (§3.3).
   - *Decisione*: H6 attuale (modo a 210 Hz) o H6 corretto sul modo a 95 Hz / modello identificato.
2. **Pesi iniziali** che riproducono I5: con $Q_f=0$ e $R=1$, cercare $q_i,q_p,q_d$ tali che $K_{c,rig}\approx[k_I\ k_P\ k_D]_{I5}$ (§4.1). È la "linea di base" equivalente al controllo attuale.
3. **Taratura** secondo i passi del §2.6:
   - $q_f$ del modo d'elevazione finché la risonanza nell'anello chiuso si spegne, senza sovrasmorzarlo;
   - $q_p$ fino al limite (corrente massima 18 A, accelerazione, banda della maglia di corrente, filtri AKD);
   - $q_i$ fino alle oscillazioni, rispettando la condizione del §3.2.
4. **Guadagni AKD** dalla mappatura del §3.2: $K_{pv}=k_D$, e $K_{pp}$, $K_{iv}$ dalle radici. Scegliere l'assegnazione delle radici; verificare $K_{vff}=1$.
5. **Stimatore**:
   - FARE continua o `dlqe` a 1 kHz, con uscite encoder (e giroscopio);
   - $V$, $W$ come al §3.3;
   - modello aumentato con il ritardo di un ciclo;
   - opzionale lo stato di disturbo (§3.5).
6. **Verifica lineare**:
   - anello chiuso (8.16–8.18) esteso alla struttura AKD reale (filtri AR1–AR4, observer OBS, maglia di corrente);
   - margini (`closed_loop.loop_frf`, `margins`);
   - autovalori del modello linearizzato completo;
   - funzioni $T_{ry}$, $T_{wy}$, $T_{ra}$, $T_{wa}$ (7.2–7.3) confrontate con I5.
7. **Verifica non lineare**:
   - modello time-domain completo (attrito Dahl+Stribeck, cavi, saturazioni, anti-windup) con EtherCAT a 1 ms ($K_f$ e stimatore nel master, guadagni rigidi nel drive);
   - metriche: errore di puntamento e jitter sopra 1/10/100 Hz sulla circonferenza di Attorrese;
   - confronto con: attuale, attuale + osservatore di Hilkert, LQG.
8. **Taratura fine**, se serve: ottimizzazione vincolata dei pesi sull'indice (8.27), adattato con $e_w$ → errore rms di puntamento e jitter.

Corrispondenze Matlab → Python (scipy):

| Matlab | Python |
|---|---|
| `lqr(A,B,Q,R)` | `S = scipy.linalg.solve_continuous_are(A, B, Q, R)`; `K = solve(R, B.T @ S)` |
| `dlqr` | `solve_discrete_are(Ad, Bd, Q, R)`; `K = solve(R + Bd.T S Bd, Bd.T S Ad)` |
| `lqe(A,G,C,V,W)` | `S = solve_continuous_are(A.T, C.T, G V G.T, W)`; `Ke = S C.T W^-1` |
| `dlqe` | `solve_discrete_are(Ad.T, C.T, V, W)`, poi il guadagno di Kalman |
| `c2d` | `scipy.linalg.expm` della matrice aumentata, oppure `scipy.signal.cont2discrete` |
| `n4sid` | nessun equivalente diretto in scipy |

---

## 6. Vincoli e punti aperti da decidere

1. **Modello per la sintesi.** Il modo d'elevazione di H6 (210 Hz) non è quello misurato (circa 95 Hz). Con il modello attuale $K_f$ ha poco da fare; con quello corretto è il caso in cui l'LQG di Gawronski dà più valore.
2. **Dove calcolare $K_f$.** L'AKD non ha una retroazione di stato configurabile, quindi va nel master a 1 kHz con un ciclo di ritardo. Lo smorzamento attivo è realistico solo per modi sotto circa 100 Hz, e lo stimatore deve modellare il ritardo.
3. **Uscite dello stimatore.** Solo encoder (come nel libro) oppure encoder + giroscopio. Il giroscopio vede l'asse ottico e quindi il modo e la torsione, ma porta il rumore MEMS già quantificato.
4. **Feedforward di velocità.** La mappatura esatta del §3.1 richiede $K_{vff}=1$, oggi spento nelle prove.
5. **Attrito.** Il libro lo tratta come disturbo (strisciamento) e propone il dither. L'elevazione di LUCY sulla circonferenza è sempre a bassissima velocità; lo stato di disturbo nello stimatore (§3.5) è il modo per dare all'LQG la reiezione che oggi dà l'osservatore di Hilkert.
6. **Filtri dell'AKD.** La PID equivalente ignora AR1–AR4 e l'observer di velocità. La verifica (passo 6 del §5) va fatta sul modello completo del drive.

---

### Riferimenti nel libro (pagine stampate)

- Cap. 1.2 (pp. 5–7): struttura delle maglie.
- 2.1 (11–12): modello rigido.
- 2.2.2–2.2.5 (14–23): modello modale e discreto.
- 3.2 (38–40): identificazione della maglia di velocità.
- 4.3–4.4 (47–50): riduzione modale.
- 6.1–6.2 (73–79): criteri e trasformazioni.
- 7.1–7.5 (81–93): PI e feedforward.
- 8.1–8.5 (95–133): LQG.
- 10.1 (145–148): maglia singola.
- 11.1.2 (162–165): anti-windup.
- 11.2 (166–175): attrito.
- 12.1 (183–): controllo del fascio.

Riferimenti citati da Gawronski utili per approfondire:
- Collins, Haddad, Ying (1994), LQG a bassa autorità (8.22);
- Athans (1971) e Johnson (1968), PID da regolatore ottimo con integrale (6.2.3);
- Gawronski & Souccar (2005), sistemi di controllo dell'LMT (8.5);
- Maneri & Gawronski (2000), taratura dell'LQG con GUI (8.4).
