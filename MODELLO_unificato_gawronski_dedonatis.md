# Modello unificato: corpo rigido (De Donatis) e modello vibrazionale multi-GdL nel formalismo modale di Gawronski

> **Nota di realismo.** Le identità di questo documento sono esatte e valgono per qualunque valore dei parametri. Tutti i **numeri** (frequenze, bande di validità, norme modali) vengono invece da parametri **segnaposto** usati solo per verificare le identità: non caratterizzano il sistema reale. Gli script di verifica sono `parte_A_simbolica.py` (SymPy) e `parte_B_numerica.py` (NumPy/SciPy).

## 0. L'idea in breve

Nel capitolo 2 di Gawronski (*Modeling and Control of Antennas and Telescopes*, 2008) il modello di una struttura nasce dal modello FEM ed è scritto in **coordinate modali**. Il **modo rigido** (frequenza nulla, $\mathbf K\boldsymbol\phi_{rb} = \mathbf 0$, eq. 2.30) è trattato come un modo qualunque, con il proprio blocco di stato (eq. 2.32). Applicato al tuo modello, questo dà tre fatti:

1. ogni catena torsionale (azimut, polare) ha **esattamente un modo rigido**, $\boldsymbol\phi_0 = \mathbf 1$, e la sua coordinata modale $\eta_0$ è un angolo di corpo rigido;
2. il modello di De Donatis è **la proiezione esatta del tuo modello sui due modi rigidi**: se si tengono solo $\eta_{\psi 0}$ e $\eta_{\theta 0}$ si ottengono le sue equazioni (3.13)–(3.18) e (4.7);
3. il resto del tuo modello è una somma di **oscillatori modali alla Gawronski** (eq. 2.14). L'unico termine non lineare, il payload, è il corpo rigido di De Donatis **agganciato ai nodi L e P** tramite le righe della matrice modale in quei nodi: gli stessi vettori che Gawronski usa per encoder e attuatori (eq. 2.21–2.22).

Ne nasce una **gerarchia** di modelli, ordinata dalle norme modali di Gawronski (cap. 4): "De Donatis" (solo modi rigidi) → "De Donatis + $n$ modi" → tuo modello ridotto completo a 7 GdL. Ogni livello è un sottoinsieme esatto del successivo.

```
FEM per asse (spec §7) ──Guyan──▶ catene 3+4 GdL (spec §8) ──Φ (Gawronski §2.2)──▶ η = [modi rigidi | modi flessibili]
                                                                             │
                     payload = corpo rigido di De Donatis su (ψ_L, θ_p) = E_η η
                                                                             │
             solo modi rigidi ⇒ De Donatis (3.13)/(4.7)      tutti i modi ⇒ tuo modello accoppiato (spec §9)
```

## 1. Cosa si prende da ciascuna fonte

| Fonte | Parte usata | Ruolo nel modello unificato |
|---|---|---|
| Gawronski §2.2.1–2.2.2 (eq. 2.4–2.22) | FEM → modale; $\mathbf B_m = \mathbf M_m^{-1}\boldsymbol\Phi^{\mathsf T}\mathbf B_o$, $\mathbf C_{mq} = \mathbf C_{oq}\boldsymbol\Phi$ | Base modale di ciascuna catena; vettori di attuatore ed encoder |
| Gawronski §2.2.3–2.2.4 (eq. 2.23–2.32) | Stato modale a blocchi 2×2, modi rigidi | Forma di stato del modello linearizzato |
| Gawronski §2.3 (eq. 2.43–2.47) | Motore e riduttore ($J_m$, $N$, $k_g$) | Il nodo motore della catena; legame con il riduttore del banco di De Donatis |
| Gawronski §3.2.4 | "Il modello di azimut dipende dall'elevazione, quello di elevazione no" | Risultato che il modello unificato **dimostra** analiticamente (§7.3) |
| Gawronski §4.3–4.4 (eq. 4.4–4.12) | Norme $H_2$, $H_\infty$, di Hankel per modo; i modi rigidi si tengono sempre | Criterio per scegliere i modi flessibili della gerarchia (§8) |
| De Donatis cap. 3 (eq. 3.1–3.27) | Cinematica e dinamica di Lagrange del gimbal rigido sbilanciato | Termini non lineari del payload $\mathbf M_p, \mathbf C_p, \mathbf g_p$ |
| De Donatis §4.2 (eq. 4.7) | Inerzia del rotore e attrito riportati con $n_{GB}^2$, coppia $\eta_{GB}\,n_{GB}\,\tau$ | Inerzie e smorzamenti al nodo motore |
| Spec multi-GdL §6–10 | FEM, Guyan, $T$, $V$, $\mathcal D$, linearizzazione | Matrici $\mathbf M_r, \mathbf K_r, \mathbf C_r$ e convenzioni |

## 2. Corrispondenza tra le convenzioni (verificata, parte A)

**Angoli.** L'azimut ψ coincide (antiorario visto dall'alto, attorno a Z). L'angolo polare di De Donatis si misura dallo zenit (linea di vista $\hat K$, eq. 3.4). Quello dello spec si misura dall'orizzonte, positivo verso il basso. Entrambe le rotazioni sono positive attorno allo stesso asse $y$:

$$\theta_{DD} = \theta + \frac{\pi}{2}$$

**Terne.** La terna del telescopio di De Donatis è la terna B ruotata di $+\pi/2$ attorno a $y_B$: $\hat I = -\mathbf z_B$, $\hat J = \mathbf y_B$, $\hat K = \mathbf x_B$ (linea di vista).

**Tensore del payload.** De Donatis dà $\mathbf I_2$ nel baricentro, con il baricentro sulla linea di vista a distanza $\rho_2$ dall'asse polare. Il tensore per O in terna B (spec) vale

$$\mathbf I_{p,O}\big|_B = \begin{bmatrix} I_{2,33} & I_{2,23} & -I_{2,13}\\ I_{2,23} & I_{2,22} + m_2\rho_2^2 & -I_{2,12}\\ -I_{2,13} & -I_{2,12} & I_{2,11} + m_2\rho_2^2 \end{bmatrix}, \qquad d_x = \rho_2,\quad d_z = 0$$

**Verifica (A5).** Con questa mappa, $\mathbf M_{DD}(\theta + \pi/2) = \mathbf M_p(\theta) + \operatorname{diag}(I_{1,33}, 0)$ e $\mathbf g_{DD}(\theta + \pi/2) = \mathbf g_p(\theta)$, identicamente.

**Corrispondenza dei parametri** (dimostrata in §6):

| De Donatis | Modello unificato |
|---|---|
| $I_{1,33} + n_{GB}^2 J_m$ (riga ψ) | $J_{\psi,s} = \mathbf 1^{\mathsf T}\mathbf M_{r,\psi}\mathbf 1 + J^{ax}$: tutta l'inerzia di azimut non dovuta al payload (forcella, albero, encoder, rotore, parte trasversale degli alberi polari) |
| $I_{2,22} + m_2\rho_2^2 + n_{GB}^2 J_m$ (riga θ) | $I_{yy,p} + J_{\theta,s}$, con $J_{\theta,s} = \mathbf 1^{\mathsf T}\mathbf M_{r,\theta}\mathbf 1$ |
| $n_{GB}^2 b_m$ | $\mathbf 1^{\mathsf T}\mathbf C_b\mathbf 1$ (attrito viscoso al nodo motore o ai cuscinetti) |
| $\eta_{GB}\,n_{GB}\,\tau_{motor}$ | $\tau$ al nodo motore (presa diretta: $n_{GB} = 1$, $\eta_{GB} = 1$) |
| $\psi, \theta$ | $\eta_{\psi 0}, \eta_{\theta 0}$: angoli degli **assi medi** delle catene (§4) |

### 2.1 Incoerenze trovate nella tesi di De Donatis (da segnalare, non da correggere in silenzio)

1. **Eq. (3.15), termine $M_{11}$.** Ricavando $M_{11}$ dall'energia cinetica (3.24) con la cinematica (3.6)–(3.8) si ottiene

   $$M_{11} = I_{1,33} + I_{2,33} - 2I_{2,31}\sin\theta\cos\theta + \left(I_{2,11} + m_2\rho_2^2 - I_{2,33}\right)\sin^2\theta$$

   La forma stampata ha $-(\ldots)\cos^2\theta$ al posto di $+(\ldots)\sin^2\theta$. Le due forme differiscono della costante $I_{2,33} - I_{2,11} - m_2\rho_2^2$. Allo zenit ($\theta_{DD} = 0$) la forma stampata dà $2I_{2,33} - I_{2,11} - m_2\rho_2^2 + I_{1,33}$, mentre l'inerzia fisica è $I_{2,33} + I_{1,33}$. $M_{12}$, $M_{22}$, $c(\theta)$ (3.17), la matrice di Coriolis (3.16) e la gravità (3.18) sono invece corretti: non dipendono dalla costante. Con i parametri segnaposto, a θ0 = −30° l'errore su $M_{11}$ è del 34 %. È probabilmente solo un refuso di stampa (il codice MATLAB può essere corretto), ma va controllato prima di riusarne i numeri.
2. **§5.2.1, punto di linearizzazione.** $\dot{\mathbf q}_0 = [0.5,\ 0.5]$ rad/s non è un moto stazionario: con $\dot\theta_0 \ne 0$ l'angolo θ cambia nel tempo, quindi quella è una linearizzazione "a tempo congelato". Lo spec usa $\dot\theta_0 = 0$, $\dot\psi_0 = \Omega_0$, che è un vero moto stazionario: va mantenuto così.
3. Refusi minori: in §3.2.1 la seconda rotazione è "di θ", non "di ψ"; nella (5.3) il termine di smorzamento è $\mathbf C\dot{\mathbf q}$.

## 3. Passo 1 — modello strutturale per asse (Gawronski eq. 2.4)

Per ciascun asse $a \in \{\psi, \theta\}$ il tuo FEM ridotto con Guyan fornisce

$$\mathbf M_{r,a}\,\ddot{\mathbf q}_a + \mathbf C_{r,a}\,\dot{\mathbf q}_a + \mathbf K_{r,a}\,\mathbf q_a = \mathbf b_{o,a}\,\tau_a, \qquad y_a = \mathbf c_{oq,a}\,\mathbf q_a$$

con $\mathbf q_\psi = [\psi_e, \psi_m, \psi_L]$, $\mathbf q_\theta = [\theta_e, \theta_p, \theta_m, \theta_r]$, $\mathbf b_{o,a} = \mathbf e_{m}$ (nodo motore) e $\mathbf c_{oq,a} = \mathbf e_{e}^{\mathsf T}$ (encoder). È esattamente la (2.4) di Gawronski, con $N = 3$ e $N = 4$.

**Inerzia di riferimento del payload.** Per costruire i modi serve l'inerzia del payload ai nodi di aggancio, valutata in una configurazione di riferimento $\theta_{ref}$:

$$\bar{\mathbf M}_\psi = \mathbf M_{r,\psi} + \left[J^{ax} + J_{\psi,p}(\theta_{ref})\right]\mathbf e_L\mathbf e_L^{\mathsf T}, \qquad \bar{\mathbf M}_\theta = \mathbf M_{r,\theta} + I_{yy,p}\,\mathbf e_P\mathbf e_P^{\mathsf T}$$

con $J_{\psi,p}(\theta) = M_{p,11}(\theta)$ (§5). Questa è la costruzione `isolated=True` dello spec (§7.3). Nel modello accoppiato la stessa inerzia viene **sottratta** nel termine del payload (§5), quindi non c'è doppio conteggio.

## 4. Passo 2 — base modale con modo rigido (Gawronski §2.2.2 e §2.2.4)

Il problema agli autovalori $(\mathbf K_{r,a} - \omega^2\bar{\mathbf M}_a)\boldsymbol\phi = \mathbf 0$ (eq. 2.6) ha una sola radice nulla, perché $\mathbf K_{r,a}\mathbf 1 = \mathbf 0$ (catena libera, eq. 2.30). Si sceglie

$$\boldsymbol\Phi_a = \begin{bmatrix}\mathbf 1 & \boldsymbol\phi_{a1} & \cdots\end{bmatrix}, \qquad \boldsymbol\Phi_a^{\mathsf T}\bar{\mathbf M}_a\boldsymbol\Phi_a = \operatorname{diag}(J_{a,tot},\ 1,\ \ldots,\ 1), \qquad \boldsymbol\Phi_a^{\mathsf T}\mathbf K_{r,a}\boldsymbol\Phi_a = \operatorname{diag}(0,\ \omega_{a1}^2,\ \ldots)$$

con il modo rigido **non normalizzato**, così che $\eta_{a0}$ sia un angolo in radianti con massa modale $J_{a,tot} = \mathbf 1^{\mathsf T}\bar{\mathbf M}_a\mathbf 1$, e i modi flessibili normalizzati sulla massa. Con $\mathbf q_a = \boldsymbol\Phi_a\boldsymbol\eta_a$ (eq. 2.13):

- **Significato di $\eta_{a0}$.** Per ortogonalità, $\eta_{a0} = \mathbf 1^{\mathsf T}\bar{\mathbf M}_a\mathbf q_a / J_{a,tot}$: è la rotazione media pesata sull'inerzia, cioè l'angolo dell'**asse medio** della catena. Gli angoli ψ, θ di De Donatis vanno letti così, non come letture degli encoder: $\psi_e = \eta_{\psi 0} + \sum_k \phi_{\psi k}(e)\,\eta_{\psi k}$.
- **Vettori di Gawronski** (eq. 2.21–2.22), righe di $\boldsymbol\Phi_a$ nei nodi:
  - encoder: $\mathbf c_{e,a} = \boldsymbol\Phi_a^{\mathsf T}\mathbf e_e$, con $c_{e,a0} = 1$;
  - motore: $\mathbf c_{m,a} = \boldsymbol\Phi_a^{\mathsf T}\mathbf e_m$, da cui $\mathbf B_m = \mathbf M_m^{-1}\mathbf c_{m,a}$ e $b_{a0} = 1/J_{a,tot}$;
  - **aggancio del payload** (nuovi rispetto a Gawronski): $\mathbf c_L = \boldsymbol\Phi_\psi^{\mathsf T}\mathbf e_L$ e $\mathbf c_P = \boldsymbol\Phi_\theta^{\mathsf T}\mathbf e_P$, entrambi con prima componente pari a 1.
- **Smorzamento** (eq. 2.12, 2.16). Con $\mathbf C = \beta\mathbf K$ (spec Q8), $\boldsymbol\Phi^{\mathsf T}\mathbf C\boldsymbol\Phi = \operatorname{diag}(0,\ \beta\omega_k^2)$, quindi $\zeta_k = \beta\omega_k/2$ esattamente e il modo rigido non è smorzato, come richiesto dallo spec. L'attrito dei cuscinetti $\mathbf C_b$ non è proporzionale: $\boldsymbol\Phi^{\mathsf T}\mathbf C_b\boldsymbol\Phi$ è piena. Il suo elemento rigido-rigido è $\mathbf 1^{\mathsf T}\mathbf C_b\mathbf 1$, cioè il $n_{GB}^2 b_m$ di De Donatis.

Si raccolgono i due assi: $\boldsymbol\Phi = \operatorname{blkdiag}(\boldsymbol\Phi_\psi, \boldsymbol\Phi_\theta)$, $\boldsymbol\eta = [\eta_{\psi 0}, \eta_{\psi 1}, \eta_{\psi 2}, \eta_{\theta 0}, \ldots, \eta_{\theta 3}]$ (7 coordinate), e la matrice d'aggancio

$$\mathbf E_\eta = \begin{bmatrix}\mathbf c_L^{\mathsf T} & \mathbf 0\\ \mathbf 0 & \mathbf c_P^{\mathsf T}\end{bmatrix}\ (2\times 7), \qquad \begin{bmatrix}\psi_L\\ \theta_p\end{bmatrix} = \mathbf E_\eta\,\boldsymbol\eta$$

## 5. Passo 3 — il payload come corpo rigido di De Donatis agganciato a L e P

Il payload è un corpo rigido la cui energia dipende solo da $(\psi_L, \theta_p)$ e dalle loro derivate (spec §9.2). Nella convenzione dello spec (verifica A2):

$$\mathbf M_p(\theta) = \begin{bmatrix} I_{xx}\sin^2\theta - I_{xz}\sin 2\theta + I_{zz}\cos^2\theta & I_{yz}\cos\theta - I_{xy}\sin\theta\\ I_{yz}\cos\theta - I_{xy}\sin\theta & I_{yy}\end{bmatrix}$$

$$\mathbf C_p\,\dot{\mathbf q}_p = \begin{bmatrix} M_{11}'(\theta)\,\dot\psi\dot\theta + M_{12}'(\theta)\,\dot\theta^2\\ -\tfrac12 M_{11}'(\theta)\,\dot\psi^2\end{bmatrix}, \qquad \mathbf g_p = \begin{bmatrix}0\\ -m_p g\left(d_x\cos\theta + d_z\sin\theta\right)\end{bmatrix}$$

con $M_{11}' = (I_{xx} - I_{zz})\sin 2\theta - 2I_{xz}\cos 2\theta$ e $M_{12}' = -I_{xy}\cos\theta - I_{yz}\sin\theta$ (tutti i momenti per O, in terna B). Sono le (3.14)–(3.18) di De Donatis riscritte con $\theta_{DD} = \theta + \pi/2$ (§2).

Le coordinate d'aggancio sono funzioni lineari di $\boldsymbol\eta$, quindi le equazioni di Lagrange si trasformano in modo covariante: la forza generalizzata del payload sulle coordinate modali è $\mathbf E_\eta^{\mathsf T}$ per le "coppie di De Donatis" calcolate in $(\psi_L, \theta_p)$.

### 5.1 Equazioni del moto non lineari in coordinate modali (risultato principale)

$$\boxed{\;\tilde{\mathbf M}\,\ddot{\boldsymbol\eta} + \tilde{\mathbf D}\,\dot{\boldsymbol\eta} + \tilde{\mathbf K}\,\boldsymbol\eta + \mathbf E_\eta^{\mathsf T}\Big\{\big[\mathbf M_p(\theta_p) - \bar{\mathbf M}_p\big]\begin{bmatrix}\ddot\psi_L\\ \ddot\theta_p\end{bmatrix} + \mathbf C_p(\theta_p, \dot\psi_L, \dot\theta_p)\begin{bmatrix}\dot\psi_L\\ \dot\theta_p\end{bmatrix} + \mathbf g_p(\theta_p)\Big\} = \boldsymbol\Phi^{\mathsf T}\mathbf B\,\boldsymbol\tau\;}$$

- $\tilde{\mathbf M} = \operatorname{diag}(J_{\psi,tot}, 1, 1, J_{\theta,tot}, 1, 1, 1)$, $\tilde{\mathbf K} = \operatorname{diag}(0, \omega_{\psi 1}^2, \omega_{\psi 2}^2, 0, \omega_{\theta 1}^2, \ldots)$ e $\tilde{\mathbf D} = \operatorname{diag}(0, 2\zeta_k\omega_k) + \boldsymbol\Phi^{\mathsf T}\mathbf C_b\boldsymbol\Phi$ formano il modello modale di Gawronski (2.14) delle due catene, ciascuna con il payload "congelato" in $\theta_{ref}$;
- $\bar{\mathbf M}_p = \operatorname{diag}(J_{\psi,p}(\theta_{ref}), I_{yy,p})$ è l'inerzia del payload già contenuta nei modi;
- $\boldsymbol\Phi^{\mathsf T}\mathbf B\boldsymbol\tau = [\mathbf c_{m,\psi}\tau_\psi;\ \mathbf c_{m,\theta}\tau_\theta]$: la coppia del motore entra in ogni modo con la sua partecipazione modale al nodo motore.

**Lettura fisica.** Il payload agisce sulla struttura flessibile come due coppie applicate in L e in P, uguali e opposte alle coppie d'inerzia, giroscopiche e di gravità del corpo rigido di De Donatis (meno l'inerzia di riferimento). Ogni modo le riceve pesate con $c_{L,k}$ e $c_{P,k}$, esattamente come riceve la coppia del motore pesata con $c_{m,k}$. L'accoppiamento tra gli assi nasce **solo** qui, in accordo con lo spec §9.2: gli alberi assialsimmetrici non danno termini giroscopici.

**Esattezza.** Tenendo tutti i 7 modi, il cambio di coordinate $\mathbf q = \boldsymbol\Phi\boldsymbol\eta$ è invertibile e il sistema coincide con quello dello spec §9.5 (verifica B9: residuo relativo $4\cdot10^{-16}$ su 20 stati casuali). Non serve alcuna ipotesi di piccole deformazioni oltre a quelle del modello torsionale: le rotazioni attorno allo stesso asse si sommano esattamente, e l'unica non linearità geometrica è $\theta_p$ dentro $\mathbf M_p$.

## 6. Risultato 1 — limite rigido = modello di De Donatis

Si tengono solo $\eta_{\psi 0}$ e $\eta_{\theta 0}$ (troncamento di Rayleigh–Ritz sui modi rigidi). Poiché $c_{L,0} = c_{P,0} = c_{m,0} = 1$, si ha $\mathbf E_\eta \to \mathbb 1_2$, $\psi_L = \eta_{\psi 0}$, $\theta_p = \eta_{\theta 0}$ e

$$\Big[\operatorname{diag}(J_{\psi,s},\ J_{\theta,s}) + \mathbf M_p(\theta)\Big]\ddot{\mathbf q} + \Big[\mathbf C_p(\theta, \dot{\mathbf q}) + \operatorname{diag}(\mathbf 1^{\mathsf T}\mathbf C_{b,\psi}\mathbf 1,\ \mathbf 1^{\mathsf T}\mathbf C_{b,\theta}\mathbf 1)\Big]\dot{\mathbf q} + \mathbf g_p(\theta) = \boldsymbol\tau, \qquad \mathbf q = [\psi, \theta]$$

perché $\mathbf 1^{\mathsf T}\mathbf K_r\mathbf 1 = 0$ (le molle spariscono), $\beta\mathbf K\mathbf 1 = \mathbf 0$ (lo smorzamento proporzionale sparisce) e $J_{a,tot} - \bar M_{p,aa} = J_{a,s}$. Con la corrispondenza della §2 questa è **esattamente** la (4.7) di De Donatis (la (3.13) più rotore e attrito riportati).

Verifiche:

- simbolica (A5): identità di $\mathbf M$ e $\mathbf g$ per ogni θ e ogni tensore;
- numerica (B3): proiezione del modello linearizzato a 7 GdL sui modi rigidi contro il modello di De Donatis linearizzato con i parametri nativi della sua tesi. Scarti $6\cdot10^{-17}$ su $\mathbf M$, $1.5\cdot10^{-15}$ su $\mathbf C$, $10^{-10}$ su $\mathbf K$ (arrotondamento su termini di $10^{5}$ N·m/rad);
- limite di rigidezza (B7, test 16 dello spec): con tutte le rigidezze ×10⁴, FRF sotto 10 Hz entro $2.5\cdot10^{-4}$ e polo lento polare entro $4\cdot10^{-6}$ da De Donatis.

**Conseguenza per la tesi.** Il modello di De Donatis non è un modello "diverso" da validare a parte: è il primo livello del tuo, ed è valido finché i modi flessibili non contano. Il §8 mostra che il limite effettivo non è la prima risonanza ma la prima **antirisonanza**.

## 7. Risultato 2 — modello linearizzato nella forma di Gawronski

### 7.1 Linearizzazione del payload (verificata in A, simbolica, con prodotti d'inerzia generici)

Nel punto di lavoro $(\theta_0, \Omega_0)$, con $\dot\theta_0 = 0$:

$$\mathbf G_{p0} = \frac{\partial(\mathbf C_p\dot{\mathbf q}_p)}{\partial\dot{\mathbf q}_p} = \begin{bmatrix}0 & \gamma\\ -\gamma & 0\end{bmatrix}, \qquad \gamma = \Omega_0\,M_{11}'(\theta_0) = \Omega_0\left[(I_{xx} - I_{zz})\sin 2\theta_0 - 2I_{xz}\cos 2\theta_0\right]$$

$$k_c = -\tfrac12\,\Omega_0^2\,M_{11}''(\theta_0), \qquad k_g = m_p g\left(d_x\sin\theta_0 - d_z\cos\theta_0\right)$$

Con prodotti nulli si ritrovano i target dello spec §10.2. $\mathbf G_{p0}$ resta **puramente antisimmetrica** anche con i prodotti; $M_{12}$ entra solo nella massa.

### 7.2 Equazioni modali linearizzate

Con $\theta_{ref} = \theta_0$ e i vettori a 7 componenti $\hat{\mathbf c}_L = [\mathbf c_L;\ \mathbf 0]$, $\hat{\mathbf c}_P = [\mathbf 0;\ \mathbf c_P]$:

$$\Big[\tilde{\mathbf M} + M_{12}(\theta_0)\big(\hat{\mathbf c}_L\hat{\mathbf c}_P^{\mathsf T} + \hat{\mathbf c}_P\hat{\mathbf c}_L^{\mathsf T}\big)\Big]\delta\ddot{\boldsymbol\eta} + \Big[\tilde{\mathbf D} + \gamma\big(\hat{\mathbf c}_L\hat{\mathbf c}_P^{\mathsf T} - \hat{\mathbf c}_P\hat{\mathbf c}_L^{\mathsf T}\big)\Big]\delta\dot{\boldsymbol\eta} + \Big[\tilde{\mathbf K} + (k_c + k_g)\,\hat{\mathbf c}_P\hat{\mathbf c}_P^{\mathsf T}\Big]\delta\boldsymbol\eta = \boldsymbol\Phi^{\mathsf T}\mathbf B\,\delta\boldsymbol\tau$$

Il modello di Gawronski diagonale è quindi perturbato **solo da tre correzioni di rango basso**, tutte costruite con i vettori d'aggancio (verifica B4, scarti a livello di arrotondamento). Scritta per componenti, per un modo flessibile $k$ di azimut e un modo flessibile $j$ polare:

$$\ddot\eta_k + 2\zeta_k\omega_k\dot\eta_k + \omega_k^2\eta_k + c_{L,k}\sum_{j'}c_{P,j'}\left(M_{12}\,\ddot\eta_{j'} + \gamma\,\dot\eta_{j'}\right) = c_{m\psi,k}\,\tau_\psi$$

$$\ddot\eta_j + 2\zeta_j\omega_j\dot\eta_j + \omega_j^2\eta_j + c_{P,j}\sum_{k'}c_{L,k'}\left(M_{12}\,\ddot\eta_{k'} - \gamma\,\dot\eta_{k'}\right) + (k_c + k_g)\,c_{P,j}\sum_{j'}c_{P,j'}\,\eta_{j'} = c_{m\theta,j}\,\tau_\theta$$

Per i modi rigidi ($k = 0$, $j = 0$) si sostituisce 1 con la massa modale $J_{a,tot}$ davanti a $\ddot\eta$ e si pone $\omega = \zeta = 0$. Sono le equazioni di De Donatis linearizzate.

**Forma di stato.** Con $\mathbf x_i = [\omega_i\eta_i,\ \dot\eta_i]$ per i modi flessibili e $[\eta_0,\ \dot\eta_0]$ per i rigidi, i blocchi diagonali sono le (2.25) e (2.32) di Gawronski. Le tre correzioni producono blocchi fuori diagonale. Nel caso disaccoppiato (sotto) la verifica B5 ritrova la FRF fisica $\mathbf P(j\omega)^{-1}$ entro $6\cdot10^{-10}$.

**Variante con diagonalizzazione completa.** Se si include $(k_c + k_g)\,\mathbf e_P\mathbf e_P^{\mathsf T}$ in $\bar{\mathbf K}_\theta$ prima del calcolo dei modi, il "modo rigido" polare diventa un modo a frequenza $\omega_{\theta 0}^2 \approx (k_c + k_g)/J_{\theta,tot}$ (pendolo dello sbilanciamento più la rigidezza centrifuga). Con valore negativo il moto è instabile senza controllo. In questa variante il blocco (2.32) diventa $\begin{bmatrix}0 & 1\\ -\omega_{\theta 0}^2 & 0\end{bmatrix}$, $\mathbf M$ e $\mathbf K$ sono diagonali, e restano solo il termine giroscopico (antisimmetrico, che nessuna trasformazione modale reale elimina) e $M_{12}$. Nell'esempio segnaposto a θ0 = −30°, Ω0 = 10 °/s, quel modo è a 94 mHz.

### 7.3 Due conseguenze che ritrovano Gawronski §3.2.4

1. **Disaccoppiamento.** Con $\Omega_0 = 0$, prodotti nulli e payload bilanciato si ha $\gamma = M_{12} = k_c = k_g = 0$: il modello è esattamente **block-diagonal**, cioè due modelli SISO di Gawronski indipendenti, e le FdT incrociate sono identicamente nulle (B5, test 15 dello spec). È il motivo per cui Gawronski può identificare azimut ed elevazione separatamente ("very weak coupling between the axes").
2. **L'azimut dipende dall'elevazione, il polare no.** L'inerzia del payload vista dal polare è $I_{yy,p}$, costante; vista dall'azimut è $J_{\psi,p}(\theta_0)$, variabile. Con payload bilanciato, quindi, le FdT polari **non cambiano** con θ0 (scarto 0 nella verifica B8), mentre quelle di azimut sì. Nell'esempio, da −90° a 0°, l'antirisonanza di τψ → ψe passa da 111 a 91 Hz e il primo modo da 521.7 a 517.8 Hz. È ciò che Gawronski misura sulle antenne da 34 m (Fig. 3.9 contro Fig. 3.10); qui ne esce anche la causa. Con lo sbilanciamento ($k_g \ne 0$) anche il polare diventa dipendente da θ0. Le antenne DSN sono controbilanciate, ed è per questo che Gawronski non lo osserva.

Figura: `fig_modi_vs_theta0.png`.

## 8. Risultato 3 — gerarchia "De Donatis + n modi" e riduzione alla Gawronski

### 8.1 Criterio di scelta dei modi (Gawronski eq. 4.4–4.6, 4.11–4.12)

Per un'uscita di posizione la matrice d'uscita del modo $i$ nello stato (2.27) vale $c_{e,i}/\omega_i$, quindi

$$\|G_i\|_\infty \simeq \frac{\|\mathbf b_i\|\,\|\mathbf c_{e,i}\|}{2\zeta_i\omega_i^2}, \qquad \|G_i\|_h \simeq \tfrac12\|G_i\|_\infty, \qquad \mathbf b_i = \frac{\boldsymbol\Phi^{\mathsf T}\mathbf B}{\tilde M_{ii}}\Big|_i$$

I modi rigidi hanno norma infinita e si tengono **sempre** (eq. 4.11–4.12). Il modello minimo di ogni gerarchia è quindi De Donatis. I modi flessibili si aggiungono in ordine di norma di Hankel decrescente.

### 8.2 Esempio (parametri segnaposto, θ0 = −30°, Ω0 = 10 °/s)

| Ordine | Asse | $f$ [Hz] | $\zeta$ | $\|G_i\|_h$ |
|---|---|---|---|---|
| 1 | azimut | 518.4 | 0.010 | 7.9·10⁻⁵ |
| 2 | polare | 254.0 | 0.010 | 2.1·10⁻⁵ |
| 3 | azimut | 1417.0 | 0.027 | 4.2·10⁻⁶ |
| 4 | polare | 603.8 | 0.024 | 1.2·10⁻⁷ |
| 5 | polare | 1156.1 | 0.046 | 1.6·10⁻⁸ |

Banda di validità di ciascun livello, definita come la frequenza fino alla quale lo scarto dal modello completo sui canali diretti resta sotto 1 dB e sotto 5°:

| Modello | Troncamento modale puro | Con correzione statica |
|---|---|---|
| De Donatis (solo modi rigidi) | 32 Hz | 86 Hz |
| rigidi + 1 modo | 61 Hz | 147 Hz |
| rigidi + 2 modi | 61 Hz | 447 Hz |
| rigidi + 3 modi | 422 Hz | 447 Hz |
| rigidi + 4 modi | 513 Hz | 831 Hz |
| rigidi + 5 modi (= tuo modello a 7 GdL) | esatto | esatto |

Figura: `fig_gerarchia_bode.png` (Bode 2×2 coppia → encoder: De Donatis, rigidi + 2 modi con correzione, modello completo).

### 8.3 Cosa insegna l'esempio (qualitativamente robusto, numeri segnaposto)

1. **De Donatis vale ben sotto la prima risonanza.** La prima risonanza è a 254 Hz, ma il modello rigido si stacca già a circa 30 Hz. Il limite è l'**antirisonanza** di azimut (circa 95 Hz), che compare in τψ → ψe a $\sqrt{k_{\psi,2}/J_L^*}$. Nella presa diretta con $J_L^* \gg J_m$, l'antirisonanza precede la risonanza del fattore $\sqrt{1 + J_L^*/J_m}$, e si sposta con θ0 (§7.3). È la versione a presa diretta del criterio di Gawronski §2.5 ("la frequenza del drive deve superare la fondamentale della struttura"): nel suo caso, con $N = 354$, l'inerzia riflessa del motore $N^2 J_m$ domina e la situazione è rovesciata.
2. **Il troncamento modale conserva i poli ma sposta gli zeri.** Aggiungere modi senza correzione migliora poco: gli zeri dipendono dalla flessibilità residua dei modi scartati. La **correzione statica** (flessibilità residua con "inertia relief", che per una catena libera è l'unica definita) ripristina la cedevolezza tra motore ed encoder:

   $$\mathbf G_{n}(s) \;\to\; \mathbf G_{n}(s) + \sum_{i\ \text{scartati}}\frac{\mathbf c_{e,i}\,\mathbf c_{m,i}^{\mathsf T}}{\tilde M_{ii}\,\omega_i^2}$$

   Il termine aggiunto è una matrice costante $\mathbf D$ nello spazio degli stati. Gawronski non la usa perché lavora con la troncatura bilanciata (§4.2), che conserva il guadagno in modo diverso. Per una gerarchia fisicamente leggibile, il troncamento modale con correzione statica è la scelta consigliata.
3. L'ordine di Hankel e l'ordine in frequenza **non coincidono**: il polare a 603.8 Hz pesa meno dell'azimut a 1417 Hz, perché quest'ultimo ha grande partecipazione all'encoder ($c_e = 20.8$).

## 9. Drive: motore e riduttore (Gawronski §2.3 ↔ De Donatis §4.1–4.2)

- **Riduttore.** Il modello di Gawronski, $T = k_g(\theta_g - \theta_p)$ con $\theta_g = \theta_m/N$ (eq. 2.46–2.47), è **un elemento della catena**: un nodo motore con inerzia $N^2 J_m$ e coordinata $\theta_m/N$, collegato alla struttura da una molla $k_g$, con coppia $N\eta\,\tau_{motor}$. La (4.7) di De Donatis è questo con $k_g \to \infty$, proiettato sul modo rigido. Nella tua presa diretta $N = 1$, e "$k_g$" sono le molle ridotte tra nodo motore e carico ($k_{\psi,2}$, $k_{\theta,2}$). Il banco di De Donatis ($n_{GB} = 10$) si modella con la stessa catena, cambiando solo il nodo motore. Il gioco (15′, §4.3) resta fuori dal modello lineare.
- **Parte elettrica** (facoltativa, se l'ingresso è la tensione e non la coppia). Gawronski (2.43)–(2.44) e De Donatis (4.1) coincidono: $L_a\,\dot i + R_a\,i + k_b\,\omega_m = v$, $\tau = k_t\,i$. Nel modello modale $\omega_m = \mathbf c_m^{\mathsf T}\dot{\boldsymbol\eta}$: la forza controelettromotrice smorza ogni modo con peso $k_t k_b c_{m,k}^2/R_a$ ed è collocata con la coppia. Con il controllo in corrente previsto dallo spec non serve.

## 10. Implementazione in `gimbal_multidof` (proposta per Claude Code)

Nuovo modulo `src/gimbal_mdof/modal.py`, che importa da `reduction.py` e `coupled.py` senza duplicarli:

| Funzione | Contenuto | Equazioni |
|---|---|---|
| `reference_inertia(axis, theta_ref)` | $\bar{\mathbf M}_a$ con il payload congelato | §3 |
| `axis_modes(Mbar, Kr)` | $\boldsymbol\Phi_a$ con modo rigido $\mathbf 1$, frequenze, masse modali, $\zeta_k = \beta\omega_k/2$ | §4 |
| `attachment_vectors(Phi)` | $\mathbf c_e, \mathbf c_m, \mathbf c_L, \mathbf c_P$; `assert` sulla prima componente pari a 1 | §4 |
| `modal_eom_nonlinear(eta, eta_d, tau)` | residuo della (5.1), per `simulate.py` | §5.1 |
| `modal_linear(theta0, Omega0)` | $\tilde{\mathbf M}, \tilde{\mathbf D}, \tilde{\mathbf K}$ più le tre correzioni di rango basso; stato a blocchi (2.25)/(2.32) | §7 |
| `hankel_norms()`, `hierarchy(n, static_correction=True)` | norme modali, modelli "De Donatis + n modi", bande di validità | §8 |
| `dedonatis_map(I2, m2, rho2)` | conversione dei dati di De Donatis nella convenzione dello spec | §2 |

**Nuovi test** (in aggiunta alla sezione 13 dello spec):

- **U1.** $\mathbf K_r\mathbf 1 = \mathbf 0$; $\boldsymbol\Phi^{\mathsf T}\bar{\mathbf M}\boldsymbol\Phi$ e $\boldsymbol\Phi^{\mathsf T}\mathbf K\boldsymbol\Phi$ diagonali; $c_{e,0} = c_{m,0} = c_{L,0} = c_{P,0} = 1$.
- **U2.** Proiezione sui modi rigidi uguale al modello di De Donatis linearizzato, con la $M_{11}$ corretta (§2.1).
- **U3.** $\boldsymbol\Phi^{\mathsf T}(\mathbf M_0, \mathbf C_0, \mathbf K_0)\boldsymbol\Phi$ uguale alla forma della §7.2.
- **U4.** Residuo della (5.1) uguale a $\boldsymbol\Phi^{\mathsf T}$ per il residuo dello spec §9.5, su stati casuali.
- **U5.** Gerarchia completa uguale al modello a 7 GdL; con payload bilanciato e $\Omega_0 = 0$, FdT polari indipendenti da θ0.
- **U6.** `dedonatis_map` più spec uguale alle (3.14)–(3.18) con $\theta_{DD} = \theta + \pi/2$ (dalla parte A).

Uscite: tabella dei modi con norme di Hankel e partecipazioni, tabella delle bande di validità della gerarchia, Bode 2×2 della gerarchia e figura "modi contro θ0". Tutte con l'avviso "risultato da parametri segnaposto".

## 11. Riepilogo delle verifiche eseguite

| ID | Verifica | Esito |
|---|---|---|
| A1 | $\boldsymbol\omega_p\vert_B$ da Poisson uguale a spec §9.1 | identità |
| A3 | $\mathbf M_{DD}$ dalla (3.24) contro la (3.14)–(3.17) stampata | $M_{12}, M_{22}, \mathbf C, \mathbf g$ corretti; $M_{11}$ con refuso (§2.1) |
| A5 | De Donatis uguale al payload dello spec più il disco, con $\theta_{DD} = \theta + \pi/2$ | identità |
| A6 | $\mathbf G_{p0}, k_c, k_g$ con prodotti generici; target dello spec §10.2 | identità |
| B3 | Modi rigidi uguali a De Donatis linearizzato | scarto ≤ 10⁻¹⁰ (arrotondamento) |
| B4 | Struttura modale della §7.2 | scarto a livello di arrotondamento |
| B5 | Stato a blocchi di Gawronski uguale alla FRF fisica; incrociate nulle nel caso disaccoppiato | 6·10⁻¹⁰; 0 |
| B6 | Gerarchia completa uguale al modello a 7 GdL | 1·10⁻⁹ |
| B7 | Rigidezze ×10⁴ tendono a De Donatis | 2.5·10⁻⁴ (FRF), 4·10⁻⁶ (polo) |
| B8 | FdT polare indipendente da θ0 (bilanciato) | scarto 0 |
| B9 | EOM non lineari modali uguali allo spec §9.5 | relativo 4·10⁻¹⁶ |

## 12. Decisioni da prendere

1. **$\theta_{ref}$ della base modale.** Default: $\theta_{ref} = \theta_0$ per ogni punto di lavoro nel modello lineare (massima diagonalità); un unico $\theta_{ref}$ (per esempio −45°, come il "modello medio" di Gawronski §3.2.4) per le simulazioni non lineari, dove la (5.1) è esatta per qualunque scelta.
2. **Correzione statica nella gerarchia.** Default: attiva per i modelli lineari destinati al progetto del controllo, disattivata nelle simulazioni non lineari, dove il modello completo a 7 GdL costa poco.
3. **Variante del modo polare** (§7.2): $k_c + k_g$ come accoppiamento (default, mantiene il legame diretto con De Donatis) oppure inclusi nei modi (diagonalizzazione completa).
4. **Dati di De Donatis.** $J_m$, $b_m$, $L_a$ identificati (Tab. 4.1) si riferiscono al banco con riduttore 10:1 e motori diversi da quelli finali. Si possono usare solo per validare il modello sul banco (catena con $N = 10$), non come parametri del gimbal finale.
