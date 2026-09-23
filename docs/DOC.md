# titrages — documentation

## La méthode

Une solution contient des espèces faibles $i$, chacune décrite comme un polyacide
$\mathrm{H}_{N_i}\mathrm{A}_i$ et introduite sous sa forme $q_i$ — le nombre de
protons qu'elle a déjà cédés —, à la concentration $C_i$, et un excès
$b = [\text{base forte}] - [\text{acide fort}]$.

Avec $h = [\mathrm{H_3O^+}]$, $w = K_e/h$, $\delta = w - h$ et les constantes de
dissociation cumulées $k_j = K_1 \cdots K_j$, les fractions sont explicites :

$$x_j = \frac{k_j/h^j}{\sum_m k_m/h^m}, \qquad Y_1 = \sum_j j\,x_j, \qquad Y_2 = \sum_j j^2 x_j.$$

L'électroneutralité s'écrit $R(\mathrm{pH}) = 0$, où le **résidu**

$$R(\mathrm{pH}) = \sum_i C_i\,(Y_{1,i} - q_i) - b + \delta$$

est strictement croissant, de dérivée le pouvoir tampon

$$\beta = \frac{\mathrm{d}R}{\mathrm{d\,pH}} = \ln 10\,\Bigl(\sum_i C_i\,(Y_{2,i} - Y_{1,i}^2) + w + h\Bigr).$$

Les ions spectateurs n'apparaissent pas : une solution préparée est neutre, ils
portent la charge $q - z_0$ qui compense celle de la forme introduite, et la charge
$z_0$ de la forme la plus acide disparaît.

**Le titrage.** Verser $V$ de burette dans $V_0$ de bécher, c'est mélanger les deux
solutions ; l'électroneutralité du mélange, multipliée par $V_0 + V$, se sépare :

$$V_0\,R_\text{bécher}(\mathrm{pH}) + V\,R_\text{burette}(\mathrm{pH}) = 0
\quad\Longrightarrow\quad
V(\mathrm{pH}) = -V_0\,\frac{R_\text{bécher}}{R_\text{burette}}.$$

Pour un acide titré par une base forte on retrouve
$n = V C_B/(C_A V_0) = (Y_1 + \delta/C_A)/(1 - \delta/C_B)$. Sans dilution
(`dilution=False`), $\delta$ disparaît du dénominateur. La pente est analytique :

$$\frac{\mathrm{d}V}{\mathrm{d\,pH}} = V_0\,\frac{R_\text{b}\,\beta_\text{bu} - \beta_\text{b}\,R_\text{bu}}{R_\text{bu}^2}.$$

**Le sens inverse** — le pH d'une solution, le pH après $V$ — se trouve par
dichotomie sur $R$ ou sur $V_0 R_\text{bécher} + V R_\text{burette}$, toutes deux
croissantes : une racine unique, sans réglage, là où un Newton divergerait au saut.

**La précision.** $x_j$ se calcule en retranchant le plus grand exposant
($\log_{10} k_j + j\,\mathrm{pH}$) avant d'élever à la puissance : aucune fraction
n'est le quotient de deux nombres démesurés.

**La dilution.** La séparation $V_0 R_\text{b} + V R_\text{bu} = 0$ est
toute la correction de dilution. Tracer $n_0 = Y_1 + \delta/C_T$ comme une
courbe de titrage — ce que fait Kalka (2021), dont la méthode est tirée — revient à
négliger l'eau qui accompagne la base versée : pour deux solutions de même
concentration, le pH est surestimé de $\log 3 \approx 0{,}48$ à deux fois
l'équivalence. `tests/test_dilution.py` confronte la formule directe à un calcul qui
n'emprunte rien au paquet : concentrations diluées à la main, formes par les
formules des manuels, électroneutralité résolue par dichotomie ; les deux
s'accordent à 10⁻⁹ près sur dix configurations.

## Construire une solution

```pycon
>>> from titrages import Couple, Solution, solution
>>> round(solution({"NH3": 0.1}).pH(), 2)
11.12
>>> round(solution({"CH3COOH": 0.1, "CH3COO-": 0.1}).pH(), 2)
4.76
>>> mon_acide = Couple("mon acide", (3.2, 8.1))
>>> round(Solution().avec(mon_acide.forme(1), 0.05).pH(), 2)
5.65

```

| Méthode | Rôle |
| --- | --- |
| `solution({nom: C, …}, pKe=14)` | par les noms de la base ; accepte aussi `Espece` et `Couple` |
| `Solution().avec(espece, C)` | ajoute une espèce (un `Couple` seul : sa forme acide) |
| `.avec_acide_fort(C)`, `.avec_base_forte(C)` | l'acide et la base forts |
| `.diluee(facteur)`, `.melangee(autre, V, V_autre)` | nouvelles solutions |
| `.pH()` | par dichotomie |
| `.concentrations(pH=None)` | chaque forme, H3O+ et HO- |
| `.pouvoir_tampon(pH=None)` | $\beta$, en mol/L |
| `.residu(pH)`, `.charge(pH)` | $R$, et $R - \delta$ |

## Un titrage

```pycon
>>> from titrages import Titrage, solution
>>> t = Titrage(solution({"H3PO4": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
>>> t.volumes_equivalents()
(10.0, 20.0, 30.0)
>>> [round(s.volume, 2) for s in t.sauts(35)]
[10.0, 20.0]
>>> round(t.tangentes(7, 18).volume, 2)
10.0

```

| Méthode | Rend |
| --- | --- |
| `volume(pH)` | le volume versé, `nan` hors du domaine |
| `derivee_volume(pH)`, `pente(pH)` | $\mathrm{d}V/\mathrm{d\,pH}$ et son inverse |
| `pH(V)` | par dichotomie, vectorisé |
| `courbe(volume_max, points=400)` | une `Courbe` : `volume`, `pH`, `pente`, points exacts à longueur d'arc constante |
| `sauts(volume_max)` | les points d'inflexion (`Saut`) |
| `volumes_equivalents()` | les équivalences stœchiométriques, titrant fort |
| `tangentes(volume_avant, volume_max)` | la construction des tangentes parallèles (`Tangentes`) |

`volumes_equivalents` compte chaque proton échangeable du bécher, rangé par pKa,
après l'acide ou la base forts que le bécher contient déjà. Pour un titrant faible,
l'équivalence n'a pas de définition aussi simple : `sauts` reste disponible.

## La conductimétrie

La loi de Kohlrausch, $\sigma = \sum_i \lambda_i^\circ c_i$ (λ° en mS·m²·mol⁻¹, c en
mol/L, σ en S/m), sur toutes les formes chargées, H₃O⁺, HO⁻ et les **ions
spectateurs**. `solution` les ajoute d'après ce qu'on nomme : le Cl⁻ de HCl, le Na⁺
de la soude, deux Na⁺ par carbonate de sodium ; un ion désigné par sa formule reçoit
Na⁺ (anion) ou Cl⁻ (cation), à changer par `contre_ions=`. Une forme chargée dont
λ° est inconnue est une erreur, pas un zéro.

```pycon
>>> from titrages import Titrage, solution
>>> t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), 10)
>>> round(float(t.conductivite(15)), 4)
0.8614

```

`Titrage.conductivite(V)` pondère la part ionique du bécher et de la burette par les
volumes, au pH du mélange ; `courbe_conductimetrique(volume_max, points, corrigee)`
donne σ(V) à pas constant, multipliée par $(V_0+V)/V_0$ si `corrigee`.

## L'oxydoréduction

Les électrons remplacent les protons, $E/s$ le pH ($s = RT\ln 10/F$). Les formes
d'un `CoupleRedox` vont de la plus réduite à la plus oxydée ; chaque étape a son E°,
son nombre d'électrons `n` et le nombre de protons `H` que consomme la réduction, pris
à pH fixé ($E^{\circ\prime} = E^\circ - s\,m/n\,\mathrm{pH}$). Le bilan d'électrons
$R(E) = \sum_i T_i (Y_i(E) - e_{q_i})$ est croissant.

Un couple **dimérique** (`stoechiometrie=(1, 2)` : I₂/I⁻, Cr₂O₇²⁻/Cr³⁺ ;
`(2, 1)` : Hg²⁺/Hg₂²⁺) n'a pas de $Y(E)$ : sa loi de Nernst porte un carré, et une
dilution change son potentiel. `TitrageRedox.volume(E)` reste explicite avec un seul
dimère — avec $W = (V_0 + V)\,u$, $u$ la forme monomère, le bilan est linéaire et la
conservation de l'élément est du second degré en $V$ — et se résout par dichotomie
au-delà. La courbe se calcule par $E(V)$, par dichotomie sur le bilan du mélange.

```pycon
>>> from titrages import TitrageRedox, solution_redox
>>> t = TitrageRedox(solution_redox({"I2": 0.05}), solution_redox({"Na2S2O3": 0.1}), 10, "ECS")
>>> t.volumes_equivalents()
(10.0,)

```

| Méthode | Rend |
| --- | --- |
| `E(V)` | le potentiel mesuré (par rapport à `reference`) |
| `volume(E)` | le volume qui y amène |
| `pente(V)` | $\mathrm{d}E/\mathrm{d}V$, par dérivation implicite |
| `courbe`, `sauts`, `volumes_equivalents` | comme en acide-base |
| `volume_initial(volume_max)` | 0, ou $V_\text{max}/1000$ si le bécher ne contient qu'une forme par couple, dont le potentiel n'est pas défini |

## Fonctions d'équilibre

`fractions(couple_ou_pKa, pH)`, `moments(…, pH, ordre)`, `log_k`, `delta(pH, pKe)`,
et `concentration_pour_pH(espece, pH)` — la concentration d'une espèce seule qui
donne ce pH, $C = -\delta/(Y_1 - q)$, le sens direct de la courbe pH = f(pC).

## La base

```pycon
>>> from titrages import couple, espece, indicateur
>>> espece("carbonate de sodium").q, couple("NH3").pKa
(2, (9.25,))
>>> indicateur("BBT").zone
(6.0, 7.6)

```

`couples()`, `indicateurs()`, `chercher(texte)`. Les données sont dans
`src/titrages/donnees/*.toml`. Casse, accents, espaces et écriture des charges ne
comptent pas ; une clé ambiguë dans la base est une erreur au chargement.

## Le tracé

`titrages.trace` (matplotlib, `pip install titrages[trace]`) : `tracer_titrage`,
`tracer_derivee`, `tracer_sauts`, `tracer_tangentes`, `tracer_fractions`,
`tracer_zone_virage`. Chacune dessine sur l'`Axes` donné, ou l'axe courant, et le
rend.

## L'export

`vers_pgfplots(courbe)` et `vers_csv(courbe, chemin)`.

## Le paquet LaTeX

Voir `latex/titrages-doc.pdf`. Son moteur Lua implémente le même algorithme ;
`tests/test_lua.py` confronte les deux.
