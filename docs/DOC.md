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

$$\frac{\mathrm{d}V}{\mathrm{d\,pH}} = -\frac{V_0\,\beta_\text{b} + V\,\beta_\text{bu}}{R_\text{bu}},$$

obtenue en dérivant $V_0 R_\text{b} + V R_\text{bu} = 0$ : le numérateur est une
somme de termes positifs, quand la forme $V_0(R_\text{b}\beta_\text{bu} -
\beta_\text{b}R_\text{bu})/R_\text{bu}^2$ soustrairait deux nombres presque égaux au
cœur du saut.

**Le sens inverse** — le pH d'une solution, le pH après $V$ — se trouve par
dichotomie sur $R$ ou sur $V_0 R_\text{bécher} + V R_\text{burette}$, toutes deux
croissantes : une racine unique, sans réglage, là où un Newton divergerait au saut.

**La précision.** $x_j$ se calcule en retranchant le plus grand exposant
($\log_{10} k_j + j\,\mathrm{pH}$) avant d'élever à la puissance : aucune fraction
n'est le quotient de deux nombres démesurés. Et rien ne soustrait deux nombres
voisins : la variance s'écrit $\sum_{j<k} x_jx_k(k-j)^2$ et non $Y_2 - Y_1^2$, qui
perd six chiffres quand une forme domine à $1 - 10^{-10}$ près — c'est-à-dire dans
un saut —, et $Y_1 - q$ s'écrit $\sum_j x_j(j-q)$.

**Les points d'inflexion**, ce qu'on mesure, se calculent à la précision de la
machine. Ce sont les minima de $|\mathrm{d}V/\mathrm{d\,pH}|$ (ou
$|\mathrm{d}V/\mathrm{d}E|$), cherchés sur la courbe directe — en $V$, le saut de
Fe²⁺ par MnO₄⁻ tient en $10^{-9}$ mL — et affinés non par une section dorée, qui
compare des valeurs et plafonne à $\sqrt\varepsilon$, mais par le zéro de
$(\ln|\mathrm{d}V/\mathrm{d\,pH}|)'$, trouvé par dichotomie
(`titrages.equilibre.minimum_lisse`). Les tests les confrontent à `mpmath` à 60
chiffres — `tests/reference_mpmath.py`, qui n'emprunte rien au paquet, régénère
les valeurs de référence — : $10^{-10}$ unité de pH, $10^{-11}$ V, couples
dimériques compris.

**La dilution.** La séparation $V_0 R_\text{b} + V R_\text{bu} = 0$ est
toute la correction de dilution. Tracer $n_0 = Y_1 + \delta/C_T$ comme une
courbe de titrage — ce que fait Kalka (2021), dont la méthode est tirée — revient à
négliger l'eau qui accompagne la base versée : pour deux solutions de même
concentration, le pH est surestimé de $\log 3 \approx 0{,}48$ à deux fois
l'équivalence. `tests/test_dilution.py` confronte la formule directe à un calcul qui
n'emprunte rien au paquet : concentrations diluées à la main, formes par les
formules des manuels, électroneutralité résolue par dichotomie ; les deux
s'accordent à 10⁻⁹ près sur dix configurations.

## Deux règles

**Un `float` pour un scalaire, un tableau pour un tableau.** Toute grandeur
calculée à un pH, un volume ou un potentiel scalaire rend un `float` ; à un
tableau (liste, tableau numpy), un tableau numpy de même forme
(`titrages.equilibre.sortie`). Les fractions font exception : à un pH scalaire elles
sont déjà un tableau, une valeur par forme.

**La conductivité est en S/m**, partout : λ° en mS·m²·mol⁻¹ et $c$ en mol/L la
donnent directement. `tracer_conductimetrie(…, unite="mS/cm")` et la clé LaTeX
`unite` convertissent pour l'affichage seulement.

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
| `.diluee(facteur)`, `.melangee(autre, V, V_autre)` | nouvelles solutions ; `SolutionRedox` a les mêmes |
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

`volumes_equivalents` fait d'abord réagir l'acide ou la base forts que le bécher
contient déjà avec ses espèces faibles — un excès de soude transforme l'acide
éthanoïque en éthanoate, qu'un acide fort titre ensuite — puis compte chaque site
restant, rangé par pKa, ce qui reste du réactif fort passant en premier :
`solution({"CH3COOH": 0.05, "NaOH": 0.1})` par HCl 0,1 donne `(5, 10)`. Deux sites
de même pKa ne font qu'une équivalence. Pour un titrant faible, l'équivalence n'a pas
de définition aussi simple : `sauts` reste disponible.

## La conductimétrie

La loi de Kohlrausch, $\sigma = \sum_i \lambda_i^\circ c_i$ (λ° en mS·m²·mol⁻¹, c en
mol/L, σ en S/m), sur toutes les formes chargées, H₃O⁺, HO⁻ et les **ions
spectateurs**. `solution` les ajoute d'après ce qu'on nomme : le Cl⁻ de HCl, le Na⁺
de la soude, deux Na⁺ par carbonate de sodium ; un ion désigné par sa formule reçoit
Na⁺ (anion) ou Cl⁻ (cation), à changer par `contre_ions=`. Une forme chargée dont
λ° est inconnue est une erreur, pas un zéro. La base couvre 34 couples sur 35
(`src/titrages/donnees/conductivites.toml`, une source par ligne : le CRC Handbook,
ou l'ion voisin qui tient lieu de valeur approchée, marqué ≈ ; l'EDTA n'a pas de
λ°). La charge d'une forme se lit dans sa formule (`charge`) : `CO3^2-`, `Fe3+`, et
les zwitterions, `H3N+CH2COO-` étant neutre.

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
| `pente(V)` | $\mathrm{d}E/\mathrm{d}V = -g_V/g_E$, dérivées partielles du bilan analytiques, dimères compris |
| `derivee_volume(E)` | $\mathrm{d}V/\mathrm{d}E$, la dérivée de la courbe directe |
| `courbe`, `sauts`, `volumes_equivalents`, `tangentes` | comme en acide-base ; `tangentes` rend une `TangentesRedox` (`E_avant`, `E_apres`, `E`) |
| `volume_initial(volume_max)` | 0, ou $V_\text{max}/1000$ si le bécher ne contient qu'une forme par couple, dont le potentiel n'est pas défini |

## La complexométrie

Le paramètre est $\mathrm{pL} = -\log \ell'$, le ligand libre sous toutes ses formes
protonées, à pH tamponné. Un `Complexe` donne les $\log\beta_j$ **cumulés** d'un métal
(complexes successifs $\mathrm{ML}_j$) ; un `Ligand` ses pKa, d'où les constantes
conditionnelles $\beta'_j = \beta_j\,\alpha^j$, $\alpha$ la fraction de la forme la plus
basique. Le bilan du ligand $R(\mathrm{pL}) = \sum_i C_i\,Y_{1,i} + \ell - C_L$ est
linéaire dans les totaux d'une solution : le volume est explicite,
$V = -V_0\,R_\text{bécher}/R_\text{burette}$, plusieurs métaux en compétition compris.
pM s'en déduit au volume obtenu.

```pycon
>>> from titrages import Complexe, Ligand, SolutionComplexation, TitrageComplexometrique
>>> edta = Ligand("Y4-", (0.0, 1.5, 2.0, 2.66, 6.16, 10.24))
>>> t = TitrageComplexometrique(
...     SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01), 20,
...     [Complexe("Ca2+", (10.7,))], ligand=edta, pH=10,
... )
>>> round(t.pM(20), 2), [round(s.volume, 3) for s in t.sauts(40)], t.volumes_equivalents()
(6.28, [20.0], (20.0,))

```

| Méthode | Rend |
| --- | --- |
| `pM(V, metal)` | pM du métal suivi (`metal`, par défaut celui du titrage, ou le seul) |
| `pL(V)`, `volume(pL)` | les deux sens du paramétrage, par dichotomie dans le premier |
| `pente(pL, metal)` | $\mathrm{d}\mathrm{pM}/\mathrm{d}V$, analytique |
| `derivee_volume(pL)` | $\mathrm{d}V/\mathrm{d}\mathrm{pL}$ |
| `courbe(volume_max, points, metal)` | une `CourbeComplexometrique` (`volume`, `pM`, `pL`, `pente`), à abscisse curviligne |
| `sauts(volume_max, points, metal)` | des `SautComplexometrique`, maxima de la pente affinés |
| `volumes_equivalents()` | ligand titrant : $n\,C_M V_0/C_L$ cumulés, métaux rangés par $\log\beta_n/n$ décroissant ; métal titrant : $V_0 C_L/(n\,C_M)$ |

`fractions_complexes(log_beta, pL)` donne les fractions de $\mathrm{M}, \mathrm{ML},
\dots, \mathrm{ML}_n$.

## La précipitation

Un `Precipite` $\mathrm{P}_\nu\mathrm{X}_\mu$ lie l'**ion commun** P (`ion_commun`,
Ag⁺ par défaut) à un partenaire X, de pKs donné, X pouvant être basique (`pKa`, pH
alors obligatoire). Le paramètre est pP. Pour un ensemble S de solides présents, le
bilan $R_S = a + \sum_{i\in S} r_i\,(C_{X_i} - x_{\text{sat},i}) - C_P$, $r = \nu/\mu$,
reste linéaire dans les totaux : chaque ensemble a son volume explicite, et l'on retient
celui qui est cohérent — chaque solide présent a une quantité positive, chaque absent
n'est pas sursaturé. La courbe a des **points anguleux** là où un solide apparaît ou
disparaît : ce sont des `Evenement`, distincts des sauts, qui ne se cherchent que dans
les morceaux lisses.

```pycon
>>> from titrages import Precipite, SolutionPrecipitation, TitrageParPrecipitation
>>> t = TitrageParPrecipitation(
...     SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
...     SolutionPrecipitation(commun=0.05), 20,
...     [Precipite("Cl-", 9.75), Precipite("CrO4^2-", 11.95, nu=2, pKa=(6.5,))], pH=7.5,
... )
>>> [(e.partenaire, round(e.volume, 3)) for e in t.evenements(30)]
[('Cl-', 0.0), ('CrO4^2-', 20.011)]

```

La méthode de Mohr : le chromate d'argent apparaît 11 µL après l'équivalence.

| Méthode | Rend |
| --- | --- |
| `p(V)`, `volume(pP)` | pP et son inverse, par dichotomie sur $\log a$ dans le premier |
| `p_partenaire(V, partenaire)` | pX du partenaire |
| `solides(pP)` | les partenaires dont le solide est présent |
| `pente(pP)`, `derivee_volume(pP)` | $\mathrm{d}\mathrm{pP}/\mathrm{d}V$ et son inverse, sur l'ensemble de solides du point |
| `courbe(volume_max, points)` | une `CourbePrecipitation` (`volume`, `p`, `pente`, `solides`) |
| `evenements(volume_max, points)` | les apparitions et disparitions de solides, affinées |
| `sauts(volume_max, points)` | des `SautPrecipitation`, dans les morceaux lisses |
| `volumes_equivalents()` | ion commun titrant : les partenaires dans l'ordre où ils précipitent ; partenaire titrant : l'ion commun restant (Charpentier-Volhard) |

Le pH n'est pas couplé : il est tamponné, et la protonation du ligand ou du partenaire
n'entre que par $\alpha$. Un pH libre demanderait un solveur à plusieurs inconnues.

## Fonctions d'équilibre

`fractions(couple_ou_pKa, pH)`, `moments(…, pH, ordre)`, `variance(…, pH)` et
`ecart_moyen(…, pH, q)` — $Y_2 - Y_1^2$ et $Y_1 - q$ écrits sans soustraction —,
`minimum_lisse(f, a, b, pas)`, `log_k`, `delta(pH, pKe)`,
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

`titrages.trace` (matplotlib, `pip install "titrages[trace]"`) : `tracer_titrage`
et `tracer_titrage_redox` (les courbes), `tracer_conductimetrie`, `tracer_derivee`,
`tracer_sauts` et `tracer_tangentes` (ces trois-là acceptent un `Titrage` comme un
`TitrageRedox`), `tracer_fractions`, `tracer_zone_virage`. Chacune dessine sur
l'`Axes` donné, ou l'axe courant, et le rend.

## L'export

`vers_pgfplots(courbe)`, pour les trois sortes de courbes (pH, σ, E) : la liste
`(V,y)` d'un `\addplot coordinates`. Pour un fichier, `numpy.savetxt` sur les
colonnes de la courbe suffit.

## Le paquet LaTeX

Voir `latex/titrages-doc.pdf`. Son moteur Lua implémente le même algorithme ;
`tests/test_lua.py` confronte les deux.
