# Projet : complexation et précipitation

État au 18 octobre 2026 : le **niveau 1 est fait** — Python, Lua, LaTeX, tracé,
base (voir *Avancement*). Ce document disait, le 3 octobre, ce que coûterait
l'ajout des équilibres de complexation et de précipitation, en deux niveaux, et
ce qu'on recommandait de faire d'abord. L'étude du niveau 2, le pH non tamponné,
est à la fin (*Niveau 2 : étude*).

## Ce qui rend le paquet simple aujourd'hui

Tout repose sur une propriété : avec **une seule inconnue chimique** — h en
acide-base, E en oxydoréduction —, chaque bilan est **linéaire en**
x = V/(V₀+V). À pH fixé, le volume est donc explicite :

```text
V(pH) = −V₀ · R_bécher(pH) / R_burette(pH)
```

C'est ce qui donne la dérivée exacte, les points d'inflexion à la précision de
la machine, la résolution du pH à volume donné par simple dichotomie, et le
moteur pgf, qui ne résout rien.

Avec un métal M et un ligand L, il y a **trois composants** (H, M, L), donc
trois concentrations libres, plus le paramètre de dilution x. En fixant l'une
des concentrations libres comme paramètre de la courbe, il reste en général un
système **non linéaire à deux inconnues** : V cesse d'être explicite.

## Niveau 1 — ce qui reste explicite

La parade est de **paramétrer par la grandeur qu'on mesure** : pM ou pAg, au
lieu du pH. Les cas des TP et des concours y entrent tous.

### Complexométrie en milieu tamponné

L'EDTA, ou tout ligand qui forme un complexe 1:1. À pH fixé, la protonation du
ligand (et une éventuelle réaction parasite du métal) se résume à une constante
conditionnelle β′. À pM donné :

- le bilan sur M, `[M] + β′[M][L′] = T_M(x)`, donne `[L′]` ;
- le bilan sur L, `[L′](1 + β′[M]) = T_L(x)`, est alors **linéaire en x**.

V(pM) est explicite, exactement comme V(pH) aujourd'hui.

### Précipitation d'un sel

Argentimétrie : AgCl, Ag₂CrO₄ (Mohr), Charpentier-Volhard. Tant que le solide
est présent, `[Cl⁻] = Ks/[Ag⁺]`. En notant s la quantité de solide par litre :

- bilan sur Ag : `[Ag⁺] + s = T_Ag(x)` ;
- bilan sur Cl : `[Cl⁻] + s = T_Cl(x)`.

Leur différence élimine s et reste linéaire en x : V(pAg) est explicite. Avant
l'apparition du solide, la solution est homogène et le calcul est trivial.

### Plusieurs solides

Cl⁻ et I⁻ titrés par Ag⁺. À pAg fixé, chaque concentration libre d'anion est
donnée par son Ks, et il reste un système **linéaire** en (x, s₁, s₂). Le solide
présent change d'un segment de courbe à l'autre ; chaque segment se borne là où
une quantité s s'annule ou bien où un produit ionique atteint son Ks.

### Ce qu'il faut écrire

- une famille « titrage à paramètre pM » à côté de l'acide-base et du rédox ;
- les segments de précipitation et la détection de leurs bornes : apparition ou
  disparition d'un solide, qui fait un **point anguleux** sur la courbe ;
- la recherche des points d'inflexion adaptée, pour qu'elle ne confonde pas un
  point anguleux avec un saut ;
- le portage en Lua, les clés LaTeX, les tests contre un calcul indépendant
  (dilution faite à la main, précision arbitraire) et la documentation, Python
  et LaTeX.

Ordre de grandeur : celui de la refonte rédox de septembre 2026, soit **800 à
1 200 lignes** en tout. Le moteur pgf peut tracer ces courbes, puisque V reste
explicite.

## Niveau 2 — le cas général, qui demande un solveur

*Estimation du 3 octobre, avant l'étude ; l'étude du 18 octobre, à la fin du
document, trouve une voie plus simple pour les cas de classe.*

Ce qui sort du niveau 1 :

- un ligand basique en milieu **non tamponné** — NH₃ titré par Cu²⁺, le pH
  bougeant au fil du titrage ;
- des complexes successifs ML₁ … ML₄ couplés à l'acide-base ;
- la précipitation d'un hydroxyde ou d'un sel d'acide faible — Mg(OH)₂, CaCO₃ —,
  où le pH et le solide sont liés.

Il faut alors la méthode standard des logiciels de spéciation (méthode du
tableau de Morel, comme MINEQL) :

- des inconnues en **logarithme des concentrations libres**, résolues par Newton
  à chaque point ;
- une **continuation** le long de la courbe, qui est une variété de dimension 1
  dans l'espace (log h, log m, log l, x) : avancer en abscisse curviligne,
  prédire, corriger ;
- la dérivée exacte par le **théorème des fonctions implicites** — un système
  linéaire à résoudre par point ;
- pour les précipités, un **problème de complémentarité** : chaque solide est
  présent ou non, sa présence change en cours de route, et il faut localiser ces
  changements.

Le paquet n'y garde son caractère « sans approximation » que dans la mesure où
le solveur converge à la précision voulue, ce qu'il faut tester sérieusement.

Ordre de grandeur : **2 500 à 4 000 lignes**, dont un moteur Lua aussi gros que
le moteur Python. Le moteur pgf en est exclu.

## Les données

Les deux niveaux demandent une base de constantes — β des complexes, Ks des
solides —, avec leurs sources, à côté de la base de couples et d'indicateurs
actuelle. Ces valeurs se prennent dans une référence, jamais de mémoire, et
c'est souvent le poste le plus long.

## Comparaison

| | Niveau 1 | Niveau 2 |
| --- | --- | --- |
| Cas couverts | EDTA en milieu tamponné, Mohr, Charpentier-Volhard, mélanges d'halogénures | en plus : ligands basiques non tamponnés, complexes successifs, hydroxydes |
| Méthode | V explicite en pM | Newton et continuation |
| Précision | celle d'aujourd'hui | celle du solveur, à établir |
| Moteur pgf | oui | non |
| Taille estimée | 800 à 1 200 lignes | 2 500 à 4 000 lignes |

## Recommandation

Faire d'abord le **niveau 1**. Il couvre ce qu'on fait en classe, garde la
méthode explicite et la précision actuelle, et tient en une session de travail
du même ordre que le rédox. Le niveau 2 ne vaut la peine que pour des systèmes
couplés non tamponnés, qui sortent largement du programme.

## Avancement

### Fait le 18 octobre 2026 : le cœur Python du niveau 1

- `complexation.py` : `TitrageComplexometrique`, paramétré par **pL**, le
  ligand libre sous toutes ses formes protonées, et non par pM comme prévu
  ci-dessus. Le bilan du ligand reste linéaire dans les totaux d'une solution
  avec **plusieurs métaux** et des **complexes successifs** ML₁ … MLₙ : le
  paramétrage par pL couvre donc, à pH tamponné, deux cas que l'étude rangeait
  au niveau 2 (Cu–NH₃ tamponné, Ca et Zn en compétition).
- `precipitation.py` : `TitrageParPrecipitation`, paramétré par pP, l'ion
  commun. Les 2ⁿ ensembles de solides présents sont énumérés, chacun avec son
  volume explicite ; on retient l'ensemble cohérent. Apparitions et
  disparitions de solides (`evenements`) sont des points anguleux, distincts
  des sauts, qui ne se cherchent que dans les morceaux lisses. Mohr (chromate
  basique), Charpentier-Volhard, I⁻ puis Cl⁻.
- Tests contre une référence indépendante en multiprécision (dilution à la
  main, bilans résolus par dichotomie) : accord à 10⁻⁹ ou mieux.
- Documentation : deux sections de `docs/DOC.md`.

### Fait le 18 octobre 2026 : Lua, LaTeX, tracé, base

- Le moteur Lua porte les deux familles (`titrage_complexation`,
  `titrage_precipitation` et leurs fonctions), confronté au moteur Python par
  `tests/test_lua.py` à 10⁻⁹ près.
- Le paquet LaTeX : `\titragecomplexationdefinir`, `\titrageprecipitationdefinir`,
  `\addplottitrageevenements`, `\titragepM`, `\titrageevenement`, moteur lua
  seulement ; `latex/tests/titrages-test-cp.tex`, comparé à Python par
  `tests/test_latex.py`. Documentation l3doc : sections 8 et 9.
- Le tracé matplotlib (`tracer_titrage_complexometrique`,
  `tracer_titrage_precipitation`) et l'export pgfplots.
- La base : 13 ligands, 92 complexes, 41 précipités, des tables de Harvey
  (*Analytical Chemistry 2.1*, § 16.10 à 16.12, d'après Martell et Smith), recopiées
  telles qu'imprimées, avec leurs conditions. Une coquille de la source (Ksp de PbCl₂)
  est signalée dans la note.
- Les équivalences de complexométrie se comptent comme celles d'un polyacide, une par
  étape : les complexes qui prédominent tour à tour, sommets de l'enveloppe concave
  des log β. Dans le dosage de Liebig, Ag(CN)₃²⁻ puis Ag(CN)₂⁻.

### Fait le 18 octobre 2026 : niveau 2, étape 1 (Python)

- `complexation_libre.py` : le moteur à pH calculé, dichotomies imbriquées à volume
  fixé et Newton à jacobienne exacte pour affiner autour d'un point connu (ajout de
  points à la courbe, recherche des sauts). Près d'une équivalence, le bilan du
  ligand ne fixe log l qu'à quelques 10⁻¹³ près : le critère d'arrêt de Newton s'y
  adapte.
- Interface : sans `pH=`, le pH se calcule, avec un `UserWarning` ;
  `protons_ligand` et `acidobasique` sur `SolutionComplexation`.
- Tests : la référence mpmath écrit l'électroneutralité avec les contre-ions, et non
  le bilan des protons ; accord à 10⁻⁹ sur le pH et le pM, cinq cas (Ca par Na₂H₂Y
  avec et sans tampon, Zn avec acide fort, Cu par NH₃, Liebig avec HCN). Les deux
  limites : ligand sans acidité (même pM qu'au niveau 1 à 10⁻¹²) et tampon de plus
  en plus concentré (l'écart au pH imposé décroît d'un facteur 5 au moins par
  décade).

### Fait le 18 octobre 2026 : niveau 2, étape 2 (Python)

- `libre.py` : ce que partagent les deux moteurs à pH calculé — l'état à volume fixé
  (pH, u, ensemble S des solides), Newton à S fixé avec reprise par la dichotomie, la
  courbe raffinée, les sauts dans les morceaux lisses.
- `precipitation_libre.py` : l'ion commun et le pH, partenaires éliminés, hydroxydes
  comme bornes de l'ion commun dont la quantité libère des protons ; événements
  affinés par dichotomie en volume sur le prolongement du morceau qui précède.
- Quinze hydroxydes et oxydes dans la base (Harvey, § 16.10).
- Tests : référence mpmath par l'électroneutralité, six cas (Mohr sans tampon et en
  milieu nitrique, Charpentier-Volhard en milieu nitrique, Mg(OH)₂ avec et sans acide,
  Ag₂O), accord à 10⁻⁹ ; les deux limites (partenaire sans acidité, tampon
  concentré).

### Décisions

- Pas de base côté LaTeX (18 octobre 2026) : les constantes s'y écrivent, comme les
  pKa et les E° ; la base est celle du paquet Python.

### Reste à faire

1. Le moteur pgf n'a pas les deux familles : leur courbe n'est pas une expression
   unique (pM se déduit de pL et du volume ; l'ensemble des solides change).
2. Le niveau 2, le pH non tamponné : étapes 1 et 2 faites en Python
   (complexométrie, précipitation, hydroxydes) ; reste l'étape 3, Lua et LaTeX.

## Niveau 2 : étude (18 octobre 2026)

### Pourquoi le faire

Le tampon n'est pas un détail du calcul : c'est ce qui rend le titrage possible.
Sans lui, l'EDTA libère ses protons en se fixant, et le pH s'effondre. Le
prototype de l'étude (ci-dessous), sur 20 mL de Ca²⁺ à 0,01 mol/L titrés par
Na₂H₂Y à 0,01 mol/L, sans tampon :

| V (mL) | 0,01 | 5 | 10 | 19,9 | 20 | 20,1 | 30 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| pH | 5,00 | 3,52 | 3,51 | 3,53 | 3,53 | 3,53 | 3,56 |
| pCa | 2,00 | 2,11 | 2,20 | 2,34 | 2,34 | 2,34 | 2,45 |

Plus de saut à l'équivalence : à pH 3,5, la constante conditionnelle du
complexe s'est effondrée. Le niveau 2 montre ce que le niveau 1 suppose, et la
même chose pour l'ammoniac qui titre Cu²⁺ en faisant monter le pH, ou pour la
soude qui précipite Mg(OH)₂. Il complète le projet sans le dénaturer, à une
condition : garder des équations exactes, résolues à la précision de la machine,
sans réaction prépondérante.

### Ce qui existe ailleurs

Les logiciels de spéciation — PHREEQC (USGS), Visual MINTEQ, ChemEQL, Reaktoro
— résolvent n'importe quel mélange, par la méthode du tableau de Morel et un
Newton à plusieurs inconnues. Aucun ne trace une courbe de titrage avec ses
points d'inflexion à la précision de la machine, ni ne tourne dans LaTeX. Les
prendre comme dépendance couperait le paquet en deux : le moteur Lua ne pourrait
plus suivre le moteur Python. On reprend leur formulation, pas leur code, et
seulement dans le cadre restreint où elle se résout par dichotomies.

### La formulation : composantes et bilans

Les inconnues sont les concentrations libres des **composantes** — H⁺, le
ligand L sous sa forme la plus basique, l'ion commun P, les métaux M, les
partenaires X —, en logarithme : u = log c. Toute espèce s'en déduit,
`c_s = K_s · 10^(ν_s · u)` : HₖY (ν_H = k, ν_L = 1), MLⱼ, OH⁻ (ν_H = −1,
K = Ke), MOH⁺, MHY⁻… Chaque composante a son bilan, `Σ_s ν_sc c_s = T_c`. Le
bilan de H est l'électroneutralité écrite en protons : un acide fort compte dans
T_H, une base forte en négatif, Na₂H₂Y apporte 2 C à T_H et C à T_L. Comme
toujours, chaque total est linéaire dans les solutions,
`T(V) = (V₀ T_bécher + V T_burette)/(V₀ + V)`.

### L'unicité : la convexité

Les bilans sont le gradient de

```text
Φ(u) = Σ_s c_s(u)/ln 10 − Σ_c T_c u_c,
```

dont la hessienne, `ln 10 · Σ_s c_s ν_s ν_sᵀ`, est définie positive dès que les
composantes sont indépendantes : Φ est strictement convexe, et l'équilibre à
volume donné est **unique**. Un solide ajoute une contrainte linéaire,
`ν_p · u ≤ log K_s,p` ; le domaine reste convexe, la solution unique, et la
quantité de solide est le multiplicateur de la contrainte. C'est le résultat
classique qui fonde la méthode du tableau de Morel.

La conséquence pratique est celle qui compte ici. Minimiser une fonction convexe
sur une partie des variables laisse une fonction convexe des autres, dont la
dérivée est le bilan restant : **le bilan de H, évalué après avoir résolu les
autres bilans à h fixé, est croissant en log h.** D'où des dichotomies
imbriquées, chacune sur une fonction monotone, sans Newton ni point de départ.
La dichotomie actuelle de `TitrageParPrecipitation.p` en est déjà un cas
particulier, à un seul niveau.

### L'algorithme : dichotomies imbriquées, à volume fixé

Une composante est **éliminable** quand son bilan se résout explicitement une
fois les autres fixées :

- un métal, si l'on connaît l et h : ses fractions (complexes MLⱼ, mais aussi
  hydroxo-complexes et complexes protonés, s'il y en a) sont celles du niveau 1 ;
- un partenaire de précipitation, si l'on connaît a et h :
  `x = min(T_X, x_sat(a))` toutes formes protonées comprises ;
- un métal qui précipite en hydroxyde, si l'on connaît h :
  `[M] = min(T_M, K_s/[OH⁻]ⁿ)`.

Restent au plus deux composantes **principales**, d'où deux niveaux :

| Cas | Dichotomie extérieure | Intérieure | Éliminées |
| --- | --- | --- | --- |
| complexométrie sans tampon (EDTA, NH₃) | log h | log l | métaux |
| précipitation, partenaire basique sans tampon (Mohr) | log h | log a | partenaires |
| hydroxydes (Mg(OH)₂ par la soude) | log h | — | métal |

Le prototype a vérifié la monotonie du bilan réduit en log h sur des grilles de
1 500 à 2 000 points, à six volumes pour Ca–EDTA et à cinq pour Cu–NH₃, de part
et d'autre des équivalences : elle tient partout, comme la convexité l'annonce.
Résultats de Cu²⁺ (CuSO₄ 0,01 mol/L, 20 mL) titré par NH₃ 0,1 mol/L sans tampon :
pH 9,33 → 10,73 et pCu 2,12 → 7,34 de 0,5 à 12 mL. Mg²⁺ 0,01 mol/L par la soude
0,1 mol/L (pKs 11,15) : la précipitation commence à pH 9,43, et le saut est
modeste, de 10,20 à 10,68 entre 3,9 et 4,1 mL.

**Le coût.** Deux dichotomies de 64 coupes font 4 096 évaluations par point,
0,3 à 0,5 ms la dichotomie intérieure dans le prototype, qui n'est pas
vectorisé. Vectorisée sur tous les points de la courbe, comme `racine_croissante`
aujourd'hui, la courbe se calcule en Python en une fraction de seconde. En Lua,
il faudra partir de l'intervalle du point précédent et le resserrer par une
méthode de Brent (dichotomie protégée), d'une dizaine d'itérations par niveau.

### Ce qu'on perd, ce qu'on garde

**On perd le volume explicite.** Le paramètre de la courbe devient V lui-même :
c'est la seule coordonnée dont on sache qu'elle est monotone. Le pH ne l'est pas
— dans l'exemple de l'EDTA sans tampon, il descend de 5,00 à 3,51 puis remonte —,
et rien ne garantit que pM le soit en général. **On perd aussi le moteur pgf**,
déjà exclu du niveau 1.

**On garde** les équations exactes, la dichotomie, l'unicité (démontrée, et non
constatée), et la précision de la machine sur les valeurs. Les dérivées restent
exactes : par le théorème des fonctions implicites, un système linéaire 2 × 2 par
point donne dh/dV et dl/dV, d'où d pM/dV et d pH/dV.

**Les points d'inflexion** se cherchent en V, comme minima de |dV/d pM|, par
`minimum_lisse` avec un pas proportionnel à la largeur du saut, estimée sur la
grille. Avec des dérivées exactes à 10⁻¹⁵ près, la position est juste à environ
10⁻¹² fois cette largeur. Il faudra vérifier que les sauts de complexométrie,
bien plus larges que le saut de Fe²⁺ par MnO₄⁻ (10⁻⁹ mL), ne posent pas le
problème qui a fait chercher les sauts rédox en E plutôt qu'en V.

**Le niveau 1 est le cas limite.** Avec h fixé, la dichotomie extérieure
disparaît et l'on retrouve le niveau 1. Les tests le vérifieront de deux façons :
exactement, en retirant H des composantes, et en limite, avec un tampon de plus
en plus concentré dans la solution.

### Ce qui est couvert, ce qui ne l'est pas

Couvert, avec deux niveaux de dichotomie au plus :

- l'EDTA sans tampon, introduit sous la forme qu'on pèse (Na₂H₂Y) ;
- un tampon décrit comme une solution (un couple acide-base du bécher), qui ne
  complexe pas le métal ;
- Cu²⁺, Ag⁺, Zn²⁺… titrés par l'ammoniac, avec la protonation du ligand ;
- le dosage de Liebig avec la protonation de HCN ;
- les hydroxo-complexes et complexes protonés d'un métal, s'ils sont tabulés ;
- les hydroxydes (Mg(OH)₂, Fe(OH)₃) précipités par une base ;
- Mohr sans tampon : le chromate prend des protons, le pH bouge.

Pas couvert :

- **deux ligands à la fois** — un tampon ammoniacal qui complexe le zinc qu'on
  titre par l'EDTA, ou un indicateur métallochrome (NET) en concurrence : trois
  composantes principales, donc trois niveaux de dichotomie, 64³ évaluations par
  point. Faisable, mais c'est là que le Newton des logiciels de spéciation
  devient le bon outil ;
- **complexation et précipitation ensemble** (AgCl qui se redissout dans
  l'ammoniac) : l'argent reste éliminable — une équation du second degré quand
  le solide est présent —, mais le cas est à étudier à part ;
- l'activité, la température : comme aujourd'hui.

### La validation

Comme pour le niveau 1 : une référence indépendante en mpmath, qui n'emprunte
rien au paquet — dilution à la main et Newton sur le système complet à 50
chiffres, à partir de la solution en double précision —, sur les cas couverts
ci-dessus. Les deux limites (h fixé, tampon concentré). L'accord des moteurs
Python et Lua par `test_lua.py`, celui des documents LaTeX par `test_latex.py`.

### L'interface

**Pas de pH, pas de tampon** (voir la décision à la fin). Aujourd'hui, un ligand ou
un partenaire qui a des pKa demande un pH, sans quoi c'est une erreur ; au niveau
2, l'absence de pH voudra dire que le pH se calcule, avec un avertissement. Les solutions gagnent ce qu'il
faut pour décrire l'acide-base — acide fort, base forte, couples (le tampon), et
la forme sous laquelle on introduit le ligand — ; la courbe gagne une colonne pH.
En LaTeX, de même : sans clé `pH`, le moteur Lua prend le niveau 2.

### La taille

De l'ordre de 500 à 700 lignes en Python, 400 à 500 en Lua, 400 de tests : un
peu plus que le niveau 1, très en dessous des 2 500 à 4 000 lignes estimées le
3 octobre pour le cas général, que ce découpage évite.

### L'ordre proposé

1. Python : la complexométrie sans tampon (H, L, métaux), avec sa référence
   mpmath et les deux limites.
2. Python : les hydroxydes et la précipitation sans tampon.
3. Lua et LaTeX, puis la documentation.

Le tampon se décrit des deux façons, dès l'étape 1 (décision du 18 octobre 2026) :
par son seul pH, idéal, comme au niveau 1, ou comme le contenu d'une solution
(NH₄Cl et NH₃ dans le bécher), dont le pH se calcule et peut bouger un peu. Le
second n'est pas un troisième calcul : c'est le niveau 2, avec un couple
acide-base de plus dans le bécher.

L'interface (décision du 18 octobre 2026) : **l'absence de pH veut dire que le pH
se calcule**, à partir du contenu des solutions. Pour qu'un pH oublié ne passe pas
inaperçu, ce choix s'annonce par un avertissement explicite chaque fois que le pH
compte, c'est-à-dire dès qu'un ligand ou un partenaire a des pKa, ou qu'une solution
contient de l'acide-base :

- en Python, un `UserWarning` à la construction du titrage — « pH non précisé :
  calculé à partir des solutions, sans tampon imposé ; préciser pH= pour un milieu
  tamponné » —, qu'on fait taire par le module `warnings` quand c'est voulu ;
- en LaTeX, un `\msg_warning:nnn` du paquet, au log et au terminal, avec le nom du
  titrage.

Sans acidité nulle part, le pH ne compte pas et rien ne s'affiche : c'est le
niveau 1 tel qu'il est aujourd'hui.
