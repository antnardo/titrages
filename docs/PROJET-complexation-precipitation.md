# Projet : complexation et précipitation

État au 18 octobre 2026 : le **niveau 1 est écrit côté Python** (voir
*Avancement*, à la fin). Ce document disait, le 3 octobre, ce que coûterait
l'ajout des équilibres de complexation et de précipitation, en deux niveaux, et
ce qu'on recommandait de faire d'abord.

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
- Le complexe qui marque l'équivalence peut se choisir (`stoechiometrie`) : dans le
  dosage de Liebig, c'est Ag(CN)₂⁻, pas le dernier complexe déclaré.

### Reste à faire

1. Le moteur pgf n'a pas les deux familles : leur courbe n'est pas une expression
   unique (pM se déduit de pL et du volume ; l'ensemble des solides change).
2. La base n'est pas accessible depuis LaTeX, comme celle des couples.
3. Le niveau 2 (pH non tamponné, complexation et précipitation couplées) reste hors
   du projet.
