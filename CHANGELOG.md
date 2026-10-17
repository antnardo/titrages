# Changelog

Le format suit [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) et le
[versionnage sémantique](https://semver.org/lang/fr/). Le paquet Python et le paquet
LaTeX partagent un numéro de version.

## [Non publié]

### Ajouté

- La complexométrie (`TitrageComplexometrique`, `SolutionComplexation`,
  `Complexe`, `Ligand`, `fractions_complexes`) : pM à pH tamponné, complexes
  successifs et métaux en compétition, paramétrés par pL, le volume restant
  explicite ; pente analytique, sauts, volumes équivalents.
- Les titrages par précipitation (`TitrageParPrecipitation`,
  `SolutionPrecipitation`, `Precipite`) : pP de l'ion commun, plusieurs solides,
  partenaire basique à pH tamponné ; apparitions et disparitions de solides
  (`evenements`), sauts dans les morceaux lisses, volumes équivalents (Mohr,
  Charpentier-Volhard, halogénures successifs).
- `volumes_equivalents(stoechiometrie)` en complexométrie : le complexe qui marque
  l'équivalence, quand ce n'est pas le dernier déclaré (Liebig, `Ag(CN)2-`).
- Une base de ligands (13), de complexes (92) et de précipités (41), tirée des tables
  de D. Harvey, *Analytical Chemistry 2.1* (Martell et Smith) : `ligand`, `ligands`,
  `complexe`, `complexes`, `precipite`, `precipites`. Champs `note` sur `Ligand`,
  `Complexe` et `Precipite`, `commun` sur `Precipite` ; l'ion commun d'un
  `TitrageParPrecipitation` se déduit des précipités.
- Le tracé : `tracer_titrage_complexometrique`, `tracer_titrage_precipitation` ;
  `tracer_derivee` et `tracer_sauts` acceptent les deux familles ; `vers_pgfplots`
  exporte pM et pP.
- Le paquet LaTeX : `\titragecomplexationdefinir`, `\titrageprecipitationdefinir`,
  `\addplottitrageevenements`, `\titragepM`, `\titrageevenement`, sous le moteur
  lua ; `\addplottitrage`, `\addplottitragederivee`, `\addplottitragesauts`,
  `\titragevolume`, `\titrageequivalence` et `\titragesaut` s'y appliquent.
- Sections « La complexométrie » et « La précipitation » de `docs/DOC.md` et de la
  documentation LaTeX.

## [0.1.4] — 2026-10-11

### Ajouté

- `abscisse=` sur toutes les fonctions de tracé qui portent un volume en abscisse
  (`tracer_sauts`, `tracer_tangentes`, `tracer_derivee`, `tracer_conductimetrie`,
  `tracer_titrage_redox`), comme `tracer_titrage` : sauts et tangentes se superposent
  à une courbe tracée en V/V_éq.
- Les exemples (`examples/`) sont sous licence CC BY-NC-SA 4.0 (`LICENSE-CC-BY-NC-SA-4.0`),
  le code reste sous MIT ; le README et les métadonnées du paquet le disent.

### Corrigé

- `09_titrage_acide_fort.py` dessinait ses tangentes en millilitres sur une courbe en
  V/V_éq.

## [0.1.3] — 2026-10-11

### Ajouté

- La galerie des exemples `examples/README.md`, la courbe de chaque exemple dans
  `examples/images/` ; les exemples qui tracent enregistrent leur figure à côté
  d'eux.

## [0.1.2] — 2026-10-11

### Modifié

- Publié sur PyPI : `pip install titrages`. Le workflow de release y envoie la roue et le sdist par « trusted publishing », sans jeton ; GitHub garde les artefacts joints à la release.
- Les liens du README sont absolus, pour la page PyPI.

## [0.1.1] — 2026-10-10

### Modifié

- La tête du paquet n'exporte plus que ce qui sert couramment (34 noms) :
  `LN10`, `PKE`, `delta`, `log_k`, `LAMBDA_H3O`, `LAMBDA_HO`, `NERNST`,
  `REFERENCES`, `Composant`, `ComposantRedox`, `Fort`, `charge`,
  `conductivite_molaire`, `electrons_cedes` et `normaliser` restent importables
  de `titrages.equilibre`, `titrages.solution`, `titrages.redox`,
  `titrages.modeles` et `titrages.base`.

### Supprimé

- `vers_csv` : `numpy.savetxt` sur les colonnes d'une courbe fait la même chose ;
  `vers_pgfplots` et la commande `titrages` restent.

## [0.1.0] — 2026-10-10

Première version publiée (GitHub `antnardo/titrages`).

### Ajouté

- `SolutionRedox.diluee` et `melangee`, `TitrageRedox.tangentes` et
  `TangentesRedox` : chaque méthode acide-base a sa symétrique rédox ;
  `tracer_sauts`, `tracer_tangentes` et `tracer_derivee` acceptent un titrage
  rédox ; en LaTeX, `\addplottitragetangentes` et `\titragetangentes` s'appliquent
  à un titrage rédox (moteur Lua, `tangentes_redox`).
- Sous le moteur pgf (pdfLaTeX), les nombres `\titragevolume`, `\titragepH`,
  `\titragepHinitial`, `\titrageE` et `\titrageconductivite`, par l3fp,
  développables ; `latex/tests/titrages-test-nombres.tex` et `tests/test_latex.py`
  confrontent les deux moteurs au paquet Python, courbes de pgfplots comprises.
- `tests/reference_mpmath.py`, qui régénère à 60 chiffres les points d'inflexion
  de référence sans rien emprunter au paquet (mpmath, dépendance de développement).
- Vingt-sept conductivités molaires de plus (CRC Handbook, 97e éd. ; dix-sept
  approchées par un ion voisin, marquées ≈) : la conductimétrie couvre 34 couples
  sur 35 au lieu de 14.
- L'intégration continue à chaque push, et la publication par release GitHub
  (roue, sdist, documentation PDF, archive CTAN joints à la release).

### Modifié

- **Règle de sortie** : toute grandeur rend un `float` pour un argument scalaire,
  un tableau numpy pour un tableau (`Titrage.pH(5)` rendait un `ndarray` de
  dimension zéro).
- **Une seule unité pour la conductivité, le S/m** : la commande `sigma`,
  `tracer_conductimetrie` et la clé LaTeX `unite` ne convertissent plus en mS/cm
  par défaut ; `unite="mS/cm"` reste disponible.
- `equilibre.changement_de_signe` remplace `titrage._changement_de_signe`.
- Le README dit vrai : installation depuis GitHub ou par copie ; PyPI et CTAN à
  venir ; `latex/titrages-doc.pdf` n'est plus suivi par git, la release et la CI le
  fournissent.

- La documentation du paquet LaTeX est réécrite en `l3doc` : une fiche de syntaxe
  par commande, toutes les clés avec leur défaut, les messages et leur remède,
  l'interface Lua, la comparaison des deux moteurs, un index.
- Les points d'inflexion (`sauts`) se calculent à la précision de la machine, en
  acide-base comme en rédox : minima de |dV/dpH| ou |dV/dE| sur la courbe directe,
  affinés par le zéro de la dérivée de leur logarithme (`equilibre.minimum_lisse`)
  et non plus par une section dorée ; confrontés à `mpmath` dans les tests. Les
  sauts rédox se cherchaient en V, où celui de Fe²⁺ par MnO₄⁻ tient en 10⁻⁹ mL
  (1,3716 V au lieu de 1,37288 V).
- La variance des formes s'écrit Σ_{j<k} x_j x_k (k−j)², et Y₁ − q Σ x_j (j−q) :
  plus de soustraction de nombres voisins dans un saut. `derivee_volume` dérive
  V₀R_b + VR_bu = 0, somme de termes positifs.
- `TitrageRedox.pente` est analytique, couples dimériques compris ; nouvelle
  méthode `TitrageRedox.derivee_volume(E)`.

### Corrigé

- `volumes_equivalents` perdait l'équivalence intermédiaire d'un titrage en
  retour, quand le réactif fort du bécher excède les sites faibles :
  `solution({"CH3COOH": 0.05, "NaOH": 0.1})` par HCl 0,1 rendait `(10,)` au lieu
  de `(5, 10)`. L'acide ou la base forts du bécher réagissent maintenant d'abord
  avec les espèces faibles, chaque proton enlevé devenant un site basique au même
  pKa, et réciproquement ; deux sites de même pKa ne font qu'une équivalence. Même
  correction dans le moteur Lua.
- `charge` lisait la glycine `H3N+CH2COO-` comme un anion : un signe entre deux
  atomes compte pour ±1, le zwitterion est neutre et l'ion glycinium vaut +1.
- Un saut sur une courbe décroissante (base titrée par un acide, titrant réducteur)
  se rabattait sur un point de la grille : l'ammoniac par l'acide chlorhydrique
  donnait pH 5,2734 au lieu de 5,27558.
- La fraction oxydée d'un couple dimérique pouvait dépasser 1 d'un ulp, ce qui
  rendait défini le potentiel d'un bécher de diiode pur.

## [0.1.0-dev] — 23 au 27 septembre 2026, non publiée

Première version.

### Ajouté

- `Solution` : pH, concentrations, pouvoir tampon d'une solution quelconque —
  espèces faibles introduites sous n'importe quelle forme, acide et base forts.
- `Titrage` : la courbe V(pH) explicite, avec ou sans dilution, pour un bécher et
  une burette quelconques ; pH à volume donné ; points d'inflexion ; équivalences
  stœchiométriques ; méthode des tangentes parallèles.
- Fractions en échelle logarithmique, sans débordement.
- Une base de 35 couples et 9 indicateurs, cherchée par nom, sel ou formule.
- `titrages.trace` (matplotlib), `vers_pgfplots`, `vers_csv`, et la commande
  `titrages`.
- La conductimétrie : `Solution.conductivite`, `Titrage.conductivite`,
  `Titrage.courbe_conductimetrique` (loi de Kohlrausch, ions spectateurs, correction
  de dilution), une base de conductivités molaires et les contre-ions des sels.
- L'oxydoréduction : `CoupleRedox`, `SolutionRedox`, `TitrageRedox` — potentiel,
  courbe E(V), sauts, équivalences, électrode de référence ; couples à plusieurs
  étapes, à pH fixé, dimériques ; `volume(E)` explicite (quotient, ou racine d'un
  trinôme pour un seul dimère) ; une base de 15 couples.
- Des tests de dilution qui confrontent le calcul direct à un calcul indépendant,
  dilution faite à la main et équilibre résolu par force brute.
- Le paquet LaTeX `titrages` : moteur Lua sous LuaLaTeX (courbes, dérivée, sauts,
  tangentes, fractions, conductivité, potentiel, nombres développables), moteur
  pgfmath ailleurs (courbes pH, σ et E des couples monomères, fractions).

### Corrigé, par rapport au module de 2024

- La pente avec dilution multipliait par `n₀` au lieu de `n`, et appelait
  `Y(h, 1)` au lieu de `Y(1, h)`.
