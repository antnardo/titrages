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
- `volumes_equivalents()` en complexométrie compte une équivalence par étape, comme
  pour un polyacide : les étapes des complexes qui prédominent tour à tour (enveloppe
  concave des log β), deux étapes de même force confondues.
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
- La complexométrie à pH calculé : sans `pH=`, et dès que le pH compte, le pH se
  calcule à partir des solutions, avec un `UserWarning`. `SolutionComplexation`
  gagne `protons_ligand` (Na₂H₂Y) et `acidobasique` (une `Solution` : acides et bases
  forts, tampon) ; `TitrageComplexometrique.pH_du_melange`, colonne `pH` de
  `CourbeComplexometrique`. Dichotomies imbriquées à volume fixé, unicité par
  convexité, dérivées exactes ; confronté à une référence mpmath qui écrit
  l'électroneutralité. Python seulement pour l'instant.
- La précipitation à pH calculé, hydroxydes compris : un précipité de partenaire
  `"OH-"` (Mg(OH)₂, Ag₂O…) ; `SolutionPrecipitation.acidobasique` ;
  `TitrageParPrecipitation.pH_du_melange`, colonne `pH` de `CourbePrecipitation` ;
  équivalences d'une base forte qui précipite un hydroxyde. Quinze hydroxydes et
  oxydes dans la base. Le moteur commun aux deux familles est `titrages.libre`.
- Le paquet LaTeX : complexométrie et précipitation à pH calculé, hydroxydes compris.
  Sans la clé `pH`, et dès que le pH compte, le pH se calcule, avec un avertissement
  du paquet ; clés de solution `protons ligand`, `acide fort`, `base forte`,
  `espece` ; partenaire `OH-` pour un hydroxyde ; `\titragepH` lit le pH de ces
  titrages. Moteur Lua confronté au moteur Python (`test_lua.py`, `test_latex.py`).
- `fractions_complexes` vit dans `titrages.equilibre`, avec `log_somme_complexes` et
  `moments_complexes`.
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

Première version publiée.

### Ajouté

- `Solution` : pH, concentrations, pouvoir tampon et conductivité d'une solution
  quelconque — espèces faibles introduites sous n'importe quelle forme, acide et base
  forts, ions spectateurs.
- `Titrage` : la courbe V(pH) explicite, avec ou sans dilution, pour un bécher et une
  burette quelconques — titrant faible, mélanges, titrage en retour — ; pH à volume
  donné ; points d'inflexion à la précision de la machine (minima de |dV/dpH|, affinés
  par le zéro de la dérivée de leur logarithme) ; équivalences stœchiométriques ;
  méthode des tangentes parallèles ; conductimétrie par la loi de Kohlrausch, avec ou
  sans correction de dilution.
- L'oxydoréduction : `CoupleRedox`, `SolutionRedox`, `TitrageRedox` — potentiel,
  courbe E(V), pente analytique, sauts, équivalences, tangentes, électrode de
  référence ; couples à plusieurs étapes, à pH fixé, dimériques.
- Des calculs qui ne perdent pas de chiffres : fractions en échelle logarithmique,
  variance des formes et pente écrites sans soustraction de nombres voisins.
- Une base de 35 couples acide-base, 15 couples rédox, 9 indicateurs et des
  conductivités molaires ioniques (CRC Handbook, 97e éd. ; dix-sept approchées par un
  ion voisin, marquées ≈), cherchée par nom, sel ou formule.
- Une règle de sortie : un `float` pour un argument scalaire, un tableau numpy pour un
  tableau ; une seule unité pour la conductivité, le S/m.
- `titrages.trace` (matplotlib), `vers_pgfplots`, `vers_csv`, et la commande `titrages`.
- Le paquet LaTeX `titrages` : moteur Lua sous LuaLaTeX (courbes, dérivée, sauts,
  équivalences, tangentes, fractions, conductivité, potentiel, nombres
  développables), moteur pgfmath et l3fp sous pdfLaTeX (courbes pH, σ et E des
  couples monomères, fractions, nombres explicites) ; documentation en `l3doc`.
- Des tests qui confrontent le calcul à des calculs indépendants — dilution faite à
  la main, points d'inflexion à 60 chiffres par mpmath — et les moteurs Lua et pgf au
  paquet Python ; l'intégration continue, et la publication par release GitHub (roue,
  sdist, documentation PDF, archive CTAN).
