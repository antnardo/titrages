# Changelog

Le format suit [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) et le
[versionnage sémantique](https://semver.org/lang/fr/). Le paquet Python et le paquet
LaTeX partagent un numéro de version.

## [Unreleased]

### Modifié

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

- Un saut sur une courbe décroissante (base titrée par un acide, titrant réducteur)
  se rabattait sur un point de la grille : l'ammoniac par l'acide chlorhydrique
  donnait pH 5,2734 au lieu de 5,27558.
- La fraction oxydée d'un couple dimérique pouvait dépasser 1 d'un ulp, ce qui
  rendait défini le potentiel d'un bécher de diiode pur.

## [0.1.0] — non publiée

Première version. Le calcul reprend un module privé de 2024
(`pyphysique/chimie/titrage.py`), réécrit.

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
