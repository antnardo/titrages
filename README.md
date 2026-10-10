# titrages

[![CI](https://github.com/antnardo/titrages/actions/workflows/ci.yml/badge.svg)](https://github.com/antnardo/titrages/actions/workflows/ci.yml)
[![Licence MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENSE)

Le pH, la conductivité ou le potentiel d'une solution quelconque, et les courbes de
titrage pH-métrique, conductimétrique et potentiométrique d'un bécher quelconque par
une burette quelconque — acides et bases forts ou faibles, polyacides, ampholytes,
mélanges, titrant faible, titrage en retour, couples rédox à plusieurs étapes ou
dimériques — avec la dilution, **sans aucune approximation**. Un paquet Python, et
un paquet LaTeX qui calcule la même chose dans `pgfplots`.

```pycon
>>> from titrages import Titrage, solution
>>> t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
>>> t.pH([0, 5, 10]).round(2).tolist()
[2.88, 4.76, 8.73]
>>> t.volumes_equivalents()
(10.0,)
>>> [(round(s.volume, 2), round(s.pH, 2)) for s in t.sauts(20)]
[(10.0, 8.73)]

```

Le paquet s'installe depuis GitHub (ci-dessous). PyPI et CTAN : à venir.

## Sommaire

- [Ce qu'il fait](#ce-quil-fait)
- [Installation](#installation)
- [Prérequis](#prérequis)
- [Prise en main](#prise-en-main)
- [La documentation](#la-documentation)
- [Les modules](#les-modules)
- [En ligne de commande](#en-ligne-de-commande)
- [Le paquet LaTeX](#le-paquet-latex)
- [Avec pychimie](#avec-pychimie)
- [Tester](#tester)
- [Antériorité](#antériorité)
- [Pour plus tard : complexation et précipitation](#pour-plus-tard--complexation-et-précipitation)
- [Limites](#limites)
- [Sources](#sources)
- [Licence](#licence)

## Ce qu'il fait

- **Le sens direct.** À pH imposé, tout est explicite : les fractions de chaque forme,
  puis l'électroneutralité, qui donne le volume versé comme une fonction directe du
  pH, `V(pH) = -V0 R_bécher(pH) / R_burette(pH)`. Pas de polynôme à résoudre, pas de
  réaction prépondérante. C'est la méthode de Kalka (2021) et de Brière (2002),
  généralisée à un bécher et une burette quelconques.
- **Le sens inverse, quand il le faut.** Le pH d'une solution ou le pH à un volume
  donné se trouvent par dichotomie sur une fonction qu'on sait croissante : une
  racine unique, trouvée à la précision de la machine, sans réglage.
- **La dilution, exactement.** Le résidu du mélange, multiplié par `V0 + V`, se
  sépare en une part du bécher et une part de la burette, chacune calculée avec ses
  concentrations d'origine : la dilution ne demande aucun traitement particulier.
  Tracer `n0 = Y1 + δ/C` comme une courbe de titrage, ainsi que le fait Kalka (2021),
  dont la méthode est tirée, c'est la négliger — jusqu'à 0,48 unité de pH d'erreur.
  Les tests confrontent la formule à un calcul indépendant qui dilue à la main.
- **Tous les cas d'un seul tenant.** Une base faible est la forme basique d'un
  couple (l'ammoniac est la forme 1 de NH₄⁺/NH₃) ; un sel s'introduit par sa forme,
  ses ions spectateurs se déduisent. Mélanges, ampholytes, titrant faible, titrage
  en retour ne demandent rien de plus.
- **Les équivalences, des deux façons.** Les volumes stœchiométriques, et les points
  d'inflexion effectivement visibles : l'acide phosphorique a trois équivalences
  stœchiométriques, et deux sauts seulement.
- **La méthode des tangentes**, construite exactement sur la courbe, pH-métrique ou
  potentiométrique, pour mesurer son erreur plutôt que la supposer.
- **Une base de couples, d'indicateurs, de conductivités et de couples rédox**,
  cherchée par ce qu'on écrit : `"NH3"`, `"ammoniac"`, `"carbonate de sodium"`,
  `"CO3^2-"`.
- **Une précision sans débordement.** Les fractions se calculent en échelle
  logarithmique : l'EDTA à pH 14 ne déborde pas. Les points d'inflexion sont
  confrontés à un calcul mpmath à 60 chiffres (`tests/reference_mpmath.py`).
- **La conductimétrie** par la loi de Kohlrausch, ions spectateurs compris (le Na⁺
  de la soude, celui du carbonate de sodium), avec ou sans correction de dilution,
  pour 34 des 35 couples de la base.
- **La potentiométrie** : les électrons remplacent les protons, E/s le pH. Couples à
  plusieurs étapes (le vanadium), couples qui échangent des protons (à pH fixé),
  couples dimériques (I₂/I⁻, Cr₂O₇²⁻/Cr³⁺, S₄O₆²⁻/S₂O₃²⁻), électrode de référence.
- **Le même calcul en LaTeX**, par un moteur Lua sous LuaLaTeX, ou par pgfmath et
  l3fp sous pdfLaTeX. Les tests confrontent les trois implémentations.

## Installation

Quatre façons, selon le poste.

**Avec pip depuis GitHub**, quand le réseau et pip sont disponibles :

```bash
pip install git+https://github.com/antnardo/titrages
pip install "titrages[trace] @ git+https://github.com/antnardo/titrages"   # avec matplotlib
```

**Avec la roue d'une release**, pour un poste qui a pip mais pas git : chaque
[release](https://github.com/antnardo/titrages/releases) porte
`titrages-<version>-py3-none-any.whl` (et le sdist `.tar.gz`, avec `docs/`,
`examples/` et `tests/`), qui s'installe par `pip install titrages-0.1.0-py3-none-any.whl`.

**Copier le dossier**, sans pip ni réseau : le dossier `src/titrages/` se pose à
côté des scripts, et `from titrages import Titrage, solution` fonctionne, à condition
d'exécuter le script depuis son dossier. Le paquet ne demande que numpy.

**Pour développer**, un clone et [uv](https://docs.astral.sh/uv/) :

```bash
git clone https://github.com/antnardo/titrages
cd titrages && uv sync && uv run pytest
```

PyPI et CTAN : à venir.

## Prérequis

Python 3.11 ou plus récent, avec numpy (1.24 ou plus) ; matplotlib (3.7 ou plus)
seulement pour `titrages.trace`. Le paquet est testé sous Python 3.11, 3.12, 3.13
et 3.14. Le paquet LaTeX demande pgfplots, et LuaLaTeX pour tout ce qui n'est pas
une courbe ou un nombre (voir [Le paquet LaTeX](#le-paquet-latex)).

## Prise en main

Deux règles valent partout. **Un `float` pour un scalaire, un tableau numpy pour un
tableau** : `t.pH(5)` est un `float`, `t.pH([5, 10])` un tableau. **Une seule unité
pour la conductivité, le S/m** — celle que donnent λ° en mS·m²·mol⁻¹ et *c* en
mol/L ; un conductimètre affiche des mS/cm, dix fois moins, et `unite="mS/cm"` le
fait pour un tracé. Les volumes s'écrivent dans l'unité qu'on veut, la même partout ;
les concentrations en mol/L.

Une courbe pH-métrique, son saut, la méthode des tangentes et la zone de virage
d'un indicateur :

```python
import matplotlib.pyplot as plt
from titrages import Titrage, indicateur, solution
from titrages.trace import tracer_sauts, tracer_tangentes, tracer_titrage, tracer_zone_virage

t = Titrage(solution({"H3PO4": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
ax = tracer_titrage(t, 35, lw=2)
tracer_sauts(t, 35, ax)
tracer_tangentes(t, 7, 35, ax)      # depuis 7 mL, avant le premier saut
tracer_zone_virage(indicateur("phénolphtaléine"), ax)
plt.show()
```

La méthode des tangentes, `t.tangentes(volume_avant, volume_max)`, part d'un point
avant le saut et rend la construction complète (`Tangentes` : les deux points de
contact, la pente commune, le point trouvé) :

```pycon
>>> t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
>>> tg = t.tangentes(7, 20)
>>> round(tg.volume, 2), round(tg.pH, 2)
(10.0, 8.65)

```

Un mélange par un titrant faible, et un titrage en retour (la soude en excès a déjà
transformé l'acide éthanoïque en éthanoate, que l'acide chlorhydrique titre après
l'excès) :

```pycon
>>> t = Titrage(solution({"HCl": 0.05, "CH3COOH": 0.05}), solution({"NH3": 0.1}), volume_becher=20)
>>> [round(s.volume, 1) for s in t.sauts(30)]
[10.0, 20.0]
>>> Titrage(solution({"CH3COOH": 0.05, "NaOH": 0.1}), solution({"HCl": 0.1}), 10).volumes_equivalents()
(5.0, 10.0)

```

Le pH d'une solution, et ce qu'elle contient :

```pycon
>>> from titrages import solution
>>> s = solution({"Na2CO3": 0.1})
>>> round(s.pH(), 2)
11.65
>>> round(s.concentrations()["HCO3-"], 4)
0.0045

```

La conductimétrie, en S/m, corrigée ou non de la dilution, et son tracé par
`tracer_conductimetrie(t, volume_max, ax, corrigee=True, unite="mS/cm")` :

```pycon
>>> t = Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), volume_becher=100)
>>> round(t.conductivite(10), 4)
0.1149
>>> t.courbe_conductimetrique(20, 3, corrigee=True).conductivite.round(4).tolist()
[0.4261, 0.1264, 0.3751]

```

La potentiométrie, mesurée contre une électrode au calomel saturé. Le pH du milieu,
supposé tamponné, se donne aux **deux** solutions, bécher et burette, puisqu'il décale
les potentiels standard apparents ; `tracer_titrage_redox(r, volume_max, ax,
sauts=True)` trace la courbe, et `tracer_sauts`, `tracer_tangentes`, `tracer_derivee`
acceptent un titrage rédox comme un titrage acide-base :

```pycon
>>> from titrages import TitrageRedox, solution_redox
>>> r = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"KMnO4": 0.02}), 10, reference="ECS")
>>> r.volumes_equivalents(), round(r.E(10), 3)
((10.0,), 1.143)
>>> v = TitrageRedox(solution_redox({"V2+": 0.1}, pH=1), solution_redox({"KMnO4": 0.02}, pH=1), 10)
>>> v.volumes_equivalents()
(10.0, 20.0, 30.0)

```

## La documentation

- [docs/DOC.md](docs/DOC.md) : l'API complète, la méthode et ses formules.
- [examples/](examples/) : six scripts à lancer, décrits dans
  [examples/README.md](examples/README.md).
- Le paquet LaTeX : [latex/README.md](latex/README.md), et sa documentation
  `titrages-doc.pdf` (index, clés, messages), jointe à chaque
  [release](https://github.com/antnardo/titrages/releases) et à l'artefact de la
  CI ; pour la recompiler, voir [latex/README.md](latex/README.md).
- [CHANGELOG.md](CHANGELOG.md).

## Les modules

| Module | Rôle |
| --- | --- |
| `titrages.base` | la base de couples, d'indicateurs, de conductivités et de couples rédox ; `solution`, `solution_redox`, `espece`, `couple`, `indicateur`, `chercher` |
| `titrages.solution` | `Solution` : pH, concentrations, pouvoir tampon, conductivité d'une solution quelconque |
| `titrages.titrage` | `Titrage` : la courbe V(pH), le pH à volume donné, les sauts, les équivalences, les tangentes, la conductimétrie |
| `titrages.redox` | `SolutionRedox`, `TitrageRedox` : les mêmes, les électrons remplaçant les protons |
| `titrages.equilibre` | les fractions, les moments, le résidu de l'électroneutralité, les dichotomies |
| `titrages.modeles` | les objets du domaine : `Couple`, `Espece`, `Courbe`, `Saut`, `Tangentes`, leurs analogues rédox |
| `titrages.trace` | les tracés matplotlib |
| `titrages.export` | `vers_pgfplots` |
| `titrages.__main__` | la commande `titrages` |

Tout ce qui sert couramment s'importe de `titrages` directement.

## En ligne de commande

```bash
titrages pH CH3COOH=0.1 CH3COO-=0.05
titrages sigma HCl=0.01                                   # en S/m
titrages E Fe2+=0.1 Fe3+=0.01
titrages titrage H3PO4=0.1 --burette NaOH=0.1 --volume 10
titrages redox V2+=0.1 --burette KMnO4=0.02 --reference ECS
titrages couples phosph
titrages indicateurs
```

## Le paquet LaTeX

Le dossier [`latex/`](latex/) contient `titrages.sty` et `titrages.lua`. Ils
s'installent par copie des deux fichiers dans un dossier que TeX trouve, par exemple
`TEXMFHOME/tex/latex/titrages/`, ou à côté du document ; chaque release porte
aussi `titrages-ctan.zip`, l'archive prête pour CTAN, avec la documentation compilée.
Le paquet demande pgfplots, et `\pgfplotsset{compat=1.18}` dans le préambule :

```latex
\usepackage{titrages}
\pgfplotsset{compat=1.18}
\titragedefinir{ethanoique}{
  volume = 10,
  becher = { espece = { pKa = 4.76, C = 0.1 } },
  burette = { base forte = 0.1 },
}
\begin{tikzpicture}
\begin{axis}[xlabel={$V$ (mL)}, ylabel=pH]
  \addplottitrage[thick]{ethanoique}
  \addplottitragetangentes[red, dashed]{ethanoique}{7}
\end{axis}
\end{tikzpicture}
Le pH à l'équivalence vaut \num{\titragepH{ethanoique}{10}}.
```

Deux moteurs derrière la même interface. Sous LuaLaTeX, `titrages.lua` calcule tout
— courbes pH, σ et E, dérivée, sauts, équivalences, tangentes, nombres
développables —, en une seconde. Sous pdfLaTeX, le moteur pgf trace les courbes par
pgfmath (quatre à cinq chiffres, invisible au tracé) et donne par l3fp les nombres
qui ne demandent qu'une évaluation ou une dichotomie : `\titragevolume`,
`\titragepH`, `\titrageE`, `\titrageconductivite` ; pas d'équivalences, de sauts
ni de tangentes, ni de couples dimériques. Son coût : 250 lignes d'expl3 à maintenir
en parallèle du moteur Lua, et une compilation deux à trente fois plus lente (une
demi-seconde par nombre). `tests/test_latex.py` vérifie que les deux moteurs donnent
les nombres du paquet Python.

## Avec pychimie

`titrages` ne dépend que de numpy et ne connaît que des concentrations en mol/L.
Pour partir d'une masse pesée, le paquet voisin
[pychimie](https://github.com/antnardo/pychimie) donne la masse molaire ; le pont
tient en une ligne, dans le script de l'utilisateur, sans que l'un importe l'autre :

```python
from pychimie import masse_molaire      # un paquet à part, installé à part
from titrages import Titrage, solution

C = 0.384 / masse_molaire("NaHCO3") / 0.100      # 384 mg dans 100 mL, en mol/L
t = Titrage(solution({"NaHCO3": C}), solution({"HCl": 0.1}), volume_becher=50)
```

Les deux paquets ont une classe `Solution` de sens différents ; quand les deux
cohabitent, importer l'une sous un alias (`from titrages import Solution as
SolutionTitrage`). Le renommage est prévu du côté de pychimie.

## Tester

```bash
uv run pytest                    # tout, moteur Lua et paquet LaTeX compris quand TeX est là
uv run pytest -m "not latex"     # sans compiler le paquet LaTeX (une quarantaine de secondes)
uv run ruff check src tests && uv run ruff format --check src tests && uv run mypy
cd latex/tests && texlua test_titrages.lua
```

Les tests marqués `lua` confrontent le moteur Lua au moteur Python par `texlua` (ou
`lua`) ; ceux marqués `latex` compilent `latex/tests/titrages-test-nombres.tex` par
`lualatex` et par `pdflatex` et comparent leurs nombres au calcul Python ; les deux
groupes sont sautés quand TeX manque. Les points d'inflexion de référence se
régénèrent par `uv run python tests/reference_mpmath.py` (mpmath, dépendance de
développement). L'intégration continue (`.github/workflows/ci.yml`) rejoue tout
cela à chaque push, sous Python 3.11 à 3.14, et compile la documentation LaTeX ; à
chaque release publiée, `publish.yml` joint à la release la roue, le sdist, la
documentation PDF et l'archive CTAN.

## Antériorité

[pHcalc](https://github.com/rnelsonchem/pHcalc) et
[chempy](https://github.com/bjodah/chempy) calculent le pH d'une composition
donnée (chempy traite aussi complexation et précipitation) ; ni l'un ni l'autre
n'a d'objet titrage, de dilution exacte en V(pH), de sauts, de tangentes, de
conductimétrie ni de rédox. Côté LaTeX, [phtitration](https://ctan.org/pkg/phtitration)
(CTAN, 2026) trace huit configurations de titrage par un titrant fort, par
Henderson-Hasselbalch par morceaux pour les monoacides et par le calcul paramétrique
exact pour les polyacides, sans mélange, sans nombres, ni conductimétrie ni rédox ;
`titrages` part du même calcul paramétrique et le généralise à un bécher et une
burette quelconques. Sans équivalent ailleurs : le couple Python + LaTeX à même
calcul, la conductimétrie avec spectateurs, le rédox dimérique, la base en français.

## Pour plus tard : complexation et précipitation

Le paquet ne traite ni la complexation ni la précipitation : un titrage par le
chlorure de baryum ou par l'EDTA se calcule à la main. L'étude
[docs/PROJET-complexation-precipitation.md](docs/PROJET-complexation-precipitation.md)
dit ce que coûterait leur ajout, en deux niveaux. Le **niveau 1** garde la méthode
explicite en paramétrant la courbe par la grandeur mesurée, pM ou pAg : la
complexométrie en milieu tamponné (EDTA, constante conditionnelle), l'argentimétrie
(Mohr, Charpentier-Volhard) et les mélanges d'halogénures y entrent, avec des
courbes à points anguleux là où un solide apparaît ou disparaît ; 800 à 1 200 lignes,
Lua et LaTeX compris, le moteur pgf pouvant suivre. Le **niveau 2** — ligand basique
en milieu non tamponné, complexes successifs, hydroxydes — demande un solveur de
spéciation (Newton sur les logarithmes des concentrations libres, continuation le
long de la courbe, problème de complémentarité pour les solides) : 2 500 à 4 000
lignes, sans moteur pgf, et une précision qui n'est plus celle du calcul explicite.
Les deux demandent une base de constantes (β, Ks) sourcées. Recommandation de
l'étude : le niveau 1 d'abord, si un TP le demande ; sinon chempy ponctuellement.

## Limites

Solutions idéales (activité = concentration, conductivités molaires limites),
température par le seul `pKe` et la pente de Nernst. Les couples rédox qui échangent
des protons sont pris à pH fixé ; les métaux et les gaz ne sont pas des couples de
solution. Les pKa, E° et λ° de la base sont ceux des tables usuelles à 25 °C ; ils
diffèrent d'une table à l'autre, et la base le signale quand l'écart est notable. Les
λ° de dix-sept ions, non tabulés dans le CRC, sont approchés par un ion voisin,
nommé dans `src/titrages/donnees/conductivites.toml` ; ceux de l'EDTA manquent.

## Sources

- H. Kalka, *Polyprotic Acids and Beyond — An Algebraic Approach*, Chemistry 2021,
  3, 454-508, [doi:10.3390/chemistry3020034](https://doi.org/10.3390/chemistry3020034).
- T. Brière, *Une méthode simple de résolution des problèmes de pH et de titrages
  acido-basiques sans approximations*, BUP 841, 2002.
- Y. Cortial, *Courbe de dosage acidobasique et méthode des tangentes*, BUP 734, 1991.
- W. M. Haynes (dir.), *CRC Handbook of Chemistry and Physics*, 97e éd., CRC Press,
  2016 : constantes d'acidité, potentiels standard, conductivités molaires ioniques.

## Licence

MIT.
