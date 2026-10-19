# titrages

[![PyPI](https://img.shields.io/pypi/v/titrages.svg)](https://pypi.org/project/titrages/)
[![CI](https://github.com/antnardo/titrages/actions/workflows/ci.yml/badge.svg)](https://github.com/antnardo/titrages/actions/workflows/ci.yml)
[![Licence MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](https://github.com/antnardo/titrages/blob/main/LICENSE)

La courbe de titrage d'un bécher quelconque par une burette quelconque — pH-métrie,
conductimétrie, potentiométrie, complexométrie, précipitation —, avec la dilution et
**sans aucune approximation**, et les nombres qu'on y lit : pH à un volume donné,
équivalences, points d'inflexion, méthode des tangentes. Un paquet Python, et un
paquet LaTeX qui calcule la même chose dans `pgfplots`.

## Sommaire

- [Le but](#le-but)
- [Les moyens](#les-moyens)
- [Installation](#installation)
- [Une solution, un titrage](#une-solution-un-titrage)
- [Tracer](#tracer)
- [Conductimétrie et potentiométrie](#conductimétrie-et-potentiométrie)
- [Complexométrie et précipitation](#complexométrie-et-précipitation)
- [Sans tampon : le pH calculé](#sans-tampon--le-ph-calculé)
- [La base de constantes](#la-base-de-constantes)
- [En ligne de commande](#en-ligne-de-commande)
- [Le paquet LaTeX](#le-paquet-latex)
- [La documentation](#la-documentation)
- [Antériorité](#antériorité)
- [Limites](#limites)
- [Développement](#développement)
- [Sources](#sources)
- [Licence](#licence)

## Le but

Une courbe de titrage se calcule d'ordinaire par morceaux — avant, à, après
l'équivalence —, chaque morceau par sa réaction prépondérante et sa formule
approchée. Cela tient pour un acide faible par la soude ; cela ne tient plus pour un
polyacide aux pKa voisins, un mélange, un titrant faible, un titrage en retour, ni
pour l'EDTA sans tampon, dont le pH s'effondre en cours de route. Résoudre
numériquement le polynôme en h ne vaut pas mieux : ses coefficients couvrent des
dizaines d'ordres de grandeur, et le saut, là où tout se joue, est son point le plus
mal conditionné.

`titrages` donne la courbe exacte, et les nombres qu'on y mesure, pour toute
composition :

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

Le paquet rend des nombres et des tableaux numpy : une courbe (`volume`, `pH`,
`pente`), des sauts, des équivalences, une construction des tangentes. Le tracé est
au choix du lecteur : matplotlib par `titrages.trace`, `pgfplots` par la liste de
coordonnées de `vers_pgfplots`, ou directement par le paquet LaTeX.

## Les moyens

1. **À pH imposé, tout est explicite.** Les fractions de chaque forme se calculent,
   puis l'électroneutralité donne le volume versé comme une fonction directe du pH,
   `V(pH) = -V0 R_bécher(pH) / R_burette(pH)` : pas de polynôme à résoudre. C'est la
   méthode de Kalka (2021) et de Brière (2002), généralisée à un bécher et une
   burette quelconques ; la dilution y est exacte, chaque résidu gardant les
   concentrations de sa solution.
2. **Le sens inverse par dichotomie.** Le pH à un volume donné annule une fonction
   croissante : une racine unique, trouvée à la précision de la machine, sans
   réglage ni point de départ.
3. **Des calculs qui ne perdent pas de chiffres.** Les fractions se calculent en
   échelle logarithmique — l'EDTA à pH 14 ne déborde pas —, et la pente sans
   soustraire deux nombres voisins. Un point d'inflexion est un minimum de
   `|dV/dpH|`, affiné jusqu'à la précision de la machine : 10⁻¹⁰ unité de pH.
4. **La même idée partout.** Le potentiel remplace le pH en rédox ; le ligand libre
   en complexométrie ; l'ion commun en précipitation, où l'on essaie chaque ensemble
   de solides et garde le seul cohérent. Quand le pH n'est pas tamponné, deux
   dichotomies imbriquées trouvent l'équilibre, dont l'unicité vient de la convexité
   des bilans.
5. **Chaque nombre est vérifié.** Les tests confrontent le paquet à des calculs
   indépendants en précision arbitraire (mpmath), qui diluent à la main, et les trois
   implémentations — Python, Lua, pgfmath — entre elles.

## Installation

```bash
pip install titrages
pip install "titrages[trace]"   # avec matplotlib, pour titrages.trace
```

Python 3.11 ou plus récent, et numpy 1.24 ou plus ; matplotlib 3.7 ou plus seulement
pour les tracés. Le paquet est testé sous Python 3.11 à 3.14.

La version de développement s'installe depuis GitHub
(`pip install git+https://github.com/antnardo/titrages`). Chaque
[release](https://github.com/antnardo/titrages/releases) porte aussi la roue
`titrages-<version>-py3-none-any.whl`, pour un poste sans git, et le sdist, avec
`docs/`, `examples/` et `tests/`. Sans pip ni réseau, le dossier `src/titrages/` posé
à côté d'un script suffit : `from titrages import Titrage, solution` fonctionne
depuis le dossier du script.

Le paquet LaTeX demande une distribution TeX avec `pgfplots` ; LuaLaTeX pour tout ce
qui n'est pas une courbe acide-base, rédox ou un nombre simple (voir
[Le paquet LaTeX](#le-paquet-latex)).

## Une solution, un titrage

Une solution se décrit par ce qu'on y dissout, en mol/L ; les ions spectateurs se
déduisent du sel :

```pycon
>>> s = solution({"Na2CO3": 0.1})
>>> round(s.pH(), 2)
11.65
>>> round(s.concentrations()["HCO3-"], 4)
0.0045

```

Un mélange par un titrant faible, et un titrage en retour — la soude en excès a déjà
transformé l'acide éthanoïque en éthanoate, que l'acide chlorhydrique titre après
l'excès :

```pycon
>>> t = Titrage(solution({"HCl": 0.05, "CH3COOH": 0.05}), solution({"NH3": 0.1}), volume_becher=20)
>>> [round(s.volume, 1) for s in t.sauts(30)]
[10.0, 20.0]
>>> Titrage(solution({"CH3COOH": 0.05, "NaOH": 0.1}), solution({"HCl": 0.1}), 10).volumes_equivalents()
(5.0, 10.0)

```

Les équivalences stœchiométriques et les points d'inflexion sont deux choses : l'acide
phosphorique a trois équivalences et deux sauts. La méthode des tangentes part d'un
point avant le saut et rend la construction complète :

```pycon
>>> t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
>>> tg = t.tangentes(7, 20)
>>> round(tg.volume, 2), round(tg.pH, 2)
(10.0, 8.65)

```

Deux règles valent partout : un `float` pour un scalaire, un tableau numpy pour un
tableau (`t.pH(5)`, `t.pH([5, 10])`) ; les volumes dans l'unité qu'on veut, la même
partout.

## Tracer

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

Le même tracé pour l'acide éthanoïque, par
[`examples/01_acide_faible.py`](https://github.com/antnardo/titrages/blob/main/examples/01_acide_faible.py) :

![Courbe pH-métrique de l'acide éthanoïque par la soude : saut à 10 mL, méthode des tangentes et zone de virage de la phénolphtaléine](https://raw.githubusercontent.com/antnardo/titrages/main/examples/images/01_acide_faible.png)

Chaque fonction de tracé dessine sur l'axe qu'on lui donne et le rend : elle se
compose avec le reste d'un script. `vers_pgfplots(t.courbe(20))` rend la même courbe
en coordonnées `pgfplots`.

## Conductimétrie et potentiométrie

La conductimétrie suit la loi de Kohlrausch, ions spectateurs compris (le Na⁺ de la
soude), en S/m — l'unité que donnent λ° en mS·m²·mol⁻¹ et c en mol/L ; un
conductimètre affiche des mS/cm, dix fois moins, et `unite="mS/cm"` le fait pour un
tracé :

```pycon
>>> t = Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), volume_becher=100)
>>> round(t.conductivite(10), 4)
0.1149
>>> t.courbe_conductimetrique(20, 3, corrigee=True).conductivite.round(4).tolist()
[0.4261, 0.1264, 0.3751]

```

En potentiométrie, les électrons remplacent les protons. Les couples à plusieurs
étapes, ceux qui échangent des protons (pris au pH du milieu, donné aux deux
solutions) et les couples dimériques (I₂/I⁻, Cr₂O₇²⁻/Cr³⁺) s'y traitent ; le
potentiel se lit contre une électrode de référence :

```pycon
>>> from titrages import TitrageRedox, solution_redox
>>> r = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"KMnO4": 0.02}), 10, reference="ECS")
>>> r.volumes_equivalents(), round(r.E(10), 3)
((10.0,), 1.143)
>>> v = TitrageRedox(solution_redox({"V2+": 0.1}, pH=1), solution_redox({"KMnO4": 0.02}, pH=1), 10)
>>> v.volumes_equivalents()
(10.0, 20.0, 30.0)

```

## Complexométrie et précipitation

Un métal titré par un ligand, à pH tamponné : la courbe pM = f(V), complexes
successifs et métaux en compétition compris ; les équivalences se comptent comme
celles d'un polyacide, une par étape de complexation.

```pycon
>>> from titrages import SolutionComplexation, TitrageComplexometrique, complexe, ligand
>>> t = TitrageComplexometrique(
...     SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01), 20,
...     [complexe("Ca2+", "EDTA")], ligand("EDTA"), pH=10,
... )
>>> round(t.pM(20), 2), [round(s.volume, 3) for s in t.sauts(40)], t.volumes_equivalents()
(6.28, [20.0], (20.0,))

```

Un titrage par précipitation suit l'ion commun, pAg = f(V). La courbe a des points
anguleux là où un solide apparaît ou disparaît — des *événements*, distincts des
sauts : dans la méthode de Mohr, le chromate d'argent apparaît juste après
l'équivalence du chlorure.

```pycon
>>> from titrages import SolutionPrecipitation, TitrageParPrecipitation, precipite
>>> m = TitrageParPrecipitation(
...     SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
...     SolutionPrecipitation(commun=0.05), 20,
...     [precipite("AgCl"), precipite("Ag2CrO4")], pH=7.5,
... )
>>> [round(v, 9) for v in m.volumes_equivalents()], [(e.partenaire, round(e.volume, 3)) for e in m.evenements(30)]
([20.0, 24.0], [('Cl-', 0.0), ('CrO4^2-', 20.012)])

```

## Sans tampon : le pH calculé

Sans `pH=`, et dès que le pH compte — un ligand ou un partenaire qui a des acidités,
un hydroxyde, une solution qui apporte de l'acide-base —, le pH se calcule à partir
du contenu des solutions, et un `UserWarning` le rappelle : un pH oublié change la
chimie. L'EDTA, pesé sous la forme Na₂H₂Y, libère ses protons en se fixant ; sans
tampon, le pH tombe et le saut disparaît :

```pycon
>>> import warnings
>>> with warnings.catch_warnings():
...     warnings.simplefilter("ignore")
...     s = TitrageComplexometrique(
...         SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01, protons_ligand=2), 20,
...         [complexe("Ca2+", "EDTA")], ligand("EDTA"),
...     )
>>> round(float(s.pH_du_melange(20)), 2), s.sauts(40)
(3.54, ())

```

Un tampon se décrit alors comme une solution (un couple acide-base dans le bécher et
la burette), et un hydroxyde comme un précipité de partenaire `"OH-"`, que précipite
une base forte.

## La base de constantes

Les espèces se cherchent par ce qu'on écrit — `"NH3"`, `"ammoniac"`,
`"carbonate de sodium"`, `"CO3^2-"` ; casse, accents et écriture des charges ne
comptent pas. La base porte 35 couples acide-base, 15 couples rédox, 9 indicateurs,
les conductivités molaires ioniques, 13 ligands, 92 complexes et 56 précipités, dont
15 hydroxydes et oxydes. Les complexes et les précipités viennent des tables de
D. Harvey, d'après Martell et Smith ; la source et les conditions de chaque valeur
sont dans les fichiers de données et dans le champ `note`.

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

Le dossier [`latex/`](https://github.com/antnardo/titrages/tree/main/latex/) contient
`titrages.sty` et `titrages.lua`, à poser ensemble dans un dossier que TeX trouve
(par exemple `TEXMFHOME/tex/latex/titrages/`) ou à côté du document ; chaque release
porte aussi `titrages-ctan.zip`, avec la documentation compilée.

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

Sous LuaLaTeX, `titrages.lua` calcule tout : courbes pH, σ, E, pM et pP, dérivée,
sauts, équivalences, tangentes, événements, pH calculé, et des nombres développables
(`\titragepH`, `\titragepM`, `\titragesaut`…). Sous pdfLaTeX, un moteur de secours
trace les courbes acide-base et rédox par pgfmath (quatre à cinq chiffres, invisible
au tracé) et donne par l3fp le pH, le potentiel, la conductivité ou le volume à un
point donné, plus lentement ; ni équivalences, ni sauts, ni tangentes, ni couples
dimériques, ni complexométrie, ni précipitation. Le détail est dans
[latex/README.md](https://github.com/antnardo/titrages/blob/main/latex/README.md) et
dans la documentation `titrages-doc.pdf`.

## La documentation

- [docs/DOC.md](https://github.com/antnardo/titrages/blob/main/docs/DOC.md) : l'API
  complète, la méthode et ses formules, les modules.
- [examples/](https://github.com/antnardo/titrages/tree/main/examples/) : onze
  scripts, décrits avec leurs courbes dans la galerie
  [`examples/README.md`](https://github.com/antnardo/titrages/blob/main/examples/README.md).
- `titrages-doc.pdf`, la documentation du paquet LaTeX (syntaxe, clés, messages,
  index), jointe à chaque release.
- [docs/technique/](https://github.com/antnardo/titrages/tree/main/docs/technique/) :
  le document technique, `titrages-technique.tex` — tous les calculs, les choix
  algorithmiques, les démonstrations d'unicité, et les écarts mesurés face aux calculs
  indépendants (mpmath), que `verifications.py` recalcule.
- [CHANGELOG.md](https://github.com/antnardo/titrages/blob/main/CHANGELOG.md).

Pour partir d'une masse pesée, le paquet [pychimie](https://github.com/antnardo/pychimie)
donne la masse molaire ; `titrages` n'en dépend pas, le pont tient dans le script :

```python
from pychimie import masse_molaire
from titrages import Titrage, solution

C = 0.384 / masse_molaire("NaHCO3") / 0.100      # 384 mg dans 100 mL, en mol/L
t = Titrage(solution({"NaHCO3": C}), solution({"HCl": 0.1}), volume_becher=50)
```

## Antériorité

[pHcalc](https://github.com/rnelsonchem/pHcalc) et
[chempy](https://github.com/bjodah/chempy) calculent l'équilibre d'une composition
donnée, chempy avec complexation et précipitation ; ni l'un ni l'autre n'a d'objet
titrage, de dilution exacte, de sauts, de tangentes, de conductimétrie ni de rédox.
Les logiciels de spéciation (PHREEQC, Visual MINTEQ) résolvent n'importe quel mélange,
sans courbe de titrage ni points d'inflexion. Côté LaTeX,
[phtitration](https://ctan.org/pkg/phtitration) trace huit configurations de titrage
par un titrant fort, sans mélange, sans nombres, ni conductimétrie ni rédox ;
`titrages` part du même calcul paramétrique exact et le généralise à un bécher et une
burette quelconques.

## Limites

- Les solutions sont idéales : l'activité est la concentration, les conductivités
  molaires sont les valeurs limites, qui surestiment σ dès quelques mmol/L. La
  température n'entre que par `pKe` et la pente de Nernst.
- Les couples rédox qui échangent des protons sont pris à pH fixé ; un métal ou un gaz
  (Ag, H₂) n'est pas un couple de solution.
- Un seul ligand par titrage : un tampon ammoniacal qui complexerait aussi le métal, ou
  un indicateur métallochrome, en serait un second. Complexation et précipitation ne
  se combinent pas dans un même titrage.
- Les constantes de la base sont celles des tables à 25 °C ; elles diffèrent d'une
  table à l'autre, et la base le signale quand l'écart est notable. Les λ° de
  dix-sept ions, non tabulés, sont approchés par un ion voisin, nommé ; ceux de l'EDTA
  manquent.

## Développement

```bash
git clone https://github.com/antnardo/titrages
cd titrages && uv sync
uv run pytest                    # tout, moteur Lua et paquet LaTeX compris quand TeX est là
uv run pytest -m "not latex"     # sans compiler le paquet LaTeX
uv run ruff check src tests && uv run ruff format --check src tests && uv run mypy
cd latex/tests && texlua test_titrages.lua
```

Les tests marqués `lua` comparent le moteur Lua au moteur Python par `texlua` ; ceux
marqués `latex` compilent les documents de `latex/tests/` par `lualatex` et
`pdflatex` et comparent leurs nombres au calcul Python. Les deux groupes sont sautés
quand TeX manque.

## Sources

- H. Kalka, *Polyprotic Acids and Beyond — An Algebraic Approach*, Chemistry 2021,
  3, 454-508, [doi:10.3390/chemistry3020034](https://doi.org/10.3390/chemistry3020034).
- T. Brière, *Une méthode simple de résolution des problèmes de pH et de titrages
  acido-basiques sans approximations*, BUP 841, 2002.
- Y. Cortial, *Courbe de dosage acidobasique et méthode des tangentes*, BUP 734, 1991.
- W. M. Haynes (dir.), *CRC Handbook of Chemistry and Physics*, 97e éd., CRC Press,
  2016 : constantes d'acidité, potentiels standard, conductivités molaires ioniques.
- D. Harvey, *Analytical Chemistry 2.1*, LibreTexts, § 16.10 à 16.12, d'après
  A. E. Martell et R. M. Smith, *Critical Stability Constants*, Plenum, 1976 :
  constantes de formation et produits de solubilité.

## Licence

Le code — le paquet Python (`src/`), le paquet LaTeX (`latex/`), les tests et les
scripts — est sous licence MIT (`LICENSE`). Les exemples (`examples/`), contenu
scientifique et pédagogique plutôt qu'algorithmes, sont sous licence CC BY-NC-SA 4.0
(`LICENSE-CC-BY-NC-SA-4.0`).
