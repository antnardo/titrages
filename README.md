# titrages

[![CI](https://github.com/antnardo/titrages/actions/workflows/ci.yml/badge.svg)](https://github.com/antnardo/titrages/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/titrages.svg)](https://pypi.org/project/titrages/)
[![Python](https://img.shields.io/pypi/pyversions/titrages.svg)](https://pypi.org/project/titrages/)
[![Licence MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENSE)

Le pH, la conductivité ou le potentiel d'une solution quelconque, et les courbes de
titrage pH-métrique, conductimétrique et potentiométrique d'un bécher quelconque par
une burette quelconque — acides et bases forts ou faibles, polyacides, ampholytes,
mélanges, titrant faible, couples rédox à plusieurs étapes ou dimériques — avec la
dilution, **sans aucune approximation**. Un paquet Python, et un paquet LaTeX qui
calcule la même chose dans `pgfplots`.

```pycon
>>> from titrages import Titrage, solution
>>> t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
>>> [round(float(p), 2) for p in t.pH([0, 5, 10])]
[2.88, 4.76, 8.73]
>>> t.volumes_equivalents()
(10.0,)
>>> [(round(s.volume, 2), round(s.pH, 2)) for s in t.sauts(20)]
[(10.0, 8.73)]

```

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
  Tracer `n0 = Y1 + δ/C` comme une courbe de titrage, ainsi que le fait
  Kalka (2021), dont la méthode est tirée, c'est la négliger — jusqu'à 0,48 unité de pH d'erreur. Les
  tests confrontent la formule à un calcul indépendant qui dilue à la main.
- **Tous les cas d'un seul tenant.** Une base faible est la forme basique d'un
  couple (l'ammoniac est la forme 1 de NH₄⁺/NH₃) ; un sel s'introduit par sa forme,
  ses ions spectateurs se déduisent. Mélanges, ampholytes, titrant faible, titrage
  en retour ne demandent rien de plus.
- **Les équivalences, des deux façons.** Les volumes stœchiométriques, et les points
  d'inflexion effectivement visibles : l'acide phosphorique a trois équivalences
  stœchiométriques, et deux sauts seulement.
- **La méthode des tangentes**, construite exactement sur la courbe, pour mesurer
  son erreur plutôt que la supposer.
- **Une base de couples, d'indicateurs, de conductivités et de couples rédox**, cherchée par ce qu'on écrit : `"NH3"`,
  `"ammoniac"`, `"carbonate de sodium"`, `"CO3^2-"`.
- **Une précision sans débordement.** Les fractions se calculent en échelle
  logarithmique : l'EDTA à pH 14 ne déborde pas.
- **La conductimétrie** par la loi de Kohlrausch, ions spectateurs compris (le Na⁺
  de la soude, celui du carbonate de sodium), avec ou sans correction de dilution.
- **La potentiométrie** : les électrons remplacent les protons, E/s le pH. Couples à
  plusieurs étapes (le vanadium), couples qui échangent des protons (à pH fixé),
  couples dimériques (I₂/I⁻, Cr₂O₇²⁻/Cr³⁺, S₄O₆²⁻/S₂O₃²⁻), électrode de référence.
- **Le même calcul en LaTeX**, par un moteur Lua sous LuaLaTeX, ou par pgfmath sous
  pdfLaTeX. Les tests confrontent les deux implémentations.

## Installer

```bash
pip install titrages              # le calcul, numpy seul
pip install "titrages[trace]"     # et les tracés matplotlib
```

Python 3.11 ou plus récent.

## Exemples

```python
import matplotlib.pyplot as plt
from titrages import Titrage, indicateur, solution
from titrages.trace import tracer_sauts, tracer_tangentes, tracer_titrage, tracer_zone_virage

t = Titrage(solution({"H3PO4": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
ax = tracer_titrage(t, 35, lw=2)
tracer_sauts(t, 35, ax)
tracer_zone_virage(indicateur("phénolphtaléine"), ax)
plt.show()
```

Un mélange, par un titrant faible :

```pycon
>>> t = Titrage(solution({"HCl": 0.05, "CH3COOH": 0.05}), solution({"NH3": 0.1}), volume_becher=20)
>>> [round(s.volume, 1) for s in t.sauts(30)]
[10.0, 20.0]

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

Une conductimétrie, et un titrage rédox mesuré contre une électrode au calomel :

```pycon
>>> from titrages import TitrageRedox, solution_redox
>>> t = Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), volume_becher=100)
>>> [round(float(x), 4) for x in t.courbe_conductimetrique(20, 3, corrigee=True).conductivite]
[0.4261, 0.1264, 0.3751]
>>> r = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"KMnO4": 0.02}), 10, reference="ECS")
>>> r.volumes_equivalents(), round(float(r.E(10)), 3)
((10.0,), 1.143)

```

En ligne de commande :

```bash
titrages pH CH3COOH=0.1 CH3COO-=0.05
titrages sigma HCl=0.01
titrages E Fe2+=0.1 Fe3+=0.01
titrages titrage H3PO4=0.1 --burette NaOH=0.1 --volume 10
titrages redox V2+=0.1 --burette KMnO4=0.02 --reference ECS
titrages couples phosph
titrages indicateurs
```

## Le paquet LaTeX

Le dossier [`latex/`](latex/) contient `titrages.sty` et `titrages.lua`, publiés sur
CTAN :

```latex
\usepackage{titrages}
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
Le pH à l'équivalence vaut \titragepH{ethanoique}{10}.
```

Voir [latex/titrages-doc.pdf](latex/titrages-doc.pdf).

## Documentation

- [docs/DOC.md](docs/DOC.md) : l'API complète, la méthode et ses formules.
- [examples/](examples/) : des scripts à lancer.
- [CHANGELOG.md](CHANGELOG.md).

## Limites

Solutions idéales (activité = concentration, conductivités molaires limites),
température par le seul `pKe` et la pente de Nernst. Les couples rédox qui échangent
des protons sont pris à pH fixé ; les métaux et les gaz ne sont pas des couples de
solution. Les pKa, E° et λ° de la base sont ceux des tables usuelles à 25 °C ; ils
diffèrent d'une table à l'autre, et la base le signale quand l'écart est notable.

## Sources

- H. Kalka, *Polyprotic Acids and Beyond — An Algebraic Approach*, Chemistry 2021,
  3, 454-508, [doi:10.3390/chemistry3020034](https://doi.org/10.3390/chemistry3020034).
- T. Brière, *Une méthode simple de résolution des problèmes de pH et de titrages
  acido-basiques sans approximations*, BUP 841, 2002.
- Y. Cortial, *Courbe de dosage acidobasique et méthode des tangentes*, BUP 734, 1991.

## Licence

MIT.
