# titrages — courbes de titrage pour pgfplots

Le paquet trace dans `pgfplots` la courbe de titrage d'un bécher quelconque par une
burette quelconque — pH, conductivité, potentiel, pM, pAg —, avec la dilution et sans
aucune approximation, et donne les nombres qu'on y lit : pH à un volume donné,
équivalences, points d'inflexion, méthode des tangentes.

## In English

This package plots, within `pgfplots`, exact titration curves — pH, conductivity,
electrode potential, pM and pAg against the added volume — for any beaker titrated by
any burette: strong or weak acids and bases, polyprotic acids, ampholytes, mixtures,
weak titrants, multi-step and dimeric redox couples, successive complexes with
competing metals, and precipitation with several solids (Mohr, Charpentier–Volhard),
in a buffered medium or with the pH computed from the solutions. Dilution is taken
into account, and no predominant-reaction approximation is made. Under LuaLaTeX, an
engine written in Lua computes in double precision and also returns expandable
numbers: pH, potential, pM or conductivity at a given volume, equivalence volumes,
inflection points, the tangent method. Under pdfLaTeX, a fallback engine plots
acid–base, conductimetric and single-step redox curves. The documentation is in
French, as are the package's commands, keys and messages. The same algorithm is
available in Python, as the `titrages` package on PyPI.

License: MIT. Bugs and suggestions: <https://github.com/antnardo/titrages/issues>.

## Le but

Tracer une courbe de titrage dans un document demande d'ordinaire de la calculer
ailleurs, puis de recopier une liste de points ; ou de l'approcher par morceaux, avec
des formules qui ne valent ni pour un polyacide, ni pour un mélange, ni pour un
titrant faible. Ici, on décrit le bécher et la burette, et le paquet trace la courbe
exacte et rend les nombres en commandes développables :

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
  \addplottitragesauts[black, mark=*]{ethanoique}
\end{axis}
\end{tikzpicture}
pH initial \num{\titragepHinitial{ethanoique}}, demi-équivalence
\num{\titragepH{ethanoique}{5}}, saut à \qty{\titragesaut{ethanoique}{1}}{mL}.
```

ce qui donne, avec `\sisetup{locale=FR}` :

```text
pH initial 2,88, demi-équivalence 4,76, saut à 10,00 mL.
```

La sortie est une courbe `\addplot` dans l'axe du lecteur, avec ses options
`pgfplots`, et des nombres que `\num`, `\qty` ou `\pgfmathsetmacro` mettent en forme ;
le style reste au document.

## Les moyens

1. **À pH imposé, tout est explicite.** Le volume versé est une fonction directe du pH,
   `V(pH) = -V0 R_bécher(pH) / R_burette(pH)`, où R est le résidu de
   l'électroneutralité de chaque solution : la courbe est une courbe paramétrée, sans
   équation à résoudre.
2. **Un moteur Lua, sous LuaLaTeX.** `titrages.lua` calcule en double précision, en
   une seconde pour une douzaine de courbes : courbes pH, σ, E, pM et pP, dérivée,
   points d'inflexion à la précision de la machine, équivalences, tangentes,
   apparitions de solides, et le pH quand aucun tampon ne l'impose.
3. **Un moteur de secours, sous pdfLaTeX.** Le paquet construit l'expression de V(pH)
   et la confie à `pgfplots` (quatre à cinq chiffres, invisible au tracé) ; l3fp en
   tire le pH, le potentiel, la conductivité ou le volume à un point donné.
4. **Le même calcul qu'en Python.** Le paquet Python `titrages`, dans le même dépôt,
   implémente le même algorithme ; les tests comparent les deux moteurs LaTeX au
   calcul Python.

## Installation

Le paquet tient en deux fichiers qui vont ensemble, `titrages.sty` et `titrages.lua`.
Ils se posent dans un même dossier du chemin de TeX, par exemple
`TEXMFHOME/tex/latex/titrages/`, ou à côté du document. On les prend dans ce dossier
du dépôt, ou dans `titrages-ctan.zip`, l'archive jointe à chaque
[release](https://github.com/antnardo/titrages/releases), qui contient aussi la
documentation compilée. Il demande `pgfplots` (1.18 conseillé) ; LuaLaTeX pour tout
ce que le moteur de secours ne fait pas (ci-dessous). Le paquet n'est pas encore sur
CTAN.

## Les familles de titrages

| Commande | Titrage | Moteur |
| --- | --- | --- |
| `\titragedefinir` | pH-métrique et conductimétrique | lua, pgf |
| `\titrageredoxdefinir` | potentiométrique ; couples dimériques sous lua seulement | lua, pgf |
| `\titragecomplexationdefinir` | complexométrique, pM = f(V) | lua |
| `\titrageprecipitationdefinir` | par précipitation, pP = f(V), points anguleux | lua |

Les deux dernières se prennent à pH imposé (clé `pH`), ou à pH calculé à partir des
solutions, sans la clé `pH` : le paquet le dit alors par un avertissement.

```latex
\titragecomplexationdefinir{calcium}{
  volume = 20, pH = 10, ligand pKa = {0, 1.5, 2, 2.66, 6.16, 10.24},
  complexe = {metal = Ca2+, log beta = 10.7},
  becher = {metal = {Ca2+, 0.01}}, burette = {ligand = 0.01},
}
\titrageprecipitationdefinir{mohr}{
  volume = 20, pH = 7.5,
  precipite = {partenaire = Cl-, pKs = 9.75},
  precipite = {partenaire = CrO4^2-, pKs = 11.95, nu = 2, pKa = 6.5},
  becher = {partenaire = {Cl-, 0.05}, partenaire = {CrO4^2-, 0.005}},
  burette = {commun = 0.05}, volume max = 30,
}
pCa à l'équivalence \num{\titragepM{calcium}{20}} ; le chromate d'argent apparaît à
\qty{\titrageevenement{mohr}{2}}{mL}.
```

```text
pCa à l'équivalence 6,28 ; le chromate d'argent apparaît à 20,01 mL.
```

## Limites

- Le moteur de secours (pdfLaTeX) ne donne ni équivalences, ni points d'inflexion, ni
  tangentes, ni couples dimériques, ni complexométrie, ni précipitation ; ses nombres
  coûtent une demi-seconde chacun, contre quelques millisecondes en Lua.
- Les solutions sont idéales. La température n'intervient pas, mais on peut changer
  `pKe` et la pente de Nernst (25 °C par défaut) ; les constantes sont celles qu'on
  écrit, le paquet LaTeX n'a pas de base de données.
- Un seul ligand par titrage complexométrique ; complexation et précipitation ne se
  combinent pas dans un même titrage.

## Développement

`titrages-doc.pdf` (syntaxe de chaque commande, clés et défauts, messages, interface
Lua, comparaison des moteurs, index) se recompile par :

```sh
lualatex titrages-doc
makeindex -s gind.ist -o titrages-doc.ind titrages-doc.idx
lualatex titrages-doc
lualatex titrages-doc
```

ou par `uv run python scripts/ctan.py` depuis la racine du dépôt, qui construit aussi
l'archive CTAN.

`tests/test_titrages.lua` éprouve le moteur Lua hors de TeX
(`texlua test_titrages.lua`). `tests/titrages-test-lua.tex` et
`tests/titrages-test-pgf.tex` compilent les tracés des deux moteurs ;
`tests/titrages-test-nombres.tex`, compilé par lualatex et par pdflatex, et
`tests/titrages-test-cp.tex`, sous lualatex, écrivent des nombres que
`tests/test_latex.py`, à la racine du dépôt, compare au calcul Python.

## Licence

MIT. Copyright (c) 2026 antnardo. Dépôt : <https://github.com/antnardo/titrages>.
