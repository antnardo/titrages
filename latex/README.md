# titrages — courbes de titrage pour pgfplots

Le paquet trace les courbes pH = f(V), σ = f(V) et E = f(V) du titrage d'un bécher
quelconque par une burette quelconque — acides et bases forts ou faibles,
polyacides, mélanges, titrant faible, couples rédox à plusieurs étapes ou
dimériques — avec la dilution, sans aucune approximation.

La méthode : à pH imposé, tout est explicite, et le volume versé est une fonction
directe du pH, V(pH) = -V0 R_bécher(pH) / R_burette(pH), où R est le résidu de
l'électroneutralité de chaque solution.

- Sous LuaLaTeX, un moteur Lua (`titrages.lua`) calcule en double précision : la
  courbe, sa dérivée, les points d'inflexion, les équivalences, la méthode des
  tangentes, et le pH, la conductivité ou le potentiel à un volume donné, en
  nombres développables.
- Sous pdfLaTeX, un moteur pgf trace les courbes et les fractions par pgfmath, et
  donne par l3fp les nombres qui ne demandent qu'une évaluation ou une dichotomie
  (`\titragevolume`, `\titragepH`, `\titrageE`, `\titrageconductivite`), plus
  lentement ; ni équivalences, ni sauts, ni tangentes, ni couples dimériques.

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
\end{axis}
\end{tikzpicture}
```

## Installation

Le paquet tient en deux fichiers qui vont ensemble, `titrages.sty` et
`titrages.lua`. Ils se posent dans un même dossier du chemin de TeX, par exemple
`TEXMFHOME/tex/latex/titrages/`, ou à côté du document. On les prend dans ce
dossier du dépôt, ou dans `titrages-ctan.zip`, l'archive jointe à chaque
[release](https://github.com/antnardo/titrages/releases), qui contient aussi la
documentation compilée. Dépend de pgfplots (1.18 conseillé). CTAN : à venir.

## Documentation

`titrages-doc.pdf` (en français : syntaxe de chaque commande, clés et défauts,
messages, interface Lua, comparaison des moteurs, index) est jointe à chaque
release et à l'artefact de l'intégration continue. Pour la recompiler :

```sh
lualatex titrages-doc
makeindex -s gind.ist -o titrages-doc.ind titrages-doc.idx
lualatex titrages-doc
lualatex titrages-doc
```

ou `uv run python scripts/ctan.py` depuis la racine du dépôt, qui compile et
construit l'archive.

## Tests

`tests/test_titrages.lua` éprouve le moteur Lua hors de TeX (`texlua
test_titrages.lua`) ; `tests/titrages-test-lua.tex` et `tests/titrages-test-pgf.tex`
compilent les tracés des deux moteurs ; `tests/titrages-test-nombres.tex`, compilé
par lualatex et par pdflatex, écrit les nombres des deux moteurs, que
`tests/test_latex.py` du paquet Python compare à son calcul. Le paquet Python
`titrages`, dans le même dépôt, implémente le même calcul avec une base de couples
et d'indicateurs.

Licence MIT. Copyright (c) 2026 antnardo.
Dépôt : <https://github.com/antnardo/titrages>.
