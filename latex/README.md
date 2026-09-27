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
  tangentes, et le pH à un volume donné, en nombres développables.
- Sous pdfLaTeX, un moteur pgfmath trace la courbe et les fractions.

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
\end{axis}
\end{tikzpicture}
```

Documentation : `titrages-doc.pdf` (en français). Pour la recompiler :

```sh
lualatex titrages-doc
makeindex -s gind.ist -o titrages-doc.ind titrages-doc.idx
lualatex titrages-doc
lualatex titrages-doc
```

Installation : `titrages.sty` et `titrages.lua` dans un même dossier du chemin de
TeX, par exemple `TEXMFHOME/tex/latex/titrages/`. Dépend de pgfplots.

Le paquet Python `titrages` (PyPI) implémente le même calcul, avec une base de
couples et d'indicateurs.

Licence MIT. Copyright (c) 2026 antnardo.
Dépôt : https://github.com/antnardo/titrages
