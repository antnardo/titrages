"""Une courbe calculée ici, composée par pgfplots : les coordonnées à coller dans un .tex."""

from titrages import Titrage, solution, vers_pgfplots

t = Titrage(solution({"Na2CO3": 0.05}), solution({"HCl": 0.1}), volume_becher=20)
print(r"\addplot coordinates {" + vers_pgfplots(t.courbe(25, 150)) + "};")
