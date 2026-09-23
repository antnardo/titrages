"""Trois polyacides titrés par la soude : toutes les équivalences ne donnent pas un saut."""

import matplotlib.pyplot as plt

from titrages import Titrage, solution
from titrages.trace import tracer_titrage

fig, ax = plt.subplots(figsize=(6, 4))
for nom in ("H2CO3", "H3PO4", "acide citrique"):
    t = Titrage(solution({nom: 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
    equivalences = t.volumes_equivalents()
    sauts = t.sauts(35)
    print(f"{nom:15s} équivalences {equivalences}  sauts {[round(s.volume, 2) for s in sauts]}")
    tracer_titrage(t, 35, ax, label=nom)
ax.legend()
ax.grid(True)
plt.show()
