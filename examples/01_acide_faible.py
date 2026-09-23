"""Le titrage de l'acide éthanoïque par la soude : courbe, saut, tangentes, indicateur."""

import matplotlib.pyplot as plt

from titrages import Titrage, indicateur, solution
from titrages.trace import tracer_sauts, tracer_tangentes, tracer_titrage, tracer_zone_virage

t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)

fig, ax = plt.subplots(figsize=(6, 4))
tracer_titrage(t, 20, ax, lw=2)
tracer_sauts(t, 20, ax)
tracer_tangentes(t, 7, 20, ax)
tracer_zone_virage(indicateur("phénolphtaléine"), ax)
ax.set_ylim(0, 14)
ax.grid(True)
ax.legend()
fig.tight_layout()

tg = t.tangentes(7, 20)
print(f"équivalence : {t.volumes_equivalents()[0]:.3f} mL, tangentes : {tg.volume:.3f} mL")
plt.show()
