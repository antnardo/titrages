"""Une base faible et une dibase titrées par l'acide chlorhydrique.

Une base faible n'est que la forme basique d'un couple : l'ammoniac est la forme 1 de
NH4+/NH3, le carbonate la forme 2 de H2CO3.
"""

import matplotlib.pyplot as plt

from titrages import Titrage, solution
from titrages.trace import tracer_sauts, tracer_titrage

fig, ax = plt.subplots(figsize=(6, 4))
for nom in ("NH3", "Na2CO3"):
    t = Titrage(solution({nom: 0.1}), solution({"HCl": 0.1}), volume_becher=10)
    equivalences = t.volumes_equivalents()
    print(
        nom,
        "équivalences",
        equivalences,
        "pH à la première demi-équivalence",
        round(float(t.pH(equivalences[0] / 2)), 2),
    )
    tracer_titrage(t, 25, ax, label=nom)
    tracer_sauts(t, 25, ax, texte=False)
ax.legend()
ax.grid(True)
plt.show()
