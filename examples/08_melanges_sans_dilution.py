"""Un acide faible, ou un mélange, additionné d'acide ou de base forts, sans dilution.

``Solution.residu(pH) / C`` est le nombre de bases fortes par acide qu'il faut ajouter,
sans dilution, pour atteindre ce pH ; négatif, c'est un acide fort qu'on ajoute.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from titrages import Solution, Titrage, couple

pH = np.linspace(0, 14, 1000)
carbonate = couple("H2CO3")
phosphate = couple("H3PO4")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
for C in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0):
    s = Solution().avec(carbonate, C)
    ax1.plot(s.residu(pH) / C, pH, label=f"$C$ = {C:.0e} mol/L")
ax1.set(xlim=(-2, 4), xlabel=r"$n = C_B/C$ ($n<0$ : acide fort)", ylabel="pH")
ax1.legend(fontsize="small")
ax1.grid(True)

melange = Solution().avec(carbonate, 0.1).avec(phosphate, 0.1)
ax2.plot(melange.residu(pH) / 0.2, pH)
ax2.set(xlim=(-3, 3), xlabel=r"$n = C_B/(C_1 + C_2)$", ylabel="pH")
ax2.set_title("carbonate + phosphate, 0,1 mol/L chacun")
ax2.grid(True)
fig.tight_layout()

# Le pH à la demi-équivalence et à l'équivalence, sans dilution
t = Titrage(Solution().avec(carbonate, 0.1), Solution().avec_base_forte(0.1), 10, dilution=False)
print("pH à n = 1/2 et n = 1, sans dilution :", [round(float(p), 3) for p in t.pH([5, 10])])
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
