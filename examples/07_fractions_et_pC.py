"""Diagramme de distribution de H2CO3, et pH = f(pC) d'un acide seul dans les deux sens.

Le sens direct donne C(pH) par ``concentration_pour_pH`` ; le sens inverse, pH(C) par
``Solution.pH`` ; les asymptotes sont celles des réactions prépondérantes du cours.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from titrages import Couple, Solution, concentration_pour_pH
from titrages.trace import tracer_fractions

carbonate = Couple("acide carbonique", (6.3, 10.3), ("H2CO3", "HCO3-", "CO3^2-"))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
tracer_fractions(carbonate, ax1, pH=(1, 13))
ax1.grid(True)

# pH = f(pC) d'un monoacide de pKa 4 : sens direct, puis sens inverse
pKa = 4.0
acide = Couple.polyacide(pKa)
pH = np.linspace(1, 13, 200)
pC = -np.log10(concentration_pour_pH(acide.acide, pH))
ax2.plot(pC, pH, label="sens direct : C(pH)")
pC_inverse = np.linspace(0.5, 8, 30)
pH_inverse = [Solution().avec(acide, 10.0**-p).pH() for p in pC_inverse]
ax2.plot(pC_inverse, pH_inverse, "x", label="sens inverse : pH(C)")
ax2.plot(pC, (pC + pKa) / 2, "r--", lw=0.5, label=r"$\frac{1}{2}(\mathrm{p}C + \mathrm{p}K_a)$")
ax2.plot(pC, pC, "r:", lw=0.5, label="acide fort")
ax2.axhline(7, color="r", lw=0.5)
ax2.set(xlim=(0, 8), ylim=(0, 8), xlabel="p$C$", ylabel="pH")
ax2.legend()
ax2.grid(True)
fig.tight_layout()
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
