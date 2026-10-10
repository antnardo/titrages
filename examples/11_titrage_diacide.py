"""Un diacide titré par la soude : dilution, demi-équivalences, tangentes, pentes.

Pour trois concentrations, la courbe exacte (avec dilution) et celle qui la néglige, le
pH à chaque demi-équivalence, la méthode des tangentes autour des deux équivalences, et
quelques tangentes tracées à partir de la pente analytique ``Courbe.pente``.
"""

import matplotlib.pyplot as plt
import numpy as np

from titrages import Couple, Solution, Titrage
from titrages.trace import tracer_tangentes, tracer_titrage

pKa1, pKa2 = 4.0, 10.0
diacide = Couple.polyacide(pKa1, pKa2)
CB, V0 = 0.1, 10.0

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
for pKa in (pKa1, pKa2):
    ax1.axhline(pKa, color="r", ls="--", lw=0.5)
for C in (0.01, 0.1, 1.0):
    becher = Solution().avec(diacide, C)
    t = Titrage(becher, Solution().avec_base_forte(CB), V0)
    Veq = t.volumes_equivalents()[0]
    (ligne,) = tracer_titrage(t, 3 * Veq, ax1, abscisse=Veq, label=f"$C$ = {C:g} mol/L").lines[-1:]
    couleur = ligne.get_color()
    sans = Titrage(becher, t.burette, V0, dilution=False)
    tracer_titrage(sans, 3 * Veq, ax1, abscisse=Veq, ls="--", lw=0.7, color=couleur)
    for n in (0.5, 1.0, 1.5, 2.0, 2.5):
        pH = float(t.pH(n * Veq))
        ax1.plot([0, n], [pH, pH], "o-", lw=0.5, ms=3, color=couleur)
    for avant in (0.6, 1.6):
        tg = t.tangentes(avant * Veq, 3 * Veq)
        n_eq = avant + 0.4
        print(f"C = {C:g}, n = {n_eq:.0f} : tangentes, n = {tg.volume / Veq:.4f}, pH = {tg.pH:.3f}")
ax1.set(xlim=(0, 3), ylim=(0, 14), xlabel=r"$n = V/V_\mathrm{éq}$")
ax1.legend()
ax1.grid(True)

# Tangentes à partir de la pente exacte, pour vérifier la dérivée
t = Titrage(Solution().avec(diacide, 0.1), Solution().avec_base_forte(CB), V0)
c = t.courbe(30, points=60)
ax2.plot(c.volume, c.pH, "k", lw=1)
for V, pH, s in zip(c.volume, c.pH, c.pente, strict=True):
    dV = 1.5 / np.hypot(1, s)
    ax2.plot([V - dV, V + dV], [pH - s * dV, pH + s * dV], color="C0", lw=0.4)
tracer_tangentes(t, 6, 15, ax2)
ax2.set(xlim=(0, 30), ylim=(0, 14), xlabel="$V$ (mL)", ylabel="pH")
ax2.grid(True)
fig.tight_layout()
plt.show()
