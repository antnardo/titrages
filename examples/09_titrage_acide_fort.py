"""Un acide fort par une base forte : calcul exact, formules du cours, méthode des tangentes.

Le calcul exact (dilution comprise) coïncide avec pH = pC_B + p((1 − n)/(y + n)) avant
l'équivalence et pH = pK_e − pC_B − p((n − 1)/(y + n)) après, où y = C_B/C_A. La méthode
des tangentes est ensuite construite pour plusieurs points de départ : son erreur dépend
de la pente choisie.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from titrages import Titrage, solution
from titrages.trace import tracer_tangentes, tracer_titrage

CA, CB, V0 = 0.1, 0.4, 10.0
y = CB / CA
t = Titrage(solution({"HCl": CA}), solution({"NaOH": CB}), V0)
(Veq,) = t.volumes_equivalents()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
tracer_titrage(t, 2 * Veq, ax1, abscisse=Veq, lw=3, label="titrages (exact)")
n1 = np.linspace(0, 0.999, 200)
n2 = np.linspace(1.001, 2, 200)
ax1.plot(n1, -np.log10(CB) - np.log10((1 - n1) / (y + n1)), "k--", label="formules du cours")
ax1.plot(n2, 14 + np.log10(CB) + np.log10((n2 - 1) / (y + n2)), "k--")
sans_dilution = Titrage(t.becher, t.burette, V0, dilution=False)
tracer_titrage(sans_dilution, 2 * Veq, ax1, abscisse=Veq, lw=0.8, label="sans dilution")
ax1.legend()
ax1.grid(True)

# Erreur de la méthode des tangentes selon le point de départ : le volume trouvé est
# juste (la courbe est presque symétrique), le pH l'est d'autant moins que la pente
# choisie est faible — y > 1 rend la courbe dissymétrique.
departs = np.linspace(0.3, 0.95, 14) * Veq
ecarts_V, ecarts_pH, pentes = [], [], []
for V1 in departs:
    tg = t.tangentes(V1, 2 * Veq)
    ecarts_V.append((tg.volume - Veq) / Veq)
    ecarts_pH.append(7 - tg.pH)
    pentes.append(tg.pente)
    print(f"V1 = {V1:5.2f} mL  pente {tg.pente:7.3f}  V = {tg.volume:.6f} mL  pH = {tg.pH:.3f}")
ax2.semilogx(pentes, ecarts_pH, "o-", label="7 - pH")
ax2.semilogx(pentes, 100 * np.array(ecarts_V), "s-", label=r"$(V - V_\mathrm{éq})/V_\mathrm{éq}$ (%)")
ax2.set(xlabel="pente de la tangente (pH/mL)", ylabel="écart du point construit")
ax2.legend()
ax2.grid(True)
tracer_tangentes(t, 0.8 * Veq, 2 * Veq, ax1, abscisse=Veq)
fig.tight_layout()
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
