"""La précipitation d'un hydroxyde : Mg²⁺ en milieu acide, titré par la soude, pH calculé.

La soude neutralise d'abord l'acide fort, puis précipite Mg(OH)₂ : le pH reste sur un
palier, fixé par le produit de solubilité, tant que du magnésium reste en solution.
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from titrages import Solution, SolutionPrecipitation, Titrage, TitrageParPrecipitation, precipite, solution

with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)  # un hydroxyde : le pH se calcule
    mg = TitrageParPrecipitation(
        SolutionPrecipitation(commun=0.01, acidobasique=Solution(acide_fort=0.002)),  # MgCl2 + HCl
        SolutionPrecipitation(acidobasique=Solution(base_forte=0.1)),  # NaOH
        20,
        [precipite("Mg(OH)2")],
    )
seul = Titrage(solution({"HCl": 0.002}), solution({"NaOH": 0.1}), 20)
print("équivalences :", [round(v, 3) for v in mg.volumes_equivalents()], "mL")
for e in mg.evenements(8):
    quoi = "apparaît" if e.apparition else "disparaît"
    print(f"Mg(OH)2 {quoi} à {e.volume:.4f} mL, pH {float(mg.pH_du_melange(e.volume)):.2f}")
print("pH du palier :", [round(float(mg.pH_du_melange(V)), 2) for V in (1, 2, 3, 4)])

fig, ax = plt.subplots(figsize=(6, 4))
c = mg.courbe(8)
ax.plot(c.volume, c.pH, lw=2, label="HCl 2 mmol/L + Mg²⁺ 10 mmol/L")
s = seul.courbe(8)
ax.plot(s.volume, s.pH, "--", lw=1.5, label="HCl 2 mmol/L seul")
for e in mg.evenements(8):
    ax.plot([e.volume], [float(mg.pH_du_melange(e.volume))], "kD", ms=5)
for v in mg.volumes_equivalents():
    ax.axvline(v, color="0.6", lw=0.8, ls=":")
ax.set_xlabel("$V$ versé (NaOH 0,1 mol/L)")
ax.set_ylabel("pH")
ax.set_ylim(0, 14)
ax.legend(loc="lower right")
ax.grid(True)
fig.tight_layout()
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
