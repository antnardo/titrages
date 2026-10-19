"""Un mélange d'iodure, de bromure et de chlorure titré par l'argent : trois précipités successifs."""

from pathlib import Path

import matplotlib.pyplot as plt

from titrages import SolutionPrecipitation, TitrageParPrecipitation, precipite
from titrages.trace import tracer_sauts, tracer_titrage_precipitation

t = TitrageParPrecipitation(
    SolutionPrecipitation(partenaires={"I-": 0.02, "Br-": 0.02, "Cl-": 0.02}),
    SolutionPrecipitation(commun=0.05),
    20,
    [precipite("AgI"), precipite("AgBr"), precipite("AgCl")],
)
print("équivalences :", [round(v, 3) for v in t.volumes_equivalents()], "mL")
for e in t.evenements(30):
    print(f"Ag{e.partenaire[:-1]:2s} apparaît à {e.volume:.4f} mL, pAg {e.p:.2f}")
print("sauts :", [(round(s.volume, 3), round(s.p, 2)) for s in t.sauts(30)])

fig, ax = plt.subplots(figsize=(6, 4))
tracer_titrage_precipitation(t, 30, ax, lw=2, label="pAg")
tracer_sauts(t, 30, ax, texte=False)
c = t.courbe(30)
for x, couleur in (("I-", "C1"), ("Br-", "C2"), ("Cl-", "C3")):
    ax.plot(c.volume, t.p_partenaire(c.volume, x), lw=1, color=couleur, label=f"p{x[:-1]}")
ax.set_ylabel("pAg, pX")
ax.set_title("I⁻, Br⁻, Cl⁻ 0,02 mol/L par Ag⁺ 0,05 mol/L")
ax.legend()
ax.grid(True)
fig.tight_layout()
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
