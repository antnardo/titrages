"""Le dosage de Liebig : les cyanures par l'argent, pH calculé.

Ag⁺ fixe deux cyanures, Ag(CN)₂⁻ ; la base ne donne pas de constante pour AgCN dissous,
qui ne se forme donc pas, et donne Ag(CN)₃²⁻, peu stable. Les équivalences se comptent
comme celles d'un polyacide, une par complexe prédominant : Ag(CN)₃²⁻ puis Ag(CN)₂⁻ ;
seule la seconde donne un saut. Au-delà, le précipité AgCN, qui marque la fin du dosage
au laboratoire, n'est pas modélisé : complexation et précipitation ne se combinent pas.
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from titrages import SolutionComplexation, TitrageComplexometrique, complexe, ligand
from titrages.trace import tracer_titrage_complexometrique

with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)  # HCN est un acide faible : le pH se calcule
    t = TitrageComplexometrique(
        SolutionComplexation(ligand=0.1),  # KCN 0,1 mol/L
        SolutionComplexation({"Ag+": 0.05}),  # AgNO3 0,05 mol/L
        20,
        [complexe("Ag+", "cyanure")],
        ligand("cyanure"),
    )
print("log beta :", complexe("Ag+", "cyanure").log_beta)
print("équivalences :", [round(v, 3) for v in t.volumes_equivalents()], "mL")
print("sauts :", [(round(s.volume, 3), round(s.pM, 2)) for s in t.sauts(30)])

fig, ax = plt.subplots(figsize=(6, 4))
tracer_titrage_complexometrique(t, 30, ax, lw=2, sauts=True, label="pAg")
ax.set_ylabel("pAg")
ax.set_ylim(0, 25)
droite = ax.twinx()
c = t.courbe(30)
droite.plot(c.volume, c.pH, "C1--", lw=1.5, label="pH calculé")
droite.set_ylabel("pH")
for v in t.volumes_equivalents():
    ax.axvline(v, color="0.6", lw=0.8, ls=":")
ax.set_title("KCN 0,1 mol/L par AgNO₃ 0,05 mol/L")
ax.legend(loc="center left")
droite.legend(loc="center right")
ax.grid(True)
fig.tight_layout()
print(f"pH : {c.pH[0]:.2f} au départ, {float(t.pH_du_melange(20)):.2f} à l'équivalence")
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
