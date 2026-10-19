"""La méthode de Charpentier-Volhard : l'argent en excès dosé en retour par le thiocyanate.

Le chlorure est précipité par un excès connu d'argent, en milieu nitrique ; l'excès est
titré par le thiocyanate, l'ion fer(III) signalant la fin par la couleur rouge du
complexe FeSCN²⁺, qui demande un peu de SCN⁻ libre. Le chlorure d'argent, plus soluble
que le thiocyanate d'argent, se redissout dès que l'argent libre manque : il retient
SCN⁻ après l'équivalence, et le virage traîne. D'où le conseil de filtrer AgCl, ou de
l'isoler (nitrobenzène), avant le titrage.
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from titrages import Solution, SolutionPrecipitation, TitrageParPrecipitation, precipite
from titrages.trace import tracer_titrage_precipitation

AGCL, AGSCN = precipite("AgCl"), precipite("AgSCN")
NITRIQUE = Solution(acide_fort=0.1)  # le milieu nitrique : SCN- y est un peu protoné (pKa 0,9)


def volhard(argent: float, chlorure: float) -> TitrageParPrecipitation:
    """20 mL d'argent et de chlorure, en milieu nitrique, par le thiocyanate à 0,1 mol/L."""
    partenaires = {"Cl-": chlorure} if chlorure else {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # le pH se calcule : c'est voulu
        return TitrageParPrecipitation(
            SolutionPrecipitation(commun=argent, partenaires=partenaires, acidobasique=NITRIQUE),
            SolutionPrecipitation(partenaires={"SCN-": 0.1}),
            20,
            [AGCL, AGSCN] if chlorure else [AGSCN],
        )


avec = volhard(0.1, 0.05)  # Cl- 0,05 mol/L précipité par Ag+ 0,1 mol/L : AgCl reste dans le bécher
sans = volhard(0.05, 0.0)  # AgCl filtré : il ne reste que l'argent en excès, 0,05 mol/L

print("équivalence (excès d'argent) :", [round(v, 3) for v in avec.volumes_equivalents()], "mL")
print(f"pH du milieu : {float(avec.pH_du_melange(10)):.2f}")
for e in avec.evenements(25):
    print(f"{e.partenaire:5s} {'apparaît' if e.apparition else 'disparaît'} à {e.volume:.4f} mL")

fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4))
tracer_titrage_precipitation(avec, 25, a, lw=2, label="AgCl laissé dans le bécher")
a.set_title("Cl⁻ 0,05 + Ag⁺ 0,1 mol/L par SCN⁻ 0,1 mol/L")
a.set_ylabel("pAg")
a.legend()
a.grid(True)

for t, nom in ((avec, "AgCl laissé dans le bécher"), (sans, "AgCl filtré")):
    c = t.courbe(25)
    b.plot(c.volume, t.p_partenaire(c.volume, "SCN-"), lw=2, label=nom)
b.axhline(5, color="C3", ls="--", lw=1)
b.text(0.5, 5.15, "SCN⁻ libre où le rouge de FeSCN²⁺ se voit (≈ 10⁻⁵ mol/L)", color="C3", fontsize="small")
b.set_xlabel("$V$ versé")
b.set_ylabel("pSCN")
b.set_title("Le thiocyanate libre, que voit l'indicateur")
b.legend(loc="lower left")
b.grid(True)
fig.tight_layout()

for t, nom in ((avec, "avec AgCl"), (sans, "AgCl filtré")):
    c = t.courbe(25, 2000)
    rouge = c.volume[t.p_partenaire(c.volume, "SCN-") <= 5][0]
    print(f"{nom:12s} : pSCN = 5 à {rouge:.3f} mL")
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
