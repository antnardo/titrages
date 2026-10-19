"""La méthode de Mohr : les chlorures par l'argent, le chromate pour indicateur."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from titrages import SolutionPrecipitation, TitrageParPrecipitation, precipite
from titrages.trace import tracer_titrage_precipitation

AGCL, AG2CRO4 = precipite("AgCl"), precipite("Ag2CrO4")


def mohr(chromate: float, pH: float) -> TitrageParPrecipitation:
    """20 mL de chlorure à 0,05 mol/L, et le chromate, par le nitrate d'argent à 0,05 mol/L."""
    return TitrageParPrecipitation(
        SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": chromate}),
        SolutionPrecipitation(commun=0.05),
        20,
        [AGCL, AG2CRO4],
        pH=pH,
    )


def ecart_indicateur(chromate: float, pH: float) -> float:
    """Le volume où apparaît le chromate d'argent, moins l'équivalence du chlorure, en µL."""
    t = mohr(chromate, pH)
    (rouge,) = [e for e in t.evenements(30) if e.partenaire == "CrO4^2-"]
    return 1000 * (rouge.volume - t.volumes_equivalents()[0])


t = mohr(0.005, 7.5)
print("équivalences :", [round(v, 3) for v in t.volumes_equivalents()])
for e in t.evenements(30):
    print(
        f"{e.partenaire:8s} {'apparaît' if e.apparition else 'disparaît'} à {e.volume:.4f} mL, pAg {e.p:.2f}"
    )

fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4))
tracer_titrage_precipitation(t, 30, a, lw=2, label="pAg")
c = t.courbe(30)
a.plot(c.volume, t.p_partenaire(c.volume, "Cl-"), lw=1.5, label="pCl")
a.set_ylabel("pAg, pCl")
a.legend()
a.grid(True)
a.set_title("Cl⁻ 0,05 mol/L + CrO₄²⁻ 0,005 mol/L, pH 7,5")

concentrations = np.logspace(-4, -1.5, 14)
for pH in (7.5, 6.5, 5.5):
    b.semilogx(
        concentrations,
        [ecart_indicateur(C, pH) for C in concentrations],
        "o-",
        ms=3,
        label=f"pH {pH}".replace(".", ","),
    )
b.axhline(0, color="k", lw=0.8)
b.set_xlabel("[CrO₄²⁻] introduit (mol/L)")
b.set_ylabel("V(Ag₂CrO₄) − V_éq (µL)")
b.legend()
b.grid(True, which="both", alpha=0.4)
b.set_title("Erreur de l'indicateur")
fig.tight_layout()

for C in (0.001, 0.005, 0.02):
    neutre, acide = ecart_indicateur(C, 7.5), ecart_indicateur(C, 5.5)
    print(f"[CrO4] = {C} mol/L : {neutre:+.1f} µL à pH 7,5, {acide:+.1f} µL à pH 5,5")
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
