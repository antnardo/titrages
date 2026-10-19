"""La complexométrie à l'EDTA : l'effet du pH, et la dureté totale d'une eau (Ca²⁺ et Mg²⁺).

Les constantes viennent de la base du paquet (Harvey, d'après Martell et Smith). L'EDTA
s'introduit sous la forme qu'on pèse, Na₂H₂Y : sans tampon, il libère ses protons en se
fixant, et le pH calculé s'effondre.
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt

from titrages import Solution, SolutionComplexation, TitrageComplexometrique, complexe, couple, ligand

EDTA = ligand("EDTA")
CA, MG = complexe("Ca2+", "EDTA"), complexe("Mg2+", "EDTA")
BECHER, BURETTE = SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01, protons_ligand=2)

fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4))
for pH in (6, 8, 10, 12):
    t = TitrageComplexometrique(BECHER, BURETTE, 20, [CA], EDTA, pH=pH)
    c = t.courbe(40)
    a.plot(c.volume, c.pM, lw=1.5, label=f"pH {pH} imposé")
    print(
        f"pH {pH:2d} : pCa à l'équivalence {t.pM(20):.2f}, sauts {[round(s.volume, 3) for s in t.sauts(40)]}"
    )

nh = couple("NH3")
tampon = Solution().avec(nh.forme(0), 0.1).avec(nh.forme(1), 0.1)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)  # le pH se calcule : c'est voulu
    sans = TitrageComplexometrique(BECHER, BURETTE, 20, [CA], EDTA)
    avec = TitrageComplexometrique(
        SolutionComplexation({"Ca2+": 0.01}, acidobasique=tampon),
        SolutionComplexation(ligand=0.01, protons_ligand=2, acidobasique=tampon),
        20,
        [CA],
        EDTA,
    )
for t, nom, style in ((sans, "sans tampon", "k:"), (avec, "tampon NH₄⁺/NH₃ 0,1 mol/L", "k--")):
    c = t.courbe(40)
    a.plot(c.volume, c.pM, style, lw=1.5, label=f"{nom}, pH calculé")
    sauts = [round(s.volume, 3) for s in t.sauts(40)]
    print(f"{nom:26s} : pH {c.pH[0]:.2f} → {float(t.pH_du_melange(20)):.2f} à l'équivalence, sauts {sauts}")
a.set_xlabel("$V$ versé")
a.set_ylabel("pCa")
a.set_title("Ca²⁺ 0,01 mol/L par Na₂H₂Y 0,01 mol/L")
a.legend(fontsize="small", loc="lower right")
a.grid(True)

eau = SolutionComplexation({"Ca2+": 0.006, "Mg2+": 0.004})
for metal, couleur in (("Ca2+", "C0"), ("Mg2+", "C1")):
    t = TitrageComplexometrique(
        eau, SolutionComplexation(ligand=0.01), 20, [CA, MG], EDTA, pH=10, metal=metal
    )
    c = t.courbe(40)
    b.plot(c.volume, c.pM, color=couleur, lw=2, label=f"p{metal[:-2]}")
    for s in t.sauts(40):
        b.plot([s.volume], [s.pM], "o", color=couleur, ms=4)
print("dureté : équivalences", [round(v, 3) for v in t.volumes_equivalents()], "mL")
b.set_xlabel("$V$ versé")
b.set_ylabel("pM")
b.set_title("Ca²⁺ 6 mmol/L + Mg²⁺ 4 mmol/L par l'EDTA, pH 10")
b.legend()
b.grid(True)
fig.tight_layout()
fig.savefig(Path(__file__).with_suffix(".pdf"))
plt.show()
