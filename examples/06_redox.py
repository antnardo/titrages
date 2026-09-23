"""Potentiométrie : le vanadium par le permanganate (trois sauts), et l'iodométrie."""

import matplotlib.pyplot as plt

from titrages import TitrageRedox, solution_redox
from titrages.trace import tracer_titrage_redox

vanadium = TitrageRedox(solution_redox({"V2+": 0.1}), solution_redox({"KMnO4": 0.02}), 10, "ECS")
iode = TitrageRedox(solution_redox({"I2": 0.05}), solution_redox({"Na2S2O3": 0.1}), 10, "ECS")
print("vanadium :", vanadium.volumes_equivalents(), [round(s.E, 3) for s in vanadium.sauts(40)])
print("iodométrie :", iode.volumes_equivalents())

fig, ax = plt.subplots(figsize=(6, 4))
tracer_titrage_redox(vanadium, 40, ax, sauts=True, label="V²⁺ par MnO₄⁻, pH 0")
tracer_titrage_redox(iode, 40, ax, label="I₂ par S₂O₃²⁻")
ax.legend()
ax.grid(True)
plt.show()
