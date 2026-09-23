"""Un mélange d'acide chlorhydrique et d'acide éthanoïque, titré par l'ammoniac."""

import matplotlib.pyplot as plt

from titrages import Titrage, solution
from titrages.trace import tracer_derivee, tracer_sauts, tracer_titrage

t = Titrage(solution({"HCl": 0.05, "CH3COOH": 0.05}), solution({"NH3": 0.1}), volume_becher=20)

fig, (haut, bas) = plt.subplots(2, 1, sharex=True, figsize=(6, 6))
tracer_titrage(t, 30, haut)
tracer_sauts(t, 30, haut)
tracer_derivee(t, 30, bas)
bas.set_yscale("log")
bas.set_ylabel("dpH/dV")
plt.show()
