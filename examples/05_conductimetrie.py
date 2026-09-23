"""Conductimétrie : un acide fort, un acide faible, et leur mélange, par la soude."""

import matplotlib.pyplot as plt

from titrages import Titrage, solution
from titrages.trace import tracer_conductimetrie

fig, ax = plt.subplots(figsize=(6, 4))
for nom, becher in [
    ("HCl", {"HCl": 0.01}),
    ("CH3COOH", {"CH3COOH": 0.01}),
    ("HCl + CH3COOH", {"HCl": 0.01, "CH3COOH": 0.01}),
]:
    t = Titrage(solution(becher), solution({"NaOH": 0.1}), volume_becher=100)
    tracer_conductimetrie(t, 30, ax, corrigee=True, label=nom)
ax.legend()
ax.grid(True)
plt.show()
