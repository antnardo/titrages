"""Écrire une courbe pour la reprendre ailleurs : pgfplots, CSV.

Le paquet LaTeX ``titrages`` calcule ses courbes lui-même. Ce module sert l'autre
cas : une figure composée en pgfplots à partir d'un calcul Python — un mélange tiré
de la base, par exemple, que le paquet LaTeX ne connaît pas.
"""

import csv
from pathlib import Path

import numpy as np
import numpy.typing as npt

from titrages.modeles import Courbe, CourbeConductimetrique, CourbeRedox

__all__ = ["vers_csv", "vers_pgfplots"]

Tableau = npt.NDArray[np.float64]
Toutes = Courbe | CourbeRedox | CourbeConductimetrique


def _colonnes(courbe: Toutes) -> dict[str, Tableau]:
    """Les colonnes d'une courbe, le volume d'abord et la grandeur mesurée ensuite."""
    if isinstance(courbe, Courbe):
        return {"volume": courbe.volume, "pH": courbe.pH, "pente": courbe.pente}
    if isinstance(courbe, CourbeRedox):
        return {"volume": courbe.volume, "E": courbe.E, "pente": courbe.pente}
    return {"volume": courbe.volume, "sigma": courbe.conductivite, "pH": courbe.pH}


def vers_pgfplots(courbe: Toutes, chiffres: int = 5) -> str:
    """La liste ``(V,y)`` à mettre dans ``\\addplot coordinates {…}`` : pH, E ou σ."""
    volume, y = list(_colonnes(courbe).values())[:2]
    return " ".join(f"({v:.{chiffres}g},{p:.{chiffres}g})" for v, p in zip(volume, y, strict=True))


def vers_csv(courbe: Toutes, chemin: str | Path) -> Path:
    """Une table ``volume,pH,pente`` (ou ``E``, ou ``sigma``), que lit aussi ``\\addplot table``."""
    chemin = Path(chemin)
    colonnes = _colonnes(courbe)
    with chemin.open("w", newline="", encoding="utf-8") as flux:
        ecrivain = csv.writer(flux)
        ecrivain.writerow(list(colonnes))
        for ligne in zip(*colonnes.values(), strict=True):
            ecrivain.writerow([f"{x:.8g}" for x in ligne])
    return chemin
