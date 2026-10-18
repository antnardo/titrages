"""Ce que partagent les titrages dont le pH se calcule : l'état à volume fixé, sa
résolution, et la courbe paramétrée par le volume.

**L'état** du mélange, à volume versé fixé, tient en deux nombres : le pH et
:math:`u`, le logarithme de la concentration libre de la composante principale — le
ligand en complexométrie, l'ion commun en précipitation —, plus l'ensemble ``S`` des
solides présents (vide en complexométrie). Les bilans en sont le gradient d'une
fonction strictement convexe : l'équilibre est unique.

**La résolution** a deux voies. La dichotomie imbriquée (pH à l'extérieur, :math:`u` à
l'intérieur) est sûre et n'a besoin d'aucun point de départ ; chaque moteur l'écrit.
Newton, à jacobienne exacte, part d'un état voisin et garde son ensemble ``S`` : il
converge en quelques pas. Un point où il n'a pas convergé, ou dont l'ensemble ``S`` ne
convient pas (un solide présent en quantité négative, un absent sursaturé), est repris
par la dichotomie.

**La courbe** est paramétrée par le volume, la seule coordonnée dont on sache qu'elle est
monotone. Des points s'ajoutent là où elle tourne, jusqu'à ce que chaque segment, en
coordonnées normalisées, soit court ; les points d'inflexion sont les minima de
:math:`|\\mathrm{d}V/\\mathrm{d}y|`, cherchés à l'intérieur des morceaux lisses.
"""

from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from titrages.equilibre import minimum_lisse

__all__ = ["Etat", "MoteurLibre"]

Tableau = npt.NDArray[np.float64]
Mesure = Callable[["Etat"], Tableau]

NEWTON = 50


@dataclass(frozen=True, slots=True)
class Etat:
    """L'équilibre du mélange en chaque volume : pH, :math:`u`, et l'ensemble ``S`` des
    solides présents (un masque de bits, 0 sans solide)."""

    volume: Tableau
    pH: Tableau
    u: Tableau
    S: npt.NDArray[np.int64]

    @classmethod
    def sans_solide(cls, volume: Tableau, pH: Tableau, u: Tableau) -> "Etat":
        return cls(volume, pH, u, np.zeros(np.shape(volume), dtype=np.int64))

    def __getitem__(self, i: "slice | npt.NDArray[np.bool_]") -> "Etat":
        return Etat(self.volume[i], self.pH[i], self.u[i], self.S[i])

    def __len__(self) -> int:
        return len(self.volume)

    def reuni(self, autre: "Etat") -> "Etat":
        """Les deux ensembles de points, rangés par volume croissant."""
        V = np.concatenate((self.volume, autre.volume))
        ordre = np.argsort(V, kind="stable")
        return Etat(
            V[ordre],
            np.concatenate((self.pH, autre.pH))[ordre],
            np.concatenate((self.u, autre.u))[ordre],
            np.concatenate((self.S, autre.S))[ordre],
        )


class MoteurLibre(ABC):
    """La résolution et la courbe ; chaque moteur fournit ses bilans."""

    volume_becher: float

    @abstractmethod
    def _dichotomie(self, V: Tableau) -> Etat:
        """L'équilibre sûr, par dichotomies imbriquées."""

    @abstractmethod
    def _pas_newton(self, e: Etat) -> tuple[Tableau, Tableau]:
        """Un pas de Newton :math:`(\\Delta\\mathrm{pH}, \\Delta u)`, à ensemble ``S`` fixé."""

    def _coherent(self, e: Etat) -> npt.NDArray[np.bool_]:
        """Où l'ensemble ``S`` convient ; sans solide, partout."""
        return np.ones(len(e), dtype=bool)

    # --- la résolution -----------------------------------------------------------------

    def continuer(self, e: Etat) -> tuple[Etat, npt.NDArray[np.bool_]]:
        """Newton depuis ``e``, à ensemble ``S`` fixé ; rend aussi où il a convergé.

        Il a convergé quand le pas tombe sous 10⁻¹³, ou quand il cesse de décroître
        alors qu'il est déjà petit : près d'une équivalence, un bilan est mal
        conditionné, et :math:`u` ne se fixe qu'à quelques 10⁻¹³ près, pour Newton comme
        pour la dichotomie.
        """
        pH, u = e.pH.copy(), e.u.copy()
        converge = np.zeros(len(e), dtype=bool)
        precedent = np.full(len(e), np.inf)
        for _ in range(NEWTON):
            dp, du = self._pas_newton(Etat(e.volume, pH, u, e.S))
            # pas bornés : un Newton qui part de loin ne sort pas du domaine
            dp, du = np.clip(dp, -2.0, 2.0), np.clip(du, -2.0, 2.0)
            libres = ~converge
            pH, u = np.where(libres, pH + dp, pH), np.where(libres, u + du, u)
            pas = np.maximum(np.abs(dp), np.abs(du))
            converge |= (pas < 1e-13) | ((pas < 1e-8) & (pas >= precedent / 2))
            precedent = pas
            if np.all(converge):
                break
        return Etat(e.volume, pH, u, e.S), converge & np.isfinite(pH) & np.isfinite(u)

    def etat(self, volume: npt.ArrayLike, depart: Etat | None = None) -> Etat:
        """L'équilibre en chaque volume. Avec ``depart`` (pH, :math:`u` et ``S`` voisins),
        Newton, puis la dichotomie là où il n'a pas convergé ou que ``S`` ne convient
        pas ; sans, la dichotomie partout."""
        V = np.atleast_1d(np.asarray(volume, dtype=np.float64))
        if np.any(V < 0):
            raise ValueError("un volume versé est positif")
        if depart is None:
            return self._dichotomie(V)
        e, ok = self.continuer(Etat(V, depart.pH, depart.u, depart.S))
        ok &= self._coherent(e)
        if np.all(ok):
            return e
        reprise = self._dichotomie(V[~ok])
        pH, u, S = e.pH.copy(), e.u.copy(), e.S.copy()
        pH[~ok], u[~ok], S[~ok] = reprise.pH, reprise.u, reprise.S
        return Etat(V, pH, u, S)

    # --- la courbe ----------------------------------------------------------------------

    def raffiner(self, volume_max: float, points: int, depart: float, mesure: Mesure) -> Etat:
        """Des points de plus là où la courbe tourne, jusqu'à ce que chaque segment, en
        coordonnées normalisées, soit plus court que ``1/points``."""
        e = self.etat(np.linspace(depart, volume_max, 200))
        for _ in range(80):
            y = mesure(e)
            etendue = max(float(np.ptp(y[np.isfinite(y)])), 1e-12)
            longueur = np.hypot(np.diff(e.volume) / volume_max, np.diff(y) / etendue)
            a_couper = (longueur > 1.0 / points) & (np.diff(e.volume) > 1e-13 * volume_max)
            if not np.any(a_couper):
                break
            i = np.flatnonzero(a_couper)
            milieu = Etat(
                (e.volume[i] + e.volume[i + 1]) / 2,
                (e.pH[i] + e.pH[i + 1]) / 2,
                (e.u[i] + e.u[i + 1]) / 2,
                e.S[i],
            )
            e = e.reuni(self.etat(milieu.volume, milieu))
        return e

    def courbe(self, volume_max: float, points: int, depart: float, mesure: Mesure) -> Etat:
        """``points`` points exacts, à longueur d'arc constante en coordonnées normalisées."""
        e = self.raffiner(volume_max, points, depart, mesure)
        y = mesure(e)
        etendue = max(float(np.ptp(y[np.isfinite(y)])), 1e-12)
        pas = np.hypot(np.diff(e.volume) / volume_max, np.diff(y) / etendue)
        s = np.concatenate(([0.0], np.cumsum(np.nan_to_num(pas))))
        cibles = np.linspace(0.0, s[-1], points)
        j = np.clip(np.searchsorted(s, cibles), 1, len(s) - 1)
        ecart = s[j] - s[j - 1]
        t = np.where(ecart > 0, (cibles - s[j - 1]) / np.where(ecart > 0, ecart, 1.0), 0.0)
        V = e.volume[j - 1] + t * (e.volume[j] - e.volume[j - 1])
        voisin = Etat(
            V, e.pH[j - 1] + t * (e.pH[j] - e.pH[j - 1]), e.u[j - 1] + t * (e.u[j] - e.u[j - 1]), e.S[j - 1]
        )
        return self.etat(V, voisin)

    def sauts(
        self, volume_max: float, points: int, depart: float, mesure: Mesure, pente: Mesure
    ) -> list[tuple[float, Etat]]:
        """Les points d'inflexion : les minima de :math:`|\\mathrm{d}V/\\mathrm{d}y|`, en
        :math:`V`, à l'intérieur des morceaux lisses (même ``S`` de part et d'autre),
        affinés par :func:`~titrages.equilibre.minimum_lisse` avec un pas proportionnel
        à la largeur locale de la grille raffinée. La recherche part d'un millionième de
        ``volume_max`` : avant, la courbe est plate à la précision de la machine."""
        e = self.raffiner(volume_max, points, depart, mesure)
        e = e[e.volume >= 1e-6 * volume_max]
        with np.errstate(divide="ignore"):
            d = np.abs(1 / pente(e))
        lisse = (e.S[1:-1] == e.S[:-2]) & (e.S[1:-1] == e.S[2:])
        resultat = []
        for i in np.flatnonzero(lisse & (d[1:-1] < d[:-2]) & (d[1:-1] <= d[2:])) + 1:
            proche = e[i : i + 1]

            def f(V: float, proche: Etat = proche) -> float:
                x, _ = self.continuer(Etat(np.array([V]), proche.pH, proche.u, proche.S))
                return float(np.abs(1 / pente(x))[0])

            largeur = min(e.volume[i + 1] - e.volume[i], e.volume[i] - e.volume[i - 1])
            V = minimum_lisse(f, e.volume[i - 1], e.volume[i + 1], 1e-3 * largeur)
            if 0 < V < volume_max:
                x, _ = self.continuer(Etat(np.array([V]), proche.pH, proche.u, proche.S))
                resultat.append((V, x))
        return resultat
