"""La complexométrie à pH calculé : le niveau 2, sans tampon imposé.

**Les inconnues** sont, à volume versé fixé, le pH et :math:`u = \\log l`, le ligand libre
sous sa forme la plus basique. Tout le reste s'en déduit : les formes protonées du
ligand, :math:`[\\mathrm{H}_k\\mathrm{L}] = K_k\\,h^k\\,l`, et les complexes de chaque
métal, dont les fractions ne dépendent que de :math:`l` (le métal est *éliminé*). Deux
bilans les fixent, écrits pour le mélange :

- celui du ligand, :math:`g_L = \\sum_m T_m Y_{1,m}(u) + l P(h) - T_L`, avec
  :math:`P = \\sum_k K_k h^k` ;
- celui des protons, :math:`g_H = R_\\text{ab}(\\mathrm{pH}) + k_0 T_L - l P_1(h)` :
  le résidu de l'électroneutralité de la partie acido-basique des solutions (acides et
  bases forts, couples, l'eau), plus les protons qu'apporte le ligand introduit sous la
  forme :math:`\\mathrm{H}_{k_0}\\mathrm{L}` (Na₂H₂Y : :math:`k_0 = 2`), moins ceux que
  portent ses formes libres, :math:`P_1 = \\sum_k k\\,K_k h^k`.

**L'unicité** vient de la convexité : les bilans sont le gradient d'une fonction
strictement convexe des :math:`\\log` des concentrations libres. À pH fixé, :math:`g_L`
croît avec :math:`u` ; une fois :math:`u` résolu, :math:`g_H` croît avec le pH. Deux
dichotomies imbriquées trouvent donc l'unique équilibre, sans point de départ. Pour
affiner autour d'un point connu — les points ajoutés à la courbe, la recherche des
sauts —, un Newton sur le système 2 × 2, à jacobienne exacte, part du point voisin ; la
racine trouvée est la seule, et la dichotomie reprend la main s'il ne converge pas.

**Les dérivées** viennent du théorème des fonctions implicites : à volume fixé,
:math:`J\\,(\\mathrm{d\\,pH}, \\mathrm{d}u)/\\mathrm{d}V = -\\partial g/\\partial V`, avec
:math:`\\partial g/\\partial V = R_\\text{burette}/(V_0 + V)` sur la courbe.

Le volume n'y est plus une fonction explicite du paramètre : la courbe est paramétrée
par :math:`V` lui-même, la seule coordonnée dont on sache qu'elle est monotone — le pH
ne l'est pas toujours.
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from titrages.equilibre import (
    LN10,
    MARGE_PH,
    log_somme_complexes,
    minimum_lisse,
    moments_complexes,
    racine_croissante,
)

if TYPE_CHECKING:
    from titrages.complexation import SolutionComplexation

__all__ = ["ComplexationLibre"]

Tableau = npt.NDArray[np.float64]

PL_MAX = 80.0
NEWTON = 50


@dataclass(frozen=True, slots=True)
class Etat:
    """L'équilibre du mélange en chaque volume : le pH et :math:`u = \\log l`."""

    volume: Tableau
    pH: Tableau
    u: Tableau


@dataclass(frozen=True, slots=True)
class ComplexationLibre:
    """Le moteur du titrage complexométrique dont le pH se calcule.

    ``log_beta`` porte, pour chaque métal, :math:`(0, \\log\\beta_1, \\dots)` **sans**
    constante conditionnelle ; ``pKa`` sont ceux du ligand.
    """

    becher: "SolutionComplexation"
    burette: "SolutionComplexation"
    volume_becher: float
    log_beta: Mapping[str, Tableau]
    pKa: tuple[float, ...]
    _log_K: Tableau = field(init=False, repr=False)
    _pKe: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        # [H_k L] = K_k h^k l : log K_k est la somme des k plus grands pKa.
        decroissants = sorted(self.pKa, reverse=True)
        log_K = np.concatenate(([0.0], np.cumsum(decroissants))) if decroissants else np.zeros(1)
        object.__setattr__(self, "_log_K", np.asarray(log_K, dtype=np.float64))
        a, b = self.becher.acidobasique.pKe, self.burette.acidobasique.pKe
        if a != b:
            raise ValueError("le bécher et la burette n'ont pas le même pKe")
        object.__setattr__(self, "_pKe", a)
        for s in (self.becher, self.burette):
            if s.protons_ligand > len(self.pKa):
                raise ValueError(f"le ligand n'a que {len(self.pKa)} protons : protons_ligand trop grand")

    # --- les briques -----------------------------------------------------------------

    def _melange(self, V: Tableau, b: float, u: float) -> Tableau:
        """La concentration dans le mélange de ce qui vaut ``b`` au bécher et ``u`` à la burette."""
        V0 = self.volume_becher
        return (V0 * b + V * u) / (V0 + V)

    def _ligand_libre(self, pH: Tableau, u: Tableau) -> tuple[Tableau, Tableau, Tableau]:
        """:math:`lP`, :math:`lP_1`, :math:`lP_2` : le ligand libre, toutes formes, et ses protons."""
        k = np.arange(len(self._log_K), dtype=np.float64).reshape((-1,) + (1,) * pH.ndim)
        termes = np.power(10.0, self._log_K.reshape(k.shape) - k * pH + u)
        return termes.sum(axis=0), (k * termes).sum(axis=0), (k * k * termes).sum(axis=0)

    def _metaux(
        self, V: Tableau, u: Tableau, s: "SolutionComplexation | None" = None
    ) -> tuple[Tableau, Tableau]:
        """:math:`\\sum_m T_m Y_{1,m}` et :math:`\\sum_m T_m \\mathrm{Var}_m`, du mélange ou de ``s``."""
        Y, W = np.zeros(u.shape), np.zeros(u.shape)
        for m, lb in self.log_beta.items():
            T = self._melange(V, self.becher.concentration(m), self.burette.concentration(m))
            if s is not None:
                T = np.full(u.shape, s.concentration(m))
            Y1, var = moments_complexes(lb, -u)
            Y, W = Y + T * Y1, W + T * var
        return Y, W

    def _acidobasique(self, V: Tableau, pH: Tableau) -> tuple[Tableau, Tableau]:
        """Le résidu acido-basique du mélange et sa dérivée en pH (pouvoir tampon)."""
        b, u = self.becher.acidobasique, self.burette.acidobasique
        V0 = self.volume_becher
        R = (V0 * np.asarray(b.residu(pH)) + V * np.asarray(u.residu(pH))) / (V0 + V)
        D = (V0 * np.asarray(b.pouvoir_tampon(pH)) + V * np.asarray(u.pouvoir_tampon(pH))) / (V0 + V)
        return R, D

    def _protons_introduits(self, V: Tableau) -> Tableau:
        b, u = self.becher, self.burette
        return self._melange(V, b.protons_ligand * b.ligand, u.protons_ligand * u.ligand)

    def _bilans(self, V: Tableau, pH: Tableau, u: Tableau) -> tuple[Tableau, Tableau]:
        """:math:`(g_H, g_L)` du mélange."""
        lP, lP1, _ = self._ligand_libre(pH, u)
        Y, _ = self._metaux(V, u)
        T_L = self._melange(V, self.becher.ligand, self.burette.ligand)
        R, _ = self._acidobasique(V, pH)
        return R + self._protons_introduits(V) - lP1, Y + lP - T_L

    def _jacobien(self, V: Tableau, pH: Tableau, u: Tableau) -> tuple[Tableau, Tableau, Tableau, Tableau]:
        """:math:`\\partial(g_H, g_L)/\\partial(\\mathrm{pH}, u)` : symétrique, définie."""
        lP, lP1, lP2 = self._ligand_libre(pH, u)
        _, W = self._metaux(V, u)
        _, D = self._acidobasique(V, pH)
        return D + LN10 * lP2, -LN10 * lP1, -LN10 * lP1, LN10 * (W + lP)

    # --- l'équilibre à volume fixé -----------------------------------------------------

    def _bornes_u(self) -> tuple[float, float]:
        concentrations = [self.becher.ligand, self.burette.ligand]
        for s in (self.becher, self.burette):
            concentrations += list(s.metaux.values())
        haut = math.log10(max(concentrations) + 1e-300) + 1
        bas = -PL_MAX - float(self._log_K[-1]) - (len(self._log_K) - 1) * MARGE_PH
        return bas, haut

    def _u_de(self, V: Tableau, pH: Tableau) -> Tableau:
        """:math:`u` à pH fixé : le bilan du ligand croît avec :math:`u`."""
        bas, haut = self._bornes_u()
        T_L = self._melange(V, self.becher.ligand, self.burette.ligand)

        def g(u: Tableau) -> Tableau:
            Y, _ = self._metaux(V, u)
            return Y + self._ligand_libre(pH, u)[0] - T_L

        return np.asarray(
            racine_croissante(g, np.full(V.shape, bas), np.full(V.shape, haut)), dtype=np.float64
        )

    def _dichotomie(self, V: Tableau) -> Etat:
        """L'équilibre par deux dichotomies imbriquées : sûr, sans point de départ."""

        def g(pH: Tableau) -> Tableau:
            return self._bilans(V, pH, self._u_de(V, pH))[0]

        pH = np.asarray(
            racine_croissante(g, np.full(V.shape, -MARGE_PH), np.full(V.shape, self._pKe + MARGE_PH)),
            dtype=np.float64,
        )
        return Etat(V, pH, self._u_de(V, pH))

    def _newton(self, V: Tableau, pH: Tableau, u: Tableau) -> tuple[Tableau, Tableau, npt.NDArray[np.bool_]]:
        """Newton depuis ``(pH, u)`` ; rend aussi où il a convergé.

        Il a convergé quand le pas tombe sous 10⁻¹³, ou quand il cesse de décroître
        alors qu'il est déjà petit : près d'une équivalence, le bilan du ligand est mal
        conditionné, et :math:`u` ne se fixe qu'à quelques 10⁻¹³ près, pour Newton comme
        pour la dichotomie.
        """
        pH, u = pH.copy(), u.copy()
        converge = np.zeros(V.shape, dtype=bool)
        precedent = np.full(V.shape, np.inf)
        for _ in range(NEWTON):
            gH, gL = self._bilans(V, pH, u)
            a, b, c, d = self._jacobien(V, pH, u)
            det = a * d - b * c
            dp, du = -(d * gH - b * gL) / det, -(a * gL - c * gH) / det
            # pas bornés : un Newton qui part de loin ne sort pas du domaine
            dp, du = np.clip(dp, -2.0, 2.0), np.clip(du, -2.0, 2.0)
            libres = ~converge
            pH, u = np.where(libres, pH + dp, pH), np.where(libres, u + du, u)
            pas = np.maximum(np.abs(dp), np.abs(du))
            converge |= (pas < 1e-13) | ((pas < 1e-8) & (pas >= precedent / 2))
            precedent = pas
            if np.all(converge):
                break
        return pH, u, converge & np.isfinite(pH) & np.isfinite(u)

    def etat(self, volume: npt.ArrayLike, depart: Etat | None = None) -> Etat:
        """L'équilibre en chaque volume. Avec ``depart``, Newton depuis cet état, puis la
        dichotomie là où il n'a pas convergé ; sans, la dichotomie partout."""
        V = np.atleast_1d(np.asarray(volume, dtype=np.float64))
        if np.any(V < 0):
            raise ValueError("un volume versé est positif")
        if depart is None:
            return self._dichotomie(V)
        pH, u, ok = self._newton(V, depart.pH, depart.u)
        if not np.all(ok):
            reprise = self._dichotomie(V[~ok])
            pH[~ok], u[~ok] = reprise.pH, reprise.u
        return Etat(V, pH, u)

    # --- les grandeurs -----------------------------------------------------------------

    def _total_metal(self, metal: str, V: Tableau) -> Tableau:
        return self._melange(V, self.becher.concentration(metal), self.burette.concentration(metal))

    def pM(self, e: Etat, metal: str) -> Tableau:
        with np.errstate(divide="ignore"):
            T = self._total_metal(metal, e.volume)
            return np.asarray(
                -np.log10(T) + log_somme_complexes(self.log_beta[metal], -e.u), dtype=np.float64
            )

    def pL(self, e: Etat) -> Tableau:
        """Le ligand libre toutes formes protonées comprises, :math:`-\\log(lP)`."""
        return np.asarray(-np.log10(self._ligand_libre(e.pH, e.u)[0]), dtype=np.float64)

    def derivees(self, e: Etat) -> tuple[Tableau, Tableau]:
        """:math:`\\mathrm{d\\,pH}/\\mathrm{d}V` et :math:`\\mathrm{d}u/\\mathrm{d}V`, exactes."""
        V, pH, u = e.volume, e.pH, e.u
        a, b, c, d = self._jacobien(V, pH, u)
        s = self.burette
        lP, lP1, _ = self._ligand_libre(pH, u)
        Y, _ = self._metaux(V, u, s)
        eH = (np.asarray(s.acidobasique.residu(pH)) + s.protons_ligand * s.ligand - lP1) / (
            self.volume_becher + V
        )
        eL = (Y + lP - s.ligand) / (self.volume_becher + V)
        det = a * d - b * c
        return -(d * eH - b * eL) / det, -(a * eL - c * eH) / det

    def pente(self, e: Etat, metal: str) -> Tableau:
        """:math:`\\mathrm{d\\,pM}/\\mathrm{d}V = -T'/(T\\ln 10) + Y_1\\,\\mathrm{d}u/\\mathrm{d}V`."""
        _, du = self.derivees(e)
        V = e.volume
        T = self._total_metal(metal, V)
        Cu = self.burette.concentration(metal)
        Y1, _ = moments_complexes(self.log_beta[metal], -e.u)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(-(Cu - T) / ((self.volume_becher + V) * T * LN10) + Y1 * du, dtype=np.float64)

    # --- la courbe ----------------------------------------------------------------------

    def _depart(self, metal: str, volume_max: float) -> float:
        """0, ou un milliardième de ``volume_max`` si le métal suivi n'est pas au bécher."""
        return 0.0 if self.becher.concentration(metal) > 0 else 1e-9 * volume_max

    def _raffiner(self, volume_max: float, points: int, metal: str) -> Etat:
        """Des points de plus là où la courbe tourne, jusqu'à ce que chaque segment, en
        coordonnées normalisées, soit plus court que ``1/points``."""
        e = self.etat(np.linspace(self._depart(metal, volume_max), volume_max, 200))
        for _ in range(80):
            y = self.pM(e, metal)
            etendue = max(float(np.ptp(y[np.isfinite(y)])), 1e-12)
            longueur = np.hypot(np.diff(e.volume) / volume_max, np.diff(y) / etendue)
            a_couper = (longueur > 1.0 / points) & (np.diff(e.volume) > 1e-13 * volume_max)
            if not np.any(a_couper):
                break
            i = np.flatnonzero(a_couper)
            milieu = Etat(
                (e.volume[i] + e.volume[i + 1]) / 2, (e.pH[i] + e.pH[i + 1]) / 2, (e.u[i] + e.u[i + 1]) / 2
            )
            nouveaux = self.etat(milieu.volume, milieu)
            ordre = np.argsort(np.concatenate((e.volume, nouveaux.volume)), kind="stable")
            e = Etat(
                *(
                    np.concatenate((x, y_))[ordre]
                    for x, y_ in zip(
                        (e.volume, e.pH, e.u), (nouveaux.volume, nouveaux.pH, nouveaux.u), strict=True
                    )
                )
            )
        return e

    def courbe(self, volume_max: float, points: int, metal: str) -> Etat:
        """``points`` points exacts, à longueur d'arc constante en coordonnées normalisées."""
        e = self._raffiner(volume_max, points, metal)
        y = self.pM(e, metal)
        etendue = max(float(np.ptp(y[np.isfinite(y)])), 1e-12)
        s = np.concatenate(([0.0], np.cumsum(np.hypot(np.diff(e.volume) / volume_max, np.diff(y) / etendue))))
        cibles = np.linspace(0.0, s[-1], points)
        j = np.clip(np.searchsorted(s, cibles), 1, len(s) - 1)
        t = np.where(
            s[j] > s[j - 1], (cibles - s[j - 1]) / np.where(s[j] > s[j - 1], s[j] - s[j - 1], 1.0), 0.0
        )
        V = e.volume[j - 1] + t * (e.volume[j] - e.volume[j - 1])
        depart = Etat(V, e.pH[j - 1] + t * (e.pH[j] - e.pH[j - 1]), e.u[j - 1] + t * (e.u[j] - e.u[j - 1]))
        return self.etat(V, depart)

    def sauts(self, volume_max: float, points: int, metal: str) -> list[tuple[float, Etat]]:
        """Les points d'inflexion : les minima de :math:`|\\mathrm{d}V/\\mathrm{d\\,pM}|`, en
        :math:`V`, affinés par :func:`~titrages.equilibre.minimum_lisse` avec un pas
        proportionnel à la largeur locale de la grille raffinée."""
        e = self._raffiner(volume_max, points, metal)
        garde = e.volume >= 1e-6 * volume_max
        e = Etat(e.volume[garde], e.pH[garde], e.u[garde])
        d = np.abs(1 / self.pente(e, metal))
        resultat = []
        for i in np.flatnonzero((d[1:-1] < d[:-2]) & (d[1:-1] <= d[2:])) + 1:
            proche = Etat(e.volume[i : i + 1], e.pH[i : i + 1], e.u[i : i + 1])

            def f(V: float, proche: Etat = proche) -> float:
                x = self.etat([V], Etat(np.array([V]), proche.pH, proche.u))
                return float(np.abs(1 / self.pente(x, metal))[0])

            largeur = min(e.volume[i + 1] - e.volume[i], e.volume[i] - e.volume[i - 1])
            V = minimum_lisse(f, e.volume[i - 1], e.volume[i + 1], 1e-3 * largeur)
            if 0 < V < volume_max:
                resultat.append((V, self.etat([V], proche)))
        return resultat
