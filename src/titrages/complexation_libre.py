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
    moments_complexes,
    racine_croissante,
)
from titrages.libre import Etat, MoteurLibre

if TYPE_CHECKING:
    from titrages.complexation import SolutionComplexation

__all__ = ["ComplexationLibre"]

Tableau = npt.NDArray[np.float64]

PL_MAX = 80.0


@dataclass(frozen=True)
class ComplexationLibre(MoteurLibre):
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
        return Etat.sans_solide(V, pH, self._u_de(V, pH))

    def _pas_newton(self, e: Etat) -> tuple[Tableau, Tableau]:
        gH, gL = self._bilans(e.volume, e.pH, e.u)
        a, b, c, d = self._jacobien(e.volume, e.pH, e.u)
        det = a * d - b * c
        return -(d * gH - b * gL) / det, -(a * gL - c * gH) / det

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

    def depart(self, metal: str, volume_max: float) -> float:
        """0, ou un milliardième de ``volume_max`` si le métal suivi n'est pas au bécher."""
        return 0.0 if self.becher.concentration(metal) > 0 else 1e-9 * volume_max
