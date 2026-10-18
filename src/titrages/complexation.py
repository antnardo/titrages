"""La complexométrie : un métal titré par un ligand, pM = f(V).

**Le ligand libre joue le rôle du proton.** Un métal M et ses complexes successifs
:math:`\\mathrm{ML}_j` se comportent comme un polyacide dont les « protons » seraient
des ligands : à concentration :math:`l` de ligand libre imposée, les fractions sont
explicites,

.. math:: x_j \\propto \\beta_j\\,l^j, \\qquad j = 0, \\dots, n,

et le nombre moyen de ligands fixés :math:`Y_1 = \\sum_j j\\,x_j` aussi. Le bilan du
ligand, :math:`\\sum_i T_i Y_{1,i} + l = T_L`, est alors **linéaire** dans les
quantités introduites : le résidu de chaque solution,

.. math:: R(\\mathrm{pL}) = \\sum_i C_i\\,Y_{1,i}(\\mathrm{pL}) + l - C_L,

se sépare entre le bécher et la burette exactement comme l'électroneutralité, et

.. math:: V(\\mathrm{pL}) = -V_0\\,\\frac{R_\\text{bécher}(\\mathrm{pL})}{R_\\text{burette}(\\mathrm{pL})}.

Le paramètre de la courbe est donc :math:`\\mathrm{pL}`, et non le :math:`\\mathrm{pM}`
qu'on mesure : celui-ci s'en déduit point par point,
:math:`[\\mathrm{M}] = T_\\mathrm{M}(V)\\,x_0(\\mathrm{pL})`. Rien n'est approché, et
cela vaut pour des complexes successifs comme pour plusieurs métaux en compétition
pour le même ligand.

**Le milieu est tamponné** quand on donne ``pH``. Si le ligand a des acidités (l'EDTA
en a six), seule sa forme la plus basique complexe ; à pH fixé, sa fraction
:math:`\\alpha` est une constante, et l'on travaille sur le ligand libre **toutes
formes protonées comprises**, :math:`l' = l/\\alpha`, avec les constantes
conditionnelles :math:`\\beta'_j = \\beta_j\\alpha^j`.

**Sans pH, le pH se calcule** à partir du contenu des solutions — la forme sous
laquelle on introduit le ligand, les acides et bases forts, un tampon décrit comme un
couple —, et le titrage le dit par un avertissement : c'est le moteur de
:mod:`titrages.complexation_libre`, par dichotomies imbriquées à volume fixé.
"""

import itertools
import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from functools import partial

import numpy as np
import numpy.typing as npt

from titrages.complexation_libre import ComplexationLibre
from titrages.equilibre import (
    LN10,
    fractions,
    log_somme_complexes,
    minimum_lisse,
    moments_complexes,
    racine_croissante,
    sortie,
)
from titrages.libre import Etat
from titrages.modeles import Complexe, CourbeComplexometrique, Ligand, SautComplexometrique
from titrages.solution import Solution

__all__ = ["SolutionComplexation", "TitrageComplexometrique"]

Tableau = npt.NDArray[np.float64]
Valeur = float | Tableau

# Le ligand libre ne descend pas sous 10^-PL_MAX : à ce niveau, même le complexe le plus
# stable des tables (log β ≈ 40) n'en fixe plus rien.
PL_MAX = 80.0


def _log_alpha(ligand: Ligand, pH: float | None) -> float:
    """:math:`\\log\\alpha`, la part de la forme la plus basique du ligand à ce pH."""
    if not ligand.pKa:
        return 0.0
    if pH is None:
        raise ValueError(f"le ligand {ligand.nom} a des acidités : préciser le pH du milieu")
    return float(np.log10(fractions(ligand.pKa, pH)[-1]))


@dataclass(frozen=True, slots=True)
class SolutionComplexation:
    """Une solution telle qu'on la prépare : des métaux et du ligand, en mol/L.

    Pour un titrage dont le pH se calcule, ``protons_ligand`` dit la forme sous laquelle on
    introduit le ligand — 2 pour Na₂H₂Y, l'EDTA qu'on pèse —, et ``acidobasique`` le reste
    de l'acide-base : acides et bases forts, et un tampon décrit comme un couple. À pH
    imposé, ni l'un ni l'autre ne compte : le tampon idéal absorbe tout.

    >>> SolutionComplexation({"Ca2+": 0.01}).concentration("Ca2+")
    0.01
    """

    metaux: Mapping[str, float] = field(default_factory=dict)
    ligand: float = 0.0
    protons_ligand: int = 0
    acidobasique: Solution = field(default_factory=Solution)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metaux", {m: float(c) for m, c in self.metaux.items()})
        if self.ligand < 0 or any(c < 0 for c in self.metaux.values()):
            raise ValueError("une concentration est positive")
        if self.protons_ligand < 0:
            raise ValueError("protons_ligand est un nombre de protons, positif")

    @property
    def a_de_l_acidobasique(self) -> bool:
        """Si la solution apporte des protons ou en prend : son pH ne va pas de soi."""
        a = self.acidobasique
        return bool(self.protons_ligand or a.composants or a.acide_fort or a.base_forte)

    def concentration(self, metal: str) -> float:
        return self.metaux.get(metal, 0.0)


@dataclass(frozen=True, slots=True)
class TitrageComplexometrique:
    """Le titrage complexométrique de ``volume_becher`` du ``becher`` par la ``burette``.

    ``complexes`` donne, pour chaque métal, ses complexes avec le ligand ; ``pH`` est
    celui du tampon. ``metal`` est celui dont on suit le pM — par défaut, celui du premier
    complexe.

    **Sans pH**, et dès que le pH compte — un ligand qui a des acidités, ou une solution
    qui en apporte —, le pH se calcule à partir des solutions, et un ``UserWarning`` le
    rappelle : un pH oublié change la chimie. La courbe est alors paramétrée par le
    volume ; :meth:`volume`, :meth:`derivee_volume` et :meth:`pente`, qui le sont par pL,
    n'ont plus de sens.

    >>> edta = Ligand("Y4-", (0.0, 1.5, 2.0, 2.66, 6.16, 10.24))
    >>> t = TitrageComplexometrique(
    ...     SolutionComplexation({"Ca2+": 0.01}),
    ...     SolutionComplexation(ligand=0.01),
    ...     volume_becher=20,
    ...     complexes=[Complexe("Ca2+", (10.7,))],
    ...     ligand=edta,
    ...     pH=10,
    ... )
    >>> round(t.pM(10), 2), round(t.pM(20), 2), round(t.pM(30), 2)
    (2.48, 6.28, 9.96)
    """

    becher: SolutionComplexation
    burette: SolutionComplexation
    volume_becher: float
    complexes: Sequence[Complexe]
    ligand: Ligand = field(default_factory=Ligand)
    pH: float | None = None
    metal: str | None = None
    _log_beta: dict[str, Tableau] = field(init=False, repr=False, compare=False)
    _libre: ComplexationLibre | None = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.volume_becher <= 0:
            raise ValueError("le volume du bécher est positif")
        complexes = tuple(self.complexes)
        object.__setattr__(self, "complexes", complexes)
        libre = self.pH is None and (
            bool(self.ligand.pKa) or self.becher.a_de_l_acidobasique or self.burette.a_de_l_acidobasique
        )
        if libre:
            warnings.warn(
                "pH non précisé : calculé à partir des solutions, sans tampon imposé ; "
                "préciser pH= pour un milieu tamponné",
                UserWarning,
                stacklevel=3,
            )
        la = 0.0 if libre else _log_alpha(self.ligand, self.pH)
        tables = {}
        for c in complexes:
            if c.metal in tables:
                raise ValueError(f"{c.metal} : deux jeux de complexes")
            j = np.arange(c.n + 1, dtype=np.float64)
            tables[c.metal] = np.concatenate(([0.0], np.asarray(c.log_beta))) + j * la
        object.__setattr__(self, "_log_beta", tables)
        for s in (self.becher, self.burette):
            for m in s.metaux:
                if m not in tables:
                    raise ValueError(f"{m} : aucun complexe déclaré")
        if self.metal is None:
            object.__setattr__(self, "metal", complexes[0].metal)
        elif self.metal not in tables:
            raise ValueError(f"{self.metal} : aucun complexe déclaré")
        moteur = (
            ComplexationLibre(self.becher, self.burette, self.volume_becher, tables, self.ligand.pKa)
            if libre
            else None
        )
        object.__setattr__(self, "_libre", moteur)

    def _parametre_pL(self) -> None:
        if self._libre is not None:
            raise ValueError("à pH calculé, la courbe est paramétrée par le volume, et non par pL")

    def _etat(self, volume: npt.ArrayLike) -> Etat:
        assert self._libre is not None
        return self._libre.etat(np.asarray(volume, dtype=np.float64).ravel())

    # --- le sens direct : V(pL) --------------------------------------------------------

    def _residu(self, s: SolutionComplexation, pL: npt.ArrayLike) -> Tableau:
        p = np.asarray(pL, dtype=np.float64)
        total = np.power(10.0, -p) - s.ligand
        for m, C in s.metaux.items():
            total = total + C * moments_complexes(self._log_beta[m], p)[0]
        return np.asarray(total, dtype=np.float64)

    def _derivee_residu(self, s: SolutionComplexation, pL: npt.ArrayLike) -> Tableau:
        """:math:`\\mathrm{d}R/\\mathrm{d\\,pL} = -\\ln 10\\,(\\sum_i C_i\\,\\mathrm{Var}_i + l)`."""
        p = np.asarray(pL, dtype=np.float64)
        total = np.power(10.0, -p)
        for m, C in s.metaux.items():
            total = total + C * moments_complexes(self._log_beta[m], p)[1]
        return np.asarray(-LN10 * total, dtype=np.float64)

    def volume(self, pL: npt.ArrayLike) -> Valeur:
        """Le volume versé qui amène le ligand libre à ce pL ; ``nan`` si aucun n'y amène."""
        self._parametre_pL()
        N = self._residu(self.becher, pL)
        D = self._residu(self.burette, pL)
        with np.errstate(divide="ignore", invalid="ignore"):
            V = -self.volume_becher * N / D
        return sortie(np.where(np.isfinite(V) & (V >= 0), V, np.nan))

    def derivee_volume(self, pL: npt.ArrayLike) -> Valeur:
        """:math:`\\mathrm{d}V/\\mathrm{d\\,pL} = -(V_0 N' + V D')/D`, somme de termes de même signe."""
        self._parametre_pL()
        D = self._residu(self.burette, pL)
        with np.errstate(divide="ignore", invalid="ignore"):
            V = -self.volume_becher * self._residu(self.becher, pL) / D
            dV = (
                -(
                    self.volume_becher * self._derivee_residu(self.becher, pL)
                    + V * self._derivee_residu(self.burette, pL)
                )
                / D
            )
        return sortie(dV)

    # --- le sens inverse : pL(V) -------------------------------------------------------

    def pL(self, volume: npt.ArrayLike) -> Valeur:
        """Le ligand libre après avoir versé ``volume``, par dichotomie.

        :math:`V_0 R_\\text{bécher} + V R_\\text{burette}` décroît avec pL pour tout
        :math:`V \\ge 0` : la racine est unique. Sans ligand du tout (le bécher au départ),
        il n'y en a pas, et la dichotomie s'arrête à ``PL_MAX``.
        """
        V = np.asarray(volume, dtype=np.float64)
        if np.any(V < 0):
            raise ValueError("un volume versé est positif")
        if self._libre is not None:
            return sortie(self._libre.pL(self._etat(V)).reshape(V.shape))
        V0 = self.volume_becher
        concentrations = [self.becher.ligand, self.burette.ligand]
        concentrations += [*self.becher.metaux.values(), *self.burette.metaux.values()]
        haut = math.log10(max(concentrations) + 1e-300) + 1

        def g(u: Tableau) -> Tableau:
            return V0 * self._residu(self.becher, -u) + V * self._residu(self.burette, -u)

        u = racine_croissante(g, np.full(V.shape, -PL_MAX), np.full(V.shape, haut))
        return sortie(-np.asarray(u))

    def _totale(self, metal: str, V: Tableau) -> Tableau:
        V0 = self.volume_becher
        return (V0 * self.becher.concentration(metal) + V * self.burette.concentration(metal)) / (V0 + V)

    def _pM(self, pL: Tableau, V: Tableau, metal: str) -> Tableau:
        with np.errstate(divide="ignore"):
            return -np.log10(self._totale(metal, V)) + log_somme_complexes(self._log_beta[metal], pL)

    def pM(self, volume: npt.ArrayLike, metal: str | None = None) -> Valeur:
        """Le pM du métal libre après avoir versé ``volume``."""
        V = np.asarray(volume, dtype=np.float64)
        m = metal or str(self.metal)
        if self._libre is not None:
            return sortie(self._libre.pM(self._etat(V), m).reshape(V.shape))
        return sortie(self._pM(np.asarray(self.pL(V)), V, m))

    def pH_du_melange(self, volume: npt.ArrayLike) -> Valeur:
        """Le pH après avoir versé ``volume`` : celui du tampon, ou celui qu'on calcule."""
        V = np.asarray(volume, dtype=np.float64)
        if self._libre is not None:
            return sortie(self._etat(V).pH.reshape(V.shape))
        if self.pH is None:
            raise ValueError("ni pH imposé, ni acide-base dans les solutions : le pH n'a pas d'objet")
        return sortie(np.full(V.shape, float(self.pH)))

    def _derivee_pM(self, pL: Tableau, metal: str) -> Tableau:
        """:math:`\\mathrm{d\\,pM}/\\mathrm{d\\,pL}` le long de la courbe.

        :math:`\\mathrm{pM} = -\\log T_\\mathrm{M}(V) + \\log\\sum_j\\beta'_jl^j`, donc
        :math:`\\mathrm{pM}' = -\\frac{C_u - T}{(V_0 + V)\\,T\\ln 10}\\,V' - Y_1`.
        """
        V = np.asarray(self.volume(pL))
        dV = np.asarray(self.derivee_volume(pL))
        T = self._totale(metal, V)
        Cu = self.burette.concentration(metal)
        Y1 = moments_complexes(self._log_beta[metal], pL)[0]
        with np.errstate(divide="ignore", invalid="ignore"):
            return -(Cu - T) / ((self.volume_becher + V) * T * LN10) * dV - Y1

    def pente(self, pL: npt.ArrayLike, metal: str | None = None) -> Valeur:
        """:math:`\\mathrm{d\\,pM}/\\mathrm{d}V` au point de la courbe de ce pL."""
        self._parametre_pL()
        p = np.asarray(pL, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            return sortie(self._derivee_pM(p, metal or str(self.metal)) / np.asarray(self.derivee_volume(p)))

    # --- la courbe ---------------------------------------------------------------------

    def _bornes(self, volume_max: float, depart: float = 1e-9) -> tuple[float, float]:
        """Les pL du départ, ``depart`` fois ``volume_max``, et de ``volume_max``."""
        if volume_max <= 0:
            raise ValueError("volume_max est positif")
        return float(self.pL(depart * volume_max)), float(self.pL(volume_max))

    def courbe(
        self, volume_max: float, points: int = 400, metal: str | None = None
    ) -> CourbeComplexometrique:
        """La courbe de 0 à ``volume_max``, en ``points`` points exacts, à longueur d'arc constante.

        Le premier point est pris à un milliardième de ``volume_max`` : à volume nul, le
        bécher n'a pas de ligand libre, et pL est infini. À pH calculé, la courbe part de
        zéro si le métal suivi est au bécher, et ses points sont ajoutés là où elle tourne.
        """
        m = metal or str(self.metal)
        if volume_max <= 0:
            raise ValueError("volume_max est positif")
        if self._libre is not None:
            libre = self._libre
            e = libre.courbe(volume_max, points, libre.depart(m, volume_max), partial(libre.pM, metal=m))
            return CourbeComplexometrique(e.volume, libre.pM(e, m), libre.pL(e), libre.pente(e, m), e.pH)
        pL0, pL1 = self._bornes(volume_max)
        grille = np.linspace(pL0, pL1, 40 * points)
        V = np.nan_to_num(np.asarray(self.volume(grille)), nan=0.0)
        pM = self._pM(grille, V, m)
        etendue = max(float(np.ptp(pM[np.isfinite(pM)])), 1e-12)
        pas = np.hypot(np.diff(V) / volume_max, np.diff(np.nan_to_num(pM)) / etendue)
        s = np.concatenate(([0.0], np.cumsum(pas)))
        pL = np.interp(np.linspace(0.0, s[-1], points), s, grille)
        volumes = np.nan_to_num(np.asarray(self.volume(pL)), nan=0.0)
        return CourbeComplexometrique(
            volumes,
            self._pM(pL, volumes, m),
            np.asarray(pL, dtype=np.float64),
            np.asarray(self.pente(pL, m)),
            np.full(points, np.nan if self.pH is None else float(self.pH)),
        )

    # --- les équivalences --------------------------------------------------------------

    def sauts(
        self, volume_max: float, points: int = 4000, metal: str | None = None
    ) -> tuple[SautComplexometrique, ...]:
        """Les points d'inflexion de :math:`\\mathrm{pM}(V)` entre 0 et ``volume_max``.

        Ce sont les minima de :math:`|\\mathrm{d}V/\\mathrm{d\\,pM}|`, cherchés en pL — le
        paramètre de la courbe — et affinés à la précision de la machine
        (:func:`~titrages.equilibre.minimum_lisse`). La recherche part d'un millionième de
        ``volume_max`` : avant, la courbe est plate à la précision de la machine, et ses
        minima n'y seraient que du bruit. À pH calculé, ils se cherchent en :math:`V`.
        """
        m = metal or str(self.metal)
        if self._libre is not None:
            libre = self._libre
            return tuple(
                SautComplexometrique(V, float(libre.pM(e, m)[0]), float(libre.pente(e, m)[0]))
                for V, e in libre.sauts(
                    volume_max,
                    points,
                    libre.depart(m, volume_max),
                    partial(libre.pM, metal=m),
                    partial(libre.pente, metal=m),
                )
            )
        pL0, pL1 = self._bornes(volume_max, depart=1e-6)

        def f(p: npt.ArrayLike) -> Tableau:
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.abs(
                    np.asarray(self.derivee_volume(p)) / self._derivee_pM(np.asarray(p, dtype=np.float64), m)
                )

        grille = np.linspace(pL0, pL1, points)
        d = f(grille)
        interieurs = np.flatnonzero((d[1:-1] < d[:-2]) & (d[1:-1] <= d[2:])) + 1
        resultat = []
        for i in interieurs:
            p = minimum_lisse(lambda x: float(f(x)), grille[i - 1], grille[i + 1], 1e-6)
            V = float(self.volume(p))
            if np.isfinite(V) and 0 < V < volume_max:
                resultat.append(
                    SautComplexometrique(
                        V, float(self._pM(np.asarray(p), np.asarray(V), m)), float(self.pente(p, m))
                    )
                )
        return tuple(resultat)

    def volumes_equivalents(self) -> tuple[float, ...]:
        """Les équivalences stœchiométriques, une par étape de complexation.

        Comme pour un polyacide, chaque étape compte, qu'elle donne un saut ou non. Les
        étapes sont celles des complexes qui prédominent tour à tour quand le ligand libre
        augmente (:meth:`_etapes`) : chacune fixe un nombre de ligands, avec une force, la
        pente :math:`\\Delta\\log\\beta/\\Delta j`.

        Le ligand titre les métaux du bécher : les étapes de tous les métaux sont rangées
        de la plus forte à la plus faible, et le :math:`k`-ième volume apporte le ligand
        que demandent les :math:`k` premières ; deux étapes de même force ne font qu'une
        équivalence. Un métal seul titre le ligand du bécher : le ligand, d'abord en
        excès, est fixé en :math:`\\mathrm{ML}_n` le plus riche, puis en complexes de moins
        en moins riches ; une équivalence par complexe prédominant.

        >>> t = TitrageComplexometrique(
        ...     SolutionComplexation({"Ca2+": 0.01, "Zn2+": 0.005}),
        ...     SolutionComplexation(ligand=0.01),
        ...     20,
        ...     [Complexe("Ca2+", (10.7,)), Complexe("Zn2+", (16.5,))],
        ... )
        >>> [round(v, 9) for v in t.volumes_equivalents()]
        [10.0, 30.0]
        """
        V0 = self.volume_becher
        if self.burette.ligand > 0 and not self.burette.metaux and self.becher.ligand == 0:
            etapes = [
                (force, nombre * C * V0)
                for m, C in self.becher.metaux.items()
                for force, nombre in self._etapes(m)
            ]
            etapes.sort(key=lambda e: -e[0])
            volumes: list[float] = []
            cumul, precedente = 0.0, None
            for force, quantite in etapes:
                cumul += quantite
                if force == precedente:
                    volumes[-1] = cumul / self.burette.ligand
                else:
                    volumes.append(cumul / self.burette.ligand)
                precedente = force
            return tuple(volumes)
        if (
            len(self.burette.metaux) == 1
            and self.burette.ligand == 0
            and self.becher.ligand > 0
            and not self.becher.metaux
        ):
            ((m, C),) = self.burette.metaux.items()
            sommets = self._sommets(m)[1:]
            return tuple(V0 * self.becher.ligand / (n * C) for n in reversed(sommets))
        raise ValueError(
            "équivalences définies pour un ligand qui titre des métaux, ou un métal seul qui titre le ligand"
        )

    def _sommets(self, metal: str) -> list[int]:
        """Les indices des complexes qui prédominent tour à tour, 0 (le métal libre) compris.

        Ce sont les sommets de l'enveloppe concave des points :math:`(j, \\log\\beta_j)` :
        un complexe dont la constante successive dépasse la précédente (Zn–NH₃) ou que les
        tables ne donnent pas (:math:`\\log\\beta = -\\infty`) ne prédomine jamais, et son
        étape se fond dans la suivante. À pente égale, le sommet le plus lointain : deux
        étapes de même force n'en font qu'une, comme deux sites de même pKa.
        """
        lb = self._log_beta[metal]
        sommets, j = [0], 0
        while j < len(lb) - 1:
            pentes = [((lb[k] - lb[j]) / (k - j), k) for k in range(j + 1, len(lb)) if np.isfinite(lb[k])]
            if not pentes:
                break
            j = max(pentes)[1]
            sommets.append(j)
        return sommets

    def _etapes(self, metal: str) -> list[tuple[float, int]]:
        """Les étapes de complexation du métal : (force, ligands fixés), de la plus forte."""
        lb, s = self._log_beta[metal], self._sommets(metal)
        return [(float(lb[b] - lb[a]) / (b - a), b - a) for a, b in itertools.pairwise(s)]
