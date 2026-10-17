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

**Le milieu est tamponné.** Si le ligand a des acidités (l'EDTA en a six), seule sa
forme la plus basique complexe ; à pH fixé, sa fraction :math:`\\alpha` est une
constante, et l'on travaille sur le ligand libre **toutes formes protonées
comprises**, :math:`l' = l/\\alpha`, avec les constantes conditionnelles
:math:`\\beta'_j = \\beta_j\\alpha^j`. Un pH qui bouge au fil du titrage — un ligand
basique en milieu non tamponné — sort de ce cadre : il faudrait un solveur.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from titrages.equilibre import LN10, fractions, minimum_lisse, racine_croissante, sortie
from titrages.modeles import Complexe, CourbeComplexometrique, Ligand, SautComplexometrique

__all__ = ["SolutionComplexation", "TitrageComplexometrique", "fractions_complexes"]

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


def _log_poids(log_beta: Tableau, pL: npt.ArrayLike) -> Tableau:
    """Les :math:`\\log(\\beta'_j l^j)`, de forme ``(n + 1, *pL.shape)``, :math:`\\beta'_0 = 1`."""
    p = np.asarray(pL, dtype=np.float64)
    forme = (-1,) + (1,) * p.ndim
    j = np.arange(len(log_beta), dtype=np.float64).reshape(forme)
    return np.asarray(log_beta.reshape(forme) - j * p, dtype=np.float64)


def fractions_complexes(log_beta: npt.ArrayLike, pL: npt.ArrayLike) -> Tableau:
    """Les fractions du métal libre (0), de ML (1), … , de forme ``(n + 1, *pL.shape)``.

    ``log_beta`` sont les constantes globales :math:`\\log\\beta_1, \\dots,
    \\log\\beta_n` ; le calcul est en échelle logarithmique, comme celui des polyacides.

    >>> x = fractions_complexes([10.7], [9.7, 10.7, 11.7])
    >>> [round(float(v), 3) for v in x[1]]
    [0.909, 0.5, 0.091]
    """
    lb = np.concatenate(([0.0], np.asarray(log_beta, dtype=np.float64)))
    t = _log_poids(lb, pL)
    t = t - t.max(axis=0)
    e = np.power(10.0, t)
    return np.asarray(e / e.sum(axis=0), dtype=np.float64)


def _log_somme(log_beta: Tableau, pL: npt.ArrayLike) -> Tableau:
    """:math:`\\log\\sum_j \\beta'_j l^j = -\\log x_0`, sans débordement."""
    t = _log_poids(log_beta, pL)
    haut = t.max(axis=0)
    return np.asarray(haut + np.log10(np.power(10.0, t - haut).sum(axis=0)), dtype=np.float64)


def _moments(log_beta: Tableau, pL: npt.ArrayLike) -> tuple[Tableau, Tableau]:
    """:math:`Y_1` et la variance :math:`\\sum_{j<k} x_jx_k(k-j)^2`, écrite sans soustraction."""
    lb = log_beta[1:]
    x = fractions_complexes(lb, pL)
    j = np.arange(x.shape[0], dtype=np.float64).reshape((-1,) + (1,) * (x.ndim - 1))
    Y1 = (j * x).sum(axis=0)
    var = np.zeros(x.shape[1:], dtype=np.float64)
    for a in range(x.shape[0]):
        for b in range(a + 1, x.shape[0]):
            var = var + x[a] * x[b] * (b - a) ** 2
    return np.asarray(Y1, dtype=np.float64), var


@dataclass(frozen=True, slots=True)
class SolutionComplexation:
    """Une solution telle qu'on la prépare : des métaux et du ligand, en mol/L.

    >>> SolutionComplexation({"Ca2+": 0.01}).concentration("Ca2+")
    0.01
    """

    metaux: Mapping[str, float] = field(default_factory=dict)
    ligand: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "metaux", {m: float(c) for m, c in self.metaux.items()})
        if self.ligand < 0 or any(c < 0 for c in self.metaux.values()):
            raise ValueError("une concentration est positive")

    def concentration(self, metal: str) -> float:
        return self.metaux.get(metal, 0.0)


@dataclass(frozen=True, slots=True)
class TitrageComplexometrique:
    """Le titrage complexométrique de ``volume_becher`` du ``becher`` par la ``burette``.

    ``complexes`` donne, pour chaque métal, ses complexes avec le ligand ; ``pH`` est
    celui du tampon, nécessaire dès que le ligand a des acidités. ``metal`` est celui dont
    on suit le pM — par défaut, celui du premier complexe.

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

    def __post_init__(self) -> None:
        if self.volume_becher <= 0:
            raise ValueError("le volume du bécher est positif")
        complexes = tuple(self.complexes)
        object.__setattr__(self, "complexes", complexes)
        la = _log_alpha(self.ligand, self.pH)
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

    # --- le sens direct : V(pL) --------------------------------------------------------

    def _residu(self, s: SolutionComplexation, pL: npt.ArrayLike) -> Tableau:
        p = np.asarray(pL, dtype=np.float64)
        total = np.power(10.0, -p) - s.ligand
        for m, C in s.metaux.items():
            total = total + C * _moments(self._log_beta[m], p)[0]
        return np.asarray(total, dtype=np.float64)

    def _derivee_residu(self, s: SolutionComplexation, pL: npt.ArrayLike) -> Tableau:
        """:math:`\\mathrm{d}R/\\mathrm{d\\,pL} = -\\ln 10\\,(\\sum_i C_i\\,\\mathrm{Var}_i + l)`."""
        p = np.asarray(pL, dtype=np.float64)
        total = np.power(10.0, -p)
        for m, C in s.metaux.items():
            total = total + C * _moments(self._log_beta[m], p)[1]
        return np.asarray(-LN10 * total, dtype=np.float64)

    def volume(self, pL: npt.ArrayLike) -> Valeur:
        """Le volume versé qui amène le ligand libre à ce pL ; ``nan`` si aucun n'y amène."""
        N = self._residu(self.becher, pL)
        D = self._residu(self.burette, pL)
        with np.errstate(divide="ignore", invalid="ignore"):
            V = -self.volume_becher * N / D
        return sortie(np.where(np.isfinite(V) & (V >= 0), V, np.nan))

    def derivee_volume(self, pL: npt.ArrayLike) -> Valeur:
        """:math:`\\mathrm{d}V/\\mathrm{d\\,pL} = -(V_0 N' + V D')/D`, somme de termes de même signe."""
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
            return -np.log10(self._totale(metal, V)) + _log_somme(self._log_beta[metal], pL)

    def pM(self, volume: npt.ArrayLike, metal: str | None = None) -> Valeur:
        """Le pM du métal libre après avoir versé ``volume``."""
        V = np.asarray(volume, dtype=np.float64)
        return sortie(self._pM(np.asarray(self.pL(V)), V, metal or str(self.metal)))

    def _derivee_pM(self, pL: Tableau, metal: str) -> Tableau:
        """:math:`\\mathrm{d\\,pM}/\\mathrm{d\\,pL}` le long de la courbe.

        :math:`\\mathrm{pM} = -\\log T_\\mathrm{M}(V) + \\log\\sum_j\\beta'_jl^j`, donc
        :math:`\\mathrm{pM}' = -\\frac{C_u - T}{(V_0 + V)\\,T\\ln 10}\\,V' - Y_1`.
        """
        V = np.asarray(self.volume(pL))
        dV = np.asarray(self.derivee_volume(pL))
        T = self._totale(metal, V)
        Cu = self.burette.concentration(metal)
        Y1 = _moments(self._log_beta[metal], pL)[0]
        with np.errstate(divide="ignore", invalid="ignore"):
            return -(Cu - T) / ((self.volume_becher + V) * T * LN10) * dV - Y1

    def pente(self, pL: npt.ArrayLike, metal: str | None = None) -> Valeur:
        """:math:`\\mathrm{d\\,pM}/\\mathrm{d}V` au point de la courbe de ce pL."""
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
        bécher n'a pas de ligand libre, et pL est infini.
        """
        m = metal or str(self.metal)
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
            volumes, self._pM(pL, volumes, m), np.asarray(pL, dtype=np.float64), np.asarray(self.pente(pL, m))
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
        minima n'y seraient que du bruit.
        """
        m = metal or str(self.metal)
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

    def volumes_equivalents(self, stoechiometrie: Mapping[str, int] | None = None) -> tuple[float, ...]:
        """Les équivalences stœchiométriques.

        Deux cas ont une définition simple. Le ligand titre les métaux du bécher : chacun
        fixe :math:`n` ligands, du complexe le plus stable (par ligand, :math:`\\log\\beta_n/n`)
        au moins stable, et le :math:`k`-ième volume apporte autant de ligand que les
        :math:`k` premiers métaux en demandent. Un métal seul titre le ligand du bécher :
        l'équivalence est celle où tout le ligand est fixé en :math:`\\mathrm{ML}_n`.

        :math:`n` est par défaut le dernier complexe déclaré ; ``stoechiometrie`` le fixe
        par métal quand ce n'est pas lui qui marque l'équivalence — :math:`\\mathrm{Ag(CN)_2^-}`
        dans le dosage de Liebig, où :math:`\\mathrm{Ag(CN)_3^{2-}}` est déclaré mais peu stable.

        >>> t = TitrageComplexometrique(
        ...     SolutionComplexation({"Ca2+": 0.01, "Zn2+": 0.005}),
        ...     SolutionComplexation(ligand=0.01),
        ...     20,
        ...     [Complexe("Ca2+", (10.7,)), Complexe("Zn2+", (16.5,))],
        ... )
        >>> [round(v, 9) for v in t.volumes_equivalents()]
        [10.0, 30.0]
        """
        choix = dict(stoechiometrie or {})
        for m, n in choix.items():
            if m not in self._log_beta or not 1 <= n <= self._n(m) or not np.isfinite(self._log_beta[m][n]):
                raise ValueError(f"{m} : pas de complexe ML{n} déclaré")

        def n_de(m: str) -> int:
            return choix.get(m, self._n(m))

        V0 = self.volume_becher
        if self.burette.ligand > 0 and not self.burette.metaux and self.becher.ligand == 0:
            ordre = sorted(self.becher.metaux, key=lambda m: -float(self._log_beta[m][n_de(m)]) / n_de(m))
            volumes, cumul = [], 0.0
            for m in ordre:
                cumul += n_de(m) * self.becher.metaux[m] * V0
                volumes.append(cumul / self.burette.ligand)
            return tuple(volumes)
        if (
            len(self.burette.metaux) == 1
            and self.burette.ligand == 0
            and self.becher.ligand > 0
            and not self.becher.metaux
        ):
            ((m, C),) = self.burette.metaux.items()
            return (V0 * self.becher.ligand / (n_de(m) * C),)
        raise ValueError(
            "équivalences définies pour un ligand qui titre des métaux, ou un métal seul qui titre le ligand"
        )

    def _n(self, metal: str) -> int:
        return len(self._log_beta[metal]) - 1
