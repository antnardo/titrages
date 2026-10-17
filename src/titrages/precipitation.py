"""Les titrages par précipitation : l'ion commun en fonction du volume, pAg = f(V).

Un ion commun P — :math:`\\mathrm{Ag^+}` en argentimétrie — forme avec des partenaires
:math:`\\mathrm{X}_i` des solides :math:`\\mathrm{P}_{\\nu}\\mathrm{X}_{\\mu}`, de produit
de solubilité :math:`K_s = a^\\nu x^\\mu`, où :math:`a = [\\mathrm{P}]`.

**À a imposé, tout est explicite.** Un solide présent fixe la concentration de son
partenaire, :math:`x_\\text{sat} = (K_s/a^\\nu)^{1/\\mu}` ; un solide absent laisse le
partenaire en solution. En notant :math:`r = \\nu/\\mu` les ions communs pris par
partenaire précipité, le bilan de l'ion commun du mélange est

.. math:: a + \\sum_{i \\in S} r_i\\,\\bigl(T_{\\mathrm{X}_i} - x_{\\text{sat},i}(a)\\bigr) = T_\\mathrm{P},

où :math:`S` est l'ensemble des solides présents. À :math:`S` fixé, c'est encore une
relation **linéaire** dans les quantités introduites : le résidu de chaque solution,

.. math::

    R_S(\\mathrm{pP}) = a + \\sum_{i \\in S} r_i\\,(C_{\\mathrm{X}_i} - x_{\\text{sat},i}) - C_\\mathrm{P},

se sépare entre le bécher et la burette, et
:math:`V = -V_0\\,R_{S,\\text{bécher}}/R_{S,\\text{burette}}`. L'ensemble :math:`S` est
celui qui est cohérent avec le volume trouvé : un solide présent a une quantité
positive, :math:`T_{\\mathrm{X}_i}(V) \\ge x_{\\text{sat},i}`, un solide absent n'est pas
sursaturé. On essaie chaque ensemble — il y en a :math:`2^n` pour :math:`n` partenaires,
quatre pour un mélange de chlorure et d'iodure — et un seul convient.

**La courbe a des points anguleux** là où un solide apparaît ou disparaît : la pente y
saute. Ce sont des :class:`~titrages.modeles.Evenement`, distincts des points
d'inflexion, cherchés à l'intérieur de chaque morceau lisse. L'apparition du chromate
d'argent est le repère de la méthode de Mohr.

**Le milieu est tamponné.** Un partenaire basique (l'ion chromate) ne précipite que sous
sa forme la plus basique ; à pH fixé, sa fraction :math:`\\alpha` est une constante, et
:math:`x_\\text{sat}` porte sur le partenaire toutes formes protonées comprises :
:math:`x_\\text{sat} = (K_s/a^\\nu)^{1/\\mu}/\\alpha`.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from functools import partial

import numpy as np
import numpy.typing as npt

from titrages.equilibre import LN10, changement_de_signe, fractions, minimum_lisse, racine_croissante, sortie
from titrages.modeles import CourbePrecipitation, Evenement, Precipite, SautPrecipitation

__all__ = ["SolutionPrecipitation", "TitrageParPrecipitation"]

Tableau = npt.NDArray[np.float64]
Valeur = float | Tableau

# L'ion commun ne descend pas sous 10^-P_MAX : là, aucun solide des tables ne se forme.
P_MAX = 60.0


@dataclass(frozen=True, slots=True)
class SolutionPrecipitation:
    """Une solution telle qu'on la prépare : l'ion commun et des partenaires, en mol/L.

    >>> SolutionPrecipitation(partenaires={"Cl-": 0.05}).concentration("Cl-")
    0.05
    """

    commun: float = 0.0
    partenaires: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "partenaires", {x: float(c) for x, c in self.partenaires.items()})
        if self.commun < 0 or any(c < 0 for c in self.partenaires.values()):
            raise ValueError("une concentration est positive")

    def concentration(self, partenaire: str) -> float:
        return self.partenaires.get(partenaire, 0.0)


@dataclass(frozen=True, slots=True)
class TitrageParPrecipitation:
    """Le titrage de ``volume_becher`` du ``becher`` par la ``burette``, suivi par pP.

    ``precipites`` donne les solides possibles ; ``pH`` est celui du tampon, nécessaire
    dès qu'un partenaire a des acidités. ``ion_commun`` ne sert qu'à nommer.

    >>> t = TitrageParPrecipitation(
    ...     SolutionPrecipitation(partenaires={"Cl-": 0.1}),
    ...     SolutionPrecipitation(commun=0.1),
    ...     volume_becher=20,
    ...     precipites=[Precipite("Cl-", 9.75)],
    ... )
    >>> round(t.p(10), 3), round(t.p(20), 3), round(t.p(30), 3)
    (8.273, 4.875, 1.699)
    """

    becher: SolutionPrecipitation
    burette: SolutionPrecipitation
    volume_becher: float
    precipites: Sequence[Precipite]
    pH: float | None = None
    ion_commun: str = "Ag+"
    _log_alpha: Tableau = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.volume_becher <= 0:
            raise ValueError("le volume du bécher est positif")
        precipites = tuple(self.precipites)
        object.__setattr__(self, "precipites", precipites)
        noms = [p.partenaire for p in precipites]
        if len(set(noms)) != len(noms):
            raise ValueError("un partenaire ne forme qu'un précipité")
        for s in (self.becher, self.burette):
            for x in s.partenaires:
                if x not in noms:
                    raise ValueError(f"{x} : aucun précipité déclaré")
        alpha = []
        for p in precipites:
            if not p.pKa:
                alpha.append(0.0)
            elif self.pH is None:
                raise ValueError(f"{p.partenaire} a des acidités : préciser le pH du milieu")
            else:
                alpha.append(float(np.log10(fractions(p.pKa, self.pH)[-1])))
        object.__setattr__(self, "_log_alpha", np.asarray(alpha, dtype=np.float64))

    # --- les briques ---------------------------------------------------------------------

    def _x_sat(self, i: int, pP: Tableau) -> Tableau:
        """La concentration du partenaire ``i`` au contact de son solide, toutes formes comprises."""
        p = self.precipites[i]
        return np.power(10.0, (-p.pKs + p.nu * pP) / p.mu - self._log_alpha[i])

    def _r(self, i: int) -> float:
        p = self.precipites[i]
        return p.nu / p.mu

    def _residu(self, s: SolutionPrecipitation, pP: Tableau, S: int) -> Tableau:
        total = np.power(10.0, -pP) - s.commun
        for i, p in enumerate(self.precipites):
            if S >> i & 1:
                total = total + self._r(i) * (s.concentration(p.partenaire) - self._x_sat(i, pP))
        return np.asarray(total, dtype=np.float64)

    def _derivee_residu(self, pP: Tableau, S: int) -> Tableau:
        """:math:`\\mathrm{d}R_S/\\mathrm{d\\,pP} = -\\ln 10\\,(a + \\sum_{i\\in S} r_i^2 x_{\\text{sat},i})`.

        Elle ne dépend que de :math:`a` et des solides présents : la même pour toute solution.
        """
        total = np.power(10.0, -pP)
        for i in range(len(self.precipites)):
            if S >> i & 1:
                total = total + self._r(i) ** 2 * self._x_sat(i, pP)
        return np.asarray(-LN10 * total, dtype=np.float64)

    def _totale(self, partenaire: str, V: Tableau) -> Tableau:
        V0 = self.volume_becher
        return (V0 * self.becher.concentration(partenaire) + V * self.burette.concentration(partenaire)) / (
            V0 + V
        )

    def _volume_S(self, pP: Tableau, S: int) -> Tableau:
        N = self._residu(self.becher, pP, S)
        D = self._residu(self.burette, pP, S)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(-self.volume_becher * N / D, dtype=np.float64)

    def _coherent(self, pP: Tableau, V: Tableau, S: int) -> npt.NDArray[np.bool_]:
        """Où l'ensemble ``S`` est le bon : ses solides existent, les autres ne sont pas sursaturés."""
        ok = np.isfinite(V) & (V >= 0)
        for i, p in enumerate(self.precipites):
            with np.errstate(invalid="ignore"):
                T = self._totale(p.partenaire, V)
                x = self._x_sat(i, pP)
                marge = 1e-12 * np.maximum(T, x)
                ok &= (T - x >= -marge) if S >> i & 1 else (T - x <= marge)
        return ok

    def _resoudre(self, pP: npt.ArrayLike) -> tuple[Tableau, npt.NDArray[np.int64]]:
        """Le volume et l'ensemble des solides présents, en chaque pP ; ``nan`` et -1 hors d'atteinte."""
        p = np.asarray(pP, dtype=np.float64)
        V = np.full(p.shape, np.nan)
        S_ = np.full(p.shape, -1, dtype=np.int64)
        for S in range(2 ** len(self.precipites)):
            libre = S_ < 0
            if not libre.any():
                break
            VS = self._volume_S(p, S)
            # Au départ, V vaut zéro à l'arrondi près, et peut sortir à -1e-16 : c'est zéro.
            with np.errstate(invalid="ignore"):
                VS = np.where((VS < 0) & (-1e-12 * self.volume_becher < VS), 0.0, VS)
            pris = libre & self._coherent(p, VS, S)
            V = np.where(pris, VS, V)
            S_ = np.where(pris, S, S_)
        return V, S_

    # --- le sens direct : V(pP) ----------------------------------------------------------

    def volume(self, pP: npt.ArrayLike) -> Valeur:
        """Le volume versé qui amène l'ion commun à ce pP ; ``nan`` si aucun n'y amène."""
        return sortie(self._resoudre(pP)[0])

    def solides(self, pP: float) -> tuple[str, ...]:
        """Les partenaires dont le solide est présent au point de la courbe de ce pP."""
        S = int(self._resoudre(np.asarray(pP))[1])
        return tuple(p.partenaire for i, p in enumerate(self.precipites) if S >= 0 and S >> i & 1)

    def _derivee_volume_S(self, pP: Tableau, S: int) -> Tableau:
        D = self._residu(self.burette, pP, S)
        V = self._volume_S(pP, S)
        dR = self._derivee_residu(pP, S)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(-(self.volume_becher + V) * dR / D, dtype=np.float64)

    def derivee_volume(self, pP: npt.ArrayLike) -> Valeur:
        """:math:`\\mathrm{d}V/\\mathrm{d\\,pP}`, sur le morceau de courbe où tombe ce pP.

        Les deux résidus ont la même dérivée — elle ne dépend que de :math:`a` et des
        solides présents —, d'où :math:`V' = -(V_0 + V)\\,R'_S/R_{S,\\text{burette}}`.
        """
        p = np.asarray(pP, dtype=np.float64)
        _, S_ = self._resoudre(p)
        dV = np.full(p.shape, np.nan)
        for S in np.unique(S_[S_ >= 0]):
            dV = np.where(S_ == S, self._derivee_volume_S(p, int(S)), dV)
        return sortie(dV)

    def pente(self, pP: npt.ArrayLike) -> Valeur:
        """:math:`\\mathrm{d\\,pP}/\\mathrm{d}V` au point de la courbe de ce pP."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return sortie(1 / np.asarray(self.derivee_volume(pP)))

    # --- le sens inverse : pP(V) ---------------------------------------------------------

    def p(self, volume: npt.ArrayLike) -> Valeur:
        """Le pP de l'ion commun après avoir versé ``volume``, par dichotomie.

        :math:`a + \\sum_i r_i \\max(0, T_{\\mathrm{X}_i} - x_{\\text{sat},i}(a)) - T_\\mathrm{P}`
        croît avec :math:`a` : la racine est unique, et les solides présents en découlent.
        """
        V = np.asarray(volume, dtype=np.float64)
        if np.any(V < 0):
            raise ValueError("un volume versé est positif")
        V0 = self.volume_becher
        T_P = (V0 * self.becher.commun + V * self.burette.commun) / (V0 + V)
        T_X = [self._totale(p.partenaire, V) for p in self.precipites]
        concentrations = [self.becher.commun, self.burette.commun]
        concentrations += [*self.becher.partenaires.values(), *self.burette.partenaires.values()]
        haut = math.log10(max(concentrations) + 1e-300) + 1

        def g(u: Tableau) -> Tableau:
            total = np.power(10.0, u) - T_P
            for i in range(len(self.precipites)):
                total = total + self._r(i) * np.maximum(0.0, T_X[i] - self._x_sat(i, -u))
            return total

        u = racine_croissante(g, np.full(V.shape, -P_MAX), np.full(V.shape, haut))
        return sortie(-np.asarray(u))

    def p_partenaire(self, volume: npt.ArrayLike, partenaire: str) -> Valeur:
        """Le cologarithme du partenaire en solution, toutes formes protonées comprises."""
        V = np.asarray(volume, dtype=np.float64)
        i = [p.partenaire for p in self.precipites].index(partenaire)
        x = np.minimum(self._totale(partenaire, V), self._x_sat(i, np.asarray(self.p(V))))
        with np.errstate(divide="ignore"):
            return sortie(-np.log10(x))

    # --- la courbe -----------------------------------------------------------------------

    def _bornes(self, volume_max: float, depart: float) -> tuple[float, float]:
        """Le pP de départ et celui de ``volume_max``. Sans ion commun dans le bécher, pP est
        infini à volume nul : on part alors de ``depart`` fois ``volume_max``."""
        if volume_max <= 0:
            raise ValueError("volume_max est positif")
        p0 = float(self.p(0.0))
        if p0 > P_MAX - 1:
            p0 = float(self.p(depart * volume_max))
        return p0, float(self.p(volume_max))

    def courbe(self, volume_max: float, points: int = 400) -> CourbePrecipitation:
        """La courbe de 0 à ``volume_max``, en ``points`` points exacts, à longueur d'arc constante.

        Sans ion commun dans le bécher, elle part d'un millionième de ``volume_max`` :
        à volume nul, pP est infini.
        """
        p0, p1 = self._bornes(volume_max, 1e-6)
        grille = np.linspace(p0, p1, 40 * points)
        V = np.nan_to_num(self._resoudre(grille)[0], nan=0.0)
        pas = np.hypot(np.diff(V) / volume_max, np.diff(grille) / max(abs(p1 - p0), 1e-12))
        s = np.concatenate(([0.0], np.cumsum(pas)))
        pP = np.interp(np.linspace(0.0, s[-1], points), s, grille)
        volumes, S_ = self._resoudre(pP)
        noms = tuple(
            tuple(p.partenaire for i, p in enumerate(self.precipites) if S >= 0 and S >> i & 1) for S in S_
        )
        return CourbePrecipitation(
            np.nan_to_num(volumes, nan=0.0),
            np.asarray(pP, dtype=np.float64),
            np.asarray(self.pente(pP)),
            noms,
        )

    def evenements(self, volume_max: float, points: int = 4000) -> tuple[Evenement, ...]:
        """Les apparitions et disparitions de solides entre 0 et ``volume_max``, dans l'ordre.

        Chacune est le point où la quantité de solide s'annule,
        :math:`T_{\\mathrm{X}}(V) = x_\\text{sat}(a)`, trouvé par dichotomie sur la courbe.
        """
        p0, p1 = self._bornes(volume_max, 1e-9)
        grille = np.linspace(p0, p1, points)
        _, S_ = self._resoudre(grille)
        resultat = []
        for k in np.flatnonzero(S_[1:] != S_[:-1]):
            avant, apres = int(S_[k]), int(S_[k + 1])
            if avant < 0 or apres < 0:
                continue
            for i, p in enumerate(self.precipites):
                if (avant ^ apres) >> i & 1:
                    pP = changement_de_signe(
                        partial(self._sursaturation, i), float(grille[k]), float(grille[k + 1])
                    )
                    resultat.append(Evenement(float(self.volume(pP)), pP, p.partenaire, bool(apres >> i & 1)))
        return tuple(resultat)

    def _sursaturation(self, i: int, pP: float) -> float:
        V = self._resoudre(np.asarray(pP))[0]
        return float(self._totale(self.precipites[i].partenaire, V) - self._x_sat(i, np.asarray(pP)))

    # --- les équivalences ----------------------------------------------------------------

    def sauts(self, volume_max: float, points: int = 4000) -> tuple[SautPrecipitation, ...]:
        """Les points d'inflexion de pP(V) entre 0 et ``volume_max``, hors des points anguleux.

        Ce sont les minima de :math:`|\\mathrm{d}V/\\mathrm{d\\,pP}|` à l'intérieur d'un
        morceau lisse — les mêmes solides de part et d'autre —, affinés à la précision de la
        machine. Un point anguleux n'est pas un saut : :meth:`evenements` les donne.
        """
        p0, p1 = self._bornes(volume_max, 1e-6)
        grille = np.linspace(p0, p1, points)
        _, S_ = self._resoudre(grille)
        d = np.abs(np.asarray(self.derivee_volume(grille)))
        resultat = []
        for i in range(1, points - 1):
            S = int(S_[i])
            if S < 0 or S_[i - 1] != S or S_[i + 1] != S or not (d[i] < d[i - 1] and d[i] <= d[i + 1]):
                continue
            q = minimum_lisse(partial(self._pente_inverse, S), grille[i - 1], grille[i + 1], 1e-6)
            V = float(self._volume_S(np.asarray(q), S))
            if np.isfinite(V) and 0 < V < volume_max:
                resultat.append(SautPrecipitation(V, q, float(1 / self._derivee_volume_S(np.asarray(q), S))))
        return tuple(resultat)

    def _pente_inverse(self, S: int, pP: float) -> float:
        """:math:`|\\mathrm{d}V/\\mathrm{d\\,pP}|` sur le morceau de courbe des solides ``S``."""
        return float(np.abs(self._derivee_volume_S(np.asarray(pP), S)))

    def volumes_equivalents(self) -> tuple[float, ...]:
        """Les équivalences stœchiométriques.

        L'ion commun titre les partenaires du bécher : ils précipitent l'un après l'autre,
        du moins soluble au plus soluble — par le pP où chacun commencerait à précipiter
        dans le bécher —, et le :math:`k`-ième volume apporte de quoi précipiter les
        :math:`k` premiers. Un partenaire titre l'ion commun du bécher : l'équivalence est
        celle où l'ion commun resté en solution, une fois précipités les partenaires du
        bécher (le titrage en retour de Charpentier-Volhard), est consommé.

        >>> t = TitrageParPrecipitation(
        ...     SolutionPrecipitation(partenaires={"Cl-": 0.05, "I-": 0.05}),
        ...     SolutionPrecipitation(commun=0.1),
        ...     20,
        ...     [Precipite("Cl-", 9.75), Precipite("I-", 16.08)],
        ... )
        >>> t.volumes_equivalents()
        (10.0, 20.0)
        """
        V0 = self.volume_becher
        if self.burette.commun > 0 and not self.burette.partenaires and self.becher.commun == 0:
            presents = [
                i for i, p in enumerate(self.precipites) if self.becher.concentration(p.partenaire) > 0
            ]
            ordre = sorted(presents, key=lambda i: -self._p_apparition(i))
            volumes, cumul = [], 0.0
            for i in ordre:
                cumul += self._r(i) * self.becher.concentration(self.precipites[i].partenaire) * V0
                volumes.append(cumul / self.burette.commun)
            return tuple(volumes)
        if len(self.burette.partenaires) == 1 and self.burette.commun == 0 and self.becher.commun > 0:
            ((x, C),) = self.burette.partenaires.items()
            j = [p.partenaire for p in self.precipites].index(x)
            reste = self.becher.commun - sum(
                self._r(i) * self.becher.concentration(p.partenaire) for i, p in enumerate(self.precipites)
            )
            if reste > 0:
                return (V0 * reste / (self._r(j) * C),)
        raise ValueError(
            "équivalences définies pour l'ion commun qui titre des partenaires, "
            "ou un partenaire seul qui titre l'ion commun restant"
        )

    def _p_apparition(self, i: int) -> float:
        """Le pP où le solide ``i`` apparaîtrait dans le bécher : :math:`x_\\text{sat} = C`."""
        p = self.precipites[i]
        C = self.becher.concentration(p.partenaire)
        return float((p.mu * (math.log10(C) + self._log_alpha[i]) + p.pKs) / p.nu)
