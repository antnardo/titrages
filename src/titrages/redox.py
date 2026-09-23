"""L'oxydoréduction : potentiel d'une solution, titrage potentiométrique E = f(V).

C'est la transposition exacte du calcul acido-basique, les électrons remplaçant les
protons et :math:`E/s` le pH, avec :math:`s = RT\\ln 10/F` (0,05916 V à 25 °C).

**À potentiel imposé, tout est explicite.** Pour un système à :math:`N+1` formes
(de la plus réduite, 0, à la plus oxydée), la loi de Nernst donne les rapports de
concentrations de proche en proche ; les fractions s'écrivent donc comme celles
d'un polyacide, en échelle logarithmique :

.. math::

    x_j \\propto 10^{\\left(e_j E - \\sum_{i\\le j} n_i E^{\\circ\\prime}_i\\right)/s},
    \\qquad E^{\\circ\\prime}_i = E^\\circ_i - s\\,\\frac{m_i}{n_i}\\,\\mathrm{pH},

où :math:`e_j` est le nombre d'électrons cédés depuis la forme 0. Les protons ne sont
pas une inconnue : le milieu est supposé tamponné, et le pH ne fait que décaler les
potentiels standard (potentiel standard apparent).

**Le bilan remplace l'électroneutralité.** Les électrons se conservent : ceux que
portent les formes présentes sont ceux qu'on a introduits. Avec :math:`Y` le nombre
moyen d'électrons cédés par unité d'élément et :math:`e_q` celui de la forme
introduite, le résidu

.. math:: R(E) = \\sum_i T_i\\,\\bigl(Y_i(E) - e_{q_i}\\bigr)

est nul au potentiel de la solution et strictement croissant.

**Les couples dimériques** (:math:`\\mathrm{I_2/I^-}`, :math:`\\mathrm{Cr_2O_7^{2-}/Cr^{3+}}`,
:math:`\\mathrm{S_4O_6^{2-}/S_2O_3^{2-}}`, :math:`\\mathrm{Hg^{2+}/Hg_2^{2+}}`) n'ont pas de
:math:`Y(E)` : leur loi de Nernst fait intervenir un carré, si bien qu'une simple
dilution, à composition fixée, change le potentiel. Leur fraction oxydée dépend de
:math:`E` **et** de la concentration — elle reste explicite, racine d'un trinôme
(:func:`electrons_cedes`), mais le résidu du mélange ne se sépare plus en une part du
bécher et une part de la burette.

Le sens direct survit pourtant à **un** dimère : en prenant pour inconnue
:math:`W = (V_0+V)\\,u`, où :math:`u` est la concentration de sa forme monomère, le
bilan d'électrons redevient linéaire, et la conservation de l'élément est une
équation **du second degré en V** (:meth:`TitrageRedox.volume`). À deux dimères
(l'iodométrie, :math:`\\mathrm{I_2}` par :math:`\\mathrm{S_2O_3^{2-}}`), le degré monte :
on résout alors par dichotomie. La courbe elle-même est tracée à partir de
:math:`E(V)`, par dichotomie sur le résidu du mélange, ce qui vaut pour tous les cas.
"""

from collections.abc import Callable
from dataclasses import dataclass, replace

import numpy as np
import numpy.typing as npt

from titrages.modeles import CoupleRedox, CourbeRedox, EspeceRedox, SautRedox

__all__ = [
    "NERNST",
    "REFERENCES",
    "ComposantRedox",
    "SolutionRedox",
    "TitrageRedox",
    "electrons_cedes",
    "fractions_redox",
]

Tableau = npt.NDArray[np.float64]

NERNST = 0.05916
"""RT ln 10 / F à 25 °C, en volts."""

REFERENCES = {"ENH": 0.0, "ECS": 0.244, "Ag/AgCl": 0.197}
"""Potentiels des électrodes de référence usuelles (V/ENH, 25 °C, KCl saturé)."""

E_MIN, E_MAX = -3.0, 4.0  # bornes de la dichotomie, en V/ENH


def _log_niveaux(couple: CoupleRedox, pH: float, s: float) -> tuple[Tableau, Tableau]:
    """Les coefficients :math:`t_j = (e_j E + L_j)/s` : ``e`` et ``L``."""
    L, cumul = [0.0], 0.0
    for E0, n, m in zip(couple.E0, couple.n, couple.H, strict=True):
        cumul -= n * (E0 - s * m / n * pH)
        L.append(cumul)
    e = np.concatenate(([0.0], np.cumsum(couple.n))).astype(np.float64)
    return e, np.asarray(L, dtype=np.float64)


def fractions_redox(couple: CoupleRedox, E: npt.ArrayLike, pH: float = 0.0, s: float = NERNST) -> Tableau:
    """Les fractions des formes d'un système sans forme dimérique, de forme ``(N + 1, *E.shape)``."""
    if couple.stoechiometrie != (1, 1):
        raise ValueError(f"{couple.nom} : fractions dépendant de la concentration, voir electrons_cedes")
    E_ = np.asarray(E, dtype=np.float64)
    e, L = _log_niveaux(couple, pH, s)
    forme = (-1,) + (1,) * E_.ndim
    t = (e.reshape(forme) * E_ + L.reshape(forme)) / s
    t = t - t.max(axis=0)
    x = np.power(10.0, t)
    return np.asarray(x / x.sum(axis=0), dtype=np.float64)


def _racine_1_plus(lx: Tableau) -> Tableau:
    """:math:`\\sqrt{1 + 10^{lx}}` sans débordement : :math:`10^{lx/2}\\sqrt{1 + 10^{-lx}}` au-delà de 0."""
    lx = np.clip(lx, -600, 600)
    grand = np.power(10.0, np.minimum(lx, 600) / 2) * np.sqrt(1 + np.power(10.0, -np.abs(lx)))
    petit = np.sqrt(1 + np.power(10.0, np.minimum(lx, 0)))
    return np.asarray(np.where(lx > 0, grand, petit), dtype=np.float64)


def electrons_cedes(
    couple: CoupleRedox,
    E: npt.ArrayLike,
    T: npt.ArrayLike = 1.0,
    pH: float = 0.0,
    s: float = NERNST,
) -> Tableau:
    """:math:`Y` : le nombre moyen d'électrons cédés par unité, à ce potentiel.

    ``T`` (en unités de l'élément) n'intervient que pour un couple dimérique.
    """
    E_ = np.asarray(E, dtype=np.float64)
    if couple.stoechiometrie == (1, 1):
        e, _ = _log_niveaux(couple, pH, s)
        x = fractions_redox(couple, E_, pH, s)
        forme = (-1,) + (1,) * E_.ndim
        return np.asarray((e.reshape(forme) * x).sum(axis=0), dtype=np.float64)
    T_ = np.maximum(np.asarray(T, dtype=np.float64), 1e-300)
    n, m, E0 = couple.n[0], couple.H[0], couple.E0[0]
    logK = n * (E_ - (E0 - s * m / n * pH)) / s
    e1 = couple.electrons[1]
    if couple.stoechiometrie == (1, 2):  # Ox dimère : o = K r², 2o + r = T
        r = 2 * T_ / (1 + _racine_1_plus(np.log10(8 * T_) + logK))
        oxyde = 1 - r / T_
    else:  # Red dimère : o² = K r, o + 2r = T
        o = 2 * T_ / (1 + _racine_1_plus(np.log10(8 * T_) - logK))
        oxyde = o / T_
    return np.asarray(e1 * oxyde, dtype=np.float64)


def _dichotomie(f: Callable[[Tableau], Tableau], forme: tuple[int, ...]) -> Tableau:
    a = np.full(forme, E_MIN)
    b = np.full(forme, E_MAX)
    for _ in range(64):
        m = (a + b) / 2
        positif = f(m) > 0
        b = np.where(positif, m, b)
        a = np.where(positif, a, m)
    return np.asarray((a + b) / 2, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class ComposantRedox:
    """Une espèce rédox et sa concentration apportée (mol/L de l'espèce, pas d'unités)."""

    espece: EspeceRedox
    concentration: float

    def __post_init__(self) -> None:
        if self.concentration < 0:
            raise ValueError(f"{self.espece} : concentration négative")

    @property
    def unites(self) -> float:
        """La concentration en unités de l'élément : 2 × [I₂] pour le diiode."""
        return self.concentration * self.espece.couple.unites[self.espece.q]

    @property
    def electrons(self) -> float:
        """Les électrons déjà cédés, en mol/L."""
        return self.unites * self.espece.couple.electrons[self.espece.q]


Bilan = dict[CoupleRedox, tuple[float, float]]


def _bilan(composants: tuple[ComposantRedox, ...]) -> Bilan:
    """Par couple : la concentration en unités, et les électrons introduits."""
    bilan: Bilan = {}
    for c in composants:
        T, intro = bilan.get(c.espece.couple, (0.0, 0.0))
        bilan[c.espece.couple] = (T + c.unites, intro + c.electrons)
    return bilan


@dataclass(frozen=True, slots=True)
class SolutionRedox:
    """Une solution d'espèces rédox, en milieu tamponné au ``pH`` donné.

    >>> from titrages import CoupleRedox
    >>> fer = CoupleRedox("fer", (0.77,), (1,), ("Fe2+", "Fe3+"))
    >>> round(SolutionRedox().avec(fer.red, 0.1).avec(fer.ox, 0.1).E(), 3)
    0.77
    """

    composants: tuple[ComposantRedox, ...] = ()
    pH: float = 0.0
    pente_nernst: float = NERNST

    def avec(self, espece: EspeceRedox, concentration: float) -> "SolutionRedox":
        return replace(self, composants=(*self.composants, ComposantRedox(espece, concentration)))

    def residu(self, E: npt.ArrayLike) -> Tableau:
        """:math:`\\sum_i T_i (Y_i - e_{q_i})`, croissant, nul au potentiel de la solution."""
        E_ = np.asarray(E, dtype=np.float64)
        total = np.zeros(E_.shape)
        for couple, (T, intro) in _bilan(self.composants).items():
            total = total + T * electrons_cedes(couple, E_, T, self.pH, self.pente_nernst) - intro
        return total

    def E(self) -> float:
        """Le potentiel de la solution (V/ENH)."""
        if not self.composants:
            raise ValueError("une solution sans espèce rédox n'a pas de potentiel défini")
        return float(_dichotomie(self.residu, ()))

    def concentrations(self, E: float | None = None) -> dict[str, float]:
        """La concentration de chaque forme (mol/L de la forme)."""
        potentiel = self.E() if E is None else E
        resultat: dict[str, float] = {}
        for couple, (T, _) in _bilan(self.composants).items():
            if couple.stoechiometrie == (1, 1):
                x = fractions_redox(couple, potentiel, self.pH, self.pente_nernst)
                valeurs = [T * float(xj) for xj in x]
            else:
                f = (
                    float(electrons_cedes(couple, potentiel, T, self.pH, self.pente_nernst))
                    / couple.electrons[1]
                )
                red, ox = couple.unites
                valeurs = [T * (1 - f) / red, T * f / ox]
            for forme, v in zip(couple.formes, valeurs, strict=True):
                resultat[forme] = resultat.get(forme, 0.0) + v
        return resultat


@dataclass(frozen=True, slots=True)
class TitrageRedox:
    """Le titrage potentiométrique de ``volume_becher`` du ``becher`` par la ``burette``.

    ``reference`` est l'électrode de référence de la mesure, par son nom
    (``"ECS"``, ``"Ag/AgCl"``, ``"ENH"``) ou son potentiel en V/ENH : les potentiels
    rendus sont ceux qu'affiche le voltmètre, :math:`E - E_\\text{réf}`.

    >>> from titrages import CoupleRedox
    >>> fer = CoupleRedox("fer", (0.77,), (1,), ("Fe2+", "Fe3+"))
    >>> cerium = CoupleRedox("cérium", (1.72,), (1,), ("Ce3+", "Ce4+"))
    >>> t = TitrageRedox(SolutionRedox().avec(fer.red, 0.1), SolutionRedox().avec(cerium.ox, 0.1), 10)
    >>> [round(float(e), 3) for e in t.E([5, 10, 20])]
    [0.77, 1.245, 1.72]
    """

    becher: SolutionRedox
    burette: SolutionRedox
    volume_becher: float
    reference: float | str = "ENH"

    def __post_init__(self) -> None:
        if self.volume_becher <= 0:
            raise ValueError("le volume du bécher est positif")
        if (self.becher.pH, self.becher.pente_nernst) != (self.burette.pH, self.burette.pente_nernst):
            raise ValueError("le bécher et la burette n'ont pas le même pH ni la même température")
        if isinstance(self.reference, str) and self.reference not in REFERENCES:
            raise ValueError(f"référence inconnue : {self.reference} ({', '.join(REFERENCES)})")

    @property
    def E_reference(self) -> float:
        return REFERENCES[self.reference] if isinstance(self.reference, str) else float(self.reference)

    def _residu(self, E: Tableau, V: Tableau) -> Tableau:
        """Le résidu du mélange, multiplié par :math:`V_0 + V`."""
        V0 = self.volume_becher
        b, bu = _bilan(self.becher.composants), _bilan(self.burette.composants)
        total = np.zeros(np.broadcast(E, V).shape)
        for couple in {**b, **bu}:
            Tb, Ib = b.get(couple, (0.0, 0.0))
            Tu, Iu = bu.get(couple, (0.0, 0.0))
            T = (V0 * Tb + V * Tu) / (V0 + V)
            Y = electrons_cedes(couple, E, T, self.becher.pH, self.becher.pente_nernst)
            total = total + (V0 * Tb + V * Tu) * Y - (V0 * Ib + V * Iu)
        return total

    def _E_enh(self, volume: npt.ArrayLike) -> Tableau:
        V = np.asarray(volume, dtype=np.float64)
        if np.any(V < 0):
            raise ValueError("un volume versé est positif")
        return _dichotomie(lambda E: self._residu(E, V), V.shape)

    def E(self, volume: npt.ArrayLike) -> Tableau:
        """Le potentiel mesuré après avoir versé ``volume`` (V, par rapport à la référence)."""
        return np.asarray(self._E_enh(volume) - self.E_reference, dtype=np.float64)

    def pente(self, volume: npt.ArrayLike) -> Tableau:
        """:math:`\\mathrm{d}E/\\mathrm{d}V`, par dérivation implicite du résidu.

        :math:`g(E, V) = 0` donne :math:`\\mathrm{d}E/\\mathrm{d}V = -g_V/g_E` ; les deux
        dérivées partielles se prennent par différences centrées, sur une fonction
        lisse et en un point où elle est exactement nulle.
        """
        V = np.asarray(volume, dtype=np.float64)
        E = self._E_enh(V)
        hE = 1e-7
        hV = 1e-7 * np.maximum(V, self.volume_becher)
        gE = (self._residu(E + hE, V) - self._residu(E - hE, V)) / (2 * hE)
        gV = (self._residu(E, V + hV) - self._residu(E, np.maximum(V - hV, 0))) / (
            V + hV - np.maximum(V - hV, 0)
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(-gV / gE, dtype=np.float64)

    # --- le sens direct : V(E) --------------------------------------------------

    def volume(self, E: npt.ArrayLike) -> Tableau:
        """Le volume versé qui amène au potentiel mesuré ``E`` ; ``nan`` si aucun.

        Trois cas, selon les couples dimériques présents :

        - **aucun** : chaque résidu est linéaire en la concentration, et
          :math:`V = -V_0 R_\\text{bécher}(E)/R_\\text{burette}(E)`, comme en acide-base ;
        - **un seul** : avec :math:`W = (V_0+V)\\,u`, où :math:`u` est la concentration de
          sa forme monomère, le bilan d'électrons est linéaire en :math:`W` et :math:`V`,
          et la conservation de l'élément,
          :math:`(V_0 T_b + V T_u - W)(V_0 + V) = 2\\kappa W^2`, est **du second degré en**
          :math:`V` : on en prend la racine physique ;
        - **plusieurs** : l'équation est de degré plus élevé ; on la résout en :math:`V`
          par dichotomie, :math:`E(V)` étant monotone.
        """
        E_ = np.asarray(E, dtype=np.float64) + self.E_reference
        V0, pH, s = self.volume_becher, self.becher.pH, self.becher.pente_nernst
        b, bu = _bilan(self.becher.composants), _bilan(self.burette.composants)
        couples = list({**b, **bu})
        dimeres = [c for c in couples if c.stoechiometrie != (1, 1)]
        if len(dimeres) > 1:
            return self._volume_dichotomie(E_)
        A = np.zeros(E_.shape)
        B = np.zeros(E_.shape)
        for c in couples:
            if c in dimeres:
                continue
            Y = electrons_cedes(c, E_, 1.0, pH, s)
            Tb, Ib = b.get(c, (0.0, 0.0))
            Tu, Iu = bu.get(c, (0.0, 0.0))
            A = A + V0 * (Tb * Y - Ib)
            B = B + (Tu * Y - Iu)
        with np.errstate(divide="ignore", invalid="ignore"):
            V = _volume_un_dimere(dimeres[0], E_, A, B, b, bu, V0, pH, s) if dimeres else -A / B
        return np.asarray(np.where(np.isfinite(V) & (V >= 0), V, np.nan), dtype=np.float64)

    def _volume_dichotomie(self, E: Tableau) -> Tableau:
        """V(E) pour plusieurs dimères : dichotomie en log V sur [10⁻¹² V₀, 10⁶ V₀]."""
        sens = 1.0 if float(self._E_enh(self.volume_becher)) > float(self._E_enh(0.0)) else -1.0
        a = np.full(E.shape, np.log10(self.volume_becher) - 12)
        c = np.full(E.shape, np.log10(self.volume_becher) + 6)
        for _ in range(80):
            m = (a + c) / 2
            # au-delà du bon volume, le potentiel du mélange dépasse E (titrant oxydant)
            trop = sens * self._residu(E, 10**m) < 0
            c = np.where(trop, m, c)
            a = np.where(trop, a, m)
        V = 10 ** ((a + c) / 2)
        ok = np.abs(self._E_enh(V) - E) < 1e-6
        return np.asarray(np.where(ok, V, np.nan), dtype=np.float64)

    def volume_initial(self, volume_max: float) -> float:
        """Le premier volume où le potentiel est défini.

        Un bécher qui ne contient qu'une forme de chaque couple — du Fe²⁺ sans Fe³⁺ —
        n'a pas de potentiel : Nernst le place à :math:`-\\infty`, et le calcul butte
        sur la borne de la dichotomie. La courbe commence alors au millième du volume
        tracé, là où la première goutte a fixé le rapport des formes.
        """
        E0 = float(self._E_enh(0.0))
        return volume_max * 1e-3 if min(E0 - E_MIN, E_MAX - E0) < 1e-6 else 0.0

    def courbe(self, volume_max: float, points: int = 400) -> CourbeRedox:
        """La courbe jusqu'à ``volume_max``, points exacts à longueur d'arc constante.

        Elle commence à :meth:`volume_initial` : 0, sauf si le potentiel initial n'est
        pas défini.
        """
        if volume_max <= 0:
            raise ValueError("volume_max est positif")
        debut = self.volume_initial(volume_max)
        grille = np.linspace(debut, volume_max, 40 * points)
        E = self._E_enh(grille)
        dE = max(float(E.max() - E.min()), 1e-12)
        pas = np.hypot(np.diff(grille) / volume_max, np.diff(E) / dE)
        s = np.concatenate(([0.0], np.cumsum(pas)))
        V = np.interp(np.linspace(0.0, s[-1], points), s, grille)
        V[0], V[-1] = debut, volume_max
        return CourbeRedox(V, self.E(V), self.pente(V))

    def sauts(self, volume_max: float, points: int = 4000) -> tuple[SautRedox, ...]:
        """Les points d'inflexion : les maxima locaux de :math:`|\\mathrm{d}E/\\mathrm{d}V|`."""
        grille = np.linspace(self.volume_initial(volume_max), volume_max, points)
        E = self._E_enh(grille)
        d = np.abs(np.diff(E)) / np.diff(grille)
        maxima = np.flatnonzero((d[1:-1] > d[:-2]) & (d[1:-1] >= d[2:])) + 1
        resultat = []
        for i in maxima.tolist():
            a, b = float(grille[max(i - 1, 0)]), float(grille[min(i + 2, points - 1)])
            V = _section_doree(lambda v: -float(abs(self.pente(v))), a, b)
            if 0 < V < volume_max:
                resultat.append(SautRedox(V, float(self.E(V)), float(self.pente(V))))
        return tuple(resultat)

    def volumes_equivalents(self) -> tuple[float, ...]:
        """Les équivalences stœchiométriques : les électrons échangés, rangés par potentiel.

        Le titrant cède ou accepte tous ses électrons, jusqu'à sa forme extrême ; les
        espèces du bécher réagissent étape par étape, la plus réductrice d'abord pour un
        titrant oxydant, la plus oxydante d'abord pour un titrant réducteur.
        """
        oxydant = self.burette.E() > self.becher.E()
        capacite = 0.0
        for c in self.burette.composants:
            e = c.espece.couple.electrons
            capacite += c.unites * (e[c.espece.q] if oxydant else e[-1] - e[c.espece.q])
        if capacite <= 0:
            raise ValueError("la burette n'apporte aucun électron à échanger")
        pH, s = self.becher.pH, self.becher.pente_nernst
        sites: list[tuple[float, float]] = []
        for c in self.becher.composants:
            couple, q = c.espece.couple, c.espece.q
            e = couple.electrons
            etapes = range(q + 1, couple.N + 1) if oxydant else range(q, 0, -1)
            for j in etapes:
                E0p = couple.E0[j - 1] - s * couple.H[j - 1] / couple.n[j - 1] * pH
                sites.append((E0p, c.unites * (e[j] - e[j - 1])))
        sites.sort(key=lambda t: t[0], reverse=not oxydant)
        volumes, cumul = [], 0.0
        for _, n in sites:
            cumul += n * self.volume_becher
            volumes.append(cumul / capacite)
        return tuple(volumes)


def _volume_un_dimere(
    d: CoupleRedox,
    E: Tableau,
    A: Tableau,
    B: Tableau,
    b: Bilan,
    bu: Bilan,
    V0: float,
    pH: float,
    s: float,
) -> Tableau:
    """La racine physique du trinôme en V, pour un seul couple dimérique.

    ``A + B V`` est la part des couples monomères du bilan d'électrons (× (V₀+V)).
    Pour le dimère, :math:`u` est la forme monomère (Red pour (1, 2), Ox pour (2, 1)),
    :math:`T = u + 2\\kappa u^2` avec :math:`\\kappa = K` ou :math:`1/K`.
    """
    Tb, Ib = b.get(d, (0.0, 0.0))
    Tu, Iu = bu.get(d, (0.0, 0.0))
    e1 = d.electrons[1]
    n, m, E0 = d.n[0], d.H[0], d.E0[0]
    logK = n * (E - (E0 - s * m / n * pH)) / s
    if d.stoechiometrie == (1, 2):  # u = Red, unités oxydées × (V0+V) = V0 Tb + V Tu - W
        w0 = (A + V0 * (e1 * Tb - Ib)) / e1
        w1 = (B + e1 * Tu - Iu) / e1
        kappa = np.power(10.0, np.clip(logK, -300, 300))
    else:  # u = Ox, unités oxydées × (V0+V) = W
        w0 = -(A - V0 * Ib) / e1
        w1 = -(B - Iu) / e1
        kappa = np.power(10.0, np.clip(-logK, -300, 300))
    d0, d1 = V0 * Tb - w0, Tu - w1
    c2 = d1 - 2 * kappa * w1**2
    c1 = d0 + V0 * d1 - 4 * kappa * w0 * w1
    c0 = V0 * d0 - 2 * kappa * w0**2
    disc = np.sqrt(np.maximum(c1**2 - 4 * c2 * c0, 0))
    q = -(c1 + np.where(c1 >= 0, disc, -disc)) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        racines = np.stack([q / c2, c0 / q, -c0 / c1])  # la dernière si c2 = 0
    W = w0 + w1 * racines
    reste = V0 * Tb + Tu * racines - W
    tol = 1e-12 * (V0 * Tb + Tu * np.abs(racines) + 1e-300)
    ok = np.isfinite(racines) & (racines >= 0) & (W + tol >= 0) & (reste + tol >= 0)
    return np.asarray(np.where(ok, racines, np.inf).min(axis=0), dtype=np.float64)


def _section_doree(f: Callable[[float], float], a: float, b: float, iterations: int = 80) -> float:
    r = (np.sqrt(5) - 1) / 2
    c, d = b - r * (b - a), a + r * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(iterations):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - r * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + r * (b - a)
            fd = f(d)
    return float((a + b) / 2)
