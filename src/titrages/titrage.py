"""Le titrage : un bécher, une burette, et la dilution entre les deux.

On verse un volume :math:`V` de la burette dans un volume :math:`V_0` du bécher.
Chaque concentration du mélange est la moyenne pondérée des deux, et
l'électroneutralité du mélange, multipliée par :math:`V_0 + V`, se sépare :

.. math::

    V_0\\,R_\\text{bécher}(\\mathrm{pH}) + V\\,R_\\text{burette}(\\mathrm{pH}) = 0
    \\quad\\Longrightarrow\\quad
    V(\\mathrm{pH}) = -V_0\\,\\frac{R_\\text{bécher}(\\mathrm{pH})}{R_\\text{burette}(\\mathrm{pH})}

où :math:`R` est le résidu de l'électroneutralité de chaque solution **telle
qu'elle est préparée**, sans dilution (:meth:`Solution.residu`). C'est la formule de
la fiche, :math:`n = (Y_1 + \\delta/C_T)/(1 - \\delta/C_B)`, rendue générale : un
bécher et une burette quelconques, fort ou faible, mélange ou non. Le volume est une
fonction **explicite** du pH, et n'a de sens qu'entre le pH initial (numérateur nul)
et le pH de la solution titrante (dénominateur nul, asymptote verticale).

Négliger la dilution revient à retirer :math:`\\delta` du dénominateur : la burette
apporte ses espèces, pas son eau (``dilution=False``).
"""

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from titrages.equilibre import MARGE_PH, minimum_lisse, racine_croissante
from titrages.modeles import Courbe, CourbeConductimetrique, Saut, Tangentes
from titrages.solution import LAMBDA_H3O, LAMBDA_HO, Solution

__all__ = ["Titrage"]

Tableau = npt.NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class Titrage:
    """Le titrage de ``volume_becher`` du ``becher`` par la ``burette``.

    Les volumes s'expriment dans l'unité qu'on veut, pourvu qu'elle soit la même
    partout : c'est un rapport de volumes qui intervient.

    >>> from titrages import Couple
    >>> becher = Solution().avec(Couple("acide éthanoïque", (4.76,)), 0.1)
    >>> t = Titrage(becher, Solution().avec_base_forte(0.1), volume_becher=10)
    >>> round(float(t.pH(5)), 2), round(float(t.pH(10)), 2)
    (4.76, 8.73)
    """

    becher: Solution
    burette: Solution
    volume_becher: float
    dilution: bool = True

    def __post_init__(self) -> None:
        if self.volume_becher <= 0:
            raise ValueError("le volume du bécher est positif")
        if self.becher.pKe != self.burette.pKe:
            raise ValueError("le bécher et la burette n'ont pas le même pKe")

    # --- le sens direct : V(pH) ------------------------------------------------

    def _denominateur(self, pH: npt.ArrayLike) -> Tableau:
        if self.dilution:
            return self.burette.residu(pH)
        return self.burette.charge(pH)

    def _derivee_denominateur(self, pH: npt.ArrayLike) -> Tableau:
        if self.dilution:
            return self.burette.pouvoir_tampon(pH)
        return self.burette.derivee_charge(pH)

    def volume(self, pH: npt.ArrayLike) -> Tableau:
        """Le volume versé qui amène à ce pH ; ``nan`` si aucun n'y amène."""
        N = self.becher.residu(pH)
        D = self._denominateur(pH)
        with np.errstate(divide="ignore", invalid="ignore"):
            V = -self.volume_becher * N / D
        return np.asarray(np.where(np.isfinite(V) & (V >= 0), V, np.nan), dtype=np.float64)

    def derivee_volume(self, pH: npt.ArrayLike) -> Tableau:
        """:math:`\\mathrm{d}V/\\mathrm{d\\,pH}`, analytiquement.

        En dérivant :math:`V_0 N + V D = 0` :
        :math:`V' = -(V_0 N' + V D')/D`, avec :math:`N' = \\beta_\\text{bécher}` et
        :math:`D' = \\beta_\\text{burette}`, les pouvoirs tampons des deux solutions
        **non diluées**, au pH du mélange. Le numérateur est une somme de termes positifs :
        il ne perd aucun chiffre, même au cœur du saut, là où la forme
        :math:`V_0(ND' - N'D)/D^2` soustrairait deux nombres presque égaux.
        """
        D = self._denominateur(pH)
        with np.errstate(divide="ignore", invalid="ignore"):
            V = -self.volume_becher * self.becher.residu(pH) / D
        dN = self.becher.pouvoir_tampon(pH)
        dD = self._derivee_denominateur(pH)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(-(self.volume_becher * dN + V * dD) / D, dtype=np.float64)

    def pente(self, pH: npt.ArrayLike) -> Tableau:
        """:math:`\\mathrm{d\\,pH}/\\mathrm{d}V` au point de la courbe de ce pH."""
        with np.errstate(divide="ignore"):
            return np.asarray(1 / self.derivee_volume(pH), dtype=np.float64)

    # --- le sens inverse : pH(V) -----------------------------------------------

    def pH(self, volume: npt.ArrayLike) -> Tableau:
        """Le pH après avoir versé ``volume``, par dichotomie.

        :math:`V_0 N + V D` est croissant en pH pour tout :math:`V \\ge 0` : la racine
        est unique.
        """
        V = np.asarray(volume, dtype=np.float64)
        if np.any(V < 0):
            raise ValueError("un volume versé est positif")

        def g(p: Tableau) -> Tableau:
            return self.volume_becher * self.becher.residu(p) + V * self._denominateur(p)

        pKe = self.becher.pKe
        return racine_croissante(g, np.full(V.shape, -MARGE_PH), np.full(V.shape, pKe + MARGE_PH))

    @property
    def pH_initial(self) -> float:
        return self.becher.pH()

    @property
    def sens(self) -> int:
        """+1 si le pH monte (titrant basique), -1 s'il descend."""
        return 1 if float(self.pH(1e-6 * self.volume_becher)) >= self.pH_initial else -1

    # --- la courbe -------------------------------------------------------------

    def courbe(self, volume_max: float, points: int = 400) -> Courbe:
        """La courbe de 0 à ``volume_max``, en ``points`` points exacts.

        Les points sont répartis à **longueur d'arc constante**, dans le repère où
        :math:`V/V_\\text{max}` et le pH normalisé ont la même échelle : un pas
        régulier en volume manquerait le saut, un pas régulier en pH laisserait les
        zones tampons presque vides. Chaque point est ensuite recalculé exactement,
        par :meth:`volume` : l'échantillonnage ne coûte aucune précision.
        """
        if volume_max <= 0:
            raise ValueError("volume_max est positif")
        pH0 = self.pH_initial
        pH1 = float(self.pH(volume_max))
        grille = np.linspace(pH0, pH1, 40 * points)
        V = np.nan_to_num(self.volume(grille), nan=0.0)
        V[0] = 0.0
        V[-1] = volume_max
        pas = np.hypot(np.diff(V) / volume_max, np.diff(grille) / max(abs(pH1 - pH0), 1e-12))
        s = np.concatenate(([0.0], np.cumsum(pas)))
        pH_ = np.interp(np.linspace(0.0, s[-1], points), s, grille)
        volumes = np.nan_to_num(self.volume(pH_), nan=0.0)
        volumes[0], volumes[-1] = 0.0, volume_max
        return Courbe(volumes, np.asarray(pH_, dtype=np.float64), self.pente(pH_))

    # --- les équivalences ------------------------------------------------------

    def sauts(self, volume_max: float, points: int = 4000) -> tuple[Saut, ...]:
        """Les points d'inflexion de la courbe entre 0 et ``volume_max``.

        Ce sont les maxima locaux de la pente, donc les minima de
        :math:`|\\mathrm{d}V/\\mathrm{d\\,pH}|`, repérés sur une grille puis affinés
        par le zéro de leur dérivée (:func:`~titrages.equilibre.minimum_lisse`), à la
        précision de la machine. Une acidité trop faible ou deux pKa trop proches ne
        donnent pas de saut : c'est justement ce que cette liste permet de voir, face à
        :meth:`volumes_equivalents`.
        """
        pH0 = self.pH_initial
        pH1 = float(self.pH(volume_max))
        grille = np.linspace(pH0, pH1, points)
        d = np.abs(self.derivee_volume(grille))
        interieurs = np.flatnonzero((d[1:-1] < d[:-2]) & (d[1:-1] <= d[2:])) + 1
        resultat = []
        for i in interieurs:
            p = minimum_lisse(
                lambda x: float(np.abs(self.derivee_volume(x))), grille[i - 1], grille[i + 1], 1e-6
            )
            V = float(self.volume(p))
            if np.isfinite(V) and 0 < V < volume_max:
                resultat.append(Saut(V, p, float(self.pente(p))))
        return tuple(resultat)

    def volumes_equivalents(self) -> tuple[float, ...]:
        """Les équivalences stœchiométriques, pour un titrant fort.

        L'acide ou la base forts que le bécher contient déjà réagissent d'abord, en
        totalité, avec les espèces faibles : un excès de base forte enlève les protons
        des sites acides du plus fort au plus faible, et chaque proton enlevé devient un
        site basique au même pKa (la soude transforme l'acide éthanoïque en éthanoate,
        qu'un acide fort titre ensuite) ; réciproquement pour un excès d'acide fort.
        Ce qui reste du réactif fort passe en premier devant le titrant. Le :math:`k`-ième
        volume est celui qui apporte autant de titrant que les :math:`k` premiers sites
        rangés par pKa — croissant pour une base forte titrante, décroissant pour un
        acide fort — ; deux sites de même pKa ne font qu'une équivalence. Avec un titrant
        faible, l'équivalence n'a plus de définition aussi simple : :meth:`sauts` reste
        disponible.

        >>> from titrages import solution
        >>> t = Titrage(solution({"CH3COOH": 0.05, "NaOH": 0.1}), solution({"HCl": 0.1}), 10)
        >>> t.volumes_equivalents()
        (5.0, 10.0)
        """
        if self.burette.composants:
            raise ValueError(
                "équivalences stœchiométriques définies pour un titrant fort ; "
                "voir sauts() pour un titrant faible"
            )
        C = self.burette.exces_base
        if C == 0:
            raise ValueError("la burette ne contient ni acide ni base")
        acides, bases = _sites(self.becher, self.volume_becher)
        sites = sorted(acides) if C > 0 else sorted(bases, reverse=True)
        volumes: list[float] = []
        cumul, precedent = 0.0, math.nan
        for pKa, n in sites:
            cumul += n
            if pKa == precedent:
                volumes[-1] = cumul / abs(C)
            else:
                volumes.append(cumul / abs(C))
            precedent = pKa
        return tuple(volumes)

    # --- la conductimétrie -------------------------------------------------------

    def conductivite(self, volume: npt.ArrayLike) -> Tableau:
        """La conductivité du mélange (S/m) après avoir versé ``volume``.

        La part des ions est linéaire en concentrations : celle du bécher et celle de
        la burette, calculées au pH du mélange, se pondèrent par les volumes. S'y
        ajoutent H3O+ et HO-, fixés par ce pH.
        """
        V = np.asarray(volume, dtype=np.float64)
        pH = self.pH(V)
        V0 = self.volume_becher
        ions = (V0 * self.becher.conductivite_ions(pH) + V * self.burette.conductivite_ions(pH)) / (V0 + V)
        eau = LAMBDA_H3O * np.power(10.0, -pH) + LAMBDA_HO * np.power(10.0, pH - self.becher.pKe)
        return np.asarray(ions + eau, dtype=np.float64)

    def courbe_conductimetrique(
        self, volume_max: float, points: int = 200, corrigee: bool = False
    ) -> CourbeConductimetrique:
        """:math:`\\sigma(V)` à pas de volume constant : la courbe est faite de segments.

        ``corrigee`` multiplie par :math:`(V_0 + V)/V_0` : la correction de dilution
        qui rend les segments rectilignes, celle qu'on applique en TP.
        """
        V = np.linspace(0.0, volume_max, points)
        sigma = self.conductivite(V)
        if corrigee:
            sigma = sigma * (self.volume_becher + V) / self.volume_becher
        return CourbeConductimetrique(V, sigma, self.pH(V))

    def tangentes(self, volume_avant: float, volume_max: float) -> Tangentes:
        """La méthode des tangentes parallèles, à partir d'un point avant le saut.

        On prend la tangente en ``volume_avant``, on cherche après le saut suivant le
        point de même pente, et la parallèle équidistante coupe la courbe au point
        cherché. Tout est calculé exactement sur la courbe : l'écart à
        :meth:`volumes_equivalents` est celui de la méthode, pas celui du calcul.
        """
        sauts = [s for s in self.sauts(volume_max) if s.volume > volume_avant]
        if not sauts:
            raise ValueError(f"aucun saut après {volume_avant} avant {volume_max}")
        saut = sauts[0]
        pH1 = float(self.pH(volume_avant))
        s = float(self.pente(pH1))
        pHmax = float(self.pH(volume_max))

        def ecart_pente(p: float) -> float:
            return abs(float(self.pente(p))) - abs(s)

        grille = np.linspace(saut.pH, pHmax, 2000)
        valeurs = np.array([ecart_pente(p) for p in grille])
        negatifs = np.flatnonzero(valeurs <= 0)
        if not len(negatifs) or negatifs[0] == 0:
            raise ValueError("pas de point de même pente après le saut : augmenter volume_max")
        i = int(negatifs[0])
        pH2 = _changement_de_signe(ecart_pente, float(grille[i - 1]), float(grille[i]))
        V1, V2 = volume_avant, float(self.volume(pH2))
        Vm, pHm = (V1 + V2) / 2, (pH1 + pH2) / 2

        def mediane(p: float) -> float:
            return (p - pHm) - s * (float(self.volume(p)) - Vm)

        pH = _changement_de_signe(mediane, pH1, pH2)
        return Tangentes(V1, pH1, V2, pH2, s, float(self.volume(pH)), pH)


def _sites(becher: Solution, V0: float) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Les protons que le bécher peut céder (sites acides) et accepter (sites basiques).

    Chaque site est ``(pKa, quantité)``, la quantité en mol pour ``V0`` dans l'unité des
    volumes. L'acide ou la base forts du bécher ont déjà réagi : un excès de base forte
    consomme les sites acides par pKa croissant et les transforme en sites basiques de
    même pKa, un excès d'acide fort consomme les sites basiques par pKa décroissant et
    les transforme en sites acides ; ce qui en reste est un site à pKa infini.
    """
    acides: list[list[float]] = []
    bases: list[list[float]] = []
    for c in becher.composants:
        pKa, q, N, n = c.espece.couple.pKa, c.espece.q, c.espece.couple.N, c.concentration * V0
        acides.extend([pKa[j], n] for j in range(q, N))
        bases.extend([pKa[j - 1], n] for j in range(q, 0, -1))
    b = becher.exces_base * V0
    source, cible, sens = (acides, bases, 1.0) if b > 0 else (bases, acides, -1.0)
    reste = abs(b)
    for site in sorted(source, key=lambda s: sens * s[0]):
        if reste <= 0:
            break
        pris = min(reste, site[1])
        site[1] -= pris
        reste -= pris
        cible.append([site[0], pris])
    if reste > 0:
        cible.append([sens * math.inf, reste])
    return [(s[0], s[1]) for s in acides if s[1] > 0], [(s[0], s[1]) for s in bases if s[1] > 0]


def _changement_de_signe(f: Callable[[float], float], a: float, b: float, iterations: int = 64) -> float:
    """Une racine de ``f`` entre ``a`` et ``b``, où ``f`` change de signe."""
    fa = f(a)
    for _ in range(iterations):
        m = (a + b) / 2
        fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2
