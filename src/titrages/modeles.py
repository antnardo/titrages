"""Les objets du domaine : couples, espèces, indicateurs, et ce que rend un calcul.

Un **couple** est décrit comme un polyacide :math:`\\mathrm{H}_N\\mathrm{A}` et ses
:math:`N+1` formes, de la plus acide (forme 0) à la plus basique (forme :math:`N`).
Une base faible n'a donc rien de particulier : l'ammoniac est la forme 1 du couple
:math:`\\mathrm{NH_4^+/NH_3}`. Ce choix unique est ce qui permet de traiter d'un seul
tenant les acides, les bases, les ampholytes et leurs mélanges.

Une **espèce** est une forme d'un couple, celle sous laquelle on l'introduit : son
indice :math:`q` est le nombre de protons qu'elle a déjà cédés. C'est la seule
information nécessaire sur les contre-ions — une solution préparée est neutre, les
ions spectateurs portent donc exactement la charge :math:`q - z_0` qui manque, et la
charge :math:`z_0` de la forme la plus acide disparaît de l'électroneutralité.

Un **couple rédox** suit la même logique, les électrons remplaçant les protons : ses
formes vont de la plus réduite (forme 0) à la plus oxydée, et une espèce est la forme
sous laquelle on l'introduit (:class:`CoupleRedox`, :class:`EspeceRedox`).

Un **complexe** (:class:`Complexe`) est la suite des complexes successifs
:math:`\\mathrm{ML}, \\mathrm{ML_2}, \\dots` d'un métal avec le ligand titrant, donnée par
ses constantes de formation globales :math:`\\beta_j`. Un **précipité**
(:class:`Precipite`) est un solide :math:`\\mathrm{P}_\\nu\\mathrm{X}_\\mu` formé par l'ion
commun P — celui qu'on titre ou qu'on verse, :math:`\\mathrm{Ag^+}` en argentimétrie — et
un partenaire X, donné par son :math:`\\mathrm{p}K_s`.

Les objets de ce module sont des données : aucun calcul d'équilibre n'y vit, ils sont
dans :mod:`titrages.equilibre` et :mod:`titrages.redox`.
"""

import re
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

__all__ = [
    "Complexe",
    "Couple",
    "CoupleRedox",
    "Courbe",
    "CourbeComplexometrique",
    "CourbeConductimetrique",
    "CourbePrecipitation",
    "CourbeRedox",
    "Espece",
    "EspeceRedox",
    "Evenement",
    "Indicateur",
    "Ligand",
    "Precipite",
    "Saut",
    "SautComplexometrique",
    "SautPrecipitation",
    "SautRedox",
    "Spectateur",
    "Tangentes",
    "TangentesRedox",
    "charge",
]

Tableau = npt.NDArray[np.float64]
Valeur = float | Tableau

_CHARGE_EXPOSANT = re.compile(r"\^\{?(\d*)([+-])\}?$")
_CHARGE_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)([+-])")
_CHARGE_INTERNE = re.compile(r"(?<=[A-Za-z0-9)])([+-])(?=[A-Za-z(])")


def charge(formule: str) -> int:
    """La charge qu'une formule écrit : à sa fin, et à l'intérieur pour un zwitterion.

    Une charge multiple s'écrit après ``^`` (``CO3^2-``) ; sans ``^``, un chiffre
    avant le signe n'est une charge qu'après un élément seul (``Fe3+``, ``Mg2+``) :
    dans ``NH4+`` ou ``HCO3-``, c'est un indice. ``CO32-`` est donc lu -1 : écrire
    ``CO3^2-``. Un signe seul entre deux atomes compte pour ±1 : la glycine
    ``H3N+CH2COO-`` est neutre, l'ion glycinium ``H3N+CH2COOH`` vaut +1.

    >>> [charge(f) for f in ("H2CO3", "HCO3-", "CO3^2-", "NH4+", "Fe3+", "PO4^{3-}", "Hg2^2+")]
    [0, -1, -2, 1, 3, -3, 2]
    >>> [charge(f) for f in ("H3N+CH2COOH", "H3N+CH2COO-", "H2NCH2COO-")]
    [1, 0, -1]
    """
    f = formule.strip()
    interne = sum(1 if m.group(1) == "+" else -1 for m in _CHARGE_INTERNE.finditer(f))
    return interne + _charge_finale(f)


def _charge_finale(f: str) -> int:
    m = _CHARGE_EXPOSANT.search(f)
    if m:
        valeur = int(m.group(1)) if m.group(1) else 1
        return valeur if m.group(2) == "+" else -valeur
    m = _CHARGE_ELEMENT.fullmatch(f)
    if m:
        valeur = int(m.group(2)) if m.group(2) else 1
        return valeur if m.group(3) == "+" else -valeur
    if f.endswith("+"):
        return 1
    if f.endswith("-"):
        return -1
    return 0


@dataclass(frozen=True, slots=True)
class Spectateur:
    """Un ion qui ne réagit pas : il ne compte que pour la conductivité.

    ``conductivite`` est la conductivité molaire ionique limite :math:`\\lambda^\\circ`,
    en mS·m²·mol⁻¹, de l'ion tel qu'il est écrit (pour :math:`\\mathrm{SO_4^{2-}}`,
    celle de l'ion, pas celle de :math:`\\frac12\\mathrm{SO_4^{2-}}`).
    """

    formule: str
    concentration: float
    conductivite: float | None = None


@dataclass(frozen=True, slots=True)
class Couple:
    """Un polyacide et ses formes successives.

    ``pKa`` est rangé par ordre croissant à la construction : l'indice 0 est la
    constante de la réaction :math:`\\mathrm{H}_N\\mathrm{A} = \\mathrm{H^+} +
    \\mathrm{H}_{N-1}\\mathrm{A^-}`. Des pKa donnés dans le désordre désignent le même
    couple ; les laisser tels quels donnerait des fractions fausses sans erreur.
    """

    nom: str
    pKa: tuple[float, ...]
    formes: tuple[str, ...] = ()
    noms_formes: tuple[tuple[str, ...], ...] = ()
    note: str = ""
    conductivites: tuple[float | None, ...] = ()
    """:math:`\\lambda^\\circ` de chaque forme (mS·m²·mol⁻¹), ``None`` si inconnue."""

    def __post_init__(self) -> None:
        pKa = tuple(sorted(float(p) for p in self.pKa))
        if not pKa:
            raise ValueError(f"{self.nom} : un couple a au moins un pKa")
        object.__setattr__(self, "pKa", pKa)
        if not self.formes:
            object.__setattr__(self, "formes", _formes_generiques(len(pKa)))
        if len(self.formes) != len(pKa) + 1:
            raise ValueError(
                f"{self.nom} : {len(pKa)} pKa appellent {len(pKa) + 1} formes, {len(self.formes)} données"
            )
        if self.noms_formes and len(self.noms_formes) != len(self.formes):
            raise ValueError(f"{self.nom} : un jeu de noms par forme, ou aucun")
        if self.conductivites and len(self.conductivites) != len(self.formes):
            raise ValueError(f"{self.nom} : une conductivité par forme, ou aucune")

    @classmethod
    def polyacide(cls, *pKa: float, nom: str = "") -> "Couple":
        """Un couple anonyme, pour un calcul qui n'a besoin que des pKa."""
        return cls(nom or "polyacide " + ", ".join(f"{p:g}" for p in pKa), tuple(pKa))

    @property
    def N(self) -> int:
        """Le nombre d'acidités."""
        return len(self.pKa)

    def forme(self, q: int) -> "Espece":
        """L'espèce qui a cédé ``q`` protons (0 : la plus acide)."""
        return Espece(self, q)

    @property
    def acide(self) -> "Espece":
        """La forme la plus acide."""
        return Espece(self, 0)

    @property
    def base(self) -> "Espece":
        """La forme la plus basique."""
        return Espece(self, self.N)


def _formes_generiques(N: int) -> tuple[str, ...]:
    """``H2A``, ``HA-``, ``A2-`` : les noms qu'on écrit au tableau."""
    formes = []
    for j in range(N + 1):
        h = N - j
        tete = "" if h == 0 else ("H" if h == 1 else f"H{h}")
        charge = "" if j == 0 else ("-" if j == 1 else f"{j}-")
        formes.append(f"{tete}A{charge}")
    return tuple(formes)


@dataclass(frozen=True, slots=True)
class Espece:
    """La forme ``q`` d'un couple, telle qu'on l'introduit en solution."""

    couple: Couple
    q: int = 0

    def __post_init__(self) -> None:
        if not 0 <= self.q <= self.couple.N:
            raise ValueError(f"{self.couple.nom} n'a pas de forme {self.q} (0 à {self.couple.N})")

    @property
    def formule(self) -> str:
        return self.couple.formes[self.q]

    def __str__(self) -> str:
        return self.formule


@dataclass(frozen=True, slots=True)
class Indicateur:
    """Un indicateur coloré : un monoacide dont les deux formes se voient.

    La zone de virage est celle que donnent les tables, pas :math:`\\mathrm{pK_a}\\pm 1` :
    elle dépend de l'œil et de l'intensité relative des deux couleurs.
    """

    nom: str
    pKa: float
    zone: tuple[float, float]
    couleur_acide: str
    couleur_base: str
    teinte_acide: str = "#ffffff"
    teinte_base: str = "#ffffff"
    alias: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class Courbe:
    """Une courbe de titrage échantillonnée : chaque point est exact.

    ``pente`` est :math:`\\mathrm{d}\\,\\mathrm{pH}/\\mathrm{d}V`, calculée
    analytiquement et non par différences finies.
    """

    volume: Tableau = field(repr=False)
    pH: Tableau = field(repr=False)
    pente: Tableau = field(repr=False)

    def __len__(self) -> int:
        return len(self.volume)


@dataclass(frozen=True, slots=True)
class Saut:
    """Un point d'inflexion de la courbe : là où la pente est localement maximale.

    C'est ce qu'on repère expérimentalement comme équivalence. Il ne coïncide avec
    l'équivalence stœchiométrique qu'à une approximation près, d'autant meilleure que
    le saut est grand.
    """

    volume: float
    pH: float
    pente: float


@dataclass(frozen=True, slots=True)
class Tangentes:
    """La construction de la méthode des tangentes parallèles.

    Deux points de même pente de part et d'autre du saut, et le point où la parallèle
    équidistante coupe la courbe.
    """

    volume_avant: float
    pH_avant: float
    volume_apres: float
    pH_apres: float
    pente: float
    volume: float
    pH: float

    def droite(self, volume: npt.ArrayLike, laquelle: str = "milieu") -> Valeur:
        """L'ordonnée d'une des trois droites : ``avant``, ``apres`` ou ``milieu``.

        Un ``float`` pour un volume scalaire, un tableau pour un tableau.
        """
        points = {
            "avant": (self.volume_avant, self.pH_avant),
            "apres": (self.volume_apres, self.pH_apres),
            "milieu": (
                (self.volume_avant + self.volume_apres) / 2,
                (self.pH_avant + self.pH_apres) / 2,
            ),
        }
        return _droite(points[laquelle], self.pente, volume)


def _droite(point: tuple[float, float], pente: float, volume: npt.ArrayLike) -> Valeur:
    V0, y0 = point
    t = np.asarray(y0 + pente * (np.asarray(volume, dtype=np.float64) - V0), dtype=np.float64)
    return float(t) if t.ndim == 0 else t


@dataclass(frozen=True, slots=True)
class CourbeConductimetrique:
    """Une courbe de titrage conductimétrique : :math:`\\sigma` en S/m selon le volume.

    ``pH`` accompagne chaque point : c'est lui qui fixe la répartition des formes.
    """

    volume: Tableau = field(repr=False)
    conductivite: Tableau = field(repr=False)
    pH: Tableau = field(repr=False)

    def __len__(self) -> int:
        return len(self.volume)


@dataclass(frozen=True, slots=True)
class CoupleRedox:
    """Un système rédox et ses formes, de la plus réduite (0) à la plus oxydée (:math:`N`).

    Chaque étape :math:`j` (de la forme :math:`j-1` à la forme :math:`j`) a son potentiel
    standard ``E0[j-1]`` (V/ENH), son nombre d'électrons ``n[j-1]`` et le nombre de
    protons ``H[j-1]`` que consomme la réduction :

    .. math:: \\mathrm{Ox} + m\\,\\mathrm{H^+} + n\\,e^- = \\mathrm{Red}.

    ``stoechiometrie = (a, b)`` décrit un couple à une étape dont une forme est
    dimérique : :math:`a\\,\\mathrm{Ox} + n\\,e^- = b\\,\\mathrm{Red}` — ``(1, 2)`` pour
    :math:`\\mathrm{I_2/I^-}` ou :math:`\\mathrm{Cr_2O_7^{2-}/Cr^{3+}}`, ``(2, 1)`` pour
    :math:`\\mathrm{Hg^{2+}/Hg_2^{2+}}`. Les concentrations se comptent alors en
    **unités** de l'élément : :math:`T = 2[\\mathrm{I_2}] + [\\mathrm{I^-}]`.
    """

    nom: str
    E0: tuple[float, ...]
    n: tuple[int, ...] = ()
    formes: tuple[str, ...] = ()
    H: tuple[int, ...] = ()
    stoechiometrie: tuple[int, int] = (1, 1)
    noms_formes: tuple[tuple[str, ...], ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        N = len(self.E0)
        if N == 0:
            raise ValueError(f"{self.nom} : au moins une étape")
        object.__setattr__(self, "E0", tuple(float(e) for e in self.E0))
        if not self.n:
            object.__setattr__(self, "n", (1,) * N)
        if not self.H:
            object.__setattr__(self, "H", (0,) * N)
        if not self.formes:
            object.__setattr__(self, "formes", tuple(f"R{j}" for j in range(N + 1)))
        if len(self.n) != N or len(self.H) != N or len(self.formes) != N + 1:
            raise ValueError(f"{self.nom} : {N} étapes appellent {N} n, {N} H et {N + 1} formes")
        if any(k <= 0 for k in self.n):
            raise ValueError(f"{self.nom} : un nombre d'électrons est positif")
        a, b = self.stoechiometrie
        if (a, b) not in {(1, 1), (1, 2), (2, 1)}:
            raise ValueError(f"{self.nom} : stœchiométrie (1, 1), (1, 2) ou (2, 1) seulement")
        if (a, b) != (1, 1) and N != 1:
            raise ValueError(f"{self.nom} : un couple dimérique n'a qu'une étape")

    @classmethod
    def simple(cls, E0: float, n: int = 1, nom: str = "", H: int = 0) -> "CoupleRedox":
        """Un couple à une étape, Ox + n e⁻ = Red."""
        return cls(nom or f"couple {E0:g} V", (E0,), (n,), ("Red", "Ox"), (H,))

    @property
    def N(self) -> int:
        return len(self.E0)

    @property
    def unites(self) -> tuple[int, ...]:
        """Le nombre d'unités de l'élément que porte chaque forme (2 pour I₂ face à I⁻)."""
        a, b = self.stoechiometrie
        if (a, b) == (1, 2):
            return (1, 2)
        if (a, b) == (2, 1):
            return (2, 1)
        return (1,) * (self.N + 1)

    @property
    def electrons(self) -> tuple[float, ...]:
        """Les électrons cédés par unité, depuis la forme la plus réduite, pour chaque forme."""
        a, b = self.stoechiometrie
        cumul, e = 0.0, [0.0]
        for k in self.n:
            cumul += k / (a * self.unites[-1]) if (a, b) != (1, 1) else k
            e.append(cumul)
        return tuple(e)

    def forme(self, q: int) -> "EspeceRedox":
        return EspeceRedox(self, q)

    @property
    def red(self) -> "EspeceRedox":
        return EspeceRedox(self, 0)

    @property
    def ox(self) -> "EspeceRedox":
        return EspeceRedox(self, self.N)


@dataclass(frozen=True, slots=True)
class EspeceRedox:
    """La forme ``q`` d'un couple rédox (0 : la plus réduite)."""

    couple: CoupleRedox
    q: int = 0

    def __post_init__(self) -> None:
        if not 0 <= self.q <= self.couple.N:
            raise ValueError(f"{self.couple.nom} n'a pas de forme {self.q} (0 à {self.couple.N})")

    @property
    def formule(self) -> str:
        return self.couple.formes[self.q]

    def __str__(self) -> str:
        return self.formule


@dataclass(frozen=True, slots=True)
class CourbeRedox:
    """Une courbe de titrage potentiométrique : :math:`E` en V selon le volume."""

    volume: Tableau = field(repr=False)
    E: Tableau = field(repr=False)
    pente: Tableau = field(repr=False)

    def __len__(self) -> int:
        return len(self.volume)


@dataclass(frozen=True, slots=True)
class SautRedox:
    """Un point d'inflexion d'une courbe E = f(V)."""

    volume: float
    E: float
    pente: float


@dataclass(frozen=True, slots=True)
class TangentesRedox:
    """La construction des tangentes parallèles sur une courbe E = f(V) (voir :class:`Tangentes`)."""

    volume_avant: float
    E_avant: float
    volume_apres: float
    E_apres: float
    pente: float
    volume: float
    E: float

    def droite(self, volume: npt.ArrayLike, laquelle: str = "milieu") -> Valeur:
        """L'ordonnée d'une des trois droites : ``avant``, ``apres`` ou ``milieu``."""
        points = {
            "avant": (self.volume_avant, self.E_avant),
            "apres": (self.volume_apres, self.E_apres),
            "milieu": ((self.volume_avant + self.volume_apres) / 2, (self.E_avant + self.E_apres) / 2),
        }
        return _droite(points[laquelle], self.pente, volume)


# --- complexation ------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Ligand:
    """Le ligand titrant, et ses acidités : l'EDTA est :math:`\\mathrm{H_6Y^{2+}}` à :math:`\\mathrm{Y^{4-}}`.

    Les ``pKa`` sont ceux de ses formes protonées ; c'est la forme la plus basique qui
    complexe. À pH fixé, seule la fraction :math:`\\alpha` de cette forme compte : les
    constantes de formation deviennent conditionnelles, :math:`\\beta'_j = \\beta_j\\alpha^j`.
    Un ligand sans acidité (``pKa`` vide) complexe tel quel, à tout pH. ``note`` dit,
    pour un ligand de la base, les conditions de la mesure.
    """

    nom: str = "L"
    pKa: tuple[float, ...] = ()
    note: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "pKa", tuple(sorted(float(p) for p in self.pKa)))


@dataclass(frozen=True, slots=True)
class Complexe:
    """Les complexes successifs :math:`\\mathrm{ML}_j` d'un ``metal`` avec le ligand.

    ``log_beta`` porte les constantes de formation **globales**,
    :math:`\\log\\beta_j` pour :math:`\\mathrm{M} + j\\,\\mathrm{L} = \\mathrm{ML}_j`, de
    :math:`j = 1` à :math:`n`. Le complexe calcium-EDTA est ``Complexe("Ca2+", (10.7,))``,
    les complexes amminecuivre(II) ``Complexe("Cu2+", (4.1, 7.6, 10.5, 12.6))``.

    Un complexe intermédiaire que les tables ne donnent pas — on ne connaît que
    :math:`\\beta_2` pour :math:`\\mathrm{Ag(CN)_2^-}` — a :math:`\\log\\beta = -\\infty` :
    il ne se forme pas. ``note`` dit, pour un complexe de la base, sa source et ses conditions.
    """

    metal: str
    log_beta: tuple[float, ...]
    nom: str = ""
    note: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "log_beta", tuple(float(b) for b in self.log_beta))
        if not self.log_beta:
            raise ValueError(f"{self.metal} : au moins une constante de formation")

    @property
    def n(self) -> int:
        """Le nombre maximal de ligands fixés."""
        return len(self.log_beta)


@dataclass(frozen=True, slots=True)
class CourbeComplexometrique:
    """Une courbe de titrage complexométrique : :math:`\\mathrm{pM}` selon le volume.

    ``pL`` est le paramètre de la courbe — le cologarithme du ligand libre, toutes formes
    protonées comprises —, ``pente`` est :math:`\\mathrm{d\\,pM}/\\mathrm{d}V`.
    """

    volume: Tableau = field(repr=False)
    pM: Tableau = field(repr=False)
    pL: Tableau = field(repr=False)
    pente: Tableau = field(repr=False)

    def __len__(self) -> int:
        return len(self.volume)


@dataclass(frozen=True, slots=True)
class SautComplexometrique:
    """Un point d'inflexion d'une courbe :math:`\\mathrm{pM} = f(V)`."""

    volume: float
    pM: float
    pente: float


# --- précipitation -----------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Precipite:
    """Le solide :math:`\\mathrm{P}_\\nu\\mathrm{X}_\\mu` d'un ``partenaire`` X avec l'ion commun P.

    :math:`K_s = [\\mathrm{P}]^\\nu[\\mathrm{X}]^\\mu`. Le chlorure d'argent est
    ``Precipite("Cl-", 9.75)``, le chromate d'argent ``Precipite("CrO4^2-", 11.95, nu=2)``.
    Les ``pKa`` sont ceux des formes protonées du partenaire (l'ion chromate est une base
    faible) : à pH fixé, seule sa forme la plus basique précipite, et le produit de
    solubilité devient conditionnel. ``commun`` nomme l'ion commun (``"Ag+"``) ; il fixe
    celui du titrage quand on ne le donne pas.
    """

    partenaire: str
    pKs: float
    nu: int = 1
    mu: int = 1
    pKa: tuple[float, ...] = ()
    nom: str = ""
    commun: str = ""
    note: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "pKa", tuple(sorted(float(p) for p in self.pKa)))
        if self.nu < 1 or self.mu < 1:
            raise ValueError(f"{self.partenaire} : stœchiométrie positive")


@dataclass(frozen=True, slots=True)
class CourbePrecipitation:
    """Une courbe de titrage par précipitation : le cologarithme de l'ion commun selon le volume.

    ``solides`` dit, en chaque point, les partenaires dont le solide est présent ;
    ``pente`` est :math:`\\mathrm{d\\,pP}/\\mathrm{d}V`, qui saute là où un solide apparaît
    ou disparaît.
    """

    volume: Tableau = field(repr=False)
    p: Tableau = field(repr=False)
    pente: Tableau = field(repr=False)
    solides: tuple[tuple[str, ...], ...] = field(repr=False)

    def __len__(self) -> int:
        return len(self.volume)


@dataclass(frozen=True, slots=True)
class SautPrecipitation:
    """Un point d'inflexion d'une courbe :math:`\\mathrm{pP} = f(V)`, hors des points anguleux."""

    volume: float
    p: float
    pente: float


@dataclass(frozen=True, slots=True)
class Evenement:
    """L'apparition ou la disparition d'un solide : un point anguleux de la courbe.

    C'est ce que repère la méthode de Mohr — l'apparition du chromate d'argent.
    """

    volume: float
    p: float
    partenaire: str
    apparition: bool
