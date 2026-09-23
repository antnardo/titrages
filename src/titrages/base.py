"""La base de couples et d'indicateurs, et la composition d'une solution par les noms.

La base vit dans deux fichiers TOML du paquet (``donnees/``), lisibles et corrigeables
sans toucher au code. On y cherche une espèce **par ce qu'on écrit** : le nom du
couple, le nom d'une forme ou du sel qui l'apporte, ou sa formule. La casse, les
accents, les espaces et l'écriture des charges ne comptent pas — ``CO3^2-``,
``CO32-`` et ``co3 2-`` désignent la même forme — parce qu'une recherche qui échoue
sur un accent est une recherche qu'on ne refait pas.

Chaque clé ne désigne qu'une chose : une ambiguïté dans la base est une erreur au
chargement, pas un choix silencieux à la recherche. Seule exception, deux formules
rédox que la normalisation des charges rend identiques (``VO^2+`` et ``VO2^+``) : la
formule est alors retirée de l'index, et la forme se cherche par son nom.

Les couples rédox ont leur propre index (:func:`espece_redox`, :func:`solution_redox`) :
un même ion peut y avoir un sens différent — :math:`\\mathrm{Cl^-}` est ici un spectateur,
là le réducteur du dichlore.
"""

import tomllib
import unicodedata
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from importlib import resources
from typing import Any

from titrages.equilibre import PKE
from titrages.modeles import Couple, CoupleRedox, Espece, EspeceRedox, Indicateur, Spectateur, charge
from titrages.redox import SolutionRedox
from titrages.solution import Solution

__all__ = [
    "Fort",
    "chercher",
    "conductivite_molaire",
    "couple",
    "couple_redox",
    "couples",
    "couples_redox",
    "espece",
    "espece_redox",
    "indicateur",
    "indicateurs",
    "normaliser",
    "solution",
    "solution_redox",
]


@dataclass(frozen=True, slots=True)
class Fort:
    """Un acide ou une base forts : une dissociation totale, pas un couple."""

    nature: str  # "acide" ou "base"
    formules: tuple[str, ...]
    noms: tuple[str, ...]
    contre_ion: str = ""

    @property
    def nom(self) -> str:
        return self.noms[0]


Entree = Espece | Fort


def normaliser(texte: str) -> str:
    """La clé de recherche : sans casse, accents, espaces ni marques de charge."""
    decompose = unicodedata.normalize("NFKD", texte)
    sans_accents = "".join(c for c in decompose if not unicodedata.combining(c))
    traduit = sans_accents.replace("−", "-").replace("’", "'").lower()
    return "".join(c for c in traduit if c not in " ^{}_$\\")


def _lire(fichier: str) -> dict[str, Any]:
    with resources.files("titrages").joinpath("donnees", fichier).open("rb") as flux:
        return tomllib.load(flux)


@cache
def _conductivites() -> dict[str, float]:
    return {k: float(v) for k, v in _lire("conductivites.toml")["lambda"].items()}


def conductivite_molaire(formule: str) -> float | None:
    """λ° de l'ion (mS·m²·mol⁻¹, 25 °C), ou ``None`` s'il n'est pas dans la base."""
    return _conductivites().get(formule)


@cache
def couples() -> tuple[Couple, ...]:
    """Tous les couples acide-base de la base, dans l'ordre du fichier."""
    return tuple(
        Couple(
            nom=c["nom"],
            pKa=tuple(c["pKa"]),
            formes=tuple(c["formes"]),
            noms_formes=tuple(tuple(n) for n in c["noms"]),
            note=c.get("note", ""),
            conductivites=tuple(conductivite_molaire(f) for f in c["formes"]),
        )
        for c in _lire("couples.toml")["couple"]
    )


@cache
def _forts() -> tuple[Fort, ...]:
    return tuple(
        Fort(f["nature"], tuple(f["formules"]), tuple(f["noms"]), f.get("contre_ion", ""))
        for f in _lire("couples.toml")["fort"]
    )


@cache
def _sels() -> dict[str, str]:
    return {normaliser(k): v for k, v in _lire("couples.toml")["sels"].items()}


@cache
def couples_redox() -> tuple[CoupleRedox, ...]:
    """Tous les couples rédox de la base."""
    return tuple(
        CoupleRedox(
            nom=c["nom"],
            E0=tuple(c["E0"]),
            n=tuple(c["n"]),
            formes=tuple(c["formes"]),
            H=tuple(c.get("H", (0,) * len(c["E0"]))),
            stoechiometrie=(int(c.get("stoechiometrie", (1, 1))[0]), int(c.get("stoechiometrie", (1, 1))[1])),
            noms_formes=tuple(tuple(n) for n in c["noms"]),
            note=c.get("note", ""),
        )
        for c in _lire("redox.toml")["couple"]
    )


@cache
def indicateurs() -> tuple[Indicateur, ...]:
    """Tous les indicateurs colorés de la base."""
    return tuple(
        Indicateur(
            nom=i["nom"],
            pKa=float(i["pKa"]),
            zone=(float(i["zone"][0]), float(i["zone"][1])),
            couleur_acide=i["couleurs"][0],
            couleur_base=i["couleurs"][1],
            teinte_acide=i["teintes"][0],
            teinte_base=i["teintes"][1],
            alias=tuple(i.get("alias", ())),
        )
        for i in _lire("indicateurs.toml")["indicateur"]
    )


def _ajouter(index: dict[str, Entree], cle: str, entree: Entree) -> None:
    k = normaliser(cle)
    deja = index.get(k)
    if deja is not None and deja != entree:
        raise ValueError(f"base ambiguë : « {cle} » désigne {deja} et {entree}")
    index[k] = entree


@cache
def _index() -> dict[str, Entree]:
    index: dict[str, Entree] = {}
    for c in couples():
        for q, forme in enumerate(c.formes):
            _ajouter(index, forme, c.forme(q))
            for nom in c.noms_formes[q] if c.noms_formes else ():
                _ajouter(index, nom, c.forme(q))
        # Le nom du couple désigne sa forme la plus acide — « ion ammonium » est
        # NH4+ — sauf si une forme porte ce nom : « glycine » est le zwitterion.
        index.setdefault(normaliser(c.nom), c.acide)
    for f in _forts():
        for cle in (*f.formules, *f.noms):
            _ajouter(index, cle, f)
    return index


@cache
def _index_redox() -> dict[str, EspeceRedox]:
    index: dict[str, EspeceRedox] = {}
    formules: dict[str, set[EspeceRedox]] = {}
    for c in couples_redox():
        for q, forme in enumerate(c.formes):
            formules.setdefault(normaliser(forme), set()).add(c.forme(q))
            for nom in c.noms_formes[q] if c.noms_formes else ():
                k = normaliser(nom)
                if k in index and index[k] != c.forme(q):
                    raise ValueError(f"base rédox ambiguë : « {nom} »")
                index[k] = c.forme(q)
        index.setdefault(normaliser(c.nom), c.ox)
    for k, especes in formules.items():
        if len(especes) == 1:
            index.setdefault(k, next(iter(especes)))
    return index


def espece_redox(nom: str) -> EspeceRedox:
    """L'espèce rédox que désigne ``nom`` : ``"Fe2+"``, ``"permanganate de potassium"``, ``"I2"``.

    Le nom du couple seul désigne sa forme oxydée.
    """
    try:
        return _index_redox()[normaliser(nom)]
    except KeyError:
        raise KeyError(f"« {nom} » n'est pas dans la base rédox") from None


def couple_redox(nom: str) -> CoupleRedox:
    return espece_redox(nom).couple


def _trouver(nom: str) -> Entree:
    try:
        return _index()[normaliser(nom)]
    except KeyError:
        proches = chercher(nom)
        suggestion = f" ; peut-être : {', '.join(proches[:5])}" if proches else ""
        raise KeyError(f"« {nom} » n'est pas dans la base{suggestion}") from None


def espece(nom: str) -> Espece:
    """L'espèce que désigne ``nom`` : ``"NH3"``, ``"ammoniac"``, ``"carbonate de sodium"``."""
    trouve = _trouver(nom)
    if isinstance(trouve, Fort):
        raise KeyError(f"« {nom} » est {'un acide fort' if trouve.nature == 'acide' else 'une base forte'}")
    return trouve


def couple(nom: str) -> Couple:
    """Le couple auquel appartient ce que désigne ``nom``."""
    return espece(nom).couple


def indicateur(nom: str) -> Indicateur:
    """L'indicateur de ce nom ou de cet alias (``"BBT"``, ``"hélianthine"``)."""
    cle = normaliser(nom)
    for i in indicateurs():
        if cle in {normaliser(n) for n in (i.nom, *i.alias)}:
            return i
    raise KeyError(f"« {nom} » n'est pas un indicateur de la base")


def chercher(texte: str) -> list[str]:
    """Les clés de la base qui contiennent ``texte``, pour une recherche approchée."""
    cle = normaliser(texte)
    noms: set[str] = set()
    for c in couples():
        for nom in (c.nom, *c.formes, *(n for groupe in c.noms_formes for n in groupe)):
            if cle in normaliser(nom):
                noms.add(nom)
    for f in _forts():
        noms.update(n for n in (*f.formules, *f.noms) if cle in normaliser(n))
    return sorted(noms, key=lambda n: (len(n), n))


def _spectateur(formule: str, concentration: float) -> Spectateur:
    return Spectateur(formule, concentration, conductivite_molaire(formule))


def solution(
    composition: Mapping[str, float] | Mapping[str | Espece | Couple, float] | None = None,
    pKe: float = PKE,
    contre_ions: tuple[str, str] = ("Na+", "Cl-"),
) -> Solution:
    """Une solution décrite par ce qu'on y a dissous, en mol/L.

    Les ions spectateurs s'en déduisent, pour la conductivité : celui de l'acide ou de
    la base forts (Cl⁻ de HCl, Na⁺ de la soude), celui du sel nommé (Na⁺ du carbonate
    de sodium). Une forme ionique désignée par sa formule ou son nom d'ion reçoit
    ``contre_ions`` : le premier pour un anion, le second pour un cation.

    >>> round(solution({"CH3COOH": 0.1, "CH3COO-": 0.1}).pH(), 2)
    4.76
    >>> round(solution({"acide chlorhydrique": 0.01}).pH(), 2)
    2.0
    >>> [(i.formule, i.concentration) for i in solution({"Na2CO3": 0.1}).spectateurs]
    [('Na+', 0.2)]
    """
    s = Solution(pKe=pKe)
    for cle, C in (composition or {}).items():
        if isinstance(cle, Espece | Couple):
            s = s.avec(cle, C)
            continue
        trouve = _trouver(cle)
        if isinstance(trouve, Fort):
            ion = _spectateur(trouve.contre_ion, C) if trouve.contre_ion else None
            s = s.avec_acide_fort(C, ion) if trouve.nature == "acide" else s.avec_base_forte(C, ion)
            continue
        s = s.avec(trouve, C)
        z = charge(trouve.formule)
        if z == 0:
            continue
        contre = _sels().get(normaliser(cle)) or (contre_ions[0] if z < 0 else contre_ions[1])
        s = s.avec_spectateur(_spectateur(contre, abs(z) / abs(charge(contre)) * C))
    return s


def solution_redox(
    composition: Mapping[str, float] | Mapping[str | EspeceRedox, float] | None = None,
    pH: float = 0.0,
) -> SolutionRedox:
    """Une solution rédox, décrite par les noms, en milieu tamponné au ``pH`` donné.

    >>> round(solution_redox({"Fe2+": 0.1, "Fe3+": 0.1}).E(), 2)
    0.77
    """
    s = SolutionRedox(pH=pH)
    for cle, C in (composition or {}).items():
        s = s.avec(cle if isinstance(cle, EspeceRedox) else espece_redox(cle), C)
    return s
