"""Une solution aqueuse : des espèces faibles, un excès d'acide ou de base forts.

Tout se ramène à une fonction, le **résidu de l'électroneutralité**

.. math::

    R(\\mathrm{pH}) = \\sum_i C_i\\,\\bigl(Y_{1,i}(h) - q_i\\bigr) - b + \\delta(h),

où :math:`b = [\\text{base forte}] - [\\text{acide fort}]` et :math:`q_i` la forme sous
laquelle l'espèce :math:`i` a été introduite. :math:`R` est nul au pH de la solution et
**strictement croissant** — :math:`Y_1` et :math:`\\delta` le sont — ce qui garantit
une racine unique et une dichotomie sans surprise. Sa dérivée est le pouvoir tampon :

.. math::

    \\beta = \\frac{\\mathrm{d}R}{\\mathrm{d\\,pH}}
    = \\ln 10\\,\\Bigl(\\sum_i C_i\\,(Y_{2,i} - Y_{1,i}^2) + w + h\\Bigr).

Un acide fort est traité à part, par :math:`b`, et non comme un couple de pKa très
négatif : c'est exact, et cela évite de choisir une valeur arbitraire.

**La conductivité** suit la loi de Kohlrausch, :math:`\\sigma = \\sum_i \\lambda_i^\\circ c_i`,
en solution idéale : toutes les concentrations étant explicites à pH imposé, elle
l'est aussi. Elle demande ce que le pH ignore, les **ions spectateurs** — le
:math:`\\mathrm{Na^+}` de la soude, le :math:`\\mathrm{Cl^-}` de l'acide chlorhydrique :
ils ne réagissent pas, mais ils conduisent.
"""

from dataclasses import dataclass, replace

import numpy as np
import numpy.typing as npt

from titrages.equilibre import (
    LN10,
    MARGE_PH,
    PKE,
    delta,
    fractions,
    moments,
    racine_croissante,
    somme_h_w,
)
from titrages.modeles import Couple, Espece, Spectateur, charge

__all__ = ["LAMBDA_H3O", "LAMBDA_HO", "Composant", "Solution"]

LAMBDA_H3O = 34.98
"""λ° de H3O+, en mS·m²·mol⁻¹ (25 °C)."""
LAMBDA_HO = 19.86
"""λ° de HO-, en mS·m²·mol⁻¹ (25 °C)."""

Tableau = npt.NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class Composant:
    """Une espèce et sa concentration apportée (mol/L)."""

    espece: Espece
    concentration: float

    def __post_init__(self) -> None:
        if self.concentration < 0:
            raise ValueError(f"{self.espece} : concentration négative ({self.concentration})")


@dataclass(frozen=True, slots=True)
class Solution:
    """Une solution aqueuse idéale, à 25 °C par défaut (:math:`\\mathrm{pK_e} = 14`).

    Immuable : :meth:`avec`, :meth:`diluee` et :meth:`melangee` rendent une nouvelle
    solution.

    >>> from titrages import Couple
    >>> ethanoique = Couple("acide éthanoïque", (4.76,))
    >>> round(Solution().avec(ethanoique.acide, 0.1).pH(), 2)
    2.88
    """

    composants: tuple[Composant, ...] = ()
    acide_fort: float = 0.0
    base_forte: float = 0.0
    pKe: float = PKE
    spectateurs: tuple[Spectateur, ...] = ()

    # --- construction --------------------------------------------------------

    @classmethod
    def eau(cls, pKe: float = PKE) -> "Solution":
        return cls(pKe=pKe)

    def avec(self, espece: Espece | Couple, concentration: float) -> "Solution":
        """Ajoute une espèce, sans dilution (un couple seul désigne sa forme acide)."""
        if isinstance(espece, Couple):
            espece = espece.acide
        return replace(self, composants=(*self.composants, Composant(espece, concentration)))

    def avec_acide_fort(self, concentration: float, contre_ion: Spectateur | None = None) -> "Solution":
        """Un acide fort ; ``contre_ion`` (sa concentration comprise) ne sert qu'à la conductivité."""
        return replace(self, acide_fort=self.acide_fort + concentration).avec_spectateur(contre_ion)

    def avec_base_forte(self, concentration: float, contre_ion: Spectateur | None = None) -> "Solution":
        return replace(self, base_forte=self.base_forte + concentration).avec_spectateur(contre_ion)

    def avec_spectateur(self, ion: Spectateur | None) -> "Solution":
        """Un ion qui ne réagit pas, pour la conductivité."""
        if ion is None:
            return self
        return replace(self, spectateurs=(*self.spectateurs, ion))

    def diluee(self, facteur: float) -> "Solution":
        """Toutes les concentrations divisées par ``facteur``."""
        if facteur <= 0:
            raise ValueError("un facteur de dilution est positif")
        return Solution(
            tuple(Composant(c.espece, c.concentration / facteur) for c in self.composants),
            self.acide_fort / facteur,
            self.base_forte / facteur,
            self.pKe,
            tuple(replace(i, concentration=i.concentration / facteur) for i in self.spectateurs),
        )

    def melangee(self, autre: "Solution", volume: float, volume_autre: float) -> "Solution":
        """Le mélange de ``volume`` de cette solution et de ``volume_autre`` de l'autre."""
        if self.pKe != autre.pKe:
            raise ValueError("deux solutions de pKe différents ne se mélangent pas ici")
        total = volume + volume_autre
        a, b = volume / total, volume_autre / total
        return Solution(
            tuple(Composant(c.espece, a * c.concentration) for c in self.composants)
            + tuple(Composant(c.espece, b * c.concentration) for c in autre.composants),
            a * self.acide_fort + b * autre.acide_fort,
            a * self.base_forte + b * autre.base_forte,
            self.pKe,
            tuple(replace(i, concentration=a * i.concentration) for i in self.spectateurs)
            + tuple(replace(i, concentration=b * i.concentration) for i in autre.spectateurs),
        )

    # --- grandeurs -----------------------------------------------------------

    @property
    def exces_base(self) -> float:
        """:math:`b = [\\text{base forte}] - [\\text{acide fort}]`."""
        return self.base_forte - self.acide_fort

    def charge(self, pH: npt.ArrayLike) -> Tableau:
        """:math:`\\sum_i C_i (Y_{1,i} - q_i) - b` : le résidu sans la part de l'eau.

        C'est la quantité de protons que les espèces ont cédés au-delà de ce qu'on a
        introduit, moins la base forte apportée.
        """
        p = np.asarray(pH, dtype=np.float64)
        total = np.full(p.shape, -self.exces_base, dtype=np.float64)
        for c in self.composants:
            (Y1,) = moments(c.espece.couple, p, ordre=1)
            total = total + c.concentration * (Y1 - c.espece.q)
        return total

    def derivee_charge(self, pH: npt.ArrayLike) -> Tableau:
        """:math:`\\mathrm{d}/\\mathrm{d\\,pH}` de :meth:`charge`, sans la part de l'eau."""
        p = np.asarray(pH, dtype=np.float64)
        total = np.zeros(p.shape, dtype=np.float64)
        for c in self.composants:
            Y1, Y2 = moments(c.espece.couple, p, ordre=2)
            total = total + c.concentration * (Y2 - Y1**2)
        return LN10 * total

    def residu(self, pH: npt.ArrayLike) -> Tableau:
        """Le résidu de l'électroneutralité, nul au pH de la solution et croissant."""
        return self.charge(pH) + delta(pH, self.pKe)

    def pouvoir_tampon(self, pH: npt.ArrayLike | None = None) -> Tableau:
        """:math:`\\beta = \\mathrm{d}C_\\text{base forte}/\\mathrm{d\\,pH}` (mol/L).

        Sans argument, au pH de la solution.
        """
        p = self.pH() if pH is None else pH
        return self.derivee_charge(p) + LN10 * somme_h_w(p, self.pKe)

    def pH(self) -> float:
        """Le pH de la solution, par dichotomie sur le résidu."""
        racine = racine_croissante(self.residu, -MARGE_PH, self.pKe + MARGE_PH)
        return float(racine)

    def concentrations(self, pH: float | None = None) -> dict[str, float]:
        """La concentration de chaque forme de chaque espèce, et de H3O+ et HO-.

        Deux composants du même couple se cumulent : c'est le même équilibre.
        """
        p = self.pH() if pH is None else pH
        resultat: dict[str, float] = {
            "H3O+": 10.0**-p,
            "HO-": 10.0 ** (p - self.pKe),
        }
        for c in self.composants:
            x = fractions(c.espece.couple, p)
            for forme, xj in zip(c.espece.couple.formes, x, strict=True):
                resultat[forme] = resultat.get(forme, 0.0) + c.concentration * float(xj)
        return resultat

    # --- conductivité -------------------------------------------------------------

    def conductivite_ions(self, pH: npt.ArrayLike) -> Tableau:
        """:math:`\\sum \\lambda_i c_i` des espèces et des spectateurs, sans H3O+ ni HO- (S/m).

        C'est la part linéaire en concentrations : celle du bécher et celle de la burette
        se mélangent donc comme les volumes, ce que fait :meth:`Titrage.conductivite`.
        """
        p = np.asarray(pH, dtype=np.float64)
        total = np.zeros(p.shape)
        manquants = []
        for c in self.composants:
            couple = c.espece.couple
            x = fractions(couple, p)
            for j, forme in enumerate(couple.formes):
                if charge(forme) == 0:
                    continue
                lam = couple.conductivites[j] if couple.conductivites else None
                if lam is None:
                    manquants.append(forme)
                    continue
                total = total + lam * c.concentration * x[j]
        for ion in self.spectateurs:
            if ion.conductivite is None:
                manquants.append(ion.formule)
                continue
            total = total + ion.conductivite * ion.concentration
        if manquants:
            raise ValueError(
                "conductivité molaire inconnue pour "
                + ", ".join(dict.fromkeys(manquants))
                + " : la donner par Couple(conductivites=...) ou Spectateur(conductivite=...)"
            )
        return np.asarray(total, dtype=np.float64)

    def conductivite(self, pH: npt.ArrayLike | None = None) -> Tableau:
        """La conductivité :math:`\\sigma` (S/m), au pH de la solution par défaut.

        λ° en mS·m²·mol⁻¹ et :math:`c` en mol/L donnent directement des S/m.

        >>> from titrages import solution
        >>> round(float(solution({"HCl": 0.01}).conductivite()), 4)
        0.4261
        """
        p = self.pH() if pH is None else pH
        h = np.power(10.0, -np.asarray(p, dtype=np.float64))
        w = np.power(10.0, np.asarray(p, dtype=np.float64) - self.pKe)
        return np.asarray(self.conductivite_ions(p) + LAMBDA_H3O * h + LAMBDA_HO * w, dtype=np.float64)
