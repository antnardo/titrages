"""Le calcul d'équilibre : fractions, moments, résidu de l'électroneutralité.

Tout part d'une remarque (Kalka 2021, Brière 2002) : **à pH imposé, tout est
explicite**. Les fractions molaires s'écrivent directement en fonction de
:math:`h`, et l'électroneutralité se lit alors comme une équation *linéaire* en la
quantité de base ajoutée. Chercher le pH d'une composition donnée, c'est au contraire
chercher la racine d'un polynôme de degré :math:`N+2` — sans formule dès
:math:`N>2`. On calcule donc dans le sens direct, et l'on n'inverse que là où c'est
indispensable, par dichotomie sur une fonction qu'on sait monotone.

**Les fractions se calculent en échelle logarithmique.** La forme naïve
:math:`x_j = x_0 k_j/h^j` déborde pour un polyacide à pH élevé (l'EDTA à pH 14 fait
intervenir :math:`10^{56}` puis :math:`10^{-56}`) : on retranche le plus grand
exposant avant d'élever à la puissance, comme dans un *log-sum-exp*, et aucune
fraction n'est jamais calculée comme le quotient de deux nombres démesurés.

Notations de la fiche : :math:`\\delta = w - h`, :math:`Y_k = \\sum_j j^k x_j`.
"""

import math
from collections.abc import Callable, Sequence

import numpy as np
import numpy.typing as npt

from titrages.modeles import Couple, Espece

__all__ = [
    "LN10",
    "PKE",
    "concentration_pour_pH",
    "delta",
    "ecart_moyen",
    "fractions",
    "log_k",
    "minimum_lisse",
    "moments",
    "racine_croissante",
    "somme_h_w",
    "variance",
]

LN10 = math.log(10)
PKE = 14.0

Tableau = npt.NDArray[np.float64]
Constantes = Couple | Sequence[float]

# Bornes de la dichotomie. Hors de [-6, pKe + 6], l'eau seule (|δ| > 10⁶ mol/L)
# l'emporte sur toute concentration réaliste : le résidu y a son signe certain.
MARGE_PH = 6.0


def _pKa(constantes: Constantes) -> tuple[float, ...]:
    if isinstance(constantes, Couple):
        return constantes.pKa
    return tuple(sorted(float(p) for p in constantes))


def log_k(constantes: Constantes) -> Tableau:
    """Les :math:`\\log_{10} k_j` des dissociations cumulées, :math:`k_0 = 1`.

    :math:`k_j = K_1 \\cdots K_j` est la constante de
    :math:`\\mathrm{H}_N\\mathrm{A} = j\\,\\mathrm{H^+} + \\mathrm{H}_{N-j}\\mathrm{A}^{j-}`.
    """
    pKa = _pKa(constantes)
    return np.concatenate(([0.0], -np.cumsum(pKa))).astype(np.float64)


def fractions(constantes: Constantes, pH: npt.ArrayLike) -> Tableau:
    """Les fractions molaires :math:`x_j`, de forme ``(N + 1, *pH.shape)``.

    >>> x = fractions([6.35, 10.33], [6.35, 8.34, 10.33])
    >>> [round(float(v), 3) for v in x[:, 0]]
    [0.5, 0.5, 0.0]
    """
    lk = log_k(constantes)
    pH_ = np.asarray(pH, dtype=np.float64)
    forme = (-1,) + (1,) * pH_.ndim
    j = np.arange(len(lk), dtype=np.float64).reshape(forme)
    t = lk.reshape(forme) + j * pH_
    t = t - t.max(axis=0)
    e = np.power(10.0, t)
    return np.asarray(e / e.sum(axis=0), dtype=np.float64)


def moments(constantes: Constantes, pH: npt.ArrayLike, ordre: int = 2) -> tuple[Tableau, ...]:
    """Les moments :math:`Y_1, \\dots, Y_{\\text{ordre}}` de la distribution des formes.

    :math:`Y_1` est le nombre moyen de protons cédés,
    :math:`Y_2 - Y_1^2` sa variance : c'est elle qui fait le pouvoir tampon.
    """
    x = fractions(constantes, pH)
    forme = (-1,) + (1,) * (x.ndim - 1)
    j = np.arange(x.shape[0], dtype=np.float64).reshape(forme)
    return tuple(np.asarray((j**k * x).sum(axis=0), dtype=np.float64) for k in range(1, ordre + 1))


def variance(constantes: Constantes, pH: npt.ArrayLike) -> Tableau:
    """:math:`Y_2 - Y_1^2`, écrite :math:`\\sum_{j<k} x_jx_k(k-j)^2` pour ne rien perdre.

    Quand une forme domine à :math:`1 - 10^{-10}`, :math:`Y_2` et :math:`Y_1^2` coïncident
    sur dix chiffres et leur différence n'en garde que six ; c'est justement le cas
    dans un saut, où l'on cherche le point d'inflexion.
    """
    x = fractions(constantes, pH)
    total = np.zeros(x.shape[1:], dtype=np.float64)
    for j in range(x.shape[0]):
        for k in range(j + 1, x.shape[0]):
            total = total + x[j] * x[k] * (k - j) ** 2
    return total


def ecart_moyen(constantes: Constantes, pH: npt.ArrayLike, q: float) -> Tableau:
    """:math:`Y_1 - q = \\sum_j x_j (j - q)`, sans soustraire deux nombres voisins."""
    x = fractions(constantes, pH)
    return np.asarray(sum(x[j] * (j - q) for j in range(x.shape[0])), dtype=np.float64)


def delta(pH: npt.ArrayLike, pKe: float = PKE) -> Tableau:
    """:math:`\\delta = w - h`, l'écart à la neutralité de l'eau."""
    p = np.asarray(pH, dtype=np.float64)
    return np.asarray(np.power(10.0, p - pKe) - np.power(10.0, -p), dtype=np.float64)


def somme_h_w(pH: npt.ArrayLike, pKe: float = PKE) -> Tableau:
    """:math:`w + h`, soit :math:`\\mathrm{d}\\delta/\\mathrm{d\\,pH}` divisé par ln 10."""
    p = np.asarray(pH, dtype=np.float64)
    return np.asarray(np.power(10.0, p - pKe) + np.power(10.0, -p), dtype=np.float64)


def concentration_pour_pH(espece: Espece, pH: npt.ArrayLike, pKe: float = PKE) -> Tableau:
    """La concentration de ``espece`` seule dans l'eau qui donne ce pH — le sens direct.

    L'électroneutralité :math:`C(Y_1 - q) + \\delta = 0` donne
    :math:`C = -\\delta/(Y_1 - q)`. Là où ce n'est pas positif, aucune concentration
    ne donne ce pH (un acide ne rend pas une solution basique) : ``nan``.
    """
    (Y1,) = moments(espece.couple, pH, ordre=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        C = -delta(pH, pKe) / (Y1 - espece.q)
    return np.asarray(np.where(C > 0, C, np.nan), dtype=np.float64)


def racine_croissante(
    f: Callable[[Tableau], Tableau],
    bas: npt.ArrayLike,
    haut: npt.ArrayLike,
    iterations: int = 64,
) -> Tableau:
    """La racine d'une fonction croissante sur ``[bas, haut]``, par dichotomie vectorisée.

    Une dichotomie plutôt qu'un Newton : les fonctions en jeu sont monotones mais
    presque plates dans les zones tampon et presque verticales aux sauts, là où un
    Newton diverge. Soixante-quatre coupes amènent l'intervalle à la précision de la
    machine, pour chaque composante à la fois, sans cas particulier.
    """
    a = np.array(bas, dtype=np.float64)
    b = np.array(haut, dtype=np.float64)
    a, b = np.broadcast_arrays(a, b)
    a, b = a.copy(), b.copy()
    for _ in range(iterations):
        m = (a + b) / 2
        positif = f(m) > 0
        b = np.where(positif, m, b)
        a = np.where(positif, a, m)
    return np.asarray((a + b) / 2, dtype=np.float64)


def minimum_lisse(f: Callable[[float], float], a: float, b: float, pas: float, iterations: int = 64) -> float:
    """Le minimum d'une fonction positive et lisse, à la précision de la machine.

    Une section dorée compare des valeurs de :math:`f` ; près du minimum elles ne
    diffèrent qu'au second ordre, et sa précision plafonne à :math:`\\sqrt\\varepsilon`
    — une dizaine de microvolts sur un saut de potentiel. On cherche plutôt le zéro de
    :math:`(\\ln f)'`, estimée par différence centrée de pas ``pas`` : une dérivée
    s'annule au premier ordre, et la dichotomie sur son signe atteint la précision de
    la machine. ``[a, b]`` doit encadrer un seul minimum.
    """

    def pente(x: float) -> float:
        return math.log(f(x + pas)) - math.log(f(x - pas))

    a, b = sorted((float(a), float(b)))
    if not (pente(a) < 0 < pente(b)):  # pas d'encadrement net : on se contente de a ou b
        return float(a if f(a) <= f(b) else b)
    for _ in range(iterations):
        m = (a + b) / 2
        if pente(m) > 0:
            b = m
        else:
            a = m
    return float((a + b) / 2)
