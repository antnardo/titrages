"""Le tracé avec matplotlib, dépendance facultative (``pip install matplotlib``, ou l'extra ``[trace]``).

Chaque fonction dessine sur un ``Axes`` qu'on lui donne, ou sur l'axe courant, et le
rend : elle ne crée ni ne montre de figure, pour se composer avec le reste d'un
script. Les réglages de style passent par ``**style`` jusqu'à ``Axes.plot``.
"""

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from titrages.equilibre import fractions
from titrages.modeles import Couple, Indicateur, Saut, Tangentes
from titrages.redox import TitrageRedox
from titrages.titrage import Titrage

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "tracer_conductimetrie",
    "tracer_derivee",
    "tracer_fractions",
    "tracer_sauts",
    "tracer_tangentes",
    "tracer_titrage",
    "tracer_titrage_redox",
    "tracer_zone_virage",
]


def _axe(ax: "Axes | None") -> "Axes":
    if ax is not None:
        return ax
    try:
        import matplotlib.pyplot as plt
    except ImportError as erreur:  # pragma: no cover — dépend de l'environnement
        raise ImportError("le tracé demande matplotlib : pip install matplotlib") from erreur
    return plt.gca()


def tracer_titrage(
    titrage: Titrage,
    volume_max: float,
    ax: "Axes | None" = None,
    *,
    points: int = 400,
    abscisse: float = 1.0,
    etiquettes: bool = True,
    **style: Any,
) -> "Axes":
    """La courbe pH = f(V), de 0 à ``volume_max``.

    ``abscisse`` divise les volumes : ``abscisse=t.volumes_equivalents()[0]`` donne
    l'abscisse normalisée :math:`n = V/V_\\text{éq}` de la fiche.
    """
    ax = _axe(ax)
    c = titrage.courbe(volume_max, points)
    ax.plot(c.volume / abscisse, c.pH, **style)
    if etiquettes:
        ax.set_xlabel("$V$ versé" if abscisse == 1.0 else r"$V/V_\mathrm{éq}$")
        ax.set_ylabel("pH")
    return ax


def tracer_derivee(
    titrage: Titrage | TitrageRedox,
    volume_max: float,
    ax: "Axes | None" = None,
    *,
    points: int = 400,
    **style: Any,
) -> "Axes":
    """La pente de la courbe selon :math:`V` : d pH/dV, ou dE/dV pour un titrage rédox."""
    ax = _axe(ax)
    c = titrage.courbe(volume_max, points)
    ax.plot(c.volume, c.pente, **style)
    return ax


def tracer_sauts(
    titrage: Titrage | TitrageRedox,
    volume_max: float,
    ax: "Axes | None" = None,
    *,
    texte: bool = True,
    **style: Any,
) -> "Axes":
    """Les points d'inflexion, marqués et reportés sur les axes, en pH ou en E."""
    ax = _axe(ax)
    style = {"marker": "o", "ls": "", "color": "k", "ms": 4, **style}
    for s in titrage.sauts(volume_max):
        y = s.pH if isinstance(s, Saut) else s.E
        ax.plot([s.volume], [y], **style)
        ax.plot([s.volume, s.volume, 0], [0, y, y], ls=":", lw=0.7, color=style["color"])
        if texte:
            ax.annotate(
                f"{s.volume:.3g} ; {y:.2f}",
                (s.volume, y),
                xytext=(4, -10),
                textcoords="offset points",
                fontsize="small",
            )
    return ax


def tracer_tangentes(
    titrage: Titrage | TitrageRedox,
    volume_avant: float,
    volume_max: float,
    ax: "Axes | None" = None,
    *,
    marge: float = 0.1,
    **style: Any,
) -> "Axes":
    """La construction des tangentes parallèles, et le point obtenu, sur pH(V) ou E(V)."""
    ax = _axe(ax)
    t = titrage.tangentes(volume_avant, volume_max)
    ordonnees = (t.pH_avant, t.pH_apres, t.pH) if isinstance(t, Tangentes) else (t.E_avant, t.E_apres, t.E)
    style = {"lw": 0.6, "ls": "--", "color": "C3", **style}
    ecart = t.volume_apres - t.volume_avant
    V = np.array([t.volume_avant - marge * ecart, t.volume_apres + marge * ecart])
    for laquelle in ("avant", "apres", "milieu"):
        ax.plot(V, t.droite(V, laquelle), **style)
    ax.plot([t.volume_avant, t.volume_apres, t.volume], ordonnees, marker="x", ls="", color=style["color"])
    return ax


def tracer_fractions(
    couple: Couple | Sequence[float],
    ax: "Axes | None" = None,
    *,
    pH: tuple[float, float] = (0.0, 14.0),
    points: int = 400,
    legende: bool = True,
    **style: Any,
) -> "Axes":
    """Le diagramme de distribution : les fractions de chaque forme selon le pH."""
    ax = _axe(ax)
    p = np.linspace(*pH, points)
    x = fractions(couple, p)
    formes = couple.formes if isinstance(couple, Couple) else Couple.polyacide(*couple).formes
    for forme, xj in zip(formes, x, strict=True):
        ax.plot(p, xj, label=forme, **style)
    ax.set_xlim(*pH)
    ax.set_ylim(0, 1)
    ax.set_xlabel("pH")
    ax.set_ylabel("fraction")
    if legende:
        ax.legend()
    return ax


def tracer_zone_virage(indicateur: Indicateur, ax: "Axes | None" = None, *, alpha: float = 0.25) -> "Axes":
    """La zone de virage d'un indicateur, en bande horizontale sur une courbe pH = f(V)."""
    ax = _axe(ax)
    a, b = indicateur.zone
    ax.axhspan(a, b, color=indicateur.teinte_base, alpha=alpha, lw=0, label=indicateur.nom)
    return ax


def tracer_conductimetrie(
    titrage: Titrage,
    volume_max: float,
    ax: "Axes | None" = None,
    *,
    points: int = 200,
    corrigee: bool = False,
    unite: str = "S/m",
    etiquettes: bool = True,
    **style: Any,
) -> "Axes":
    """La courbe :math:`\\sigma = f(V)`, en S/m, l'unité du paquet.

    ``unite="mS/cm"`` affiche ce que lit un conductimètre (dix fois moins).
    """
    ax = _axe(ax)
    c = titrage.courbe_conductimetrique(volume_max, points, corrigee)
    facteur = {"S/m": 1.0, "mS/cm": 10.0, "mS/m": 1000.0}[unite]
    ax.plot(c.volume, facteur * c.conductivite, **style)
    if etiquettes:
        ax.set_xlabel("$V$ versé")
        nom = r"$\sigma\,(V_0+V)/V_0$" if corrigee else r"$\sigma$"
        ax.set_ylabel(f"{nom} ({unite})")
    return ax


def tracer_titrage_redox(
    titrage: TitrageRedox,
    volume_max: float,
    ax: "Axes | None" = None,
    *,
    points: int = 400,
    sauts: bool = False,
    etiquettes: bool = True,
    **style: Any,
) -> "Axes":
    """La courbe :math:`E = f(V)`, par rapport à l'électrode de référence du titrage."""
    ax = _axe(ax)
    c = titrage.courbe(volume_max, points)
    ax.plot(c.volume, c.E, **style)
    if sauts:
        for s in titrage.sauts(volume_max):
            ax.plot([s.volume], [s.E], "ko", ms=3)
    if etiquettes:
        ref = titrage.reference if isinstance(titrage.reference, str) else f"{titrage.reference:g} V"
        ax.set_xlabel("$V$ versé")
        ax.set_ylabel(f"$E$ (V/{ref})")
    return ax
