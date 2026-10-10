"""Les points d'inflexion de référence, calculés par mpmath à 60 chiffres, sans le paquet.

    uv run python tests/reference_mpmath.py

Le script écrit les constantes que ``tests/test_titrage.py`` et ``tests/test_redox.py``
comparent aux sauts du paquet (``INFLEXIONS`` et ``test_inflexion_a_la_precision_de_la_machine``).
Il n'emprunte rien à ``titrages`` : les fractions sont écrites par les formules des manuels,
la courbe directe :math:`V(\\mathrm{pH}) = -V_0 R_\\text{bécher}/R_\\text{burette}` (ou
:math:`V(E)`) l'est en arithmétique multiprécision, et le point d'inflexion est le zéro de sa
dérivée seconde, trouvé par ``mp.findroot`` depuis le minimum de :math:`|V'|` sur une grille.
Pour le dichromate, dimérique, :math:`V(E)` est la racine du bilan d'électrons du mélange,
résolue par ``mp.findroot`` à chaque potentiel. ``mpmath`` n'est qu'une dépendance de
développement : le paquet n'en a pas besoin.
"""

from collections.abc import Callable

from mpmath import mp, mpf

mp.dps = 60

PKE = mpf(14)
S = mpf("0.05916")
V0 = mpf(10)


def _dix(x: mpf) -> mpf:
    return mp.power(10, x)


# --- acide-base ---------------------------------------------------------------------------


def residu(especes: list[tuple[list[float], float, int]], acide: float, base: float) -> Callable[[mpf], mpf]:
    """R(pH) = Σ C (Y1 - q) - (base - acide) + w - h, pour des espèces (pKa, C, forme q)."""

    def R(pH: mpf) -> mpf:
        h = _dix(-pH)
        total = mpf(acide) - mpf(base) + _dix(pH - PKE) - h
        for pKa, C, q in especes:
            k, termes = mpf(1), [mpf(1)]
            for p in sorted(pKa):
                k = k * _dix(-mpf(p))
                termes.append(k / h ** len(termes))
            somme = sum(termes)
            Y1 = sum(j * t for j, t in enumerate(termes)) / somme
            total += mpf(C) * (Y1 - q)
        return total

    return R


def volume_acide_base(R_b: Callable[[mpf], mpf], R_bu: Callable[[mpf], mpf]) -> Callable[[mpf], mpf]:
    """V(pH) = -V0 R_bécher / R_burette."""
    return lambda pH: -V0 * R_b(pH) / R_bu(pH)


def inflexion(V: Callable[[mpf], mpf], a: float, b: float, pas: float = 0.01) -> mpf:
    """Le zéro de V'' entre a et b, depuis le minimum de |V'| sur une grille de pas ``pas``."""
    grille = mp.arange(a, b, pas)
    depart = min(grille, key=lambda x: abs(mp.diff(V, x)))
    return mp.findroot(lambda x: mp.diff(V, x, 2), depart)


ACIDE_BASE = {
    "CH3COOH par NaOH": (([([4.76], 0.1, 0)], 0, 0), ([], 0, 0.1), 7, 11),
    "HCl par NaOH": (([], 0.1, 0), ([], 0, 0.1), 5, 9),
    "NH3 par HCl": (([([9.25], 0.1, 1)], 0, 0), ([], 0.1, 0), 3, 8),
}


# --- oxydoréduction ---------------------------------------------------------------------------


def fraction_oxydee(E: mpf, E0: float, n: int) -> mpf:
    """x_ox = 1 / (1 + 10^(-n (E - E°)/s)), pour un couple monomère à pH 0."""
    return 1 / (1 + _dix(-n * (E - mpf(E0)) / S))


def volume_monomere(
    E0_becher: float, n_b: int, C_b: float, E0_burette: float, n_u: int, C_u: float
) -> Callable[[mpf], mpf]:
    """V(E) = -V0 R_b / R_bu : le bécher est réduit, la burette oxydée."""

    def V(E: mpf) -> mpf:
        R_b = mpf(C_b) * n_b * fraction_oxydee(E, E0_becher, n_b)
        R_bu = mpf(C_u) * n_u * (fraction_oxydee(E, E0_burette, n_u) - 1)
        return -V0 * R_b / R_bu

    return V


def volume_dichromate(C_fe: float, C_cr: float) -> Callable[[mpf], mpf]:
    """Fe2+ par Cr2O7^2- : V(E) est la racine en V du bilan d'électrons [Fe3+] = 3 [Cr3+]."""

    def V(E: mpf) -> mpf:
        K = _dix(6 * (E - mpf("1.33")) / S)  # [Cr2O7^2-] / [Cr3+]^2

        def bilan(v: mpf) -> mpf:
            fe3 = mpf(C_fe) * V0 / (V0 + v) * fraction_oxydee(E, 0.77, 1)
            T = 2 * mpf(C_cr) * v / (V0 + v)  # unités de chrome
            cr3 = (
                2 * T / (1 + mp.sqrt(1 + 8 * K * T))
            )  # la racine positive de 2K u² + u = T, sans soustraction
            return fe3 - 3 * cr3

        return mp.findroot(
            bilan, (mpf("1e-9"), mpf("1e6")), solver="illinois", tol=mpf(10) ** -80, maxsteps=200
        )

    return V


REDOX = {
    "Fe2+ par Ce4+": (volume_monomere(0.77, 1, 0.1, 1.72, 1, 0.1), 0.9, 1.6),
    "Fe2+ par MnO4-": (volume_monomere(0.77, 1, 0.1, 1.51, 5, 0.02), 0.9, 1.5),
    "Fe2+ par Cr2O7 2-": (volume_dichromate(0.1, 0.1 / 6), 0.9, 1.4),
}


def main() -> None:
    print("# acide-base : pH du point d'inflexion (tests/test_titrage.py)")
    for nom, (becher, burette, a, b) in ACIDE_BASE.items():
        pH = inflexion(volume_acide_base(residu(*becher), residu(*burette)), a, b)
        print(f"{nom:22s} {mp.nstr(pH, 13)}")
    print("# rédox : potentiel du point d'inflexion, V/ENH (tests/test_redox.py)")
    for nom, (V, a, b) in REDOX.items():
        E = inflexion(V, a, b, pas=0.005)
        print(f"{nom:22s} {mp.nstr(E, 12)}")


if __name__ == "__main__":
    main()
