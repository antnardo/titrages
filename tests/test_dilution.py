"""La dilution : le calcul direct face à un calcul indirect indépendant.

Le calcul direct est la formule du paquet, V(pH) = -V0 R_bécher(pH) / R_burette(pH),
qui repose sur une seule idée : l'électroneutralité du mélange, multipliée par
V0 + V, se sépare en une part du bécher et une part de la burette, chacune calculée
avec les concentrations d'avant le mélange. Cette correction de dilution est propre à
ce paquet : l'article d'origine (Kalka 2021) trace n0(h) = Y1 + δ/C_T comme une courbe de
titrage, ce qui revient à la négliger.

Le calcul indirect n'emprunte rien au paquet. Pour un volume versé V, il dilue à la
main chaque concentration — C V0/(V0+V) pour le bécher, C V/(V0+V) pour la burette —,
écrit les concentrations de chaque forme par les formules des manuels
(x_j = k_j h^-j / sum_m k_m h^-m, loi de Nernst écrite telle quelle), et cherche le pH
(ou le potentiel) qui annule l'électroneutralité (ou le bilan d'électrons) par une
dichotomie écrite ici. Si la séparation était fausse, les deux calculs divergeraient.
"""

import math

import numpy as np
import pytest

from titrages import Couple, CoupleRedox, Solution, SolutionRedox, Titrage, TitrageRedox

# --- le calcul indirect, acide-base -----------------------------------------------


def _formes(pKa, C, h):
    """[A_j] par la formule des manuels : C k_j h^-j / sum_m k_m h^-m."""
    k = [1.0]
    for p in sorted(pKa):
        k.append(k[-1] * 10.0**-p)
    termes = [kj / h**j for j, kj in enumerate(k)]
    total = sum(termes)
    return [C * t / total for t in termes]


def _electroneutralite(h, especes, Ca, Cb, pKe):
    """Cations moins anions, chaque forme j de charge z0 - j.

    Les contre-ions sont écrits explicitement : une espèce introduite sous sa forme q,
    de charge z0 - q, apporte un contre-ion qui porte la charge opposée.
    """
    w = 10.0**-pKe / h
    bilan = h - w + Cb - Ca
    for pKa, C, q, z0 in especes:
        for j, c in enumerate(_formes(pKa, C, h)):
            bilan += (z0 - j) * c
        bilan -= (z0 - q) * C  # le contre-ion
    return bilan


def ph_force_brute(especes, Ca, Cb, pKe=14.0):
    """Le pH qui annule l'électroneutralité, par dichotomie sur log h."""
    a, b = -pKe - 4, 4.0  # log10 h ; l'électroneutralité décroît quand h diminue
    for _ in range(200):
        m = (a + b) / 2
        if _electroneutralite(10.0**m, especes, Ca, Cb, pKe) > 0:
            b = m
        else:
            a = m
    return -(a + b) / 2


def melange(becher, burette, V0, V):
    """Les concentrations du mélange, diluées à la main."""
    fb, fu = V0 / (V0 + V), V / (V0 + V)
    especes = [(p, C * fb, q, z) for p, C, q, z in becher["especes"]]
    especes += [(p, C * fu, q, z) for p, C, q, z in burette["especes"]]
    Ca = becher.get("Ca", 0) * fb + burette.get("Ca", 0) * fu
    Cb = becher.get("Cb", 0) * fb + burette.get("Cb", 0) * fu
    return especes, Ca, Cb


def construire(desc):
    """La même solution, décrite pour le paquet."""
    s = Solution()
    for pKa, C, q, z0 in desc["especes"]:
        formes = tuple(f"F{j}{'+' * (z0 - j) if z0 > j else '-' * (j - z0)}" for j in range(len(pKa) + 1))
        s = s.avec(Couple("x", tuple(pKa), formes).forme(q), C)
    return s.avec_acide_fort(desc.get("Ca", 0)).avec_base_forte(desc.get("Cb", 0))


# (pKa, C, forme introduite q, charge z0 de la forme la plus acide)
CAS = {
    "acide fort par base forte": ({"Ca": 0.1, "especes": []}, {"Cb": 0.1, "especes": []}),
    "acide faible par base forte": ({"especes": [([4.76], 0.1, 0, 0)]}, {"Cb": 0.1, "especes": []}),
    "acide faible dilué par base concentrée": (
        {"especes": [([4.76], 0.001, 0, 0)]},
        {"Cb": 1.0, "especes": []},
    ),
    "acide faible concentré par base diluée": (
        {"especes": [([4.76], 1.0, 0, 0)]},
        {"Cb": 0.05, "especes": []},
    ),
    "triacide par base forte": ({"especes": [([2.15, 7.2, 12.35], 0.05, 0, 0)]}, {"Cb": 0.1, "especes": []}),
    "base faible par acide fort": ({"especes": [([9.25], 0.1, 1, 1)]}, {"Ca": 0.1, "especes": []}),
    "sel de diacide par acide fort": ({"especes": [([6.35, 10.33], 0.05, 2, 0)]}, {"Ca": 0.1, "especes": []}),
    "mélange par base forte": (
        {"Ca": 0.02, "especes": [([4.76], 0.05, 0, 0), ([9.25], 0.03, 0, 1)]},
        {"Cb": 0.1, "especes": []},
    ),
    "acide faible par base faible": ({"especes": [([4.76], 0.1, 0, 0)]}, {"especes": [([9.25], 0.1, 1, 1)]}),
    "tampon par tampon": (
        {"especes": [([4.76], 0.1, 0, 0), ([4.76], 0.05, 1, 0)]},
        {"Cb": 0.02, "especes": [([6.35, 10.33], 0.1, 1, 0)]},
    ),
}
V0 = 10.0


@pytest.fixture(params=list(CAS))
def cas(request):
    becher, burette = CAS[request.param]
    return becher, burette, Titrage(construire(becher), construire(burette), V0)


def test_independance_du_calcul_indirect():
    # garde-fou : la force brute retrouve les pH des manuels
    assert ph_force_brute([([4.76], 0.1, 0, 0)], 0, 0) == pytest.approx(2.88, abs=0.005)
    assert ph_force_brute([], 0.01, 0) == pytest.approx(2.0, abs=1e-6)
    assert ph_force_brute([([9.25], 0.1, 1, 1)], 0, 0) == pytest.approx(11.12, abs=0.005)


def test_direct_puis_indirect(cas):
    """On choisit des pH, le calcul direct donne V, la force brute redonne le pH."""
    becher, burette, t = cas
    pH0, pH1 = t.pH_initial, float(t.pH(3 * V0))
    for p in np.linspace(pH0, pH1, 25)[1:-1]:
        V = float(t.volume(p))
        assert math.isfinite(V) and V > 0
        assert ph_force_brute(*melange(becher, burette, V0, V)) == pytest.approx(p, abs=1e-9)


@pytest.mark.parametrize("V", [0.0, 0.5, 2.5, 5.0, 9.0, 9.99, 10.0, 10.01, 11.0, 15.0, 20.0, 29.0])
def test_indirect_puis_direct(cas, V):
    """On choisit des volumes, la force brute donne le pH, le calcul direct redonne V."""
    becher, burette, t = cas
    p = ph_force_brute(*melange(becher, burette, V0, V))
    assert float(t.pH(V)) == pytest.approx(p, abs=1e-9)
    if V > 0:
        assert float(t.volume(p)) == pytest.approx(V, rel=1e-7)


def test_le_melange_du_paquet_est_la_dilution_a_la_main(cas):
    """Solution.melangee, qui dilue, rend le même pH que le calcul direct."""
    _, _, t = cas
    for V in (3.0, 10.0, 17.0):
        assert t.becher.melangee(t.burette, V0, V).pH() == pytest.approx(float(t.pH(V)), abs=1e-9)


def test_sans_dilution_l_erreur_de_la_formule_d_origine():
    """n0 = Y1 + δ/C_T, tracé comme courbe de titrage, se trompe loin de l'équivalence.

    Pour l'acide éthanoïque 0,1 M par la soude 0,1 M, à deux fois l'équivalence, la
    formule sans dilution donne un pH trop élevé de log(3) ≈ 0,48 : elle compte la
    base versée sans l'eau qui l'accompagne.
    """
    becher, burette = CAS["acide faible par base forte"]
    avec = Titrage(construire(becher), construire(burette), V0)
    sans = Titrage(construire(becher), construire(burette), V0, dilution=False)
    V = 20.0
    exact = ph_force_brute(*melange(becher, burette, V0, V))
    assert float(avec.pH(V)) == pytest.approx(exact, abs=1e-9)
    assert float(sans.pH(V)) - exact == pytest.approx(math.log10(3), abs=0.01)
    # à la demi-équivalence, la dilution ne change presque rien
    assert float(sans.pH(5.0)) == pytest.approx(float(avec.pH(5.0)), abs=0.01)


def test_conductivite_dilution_a_la_main():
    """σ(V) du paquet, face à Kohlrausch sur les concentrations diluées à la main."""
    from titrages import solution

    t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), V0)
    for V in (2.0, 10.0, 15.0):
        C_ac, Cb = 0.1 * V0 / (V0 + V), 0.1 * V / (V0 + V)
        p = ph_force_brute([([4.76], C_ac, 0, 0)], 0, Cb)
        h, w = 10**-p, 10 ** (p - 14)
        acetate = _formes([4.76], C_ac, h)[1]
        sigma = 34.98 * h + 19.86 * w + 4.09 * acetate + 5.01 * Cb
        assert float(t.conductivite(V)) == pytest.approx(sigma, rel=1e-9)


# --- le calcul indirect, oxydoréduction -------------------------------------------

S = 0.05916


def _bilan_electrons(E, couples):
    """Électrons portés par les formes oxydées, moins ceux introduits, par force brute.

    Monomère : [Ox]/[Red] = 10^(n (E - E0)/s). Dimère (1, 2) : [Ox]/[Red]² = K, la
    concentration de Red tirée de T = [Red] + 2 K [Red]² par dichotomie, sans formule.
    """
    total = 0.0
    for E0, n, st, T, introduits in couples:
        K = 10 ** (n * (E - E0) / S)
        if st == (1, 1):
            oxyde = T * K / (1 + K)
            total += n * oxyde - introduits
        else:
            a, b = 0.0, T
            for _ in range(200):
                r = (a + b) / 2
                if r + 2 * K * r * r > T:
                    b = r
                else:
                    a = r
            total += n / 2 * (T - r) - introduits  # Ox dimère : n électrons pour deux unités
    return total


def e_force_brute(couples):
    a, b = -1.0, 3.0
    for _ in range(200):
        m = (a + b) / 2
        if _bilan_electrons(m, couples) > 0:
            b = m
        else:
            a = m
    return (a + b) / 2


FER = CoupleRedox("fer", (0.77,), (1,), ("Fe2+", "Fe3+"))
CERIUM = CoupleRedox("cérium", (1.72,), (1,), ("Ce3+", "Ce4+"))
DICHROMATE = CoupleRedox("dichromate", (1.33,), (6,), ("Cr3+", "Cr2O7^2-"), stoechiometrie=(1, 2))


@pytest.mark.parametrize("V", [1.0, 5.0, 9.5, 10.5, 15.0, 19.0])
def test_redox_monomeres(V):
    t = TitrageRedox(SolutionRedox().avec(FER.red, 0.1), SolutionRedox().avec(CERIUM.ox, 0.1), V0)
    fb, fu = V0 / (V0 + V), V / (V0 + V)
    # Fe2+ introduit : aucun électron cédé ; Ce4+ introduit : 1 électron cédé par unité
    E = e_force_brute([(0.77, 1, (1, 1), 0.1 * fb, 0.0), (1.72, 1, (1, 1), 0.1 * fu, 0.1 * fu)])
    assert float(t.E(V)) == pytest.approx(E, abs=1e-9)
    assert float(t.volume(E)) == pytest.approx(V, rel=1e-7)


@pytest.mark.parametrize("V", [1.0, 5.0, 9.5, 10.5, 15.0, 19.0])
def test_redox_un_dimere_trinome(V):
    """La racine du trinôme (calcul direct) face au bilan d'électrons résolu pas à pas."""
    C = 1 / 60
    t = TitrageRedox(SolutionRedox().avec(FER.red, 0.1), SolutionRedox().avec(DICHROMATE.ox, C), V0)
    fb, fu = V0 / (V0 + V), V / (V0 + V)
    T = 2 * C * fu  # unités de chrome
    E = e_force_brute([(0.77, 1, (1, 1), 0.1 * fb, 0.0), (1.33, 6, (1, 2), T, 3 * T)])
    assert float(t.E(V)) == pytest.approx(E, abs=1e-9)
    assert float(t.volume(E)) == pytest.approx(V, rel=1e-7)
