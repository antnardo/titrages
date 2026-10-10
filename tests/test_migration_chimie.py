"""L'ancien ``chimie/titrage.py`` de pychimie est remplacé par ``titrages`` : rien ne s'est perdu.

Les valeurs attendues ont été calculées le 8 octobre 2026 avec l'ancien module, avant
sa suppression. Chaque test refait le même calcul avec ``titrages`` et indique l'API
qui remplace l'ancienne. Là où l'ancien code interpolait sur une grille, l'accord est
de l'ordre de la finesse de la grille ; là où il se trompait (``beta_dilue``), le test
vérifie que ``titrages`` donne la valeur juste, recalculée par différence finie. Les
fixtures sont propres à ce fichier : celles de ``conftest.py`` (acétique, phosphorique,
ammoniac) ne sont pas les solutions de l'ancien module.
"""

import numpy as np
import pytest

from titrages import (
    Couple,
    Solution,
    Titrage,
    concentration_pour_pH,
    fractions,
    solution,
)

PH = np.array([5.0, 7.0, 8.3, 9.0, 11.0])
CARBONATE = Couple.polyacide(6.3, 10.3)
DIACIDE = Couple.polyacide(4.0, 10.0)


@pytest.fixture
def carbonate() -> Solution:
    return Solution().avec(CARBONATE, 0.1)


@pytest.fixture
def titrage_carbonate(carbonate) -> Titrage:
    return Titrage(carbonate, Solution().avec_base_forte(0.1), volume_becher=10)


def test_fractions_reproduisent_polyacide_fractions():
    attendu = [
        [
            0.9522730511818606,
            0.16626806053145238,
            0.009803921568627432,
            0.0018964317004907742,
            3.318859053170638e-06,
        ],
        [
            0.047726709617963874,
            0.8333142929833052,
            0.9803921568627452,
            0.9504673577320353,
            0.1663369787657421,
        ],
        [
            2.39200175635594e-07,
            0.0004176464852423683,
            0.00980392156862745,
            0.047636210567474094,
            0.8336597023752047,
        ],
    ]
    assert fractions(CARBONATE, PH) == pytest.approx(np.array(attendu), rel=1e-12)


def test_residu_sur_C_reproduit_n_sans_dilution(carbonate):
    attendu = [
        0.047627198018315144,
        0.83414958595379,
        1.0000199025044263,
        1.0458397688669836,
        1.8436563834161515,
    ]
    assert carbonate.residu(PH) / 0.1 == pytest.approx(attendu, rel=1e-12)


def test_volume_reproduit_n_dilue(titrage_carbonate):
    attendu = [
        0.04762243625091441,
        0.83414958595379,
        1.0000398058010878,
        1.0459443528428245,
        1.8622791749797207,
    ]
    assert titrage_carbonate.volume(PH) / 10 == pytest.approx(attendu, rel=1e-12)


def test_pouvoir_tampon_reproduit_beta(carbonate):
    attendu = [
        0.10488227191585174,
        0.3204768952334999,
        0.045194785328864985,
        0.10946610394132199,
        0.3423484795511438,
    ]
    assert carbonate.pouvoir_tampon(PH) / 0.1 == pytest.approx(attendu, rel=1e-12)


def test_derivee_volume_corrige_beta_dilue(titrage_carbonate):
    """``beta_dilue`` employait ``n`` au lieu de ``n_dilue`` : 1,2 % d'erreur à pH 11."""
    ancien_faux = 0.38868710806804657
    h = 1e-6
    juste = (titrage_carbonate.volume(11 + h) - titrage_carbonate.volume(11 - h)) / (2 * h) / 10
    assert titrage_carbonate.derivee_volume(11.0) / 10 == pytest.approx(juste, rel=1e-6)
    assert abs(juste / ancien_faux - 1) > 1e-3


def test_concentration_pour_pH_reproduit_Ctot():
    attendu = [0.020052583046987066, 0.00020950322898057432, 2.9650620991011235e-06]
    assert concentration_pour_pH(CARBONATE.acide, [4.0, 5.0, 6.0]) == pytest.approx(attendu, rel=1e-12)


@pytest.mark.parametrize(
    ("C", "pH_attendu"),
    [(1e-6, 6.289470906927095), (1e-3, 4.654855756208415), (0.1, 3.6504859912094636)],
)
def test_pH_reproduit_pH_inverse(C, pH_attendu):
    assert Solution().avec(CARBONATE, C).pH() == pytest.approx(pH_attendu, abs=1e-9)


def test_pH_du_titrage_reproduit_pH_EP():
    """L'ancien ``pH_EP`` interpolait sur 200 points : accord à 1e-4 près."""
    t = Titrage(Solution().avec(DIACIDE, 0.1), Solution().avec_base_forte(0.1), 10)
    attendu = [4.002586669525425, 7.000000007360507, 9.995696214870378, 11.24949645910876]
    assert t.pH([5, 10, 15, 20]) == pytest.approx(attendu, abs=1e-4)


def test_tangentes_reproduisent_methode_tangentes():
    t = Titrage(Solution().avec(DIACIDE, 0.1), Solution().avec_base_forte(0.1), 10)
    tg = t.tangentes(6, 15)
    assert (tg.volume / 10, tg.pH) == pytest.approx((0.999997878871766, 6.999539672051506), abs=1e-4)


def test_pH_avec_acide_et_base_forts_reproduit_pH_melange():
    s = Solution().avec(CARBONATE, 0.1).avec_acide_fort(0.02).avec_base_forte(0.05)
    assert s.pH() == pytest.approx(5.9315189152510595, abs=1e-3)


def test_solution_a_deux_couples_reproduit_melange():
    s = Solution().avec(CARBONATE, 0.1).avec(Couple.polyacide(2.1, 7.2, 12.7), 0.1)
    attendu = [
        0.5263199538078107,
        1.110501652971882,
        1.4632359253102323,
        1.51521877060999,
        1.9315301416173032,
    ]
    assert s.residu(PH) / 0.2 == pytest.approx(attendu, rel=1e-12)


def test_acide_fort_exact_reproduit_acidefort_approche():
    """L'ancien ``AcideFort()`` était un pKa de -5 : accord à 1e-6 près."""
    t = Titrage(solution({"HCl": 0.1}), solution({"NaOH": 0.4}), 10)
    attendu = [
        0.6336754729606564,
        0.8780486829387368,
        0.9875311622191405,
        1.0125313281951747,
        1.1282051281919792,
    ]
    assert t.volume([1.5, 2.0, 3.0, 11.0, 12.0]) * 0.4 / 1.0 == pytest.approx(attendu, rel=1e-6)


def test_titrage_sans_dilution_reproduit_le_script_minimal(carbonate):
    t = Titrage(carbonate, Solution().avec_base_forte(0.1), 10, dilution=False)
    assert t.pH([5, 10]) == pytest.approx([6.2998704946847734, 8.29955170090878], abs=1e-4)
