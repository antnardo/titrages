import numpy as np
import pytest

from titrages import Couple, Solution, solution


@pytest.mark.parametrize(
    ("composition", "attendu"),
    [
        ({}, 7.00),
        ({"HCl": 0.01}, 2.00),
        ({"NaOH": 0.01}, 12.00),
        ({"HCl": 1e-8}, 6.98),  # l'eau l'emporte : pH < 7 mais pas 8
        ({"CH3COOH": 0.1}, 2.88),
        ({"CH3COOH": 0.1, "CH3COO-": 0.1}, 4.76),
        ({"NH3": 0.1}, 11.12),
        ({"NH4Cl": 0.1}, 5.13),
        ({"NaHCO3": 0.1}, 8.34),
        ({"Na2CO3": 0.1}, 11.66),
        ({"H3PO4": 0.1}, 1.63),
    ],
)
def test_pH_de_reference(composition, attendu):
    assert solution(composition).pH() == pytest.approx(attendu, abs=0.01)


def test_pKe_autre_temperature():
    assert Solution(pKe=13.0).pH() == pytest.approx(6.5)


def test_residu_croissant():
    s = solution({"H3PO4": 0.05, "NH3": 0.02, "HCl": 0.01})
    assert np.all(np.diff(s.residu(np.linspace(-2, 16, 500))) > 0)


def test_pouvoir_tampon_est_la_derivee_du_residu():
    s = solution({"CH3COOH": 0.1, "CH3COO-": 0.05})
    p, h = 4.5, 1e-6
    fd = (s.residu(p + h) - s.residu(p - h)) / (2 * h)
    assert s.pouvoir_tampon(p) == pytest.approx(fd, rel=1e-6)


def test_pouvoir_tampon_maximal_au_pKa():
    s = solution({"CH3COOH": 0.1, "CH3COO-": 0.1})
    pH = np.linspace(3, 6.5, 351)
    assert pH[np.argmax(s.pouvoir_tampon(pH))] == pytest.approx(4.76, abs=0.02)


def test_diluee_et_melangee():
    s = solution({"CH3COOH": 0.1, "HCl": 0.02})
    assert s.diluee(10).composants[0].concentration == pytest.approx(0.01)
    m = s.melangee(solution({"NaOH": 0.1}), 10, 10)
    assert m.acide_fort == pytest.approx(0.01)
    assert m.base_forte == pytest.approx(0.05)


def test_concentrations_bilan_de_matiere():
    s = solution({"H3PO4": 0.1})
    c = s.concentrations()
    total = sum(c[f] for f in ("H3PO4", "H2PO4-", "HPO4^2-", "PO4^3-"))
    assert total == pytest.approx(0.1)
    assert c["H3O+"] == pytest.approx(10 ** -s.pH())


def test_avec_couple_designe_la_forme_acide():
    c = Couple.polyacide(4.76)
    assert Solution().avec(c, 0.1).composants[0].espece.q == 0


def test_concentration_negative_refusee():
    with pytest.raises(ValueError):
        solution({"NH3": -1})
