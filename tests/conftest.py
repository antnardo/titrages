import matplotlib
import pytest

from titrages import Titrage, solution

matplotlib.use("Agg")


@pytest.fixture
def acetique():
    """Le titrage de la fiche : 10 mL d'acide éthanoïque 0,1 M par la soude 0,1 M."""
    return Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)


@pytest.fixture
def phosphorique():
    return Titrage(solution({"H3PO4": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)


@pytest.fixture
def ammoniac():
    return Titrage(solution({"NH3": 0.1}), solution({"HCl": 0.1}), volume_becher=10)
