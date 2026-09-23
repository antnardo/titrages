import pytest

from titrages import chercher, couple, couples, espece, indicateur, indicateurs, normaliser, solution
from titrages.base import _index


def test_base_sans_ambiguite():
    assert len(_index()) > 200


@pytest.mark.parametrize(
    ("nom", "couple_attendu", "q"),
    [
        ("NH3", "ion ammonium", 1),
        ("ammoniac", "ion ammonium", 1),
        ("ion ammonium", "ion ammonium", 0),
        ("CO3^2-", "acide carbonique", 2),
        ("CO32-", "acide carbonique", 2),
        ("co3 2-", "acide carbonique", 2),
        ("Carbonate de sodium", "acide carbonique", 2),
        ("acide ethanoique", "acide éthanoïque", 0),
        ("acide acétique", "acide éthanoïque", 0),
        ("glycine", "glycine", 1),
        ("Na2H2Y", "EDTA", 2),
    ],
)
def test_recherche(nom, couple_attendu, q):
    e = espece(nom)
    assert (e.couple.nom, e.q) == (couple_attendu, q)


def test_couple_depuis_une_forme():
    assert couple("HPO4^2-").pKa == (2.15, 7.20, 12.35)


def test_inconnu_avec_suggestion():
    with pytest.raises(KeyError, match="peut-être"):
        espece("carbonate")


def test_fort_n_est_pas_une_espece():
    with pytest.raises(KeyError, match="acide fort"):
        espece("HCl")


def test_forts_dans_une_solution():
    s = solution({"soude": 0.1, "acide nitrique": 0.03})
    assert (s.base_forte, s.acide_fort) == (0.1, 0.03)


def test_pKa_ranges_et_formes_coherentes():
    for c in couples():
        assert list(c.pKa) == sorted(c.pKa)
        assert len(c.formes) == c.N + 1


@pytest.mark.parametrize("nom", ["BBT", "bleu de bromothymol", "phenolphtaleine", "méthylorange"])
def test_indicateur(nom):
    i = indicateur(nom)
    assert i.zone[0] < i.pKa < i.zone[1]


def test_indicateurs_zones_ordonnees():
    assert all(a < b for a, b in (i.zone for i in indicateurs()))


def test_indicateur_inconnu():
    with pytest.raises(KeyError):
        indicateur("tournesol mauve")


def test_normaliser():
    assert normaliser("CO_3^{2-}") == normaliser("co3 2-") == "co32-"
    assert normaliser("Éthanoïque") == "ethanoique"


def test_chercher():
    assert "HCO3-" in chercher("hco3")
