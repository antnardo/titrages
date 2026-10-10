import numpy as np
import pytest

from titrages import Couple, Espece, Tangentes


@pytest.mark.parametrize(
    ("formule", "attendu"),
    [
        ("H3N+CH2COOH", 1),  # ion glycinium : la charge est à l'intérieur
        ("H3N+CH2COO-", 0),  # la glycine, zwitterion
        ("H2NCH2COO-", -1),
        ("(C2H5)3NH+", 1),
        ("B(OH)4-", -1),
        ("CH2(COO)2^2-", -2),
        ("Fe(CN)6^3-", -3),
        ("VO2^+", 1),
    ],
)
def test_charge_zwitterions_et_parentheses(formule, attendu):
    from titrages.modeles import charge

    assert charge(formule) == attendu


def test_couple_range_les_pKa():
    assert Couple("x", (10.33, 6.35)).pKa == (6.35, 10.33)


def test_couple_formes_generiques():
    assert Couple.polyacide(2, 7, 12).formes == ("H3A", "H2A-", "HA2-", "A3-")


@pytest.mark.parametrize(
    ("pKa", "formes"),
    [((), ()), ((4.0,), ("HA",)), ((4.0,), ("HA", "A-", "A2-"))],
)
def test_couple_incoherent_refuse(pKa, formes):
    with pytest.raises(ValueError):
        Couple("x", pKa, formes)


def test_espece_hors_des_formes_refusee():
    with pytest.raises(ValueError):
        Espece(Couple.polyacide(4.0), 2)


def test_acide_et_base():
    c = Couple.polyacide(2, 7)
    assert (c.acide.q, c.base.q, c.N) == (0, 2, 2)
    assert str(c.forme(1)) == "HA-"


def test_tangentes_droites():
    t = Tangentes(1.0, 4.0, 3.0, 10.0, 0.5, 2.0, 7.0)
    assert t.droite([1.0], "avant")[0] == pytest.approx(4)
    assert t.droite(np.array([2.0]), "milieu")[0] == pytest.approx(7)
