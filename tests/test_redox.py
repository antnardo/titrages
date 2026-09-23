import numpy as np
import pytest

from titrages import (
    CoupleRedox,
    SolutionRedox,
    TitrageRedox,
    couple_redox,
    electrons_cedes,
    espece_redox,
    fractions_redox,
    solution_redox,
)

FER = CoupleRedox("fer", (0.77,), (1,), ("Fe2+", "Fe3+"))
CERIUM = CoupleRedox("cérium", (1.72,), (1,), ("Ce3+", "Ce4+"))
VOLUMES = np.array([0.5, 2, 5, 9, 9.9, 10, 10.1, 11, 15, 19])


@pytest.fixture
def fer_cerium():
    return TitrageRedox(SolutionRedox().avec(FER.red, 0.1), SolutionRedox().avec(CERIUM.ox, 0.1), 10)


def test_fractions_egales_au_potentiel_standard():
    x = fractions_redox(FER, 0.77)
    assert x[0] == pytest.approx(0.5)


def test_fractions_potentiel_apparent_selon_le_pH():
    mn = CoupleRedox("Mn", (1.51,), (5,), ("Mn2+", "MnO4-"), (8,))
    E_app = 1.51 - 0.05916 * 8 / 5 * 2
    assert fractions_redox(mn, E_app, pH=2)[0] == pytest.approx(0.5)


def test_fractions_sans_debordement():
    v = couple_redox("vanadium")
    x = fractions_redox(v, np.array([-3.0, 4.0]))
    assert np.all(np.isfinite(x)) and x[0, 0] == pytest.approx(1) and x[-1, 1] == pytest.approx(1)


def test_dimere_fraction_depend_de_la_concentration():
    iode = couple_redox("I2")
    Y1 = electrons_cedes(iode, 0.54, T=0.1)
    Y2 = electrons_cedes(iode, 0.54, T=0.001)
    assert Y1 != pytest.approx(Y2)


@pytest.mark.parametrize(
    ("composition", "attendu"),
    [
        ({"Fe2+": 0.1, "Fe3+": 0.1}, 0.77),
        ({"Fe2+": 0.1, "Fe3+": 0.01}, 0.77 - 0.05916),
        ({"I2": 0.01, "I-": 0.1}, 0.54 + 0.05916 / 2 * np.log10(0.01 / 0.1**2)),
    ],
)
def test_potentiel_nernst(composition, attendu):
    assert solution_redox(composition).E() == pytest.approx(attendu, abs=1e-6)


def test_concentrations_dimere():
    c = solution_redox({"I2": 0.01, "I-": 0.1}).concentrations()
    assert c["I2"] == pytest.approx(0.01) and c["I-"] == pytest.approx(0.1)


class TestTitrage:
    def test_points_remarquables(self, fer_cerium):
        E = fer_cerium.E([5, 10, 20])
        assert E == pytest.approx([0.77, (0.77 + 1.72) / 2, 1.72], abs=1e-4)

    def test_reference(self):
        t = TitrageRedox(SolutionRedox().avec(FER.red, 0.1), SolutionRedox().avec(CERIUM.ox, 0.1), 10, "ECS")
        assert float(t.E(5)) == pytest.approx(0.77 - 0.244, abs=1e-4)

    def test_reference_inconnue(self):
        with pytest.raises(ValueError):
            TitrageRedox(SolutionRedox(), SolutionRedox(), 10, "ESM")

    def test_permanganate_equivalence(self):
        t = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"KMnO4": 0.02}), 10)
        assert t.volumes_equivalents() == pytest.approx((10,))
        assert float(t.E(10)) == pytest.approx((0.77 + 5 * 1.51) / 6, abs=2e-3)

    def test_vanadium_trois_sauts(self):
        t = TitrageRedox(solution_redox({"V2+": 0.1}), solution_redox({"KMnO4": 0.02}), 10)
        assert t.volumes_equivalents() == pytest.approx((10, 20, 30))
        assert [round(s.volume, 3) for s in t.sauts(40)] == [10, 20, 30]

    def test_titrant_reducteur(self):
        t = TitrageRedox(solution_redox({"Ce4+": 0.1}), solution_redox({"Fe2+": 0.1}), 10)
        assert t.volumes_equivalents() == pytest.approx((10,))
        assert np.all(np.diff(t.courbe(20).E) < 0)

    @pytest.mark.parametrize(
        ("becher", "burette", "equivalence"),
        [
            ({"Fe2+": 0.1}, {"Ce4+": 0.1}, 10),
            ({"Fe2+": 0.1}, {"K2Cr2O7": 1 / 60}, 10),  # un dimère, dans la burette
            ({"K2Cr2O7": 1 / 60}, {"Fe2+": 0.1}, 10),  # un dimère, dans le bécher
            ({"Hg2^2+": 0.05}, {"Ce4+": 0.1}, 10),  # dimère réduit (2, 1)
            ({"I2": 0.05}, {"Na2S2O3": 0.1}, 10),  # deux dimères : l'iodométrie
        ],
    )
    def test_volume_et_potentiel_reciproques(self, becher, burette, equivalence):
        t = TitrageRedox(solution_redox(becher), solution_redox(burette), 10)
        assert t.volumes_equivalents()[0] == pytest.approx(equivalence)
        assert t.volume(t.E(VOLUMES)) == pytest.approx(VOLUMES, rel=1e-9)

    def test_volume_hors_domaine(self):
        # une burette qui contient les deux formes a un potentiel : on ne le dépasse pas
        t = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"Ce4+": 0.1, "Ce3+": 0.1}), 10)
        assert np.isnan(t.volume(1.8))
        assert np.isfinite(t.volume(1.7))

    def test_potentiel_initial_indefini(self, fer_cerium):
        assert fer_cerium.volume_initial(20) == pytest.approx(0.02)
        melange = TitrageRedox(
            solution_redox({"Fe2+": 0.1, "Fe3+": 0.001}), solution_redox({"Ce4+": 0.1}), 10
        )
        assert melange.volume_initial(20) == 0

    def test_pente_implicite(self, fer_cerium):
        V, h = np.array([3.0, 7.0, 13.0]), 1e-4
        fd = (fer_cerium.E(V + h) - fer_cerium.E(V - h)) / (2 * h)
        assert fer_cerium.pente(V) == pytest.approx(fd, rel=1e-4)

    def test_courbe(self, fer_cerium):
        c = fer_cerium.courbe(20, 100)
        assert len(c) == 100 and c.volume[0] == pytest.approx(0.02) and c.volume[-1] == 20
        assert np.all(np.diff(c.E) > 0)
        assert ((c.E > 0.9) & (c.E < 1.6)).sum() >= 10

    def test_saut(self, fer_cerium):
        (s,) = fer_cerium.sauts(20)
        assert s.volume == pytest.approx(10, abs=1e-6)
        assert s.E == pytest.approx(1.245, abs=1e-3)

    def test_pH_different_refuse(self):
        with pytest.raises(ValueError):
            TitrageRedox(SolutionRedox(pH=0), SolutionRedox(pH=1), 10)


def test_base_redox():
    assert espece_redox("permanganate de potassium").q == 1
    assert espece_redox("ion vanadyle").formule == "VO^2+"
    assert couple_redox("I-").stoechiometrie == (1, 2)
    with pytest.raises(KeyError):
        espece_redox("VO2+")  # ambigu : retiré de l'index


def test_couple_incoherent():
    with pytest.raises(ValueError):
        CoupleRedox("x", (0.5, 1.0), (1, 1), ("A", "B", "C"), stoechiometrie=(1, 2))
