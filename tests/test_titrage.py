import numpy as np
import pytest

from titrages import Titrage, solution


def test_pH_aux_points_remarquables(acetique):
    pH = acetique.pH([0, 5, 10])
    assert pH == pytest.approx([2.88, 4.76, 8.73], abs=0.01)


def test_asymptote_de_la_burette(acetique):
    assert float(acetique.pH(1e6)) == pytest.approx(13, abs=1e-3)


def test_volume_et_pH_sont_reciproques(phosphorique):
    V = np.array([0.5, 5, 10, 15, 22, 29])
    assert phosphorique.volume(phosphorique.pH(V)) == pytest.approx(V, rel=1e-9)


def test_volume_hors_domaine_est_nan(acetique):
    assert np.isnan(acetique.volume(1.0))  # sous le pH initial
    assert np.isnan(acetique.volume(13.5))  # au-delà de la burette


def test_formule_de_la_fiche():
    # n = (Y1 + δ/CT) / (1 - δ/CB), ici pour un diacide
    from titrages import delta, moments

    CT, CB, V0 = 0.01, 0.1, 10
    t = Titrage(solution({"H2CO3": CT}), solution({"NaOH": CB}), V0)
    pH = np.linspace(5, 12, 20)
    (Y1,) = moments([6.35, 10.33], pH, 1)
    d = delta(pH)
    n = (Y1 + d / CT) / (1 - d / CB)
    assert t.volume(pH) * CB / (CT * V0) == pytest.approx(n)


def test_sans_dilution_formule_de_la_fiche():
    from titrages import delta, moments

    CT, CB, V0 = 0.01, 0.1, 10
    t = Titrage(solution({"H2CO3": CT}), solution({"NaOH": CB}), V0, dilution=False)
    pH = np.linspace(5, 11, 20)
    (Y1,) = moments([6.35, 10.33], pH, 1)
    assert t.volume(pH) * CB / (CT * V0) == pytest.approx(Y1 + delta(pH) / CT)


def test_derivee_volume_analytique(phosphorique):
    pH, h = np.array([3.0, 5.0, 8.0, 11.0]), 1e-6
    fd = (phosphorique.volume(pH + h) - phosphorique.volume(pH - h)) / (2 * h)
    assert phosphorique.derivee_volume(pH) == pytest.approx(fd, rel=1e-5)


def test_derivee_volume_titrant_faible():
    t = Titrage(solution({"CH3COOH": 0.1}), solution({"NH3": 0.1}), 10)
    pH, h = np.array([4.0, 7.0, 8.5]), 1e-6
    fd = (t.volume(pH + h) - t.volume(pH - h)) / (2 * h)
    assert t.derivee_volume(pH) == pytest.approx(fd, rel=1e-5)


@pytest.mark.parametrize("volume", [-1.0])
def test_volume_negatif_refuse(acetique, volume):
    with pytest.raises(ValueError):
        acetique.pH(volume)


class TestCourbe:
    def test_points_exacts(self, acetique):
        c = acetique.courbe(20, 100)
        assert len(c) == 100
        assert c.volume[0] == 0 and c.volume[-1] == 20
        assert acetique.pH(c.volume[1:-1]) == pytest.approx(c.pH[1:-1], abs=1e-9)

    def test_croissante(self, acetique):
        c = acetique.courbe(20)
        assert np.all(np.diff(c.volume) > 0) and np.all(np.diff(c.pH) > 0)

    def test_le_saut_est_echantillonne(self, acetique):
        c = acetique.courbe(20, 200)
        dans_le_saut = (c.pH > 6) & (c.pH < 11)
        assert dans_le_saut.sum() >= 20

    def test_decroissante_pour_un_titrant_acide(self, ammoniac):
        assert ammoniac.sens == -1
        assert np.all(np.diff(ammoniac.courbe(20).pH) < 0)

    def test_volume_max_positif(self, acetique):
        with pytest.raises(ValueError):
            acetique.courbe(0)


class TestEquivalences:
    def test_monoacide(self, acetique):
        assert acetique.volumes_equivalents() == pytest.approx((10,))
        (s,) = acetique.sauts(20)
        assert s.volume == pytest.approx(10, abs=1e-3)
        assert s.pH == pytest.approx(8.73, abs=0.01)

    def test_troisieme_acidite_phosphorique_invisible(self, phosphorique):
        assert phosphorique.volumes_equivalents() == pytest.approx((10, 20, 30))
        sauts = phosphorique.sauts(35)
        assert [round(s.volume) for s in sauts] == [10, 20]

    def test_base_faible_par_acide_fort(self, ammoniac):
        assert ammoniac.volumes_equivalents() == pytest.approx((10,))
        (s,) = ammoniac.sauts(20)
        assert s.pH == pytest.approx(5.27, abs=0.01)

    def test_melange_acide_fort_acide_faible(self):
        t = Titrage(solution({"HCl": 0.1, "CH3COOH": 0.1}), solution({"NaOH": 0.1}), 10)
        assert t.volumes_equivalents() == pytest.approx((10, 20))
        assert [round(s.volume, 1) for s in t.sauts(30)] == [10.0, 20.0]

    def test_base_deja_presente_passe_en_premier(self):
        # 0,05 M de soude neutralise la moitié du premier proton avant la burette
        t = Titrage(solution({"H2CO3": 0.1, "NaOH": 0.05}), solution({"NaOH": 0.1}), 10)
        assert t.volumes_equivalents() == pytest.approx((5, 15))

    def test_forme_introduite(self):
        t = Titrage(solution({"Na2CO3": 0.1}), solution({"HCl": 0.1}), 10)
        assert t.volumes_equivalents() == pytest.approx((10, 20))

    def test_titrant_faible_sans_equivalence_stoechiometrique(self):
        t = Titrage(solution({"CH3COOH": 0.1}), solution({"NH3": 0.1}), 10)
        with pytest.raises(ValueError):
            t.volumes_equivalents()
        (s,) = t.sauts(20)
        assert s.volume == pytest.approx(10, abs=0.01)
        assert s.pH == pytest.approx(7.0, abs=0.01)

    def test_burette_vide(self):
        with pytest.raises(ValueError):
            Titrage(solution({"HCl": 0.1}), solution(), 10).volumes_equivalents()


class TestTangentes:
    def test_proche_de_l_equivalence(self, acetique):
        t = acetique.tangentes(7, 20)
        assert t.volume == pytest.approx(10, abs=0.01)
        assert t.pente == pytest.approx(float(acetique.pente(t.pH_apres)), rel=1e-6)
        assert t.volume_avant < t.volume < t.volume_apres

    def test_titrant_acide(self, ammoniac):
        t = ammoniac.tangentes(7, 20)
        assert t.volume == pytest.approx(10, abs=0.01)
        assert t.pH_apres < t.pH < t.pH_avant

    def test_aucun_saut_apres(self, acetique):
        with pytest.raises(ValueError):
            acetique.tangentes(15, 20)


def test_pKe_different_refuse():
    from titrages import Solution

    with pytest.raises(ValueError):
        Titrage(Solution(pKe=14), Solution(pKe=13), 10)


@pytest.mark.parametrize(
    "becher, burette, attendu",
    [
        ({"CH3COOH": 0.1}, {"NaOH": 0.1}, 8.729420699201),
        ({"HCl": 0.1}, {"NaOH": 0.1}, 6.999998262822),
        ({"NH3": 0.1}, {"HCl": 0.1}, 5.275575942902),  # pH décroissant
    ],
)
def test_inflexion_a_la_precision_de_la_machine(becher, burette, attendu):
    """Face à mpmath : le zéro de d²V/dpH² de V(pH) = -V0 R_bécher/R_burette, 60 chiffres."""
    t = Titrage(solution(becher), solution(burette), 10)
    (s,) = t.sauts(20)
    assert s.pH == pytest.approx(attendu, abs=1e-9)
