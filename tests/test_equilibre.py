import numpy as np
import pytest

from titrages import Couple, concentration_pour_pH, fractions, moments
from titrages.equilibre import delta, log_k, racine_croissante


def test_log_k_cumule_les_pKa():
    assert log_k([6.35, 10.33]).tolist() == pytest.approx([0, -6.35, -16.68])


@pytest.mark.parametrize("pKa", [[4.76], [6.35, 10.33], [2.15, 7.20, 12.35], [2.0, 2.67, 6.16, 10.26]])
def test_fractions_somment_a_un(pKa):
    x = fractions(pKa, np.linspace(-2, 16, 50))
    assert x.shape == (len(pKa) + 1, 50)
    assert np.allclose(x.sum(axis=0), 1)


def test_fractions_egales_au_pKa():
    x = fractions([6.35, 10.33], [6.35, 10.33])
    assert x[0, 0] == pytest.approx(x[1, 0], rel=1e-3)
    assert x[1, 1] == pytest.approx(x[2, 1], rel=1e-3)


def test_fractions_sans_debordement_aux_extremes():
    # la forme naïve déborde : 10^(4 × 30) pour l'EDTA à pH 30
    x = fractions([2.0, 2.67, 6.16, 10.26], [-30.0, 30.0])
    assert np.all(np.isfinite(x))
    assert x[0, 0] == pytest.approx(1) and x[-1, 1] == pytest.approx(1)


def test_fractions_ordre_des_pKa_indifferent():
    assert np.allclose(fractions([10.33, 6.35], 8), fractions([6.35, 10.33], 8))


def test_fractions_conservent_la_forme_du_pH():
    assert fractions([4.76], np.zeros((3, 4))).shape == (2, 3, 4)


def test_moments_ampholyte_au_point_isoelectrique():
    (Y1,) = moments([6.35, 10.33], (6.35 + 10.33) / 2, ordre=1)
    assert Y1 == pytest.approx(1)


def test_moments_derivee_de_Y1_est_la_variance():
    p, h = 7.0, 1e-6
    (Ya,), (Yb,) = moments([6.35, 10.33], p - h, 1), moments([6.35, 10.33], p + h, 1)
    Y1, Y2 = moments([6.35, 10.33], p)
    assert (Yb - Ya) / (2 * h) == pytest.approx(np.log(10) * (Y2 - Y1**2), rel=1e-6)


def test_scalaire_rend_float_tableau_rend_tableau():
    assert isinstance(delta(7.0), float) and isinstance(delta([7.0, 8.0]), np.ndarray)
    (Y1,) = moments([4.76], 4.76, 1)
    assert isinstance(Y1, float)
    assert isinstance(concentration_pour_pH(Couple("acide", (4.0,)).acide, 3.0), float)
    assert fractions([4.76], 4.76).shape == (2,)  # les fractions : un tableau, une valeur par forme


def test_delta_nul_a_la_neutralite():
    assert delta(7.0) == pytest.approx(0, abs=1e-20)
    assert delta(6.5, pKe=13.0) == pytest.approx(0, abs=1e-20)


def test_concentration_pour_pH_retrouve_le_pH():
    from titrages import Solution

    acide = Couple("acide", (4.0,)).acide
    C = float(concentration_pour_pH(acide, 3.0))
    assert Solution().avec(acide, C).pH() == pytest.approx(3.0)


def test_concentration_pour_pH_impossible_donne_nan():
    assert np.isnan(concentration_pour_pH(Couple("acide", (4.0,)).acide, 9.0))


def test_racine_croissante_vectorisee():
    r = racine_croissante(lambda x: x**3 - np.array([1.0, 8.0, 27.0]), 0.0, 10.0)
    assert r == pytest.approx([1, 2, 3])
