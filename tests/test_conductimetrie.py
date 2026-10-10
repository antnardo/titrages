import numpy as np
import pytest

from titrages import Couple, Solution, Spectateur, Titrage, solution


def test_acide_fort_kohlrausch():
    assert float(solution({"HCl": 0.01}).conductivite()) == pytest.approx((34.98 + 7.63) * 0.01, rel=1e-6)


def test_contre_ions_des_sels():
    s = solution({"Na2CO3": 0.1, "NH4Cl": 0.05, "CH3COO-": 0.1})
    ions: dict[str, float] = {}
    for i in s.spectateurs:
        ions[i.formule] = ions.get(i.formule, 0) + i.concentration
    assert ions == pytest.approx({"Na+": 0.3, "Cl-": 0.05})


def test_couverture_de_la_base():
    """Toute forme chargée de la base a son λ°, sauf les quatre de l'EDTA (complexométrie)."""
    from titrages import charge, couples

    manquants = [
        f
        for c in couples()
        for f, lam in zip(c.formes, c.conductivites, strict=True)
        if charge(f) and lam is None
    ]
    assert manquants == ["H3Y-", "H2Y^2-", "HY^3-", "Y^4-"]
    assert all(
        lam is None
        for c in couples()
        for f, lam in zip(c.formes, c.conductivites, strict=True)
        if not charge(f)
    )


def test_glycine_zwitterion_neutre():
    """La glycine n'apporte ni contre-ion ni λ° : elle est neutre ; ses deux autres formes conduisent."""
    s = solution({"glycine": 0.1})
    assert s.spectateurs == ()
    assert s.conductivite() > 0
    assert solution({"ion glycinium": 0.1}).spectateurs[0].formule == "Cl-"


def test_conductivite_inconnue_signalee():
    s = Solution().avec(Couple.polyacide(4.0), 0.1)
    with pytest.raises(ValueError, match="A-"):
        s.conductivite()


def test_couple_avec_conductivites():
    c = Couple("x", (4.0,), ("HA", "A-"), conductivites=(None, 4.0))
    s = Solution().avec(c, 0.1).avec_spectateur(Spectateur("Na+", 0.0, 5.01))
    assert float(s.conductivite()) > 0


class TestTitrageConductimetrique:
    @pytest.fixture
    def fort(self):
        return Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), 100)

    def test_valeurs_corrigees(self, fort):
        c = fort.courbe_conductimetrique(20, 3, corrigee=True)
        attendu = [(34.98 + 7.63) * 0.01, (5.01 + 7.63) * 0.01, 5.01 * 0.02 + 7.63 * 0.01 + 19.86 * 0.01]
        assert c.conductivite == pytest.approx(attendu, rel=1e-3)

    def test_segments_rectilignes_une_fois_corriges(self, fort):
        c = fort.courbe_conductimetrique(8, 9, corrigee=True)
        assert np.allclose(np.diff(c.conductivite, 2), 0, atol=1e-5)

    def test_minimum_a_l_equivalence(self, fort):
        c = fort.courbe_conductimetrique(20, 201)
        assert c.volume[np.argmin(c.conductivite)] == pytest.approx(10, abs=0.1)

    def test_acide_faible_croissant(self):
        t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), 10)
        c = t.courbe_conductimetrique(20, 41, corrigee=True)
        pente_avant = np.diff(c.conductivite[5:15]).mean()
        pente_apres = np.diff(c.conductivite[25:35]).mean()
        assert 0 < pente_avant < pente_apres

    def test_coherent_avec_le_melange(self):
        t = Titrage(solution({"H3PO4": 0.05}), solution({"NaOH": 0.1}), 10)
        V = 7.0
        melange = t.becher.melangee(t.burette, 10, V)
        assert float(t.conductivite(V)) == pytest.approx(float(melange.conductivite()), rel=1e-9)
