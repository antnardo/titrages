"""Les titrages par précipitation, confrontés à un calcul indépendant.

La référence dilue les concentrations à la main, au volume versé, et résout le bilan de
l'ion commun du mélange par dichotomie en multiprécision, chaque solide présent ou non
selon qu'il est sursaturé : elle ne connaît ni le paramétrage par pP, ni l'énumération des
ensembles de solides du paquet.
"""

import numpy as np
import pytest
from mpmath import mp, mpf

from titrages import Precipite, SolutionPrecipitation, TitrageParPrecipitation

mp.dps = 50

CL = Precipite("Cl-", 9.75)
IODURE = Precipite("I-", 16.08)
SCN = Precipite("SCN-", 12.0)
CRO4 = Precipite("CrO4^2-", 11.95, nu=2, pKa=(6.5,))


def reference_p(t: TitrageParPrecipitation, V: float) -> float:
    """pP après avoir versé V, sans rien emprunter au paquet."""
    V, V0 = mpf(V), mpf(t.volume_becher)
    TP = (V0 * mpf(t.becher.commun) + V * mpf(t.burette.commun)) / (V0 + V)
    termes = []
    for p in t.precipites:
        TX = (
            V0 * mpf(t.becher.concentration(p.partenaire)) + V * mpf(t.burette.concentration(p.partenaire))
        ) / (V0 + V)
        alpha = mpf(1)
        if p.pKa:
            h = mp.power(10, -mpf(t.pH))
            alpha = 1 / (1 + h / mp.power(10, -mpf(p.pKa[-1])))
        termes.append((TX, mpf(p.nu) / p.mu, p, alpha))

    def bilan(log_a):
        a = mp.power(10, log_a)
        total = a - TP
        for TX, r, p, alpha in termes:
            x_sat = mp.power(mp.power(10, -mpf(p.pKs)) / a**p.nu, mpf(1) / p.mu) / alpha
            total += r * max(mpf(0), TX - x_sat)
        return total

    bas, haut = mpf(-60), mpf(1)
    for _ in range(400):
        milieu = (bas + haut) / 2
        if bilan(milieu) > 0:
            haut = milieu
        else:
            bas = milieu
    return float(-(bas + haut) / 2)


CAS = {
    "chlorure par l'argent": TitrageParPrecipitation(
        SolutionPrecipitation(partenaires={"Cl-": 0.1}), SolutionPrecipitation(commun=0.1), 20, [CL]
    ),
    "chlorure et iodure": TitrageParPrecipitation(
        SolutionPrecipitation(partenaires={"Cl-": 0.05, "I-": 0.05}),
        SolutionPrecipitation(commun=0.1),
        20,
        [CL, IODURE],
    ),
    "Mohr, pH 7,5": TitrageParPrecipitation(
        SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
        SolutionPrecipitation(commun=0.05),
        20,
        [CL, CRO4],
        pH=7.5,
    ),
    "Charpentier-Volhard": TitrageParPrecipitation(
        SolutionPrecipitation(commun=0.1, partenaires={"Cl-": 0.05}),
        SolutionPrecipitation(partenaires={"SCN-": 0.1}),
        20,
        [CL, SCN],
    ),
    "argent par le chlorure": TitrageParPrecipitation(
        SolutionPrecipitation(commun=0.05), SolutionPrecipitation(partenaires={"Cl-": 0.1}), 20, [CL]
    ),
}


class TestTitrageParPrecipitation:
    @pytest.mark.parametrize("nom", CAS)
    @pytest.mark.parametrize("V", [0.5, 5.0, 9.99, 10.0, 10.01, 15.0, 19.99, 20.0, 20.01, 24.0, 30.0])
    def test_p_identique_a_la_reference(self, nom, V):
        t = CAS[nom]
        assert t.p(V) == pytest.approx(reference_p(t, V), abs=1e-9)

    @pytest.mark.parametrize("nom", CAS)
    def test_volume_et_p_reciproques(self, nom):
        t = CAS[nom]
        V = np.array([1.0, 7.0, 13.0, 19.0, 22.0, 28.0])
        assert t.volume(t.p(V)) == pytest.approx(V, rel=1e-10)

    @pytest.mark.parametrize("nom", CAS)
    def test_pente_egale_aux_differences_finies(self, nom):
        t = CAS[nom]
        for V in (3.0, 12.0, 26.0):
            h = 1e-5
            attendu = (t.p(V + h) - t.p(V - h)) / (2 * h)
            assert t.pente(t.p(V)) == pytest.approx(attendu, rel=1e-5)

    def test_equivalence_du_chlorure_d_argent(self):
        # À l'équivalence, [Ag+] = [Cl-] et le solide est présent : pAg = pKs/2.
        assert CAS["chlorure par l'argent"].p(20) == pytest.approx(9.75 / 2, abs=1e-12)

    def test_saut_a_l_equivalence(self):
        (saut,) = CAS["chlorure par l'argent"].sauts(40)
        assert saut.volume == pytest.approx(20, abs=1e-4)


class TestEvenements:
    def test_mohr_le_chromate_apparait_apres_l_equivalence(self):
        t = CAS["Mohr, pH 7,5"]
        chromate = [e for e in t.evenements(30) if e.partenaire == "CrO4^2-"]
        assert len(chromate) == 1 and chromate[0].apparition and 20 < chromate[0].volume < 20.05

    def test_evenement_sur_la_courbe(self):
        t = CAS["Mohr, pH 7,5"]
        for e in t.evenements(30):
            assert t.p(e.volume) == pytest.approx(e.p, abs=1e-8)

    def test_iodure_puis_chlorure(self):
        (e,) = CAS["chlorure et iodure"].evenements(30)
        assert (e.partenaire, e.apparition, round(e.volume, 4)) == ("Cl-", True, 10.0)

    def test_volhard_le_chlorure_d_argent_se_redissout(self):
        t = CAS["Charpentier-Volhard"]
        noms = [(e.partenaire, e.apparition) for e in t.evenements(25)]
        assert noms == [("SCN-", True), ("Cl-", False)]

    def test_solides_presents(self):
        t = CAS["chlorure et iodure"]
        assert (t.solides(t.p(5.0)), t.solides(t.p(15.0))) == (("I-",), ("Cl-", "I-"))


class TestEquivalences:
    def test_l_iodure_precipite_d_abord(self):
        assert CAS["chlorure et iodure"].volumes_equivalents() == pytest.approx((10.0, 20.0))

    def test_chromate_compte_deux_argent(self):
        assert CAS["Mohr, pH 7,5"].volumes_equivalents() == pytest.approx((20.0, 24.0))

    def test_volhard_titre_l_argent_restant(self):
        assert CAS["Charpentier-Volhard"].volumes_equivalents() == pytest.approx((10.0,))

    def test_partenaire_titrant(self):
        assert CAS["argent par le chlorure"].volumes_equivalents() == pytest.approx((10.0,))


class TestErreurs:
    def test_partenaire_acide_sans_pH_le_pH_se_calcule(self):
        # ce n'est plus une erreur : le pH se calcule (niveau 2), et un avertissement le dit
        with pytest.warns(UserWarning, match="pH non précisé"):
            TitrageParPrecipitation(
                SolutionPrecipitation(partenaires={"CrO4^2-": 0.01}),
                SolutionPrecipitation(commun=0.1),
                20,
                [CRO4],
            )

    def test_partenaire_sans_precipite(self):
        with pytest.raises(ValueError, match="Br-"):
            TitrageParPrecipitation(
                SolutionPrecipitation(partenaires={"Br-": 0.01}), SolutionPrecipitation(commun=0.1), 20, [CL]
            )

    def test_p_partenaire_au_contact_du_solide(self):
        t = CAS["chlorure par l'argent"]
        assert t.p(10.0) + t.p_partenaire(10.0, "Cl-") == pytest.approx(9.75)
