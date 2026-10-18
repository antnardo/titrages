"""La base des ligands, des complexes et des précipités."""

import math

import pytest

from titrages import (
    SolutionComplexation,
    SolutionPrecipitation,
    TitrageComplexometrique,
    TitrageParPrecipitation,
    complexe,
    complexes,
    ligand,
    ligands,
    precipite,
    precipites,
)


class TestLigands:
    def test_recherche_par_nom_formule_ou_alias(self):
        assert ligand("EDTA") is ligand("Y4-") is ligand("acide éthylènediaminetétraacétique")

    def test_pKa_de_la_source(self):
        assert ligand("ammoniac").pKa == (9.244,)

    def test_inconnu(self):
        with pytest.raises(KeyError, match="glucose"):
            ligand("glucose")

    def test_tous_ont_des_complexes(self):
        assert all(complexes(lig.nom) for lig in ligands())


class TestComplexes:
    def test_constantes_successives_cumulees(self):
        assert complexe("Zn2+", "NH3").log_beta == pytest.approx((2.21, 4.50, 6.86, 8.89))

    def test_constante_globale_seule(self):
        # la source ne donne que log beta_6 pour l'hexacyanoferrate(III)
        c = complexe("Fe3+", "CN-")
        assert c.log_beta[:5] == (-math.inf,) * 5 and c.log_beta[5] == 43.6

    def test_constantes_successives_croissantes_une_seule_etape(self):
        """Zn–NH3 : log K3 > log K2 > log K1, ZnL et ZnL2 ne prédominent jamais."""
        t = TitrageComplexometrique(
            SolutionComplexation({"Zn2+": 0.01}),
            SolutionComplexation(ligand=0.1),
            20,
            [complexe("Zn2+", "NH3")],
            ligand("NH3"),
            pH=11,
        )
        # trois ligands d'un coup (6 mL), puis le quatrième (8 mL)
        assert t.volumes_equivalents() == pytest.approx((6.0, 8.0))

    def test_metal_sans_complexe(self):
        with pytest.raises(KeyError, match="Ba2"):
            complexe("Ba2+", "NH3")

    def test_note_de_conditions(self):
        assert "5,0" in complexe("Ag+", "Cl-").note

    def test_dosage_de_liebig(self):
        """Ag+ dans les cyanures : Ag(CN)2-, deux cyanures par argent, ML1 absent de la source."""
        t = TitrageComplexometrique(
            SolutionComplexation(ligand=0.1),
            SolutionComplexation({"Ag+": 0.05}),
            20,
            [complexe("Ag+", "cyanure")],
            ligand("cyanure"),
            pH=12,
        )
        # Ag(CN)3^2- puis Ag(CN)2- ; AgCN, que la source ne donne pas, ne se forme pas.
        # Seule la seconde a un saut, comme la troisième acidité de H3PO4.
        assert t.volumes_equivalents() == pytest.approx((40 / 3, 20.0))
        (saut,) = t.sauts(40)
        assert saut.volume == pytest.approx(20, abs=1e-3)


class TestPrecipites:
    def test_recherche(self):
        assert precipite("AgCl") is precipite("chlorure d'argent")

    def test_stoechiometrie(self):
        p = precipite("PbI2")
        assert (p.commun, p.partenaire, p.nu, p.mu, p.pKs) == ("Pb2+", "I-", 1, 2, 8.1)

    def test_formules_uniques(self):
        noms = [p.nom for p in precipites()]
        assert len(noms) == len(set(noms)) == 56

    def test_ion_commun_deduit(self):
        t = TitrageParPrecipitation(
            SolutionPrecipitation(partenaires={"SO4^2-": 0.01}),
            SolutionPrecipitation(commun=0.05),
            20,
            [precipite("BaSO4")],
            pH=7,
        )
        assert t.ion_commun == "Ba2+" and t.volumes_equivalents() == pytest.approx((4.0,))

    def test_ions_communs_differents(self):
        with pytest.raises(ValueError, match="ion commun"):
            TitrageParPrecipitation(
                SolutionPrecipitation(partenaires={"Cl-": 0.01}),
                SolutionPrecipitation(commun=0.05),
                20,
                [precipite("AgCl"), precipite("PbI2")],
            )

    def test_mohr_depuis_la_base(self):
        t = TitrageParPrecipitation(
            SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
            SolutionPrecipitation(commun=0.05),
            20,
            [precipite("AgCl"), precipite("Ag2CrO4")],
            pH=7.5,
        )
        chromate = [e for e in t.evenements(30) if e.partenaire == "CrO4^2-"]
        assert len(chromate) == 1 and 20 < chromate[0].volume < 20.1
