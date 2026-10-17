from titrages import Titrage, TitrageRedox, solution, solution_redox, vers_pgfplots


def test_pgfplots(acetique):
    texte = vers_pgfplots(acetique.courbe(20, 5), chiffres=3)
    assert texte.startswith("(0,2.88) ")
    assert texte.count("(") == 5


def test_export_redox_et_conductimetrie():
    r = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"Ce4+": 0.1}), 10).courbe(20, 5)
    assert vers_pgfplots(r, 3).startswith("(0.02,")  # Fe2+ seul : pas de potentiel à V = 0
    c = Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), 100).courbe_conductimetrique(20, 5)
    texte = vers_pgfplots(c, 3)
    assert texte.count("(") == 5
    assert texte.startswith("(0,")  # le volume d'abord, puis σ


def test_export_complexometrie_et_precipitation():
    from titrages import (
        Complexe,
        Precipite,
        SolutionComplexation,
        SolutionPrecipitation,
        TitrageComplexometrique,
        TitrageParPrecipitation,
    )

    c = TitrageComplexometrique(
        SolutionComplexation({"Ca2+": 0.01}),
        SolutionComplexation(ligand=0.01),
        20,
        [Complexe("Ca2+", (8.0,))],
    ).courbe(40, 5)
    assert vers_pgfplots(c, 3).startswith("(4e-08,2)")
    p = TitrageParPrecipitation(
        SolutionPrecipitation(partenaires={"Cl-": 0.1}),
        SolutionPrecipitation(commun=0.1),
        20,
        [Precipite("Cl-", 9.75)],
    ).courbe(40, 5)
    texte = vers_pgfplots(p, 3)
    assert texte.count("(") == 5 and texte.endswith(",1.48)")
