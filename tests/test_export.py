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
