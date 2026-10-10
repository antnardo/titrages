from titrages.__main__ import main


def test_pH(capsys):
    assert main(["pH", "CH3COOH=0.1"]) == 0
    assert capsys.readouterr().out.strip() == "2.88"


def test_titrage(capsys):
    assert main(["titrage", "H3PO4=0.1", "--burette", "NaOH=0.1"]) == 0
    sortie = capsys.readouterr().out
    assert "10, 20, 30" in sortie and sortie.count("saut") == 2


def test_couples_filtre(capsys):
    main(["couples", "phosph"])
    assert capsys.readouterr().out.startswith("acide phosphorique")


def test_indicateurs(capsys):
    main(["indicateurs"])
    assert "phénolphtaléine" in capsys.readouterr().out


def test_inconnu(capsys):
    assert main(["pH", "XYZ=0.1"]) == 1


def test_sigma(capsys):
    assert main(["sigma", "HCl=0.01"]) == 0
    assert capsys.readouterr().out.strip() == "0.4261"


def test_potentiel(capsys):
    assert main(["E", "Fe2+=0.1", "Fe3+=0.1"]) == 0
    assert capsys.readouterr().out.strip() == "0.770"


def test_redox(capsys):
    assert main(["redox", "V2+=0.1", "--burette", "KMnO4=0.02", "--reference", "ECS"]) == 0
    sortie = capsys.readouterr().out
    assert "10, 20, 30" in sortie and sortie.count("saut") == 3


def test_redox_couples(capsys):
    main(["redox-couples"])
    assert "vanadium" in capsys.readouterr().out
