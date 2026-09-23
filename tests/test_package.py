import titrages


def test_all_exporte():
    for nom in titrages.__all__:
        assert hasattr(titrages, nom)


def test_version():
    assert titrages.__version__.count(".") == 2
