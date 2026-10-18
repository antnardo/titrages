"""Les exemples des docstrings et de la documentation sont exécutés : ils ne mentent pas."""

import doctest
import importlib
from pathlib import Path

import pytest

MODULES = [
    "titrages",
    "titrages.base",
    "titrages.complexation",
    "titrages.complexation_libre",
    "titrages.equilibre",
    "titrages.modeles",
    "titrages.precipitation",
    "titrages.redox",
    "titrages.solution",
    "titrages.titrage",
]
RACINE = Path(__file__).parent.parent


@pytest.mark.parametrize("nom", MODULES)
def test_docstrings(nom):
    resultat = doctest.testmod(importlib.import_module(nom), optionflags=doctest.NORMALIZE_WHITESPACE)
    assert resultat.failed == 0


@pytest.mark.parametrize("fichier", ["README.md", "docs/DOC.md"])
def test_documentation(fichier):
    chemin = RACINE / fichier
    if not chemin.exists():
        pytest.skip(f"{fichier} absent")
    resultat = doctest.testfile(str(chemin), module_relative=False, optionflags=doctest.NORMALIZE_WHITESPACE)
    assert resultat.failed == 0
