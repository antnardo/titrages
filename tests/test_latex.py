"""Les deux moteurs du paquet LaTeX donnent les nombres du paquet Python.

``latex/tests/titrages-test-nombres.tex``, compilé par lualatex (moteur lua) et par
pdflatex (moteur pgf), écrit un fichier de nombres — V(pH), pH(V), σ(V), E(V), V(E) —
que l'on compare ici au calcul Python des mêmes titrages. Sous pgf, les lignes
« pgfmath » évaluent l'expression même que pgfplots trace, avec son arithmétique en
virgule flottante (fpu, quatre à cinq chiffres) : c'est la vérification que les
courbes du moteur pgf sont celles du moteur lua, à l'épaisseur du trait près. Ces tests
demandent une distribution TeX ; ``uv run pytest -m "not latex"`` les laisse de côté.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from titrages import Titrage, TitrageRedox, solution, solution_redox

LATEX = Path(__file__).parent.parent / "latex"
DOCUMENT = LATEX / "tests" / "titrages-test-nombres.tex"
COMPILATEURS = {"lua": "lualatex", "pgf": "pdflatex"}

pytestmark = pytest.mark.latex


def _titrages() -> dict[str, Titrage | TitrageRedox]:
    """Les titrages définis dans le document, écrits ici avec le paquet Python."""
    return {
        "acetique": Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), 10),
        "phosphorique": Titrage(solution({"H3PO4": 0.1}), solution({"NaOH": 0.1}), 10),
        "ammoniac": Titrage(solution({"NH3": 0.1}), solution({"HCl": 0.1}), 10),
        "fort": Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), 100),
        "faible": Titrage(solution({"HCl": 0.05, "CH3COOH": 0.05}), solution({"NH3": 0.1}), 20),
        "sansdilution": Titrage(solution({"H2CO3": 0.01}), solution({"NaOH": 0.1}), 10, dilution=False),
        "cerium": TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"Ce4+": 0.1}), 10, "ECS"),
        "vanadium": TitrageRedox(
            solution_redox({"V2+": 0.1}, pH=1), solution_redox({"KMnO4": 0.02}, pH=1), 10
        ),
    }


def _attendu(t: Titrage | TitrageRedox, nom: str, grandeur: str, argument: float) -> float:
    if grandeur in ("volume", "pgfmath"):
        return float(t.volume(argument))
    if grandeur == "E":
        assert isinstance(t, TitrageRedox)
        return float(t.E(argument))
    assert isinstance(t, Titrage)
    if grandeur == "pH":
        return float(t.pH(argument))
    if grandeur == "pHinitial":
        return t.pH_initial
    assert grandeur == "conductivite"
    sigma = float(t.conductivite(argument))
    if nom == "fort":  # conductivite corrigee = true dans le document
        sigma *= (t.volume_becher + argument) / t.volume_becher
    return sigma


@pytest.fixture(scope="module", params=list(COMPILATEURS))
def nombres(request, tmp_path_factory):
    """Le fichier de nombres écrit par le document, compilé par le moteur demandé."""
    moteur = request.param
    compilateur = shutil.which(COMPILATEURS[moteur])
    if compilateur is None:
        pytest.skip(f"{COMPILATEURS[moteur]} absent")
    dossier = tmp_path_factory.mktemp(moteur)
    shutil.copy(DOCUMENT, dossier)
    subprocess.run(
        [compilateur, "-interaction=nonstopmode", "-halt-on-error", DOCUMENT.name],
        cwd=dossier,
        env={**os.environ, "TEXINPUTS": f"{LATEX}:"},
        capture_output=True,
        check=True,
    )
    lignes = (dossier / f"{DOCUMENT.stem}-{moteur}.dat").read_text(encoding="utf-8").split("\n")
    return moteur, [ligne.split() for ligne in lignes if ligne.strip()]


def test_memes_nombres_que_python(nombres):
    moteur, lignes = nombres
    titrages = _titrages()
    assert len(lignes) >= 50
    for nom, grandeur, argument, valeur in lignes:
        attendu = _attendu(titrages[nom], nom, grandeur, float(argument))
        if grandeur == "pgfmath":  # l'arithmétique de pgfplots : quatre à cinq chiffres
            assert float(valeur) == pytest.approx(attendu, rel=2e-3), (moteur, nom, grandeur, argument)
        else:
            assert float(valeur) == pytest.approx(attendu, abs=1e-9), (moteur, nom, grandeur, argument)


def test_pgfmath_seulement_sous_pgf(nombres):
    """Les courbes de pgfplots ne sont vérifiées que là où pgfplots les calcule lui-même."""
    moteur, lignes = nombres
    assert any(grandeur == "pgfmath" for _, grandeur, _, _ in lignes) == (moteur == "pgf")
