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


# --- complexométrie et précipitation, moteur lua seulement ---------------------------

DOCUMENT_CP = LATEX / "tests" / "titrages-test-cp.tex"


def _titrages_cp() -> dict:
    from titrages import (
        Complexe,
        Ligand,
        Precipite,
        SolutionComplexation,
        SolutionPrecipitation,
        TitrageComplexometrique,
        TitrageParPrecipitation,
    )

    edta = Ligand("Y4-", (0, 1.5, 2, 2.66, 6.16, 10.24))
    ca, zn = Complexe("Ca2+", (10.7,)), Complexe("Zn2+", (16.5,))
    cl = Precipite("Cl-", 9.75)
    return {
        "edta": TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01), 20, [ca], edta, 10
        ),
        "competition": TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01, "Zn2+": 0.005}),
            SolutionComplexation(ligand=0.01),
            20,
            [ca, zn],
            edta,
            10,
            "Zn2+",
        ),
        "cuivre": TitrageComplexometrique(
            SolutionComplexation(ligand=0.1),
            SolutionComplexation({"Cu2+": 0.05}),
            20,
            [Complexe("Cu2+", (4.1, 7.6, 10.5, 12.6))],
            Ligand("NH3", (9.25,)),
            11,
        ),
        "liebig": TitrageComplexometrique(
            SolutionComplexation(ligand=0.1),
            SolutionComplexation({"Ag+": 0.05}),
            20,
            [Complexe("Ag+", (-float("inf"), 20.48, 21.40))],
            Ligand("CN-", (9.21,)),
            12,
        ),
        "mohr": TitrageParPrecipitation(
            SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
            SolutionPrecipitation(commun=0.05),
            20,
            [cl, Precipite("CrO4^2-", 11.95, nu=2, pKa=(6.5,))],
            pH=7.5,
        ),
        "volhard": TitrageParPrecipitation(
            SolutionPrecipitation(commun=0.1, partenaires={"Cl-": 0.05}),
            SolutionPrecipitation(partenaires={"SCN-": 0.1}),
            20,
            [cl, Precipite("SCN-", 12.0)],
        ),
    }


# Le titrage de Volhard est suivi par pSCN (clé suivi) ; les autres par leur grandeur.
SUIVI_CP = {"volhard": "SCN-"}
VOLUME_MAX_CP = {"edta": 40, "competition": 40, "cuivre": 20, "liebig": 40, "mohr": 30, "volhard": 25}
STOECHIOMETRIE_CP = {"liebig": {"Ag+": 2}}


def _mesure_cp(t, nom: str, V: float, espece: str | None = None) -> float:
    from titrages import TitrageComplexometrique

    if isinstance(t, TitrageComplexometrique):
        return float(t.pM(V, espece))
    espece = espece or SUIVI_CP.get(nom)
    if espece and espece != t.ion_commun:
        return float(t.p_partenaire(V, espece))
    return float(t.p(V))


@pytest.fixture(scope="module")
def nombres_cp(tmp_path_factory):
    compilateur = shutil.which("lualatex")
    if compilateur is None:
        pytest.skip("lualatex absent")
    dossier = tmp_path_factory.mktemp("cp")
    shutil.copy(DOCUMENT_CP, dossier)
    subprocess.run(
        [compilateur, "-interaction=nonstopmode", "-halt-on-error", DOCUMENT_CP.name],
        cwd=dossier,
        env={**os.environ, "TEXINPUTS": f"{LATEX}:"},
        capture_output=True,
        check=True,
    )
    lignes = (dossier / f"{DOCUMENT_CP.stem}.dat").read_text(encoding="utf-8").split("\n")
    return [ligne.split() for ligne in lignes if ligne.strip()]


def test_complexation_precipitation_memes_nombres_que_python(nombres_cp):
    titrages = _titrages_cp()
    assert len(nombres_cp) >= 35
    for nom, grandeur, argument, valeur in nombres_cp:
        t, x, a = titrages[nom], float(valeur), float(argument)
        Vmax = VOLUME_MAX_CP[nom]
        if grandeur == "pM":
            attendu = _mesure_cp(t, nom, a)
        elif grandeur.startswith("pM"):  # pM<espèce> : une autre espèce que celle suivie
            attendu = _mesure_cp(t, nom, a, grandeur[2:])
        elif grandeur == "volume":  # le volume qui amène à ce pM : on y relit le pM
            x, attendu = _mesure_cp(t, nom, x), a
        elif grandeur == "equivalence":
            if nom in STOECHIOMETRIE_CP:
                attendu = t.volumes_equivalents(STOECHIOMETRIE_CP[nom])[int(a) - 1]
            else:
                attendu = t.volumes_equivalents()[int(a) - 1]
        elif grandeur.startswith("saut-"):
            s = t.sauts(Vmax)[int(a) - 1]
            attendu = {"volume": s.volume, "pM": getattr(s, "pM", None), "pente": s.pente}[grandeur[5:]]
        else:
            e = t.evenements(Vmax)[int(a) - 1]
            attendu = e.volume if grandeur == "evenement-volume" else _mesure_cp(t, nom, e.volume)
        tolerance = 1e-6 if grandeur.startswith(("saut", "evenement")) else 1e-9
        assert x == pytest.approx(attendu, abs=tolerance, rel=1e-9), (nom, grandeur, argument)
