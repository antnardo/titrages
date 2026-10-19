"""Le moteur Lua du paquet LaTeX et le moteur Python donnent les mêmes nombres.

Les deux implémentent le même algorithme ; ce test empêche qu'ils divergent. Il
demande ``texlua`` (TeX Live) ou ``lua``, et se passe de LaTeX.
"""

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from titrages import Titrage, solution

LATEX = Path(__file__).parent.parent / "latex"
INTERPRETE = shutil.which("texlua") or shutil.which("lua5.4") or shutil.which("lua5.3") or shutil.which("lua")

pytestmark = [
    pytest.mark.lua,
    pytest.mark.skipif(INTERPRETE is None, reason="ni texlua ni lua"),
]

CAS = {
    "acetique": ({"CH3COOH": 0.1}, {"NaOH": 0.1}, "{T.espece({4.76},0.1,0)}", 0, "{}", 0.1, 20),
    "phosphorique": ({"H3PO4": 0.1}, {"NaOH": 0.1}, "{T.espece({2.15,7.20,12.35},0.1,0)}", 0, "{}", 0.1, 35),
    "ammoniac": ({"NH3": 0.1}, {"HCl": 0.1}, "{T.espece({9.25},0.1,1)}", 0, "{}", -0.1, 20),
    "faible": (
        {"HCl": 0.05, "CH3COOH": 0.05},
        {"NH3": 0.1},
        "{T.espece({4.76},0.05,0)}",
        -0.05,
        "{T.espece({9.25},0.1,1)}",
        0,
        30,
    ),
}
VOLUMES = [0, 2.5, 5, 9.9, 10, 10.1, 15, 19]


def _lua(becher, b0, burette, b1, vmax, dossier):
    script = f"""
package.path = "{LATEX.as_posix()}/?.lua;" .. package.path
local T = require("titrages")
local t = T.titrage(T.solution({becher}, {b0}), T.solution({burette}, {b1}), 10)
local pH = {{}}
for i, V in ipairs({{{", ".join(map(str, VOLUMES))}}}) do pH[i] = string.format("%.17g", T.pH(t, V)) end
local sauts = {{}}
for i, s in ipairs(T.sauts(t, {vmax})) do sauts[i] = string.format("[%.17g,%.17g]", s[1], s[2]) end
io.write('{{"pH":[' .. table.concat(pH, ",") .. '],"sauts":[' .. table.concat(sauts, ",") .. ']}}')
"""
    fichier = dossier / "cas.lua"  # texlua ne lit pas l'entrée standard
    fichier.write_text(script, encoding="utf-8")
    sortie = subprocess.run([INTERPRETE, str(fichier)], capture_output=True, text=True, check=True)
    return json.loads(sortie.stdout)


@pytest.mark.parametrize("nom", CAS)
def test_memes_nombres(nom, tmp_path):
    becher, burette, lb, b0, lu, b1, vmax = CAS[nom]
    t = Titrage(solution(becher), solution(burette), 10)
    lua = _lua(lb, b0, lu, b1, vmax, tmp_path)
    assert lua["pH"] == pytest.approx(t.pH(VOLUMES).tolist(), abs=1e-9)
    sauts = t.sauts(vmax)
    assert len(lua["sauts"]) == len(sauts)
    for (V, pH), s in zip(lua["sauts"], sauts, strict=True):
        assert V == pytest.approx(s.volume, abs=1e-6)
        assert pH == pytest.approx(s.pH, abs=1e-6)


EQUIVALENCES = {
    "phosphorique": ({"H3PO4": 0.1}, {"NaOH": 0.1}, "{T.espece({2.15,7.20,12.35},0.1,0)}", 0, 0.1),
    "carbonate": ({"Na2CO3": 0.1}, {"HCl": 0.1}, "{T.espece({6.35,10.33},0.1,2)}", 0, -0.1),
    "retour acide": ({"CH3COOH": 0.05, "NaOH": 0.1}, {"HCl": 0.1}, "{T.espece({4.76},0.05,0)}", 0.1, -0.1),
    "retour base": ({"NH3": 0.05, "HCl": 0.1}, {"NaOH": 0.1}, "{T.espece({9.25},0.05,1)}", -0.1, 0.1),
    "base déjà présente": (
        {"H2CO3": 0.1, "NaOH": 0.05},
        {"NaOH": 0.1},
        "{T.espece({6.35,10.33},0.1,0)}",
        0.05,
        0.1,
    ),
    "polyacide en retour": (
        {"H3PO4": 0.05, "NaOH": 0.2},
        {"HCl": 0.1},
        "{T.espece({2.15,7.20,12.35},0.05,0)}",
        0.2,
        -0.1,
    ),
}


@pytest.mark.parametrize("nom", EQUIVALENCES)
def test_memes_equivalences(nom, tmp_path):
    becher, burette, lb, b0, b1 = EQUIVALENCES[nom]
    t = Titrage(solution(becher), solution(burette), 10)
    script = f"""
package.path = "{LATEX.as_posix()}/?.lua;" .. package.path
local T = require("titrages")
local t = T.titrage(T.solution({lb}, {b0}), T.solution({{}}, {b1}), 10)
local eq = {{}}
for i, V in ipairs(T.equivalences(t)) do eq[i] = string.format("%.17g", V) end
io.write("[" .. table.concat(eq, ",") .. "]")
"""
    fichier = tmp_path / "equivalences.lua"
    fichier.write_text(script, encoding="utf-8")
    sortie = subprocess.run([INTERPRETE, str(fichier)], capture_output=True, text=True, check=True)
    assert json.loads(sortie.stdout) == pytest.approx(list(t.volumes_equivalents()), abs=1e-9)


REDOX = {
    "fer/cérium": (
        {"Fe2+": 0.1},
        {"Ce4+": 0.1},
        "{T.espece_redox(T.couple_redox({0.77},{1}),0.1,0)}",
        "{T.espece_redox(T.couple_redox({1.72},{1}),0.1,1)}",
    ),
    "fer/dichromate": (
        {"Fe2+": 0.1},
        {"K2Cr2O7": 1 / 60},
        "{T.espece_redox(T.couple_redox({0.77},{1}),0.1,0)}",
        "{T.espece_redox(T.couple_redox({1.33},{6},{14},{1,2}),1/60,1)}",
    ),
    "iode/thiosulfate": (
        {"I2": 0.05},
        {"Na2S2O3": 0.1},
        "{T.espece_redox(T.couple_redox({0.54},{2},{0},{1,2}),0.05,1)}",
        "{T.espece_redox(T.couple_redox({0.08},{2},{0},{1,2}),0.1,0)}",
    ),
}


def lua_redox(nom: str, dossier: Path) -> dict:
    """Le potentiel, les sauts et les tangentes du moteur Lua, pour un cas de REDOX."""
    _, _, lb, lu = REDOX[nom]
    script = f"""
package.path = "{LATEX.as_posix()}/?.lua;" .. package.path
local T = require("titrages")
local t = T.titrage(T.solution_redox({lb}), T.solution_redox({lu}), 10)
local E = {{}}
for i, V in ipairs({{{", ".join(map(str, VOLUMES[1:]))}}}) do E[i] = string.format("%.17g", T.E_enh(t, V)) end
local S = {{}}
for i, s in ipairs(T.sauts_redox(t, 20)) do S[i] = string.format("[%.17g,%.17g,%.17g]", s[1], s[2], s[3]) end
local tg = T.tangentes_redox(t, 7, 20)
io.write('{{"E":[' .. table.concat(E, ",") .. '],"sauts":[' .. table.concat(S, ",") .. '],')
io.write(string.format('"tangentes":[%.17g,%.17g]}}', tg.V, tg.E))
"""
    fichier = dossier / "redox.lua"
    fichier.write_text(script, encoding="utf-8")
    sortie = subprocess.run([INTERPRETE, str(fichier)], capture_output=True, text=True, check=True)
    return json.loads(sortie.stdout)


def python_redox(nom: str):
    from titrages import TitrageRedox, solution_redox

    becher, burette, _, _ = REDOX[nom]
    return TitrageRedox(solution_redox(becher), solution_redox(burette), 10)


@pytest.mark.parametrize("nom", REDOX)
def test_redox_memes_nombres(nom, tmp_path):
    lua, t = lua_redox(nom, tmp_path), python_redox(nom)
    assert lua["E"] == pytest.approx(t.E(VOLUMES[1:]).tolist(), abs=1e-9)
    tangentes = t.tangentes(7, 20)
    assert lua["tangentes"] == pytest.approx([tangentes.volume, tangentes.E], abs=1e-6)
    sauts = t.sauts(20)
    assert len(lua["sauts"]) == len(sauts)
    for (V, E, pente), s in zip(lua["sauts"], sauts, strict=True):
        assert V == pytest.approx(s.volume, rel=1e-9)
        assert E == pytest.approx(s.E, abs=1e-9)
        assert pente == pytest.approx(s.pente, rel=1e-6)


# --- complexométrie et précipitation -------------------------------------------------

CP_VOLUMES = [0.5, 5, 9.99, 10, 10.01, 15, 19.99, 20, 20.01, 24, 30]


def _executer(script: str, dossier: Path) -> dict:
    fichier = dossier / "cp.lua"
    fichier.write_text(
        f'package.path = "{LATEX.as_posix()}/?.lua;" .. package.path\nlocal T = require("titrages")\n'
        + script,
        encoding="utf-8",
    )
    sortie = subprocess.run([INTERPRETE, str(fichier)], capture_output=True, text=True, check=True)
    return json.loads(sortie.stdout)


_EXPORT = """
local function liste(t) local r = {} for i, x in ipairs(t) do r[i] = string.format("%.17g", x) end return "[" .. table.concat(r, ",") .. "]" end
local function triplets(l) local r = {} for i, s in ipairs(l) do r[i] = string.format("[%.17g,%.17g,%.17g]", s[1], s[2], s[3]) end return "[" .. table.concat(r, ",") .. "]" end
"""

COMPLEXATION = {
    "Ca-EDTA": (
        '{{"Ca2+",0.01}}, 0',
        "{}, 0.01",
        '{{"Ca2+",{10.7}}}',
        "{0.0,1.5,2.0,2.66,6.16,10.24}",
        "10",
    ),
    "Ca et Zn": (
        '{{"Ca2+",0.01},{"Zn2+",0.005}}, 0',
        "{}, 0.01",
        '{{"Ca2+",{10.7}},{"Zn2+",{16.5}}}',
        "{0.0,1.5,2.0,2.66,6.16,10.24}",
        "10",
    ),
    "Cu-NH3, métal titrant": (
        "{}, 0.1",
        '{{"Cu2+",0.05}}, 0',
        '{{"Cu2+",{4.1,7.6,10.5,12.6}}}',
        "{9.25}",
        "11",
    ),
}


def _complexation_python(nom: str):
    from titrages import Complexe, Ligand, SolutionComplexation, TitrageComplexometrique

    edta = Ligand("Y4-", (0.0, 1.5, 2.0, 2.66, 6.16, 10.24))
    ca, zn = Complexe("Ca2+", (10.7,)), Complexe("Zn2+", (16.5,))
    if nom == "Ca-EDTA":
        return TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01), 20, [ca], edta, 10
        )
    if nom == "Ca et Zn":
        return TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01, "Zn2+": 0.005}),
            SolutionComplexation(ligand=0.01),
            20,
            [ca, zn],
            edta,
            10,
        )
    return TitrageComplexometrique(
        SolutionComplexation(ligand=0.1),
        SolutionComplexation({"Cu2+": 0.05}),
        20,
        [Complexe("Cu2+", (4.1, 7.6, 10.5, 12.6))],
        Ligand("NH3", (9.25,)),
        11,
    )


def lua_complexation(nom: str, dossier: Path) -> dict:
    """Le moteur Lua sur un cas : ses nombres, en JSON."""
    becher, burette, complexes, pKa, pH = COMPLEXATION[nom]
    return _executer(
        _EXPORT
        + f"""
local t = T.titrage_complexation(T.solution_complexation({becher}), T.solution_complexation({burette}),
  20, {complexes}, {pKa}, {pH})
local pM, pL, pente = {{}}, {{}}, {{}}
for i, V in ipairs({{{", ".join(map(str, CP_VOLUMES))}}}) do
  pM[i], pL[i] = T.pM(t, V), T.pL(t, V)
  pente[i] = T.pente_complexation(t, pL[i])
end
local ok, eq = pcall(T.equivalences_complexation, t)
io.write('{{"pM":' .. liste(pM) .. ',"pente":' .. liste(pente) .. ',"sauts":' .. triplets(T.sauts_complexation(t, 40))
  .. ',"equivalences":' .. (ok and liste(eq) or "null") .. '}}')
""",
        dossier,
    )


@pytest.mark.parametrize("nom", COMPLEXATION)
def test_complexation_memes_nombres(nom, tmp_path):
    lua = lua_complexation(nom, tmp_path)
    t = _complexation_python(nom)
    assert lua["pM"] == pytest.approx(t.pM(CP_VOLUMES).tolist(), abs=1e-9)
    assert lua["pente"] == pytest.approx(t.pente(t.pL(CP_VOLUMES)).tolist(), rel=1e-9)
    sauts = t.sauts(40)
    assert len(lua["sauts"]) == len(sauts)
    for (V, pM, pente), s in zip(lua["sauts"], sauts, strict=True):
        assert (V, pM, pente) == pytest.approx((s.volume, s.pM, s.pente), rel=1e-9)
    assert lua["equivalences"] == pytest.approx(list(t.volumes_equivalents()), rel=1e-12)


PRECIPITATION = {
    "Mohr": (
        '0, {{"Cl-",0.05},{"CrO4^2-",0.005}}',
        "0.05",
        '{{partenaire="Cl-",pKs=9.75},{partenaire="CrO4^2-",pKs=11.95,nu=2,pKa={6.5}}}',
        "7.5",
    ),
    "iodure et chlorure": (
        '0, {{"Cl-",0.05},{"I-",0.05}}',
        "0.1",
        '{{partenaire="Cl-",pKs=9.75},{partenaire="I-",pKs=16.08}}',
        "nil",
    ),
    "Charpentier-Volhard": (
        '0.1, {{"Cl-",0.05}}',
        '0, {{"SCN-",0.1}}',
        '{{partenaire="Cl-",pKs=9.75},{partenaire="SCN-",pKs=12}}',
        "nil",
    ),
}


def _precipitation_python(nom: str):
    from titrages import Precipite, SolutionPrecipitation, TitrageParPrecipitation

    cl = Precipite("Cl-", 9.75)
    if nom == "Mohr":
        return TitrageParPrecipitation(
            SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
            SolutionPrecipitation(commun=0.05),
            20,
            [cl, Precipite("CrO4^2-", 11.95, nu=2, pKa=(6.5,))],
            pH=7.5,
        )
    if nom == "iodure et chlorure":
        return TitrageParPrecipitation(
            SolutionPrecipitation(partenaires={"Cl-": 0.05, "I-": 0.05}),
            SolutionPrecipitation(commun=0.1),
            20,
            [cl, Precipite("I-", 16.08)],
        )
    return TitrageParPrecipitation(
        SolutionPrecipitation(commun=0.1, partenaires={"Cl-": 0.05}),
        SolutionPrecipitation(partenaires={"SCN-": 0.1}),
        20,
        [cl, Precipite("SCN-", 12.0)],
    )


def lua_precipitation(nom: str, dossier: Path) -> dict:
    """Le moteur Lua sur un cas : ses nombres, en JSON."""
    becher, burette, precipites, pH = PRECIPITATION[nom]
    return _executer(
        _EXPORT
        + f"""
local t = T.titrage_precipitation(T.solution_precipitation({becher}), T.solution_precipitation({burette}),
  20, {precipites}, {pH})
local p, pente = {{}}, {{}}
for i, V in ipairs({{{", ".join(map(str, CP_VOLUMES))}}}) do
  p[i] = T.pP(t, V)
  pente[i] = T.pente_precipitation(t, p[i])
end
local ev = {{}}
for i, e in ipairs(T.evenements(t, 30)) do
  ev[i] = string.format('[%.17g,%.17g,"%s",%s]', e[1], e[2], e[3], tostring(e[4]))
end
io.write('{{"p":' .. liste(p) .. ',"pente":' .. liste(pente) .. ',"sauts":' .. triplets(T.sauts_precipitation(t, 30))
  .. ',"evenements":[' .. table.concat(ev, ",") .. '],"equivalences":' .. liste(T.equivalences_precipitation(t)) .. '}}')
""",
        dossier,
    )


@pytest.mark.parametrize("nom", PRECIPITATION)
def test_precipitation_memes_nombres(nom, tmp_path):
    lua = lua_precipitation(nom, tmp_path)
    t = _precipitation_python(nom)
    assert lua["p"] == pytest.approx(t.p(CP_VOLUMES).tolist(), abs=1e-9)
    assert lua["pente"] == pytest.approx(t.pente(t.p(CP_VOLUMES)).tolist(), rel=1e-9)
    sauts = t.sauts(30)
    assert len(lua["sauts"]) == len(sauts)
    for (V, p, pente), s in zip(lua["sauts"], sauts, strict=True):
        assert (V, p, pente) == pytest.approx((s.volume, s.p, s.pente), rel=1e-9)
    evenements = t.evenements(30)
    assert [(x, a) for _, _, x, a in lua["evenements"]] == [(e.partenaire, e.apparition) for e in evenements]
    for (V, p, _, _), e in zip(lua["evenements"], evenements, strict=True):
        assert (V, p) == pytest.approx((e.volume, e.p), rel=1e-9, abs=1e-12)
    assert lua["equivalences"] == pytest.approx(list(t.volumes_equivalents()), rel=1e-12)


# --- le pH calculé (niveau 2) ----------------------------------------------------------

_LIBRE_EXPORT = (
    _EXPORT
    + """
local function etats(t, volumes, mesure)
  local pH, y = {}, {}
  for i, V in ipairs(volumes) do pH[i], y[i] = T.pH_du_melange(t, V), mesure(t, V) end
  return liste(pH), liste(y)
end
"""
)

EDTA_LUA = "{0.0,1.5,2.0,2.66,6.16,10.24}"
TAMPON_LUA = "T.solution({T.espece({9.25},0.1,0), T.espece({9.25},0.1,1)},0,14)"

COMPLEXATION_LIBRE = {
    "Ca par Na2H2Y": (
        'T.solution_complexation({{"Ca2+",0.01}})',
        "T.solution_complexation({}, 0.01, 2)",
        '{{"Ca2+",{10.7}}}',
        EDTA_LUA,
    ),
    "Ca, tampon en solution": (
        f'T.solution_complexation({{{{"Ca2+",0.01}}}}, 0, 0, {TAMPON_LUA})',
        f"T.solution_complexation({{}}, 0.01, 2, {TAMPON_LUA})",
        '{{"Ca2+",{10.7}}}',
        EDTA_LUA,
    ),
    "Cu par NH3": (
        'T.solution_complexation({{"Cu2+",0.01}})',
        "T.solution_complexation({}, 0.1)",
        '{{"Cu2+",{4.04,7.47,10.27,11.75}}}',
        "{9.244}",
    ),
    "Liebig": (
        "T.solution_complexation({}, 0.1)",
        'T.solution_complexation({{"Ag+",0.05}})',
        '{{"Ag+",{-math.huge,20.48,21.40}}}',
        "{9.21}",
    ),
}


def _complexation_libre_python(nom: str):
    import warnings

    from titrages import Complexe, Ligand, Solution, SolutionComplexation, TitrageComplexometrique, couple

    edta = Ligand("Y4-", (0.0, 1.5, 2.0, 2.66, 6.16, 10.24))
    ca = Complexe("Ca2+", (10.7,))
    nh = couple("NH3")
    tampon = Solution().avec(nh.forme(0), 0.1).avec(nh.forme(1), 0.1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        if nom == "Ca par Na2H2Y":
            return TitrageComplexometrique(
                SolutionComplexation({"Ca2+": 0.01}),
                SolutionComplexation(ligand=0.01, protons_ligand=2),
                20,
                [ca],
                edta,
            )
        if nom == "Ca, tampon en solution":
            return TitrageComplexometrique(
                SolutionComplexation({"Ca2+": 0.01}, acidobasique=tampon),
                SolutionComplexation(ligand=0.01, protons_ligand=2, acidobasique=tampon),
                20,
                [ca],
                edta,
            )
        if nom == "Cu par NH3":
            return TitrageComplexometrique(
                SolutionComplexation({"Cu2+": 0.01}),
                SolutionComplexation(ligand=0.1),
                20,
                [Complexe("Cu2+", (4.04, 7.47, 10.27, 11.75))],
                Ligand("NH3", (9.244,)),
            )
        return TitrageComplexometrique(
            SolutionComplexation(ligand=0.1),
            SolutionComplexation({"Ag+": 0.05}),
            20,
            [Complexe("Ag+", (-float("inf"), 20.48, 21.40))],
            Ligand("CN-", (9.21,)),
        )


LIBRE_VOLUMES = [0.5, 5.0, 10.0, 15.0, 19.9, 20.0, 20.1, 25.0, 39.0]


def lua_complexation_libre(nom: str, dossier: Path) -> dict:
    """Le moteur Lua sur un cas : ses nombres, en JSON."""
    becher, burette, complexes, pKa = COMPLEXATION_LIBRE[nom]
    return _executer(
        _LIBRE_EXPORT
        + f"""
local t = T.titrage_complexation({becher}, {burette}, 20, {complexes}, {pKa})
local pH, pM = etats(t, {{{", ".join(map(str, LIBRE_VOLUMES))}}}, function(t, V) return T.pM(t, V) end)
io.write('{{"libre":' .. tostring(t.libre) .. ',"pH":' .. pH .. ',"pM":' .. pM
  .. ',"sauts":' .. triplets(T.sauts_complexation(t, 40)) .. '}}')
""",
        dossier,
    )


@pytest.mark.parametrize("nom", COMPLEXATION_LIBRE)
def test_complexation_libre_memes_nombres(nom, tmp_path):
    lua = lua_complexation_libre(nom, tmp_path)
    t = _complexation_libre_python(nom)
    assert lua["libre"] is True
    assert lua["pH"] == pytest.approx(t.pH_du_melange(LIBRE_VOLUMES).tolist(), abs=1e-9)
    assert lua["pM"] == pytest.approx(t.pM(LIBRE_VOLUMES).tolist(), abs=1e-9)
    sauts = t.sauts(40)
    assert len(lua["sauts"]) == len(sauts)
    for (V, pM, pente), s in zip(lua["sauts"], sauts, strict=True):
        assert (V, pM, pente) == pytest.approx((s.volume, s.pM, s.pente), rel=1e-8)


PRECIPITATION_LIBRE = {
    "Mohr sans tampon": (
        'T.solution_precipitation(0, {{"Cl-",0.05},{"CrO4^2-",0.005}})',
        "T.solution_precipitation(0.05)",
        '{{partenaire="Cl-",pKs=9.74},{partenaire="CrO4^2-",pKs=11.92,nu=2,pKa={-0.2,6.51}}}',
        30,
        [0.5, 5.0, 10.0, 19.9, 20.0, 20.1, 25.0],
    ),
    "Mohr en milieu acide": (
        'T.solution_precipitation(0, {{"Cl-",0.05},{"CrO4^2-",0.005}}, T.solution({}, -0.001, 14))',
        "T.solution_precipitation(0.05)",
        '{{partenaire="Cl-",pKs=9.74},{partenaire="CrO4^2-",pKs=11.92,nu=2,pKa={-0.2,6.51}}}',
        30,
        [0.5, 5.0, 10.0, 19.9, 20.0, 20.1, 25.0],
    ),
    "Volhard en milieu nitrique": (
        'T.solution_precipitation(0.1, {{"Cl-",0.05}}, T.solution({}, -0.1, 14))',
        'T.solution_precipitation(0, {{"SCN-",0.1}})',
        '{{partenaire="Cl-",pKs=9.74},{partenaire="SCN-",pKs=11.97,pKa={0.9}}}',
        25,
        [0.5, 5.0, 9.9, 10.0, 10.1, 15.0, 24.0],
    ),
    "Mg(OH)2 après un acide": (
        "T.solution_precipitation(0.01, {}, T.solution({}, -0.002, 14))",
        "T.solution_precipitation(0, {}, T.solution({}, 0.1, 14))",
        '{{partenaire="OH-",pKs=11.15,nu=1,mu=2}}',
        8,
        [0.1, 0.39, 0.4, 0.41, 2.0, 4.39, 4.4, 4.41, 6.0],
    ),
    "Ag2O": (
        "T.solution_precipitation(0.01)",
        "T.solution_precipitation(0, {}, T.solution({}, 0.1, 14))",
        '{{partenaire="OH-",pKs=15.42,nu=2,mu=2}}',
        4,
        [0.1, 1.0, 1.9, 2.0, 2.1, 3.0],
    ),
}


def _precipitation_libre_python(nom: str):
    import warnings

    from titrages import Precipite, Solution, SolutionPrecipitation, TitrageParPrecipitation

    cl, cro4 = Precipite("Cl-", 9.74), Precipite("CrO4^2-", 11.92, nu=2, pKa=(-0.2, 6.51))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        if nom.startswith("Mohr"):
            acide = Solution(acide_fort=0.001) if "acide" in nom else Solution()
            return TitrageParPrecipitation(
                SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}, acidobasique=acide),
                SolutionPrecipitation(commun=0.05),
                20,
                [cl, cro4],
            )
        if nom.startswith("Volhard"):
            return TitrageParPrecipitation(
                SolutionPrecipitation(
                    commun=0.1, partenaires={"Cl-": 0.05}, acidobasique=Solution(acide_fort=0.1)
                ),
                SolutionPrecipitation(partenaires={"SCN-": 0.1}),
                20,
                [cl, Precipite("SCN-", 11.97, pKa=(0.9,))],
            )
        soude = SolutionPrecipitation(acidobasique=Solution(base_forte=0.1))
        if nom.startswith("Mg"):
            return TitrageParPrecipitation(
                SolutionPrecipitation(commun=0.01, acidobasique=Solution(acide_fort=0.002)),
                soude,
                20,
                [Precipite("OH-", 11.15, mu=2)],
            )
        return TitrageParPrecipitation(
            SolutionPrecipitation(commun=0.01), soude, 20, [Precipite("OH-", 15.42, nu=2, mu=2)]
        )


def lua_precipitation_libre(nom: str, dossier: Path) -> dict:
    """Le moteur Lua sur un cas : ses nombres, en JSON."""
    becher, burette, precipites, Vmax, volumes = PRECIPITATION_LIBRE[nom]
    return _executer(
        _LIBRE_EXPORT
        + f"""
local t = T.titrage_precipitation({becher}, {burette}, 20, {precipites})
local pH, p = etats(t, {{{", ".join(map(str, volumes))}}}, T.pP)
local ev = {{}}
for i, e in ipairs(T.evenements(t, {Vmax})) do
  ev[i] = string.format('[%.17g,%.17g,"%s",%s]', e[1], e[2], e[3], tostring(e[4]))
end
io.write('{{"libre":' .. tostring(t.libre) .. ',"pH":' .. pH .. ',"p":' .. p
  .. ',"sauts":' .. triplets(T.sauts_precipitation(t, {Vmax}))
  .. ',"evenements":[' .. table.concat(ev, ",") .. '],"equivalences":' .. liste(T.equivalences_precipitation(t)) .. '}}')
""",
        dossier,
    )


@pytest.mark.parametrize("nom", PRECIPITATION_LIBRE)
def test_precipitation_libre_memes_nombres(nom, tmp_path):
    _, _, _, Vmax, volumes = PRECIPITATION_LIBRE[nom]
    lua = lua_precipitation_libre(nom, tmp_path)
    t = _precipitation_libre_python(nom)
    assert lua["libre"] is True
    assert lua["pH"] == pytest.approx(t.pH_du_melange(volumes).tolist(), abs=1e-9)
    assert lua["p"] == pytest.approx(t.p(volumes).tolist(), abs=1e-9)
    sauts = t.sauts(Vmax)
    assert len(lua["sauts"]) == len(sauts)
    for (V, p, pente), s in zip(lua["sauts"], sauts, strict=True):
        assert (V, p, pente) == pytest.approx((s.volume, s.p, s.pente), rel=1e-8)
    evenements = t.evenements(Vmax)
    assert [(x, a) for _, _, x, a in lua["evenements"]] == [(e.partenaire, e.apparition) for e in evenements]
    for (V, p, _, _), e in zip(lua["evenements"], evenements, strict=True):
        assert (V, p) == pytest.approx((e.volume, e.p), rel=1e-9, abs=1e-12)
    assert lua["equivalences"] == pytest.approx(list(t.volumes_equivalents()), rel=1e-12)
