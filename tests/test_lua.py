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
for i, V in ipairs({{{", ".join(map(str, VOLUMES))}}}) do pH[i] = string.format("%.12g", T.pH(t, V)) end
local sauts = {{}}
for i, s in ipairs(T.sauts(t, {vmax})) do sauts[i] = string.format("[%.12g,%.12g]", s[1], s[2]) end
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
for i, V in ipairs(T.equivalences(t)) do eq[i] = string.format("%.12g", V) end
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


@pytest.mark.parametrize("nom", REDOX)
def test_redox_memes_nombres(nom, tmp_path):
    from titrages import TitrageRedox, solution_redox

    becher, burette, lb, lu = REDOX[nom]
    t = TitrageRedox(solution_redox(becher), solution_redox(burette), 10)
    script = f"""
package.path = "{LATEX.as_posix()}/?.lua;" .. package.path
local T = require("titrages")
local t = T.titrage(T.solution_redox({lb}), T.solution_redox({lu}), 10)
local E = {{}}
for i, V in ipairs({{{", ".join(map(str, VOLUMES[1:]))}}}) do E[i] = string.format("%.12g", T.E_enh(t, V)) end
local S = {{}}
for i, s in ipairs(T.sauts_redox(t, 20)) do S[i] = string.format("[%.15g,%.15g,%.12g]", s[1], s[2], s[3]) end
io.write('{{"E":[' .. table.concat(E, ",") .. '],"sauts":[' .. table.concat(S, ",") .. ']}}')
"""
    fichier = tmp_path / "redox.lua"
    fichier.write_text(script, encoding="utf-8")
    sortie = subprocess.run([INTERPRETE, str(fichier)], capture_output=True, text=True, check=True)
    lua = json.loads(sortie.stdout)
    assert lua["E"] == pytest.approx(t.E(VOLUMES[1:]).tolist(), abs=1e-9)
    sauts = t.sauts(20)
    assert len(lua["sauts"]) == len(sauts)
    for (V, E, pente), s in zip(lua["sauts"], sauts, strict=True):
        assert V == pytest.approx(s.volume, rel=1e-9)
        assert E == pytest.approx(s.E, abs=1e-9)
        assert pente == pytest.approx(s.pente, rel=1e-6)
