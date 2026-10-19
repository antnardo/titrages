"""Les tables de vérification du document technique, recalculées.

    uv run python docs/technique/verifications.py

Chaque table confronte le paquet à un calcul qui ne lui emprunte rien — les références
des tests, en double précision par force brute ou en précision arbitraire par mpmath —
et en donne le plus grand écart mesuré, et non la tolérance des tests. Les tables
s'écrivent dans ``docs/technique/tables/``, que ``titrages-technique.tex`` inclut.
Il faut les dépendances de développement (``uv sync``) et, pour la dernière, ``texlua``.
"""

import re
import sys
import tempfile
import warnings
from collections.abc import Callable, Iterable
from pathlib import Path

import numpy as np

RACINE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RACINE / "tests"))
sys.path.insert(0, str(RACINE / "src"))

from mpmath import mp  # noqa: E402

import reference_mpmath as ref  # noqa: E402
import test_complexation as tc  # noqa: E402
import test_complexation_libre as tcl  # noqa: E402
import test_dilution as td  # noqa: E402
import test_lua as tl  # noqa: E402
import test_precipitation as tp  # noqa: E402
import test_precipitation_libre as tpl  # noqa: E402
from titrages import Titrage, TitrageRedox, solution, solution_redox  # noqa: E402

TABLES = Path(__file__).parent / "tables"
warnings.simplefilter("ignore", UserWarning)


# Les noms des cas, tels que les tests les écrivent, et leur version composée.
NOMS_LUA = {
    "acetique": "CH3COOH par NaOH",
    "phosphorique": "H3PO4 par NaOH",
    "ammoniac": "NH3 par HCl",
    "faible": "HCl et CH3COOH par NH3",
    "fer/cérium": "Fe2+ par Ce4+",
    "fer/dichromate": "Fe2+ par Cr2O7^2-",
    "iode/thiosulfate": "I2 par S2O3^2-",
}
FORMULES = re.compile(
    r"(?<![\\w])(Cr2O7(?: |\^)2-|S2O3\^2-|CrO4\^2-|Mg\(OH\)2|Na2H2Y|CH3COOH|H3PO4|NaOH|HCl|NH3|HCN|"
    r"Ag2O|Fe2\+|Ce4\+|MnO4-|I2|Ca|Zn|Cu|Ag|Fe)(?![\w])"
)


def libelle(nom: str) -> str:
    """Le nom d'un cas, ses formules chimiques composées par mhchem."""
    nom = NOMS_LUA.get(nom, nom)
    return FORMULES.sub(lambda m: "\\ce{" + m.group(1).replace(" 2-", "^2-") + "}", nom)


def nombre(x: float) -> str:
    """Un écart, en notation scientifique pour siunitx ; 0 quand il est nul."""
    return "0" if x == 0 else f"\\num{{{x:.1e}}}"


def ecrire(nom: str, entete: str, colonnes: str, lignes: Iterable[Iterable[str]]) -> None:
    corps = "\n".join(" & ".join(libelle(cellule) for cellule in ligne) + r" \\" for ligne in lignes)
    texte = (
        f"% écrit par docs/technique/verifications.py : ne pas éditer\n"
        f"\\begin{{tabular}}{{{colonnes}}}\n\\toprule\n{entete} \\\\\n\\midrule\n"
        f"{corps}\n\\bottomrule\n\\end{{tabular}}\n"
    )
    (TABLES / f"{nom}.tex").write_text(texte, encoding="utf-8")
    print(f"tables/{nom}.tex")


def ecart_max(f: Callable[[float], float], g: Callable[[float], float], volumes: Iterable[float]) -> float:
    return max(abs(f(V) - g(V)) for V in volumes)


# --- acide-base et rédox : la dilution -----------------------------------------------------


def dilution() -> None:
    volumes = [0.0, 0.5, 2.5, 5.0, 9.0, 9.99, 10.0, 10.01, 11.0, 15.0, 20.0, 29.0]
    lignes = []
    for nom, (becher, burette) in td.CAS.items():
        t = Titrage(td.construire(becher), td.construire(burette), td.V0)

        def brute(V: float, becher: dict = becher, burette: dict = burette) -> float:
            return td.ph_force_brute(*td.melange(becher, burette, td.V0, V))

        lignes.append(
            [nom, str(len(volumes)), nombre(ecart_max(lambda V, t=t: float(t.pH(V)), brute, volumes))]
        )
    V = [1.0, 5.0, 9.5, 10.5, 15.0, 19.0]
    fer, ce, cr = td.FER, td.CERIUM, td.DICHROMATE
    t = TitrageRedox(td.SolutionRedox().avec(fer.red, 0.1), td.SolutionRedox().avec(ce.ox, 0.1), td.V0)

    def brute_ce(v: float) -> float:
        fb, fu = td.V0 / (td.V0 + v), v / (td.V0 + v)
        return td.e_force_brute([(0.77, 1, (1, 1), 0.1 * fb, 0.0), (1.72, 1, (1, 1), 0.1 * fu, 0.1 * fu)])

    lignes.append(
        ["Fe2+ par Ce4+ ($E$)", str(len(V)), nombre(ecart_max(lambda v: float(t.E(v)), brute_ce, V))]
    )
    C = 1 / 60
    t2 = TitrageRedox(td.SolutionRedox().avec(fer.red, 0.1), td.SolutionRedox().avec(cr.ox, C), td.V0)

    def brute_cr(v: float) -> float:
        fb, fu = td.V0 / (td.V0 + v), v / (td.V0 + v)
        T = 2 * C * fu
        return td.e_force_brute([(0.77, 1, (1, 1), 0.1 * fb, 0.0), (1.33, 6, (1, 2), T, 3 * T)])

    lignes.append(
        [
            "Fe2+ par Cr2O7^2-, dimère ($E$)",
            str(len(V)),
            nombre(ecart_max(lambda v: float(t2.E(v)), brute_cr, V)),
        ]
    )
    ecrire("dilution", "Titrage & volumes & écart maximal (pH ou V)", "lrr", lignes)


# --- les points d'inflexion, face à mpmath à 60 chiffres -------------------------------------


def inflexions() -> None:
    titrages = {
        "CH3COOH par NaOH": Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), 10),
        "HCl par NaOH": Titrage(solution({"HCl": 0.1}), solution({"NaOH": 0.1}), 10),
        "NH3 par HCl": Titrage(solution({"NH3": 0.1}), solution({"HCl": 0.1}), 10),
    }
    lignes = []
    for nom, (becher, burette, a, b) in ref.ACIDE_BASE.items():
        attendu = ref.inflexion(ref.volume_acide_base(ref.residu(*becher), ref.residu(*burette)), a, b)
        (s,) = titrages[nom].sauts(20)
        lignes.append(
            [
                nom,
                "\\num{" + mp.nstr(attendu, 13) + "}",
                f"\\num{{{s.pH:.13f}}}",
                nombre(abs(s.pH - float(attendu))),
            ]
        )
    redox = {
        "Fe2+ par Ce4+": {"Ce4+": 0.1},
        "Fe2+ par MnO4-": {"KMnO4": 0.02},
        "Fe2+ par Cr2O7 2-": {"K2Cr2O7": 1 / 60},
    }
    for nom, (V, a, b) in ref.REDOX.items():
        attendu = ref.inflexion(V, a, b, pas=0.005)
        (s,) = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox(redox[nom]), 10).sauts(20)
        lignes.append(
            [
                nom + " (V/ENH)",
                "\\num{" + mp.nstr(attendu, 13) + "}",
                f"\\num{{{s.E:.13f}}}",
                nombre(abs(s.E - float(attendu))),
            ]
        )
    ecrire("inflexions", "Titrage & mpmath & paquet & écart", "lrrr", lignes)


# --- complexométrie et précipitation à pH imposé ----------------------------------------------


def complexation() -> None:
    volumes = [0.5, 5.0, 9.99, 10.0, 15.0, 19.9, 20.0, 20.1, 25.0, 39.0]
    lignes = []
    for nom, cas in tc.CAS.items():
        becher, burette, V0, complexes, ligand, pH, metal = cas
        t = tc._titrage(cas)
        ecarts = []
        for V in volumes:
            if burette.metaux and V == 0.5:  # le métal titrant à l'état de traces, comme dans les tests
                continue
            attendu = tc.reference_pM(
                becher, burette, V0, complexes, ligand, pH if pH is not None else 7.0, V, metal
            )
            ecarts.append(abs(float(t.pM(V)) - attendu))
        lignes.append([nom, str(len(ecarts)), nombre(max(ecarts))])
    ecrire("complexation", "Titrage & volumes & écart maximal sur pM", "lrr", lignes)


def precipitation() -> None:
    volumes = [0.5, 5.0, 9.99, 10.0, 10.01, 15.0, 19.99, 20.0, 20.01, 24.0, 30.0]
    lignes = [
        [
            nom,
            str(len(volumes)),
            nombre(ecart_max(lambda V, t=t: float(t.p(V)), lambda V, t=t: tp.reference_p(t, V), volumes)),
        ]
        for nom, t in tp.CAS.items()
    ]
    ecrire("precipitation", "Titrage & volumes & écart maximal sur pP", "lrr", lignes)


# --- le pH calculé ------------------------------------------------------------------------------


def complexation_libre() -> None:
    lignes = []
    for nom, (cas, _, metal) in tcl.CAS.items():
        t = tcl._titrage(nom)
        ecarts_pH, ecarts_pM = [], []
        for V in tcl.VOLUMES:
            pH, pM = tcl._reference(cas, V, metal)
            ecarts_pH.append(abs(float(t.pH_du_melange(V)) - pH))
            ecarts_pM.append(abs(float(t.pM(V)) - pM))
        lignes.append([nom, str(len(tcl.VOLUMES)), nombre(max(ecarts_pH)), nombre(max(ecarts_pM))])
    ecrire("complexation_libre", "Titrage & volumes & écart sur pH & écart sur pM", "lrrr", lignes)


def precipitation_libre() -> None:
    lignes = []
    for nom, (cas, _, volumes) in tpl.CAS.items():
        t = tpl._titrage(nom)
        ecarts_pH, ecarts_p = [], []
        for V in volumes:
            pH, p = tpl._reference(cas, V)
            ecarts_pH.append(abs(float(t.pH_du_melange(V)) - pH))
            ecarts_p.append(abs(float(t.p(V)) - p))
        lignes.append([nom, str(len(volumes)), nombre(max(ecarts_pH)), nombre(max(ecarts_p))])
    ecrire("precipitation_libre", "Titrage & volumes & écart sur pH & écart sur pP", "lrrr", lignes)


# --- Lua face à Python ------------------------------------------------------------------------


def ecart_listes(a: Iterable[float], b: Iterable[float]) -> float:
    return float(np.max(np.abs(np.asarray(list(a), dtype=float) - np.asarray(list(b), dtype=float))))


def ecart_sauts(lua: list, python: list, grandeur: str) -> tuple[float, float]:
    """Le plus grand écart relatif sur le volume, et absolu sur la grandeur, des sauts."""
    if not python:
        return 0.0, 0.0
    dV = max(abs(x[0] - s.volume) / s.volume for x, s in zip(lua, python, strict=True))
    dy = max(abs(x[1] - getattr(s, grandeur)) for x, s in zip(lua, python, strict=True))
    return dV, dy


def lua() -> None:
    lignes = []
    with tempfile.TemporaryDirectory() as dossier:
        d = Path(dossier)
        for nom, (becher, burette, lb, b0, lu, b1, vmax) in tl.CAS.items():
            r = tl._lua(lb, b0, lu, b1, vmax, d)
            t = Titrage(solution(becher), solution(burette), 10)
            sauts = [[V, pH] for V, pH in r["sauts"]]
            dV, dy = ecart_sauts(sauts, list(t.sauts(vmax)), "pH")
            lignes.append(
                ["acide-base", nom, nombre(ecart_listes(r["pH"], t.pH(tl.VOLUMES))), nombre(dV), nombre(dy)]
            )
        for nom in tl.REDOX:
            r, t = tl.lua_redox(nom, d), tl.python_redox(nom)
            dV, dy = ecart_sauts(r["sauts"], list(t.sauts(20)), "E")
            lignes.append(
                ["rédox", nom, nombre(ecart_listes(r["E"], t.E(tl.VOLUMES[1:]))), nombre(dV), nombre(dy)]
            )
        for nom in tl.COMPLEXATION:
            r, t = tl.lua_complexation(nom, d), tl._complexation_python(nom)
            dV, dy = ecart_sauts(r["sauts"], list(t.sauts(40)), "pM")
            lignes.append(
                [
                    "complexation",
                    nom,
                    nombre(ecart_listes(r["pM"], t.pM(tl.CP_VOLUMES))),
                    nombre(dV),
                    nombre(dy),
                ]
            )
        for nom in tl.PRECIPITATION:
            r, t = tl.lua_precipitation(nom, d), tl._precipitation_python(nom)
            dV, dy = ecart_sauts(r["sauts"], list(t.sauts(30)), "p")
            lignes.append(
                [
                    "précipitation",
                    nom,
                    nombre(ecart_listes(r["p"], t.p(tl.CP_VOLUMES))),
                    nombre(dV),
                    nombre(dy),
                ]
            )
        for nom in tl.COMPLEXATION_LIBRE:
            r, t = tl.lua_complexation_libre(nom, d), tl._complexation_libre_python(nom)
            dV, dy = ecart_sauts(r["sauts"], list(t.sauts(40)), "pM")
            ecart = max(
                ecart_listes(r["pH"], t.pH_du_melange(tl.LIBRE_VOLUMES)),
                ecart_listes(r["pM"], t.pM(tl.LIBRE_VOLUMES)),
            )
            lignes.append(["complexation, pH calculé", nom, nombre(ecart), nombre(dV), nombre(dy)])
        for nom, (_, _, _, vmax, volumes) in tl.PRECIPITATION_LIBRE.items():
            r, t = tl.lua_precipitation_libre(nom, d), tl._precipitation_libre_python(nom)
            dV, dy = ecart_sauts(r["sauts"], list(t.sauts(vmax)), "p")
            ecart = max(ecart_listes(r["pH"], t.pH_du_melange(volumes)), ecart_listes(r["p"], t.p(volumes)))
            lignes.append(["précipitation, pH calculé", nom, nombre(ecart), nombre(dV), nombre(dy)])
    ecrire(
        "lua",
        "Famille & Titrage & écart sur la courbe & sauts : $V$ (relatif) & sauts : grandeur",
        "llrrr",
        lignes,
    )


if __name__ == "__main__":
    TABLES.mkdir(exist_ok=True)
    for etape in (
        dilution,
        inflexions,
        complexation,
        precipitation,
        complexation_libre,
        precipitation_libre,
        lua,
    ):
        etape()
