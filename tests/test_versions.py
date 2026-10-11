"""Le paquet Python, le paquet LaTeX et le moteur Lua partagent un numéro de version :
``scripts/ctan.py`` refuse de construire l'archive CTAN sinon, et la publication
d'une release échoue. Ce test le dit avant."""

import re
from pathlib import Path

RACINE = Path(__file__).resolve().parents[1]


def _version(chemin: str, motif: str) -> str:
    texte = (RACINE / chemin).read_text(encoding="utf-8")
    trouve = re.search(motif, texte)
    assert trouve is not None, f"version introuvable dans {chemin}"
    return trouve.group(1)


def test_python_sty_lua_et_doc_ont_la_meme_version() -> None:
    versions = {
        "python": _version("src/titrages/__init__.py", r'__version__ = "([^"]+)"'),
        "sty": _version("latex/titrages.sty", r"\{titrages\}\{[^}]*\}\{([^}]+)\}"),
        "lua": _version("latex/titrages.lua", r'M\.version = "([^"]+)"'),
        "doc": _version("latex/titrages-doc.tex", r"\\date\{Version ([^ ]+) "),
    }
    assert len(set(versions.values())) == 1, versions
