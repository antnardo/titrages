"""Construire l'archive CTAN du paquet LaTeX, sans rien envoyer.

    uv run python scripts/ctan.py
    uv run python scripts/ctan.py --sans-compilation   # le PDF est déjà là

Compile la documentation (lualatex, makeindex, puis deux passes) — sauf avec
``--sans-compilation``, où ``latex/titrages-doc.pdf`` doit exister —, vérifie que
les versions du .sty, du .lua et du paquet Python concordent, puis écrit
``dist/titrages-ctan.zip`` au format que CTAN attend : un dossier ``titrages/``
avec les sources, la documentation et son PDF, le README et la licence. L'envoi à
CTAN reste un geste manuel ; la publication GitHub joint l'archive à la release.
"""

import re
import subprocess
import sys
import zipfile
from pathlib import Path

RACINE = Path(__file__).resolve().parent.parent
LATEX = RACINE / "latex"
FICHIERS = ["titrages.sty", "titrages.lua", "titrages-doc.tex", "titrages-doc.pdf", "README.md"]


def versions() -> dict[str, str]:
    lire = lambda chemin, motif: re.search(motif, (RACINE / chemin).read_text(encoding="utf-8")).group(1)  # noqa: E731
    return {
        "python": lire("src/titrages/__init__.py", r'__version__ = "([^"]+)"'),
        "sty": lire("latex/titrages.sty", r"\{titrages\}\{[^}]*\}\{([^}]+)\}"),
        "lua": lire("latex/titrages.lua", r'M\.version = "([^"]+)"'),
    }


def compiler() -> None:
    lualatex = ["lualatex", "-interaction=nonstopmode", "-halt-on-error", "titrages-doc.tex"]
    index = ["makeindex", "-q", "-s", "gind.ist", "-o", "titrages-doc.ind", "titrages-doc.idx"]
    for commande in (lualatex, index, lualatex, lualatex):
        subprocess.run(commande, cwd=LATEX, check=True, stdout=subprocess.DEVNULL)
    for ext in ("aux", "log", "out", "toc", "listing", "idx", "ind", "ilg", "glo", "hd"):
        (LATEX / f"titrages-doc.{ext}").unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    v = versions()
    if len(set(v.values())) != 1:
        print(f"versions discordantes : {v}", file=sys.stderr)
        return 1
    if "--sans-compilation" in (sys.argv[1:] if argv is None else argv):
        if not (LATEX / "titrages-doc.pdf").exists():
            print("latex/titrages-doc.pdf manque : compiler, ou retirer --sans-compilation", file=sys.stderr)
            return 1
    else:
        compiler()
    dist = RACINE / "dist"
    dist.mkdir(exist_ok=True)
    archive = dist / "titrages-ctan.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for nom in FICHIERS:
            z.write(LATEX / nom, f"titrages/{nom}")
        z.write(RACINE / "LICENSE", "titrages/LICENSE")
    print(f"{archive} — version {v['sty']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
