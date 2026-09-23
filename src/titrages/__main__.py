"""``python -m titrages`` : la base, le pH, la conductivité, le potentiel, un titrage.

.. code-block:: console

    $ python -m titrages pH CH3COOH=0.1 CH3COO-=0.05
    $ python -m titrages sigma HCl=0.01
    $ python -m titrages E Fe2+=0.1 Fe3+=0.01
    $ python -m titrages titrage "H3PO4=0.1" --burette NaOH=0.1 --volume 10
    $ python -m titrages redox Fe2+=0.1 --burette KMnO4=0.02 --reference ECS
    $ python -m titrages couples
"""

import argparse
import sys
from collections.abc import Sequence

from titrages.base import chercher, couple, couples, couples_redox, indicateurs, solution, solution_redox
from titrages.redox import TitrageRedox
from titrages.titrage import Titrage

__all__ = ["main"]


def _composition(termes: Sequence[str]) -> dict[str, float]:
    resultat: dict[str, float] = {}
    for terme in termes:
        nom, _, valeur = terme.rpartition("=")
        if not nom:
            raise SystemExit(f"« {terme} » : écrire espèce=concentration, par exemple NH3=0.1")
        resultat[nom.strip()] = float(valeur)
    return resultat


def _couples(motif: str | None) -> None:
    retenus = _groupe(motif) if motif else None
    for c in couples():
        if retenus is not None and c.nom not in retenus:
            continue
        pKa = ", ".join(f"{p:g}" for p in c.pKa)
        print(f"{c.nom:32s} pKa {pKa:22s} {' / '.join(c.formes)}")


def _groupe(motif: str) -> set[str]:
    """Les couples dont un nom ou une forme contient ``motif``."""
    noms = set()
    for cle in chercher(motif):
        try:
            noms.add(couple(cle).nom)
        except KeyError:  # un acide ou une base forts
            continue
    return noms


def _titrage(args: argparse.Namespace) -> None:
    t = Titrage(solution(_composition(args.becher)), solution(_composition(args.burette)), args.volume)
    print(f"pH initial : {t.pH_initial:.2f}")
    try:
        equivalences = t.volumes_equivalents()
        print("équivalences stœchiométriques : " + ", ".join(f"{v:.4g}" for v in equivalences))
        volume_max = args.volume_max or 1.5 * equivalences[-1]
    except ValueError:
        volume_max = args.volume_max or 3 * args.volume
    for s in t.sauts(volume_max):
        print(f"saut : V = {s.volume:.4g}, pH = {s.pH:.2f}, pente = {s.pente:.3g}")


def _redox(args: argparse.Namespace) -> None:
    t = TitrageRedox(
        solution_redox(_composition(args.becher), pH=args.pH),
        solution_redox(_composition(args.burette), pH=args.pH),
        args.volume,
        args.reference,
    )
    equivalences = t.volumes_equivalents()
    print("équivalences stœchiométriques : " + ", ".join(f"{v:.4g}" for v in equivalences))
    for s in t.sauts(args.volume_max or 1.5 * equivalences[-1]):
        print(f"saut : V = {s.volume:.4g}, E = {s.E:.3f} V/{args.reference}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="titrages", description=__doc__.splitlines()[0])
    sous = parser.add_subparsers(dest="commande", required=True)

    p_couples = sous.add_parser("couples", help="les couples de la base")
    p_couples.add_argument("motif", nargs="?", help="ne garder que ceux qui contiennent ce texte")

    sous.add_parser("indicateurs", help="les indicateurs colorés de la base")

    p_pH = sous.add_parser("pH", help="le pH d'une solution")
    p_pH.add_argument("composition", nargs="+", metavar="espèce=C")

    p_sigma = sous.add_parser("sigma", help="la conductivité d'une solution (mS/cm)")
    p_sigma.add_argument("composition", nargs="+", metavar="espèce=C")

    p_E = sous.add_parser("E", help="le potentiel d'une solution rédox (V/ENH)")
    p_E.add_argument("composition", nargs="+", metavar="espèce=C")
    p_E.add_argument("--pH", type=float, default=0.0)

    sous.add_parser("redox-couples", help="les couples rédox de la base")

    p_r = sous.add_parser("redox", help="les équivalences d'un titrage potentiométrique")
    p_r.add_argument("becher", nargs="+", metavar="espèce=C")
    p_r.add_argument("--burette", nargs="+", required=True, metavar="espèce=C")
    p_r.add_argument("--volume", type=float, default=10.0)
    p_r.add_argument("--volume-max", type=float)
    p_r.add_argument("--pH", type=float, default=0.0, help="pH du milieu, supposé tamponné")
    p_r.add_argument("--reference", default="ENH", help="ENH, ECS ou Ag/AgCl")

    p_t = sous.add_parser("titrage", help="les équivalences d'un titrage")
    p_t.add_argument("becher", nargs="+", metavar="espèce=C")
    p_t.add_argument("--burette", nargs="+", required=True, metavar="espèce=C")
    p_t.add_argument("--volume", type=float, default=10.0, help="volume du bécher (10 par défaut)")
    p_t.add_argument("--volume-max", type=float, help="volume versé maximal")

    args = parser.parse_args(argv)
    try:
        if args.commande == "couples":
            _couples(args.motif)
        elif args.commande == "indicateurs":
            for i in indicateurs():
                print(
                    f"{i.nom:22s} pKa {i.pKa:<5g} virage {i.zone[0]:g}–{i.zone[1]:g}  "
                    f"{i.couleur_acide} → {i.couleur_base}"
                )
        elif args.commande == "pH":
            print(f"{solution(_composition(args.composition)).pH():.2f}")
        elif args.commande == "sigma":
            print(f"{10 * float(solution(_composition(args.composition)).conductivite()):.4g}")
        elif args.commande == "E":
            print(f"{solution_redox(_composition(args.composition), pH=args.pH).E():.3f}")
        elif args.commande == "redox-couples":
            for c in couples_redox():
                E0 = ", ".join(f"{e:g}" for e in c.E0)
                print(f"{c.nom:22s} E° {E0:18s} {' / '.join(c.formes)}")
        elif args.commande == "redox":
            _redox(args)
        else:
            _titrage(args)
    except (KeyError, ValueError) as erreur:
        print(erreur.args[0], file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
