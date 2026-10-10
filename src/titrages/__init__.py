"""Courbes de titrage sans approximation : pH-métrie, conductimétrie, potentiométrie.

Le pH, la conductivité ou le potentiel d'une solution quelconque, et la courbe de
titrage d'un bécher quelconque par une burette quelconque — acides et bases forts ou
faibles, polyacides, mélanges, couples rédox — avec la dilution. La méthode est celle
de Kalka (2021) et Brière (2002) : à pH (ou à potentiel) imposé, tout est explicite,
et le volume versé est une fonction directe de la grandeur mesurée.

>>> from titrages import Titrage, solution
>>> t = Titrage(solution({"CH3COOH": 0.1}), solution({"NaOH": 0.1}), volume_becher=10)
>>> t.pH([0, 5, 10]).round(2).tolist()
[2.88, 4.76, 8.73]
"""

from titrages.base import (
    Fort,
    chercher,
    conductivite_molaire,
    couple,
    couple_redox,
    couples,
    couples_redox,
    espece,
    espece_redox,
    indicateur,
    indicateurs,
    normaliser,
    solution,
    solution_redox,
)
from titrages.equilibre import (
    LN10,
    PKE,
    concentration_pour_pH,
    delta,
    fractions,
    log_k,
    moments,
)
from titrages.export import vers_csv, vers_pgfplots
from titrages.modeles import (
    Couple,
    CoupleRedox,
    Courbe,
    CourbeConductimetrique,
    CourbeRedox,
    Espece,
    EspeceRedox,
    Indicateur,
    Saut,
    SautRedox,
    Spectateur,
    Tangentes,
    TangentesRedox,
    charge,
)
from titrages.redox import (
    NERNST,
    REFERENCES,
    ComposantRedox,
    SolutionRedox,
    TitrageRedox,
    electrons_cedes,
    fractions_redox,
)
from titrages.solution import LAMBDA_H3O, LAMBDA_HO, Composant, Solution
from titrages.titrage import Titrage

__version__ = "0.1.0"

__all__ = [
    "LAMBDA_H3O",
    "LAMBDA_HO",
    "LN10",
    "NERNST",
    "PKE",
    "REFERENCES",
    "Composant",
    "ComposantRedox",
    "Couple",
    "CoupleRedox",
    "Courbe",
    "CourbeConductimetrique",
    "CourbeRedox",
    "Espece",
    "EspeceRedox",
    "Fort",
    "Indicateur",
    "Saut",
    "SautRedox",
    "Solution",
    "SolutionRedox",
    "Spectateur",
    "Tangentes",
    "TangentesRedox",
    "Titrage",
    "TitrageRedox",
    "__version__",
    "charge",
    "chercher",
    "concentration_pour_pH",
    "conductivite_molaire",
    "couple",
    "couple_redox",
    "couples",
    "couples_redox",
    "delta",
    "electrons_cedes",
    "espece",
    "espece_redox",
    "fractions",
    "fractions_redox",
    "indicateur",
    "indicateurs",
    "log_k",
    "moments",
    "normaliser",
    "solution",
    "solution_redox",
    "vers_csv",
    "vers_pgfplots",
]
