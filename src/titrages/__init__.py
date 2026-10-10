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
    chercher,
    couple,
    couple_redox,
    couples,
    couples_redox,
    espece,
    espece_redox,
    indicateur,
    indicateurs,
    solution,
    solution_redox,
)
from titrages.equilibre import (
    concentration_pour_pH,
    fractions,
    moments,
)
from titrages.export import vers_pgfplots
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
)
from titrages.redox import (
    SolutionRedox,
    TitrageRedox,
    fractions_redox,
)
from titrages.solution import Solution
from titrages.titrage import Titrage

__version__ = "0.1.1"

__all__ = [
    "Couple",
    "CoupleRedox",
    "Courbe",
    "CourbeConductimetrique",
    "CourbeRedox",
    "Espece",
    "EspeceRedox",
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
    "chercher",
    "concentration_pour_pH",
    "couple",
    "couple_redox",
    "couples",
    "couples_redox",
    "espece",
    "espece_redox",
    "fractions",
    "fractions_redox",
    "indicateur",
    "indicateurs",
    "moments",
    "solution",
    "solution_redox",
    "vers_pgfplots",
]
