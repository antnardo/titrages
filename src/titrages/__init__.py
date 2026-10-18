"""Courbes de titrage sans approximation : pH-métrie, conductimétrie, potentiométrie,
complexométrie, précipitation.

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
    complexe,
    complexes,
    couple,
    couple_redox,
    couples,
    couples_redox,
    espece,
    espece_redox,
    indicateur,
    indicateurs,
    ligand,
    ligands,
    precipite,
    precipites,
    solution,
    solution_redox,
)
from titrages.complexation import (
    SolutionComplexation,
    TitrageComplexometrique,
)
from titrages.equilibre import (
    concentration_pour_pH,
    fractions,
    fractions_complexes,
    moments,
)
from titrages.export import vers_pgfplots
from titrages.modeles import (
    Complexe,
    Couple,
    CoupleRedox,
    Courbe,
    CourbeComplexometrique,
    CourbeConductimetrique,
    CourbePrecipitation,
    CourbeRedox,
    Espece,
    EspeceRedox,
    Evenement,
    Indicateur,
    Ligand,
    Precipite,
    Saut,
    SautComplexometrique,
    SautPrecipitation,
    SautRedox,
    Spectateur,
    Tangentes,
    TangentesRedox,
)
from titrages.precipitation import SolutionPrecipitation, TitrageParPrecipitation
from titrages.redox import (
    SolutionRedox,
    TitrageRedox,
    fractions_redox,
)
from titrages.solution import Solution
from titrages.titrage import Titrage

__version__ = "0.1.4"

__all__ = [
    "Complexe",
    "Couple",
    "CoupleRedox",
    "Courbe",
    "CourbeComplexometrique",
    "CourbeConductimetrique",
    "CourbePrecipitation",
    "CourbeRedox",
    "Espece",
    "EspeceRedox",
    "Evenement",
    "Indicateur",
    "Ligand",
    "Precipite",
    "Saut",
    "SautComplexometrique",
    "SautPrecipitation",
    "SautRedox",
    "Solution",
    "SolutionComplexation",
    "SolutionPrecipitation",
    "SolutionRedox",
    "Spectateur",
    "Tangentes",
    "TangentesRedox",
    "Titrage",
    "TitrageComplexometrique",
    "TitrageParPrecipitation",
    "TitrageRedox",
    "__version__",
    "chercher",
    "complexe",
    "complexes",
    "concentration_pour_pH",
    "couple",
    "couple_redox",
    "couples",
    "couples_redox",
    "espece",
    "espece_redox",
    "fractions",
    "fractions_complexes",
    "fractions_redox",
    "indicateur",
    "indicateurs",
    "ligand",
    "ligands",
    "moments",
    "precipite",
    "precipites",
    "solution",
    "solution_redox",
    "vers_pgfplots",
]
