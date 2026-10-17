import matplotlib.pyplot as plt
import pytest

from titrages import Titrage, TitrageRedox, couple, indicateur, solution, solution_redox
from titrages.trace import (
    tracer_conductimetrie,
    tracer_derivee,
    tracer_fractions,
    tracer_sauts,
    tracer_tangentes,
    tracer_titrage,
    tracer_titrage_redox,
    tracer_zone_virage,
)


def test_tout_se_trace(acetique):
    fig, ax = plt.subplots()
    tracer_titrage(acetique, 20, ax)
    tracer_sauts(acetique, 20, ax)
    tracer_tangentes(acetique, 7, 20, ax)
    tracer_zone_virage(indicateur("phénolphtaléine"), ax)
    tracer_derivee(acetique, 20, ax.twinx())
    assert len(ax.lines) >= 6
    plt.close(fig)


def test_abscisse_normalisee(acetique):
    fig, ax = plt.subplots()
    tracer_titrage(acetique, 20, ax, abscisse=10)
    assert ax.lines[0].get_xdata().max() == 2
    plt.close(fig)


def test_abscisse_normalisee_partout(acetique):
    """Sauts, tangentes et dérivée se superposent à la courbe tracée en V/V_éq."""
    fig, ax = plt.subplots()
    tracer_titrage(acetique, 20, ax, abscisse=10)
    tracer_sauts(acetique, 20, ax, abscisse=10)
    tracer_tangentes(acetique, 7, 20, ax, abscisse=10)
    tracer_derivee(acetique, 20, ax.twinx(), abscisse=10)
    for ligne in ax.lines + ax.figure.axes[1].lines:
        assert ligne.get_xdata().max() <= 2.01
    plt.close(fig)


def test_fractions():
    fig, ax = plt.subplots()
    tracer_fractions(couple("H3PO4"), ax)
    assert [line.get_label() for line in ax.lines] == ["H3PO4", "H2PO4-", "HPO4^2-", "PO4^3-"]
    tracer_fractions([4.76], ax, legende=False)
    plt.close(fig)


def test_conductimetrie_et_redox():
    fig, (a, b) = plt.subplots(1, 2)
    tracer_conductimetrie(
        Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), 100), 20, a, corrigee=True
    )
    t = TitrageRedox(solution_redox({"Fe2+": 0.1}), solution_redox({"Ce4+": 0.1}), 10, "ECS")
    tracer_titrage_redox(t, 20, b, sauts=True)
    tracer_sauts(t, 20, b)
    tracer_tangentes(t, 7, 20, b)
    tracer_derivee(t, 20, b.twinx())
    assert len(b.lines) >= 8
    assert "ECS" in b.get_ylabel() and "S/m" in a.get_ylabel()
    tracer_conductimetrie(
        Titrage(solution({"HCl": 0.01}), solution({"NaOH": 0.1}), 100), 20, a, unite="mS/cm"
    )
    assert a.lines[-1].get_ydata()[0] == pytest.approx(10 * a.lines[0].get_ydata()[0] / (1 + 0))
    plt.close(fig)


def test_complexometrie_et_precipitation():
    from titrages import (
        Complexe,
        Ligand,
        Precipite,
        SolutionComplexation,
        SolutionPrecipitation,
        TitrageComplexometrique,
        TitrageParPrecipitation,
    )
    from titrages.trace import tracer_titrage_complexometrique, tracer_titrage_precipitation

    edta = Ligand("Y4-", (0.0, 1.5, 2.0, 2.66, 6.16, 10.24))
    c = TitrageComplexometrique(
        SolutionComplexation({"Ca2+": 0.01}),
        SolutionComplexation(ligand=0.01),
        20,
        [Complexe("Ca2+", (10.7,))],
        edta,
        10,
    )
    p = TitrageParPrecipitation(
        SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
        SolutionPrecipitation(commun=0.05),
        20,
        [Precipite("Cl-", 9.75), Precipite("CrO4^2-", 11.95, nu=2, pKa=(6.5,))],
        pH=7.5,
    )
    fig, (a, b) = plt.subplots(1, 2)
    tracer_titrage_complexometrique(c, 40, a, sauts=True)
    tracer_sauts(c, 40, a)
    tracer_derivee(c, 40, a.twinx())
    tracer_titrage_precipitation(p, 30, b, abscisse=20)
    tracer_sauts(p, 30, b, abscisse=20)
    assert a.get_ylabel() == "pCa2+" and b.get_ylabel() == "pAg+"
    # la courbe, puis un losange par événement : Cl- au départ, CrO4^2- après l'équivalence
    assert len(b.lines) >= 3
    assert all(ligne.get_xdata().max() <= 1.5 + 1e-9 for ligne in b.lines)
    plt.close(fig)
