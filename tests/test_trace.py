import matplotlib.pyplot as plt

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
    assert "ECS" in b.get_ylabel() and "mS/cm" in a.get_ylabel()
    plt.close(fig)
