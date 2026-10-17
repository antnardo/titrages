"""La complexométrie, confrontée à un calcul indépendant.

La référence n'emprunte rien au paquet : elle dilue les concentrations à la main, au
volume versé, et résout les bilans du mélange — celui de chaque métal et celui du ligand
— par dichotomie en multiprécision, sur la concentration de ligand libre. Le paquet, lui,
part du ligand libre et en déduit le volume : les deux chemins ne se croisent qu'au
résultat.
"""

import math

import numpy as np
import pytest
from mpmath import mp, mpf

from titrages import Complexe, Ligand, SolutionComplexation, TitrageComplexometrique, fractions_complexes

mp.dps = 50

EDTA = Ligand("Y4-", (0.0, 1.5, 2.0, 2.66, 6.16, 10.24))
NH3 = Ligand("NH3", (9.25,))
CA = Complexe("Ca2+", (10.7,))
ZN = Complexe("Zn2+", (16.5,))
CU = Complexe("Cu2+", (4.1, 7.6, 10.5, 12.6))


def _alpha(ligand: Ligand, pH: float) -> mpf:
    """La fraction de la forme la plus basique, par la formule des manuels."""
    h = mp.power(10, -mpf(pH))
    somme, terme = mpf(1), mpf(1)
    for pKa in sorted(ligand.pKa, reverse=True):
        terme = terme * h / mp.power(10, -mpf(pKa))
        somme += terme
    return 1 / somme


def reference_pM(becher, burette, V0, complexes, ligand, pH, V, metal):
    """pM après avoir versé V : dilution à la main, bilans du mélange, dichotomie sur log l."""
    V, V0 = mpf(V), mpf(V0)
    a = _alpha(ligand, pH) if ligand.pKa else mpf(1)
    betas = {
        c.metal: [mp.power(10, mpf(b)) * a ** (j + 1) for j, b in enumerate(c.log_beta)] for c in complexes
    }
    T = {m: (V0 * becher.concentration(m) + V * burette.concentration(m)) / (V0 + V) for m in betas}
    TL = (V0 * mpf(becher.ligand) + V * mpf(burette.ligand)) / (V0 + V)

    def lie(m, lib):
        poids = [b * lib ** (j + 1) for j, b in enumerate(betas[m])]
        return sum((j + 1) * p for j, p in enumerate(poids)) / (1 + sum(poids))

    def bilan(log_l):
        lib = mp.power(10, log_l)
        return lib + sum(T[m] * lie(m, lib) for m in betas) - TL

    bas, haut = mpf(-80), mp.log10(TL + max(T.values()) + 1)
    for _ in range(400):
        milieu = (bas + haut) / 2
        if bilan(milieu) > 0:
            haut = milieu
        else:
            bas = milieu
    lib = mp.power(10, (bas + haut) / 2)
    m = T[metal] / (1 + sum(b * lib ** (j + 1) for j, b in enumerate(betas[metal])))
    return float(-mp.log10(m))


CAS = {
    "Ca-EDTA, pH 10": (
        SolutionComplexation({"Ca2+": 0.01}),
        SolutionComplexation(ligand=0.01),
        20,
        [CA],
        EDTA,
        10.0,
        "Ca2+",
    ),
    "Ca-EDTA, pH 6": (
        SolutionComplexation({"Ca2+": 0.01}),
        SolutionComplexation(ligand=0.01),
        20,
        [CA],
        EDTA,
        6.0,
        "Ca2+",
    ),
    "Ca et Zn en compétition, pZn": (
        SolutionComplexation({"Ca2+": 0.01, "Zn2+": 0.005}),
        SolutionComplexation(ligand=0.01),
        20,
        [CA, ZN],
        EDTA,
        10.0,
        "Zn2+",
    ),
    "Ca et Zn en compétition, pCa": (
        SolutionComplexation({"Ca2+": 0.01, "Zn2+": 0.005}),
        SolutionComplexation(ligand=0.01),
        20,
        [CA, ZN],
        EDTA,
        10.0,
        "Ca2+",
    ),
    "complexes successifs Cu-NH3, ligand titrant": (
        SolutionComplexation({"Cu2+": 0.01}),
        SolutionComplexation(ligand=0.1),
        20,
        [CU],
        NH3,
        11.0,
        "Cu2+",
    ),
    "complexes successifs Cu-NH3, métal titrant": (
        SolutionComplexation(ligand=0.1),
        SolutionComplexation({"Cu2+": 0.05}),
        20,
        [CU],
        NH3,
        11.0,
        "Cu2+",
    ),
    "ligand sans acidité": (
        SolutionComplexation({"Ca2+": 0.01}),
        SolutionComplexation(ligand=0.02),
        10,
        [Complexe("Ca2+", (8.0,))],
        Ligand(),
        None,
        "Ca2+",
    ),
}


def _titrage(cas):
    becher, burette, V0, complexes, ligand, pH, metal = cas
    return TitrageComplexometrique(becher, burette, V0, complexes, ligand, pH, metal)


class TestTitrageComplexometrique:
    @pytest.mark.parametrize("nom", CAS)
    @pytest.mark.parametrize("V", [0.5, 5.0, 9.99, 10.0, 15.0, 19.9, 20.0, 20.1, 25.0, 39.0])
    def test_pM_identique_a_la_reference(self, nom, V):
        cas = CAS[nom]
        becher, burette, V0, complexes, ligand, pH, metal = cas
        if burette.metaux and V == 0.5:
            pytest.skip("le métal titrant est à l'état de traces : pM y est grand mais exact")
        attendu = reference_pM(
            becher, burette, V0, complexes, ligand, pH if pH is not None else 7.0, V, metal
        )
        assert _titrage(cas).pM(V) == pytest.approx(attendu, abs=1e-9)

    @pytest.mark.parametrize("nom", CAS)
    def test_volume_et_pL_reciproques(self, nom):
        t = _titrage(CAS[nom])
        V = np.array([1.0, 7.0, 13.0, 22.0])
        assert t.volume(t.pL(V)) == pytest.approx(V, rel=1e-10)

    @pytest.mark.parametrize("nom", CAS)
    def test_pente_egale_aux_differences_finies(self, nom):
        t = _titrage(CAS[nom])
        for V in (3.0, 12.0, 26.0):
            h = 1e-5
            attendu = (t.pM(V + h) - t.pM(V - h)) / (2 * h)
            assert t.pente(t.pL(V)) == pytest.approx(attendu, rel=1e-5)

    def test_equivalence_monometallique_exacte(self):
        # À l'équivalence d'un complexe 1:1, [M] = [L'] : β' m² + m - T = 0.
        t = _titrage(CAS["Ca-EDTA, pH 10"])
        b = 10 ** (10.7 + math.log10(float(_alpha(EDTA, 10.0))))
        T = 0.01 * 20 / 40
        m = (-1 + math.sqrt(1 + 4 * b * T)) / (2 * b)
        assert t.pM(20) == pytest.approx(-math.log10(m), abs=1e-10)

    def test_saut_a_l_equivalence_pour_un_complexe_stable(self):
        t = _titrage(CAS["Ca-EDTA, pH 10"])
        (saut,) = t.sauts(40)
        assert saut.volume == pytest.approx(20, abs=1e-5)

    def test_saut_est_le_maximum_de_la_pente(self):
        t = _titrage(CAS["Ca-EDTA, pH 6"])
        (saut,) = t.sauts(40)
        pentes = [abs(t.pente(t.pL(saut.volume + d))) for d in (-1e-3, 1e-3)]
        assert abs(saut.pente) > max(pentes)

    def test_milieu_trop_acide_sans_saut(self):
        cas = list(CAS["Ca-EDTA, pH 10"])
        cas[5] = 4.0
        assert _titrage(cas).sauts(40) == ()

    def test_competition_deux_sauts_sur_pZn(self):
        t = _titrage(CAS["Ca et Zn en compétition, pZn"])
        assert [round(s.volume, 3) for s in t.sauts(40)] == [10.0, 30.0]

    def test_courbe_part_de_pM_initial(self):
        c = _titrage(CAS["Ca-EDTA, pH 10"]).courbe(40, 50)
        assert (c.volume[0], c.pM[0]) == pytest.approx((0.0, 2.0), abs=1e-6)

    def test_courbe_finit_a_volume_max(self):
        assert _titrage(CAS["Ca-EDTA, pH 10"]).courbe(40, 50).volume[-1] == pytest.approx(40)


class TestEquivalences:
    def test_ligand_titrant_range_par_stabilite(self):
        t = _titrage(CAS["Ca et Zn en compétition, pZn"])
        assert t.volumes_equivalents() == pytest.approx((10.0, 30.0))

    def test_complexes_successifs_comptent_n_ligands(self):
        t = _titrage(CAS["complexes successifs Cu-NH3, ligand titrant"])
        assert t.volumes_equivalents() == pytest.approx((8.0,))

    def test_metal_titrant(self):
        t = _titrage(CAS["complexes successifs Cu-NH3, métal titrant"])
        assert t.volumes_equivalents() == pytest.approx((10.0,))


class TestErreurs:
    def test_ligand_acide_sans_pH(self):
        with pytest.raises(ValueError, match="pH"):
            TitrageComplexometrique(
                SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01), 20, [CA], EDTA
            )

    def test_metal_sans_complexe(self):
        with pytest.raises(ValueError, match="Mg2"):
            TitrageComplexometrique(
                SolutionComplexation({"Mg2+": 0.01}), SolutionComplexation(ligand=0.01), 20, [CA]
            )

    def test_complexe_sans_constante(self):
        with pytest.raises(ValueError):
            Complexe("Ca2+", ())


def test_fractions_complexes_somment_a_un():
    x = fractions_complexes(CU.log_beta, np.linspace(0, 8, 9))
    assert x.sum(axis=0) == pytest.approx(np.ones(9))
