"""La complexométrie à pH calculé, confrontée à un calcul indépendant.

La référence n'emprunte rien au paquet, ni sa formulation : là où le paquet écrit le
bilan des **protons**, elle écrit l'**électroneutralité**, avec les contre-ions écrits à
la main (Cl⁻ de CaCl₂, Na⁺ de Na₂H₂Y, K⁺ de KCN…), les charges de chaque forme du
ligand, de chaque complexe et de chaque couple. Elle dilue les concentrations au volume
versé, trouve un point de départ par deux dichotomies en double précision, puis résout
le système complet par Newton en multiprécision (``mpmath.findroot``, 40 chiffres).
"""

import math
import warnings

import numpy as np
import pytest
from mpmath import mp, mpf

from titrages import (
    Complexe,
    Ligand,
    Solution,
    SolutionComplexation,
    TitrageComplexometrique,
    complexe,
    couple,
    ligand,
)

mp.dps = 40

EDTA = (0.0, 1.5, 2.0, 2.66, 6.16, 10.24)
NH4 = couple("NH3")


def _tampon(C: float) -> Solution:
    return Solution().avec(NH4.forme(0), C).avec(NH4.forme(1), C)


# --- la référence ---------------------------------------------------------------------
#
# Un cas : V0, le ligand (pKa, charge de sa forme la plus basique), les complexes
# (log beta globaux, -inf permis), et pour chaque solution : métaux {nom: (C, charge)},
# ligand C, charge totale des contre-ions, couples [(pKa, C total, charge de la forme
# la plus acide)].


def _reference(cas: dict, V: float, metal: str) -> tuple[float, float]:
    V0, V = mpf(cas["V0"]), mpf(V)

    def melange(b, u):
        return (V0 * mpf(b) + V * mpf(u)) / (V0 + V)

    b, u = cas["becher"], cas["burette"]
    metaux = set(b["metaux"]) | set(u["metaux"])
    T = {m: melange(b["metaux"].get(m, (0, 0))[0], u["metaux"].get(m, (0, 0))[0]) for m in metaux}
    charge = {m: (b["metaux"].get(m) or u["metaux"].get(m))[1] for m in metaux}
    TL = melange(b["ligand"], u["ligand"])
    Q = melange(b["contre_ions"], u["contre_ions"])
    couples = [
        (pKa, melange(Cb, Cu), z) for (pKa, Cb, z), (_, Cu, _) in zip(b["couples"], u["couples"], strict=True)
    ]
    pKa_L, zL = cas["ligand"]
    Kw = mp.power(10, -14)

    def log_K(pKa):  # [H_k A] = K_k h^k [A]
        d = sorted(pKa, reverse=True)
        return [sum(d[:k]) for k in range(len(d) + 1)]

    def formes(pKa, h):  # fractions de A, HA, H2A… (de la plus basique à la plus acide)
        termes = [mp.power(10, mpf(lk)) * h**k for k, lk in enumerate(log_K(pKa))]
        s = sum(termes)
        return [t / s for t in termes]

    betas = {
        m: [mpf(1)] + [mpf(0) if lb == -math.inf else mp.power(10, mpf(lb)) for lb in cas["complexes"][m]]
        for m in metaux
    }

    def fractions_metal(m, lib):
        poids = [bj * lib**j for j, bj in enumerate(betas[m])]
        s = sum(poids)
        return [p / s for p in poids]

    def bilans(x, y):  # x = log h, y = log l
        h, lib = mp.power(10, x), mp.power(10, y)
        lk = log_K(pKa_L)
        libres = [mp.power(10, mpf(c)) * h**k * lib for k, c in enumerate(lk)]
        lie = sum(T[m] * sum(j * f for j, f in enumerate(fractions_metal(m, lib))) for m in metaux)
        g_L = sum(libres) + lie - TL
        z = h - Kw / h + Q
        z += sum((zL + k) * c for k, c in enumerate(libres))
        for m in metaux:
            z += T[m] * sum((charge[m] + j * zL) * f for j, f in enumerate(fractions_metal(m, lib)))
        for pKa, C, z_acide in couples:  # z de la forme H_N A, la plus acide
            N = len(pKa)
            z += C * sum((z_acide - N + k) * f for k, f in enumerate(formes(pKa, h)))
        return z, g_L

    def l_de(x):
        a, c = -120.0, 2.0
        for _ in range(70):
            m = (a + c) / 2
            a, c = (a, m) if float(bilans(mpf(x), mpf(m))[1]) > 0 else (m, c)
        return (a + c) / 2

    a, c = -20.0, 3.0
    signe = float(bilans(mpf(a), mpf(l_de(a)))[0]) > 0
    for _ in range(60):
        m = (a + c) / 2
        if (float(bilans(mpf(m), mpf(l_de(m)))[0]) > 0) == signe:
            a = m
        else:
            c = m
    x0 = (a + c) / 2
    x, y = mp.findroot(lambda x, y: bilans(x, y), (mpf(x0), mpf(l_de(x0))))
    lib = mp.power(10, y)
    libre = T[metal] * fractions_metal(metal, lib)[0]
    return float(-x), float(-mp.log10(libre))


def _solution(d: dict, acidobasique: Solution, protons: int) -> SolutionComplexation:
    return SolutionComplexation(
        {m: c for m, (c, _) in d["metaux"].items()}, d["ligand"], protons, acidobasique
    )


CA = Complexe("Ca2+", (10.7,))
ZN = Complexe("Zn2+", (16.5,))

# nom : (cas de la référence, titrage du paquet, métal suivi)
CAS = {
    "Ca par Na2H2Y, sans tampon": (
        {
            "V0": 20,
            "ligand": (EDTA, -4),
            "complexes": {"Ca2+": (10.7,)},
            "becher": {"metaux": {"Ca2+": (0.01, 2)}, "ligand": 0, "contre_ions": -0.02, "couples": []},
            "burette": {"metaux": {}, "ligand": 0.01, "contre_ions": 0.02, "couples": []},
        },
        lambda: TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01}),
            SolutionComplexation(ligand=0.01, protons_ligand=2),
            20,
            [CA],
            Ligand("Y4-", EDTA),
        ),
        "Ca2+",
    ),
    "Ca par Na2H2Y, tampon ammoniacal en solution": (
        {
            "V0": 20,
            "ligand": (EDTA, -4),
            "complexes": {"Ca2+": (10.7,)},
            "becher": {
                "metaux": {"Ca2+": (0.01, 2)},
                "ligand": 0,
                "contre_ions": -0.12,
                "couples": [((9.25,), 0.2, 1)],
            },
            "burette": {
                "metaux": {},
                "ligand": 0.01,
                "contre_ions": 0.02 - 0.1,
                "couples": [((9.25,), 0.2, 1)],
            },
        },
        lambda: TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01}, acidobasique=_tampon(0.1)),
            SolutionComplexation(ligand=0.01, protons_ligand=2, acidobasique=_tampon(0.1)),
            20,
            [CA],
            Ligand("Y4-", EDTA),
        ),
        "Ca2+",
    ),
    "Zn par Na2H2Y, avec de l'acide fort": (
        {
            "V0": 20,
            "ligand": (EDTA, -4),
            "complexes": {"Zn2+": (16.5,)},
            "becher": {"metaux": {"Zn2+": (0.005, 2)}, "ligand": 0, "contre_ions": -0.012, "couples": []},
            "burette": {"metaux": {}, "ligand": 0.01, "contre_ions": 0.02, "couples": []},
        },
        lambda: TitrageComplexometrique(
            SolutionComplexation({"Zn2+": 0.005}, acidobasique=Solution(acide_fort=0.002)),
            SolutionComplexation(ligand=0.01, protons_ligand=2),
            20,
            [ZN],
            Ligand("Y4-", EDTA),
        ),
        "Zn2+",
    ),
    "Cu par l'ammoniac, sans tampon": (
        {
            "V0": 20,
            "ligand": ((9.244,), 0),
            "complexes": {"Cu2+": (4.04, 7.47, 10.27, 11.75)},
            "becher": {"metaux": {"Cu2+": (0.01, 2)}, "ligand": 0, "contre_ions": -0.02, "couples": []},
            "burette": {"metaux": {}, "ligand": 0.1, "contre_ions": 0, "couples": []},
        },
        lambda: TitrageComplexometrique(
            SolutionComplexation({"Cu2+": 0.01}),
            SolutionComplexation(ligand=0.1),
            20,
            [complexe("Cu2+", "NH3")],
            ligand("NH3"),
        ),
        "Cu2+",
    ),
    "Liebig, avec la protonation de HCN": (
        {
            "V0": 20,
            "ligand": ((9.21,), -1),
            "complexes": {"Ag+": (-math.inf, 20.48, 21.40)},
            "becher": {"metaux": {}, "ligand": 0.1, "contre_ions": 0.1, "couples": []},
            "burette": {"metaux": {"Ag+": (0.05, 1)}, "ligand": 0, "contre_ions": -0.05, "couples": []},
        },
        lambda: TitrageComplexometrique(
            SolutionComplexation(ligand=0.1),
            SolutionComplexation({"Ag+": 0.05}),
            20,
            [complexe("Ag+", "cyanure")],
            ligand("cyanure"),
        ),
        "Ag+",
    ),
}
VOLUMES = [0.5, 5.0, 10.0, 15.0, 19.9, 20.0, 20.1, 25.0, 39.0]


def _titrage(nom: str) -> TitrageComplexometrique:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return CAS[nom][1]()


class TestReference:
    @pytest.mark.parametrize("nom", CAS)
    @pytest.mark.parametrize("V", VOLUMES)
    def test_pH_et_pM_identiques_a_la_reference(self, nom, V):
        cas, _, metal = CAS[nom]
        t = _titrage(nom)
        pH, pM = _reference(cas, V, metal)
        assert (t.pH_du_melange(V), t.pM(V)) == pytest.approx((pH, pM), abs=1e-9)

    @pytest.mark.parametrize("nom", CAS)
    def test_pente_egale_aux_differences_finies(self, nom):
        t = _titrage(nom)
        c = t.courbe(40, 60)
        for V, pente in zip(c.volume[5::10], c.pente[5::10], strict=True):
            h = 1e-8 * max(V, 1.0)  # petit : dans un saut, la courbure fausse un pas plus grand
            attendu = (t.pM(V + h) - t.pM(V - h)) / (2 * h)
            assert pente == pytest.approx(attendu, rel=1e-5, abs=1e-9)


class TestAvertissement:
    def test_sans_pH_le_pH_se_calcule_et_le_dit(self):
        with pytest.warns(UserWarning, match="pH non précisé"):
            CAS["Ca par Na2H2Y, sans tampon"][1]()

    def test_avec_pH_rien(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            TitrageComplexometrique(
                SolutionComplexation({"Ca2+": 0.01}),
                SolutionComplexation(ligand=0.01, protons_ligand=2),
                20,
                [CA],
                Ligand("Y4-", EDTA),
                pH=10,
            )

    def test_sans_acidite_rien(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            TitrageComplexometrique(
                SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01), 20, [CA]
            )


class TestLimites:
    def test_ligand_sans_acidite_meme_pM_qu_au_niveau_1(self):
        """Le pH se calcule (la soude du bécher), mais il ne compte pas : pM est celui du niveau 1."""
        sans = TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01}),
            SolutionComplexation(ligand=0.02),
            10,
            [Complexe("Ca2+", (8.0,))],
        )
        with pytest.warns(UserWarning):
            avec = TitrageComplexometrique(
                SolutionComplexation({"Ca2+": 0.01}, acidobasique=Solution(base_forte=0.001)),
                SolutionComplexation(ligand=0.02),
                10,
                [Complexe("Ca2+", (8.0,))],
            )
        V = [0.5, 4.0, 5.0, 6.0, 9.0]
        assert avec.pM(V) == pytest.approx(sans.pM(V), abs=1e-12)

    def test_tampon_concentre_tend_vers_le_niveau_1(self):
        """Un tampon de plus en plus concentré rapproche pM de celui du pH imposé."""
        V = np.array([10.0, 19.9, 20.0, 20.1, 30.0])
        ecarts = []
        for C in (0.1, 1.0, 10.0):
            tampon = _tampon(C)
            impose = TitrageComplexometrique(
                SolutionComplexation({"Ca2+": 0.01}),
                SolutionComplexation(ligand=0.01),
                20,
                [CA],
                Ligand("Y4-", EDTA),
                pH=tampon.pH(),
            )
            with pytest.warns(UserWarning):
                calcule = TitrageComplexometrique(
                    SolutionComplexation({"Ca2+": 0.01}, acidobasique=tampon),
                    SolutionComplexation(ligand=0.01, protons_ligand=2, acidobasique=tampon),
                    20,
                    [CA],
                    Ligand("Y4-", EDTA),
                )
            ecarts.append(float(np.max(np.abs(calcule.pM(V) - impose.pM(V)))))
        assert ecarts[1] < ecarts[0] / 5 and ecarts[2] < ecarts[1] / 5 and ecarts[2] < 1e-3


class TestCourbe:
    def test_sans_tampon_pas_de_saut(self):
        """Le pH tombe à 3,5 : la constante conditionnelle s'effondre, le saut disparaît."""
        t = _titrage("Ca par Na2H2Y, sans tampon")
        assert t.sauts(40) == ()
        assert float(t.pH_du_melange(20)) == pytest.approx(3.533, abs=1e-3)

    def test_tampon_en_solution_saut_a_l_equivalence(self):
        t = _titrage("Ca par Na2H2Y, tampon ammoniacal en solution")
        (saut,) = t.sauts(40)
        assert saut.volume == pytest.approx(20, abs=1e-4)

    def test_saut_est_le_maximum_de_la_pente(self):
        t = _titrage("Ca par Na2H2Y, tampon ammoniacal en solution")
        (saut,) = t.sauts(40)
        c = t.courbe(40, 400)
        assert abs(saut.pente) >= np.max(np.abs(c.pente)) * (1 - 1e-9)

    def test_liebig_saut_a_Ag_CN_2(self):
        t = _titrage("Liebig, avec la protonation de HCN")
        assert [round(s.volume, 3) for s in t.sauts(40)] == [20.0]

    def test_courbe_part_de_zero_et_porte_le_pH(self):
        c = _titrage("Ca par Na2H2Y, sans tampon").courbe(40, 50)
        assert len(c) == 50 and c.volume[0] == 0 and c.volume[-1] == pytest.approx(40)
        assert np.all(np.isfinite(c.pH)) and c.pH[0] == pytest.approx(7, abs=1e-6)

    def test_metal_titrant_part_d_un_milliardieme(self):
        c = _titrage("Liebig, avec la protonation de HCN").courbe(40, 50)
        assert c.volume[0] == pytest.approx(4e-8)


class TestErreurs:
    def test_parametre_pL_sans_objet(self):
        t = _titrage("Ca par Na2H2Y, sans tampon")
        for methode in (t.volume, t.derivee_volume, t.pente):
            with pytest.raises(ValueError, match="volume"):
                methode(5.0)

    def test_trop_de_protons(self):
        with pytest.raises(ValueError, match="protons"), pytest.warns(UserWarning):
            TitrageComplexometrique(
                SolutionComplexation({"Cu2+": 0.01}),
                SolutionComplexation(ligand=0.1, protons_ligand=2),
                20,
                [complexe("Cu2+", "NH3")],
                ligand("NH3"),
            )

    def test_pH_sans_objet(self):
        t = TitrageComplexometrique(
            SolutionComplexation({"Ca2+": 0.01}), SolutionComplexation(ligand=0.01), 20, [CA]
        )
        with pytest.raises(ValueError, match="objet"):
            t.pH_du_melange(5)
