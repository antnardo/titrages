"""La précipitation à pH calculé, hydroxydes compris, confrontée à un calcul indépendant.

Comme pour la complexométrie, la référence écrit l'**électroneutralité** — contre-ions
écrits à la main, charge de chaque forme protonée des partenaires — là où le paquet écrit
le bilan des protons. Elle trouve un point de départ et l'ensemble des solides présents
par deux dichotomies en double précision, puis résout le système de cet ensemble par
Newton en multiprécision (``mpmath.findroot``, 40 chiffres).
"""

import warnings

import numpy as np
import pytest
from mpmath import mp, mpf

from titrages import (
    Couple,
    Precipite,
    Solution,
    SolutionPrecipitation,
    TitrageParPrecipitation,
    precipite,
)

mp.dps = 40

CL = Precipite("Cl-", 9.74)
CRO4 = Precipite("CrO4^2-", 11.92, nu=2, pKa=(-0.2, 6.51))
SCN = Precipite("SCN-", 11.97, pKa=(0.9,))
MGOH2 = Precipite("OH-", 11.15, nu=1, mu=2, commun="Mg2+")
AG2O = Precipite("OH-", 15.42, nu=2, mu=2, commun="Ag+")

# --- la référence ---------------------------------------------------------------------
#
# Un cas : V0, la charge de l'ion commun, les précipités (partenaire, pKs, nu, mu, pKa,
# charge de la forme la plus basique ; partenaire "OH-" pour un hydroxyde), et pour
# chaque solution : ion commun, partenaires, charge totale des contre-ions.


def _reference(cas: dict, V: float) -> tuple[float, float]:
    V0, V = mpf(cas["V0"]), mpf(V)
    b, u = cas["becher"], cas["burette"]

    def melange(x, y):
        return (V0 * mpf(x) + V * mpf(y)) / (V0 + V)

    TP = melange(b["commun"], u["commun"])
    Q = melange(b["contre_ions"], u["contre_ions"])
    zP = cas["charge_commun"]
    Kw = mp.power(10, -14)
    solides = [p for p in cas["precipites"] if p[0] != "OH-"]
    hydroxydes = [p for p in cas["precipites"] if p[0] == "OH-"]
    TX = {p[0]: melange(b["partenaires"].get(p[0], 0), u["partenaires"].get(p[0], 0)) for p in solides}

    def P(pKa, h):  # somme des K_k h^k, et les fractions de A, HA, H2A…
        d = sorted(pKa, reverse=True)
        termes = [mp.power(10, mpf(sum(d[:k]))) * h**k for k in range(len(d) + 1)]
        return sum(termes), [t / sum(termes) for t in termes]

    def x_sat(p, a, h):
        _, pKs, nu, mu, pKa, _ = p
        return mp.power(mp.power(10, -mpf(pKs)) / a**nu, mpf(1) / mu) * P(pKa, h)[0]

    def a_sat(p, h):
        _, pKs, nu, mu, _, _ = p
        return mp.power(mp.power(10, -mpf(pKs)) / (Kw / h) ** mu, mpf(1) / nu)

    def en_solution(p, a, h, present):
        return x_sat(p, a, h) if present else TX[p[0]]

    def charge(h, a, S):
        z = h - Kw / h + zP * a + Q
        for p in solides:
            x = en_solution(p, a, h, S[p[0]])
            _, f = P(p[4], h)
            z += sum((p[5] + k) * x * fk for k, fk in enumerate(f))
        return z

    def bilan_commun(h, a, S):
        g = a - TP
        for p in solides:
            if S[p[0]]:
                g += mpf(p[2]) / p[3] * (TX[p[0]] - x_sat(p, a, h))
        return g

    def interieur(x):  # log a, et les solides, au pH donné
        h = mp.power(10, mpf(x))

        def g(y):
            a = mp.power(10, mpf(y))
            return float(a - TP + sum(mpf(p[2]) / p[3] * max(0, TX[p[0]] - x_sat(p, a, h)) for p in solides))

        lo, hi = -60.0, 2.0
        for _ in range(70):
            m = (lo + hi) / 2
            lo, hi = (lo, m) if g(m) > 0 else (m, hi)
        y = (lo + hi) / 2
        hyd = None
        for p in hydroxydes:
            borne = float(mp.log10(a_sat(p, h)))
            if y > borne:
                y, hyd = borne, p
        a = mp.power(10, mpf(y))
        S = {p[0]: TX[p[0]] > x_sat(p, a, h) for p in solides}
        return y, S, hyd

    def exterieur(x):
        y, S, _ = interieur(x)
        return float(charge(mp.power(10, mpf(x)), mp.power(10, mpf(y)), S))

    lo, hi = -20.0, 2.0
    signe = exterieur(lo) > 0
    for _ in range(60):
        m = (lo + hi) / 2
        if (exterieur(m) > 0) == signe:
            lo = m
        else:
            hi = m
    x0 = (lo + hi) / 2
    y0, S, hyd = interieur(x0)
    if hyd is not None:
        x = mp.findroot(lambda x: charge(mp.power(10, x), a_sat(hyd, mp.power(10, x)), S), mpf(x0))
        a = a_sat(hyd, mp.power(10, x))
    else:
        x, y = mp.findroot(
            lambda x, y: (
                charge(mp.power(10, x), mp.power(10, y), S),
                bilan_commun(mp.power(10, x), mp.power(10, y), S),
            ),
            (mpf(x0), mpf(y0)),
        )
        a = mp.power(10, y)
    return float(-x), float(-mp.log10(a))


def _p(precipite: Precipite, charge: int) -> tuple:
    return (precipite.partenaire, precipite.pKs, precipite.nu, precipite.mu, precipite.pKa, charge)


MOHR_VOLUMES = [0.5, 5.0, 10.0, 15.0, 19.9, 20.0, 20.1, 25.0]
MG_VOLUMES = [0.1, 1.0, 2.0, 3.9, 4.0, 4.1, 6.0]

# nom : (cas de la référence, titrage du paquet, volumes)
CAS = {
    "Mohr sans tampon": (
        {
            "V0": 20,
            "charge_commun": 1,
            "precipites": [_p(CL, -1), _p(CRO4, -2)],
            "becher": {"commun": 0, "partenaires": {"Cl-": 0.05, "CrO4^2-": 0.005}, "contre_ions": 0.06},
            "burette": {"commun": 0.05, "partenaires": {}, "contre_ions": -0.05},
        },
        lambda: TitrageParPrecipitation(
            SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
            SolutionPrecipitation(commun=0.05),
            20,
            [CL, CRO4],
        ),
        MOHR_VOLUMES,
    ),
    "Mohr en milieu acide nitrique": (
        {
            "V0": 20,
            "charge_commun": 1,
            "precipites": [_p(CL, -1), _p(CRO4, -2)],
            "becher": {
                "commun": 0,
                "partenaires": {"Cl-": 0.05, "CrO4^2-": 0.005},
                "contre_ions": 0.06 - 0.001,
            },
            "burette": {"commun": 0.05, "partenaires": {}, "contre_ions": -0.05},
        },
        lambda: TitrageParPrecipitation(
            SolutionPrecipitation(
                partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}, acidobasique=Solution(acide_fort=0.001)
            ),
            SolutionPrecipitation(commun=0.05),
            20,
            [CL, CRO4],
        ),
        MOHR_VOLUMES,
    ),
    "Charpentier-Volhard en milieu nitrique": (
        {
            "V0": 20,
            "charge_commun": 1,
            "precipites": [_p(CL, -1), _p(SCN, -1)],
            "becher": {"commun": 0.1, "partenaires": {"Cl-": 0.05}, "contre_ions": -0.1 + 0.05 - 0.1},
            "burette": {"commun": 0, "partenaires": {"SCN-": 0.1}, "contre_ions": 0.1},
        },
        lambda: TitrageParPrecipitation(
            SolutionPrecipitation(
                commun=0.1, partenaires={"Cl-": 0.05}, acidobasique=Solution(acide_fort=0.1)
            ),
            SolutionPrecipitation(partenaires={"SCN-": 0.1}),
            20,
            [CL, SCN],
        ),
        [0.5, 5.0, 9.9, 10.0, 10.1, 15.0, 25.0],
    ),
    "Mg(OH)2 par la soude": (
        {
            "V0": 20,
            "charge_commun": 2,
            "precipites": [_p(MGOH2, 0)],
            "becher": {"commun": 0.01, "partenaires": {}, "contre_ions": -0.02},
            "burette": {"commun": 0, "partenaires": {}, "contre_ions": 0.1},
        },
        lambda: TitrageParPrecipitation(
            SolutionPrecipitation(commun=0.01),
            SolutionPrecipitation(acidobasique=Solution(base_forte=0.1)),
            20,
            [MGOH2],
        ),
        MG_VOLUMES,
    ),
    "Mg(OH)2 après un acide fort": (
        {
            "V0": 20,
            "charge_commun": 2,
            "precipites": [_p(MGOH2, 0)],
            "becher": {"commun": 0.01, "partenaires": {}, "contre_ions": -0.02 - 0.002},
            "burette": {"commun": 0, "partenaires": {}, "contre_ions": 0.1},
        },
        lambda: TitrageParPrecipitation(
            SolutionPrecipitation(commun=0.01, acidobasique=Solution(acide_fort=0.002)),
            SolutionPrecipitation(acidobasique=Solution(base_forte=0.1)),
            20,
            [MGOH2],
        ),
        [0.1, 0.39, 0.4, 0.41, 2.0, 4.39, 4.4, 4.41, 6.0],
    ),
    "Ag2O par la soude": (
        {
            "V0": 20,
            "charge_commun": 1,
            "precipites": [_p(AG2O, 0)],
            "becher": {"commun": 0.01, "partenaires": {}, "contre_ions": -0.01},
            "burette": {"commun": 0, "partenaires": {}, "contre_ions": 0.1},
        },
        lambda: TitrageParPrecipitation(
            SolutionPrecipitation(commun=0.01),
            SolutionPrecipitation(acidobasique=Solution(base_forte=0.1)),
            20,
            [AG2O],
        ),
        [0.1, 1.0, 1.9, 2.0, 2.1, 3.0],
    ),
}


def _titrage(nom: str) -> TitrageParPrecipitation:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return CAS[nom][1]()


POINTS = [(nom, V) for nom, (_, _, volumes) in CAS.items() for V in volumes]


class TestReference:
    @pytest.mark.parametrize(("nom", "V"), POINTS)
    def test_pH_et_p_identiques_a_la_reference(self, nom, V):
        t = _titrage(nom)
        pH, p = _reference(CAS[nom][0], V)
        assert (t.pH_du_melange(V), t.p(V)) == pytest.approx((pH, p), abs=1e-9)

    @pytest.mark.parametrize("nom", CAS)
    def test_pente_egale_aux_differences_finies(self, nom):
        t = _titrage(nom)
        Vmax = 30 if "Mohr" in nom or "Volhard" in nom else 8
        c = t.courbe(Vmax, 60)
        anguleux = [e.volume for e in t.evenements(Vmax)]
        for V, pente in zip(c.volume[5::10], c.pente[5::10], strict=True):
            h = 1e-8 * max(V, 1.0)
            if any(abs(V - a) < 1e3 * h for a in anguleux):
                continue
            attendu = (t.p(V + h) - t.p(V - h)) / (2 * h)
            assert pente == pytest.approx(attendu, rel=1e-5, abs=1e-9)


class TestEvenements:
    def test_mohr_les_deux_solides(self):
        t = _titrage("Mohr sans tampon")
        noms = [(e.partenaire, e.apparition) for e in t.evenements(30)]
        assert noms == [("Cl-", True), ("CrO4^2-", True)]

    def test_le_milieu_acide_retarde_le_chromate(self):
        sans = [e.volume for e in _titrage("Mohr sans tampon").evenements(30)][1]
        acide = [e.volume for e in _titrage("Mohr en milieu acide nitrique").evenements(30)][1]
        assert 20 < sans < acide < 20.05

    def test_evenement_sur_la_courbe(self):
        t = _titrage("Mohr en milieu acide nitrique")
        for e in t.evenements(30):
            assert t.p(e.volume) == pytest.approx(e.p, abs=1e-8)

    def test_l_hydroxyde_apparait_a_la_premiere_goutte(self):
        (e,) = _titrage("Mg(OH)2 par la soude").evenements(8)
        assert (e.partenaire, e.apparition) == ("OH-", True) and 0 < e.volume < 0.01

    def test_l_acide_retarde_l_hydroxyde(self):
        (e,) = _titrage("Mg(OH)2 après un acide fort").evenements(8)
        assert 0.4 < e.volume < 0.41


class TestEquivalences:
    def test_hydroxyde(self):
        assert _titrage("Mg(OH)2 par la soude").volumes_equivalents() == pytest.approx((4.0,))

    def test_acide_puis_hydroxyde(self):
        assert _titrage("Mg(OH)2 après un acide fort").volumes_equivalents() == pytest.approx((0.4, 4.4))

    def test_oxyde_d_argent(self):
        assert _titrage("Ag2O par la soude").volumes_equivalents() == pytest.approx((2.0,))

    def test_mohr_inchangees(self):
        assert _titrage("Mohr sans tampon").volumes_equivalents() == pytest.approx((20.0, 24.0))


class TestLimites:
    def test_partenaire_sans_acidite_meme_p_qu_au_niveau_1(self):
        """Le pH se calcule (l'acide du bécher), mais il ne compte pas : pP est celui du niveau 1."""
        sans = TitrageParPrecipitation(
            SolutionPrecipitation(partenaires={"Cl-": 0.1}), SolutionPrecipitation(commun=0.1), 20, [CL]
        )
        with pytest.warns(UserWarning):
            avec = TitrageParPrecipitation(
                SolutionPrecipitation(partenaires={"Cl-": 0.1}, acidobasique=Solution(acide_fort=0.01)),
                SolutionPrecipitation(commun=0.1),
                20,
                [CL],
            )
        V = [1.0, 10.0, 19.9, 20.0, 20.1, 30.0]
        assert avec.p(V) == pytest.approx(sans.p(V), abs=1e-12)

    def test_tampon_concentre_tend_vers_le_niveau_1(self):
        """Mohr dans un tampon de plus en plus concentré : pAg tend vers celui du pH imposé."""
        tampon_couple = Couple("tampon", (7.5,))
        V = np.array([10.0, 19.9, 20.0, 20.1, 25.0])
        ecarts = []
        for C in (0.01, 0.1, 1.0):
            tampon = Solution().avec(tampon_couple.forme(0), C).avec(tampon_couple.forme(1), C)
            impose = TitrageParPrecipitation(
                SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}),
                SolutionPrecipitation(commun=0.05),
                20,
                [CL, CRO4],
                pH=tampon.pH(),
            )
            with pytest.warns(UserWarning):
                calcule = TitrageParPrecipitation(
                    SolutionPrecipitation(partenaires={"Cl-": 0.05, "CrO4^2-": 0.005}, acidobasique=tampon),
                    SolutionPrecipitation(commun=0.05, acidobasique=tampon),
                    20,
                    [CL, CRO4],
                )
            ecarts.append(float(np.max(np.abs(calcule.p(V) - impose.p(V)))))
        assert ecarts[1] < ecarts[0] / 5 and ecarts[2] < ecarts[1] / 5 and ecarts[2] < 1e-4


class TestCourbe:
    def test_courbe_porte_le_pH_et_les_solides(self):
        c = _titrage("Mg(OH)2 par la soude").courbe(8, 50)
        assert c.volume[0] == 0 and c.volume[-1] == pytest.approx(8)
        assert c.solides[0] == () and c.solides[-1] == ("OH-",)
        assert c.pH[0] == pytest.approx(7, abs=1e-6)

    def test_sauts_dans_les_morceaux_lisses(self):
        t = _titrage("Mohr sans tampon")
        assert round(t.sauts(30)[0].volume, 3) == 20.0

    def test_base_hydroxyde(self):
        p = precipite("hydroxyde de magnésium")
        assert (p.commun, p.partenaire, p.nu, p.mu, p.pKs) == ("Mg2+", "OH-", 1, 2, 11.15)


class TestErreurs:
    def test_hydroxyde_a_pH_impose(self):
        with pytest.raises(ValueError, match="hydroxyde"):
            TitrageParPrecipitation(
                SolutionPrecipitation(commun=0.01),
                SolutionPrecipitation(acidobasique=Solution(base_forte=0.1)),
                20,
                [MGOH2],
                pH=10,
            )

    def test_hydroxyde_dans_les_partenaires(self):
        with pytest.raises(ValueError, match="hydroxyde"):
            SolutionPrecipitation(partenaires={"OH-": 0.1})

    def test_parametre_pP_sans_objet(self):
        t = _titrage("Mohr sans tampon")
        for methode in (t.volume, t.derivee_volume, t.pente, t.solides):
            with pytest.raises(ValueError, match="volume"):
                methode(5.0)

    def test_avertissement(self):
        with pytest.warns(UserWarning, match="pH non précisé"):
            CAS["Mg(OH)2 par la soude"][1]()
