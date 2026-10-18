"""La précipitation à pH calculé : le niveau 2, sans tampon imposé, hydroxydes compris.

**Les inconnues** sont, à volume versé fixé, le pH et :math:`u = \\log a`, l'ion commun
libre. Les partenaires sont *éliminés* : au contact de son solide, un partenaire a sa
concentration de saturation, toutes formes protonées comprises,
:math:`x_\\text{sat} = (K_s/a^\\nu)^{1/\\mu}\\,P(h)`, avec :math:`P = \\sum_k K_k h^k` ;
sans solide, il est en solution tout entier. Deux bilans :

- celui de l'ion commun, :math:`g_P = a + \\sum_{i\\in S} r_i\\,(T_{X_i} - x_{\\text{sat},i}) - T_P` ;
- celui des protons, :math:`g_H = R_\\text{ab}(\\mathrm{pH}) - \\sum_i x'_i\\,m_i`, où
  :math:`R_\\text{ab}` est le résidu de l'électroneutralité de la partie acido-basique des
  solutions, :math:`x'_i` le partenaire en solution et :math:`m_i` le nombre moyen de
  protons qu'il porte : un partenaire s'introduit sous sa forme la plus basique (NaCl,
  K₂CrO₄), et ce qu'il fixe de protons en solution manque au bilan.

**Un hydroxyde** :math:`\\mathrm{P}_\\nu(\\mathrm{OH})_\\mu` — un précipité dont le partenaire
est ``"OH-"`` — borne l'ion commun : :math:`a \\le a_\\text{sat}(h)`, avec
:math:`K_s = a^\\nu\\,(K_e/h)^\\mu`. Présent, il fixe :math:`u` au pH, et sa quantité
:math:`s = -g_P/\\nu` libère :math:`\\mu s` protons : l'équation du pH devient
:math:`g_H - (\\mu/\\nu)\\,g_P = 0`.

**L'unicité** vient de la convexité, solides compris : la dichotomie sur :math:`u` à pH
fixé (le bilan de l'ion commun, avec :math:`\\max(0, T_X - x_\\text{sat})` pour chaque
partenaire, puis la borne des hydroxydes), et la dichotomie sur le pH autour d'elle,
trouvent l'unique équilibre. Newton, à ensemble de solides fixé, affine autour d'un point
connu (:mod:`titrages.libre`).
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from titrages.equilibre import LN10, MARGE_PH, changement_de_signe, racine_croissante
from titrages.libre import Etat, MoteurLibre
from titrages.modeles import Evenement, Precipite

if TYPE_CHECKING:
    from titrages.precipitation import SolutionPrecipitation

__all__ = ["PrecipitationLibre"]

Tableau = npt.NDArray[np.float64]

P_MAX = 60.0


@dataclass(frozen=True)
class PrecipitationLibre(MoteurLibre):
    """Le moteur du titrage par précipitation dont le pH se calcule."""

    becher: "SolutionPrecipitation"
    burette: "SolutionPrecipitation"
    volume_becher: float
    precipites: Sequence[Precipite]
    _log_K: tuple[Tableau, ...] = field(init=False, repr=False)
    _pKe: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        tables = []
        for p in self.precipites:
            decroissants = sorted(p.pKa, reverse=True)
            tables.append(np.concatenate(([0.0], np.cumsum(decroissants))) if decroissants else np.zeros(1))
        object.__setattr__(self, "_log_K", tuple(np.asarray(t, dtype=np.float64) for t in tables))
        a, b = self.becher.acidobasique.pKe, self.burette.acidobasique.pKe
        if a != b:
            raise ValueError("le bécher et la burette n'ont pas le même pKe")
        object.__setattr__(self, "_pKe", a)

    # --- les briques -------------------------------------------------------------------

    def _hydroxyde(self, i: int) -> bool:
        return self.precipites[i].partenaire == "OH-"

    @property
    def _partenaires(self) -> list[int]:
        return [i for i in range(len(self.precipites)) if not self._hydroxyde(i)]

    @property
    def _hydroxydes(self) -> list[int]:
        return [i for i in range(len(self.precipites)) if self._hydroxyde(i)]

    def _melange(self, V: Tableau, b: float, u: float) -> Tableau:
        V0 = self.volume_becher
        return (V0 * b + V * u) / (V0 + V)

    def _total(self, i: int, V: Tableau) -> Tableau:
        x = self.precipites[i].partenaire
        return self._melange(V, self.becher.concentration(x), self.burette.concentration(x))

    def _protons(self, i: int, pH: Tableau) -> tuple[Tableau, Tableau, Tableau, Tableau]:
        """:math:`\\log P`, et les protons du partenaire : moyenne, moyenne des carrés, variance."""
        lk = self._log_K[i]
        k = np.arange(len(lk), dtype=np.float64).reshape((-1,) + (1,) * pH.ndim)
        t = lk.reshape(k.shape) - k * pH
        haut = t.max(axis=0)
        e = np.power(10.0, t - haut)
        somme = e.sum(axis=0)
        x = e / somme
        moyenne, carre = (k * x).sum(axis=0), (k * k * x).sum(axis=0)
        var = np.zeros(pH.shape)
        for a in range(len(lk)):
            for b in range(a + 1, len(lk)):
                var = var + x[a] * x[b] * (b - a) ** 2
        return haut + np.log10(somme), moyenne, carre, var

    def _x_sat(self, i: int, pH: Tableau, u: Tableau) -> Tableau:
        """Le partenaire au contact de son solide, toutes formes protonées comprises."""
        p = self.precipites[i]
        return np.power(10.0, (-p.pKs - p.nu * u) / p.mu + self._protons(i, pH)[0])

    def _u_sat(self, i: int, pH: Tableau) -> Tableau:
        """:math:`\\log a` au contact de l'hydroxyde ``i`` : :math:`K_s = a^\\nu (K_e/h)^\\mu`."""
        p = self.precipites[i]
        return (-p.pKs - p.mu * (pH - self._pKe)) / p.nu

    def _present(self, S: npt.NDArray[np.int64], i: int) -> npt.NDArray[np.bool_]:
        return ((S >> i) & 1).astype(bool)

    def _hydroxyde_present(self, S: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
        """L'indice de l'hydroxyde présent en chaque point, ou -1."""
        k = np.full(S.shape, -1, dtype=np.int64)
        for i in reversed(self._hydroxydes):
            k = np.where(self._present(S, i), i, k)
        return k

    def _en_solution(self, e: Etat, totaux: dict[int, Tableau]) -> dict[int, Tableau]:
        """:math:`x'_i` : la saturation si le solide est présent, sinon le total donné."""
        return {
            i: np.where(self._present(e.S, i), self._x_sat(i, e.pH, e.u), totaux[i])
            for i in self._partenaires
        }

    def _acidobasique(self, V: Tableau, pH: Tableau) -> tuple[Tableau, Tableau]:
        b, u = self.becher.acidobasique, self.burette.acidobasique
        V0 = self.volume_becher
        R = (V0 * np.asarray(b.residu(pH)) + V * np.asarray(u.residu(pH))) / (V0 + V)
        D = (V0 * np.asarray(b.pouvoir_tampon(pH)) + V * np.asarray(u.pouvoir_tampon(pH))) / (V0 + V)
        return R, D

    def _bilans(self, e: Etat) -> tuple[Tableau, Tableau]:
        """:math:`(g_H, g_P)` du mélange, à ensemble ``S`` fixé, sans les hydroxydes."""
        V = e.volume
        T = {i: self._total(i, V) for i in self._partenaires}
        x = self._en_solution(e, T)
        T_P = self._melange(V, self.becher.commun, self.burette.commun)
        gP = np.power(10.0, e.u) - T_P
        gH, _ = self._acidobasique(V, e.pH)
        for i in self._partenaires:
            gP = gP + self.precipites[i].nu / self.precipites[i].mu * (T[i] - x[i])
            gH = gH - x[i] * self._protons(i, e.pH)[1]
        return gH, gP

    def _jacobien(self, e: Etat) -> tuple[Tableau, Tableau, Tableau]:
        """Les dérivées des bilans : :math:`\\partial g_H/\\partial\\mathrm{pH}`, puis
        :math:`\\partial g_H/\\partial u = \\partial g_P/\\partial\\mathrm{pH}`, puis
        :math:`\\partial g_P/\\partial u`."""
        _, D = self._acidobasique(e.volume, e.pH)
        a11, a12, a22 = D, np.zeros(len(e)), LN10 * np.power(10.0, e.u)
        for i in self._partenaires:
            present = self._present(e.S, i)
            _, m, carre, var = self._protons(i, e.pH)
            x = self._x_sat(i, e.pH, e.u)
            r = self.precipites[i].nu / self.precipites[i].mu
            a11 = a11 + LN10 * np.where(present, x * carre, self._total(i, e.volume) * var)
            a12 = a12 + LN10 * np.where(present, r * x * m, 0.0)
            a22 = a22 + LN10 * np.where(present, r * r * x, 0.0)
        return a11, a12, a22

    # --- la résolution -----------------------------------------------------------------

    def _pas_newton(self, e: Etat) -> tuple[Tableau, Tableau]:
        k = self._hydroxyde_present(e.S)
        # avec un hydroxyde, u suit le pH : u = u_sat(pH), et l'équation est g_H - c g_P
        c = np.ones(len(e))
        u_sat = e.u.copy()
        for i in self._hydroxydes:
            p = self.precipites[i]
            c = np.where(k == i, p.mu / p.nu, c)
            u_sat = np.where(k == i, self._u_sat(i, e.pH), u_sat)
        f = Etat(e.volume, e.pH, u_sat, e.S)
        gH, gP = self._bilans(f)
        a11, a12, a22 = self._jacobien(f)
        det = a11 * a22 - a12 * a12
        dp2, du2 = -(a22 * gH - a12 * gP) / det, -(a11 * gP - a12 * gH) / det
        dp1 = -(gH - c * gP) / (a11 - 2 * c * a12 + c * c * a22)
        avec = k >= 0
        dp = np.where(avec, dp1, dp2)
        du = np.where(avec, u_sat - c * dp - e.u, du2)
        return dp, du

    def _coherent(self, e: Etat) -> npt.NDArray[np.bool_]:
        """Les solides présents ont une quantité positive, les absents ne sont pas sursaturés."""
        ok = np.isfinite(e.pH) & np.isfinite(e.u)
        for i in self._partenaires:
            T, x = self._total(i, e.volume), self._x_sat(i, e.pH, e.u)
            marge = 1e-12 * np.maximum(T, x)
            ok &= np.where(self._present(e.S, i), T - x >= -marge, T - x <= marge)
        k = self._hydroxyde_present(e.S)
        _, gP = self._bilans(e)
        T_P = self._melange(e.volume, self.becher.commun, self.burette.commun)
        for i in self._hydroxydes:
            ecart = e.u - self._u_sat(i, e.pH)
            ok &= np.where(k == i, (np.abs(ecart) < 1e-9) & (-gP >= -1e-12 * T_P), ecart <= 1e-12)
        return ok

    def _u_de(self, V: Tableau, pH: Tableau) -> tuple[Tableau, Tableau, npt.NDArray[np.int64]]:
        """:math:`u` à pH fixé, la quantité d'hydroxyde (en ion commun) et ``S``."""
        T = {i: self._total(i, V) for i in self._partenaires}
        T_P = self._melange(V, self.becher.commun, self.burette.commun)

        def g(u: Tableau) -> Tableau:
            total = np.power(10.0, u) - T_P
            for i in self._partenaires:
                p = self.precipites[i]
                total = total + p.nu / p.mu * np.maximum(0.0, T[i] - self._x_sat(i, pH, u))
            return total

        concentrations = [self.becher.commun, self.burette.commun]
        concentrations += [*self.becher.partenaires.values(), *self.burette.partenaires.values()]
        haut = math.log10(max(concentrations) + 1e-300) + 1
        u = np.asarray(
            racine_croissante(g, np.full(V.shape, -P_MAX), np.full(V.shape, haut)), dtype=np.float64
        )
        S = np.zeros(V.shape, dtype=np.int64)
        borne, laquelle = np.full(V.shape, np.inf), np.full(V.shape, -1)
        for i in self._hydroxydes:
            ui = self._u_sat(i, pH)
            laquelle = np.where(ui < borne, i, laquelle)
            borne = np.minimum(borne, ui)
        bornee = u > borne
        u = np.where(bornee, borne, u)
        exces = np.where(bornee, -g(u), 0.0)  # l'ion commun pris par l'hydroxyde
        for i in self._hydroxydes:
            S = S | np.where(bornee & (laquelle == i), 1 << i, 0)
        for i in self._partenaires:
            S = S | np.where(T[i] > self._x_sat(i, pH, u), 1 << i, 0)
        return u, exces, S

    def _bilan_protons(self, V: Tableau, pH: Tableau) -> Tableau:
        u, exces, S = self._u_de(V, pH)
        gH, _ = self._bilans(Etat(V, pH, u, S))
        k = self._hydroxyde_present(S)
        for i in self._hydroxydes:
            p = self.precipites[i]
            gH = gH + np.where(k == i, p.mu / p.nu * exces, 0.0)
        return gH

    def _dichotomie(self, V: Tableau) -> Etat:
        pH = np.asarray(
            racine_croissante(
                lambda p: self._bilan_protons(V, p),
                np.full(V.shape, -MARGE_PH),
                np.full(V.shape, self._pKe + MARGE_PH),
            ),
            dtype=np.float64,
        )
        u, _, S = self._u_de(V, pH)
        return Etat(V, pH, u, S)

    # --- les grandeurs -----------------------------------------------------------------

    def depart(self, volume_max: float, fraction: float = 1e-6) -> float:
        """0, ou ``fraction`` de ``volume_max`` sans ion commun au bécher (pP y est infini)."""
        return 0.0 if self.becher.commun > 0 else fraction * volume_max

    def p(self, e: Etat) -> Tableau:
        return np.asarray(-e.u, dtype=np.float64)

    def p_partenaire(self, e: Etat, partenaire: str) -> Tableau:
        i = [p.partenaire for p in self.precipites].index(partenaire)
        x = self._en_solution(e, {j: self._total(j, e.volume) for j in self._partenaires})[i]
        with np.errstate(divide="ignore"):
            return np.asarray(-np.log10(x), dtype=np.float64)

    def solides(self, e: Etat) -> tuple[tuple[str, ...], ...]:
        return tuple(
            tuple(p.partenaire for i, p in enumerate(self.precipites) if (int(S) >> i) & 1) for S in e.S
        )

    def derivees(self, e: Etat) -> tuple[Tableau, Tableau]:
        """:math:`\\mathrm{d\\,pH}/\\mathrm{d}V` et :math:`\\mathrm{d}u/\\mathrm{d}V`, à ``S`` fixé."""
        V, s = e.volume, self.burette
        a11, a12, a22 = self._jacobien(e)
        C = {i: np.full(len(e), s.concentration(self.precipites[i].partenaire)) for i in self._partenaires}
        x = self._en_solution(e, C)
        eP = np.power(10.0, e.u) - s.commun
        eH = np.asarray(s.acidobasique.residu(e.pH), dtype=np.float64)
        for i in self._partenaires:
            p = self.precipites[i]
            eP = eP + p.nu / p.mu * (C[i] - x[i])
            eH = eH - x[i] * self._protons(i, e.pH)[1]
        eP, eH = eP / (self.volume_becher + V), eH / (self.volume_becher + V)
        det = a11 * a22 - a12 * a12
        dp, du = -(a22 * eH - a12 * eP) / det, -(a11 * eP - a12 * eH) / det
        k = self._hydroxyde_present(e.S)
        for i in self._hydroxydes:
            c = self.precipites[i].mu / self.precipites[i].nu
            dp1 = -(eH - c * eP) / (a11 - 2 * c * a12 + c * c * a22)
            dp, du = np.where(k == i, dp1, dp), np.where(k == i, -c * dp1, du)
        return dp, du

    def pente(self, e: Etat) -> Tableau:
        """:math:`\\mathrm{d\\,pP}/\\mathrm{d}V`."""
        return np.asarray(-self.derivees(e)[1], dtype=np.float64)

    # --- les événements ----------------------------------------------------------------

    def _sursaturation(self, i: int, e: Etat) -> float:
        """Nulle là où le solide ``i`` apparaît ou disparaît, le long du morceau de ``e``."""
        if not self._hydroxyde(i):
            return float((self._total(i, e.volume) - self._x_sat(i, e.pH, e.u))[0])
        if self._present(e.S, i)[0]:
            return float(-self._bilans(e)[1][0])  # la quantité d'hydroxyde
        return float((e.u - self._u_sat(i, e.pH))[0])

    def evenements(self, volume_max: float, points: int) -> tuple[Evenement, ...]:
        """Les apparitions et disparitions de solides, dans l'ordre, chacune affinée par
        dichotomie en :math:`V` sur le prolongement du morceau de courbe qui la précède."""
        # dès un milliardième de volume_max : le premier solide apparaît souvent à la première goutte
        e = self.raffiner(volume_max, points, self.depart(volume_max, 1e-9), self.p)
        resultat = []
        for k in np.flatnonzero(e.S[1:] != e.S[:-1]):
            avant, apres = int(e.S[k]), int(e.S[k + 1])
            gauche = e[k : k + 1]

            def prolonge(V: float, gauche: Etat = gauche) -> Etat:
                return self.continuer(Etat(np.array([V]), gauche.pH, gauche.u, gauche.S))[0]

            for i, p in enumerate(self.precipites):
                if (avant ^ apres) >> i & 1:

                    def ecart(v: float, i: int = i) -> float:
                        return self._sursaturation(i, prolonge(v))

                    V = changement_de_signe(ecart, float(e.volume[k]), float(e.volume[k + 1]))
                    resultat.append(
                        Evenement(V, float(-prolonge(V).u[0]), p.partenaire, bool(apres >> i & 1))
                    )
        return tuple(resultat)
