-- Les tests du moteur Lua, hors de TeX : texlua tests/test_titrages.lua
-- (ou lua5.3). Les valeurs de référence sont celles des tests Python.

package.path = "./?.lua;../?.lua;" .. package.path
local T = require("titrages")

local echecs, total = 0, 0
local function proche(nom, obtenu, attendu, tol)
  total = total + 1
  if math.abs(obtenu - attendu) > (tol or 0.01) then
    echecs = echecs + 1
    print(("ÉCHEC %s : %.6g au lieu de %.6g"):format(nom, obtenu, attendu))
  end
end

local function sol(especes, b) return T.solution(especes, b or 0, 14) end
local ethanoique = T.espece({ 4.76 }, 0.1, 0)

proche("eau", T.pH_solution(sol({})), 7)
proche("HCl 0,01", T.pH_solution(sol({}, -0.01)), 2)
proche("CH3COOH 0,1", T.pH_solution(sol({ ethanoique })), 2.88)
proche("NH3 0,1", T.pH_solution(sol({ T.espece({ 9.25 }, 0.1, 1) })), 11.12)
proche("Na2CO3 0,1", T.pH_solution(sol({ T.espece({ 10.33, 6.35 }, 0.1, 2) })), 11.66)
proche("H3PO4 0,1", T.pH_solution(sol({ T.espece({ 2.15, 7.20, 12.35 }, 0.1, 0) })), 1.63)

local x = T.fractions({ 6.35, 10.33 }, 6.35)
proche("fractions au pKa", x[1] - x[2], 0, 1e-3)
local x2 = T.fractions({ 2.0, 2.67, 6.16, 10.26 }, 40)
proche("fractions sans débordement", x2[5], 1, 1e-12)

local t = T.titrage(sol({ ethanoique }), sol({}, 0.1), 10)
proche("pH(0)", T.pH(t, 0), 2.88)
proche("demi-équivalence", T.pH(t, 5), 4.76)
proche("équivalence", T.pH(t, 10), 8.73)
proche("volume(pH(7))", T.volume(t, T.pH(t, 7)), 7, 1e-9)

local h = 1e-6
local fd = (T.volume(t, 5 + h) - T.volume(t, 5 - h)) / (2 * h)
proche("dérivée analytique", T.derivee_volume(t, 5) / fd, 1, 1e-5)

local sauts = T.sauts(t, 20)
proche("un saut", #sauts, 1, 0)
proche("saut en V", sauts[1][1], 10, 1e-3)

local phos = T.titrage(sol({ T.espece({ 2.15, 7.20, 12.35 }, 0.1, 0) }), sol({}, 0.1), 10)
local eq = T.equivalences(phos)
proche("3 équivalences", #eq, 3, 0)
proche("éq. 3", eq[3], 30, 1e-12)
proche("2 sauts seulement", #T.sauts(phos, 35), 2, 0)

local amm = T.titrage(sol({ T.espece({ 9.25 }, 0.1, 1) }), sol({}, -0.1), 10)
proche("éq. ammoniac", T.equivalences(amm)[1], 10, 1e-12)
proche("saut ammoniac", T.sauts(amm, 20)[1][2], 5.27)

-- le réactif fort du bécher excède les sites faibles : titrage en retour
local retour = T.titrage(sol({ T.espece({ 4.76 }, 0.05, 0) }, 0.1), sol({}, -0.1), 10)
proche("retour acide : 2 équivalences", #T.equivalences(retour), 2, 0)
proche("retour acide : éq. 1", T.equivalences(retour)[1], 5, 1e-12)
proche("retour acide : éq. 2", T.equivalences(retour)[2], 10, 1e-12)
local retour_base = T.titrage(sol({ T.espece({ 9.25 }, 0.05, 1) }, -0.1), sol({}, 0.1), 10)
proche("retour base : 2 équivalences", #T.equivalences(retour_base), 2, 0)
proche("retour base : éq. 1", T.equivalences(retour_base)[1], 5, 1e-12)
proche("retour base : éq. 2", T.equivalences(retour_base)[2], 10, 1e-12)
local deja = T.titrage(sol({ T.espece({ 6.35, 10.33 }, 0.1, 0) }, 0.05), sol({}, 0.1), 10)
proche("base déjà présente : éq. 1", T.equivalences(deja)[1], 5, 1e-12)
proche("base déjà présente : éq. 2", T.equivalences(deja)[2], 15, 1e-12)

local tg = T.tangentes(t, 7, 20)
proche("tangentes", tg.V, 10, 0.01)

local c = T.courbe(t, 20, 100)
proche("courbe : points", #c, 100, 0)
proche("courbe : fin", c[100][1], 20, 0)
local dans = 0
for _, p in ipairs(c) do if p[2] > 6 and p[2] < 11 then dans = dans + 1 end end
proche("courbe : saut échantillonné", dans >= 20 and 1 or 0, 1, 0)

local sans = T.titrage(sol({ T.espece({ 6.35, 10.33 }, 0.01, 0) }), sol({}, 0.1), 10, false)
local Y1 = T.moments(T.espece({ 6.35, 10.33 }), 8)
proche("sans dilution", T.volume(sans, 8) * 0.1 / (0.01 * 10), Y1 + (1e-6 - 1e-8) / 0.01, 1e-9)

-- conductimétrie
local ac = T.espece({ 4.76 }, 0.1, 0); ac.z0, ac.lam = 0, { 0, 4.09 }
local tc = T.titrage(T.solution({ ac }, 0, 14), T.solution({}, 0.1, 14), 10)
tc.burette.base = 0.1
proche("σ(2 mL)", T.conductivite(tc, 2), 0.15436848, 1e-7)
proche("σ(15 mL)", T.conductivite(tc, 15), 0.86140002, 1e-7)
local fort = T.titrage(T.solution({}, -0.01, 14), T.solution({}, 0.1, 14), 100)
fort.becher.acide, fort.burette.base = 0.01, 0.1
proche("σ fort à l'équivalence (corrigée)", T.conductivite(fort, 10) * 110 / 100, (5.01 + 7.63) * 0.01, 1e-5)  -- plus H3O+ et HO- : 5,5e-6

-- oxydoréduction
local fer, ce = T.couple_redox({ 0.77 }, { 1 }), T.couple_redox({ 1.72 }, { 1 })
local cr = T.couple_redox({ 1.33 }, { 6 }, { 14 }, { 1, 2 })
local tr = T.titrage(T.solution_redox({ T.espece_redox(fer, 0.1, 0) }), T.solution_redox({ T.espece_redox(ce, 0.1, 1) }), 10)
proche("E demi-équivalence", T.E_enh(tr, 5), 0.77, 1e-6)
proche("E équivalence", T.E_enh(tr, 10), (0.77 + 1.72) / 2, 1e-4)
proche("équivalence rédox", T.equivalences_redox(tr)[1], 10, 1e-12)
proche("saut rédox", T.sauts_redox(tr, 20)[1][1], 10, 1e-4)
proche("tangentes rédox", T.tangentes_redox(tr, 7, 20).V, 10, 1e-6)
local td = T.titrage(T.solution_redox({ T.espece_redox(fer, 0.1, 0) }), T.solution_redox({ T.espece_redox(cr, 1 / 60, 1) }), 10)
proche("dichromate E(10)", T.E_enh(td, 10), 1.26248379, 1e-6)
proche("volume_redox", T.volume_redox(td, T.E_enh(td, 7)), 7, 1e-6)
local v = T.couple_redox({ -0.26, 0.34, 1.00 }, { 1, 1, 1 }, { 0, 2, 2 })
local mn = T.couple_redox({ 1.51 }, { 5 }, { 8 })
local tv = T.titrage(T.solution_redox({ T.espece_redox(v, 0.1, 0) }), T.solution_redox({ T.espece_redox(mn, 0.02, 1) }), 10)
proche("vanadium : 3 sauts", #T.sauts_redox(tv, 40), 3, 0)

-- complexométrie
local edta = { 0.0, 1.5, 2.0, 2.66, 6.16, 10.24 }
local tca = T.titrage_complexation(T.solution_complexation({ { "Ca2+", 0.01 } }, 0),
  T.solution_complexation({}, 0.01), 20, { { "Ca2+", { 10.7 } } }, edta, 10)
proche("pCa(10)", T.pM(tca, 10), 2.4771212476, 1e-9)
proche("pCa(20)", T.pM(tca, 20), 6.2818168487, 1e-9)
proche("volume(pL(7))", T.volume_complexation(tca, T.pL(tca, 7)), 7, 1e-9)
proche("saut complexométrique", T.sauts_complexation(tca, 40)[1][1], 20, 1e-5)
proche("équivalence EDTA", T.equivalences_complexation(tca)[1], 20, 1e-12)
local x = T.fractions_complexes({ 10.7 }, 10.7)
proche("fractions complexes", x[1], 0.5, 1e-12)

-- précipitation
local tag = T.titrage_precipitation(T.solution_precipitation(0, { { "Cl-", 0.1 } }),
  T.solution_precipitation(0.1), 20, { { partenaire = "Cl-", pKs = 9.75 } })
proche("pAg à l'équivalence", T.pP(tag, 20), 9.75 / 2, 1e-12)
proche("pAg + pCl", T.pP(tag, 10) + T.p_partenaire(tag, 10, "Cl-"), 9.75, 1e-9)
proche("volume(pAg(13))", T.volume_precipitation(tag, T.pP(tag, 13)), 13, 1e-9)
proche("saut AgCl", T.sauts_precipitation(tag, 40)[1][1], 20, 1e-4)
local mohr = T.titrage_precipitation(T.solution_precipitation(0, { { "Cl-", 0.05 }, { "CrO4^2-", 0.005 } }),
  T.solution_precipitation(0.05), 20,
  { { partenaire = "Cl-", pKs = 9.75 }, { partenaire = "CrO4^2-", pKs = 11.95, nu = 2, pKa = { 6.5 } } }, 7.5)
local ev = T.evenements(mohr, 30)
proche("Mohr : deux événements", #ev, 2, 0)
proche("Mohr : le chromate apparaît", ev[2][1], 20.0113792425, 1e-8)
proche("Mohr : équivalence du chromate", T.equivalences_precipitation(mohr)[2], 24, 1e-12)

-- le pH calculé (niveau 2)
local libre = T.titrage_complexation(T.solution_complexation({ { "Ca2+", 0.01 } }),
  T.solution_complexation({}, 0.01, 2), 20, { { "Ca2+", { 10.7 } } }, edta)
proche("sans tampon : le pH se calcule", libre.libre and 1 or 0, 1, 0)
proche("sans tampon : pH à 10 mL", T.pH_du_melange(libre, 10), 3.5113151725, 1e-9)
proche("sans tampon : pas de saut", #T.sauts_complexation(libre, 40), 0, 0)
local mg = T.titrage_precipitation(T.solution_precipitation(0.01),
  T.solution_precipitation(0, {}, T.solution({}, 0.1, 14)), 20, { { partenaire = "OH-", pKs = 11.15, mu = 2 } })
proche("Mg(OH)2 : pH à l'équivalence", T.pH_du_melange(mg, 4), 10.383676690, 1e-8)
proche("Mg(OH)2 : équivalence", T.equivalences_precipitation(mg)[1], 4, 1e-12)
proche("Mg(OH)2 : l'hydroxyde apparaît", T.evenements(mg, 8)[1][1], 0.0053234997, 1e-9)
local ok = pcall(T.titrage_precipitation, T.solution_precipitation(0.01), T.solution_precipitation(), 20,
  { { partenaire = "OH-", pKs = 11.15, mu = 2 } }, 10)
proche("hydroxyde à pH imposé : erreur", ok and 0 or 1, 1, 0)

print(("%d/%d vérifications passées"):format(total - echecs, total))
os.exit(echecs == 0 and 0 or 1)
