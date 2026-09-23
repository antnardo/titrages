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
local td = T.titrage(T.solution_redox({ T.espece_redox(fer, 0.1, 0) }), T.solution_redox({ T.espece_redox(cr, 1 / 60, 1) }), 10)
proche("dichromate E(10)", T.E_enh(td, 10), 1.26248379, 1e-6)
proche("volume_redox", T.volume_redox(td, T.E_enh(td, 7)), 7, 1e-6)
local v = T.couple_redox({ -0.26, 0.34, 1.00 }, { 1, 1, 1 }, { 0, 2, 2 })
local mn = T.couple_redox({ 1.51 }, { 5 }, { 8 })
local tv = T.titrage(T.solution_redox({ T.espece_redox(v, 0.1, 0) }), T.solution_redox({ T.espece_redox(mn, 0.02, 1) }), 10)
proche("vanadium : 3 sauts", #T.sauts_redox(tv, 40), 3, 0)

print(("%d/%d vérifications passées"):format(total - echecs, total))
os.exit(echecs == 0 and 0 or 1)
