-- titrages.lua — le moteur de calcul du paquet LaTeX titrages.
--
-- Le même algorithme que le paquet Python titrages, du même dépôt
-- (https://github.com/antnardo/titrages), en Lua 5.3 pur : il tourne sous LuaTeX
-- comme sous un lua ordinaire, ce qui permet de le tester hors de TeX
-- (tests/test_titrages.lua) et de le confronter au moteur Python (tests/test_lua.py).
--
-- À pH imposé, tout est explicite. Pour une solution,
--
--   R(pH) = sum_i C_i (Y1_i - q_i) - b + delta,    delta = w - h,
--
-- est le résidu de l'électroneutralité : nul au pH de la solution, strictement
-- croissant. Pour un titrage, V0 R_becher + V R_burette = 0 donne le volume versé
-- comme une fonction directe du pH :
--
--   V(pH) = -V0 R_becher(pH) / R_burette(pH).
--
-- Les fractions se calculent en échelle logarithmique (log-sum-exp) : la forme
-- naïve x_j = x_0 k_j / h^j déborde pour un polyacide loin de ses pKa.
--
-- La même idée vaut pour l'oxydoréduction (le potentiel), la complexométrie (pL, le
-- ligand libre) et la précipitation (pP, l'ion commun, pour chaque ensemble de solides
-- présents) : voir l'en-tête de chaque section.
--
-- Copyright (c) 2026 antnardo — licence MIT (voir LICENSE).

local M = {}
M.version = "0.1.4"
M.titrages = {}

local LN10 = math.log(10)
local MARGE = 6          -- bornes de la dichotomie : [-6, pKe + 6]
local ITERATIONS = 64

-- --- équilibre ----------------------------------------------------------------

--- Une espèce normalisée : pKa triés, log10 k_j cumulés.
function M.espece(pKa, C, q)
  local p = {}
  for i, v in ipairs(pKa) do p[i] = v end
  table.sort(p)
  local L, s = { [0] = 0 }, 0
  for j = 1, #p do
    s = s - p[j]
    L[j] = s
  end
  q = q or 0
  if q < 0 or q > #p then
    error(("titrages : forme %d hors de 0..%d"):format(q, #p))
  end
  return { pKa = p, L = L, N = #p, C = C or 0, q = q }
end

--- Les moments Y1 et Y2 de la distribution des formes.
local function moments(e, pH)
  local tmax, t = -math.huge, {}
  for j = 0, e.N do
    t[j] = e.L[j] + j * pH
    if t[j] > tmax then tmax = t[j] end
  end
  local S, S1, S2 = 0, 0, 0
  for j = 0, e.N do
    local x = 10 ^ (t[j] - tmax)
    S, S1, S2 = S + x, S1 + j * x, S2 + j * j * x
  end
  return S1 / S, S2 / S
end
M.moments = moments

--- Les fractions x_0 … x_N d'une espèce (indices Lua 1 … N+1).
local function fractions_espece(e, pH)
  local tmax, t = -math.huge, {}
  for j = 0, e.N do
    t[j] = e.L[j] + j * pH
    if t[j] > tmax then tmax = t[j] end
  end
  local S, x = 0, {}
  for j = 0, e.N do
    x[j + 1] = 10 ^ (t[j] - tmax)
    S = S + x[j + 1]
  end
  for j = 1, e.N + 1 do x[j] = x[j] / S end
  return x
end

--- Y1 - q et la variance Y2 - Y1^2, sans soustraire deux nombres voisins :
--- sum_j x_j (j - q) et sum_{j<k} x_j x_k (k - j)^2. Dans un saut, une forme domine
--- à 1 - 1e-10 près, et Y2 - Y1*Y1 n'y garderait que six chiffres.
local function ecart_et_variance(e, pH)
  local x = fractions_espece(e, pH)
  local ecart, var = 0, 0
  for j = 0, e.N do
    ecart = ecart + x[j + 1] * (j - e.q)
    for k = j + 1, e.N do var = var + x[j + 1] * x[k + 1] * (k - j) ^ 2 end
  end
  return ecart, var
end

--- Les fractions x_0 … x_N (indices Lua 1 … N+1).
function M.fractions(pKa, pH)
  local e = M.espece(pKa)
  local tmax, t = -math.huge, {}
  for j = 0, e.N do
    t[j] = e.L[j] + j * pH
    if t[j] > tmax then tmax = t[j] end
  end
  local S, x = 0, {}
  for j = 0, e.N do
    x[j + 1] = 10 ^ (t[j] - tmax)
    S = S + x[j + 1]
  end
  for j = 1, e.N + 1 do x[j] = x[j] / S end
  return x
end

local function delta(pH, pKe) return 10 ^ (pH - pKe) - 10 ^ (-pH) end
local function somme(pH, pKe) return 10 ^ (pH - pKe) + 10 ^ (-pH) end

--- Une solution : { especes = {...}, b = [base forte] - [acide fort], pKe = 14 }.
function M.solution(especes, b, pKe)
  return { especes = especes or {}, b = b or 0, pKe = pKe or 14 }
end

function M.charge(s, pH)
  local r = -s.b
  for _, e in ipairs(s.especes) do
    local ecart = ecart_et_variance(e, pH)
    r = r + e.C * ecart
  end
  return r
end

function M.derivee_charge(s, pH)
  local r = 0
  for _, e in ipairs(s.especes) do
    local _, var = ecart_et_variance(e, pH)
    r = r + e.C * var
  end
  return LN10 * r
end

function M.residu(s, pH) return M.charge(s, pH) + delta(pH, s.pKe) end

--- Le pouvoir tampon dC_base/dpH.
function M.tampon(s, pH) return M.derivee_charge(s, pH) + LN10 * somme(pH, s.pKe) end

--- La racine d'une fonction croissante, par dichotomie.
local function dichotomie(f, a, b)
  for _ = 1, ITERATIONS do
    local m = (a + b) / 2
    if f(m) > 0 then b = m else a = m end
  end
  return (a + b) / 2
end

--- La racine d'une fonction qui change de signe entre a et b.
local function changement(f, a, b)
  local fa = f(a)
  for _ = 1, ITERATIONS do
    local m = (a + b) / 2
    local fm = f(m)
    if (fm > 0) == (fa > 0) then a, fa = m, fm else b = m end
  end
  return (a + b) / 2
end

function M.pH_solution(s)
  return dichotomie(function(p) return M.residu(s, p) end, -MARGE, s.pKe + MARGE)
end

-- --- titrage ------------------------------------------------------------------

--- Un titrage : { becher, burette, V0, dilution = true }.
function M.titrage(becher, burette, V0, dilution)
  if V0 <= 0 then error("titrages : le volume du bécher est positif") end
  if dilution == nil then dilution = true end
  return { becher = becher, burette = burette, V0 = V0, dilution = dilution }
end

local function denominateur(t, p)
  if t.dilution then return M.residu(t.burette, p) end
  return M.charge(t.burette, p)
end

local function derivee_denominateur(t, p)
  if t.dilution then return M.tampon(t.burette, p) end
  return M.derivee_charge(t.burette, p)
end

--- Le volume versé qui amène à ce pH, ou nil.
function M.volume(t, p)
  local V = -t.V0 * M.residu(t.becher, p) / denominateur(t, p)
  if V ~= V or V < 0 or V == math.huge then return nil end
  return V
end

--- dV/dpH = -(V0 N' + V D')/D, en dérivant V0 N + V D = 0 : une somme de termes
--- positifs, qui ne perd aucun chiffre au cœur du saut.
function M.derivee_volume(t, p)
  local D = denominateur(t, p)
  local V = -t.V0 * M.residu(t.becher, p) / D
  return -(t.V0 * M.tampon(t.becher, p) + V * derivee_denominateur(t, p)) / D
end

function M.pente(t, p) return 1 / M.derivee_volume(t, p) end

--- Le pH après avoir versé V.
function M.pH(t, V)
  return dichotomie(function(p)
    return t.V0 * M.residu(t.becher, p) + V * denominateur(t, p)
  end, -MARGE, t.becher.pKe + MARGE)
end

--- La courbe de 0 à Vmax en n points, à longueur d'arc constante, chaque point exact.
function M.courbe(t, Vmax, n)
  n = n or 300
  local pH0, pH1 = M.pH_solution(t.becher), M.pH(t, Vmax)
  local m = 40 * n
  local g, V, s = {}, {}, { 0 }
  local dp = math.max(math.abs(pH1 - pH0), 1e-12)
  for i = 1, m do
    g[i] = pH0 + (pH1 - pH0) * (i - 1) / (m - 1)
    V[i] = M.volume(t, g[i]) or 0
  end
  V[1], V[m] = 0, Vmax
  for i = 2, m do
    local dx, dy = (V[i] - V[i - 1]) / Vmax, (g[i] - g[i - 1]) / dp
    s[i] = s[i - 1] + math.sqrt(dx * dx + dy * dy)
  end
  local points, k = {}, 2
  for i = 1, n do
    local cible = s[m] * (i - 1) / (n - 1)
    while k < m and s[k] < cible do k = k + 1 end
    local a = s[k] > s[k - 1] and (cible - s[k - 1]) / (s[k] - s[k - 1]) or 0
    local p = g[k - 1] + a * (g[k] - g[k - 1])
    local v = M.volume(t, p) or 0
    if i == 1 then v, p = 0, pH0 elseif i == n then v, p = Vmax, pH1 end
    points[i] = { v, p }
  end
  return points
end

--- Le minimum d'une fonction positive et lisse, à la précision de la machine : le zéro
--- de (ln f)', estimée par différence centrée de pas `pas`, trouvé par dichotomie.
--- Une section dorée, qui compare des valeurs de f, plafonnerait à sqrt(epsilon).
local function minimum_lisse(f, a, b, pas)
  local function pente(x) return math.log(f(x + pas)) - math.log(f(x - pas)) end
  if a > b then a, b = b, a end
  if not (pente(a) < 0 and pente(b) > 0) then
    if f(a) <= f(b) then return a end
    return b
  end
  for _ = 1, ITERATIONS do
    local m = (a + b) / 2
    if pente(m) > 0 then b = m else a = m end
  end
  return (a + b) / 2
end

--- Les points d'inflexion : { {V, pH, pente}, ... }.
function M.sauts(t, Vmax, n)
  n = n or 4000
  local pH0, pH1 = M.pH_solution(t.becher), M.pH(t, Vmax)
  local g, d = {}, {}
  local function absd(p) return math.abs(M.derivee_volume(t, p)) end
  for i = 1, n do
    g[i] = pH0 + (pH1 - pH0) * (i - 1) / (n - 1)
    d[i] = absd(g[i])
  end
  local resultat = {}
  for i = 2, n - 1 do
    if d[i] < d[i - 1] and d[i] <= d[i + 1] then
      local p = minimum_lisse(absd, g[i - 1], g[i + 1], 1e-6)
      local V = M.volume(t, p)
      if V and V > 0 and V < Vmax then
        resultat[#resultat + 1] = { V, p, M.pente(t, p) }
      end
    end
  end
  return resultat
end

--- Les sites acides et basiques du bécher, { {pKa, quantité}, ... } chacun, une fois
--- que l'acide ou la base forts qu'il contient ont réagi : un excès de base forte
--- consomme les sites acides par pKa croissant et en fait des sites basiques de même
--- pKa, un excès d'acide fort fait l'inverse ; ce qui en reste est un site à pKa infini.
local function sites_becher(becher, V0)
  local acides, bases = {}, {}
  for _, e in ipairs(becher.especes) do
    local n = e.C * V0
    for j = e.q + 1, e.N do acides[#acides + 1] = { e.pKa[j], n } end
    for j = e.q, 1, -1 do bases[#bases + 1] = { e.pKa[j], n } end
  end
  local b = becher.b * V0
  local source, cible, sens = acides, bases, 1
  if b < 0 then source, cible, sens = bases, acides, -1 end
  local tri = {}
  for i, s in ipairs(source) do tri[i] = s end
  table.sort(tri, function(x, y) return sens * x[1] < sens * y[1] end)
  local reste = math.abs(b)
  for _, s in ipairs(tri) do
    if reste <= 0 then break end
    local pris = math.min(reste, s[2])
    s[2] = s[2] - pris
    reste = reste - pris
    cible[#cible + 1] = { s[1], pris }
  end
  if reste > 0 then cible[#cible + 1] = { sens * math.huge, reste } end
  return acides, bases
end

--- Les équivalences stœchiométriques, pour un titrant fort : les sites du bécher
--- rangés par pKa, croissant pour une base forte titrante, décroissant pour un acide
--- fort ; deux sites de même pKa ne font qu'une équivalence.
function M.equivalences(t)
  if #t.burette.especes > 0 then
    error("titrages : équivalences stœchiométriques définies pour un titrant fort seulement")
  end
  local C = t.burette.b
  if C == 0 then error("titrages : la burette ne contient ni acide ni base") end
  local acides, bases = sites_becher(t.becher, t.V0)
  local sites = (C > 0) and acides or bases
  table.sort(sites, function(x, y)
    if C > 0 then return x[1] < y[1] end
    return x[1] > y[1]
  end)
  local volumes, cumul, precedent = {}, 0, nil
  for _, s in ipairs(sites) do
    if s[2] > 0 then
      cumul = cumul + s[2]
      if s[1] == precedent then
        volumes[#volumes] = cumul / math.abs(C)
      else
        volumes[#volumes + 1] = cumul / math.abs(C)
      end
      precedent = s[1]
    end
  end
  return volumes
end

--- La méthode des tangentes parallèles, depuis V1 avant le saut.
function M.tangentes(t, V1, Vmax)
  local saut
  for _, s in ipairs(M.sauts(t, Vmax)) do
    if s[1] > V1 then saut = s break end
  end
  if not saut then error(("titrages : aucun saut après %g avant %g"):format(V1, Vmax)) end
  local pH1 = M.pH(t, V1)
  local p1 = M.pente(t, pH1)
  local pHmax = M.pH(t, Vmax)
  local function ecart(p) return math.abs(M.pente(t, p)) - math.abs(p1) end
  local m, pH2 = 2000
  local prec = saut[2]
  for i = 1, m do
    local p = saut[2] + (pHmax - saut[2]) * i / m
    if ecart(p) <= 0 then
      pH2 = changement(ecart, prec, p)
      break
    end
    prec = p
  end
  if not pH2 then error("titrages : pas de point de même pente après le saut ; augmenter le volume maximal") end
  local V2 = M.volume(t, pH2)
  local Vm, pHm = (V1 + V2) / 2, (pH1 + pH2) / 2
  local pH = changement(function(p) return (p - pHm) - p1 * (M.volume(t, p) - Vm) end, pH1, pH2)
  return { V1 = V1, pH1 = pH1, V2 = V2, pH2 = pH2, pente = p1, V = M.volume(t, pH), pH = pH }
end

-- --- conductimétrie -------------------------------------------------------------
--
-- Loi de Kohlrausch, solution idéale : sigma = sum lambda_i c_i. À pH imposé, toutes
-- les concentrations sont explicites. Les lambda sont en mS m2 mol-1 et les
-- concentrations en mol/L : la somme est en S/m, l'unité du paquet (la clé unite
-- ne change que l'affichage). Une espèce porte lam (lambda de
-- chaque forme, indices 1..N+1), z0 (charge de sa forme la plus acide) et lamci
-- (lambda de son contre-ion, supposé monovalent) ; une solution, lamA et lamB, ceux
-- des contre-ions de l'acide et de la base forts.

M.LAMBDA_H3O, M.LAMBDA_HO = 34.98, 19.86

--- La part des ions, sans H3O+ ni HO- : linéaire en concentrations.
function M.conductivite_ions(s, pH)
  local total = 0
  for i, e in ipairs(s.especes) do
    local x = fractions_espece(e, pH)
    for j = 0, e.N do
      if (e.z0 or 0) - j ~= 0 then
        local lam = e.lam and e.lam[j + 1]
        if not lam then
          error(("titrages : lambda inconnu pour la forme %d de l'espèce %d (clé lambda)"):format(j, i))
        end
        total = total + lam * e.C * x[j + 1]
      end
    end
    local zq = math.abs((e.z0 or 0) - e.q)
    if zq > 0 then
      if not e.lamci then error(("titrages : lambda du contre-ion de l'espèce %d inconnu (clé contre ion)"):format(i)) end
      total = total + e.lamci * zq * e.C
    end
  end
  return total + (s.acide or 0) * (s.lamA or 7.63) + (s.base or 0) * (s.lamB or 5.01)
end

--- La conductivité du mélange après avoir versé V (S/m), et son pH.
function M.conductivite(t, V)
  local p = M.pH(t, V)
  local ions = (t.V0 * M.conductivite_ions(t.becher, p) + V * M.conductivite_ions(t.burette, p)) / (t.V0 + V)
  return ions + M.LAMBDA_H3O * 10 ^ (-p) + M.LAMBDA_HO * 10 ^ (p - t.becher.pKe), p
end

-- --- oxydoréduction ---------------------------------------------------------------
--
-- Les électrons remplacent les protons, E/s le pH (s = RT ln10/F). Un système a ses
-- formes de la plus réduite (0) à la plus oxydée (N) ; chaque étape j a E0, n et H
-- (protons consommés par la réduction), ce qui donne à pH fixé le potentiel apparent
-- E0 - s H/n pH. Le résidu R(E) = sum T (Y(E) - e_q) est croissant. Un couple
-- dimérique (st = {a, b}, a Ox + n e- = b Red) a une fraction oxydée qui dépend de la
-- concentration : racine d'un trinôme. Le potentiel se cherche donc à volume donné,
-- E(V), par dichotomie sur le résidu du mélange : exact pour tous les couples.

M.NERNST = 0.05916
M.REFERENCES = { ENH = 0, ECS = 0.244, ["Ag/AgCl"] = 0.197 }
local EMIN, EMAX = -3, 4

--- Un couple rédox normalisé.
function M.couple_redox(E0, n, H, st)
  local N = #E0
  n, H, st = n or {}, H or {}, st or { 1, 1 }
  local c = { E0 = E0, n = {}, H = {}, N = N, a = st[1], b = st[2] }
  for j = 1, N do c.n[j], c.H[j] = n[j] or 1, H[j] or 0 end
  if c.a ~= 1 or c.b ~= 1 then
    if N ~= 1 then error("titrages : un couple dimérique n'a qu'une étape") end
    c.unites = (c.a == 1) and { 1, 2 } or { 2, 1 }
    c.e = { [0] = 0, c.n[1] / (c.a * c.unites[2]) }
  else
    c.unites = {}
    c.e = { [0] = 0 }
    for j = 0, N do c.unites[j + 1] = 1 end
    for j = 1, N do c.e[j] = c.e[j - 1] + c.n[j] end
  end
  c.cle = table.concat(E0, ",") .. "|" .. table.concat(c.n, ",") .. "|" .. table.concat(c.H, ",")
    .. "|" .. c.a .. "," .. c.b
  return c
end

local function racine_1_plus(lx)
  lx = math.max(math.min(lx, 600), -600)
  if lx > 0 then return 10 ^ (lx / 2) * math.sqrt(1 + 10 ^ (-lx)) end
  return math.sqrt(1 + 10 ^ lx)
end

--- Y - eref, dY/dE et dY/dT, exactement, à ce potentiel (V/ENH).
--- Monomère : dY/dE = ln10/s * sum_{j<k} x_j x_k (e_k - e_j)^2, et Y - eref =
--- sum_j x_j (e_j - eref) : rien ne s'y soustrait. Dimère : u la forme monomère,
--- D = T - u, T = u + 2 kappa u^2, u = 2T/(1+q), D = T (q-1)/(q+1), q = sqrt(1+8 kappa T) ;
--- dY/dE = e1 n ln10/s u D/(T(2T-u)), dY/dT = ±e1 u D/(T^2(2T-u)).
local function electrons_et_derivees(c, E, T, pH, s, eref)
  eref = eref or 0
  if c.a == 1 and c.b == 1 then
    local tmax, t, L = -math.huge, {}, 0
    for j = 0, c.N do
      if j > 0 then L = L - c.n[j] * (c.E0[j] - s * c.H[j] / c.n[j] * pH) end
      t[j] = (c.e[j] * E + L) / s
      if t[j] > tmax then tmax = t[j] end
    end
    local S, x = 0, {}
    for j = 0, c.N do
      x[j] = 10 ^ (t[j] - tmax)
      S = S + x[j]
    end
    local Y, var = 0, 0
    for j = 0, c.N do
      x[j] = x[j] / S
      Y = Y + x[j] * (c.e[j] - eref)
    end
    for j = 0, c.N do
      for k = j + 1, c.N do var = var + x[j] * x[k] * (c.e[k] - c.e[j]) ^ 2 end
    end
    return Y, LN10 / s * var, 0
  end
  T = math.max(T, 1e-300)
  local n = c.n[1]
  local logK = n * (E - (c.E0[1] - s * c.H[1] / n * pH)) / s
  local lx = math.log(8 * T, 10) + ((c.a == 1) and logK or -logK)
  local q = racine_1_plus(lx)
  local q_1 = (lx < 0) and 10 ^ lx / (q + 1) or q - 1
  local u, D = 2 * T / (1 + q), T * q_1 / (q + 1)
  local petit = math.min(u, D) / T
  local mono, dim = 1 - petit, petit
  if u < D then mono, dim = petit, 1 - petit end
  local ox, red = dim, mono
  if c.a ~= 1 then ox, red = mono, dim end
  local e1 = c.e[1]
  local Y = (eref == e1) and -e1 * red or e1 * ox - eref
  local den = T * (2 * T - u)
  local commun = (den > 0) and e1 * u * D / den or 0
  local signe = (c.a == 1) and 1 or -1
  return Y, n * LN10 / s * commun, signe * commun / T
end

--- Les électrons cédés par unité, à ce potentiel (V/ENH).
function M.electrons_cedes(c, E, T, pH, s)
  return (electrons_et_derivees(c, E, T, pH, s, 0))
end

--- Le bilan d'une solution rédox : par couple, T (unités) et I (électrons introduits).
local function bilan(s)
  local b, ordre = {}, {}
  for _, e in ipairs(s.especes) do
    local c = e.couple
    local T = e.C * c.unites[e.q + 1]
    if not b[c.cle] then
      b[c.cle] = { couple = c, T = 0, I = 0 }
      ordre[#ordre + 1] = c.cle
    end
    b[c.cle].T = b[c.cle].T + T
    b[c.cle].I = b[c.cle].I + T * c.e[e.q]
  end
  return b, ordre
end

function M.solution_redox(especes, pH, s)
  return { especes = especes or {}, pH = pH or 0, s = s or M.NERNST, redox = true }
end

function M.espece_redox(couple, C, q)
  if q < 0 or q > couple.N then error(("titrages : forme %d hors de 0..%d"):format(q, couple.N)) end
  return { couple = couple, C = C, q = q }
end

function M.E_solution(sol)
  local b, ordre = bilan(sol)
  return dichotomie(function(E)
    local r = 0
    for _, k in ipairs(ordre) do
      local x = b[k]
      r = r + x.T * M.electrons_cedes(x.couple, E, x.T, sol.pH, sol.s) - x.I
    end
    return r
  end, EMIN, EMAX)
end

--- Le résidu du mélange, multiplié par V0 + V.
local function residu_redox(t, E, V)
  local bb, ob = bilan(t.becher)
  local bu, ou = bilan(t.burette)
  local cles, vu = {}, {}
  for _, k in ipairs(ob) do cles[#cles + 1], vu[k] = k, true end
  for _, k in ipairs(ou) do if not vu[k] then cles[#cles + 1] = k end end
  local r = 0
  for _, k in ipairs(cles) do
    local x, y = bb[k] or { T = 0, I = 0 }, bu[k] or { T = 0, I = 0 }
    local c = (bb[k] or bu[k]).couple
    local n = t.V0 * x.T + V * y.T
    r = r + n * M.electrons_cedes(c, E, n / (t.V0 + V), t.becher.pH, t.becher.s) - (t.V0 * x.I + V * y.I)
  end
  return r
end

--- Le potentiel (V/ENH) après avoir versé V.
function M.E_enh(t, V)
  return dichotomie(function(E) return residu_redox(t, E, V) end, EMIN, EMAX)
end

--- dg/dE et dg/dV, exactement, g étant le résidu du mélange multiplié par V0 + V.
local function derivees_redox(t, E, V)
  local bb, ob = bilan(t.becher)
  local bu, ou = bilan(t.burette)
  local cles, vu = {}, {}
  for _, k in ipairs(ob) do cles[#cles + 1], vu[k] = k, true end
  for _, k in ipairs(ou) do if not vu[k] then cles[#cles + 1] = k end end
  local gE, gV = 0, 0
  for _, k in ipairs(cles) do
    local x, y = bb[k] or { T = 0, I = 0 }, bu[k] or { T = 0, I = 0 }
    local c = (bb[k] or bu[k]).couple
    local S = t.V0 * x.T + V * y.T
    local T = S / (t.V0 + V)
    local eref = (y.T > 0) and y.I / y.T or 0
    local Yr, dE, dT = electrons_et_derivees(c, E, T, t.becher.pH, t.becher.s, eref)
    gE = gE + S * dE
    gV = gV + y.T * Yr + S * dT * (y.T - T) / (t.V0 + V)
  end
  return gE, gV
end

--- dE/dV = -g_V/g_E, par dérivation implicite exacte.
function M.pente_redox(t, V)
  local gE, gV = derivees_redox(t, M.E_enh(t, V), V)
  return -gV / gE
end

--- Le premier volume où le potentiel est défini (0, ou Vmax/1000 si une seule forme).
function M.volume_initial(t, Vmax)
  local E0 = M.E_enh(t, 0)
  if math.min(E0 - EMIN, EMAX - E0) < 1e-6 then return Vmax * 1e-3 end
  return 0
end

--- La courbe E(V) (V/ENH), à longueur d'arc constante.
function M.courbe_redox(t, Vmax, n)
  n = n or 300
  local debut = M.volume_initial(t, Vmax)
  local m = 20 * n
  local g, E, s = {}, {}, { 0 }
  local Emin, Emax = math.huge, -math.huge
  for i = 1, m do
    g[i] = debut + (Vmax - debut) * (i - 1) / (m - 1)
    E[i] = M.E_enh(t, g[i])
    Emin, Emax = math.min(Emin, E[i]), math.max(Emax, E[i])
  end
  local dE = math.max(Emax - Emin, 1e-12)
  for i = 2, m do
    local dx, dy = (g[i] - g[i - 1]) / Vmax, (E[i] - E[i - 1]) / dE
    s[i] = s[i - 1] + math.sqrt(dx * dx + dy * dy)
  end
  local points, k = {}, 2
  for i = 1, n do
    local cible = s[m] * (i - 1) / (n - 1)
    while k < m and s[k] < cible do k = k + 1 end
    local a = s[k] > s[k - 1] and (cible - s[k - 1]) / (s[k] - s[k - 1]) or 0
    local V = g[k - 1] + a * (g[k] - g[k - 1])
    if i == 1 then V = debut elseif i == n then V = Vmax end
    points[i] = { V, M.E_enh(t, V) }
  end
  return points
end

--- Le volume versé qui amène au potentiel E (V/ENH), par dichotomie en log V, sans contrôle.
local function volume_brut(t, E)
  if not t.sens_redox then
    t.sens_redox = (M.E_enh(t, t.V0) > M.E_enh(t, M.volume_initial(t, t.V0))) and 1 or -1
  end
  local sens = t.sens_redox
  local a, b = math.log(t.V0, 10) - 12, math.log(t.V0, 10) + 6
  for _ = 1, 80 do
    local m = (a + b) / 2
    if sens * residu_redox(t, E, 10 ^ m) < 0 then b = m else a = m end
  end
  return 10 ^ ((a + b) / 2)
end

--- dV/dE (V/ENH) le long de la courbe directe : -g_E/g_V, au volume V(E).
function M.derivee_volume_redox(t, E)
  local V = volume_brut(t, E)
  local gE, gV = derivees_redox(t, E, V)
  return -gE / gV, V
end

--- Les points d'inflexion de E(V) : { {V, E, pente}, ... } (E en V/ENH).
--- Comme en acide-base, les minima de |dV/dE|, fonction lisse du potentiel : en V,
--- un saut tient parfois en 1e-9 mL.
function M.sauts_redox(t, Vmax, n)
  n = n or 2000
  local E0, E1 = M.E_enh(t, M.volume_initial(t, Vmax)), M.E_enh(t, Vmax)
  local g, d = {}, {}
  local function absd(E) return math.abs((M.derivee_volume_redox(t, E))) end
  for i = 1, n do
    g[i] = E0 + (E1 - E0) * (i - 1) / (n - 1)
    d[i] = absd(g[i])
  end
  local resultat = {}
  for i = 2, n - 1 do
    if d[i] < d[i - 1] and d[i] <= d[i + 1] then
      local E = minimum_lisse(absd, g[i - 1], g[i + 1], 1e-6 * t.becher.s)
      local dV, V = M.derivee_volume_redox(t, E)
      if V > 0 and V < Vmax then resultat[#resultat + 1] = { V, E, 1 / dV } end
    end
  end
  return resultat
end

--- La méthode des tangentes parallèles sur E(V), depuis V1 avant le saut (E en V/ENH).
function M.tangentes_redox(t, V1, Vmax)
  local saut
  for _, s in ipairs(M.sauts_redox(t, Vmax)) do
    if s[1] > V1 then saut = s break end
  end
  if not saut then error(("titrages : aucun saut après %g avant %g"):format(V1, Vmax)) end
  local E1 = M.E_enh(t, V1)
  local p1 = M.pente_redox(t, V1)
  local Emax = M.E_enh(t, Vmax)
  local function ecart(E) return math.abs(1 / (M.derivee_volume_redox(t, E))) - math.abs(p1) end
  local m, E2 = 2000
  local prec = saut[2]
  for i = 1, m do
    local E = saut[2] + (Emax - saut[2]) * i / m
    if ecart(E) <= 0 then
      E2 = changement(ecart, prec, E)
      break
    end
    prec = E
  end
  if not E2 then error("titrages : pas de point de même pente après le saut ; augmenter le volume maximal") end
  local V2 = volume_brut(t, E2)
  local Vm, Em = (V1 + V2) / 2, (E1 + E2) / 2
  local E = changement(function(e) return (e - Em) - p1 * (volume_brut(t, e) - Vm) end, E1, E2)
  return { V1 = V1, E1 = E1, V2 = V2, E2 = E2, pente = p1, V = volume_brut(t, E), E = E }
end

--- Les équivalences stœchiométriques : les électrons, rangés par potentiel apparent.
function M.equivalences_redox(t)
  local oxydant = M.E_solution(t.burette) > M.E_solution(t.becher)
  local capacite = 0
  for _, e in ipairs(t.burette.especes) do
    local c, T = e.couple, e.C * e.couple.unites[e.q + 1]
    capacite = capacite + T * (oxydant and c.e[e.q] or (c.e[c.N] - c.e[e.q]))
  end
  if capacite <= 0 then error("titrages : la burette n'apporte aucun électron") end
  local pH, s, sites = t.becher.pH, t.becher.s, {}
  for _, e in ipairs(t.becher.especes) do
    local c, T = e.couple, e.C * e.couple.unites[e.q + 1]
    local a, b, pas = e.q + 1, c.N, 1
    if not oxydant then a, b, pas = e.q, 1, -1 end
    for j = a, b, pas do
      sites[#sites + 1] = { c.E0[j] - s * c.H[j] / c.n[j] * pH, T * (c.e[j] - c.e[j - 1]) }
    end
  end
  table.sort(sites, function(x, y) if oxydant then return x[1] < y[1] end return x[1] > y[1] end)
  local volumes, cumul = {}, 0
  for _, x in ipairs(sites) do
    cumul = cumul + x[2] * t.V0
    volumes[#volumes + 1] = cumul / capacite
  end
  return volumes
end

-- --- complexométrie ---------------------------------------------------------------
--
-- Le ligand libre joue le rôle du proton : à pL imposé (l' = 10^-pL, le ligand libre
-- toutes formes protonées comprises, à pH tamponné), les fractions des complexes
-- successifs d'un métal sont explicites, x_j ~ beta'_j l^j, beta'_j = beta_j alpha^j.
-- Le bilan du ligand, R(pL) = sum_i C_i Y1_i + l - C_L, est linéaire dans les totaux
-- d'une solution : V(pL) = -V0 R_becher / R_burette, comme en acide-base. pM s'en déduit,
-- pM = -log T_M(V) + log sum_j beta'_j l^j.

M.PL_MAX = 80

--- log alpha, la part de la forme la plus basique d'une espèce de pKa donnés.
local function log_alpha(pKa, pH, nom)
  if not pKa or #pKa == 0 then return 0 end
  if not pH then error(("titrages : %s a des acidités : préciser le pH du milieu"):format(nom)) end
  local x = M.fractions(pKa, pH)
  return math.log(x[#x], 10)
end

--- Les fractions x_0 … x_n des complexes (indices 0 … n) et log sum_j beta'_j l^j.
--- lb[0] = 0, lb[j] = log beta'_j.
local function fractions_lb(lb, pL)
  local tmax, t = -math.huge, {}
  for j = 0, #lb do
    t[j] = lb[j] - j * pL
    if t[j] > tmax then tmax = t[j] end
  end
  local S, x = 0, {}
  for j = 0, #lb do
    x[j] = 10 ^ (t[j] - tmax)
    S = S + x[j]
  end
  for j = 0, #lb do x[j] = x[j] / S end
  return x, tmax + math.log(S, 10)
end

--- Y1 et la variance sum_{j<k} x_j x_k (k - j)^2, écrite sans soustraction.
local function moments_lb(lb, pL)
  local x = fractions_lb(lb, pL)
  local Y1, var = 0, 0
  for j = 0, #lb do
    Y1 = Y1 + j * x[j]
    for k = j + 1, #lb do var = var + x[j] * x[k] * (k - j) ^ 2 end
  end
  return Y1, var
end

--- Les fractions du métal libre et de ses complexes (indices Lua 1 … n+1).
function M.fractions_complexes(log_beta, pL)
  local lb = { [0] = 0 }
  for j, b in ipairs(log_beta) do lb[j] = b end
  local x, r = fractions_lb(lb, pL), {}
  for j = 0, #lb do r[j + 1] = x[j] end
  return r
end

--- Une solution : { metaux = { {nom, C}, ... }, ligand = C }.
function M.solution_complexation(metaux, ligand)
  return { metaux = metaux or {}, ligand = ligand or 0 }
end

local function concentration(liste, nom)
  for _, m in ipairs(liste) do if m[1] == nom then return m[2] end end
  return 0
end

--- Un titrage complexométrique. complexes = { {metal, {log beta_1, …}, n}, … } (cumulés ;
--- n, facultatif, le complexe ML_n qui marque l'équivalence, le dernier par défaut),
--- pKa ceux du ligand, pH celui du tampon, metal celui qu'on suit.
function M.titrage_complexation(becher, burette, V0, complexes, pKa, pH, metal)
  if V0 <= 0 then error("titrages : le volume du bécher est positif") end
  local la = log_alpha(pKa, pH, "le ligand")
  local lb, ordre, stoechiometrie = {}, {}, {}
  for _, c in ipairs(complexes) do
    if lb[c[1]] then error(("titrages : %s : deux jeux de complexes"):format(c[1])) end
    if #c[2] == 0 then error(("titrages : %s : aucune constante"):format(c[1])) end
    local t = { [0] = 0 }
    for j, b in ipairs(c[2]) do t[j] = b + j * la end
    lb[c[1]] = t
    ordre[#ordre + 1] = c[1]
    if c[3] then
      if c[3] < 1 or c[3] > #t then error(("titrages : %s : pas de complexe ML%d déclaré"):format(c[1], c[3])) end
      stoechiometrie[c[1]] = c[3]
    end
  end
  for _, s in ipairs({ becher, burette }) do
    for _, m in ipairs(s.metaux) do
      if not lb[m[1]] then error(("titrages : %s : aucun complexe déclaré"):format(m[1])) end
    end
  end
  metal = metal or ordre[1]
  if not lb[metal] then error(("titrages : %s : aucun complexe déclaré"):format(metal)) end
  return { type = "complexation", becher = becher, burette = burette, V0 = V0, lb = lb, metal = metal,
           stoechiometrie = stoechiometrie }
end

local function residu_c(t, s, pL)
  local r = 10 ^ (-pL) - s.ligand
  for _, m in ipairs(s.metaux) do r = r + m[2] * (moments_lb(t.lb[m[1]], pL)) end
  return r
end

local function derivee_residu_c(t, s, pL)
  local r = 10 ^ (-pL)
  for _, m in ipairs(s.metaux) do
    local _, var = moments_lb(t.lb[m[1]], pL)
    r = r + m[2] * var
  end
  return -LN10 * r
end

local function volume_c(t, pL)
  return -t.V0 * residu_c(t, t.becher, pL) / residu_c(t, t.burette, pL)
end

--- Le volume versé qui amène le ligand libre à ce pL, ou nil.
function M.volume_complexation(t, pL)
  local V = volume_c(t, pL)
  if V ~= V or V < 0 or V == math.huge then return nil end
  return V
end

--- dV/dpL = -(V0 N' + V D')/D.
function M.derivee_volume_complexation(t, pL)
  local D = residu_c(t, t.burette, pL)
  local V = -t.V0 * residu_c(t, t.becher, pL) / D
  return -(t.V0 * derivee_residu_c(t, t.becher, pL) + V * derivee_residu_c(t, t.burette, pL)) / D
end

--- Le pL après avoir versé V : V0 R_b + V R_bu décroît avec pL.
function M.pL(t, V)
  local haut = 0
  for _, s in ipairs({ t.becher, t.burette }) do
    haut = math.max(haut, s.ligand)
    for _, m in ipairs(s.metaux) do haut = math.max(haut, m[2]) end
  end
  haut = math.log(haut + 1e-300, 10) + 1
  return -dichotomie(function(u)
    return t.V0 * residu_c(t, t.becher, -u) + V * residu_c(t, t.burette, -u)
  end, -M.PL_MAX, haut)
end

local function totale_c(t, metal, V)
  return (t.V0 * concentration(t.becher.metaux, metal) + V * concentration(t.burette.metaux, metal)) / (t.V0 + V)
end

local function pM_c(t, pL, V, metal)
  local _, ls = fractions_lb(t.lb[metal], pL)
  return -math.log(totale_c(t, metal, V), 10) + ls
end

--- pM du métal (par défaut, celui du titrage) après avoir versé V.
function M.pM(t, V, metal)
  return pM_c(t, M.pL(t, V), V, metal or t.metal)
end

--- dpM/dpL le long de la courbe : -(Cu - T)/((V0+V) T ln10) V' - Y1.
local function derivee_pM_c(t, pL, metal)
  local V = volume_c(t, pL)
  local dV = M.derivee_volume_complexation(t, pL)
  local T = totale_c(t, metal, V)
  local Cu = concentration(t.burette.metaux, metal)
  return -(Cu - T) / ((t.V0 + V) * T * LN10) * dV - (moments_lb(t.lb[metal], pL))
end

--- dpM/dV au point de la courbe de ce pL.
function M.pente_complexation(t, pL, metal)
  return derivee_pM_c(t, pL, metal or t.metal) / M.derivee_volume_complexation(t, pL)
end

--- La courbe { {V, pM}, … } à longueur d'arc constante, depuis 1e-9 Vmax : à volume nul,
--- le bécher n'a pas de ligand libre.
function M.courbe_complexation(t, Vmax, n, metal)
  n, metal = n or 300, metal or t.metal
  local p0, p1 = M.pL(t, 1e-9 * Vmax), M.pL(t, Vmax)
  local m = 40 * n
  local g, V, y = {}, {}, {}
  local ymin, ymax = math.huge, -math.huge
  for i = 1, m do
    g[i] = p0 + (p1 - p0) * (i - 1) / (m - 1)
    V[i] = M.volume_complexation(t, g[i]) or 0
    y[i] = pM_c(t, g[i], V[i], metal)
    if y[i] == y[i] and math.abs(y[i]) < math.huge then
      ymin, ymax = math.min(ymin, y[i]), math.max(ymax, y[i])
    end
  end
  local dy, s = math.max(ymax - ymin, 1e-12), { 0 }
  for i = 2, m do
    local a, b = (V[i] - V[i - 1]) / Vmax, (y[i] - y[i - 1]) / dy
    if b ~= b then b = 0 end
    s[i] = s[i - 1] + math.sqrt(a * a + b * b)
  end
  local points, k = {}, 2
  for i = 1, n do
    local cible = s[m] * (i - 1) / (n - 1)
    while k < m and s[k] < cible do k = k + 1 end
    local a = s[k] > s[k - 1] and (cible - s[k - 1]) / (s[k] - s[k - 1]) or 0
    local p = g[k - 1] + a * (g[k] - g[k - 1])
    local v = M.volume_complexation(t, p) or 0
    points[i] = { v, pM_c(t, p, v, metal), p }
  end
  return points
end

--- Les points d'inflexion { {V, pM, pente}, … } : minima de |dV/dpM|, cherchés en pL
--- depuis 1e-6 Vmax (avant, la courbe est plate à la précision de la machine).
function M.sauts_complexation(t, Vmax, n, metal)
  n, metal = n or 4000, metal or t.metal
  local p0, p1 = M.pL(t, 1e-6 * Vmax), M.pL(t, Vmax)
  local function f(p) return math.abs(M.derivee_volume_complexation(t, p) / derivee_pM_c(t, p, metal)) end
  local g, d = {}, {}
  for i = 1, n do
    g[i] = p0 + (p1 - p0) * (i - 1) / (n - 1)
    d[i] = f(g[i])
  end
  local resultat = {}
  for i = 2, n - 1 do
    if d[i] < d[i - 1] and d[i] <= d[i + 1] then
      local p = minimum_lisse(f, g[i - 1], g[i + 1], 1e-6)
      local V = M.volume_complexation(t, p)
      if V and V > 0 and V < Vmax then
        resultat[#resultat + 1] = { V, pM_c(t, p, V, metal), M.pente_complexation(t, p, metal) }
      end
    end
  end
  return resultat
end

local function n_de(t, metal) return t.stoechiometrie[metal] or #t.lb[metal] end

--- Les équivalences : le ligand titre les métaux du bécher, rangés par log beta_n/n
--- décroissant, chacun fixant n ligands ; ou un métal seul titre le ligand du bécher.
function M.equivalences_complexation(t)
  local b, u = t.becher, t.burette
  if u.ligand > 0 and #u.metaux == 0 and b.ligand == 0 then
    local ordre = {}
    for i, m in ipairs(b.metaux) do ordre[i] = m end
    local function stabilite(m) return t.lb[m[1]][n_de(t, m[1])] / n_de(t, m[1]) end
    table.sort(ordre, function(x, y) return stabilite(x) > stabilite(y) end)
    local volumes, cumul = {}, 0
    for _, m in ipairs(ordre) do
      cumul = cumul + n_de(t, m[1]) * m[2] * t.V0
      volumes[#volumes + 1] = cumul / u.ligand
    end
    return volumes
  end
  if #u.metaux == 1 and u.ligand == 0 and b.ligand > 0 and #b.metaux == 0 then
    local m = u.metaux[1]
    return { t.V0 * b.ligand / (n_de(t, m[1]) * m[2]) }
  end
  error("titrages : équivalences définies pour un ligand qui titre des métaux, ou un métal seul qui titre le ligand")
end

-- --- précipitation -------------------------------------------------------------------
--
-- Un ion commun P (Ag+) forme avec des partenaires X_i des solides P_nu X_mu,
-- Ks = a^nu x^mu, a = [P]. À pP imposé, un solide présent fixe son partenaire,
-- x_sat = (Ks/a^nu)^(1/mu)/alpha (alpha : la forme la plus basique, à pH tamponné).
-- Pour un ensemble S de solides présents, R_S = a + sum_{i in S} r_i (C_Xi - x_sat,i) - C_P,
-- r = nu/mu, est linéaire dans les totaux : V = -V0 R_S,becher / R_S,burette. On essaie
-- les 2^n ensembles et l'on garde celui qui est cohérent. Là où un solide apparaît ou
-- disparaît, la courbe a un point anguleux : un événement, distinct des sauts.

M.P_MAX = 60

--- Une solution : { commun = C, partenaires = { {nom, C}, … } }.
function M.solution_precipitation(commun, partenaires)
  return { commun = commun or 0, partenaires = partenaires or {} }
end

--- Un titrage par précipitation. precipites = { {partenaire=, pKs=, nu=, mu=, pKa=}, … }.
function M.titrage_precipitation(becher, burette, V0, precipites, pH)
  if V0 <= 0 then error("titrages : le volume du bécher est positif") end
  local P, vus = {}, {}
  for i, p in ipairs(precipites) do
    if vus[p.partenaire] then error("titrages : un partenaire ne forme qu'un précipité") end
    vus[p.partenaire] = true
    local nu, mu = p.nu or 1, p.mu or 1
    if nu < 1 or mu < 1 then error(("titrages : %s : stœchiométrie positive"):format(p.partenaire)) end
    P[i] = { partenaire = p.partenaire, pKs = p.pKs, nu = nu, mu = mu, r = nu / mu,
             la = log_alpha(p.pKa, pH, p.partenaire) }
  end
  for _, s in ipairs({ becher, burette }) do
    for _, x in ipairs(s.partenaires) do
      if not vus[x[1]] then error(("titrages : %s : aucun précipité déclaré"):format(x[1])) end
    end
  end
  return { type = "precipitation", becher = becher, burette = burette, V0 = V0, precipites = P }
end

local function x_sat(p, pP) return 10 ^ ((-p.pKs + p.nu * pP) / p.mu - p.la) end

local function present(S, i) return (S >> (i - 1)) & 1 == 1 end

local function residu_p(t, s, pP, S)
  local r = 10 ^ (-pP) - s.commun
  for i, p in ipairs(t.precipites) do
    if present(S, i) then r = r + p.r * (concentration(s.partenaires, p.partenaire) - x_sat(p, pP)) end
  end
  return r
end

--- dR_S/dpP = -ln10 (a + sum_{i in S} r_i^2 x_sat,i), la même pour toute solution.
local function derivee_residu_p(t, pP, S)
  local r = 10 ^ (-pP)
  for i, p in ipairs(t.precipites) do
    if present(S, i) then r = r + p.r ^ 2 * x_sat(p, pP) end
  end
  return -LN10 * r
end

local function totale_p(t, nom, V)
  return (t.V0 * concentration(t.becher.partenaires, nom) + V * concentration(t.burette.partenaires, nom)) / (t.V0 + V)
end

local function volume_S(t, pP, S)
  return -t.V0 * residu_p(t, t.becher, pP, S) / residu_p(t, t.burette, pP, S)
end

local function coherent(t, pP, V, S)
  if V ~= V or V < 0 or V == math.huge then return false end
  for i, p in ipairs(t.precipites) do
    local T, x = totale_p(t, p.partenaire, V), x_sat(p, pP)
    local marge = 1e-12 * math.max(T, x)
    if present(S, i) then
      if T - x < -marge then return false end
    elseif T - x > marge then
      return false
    end
  end
  return true
end

--- Le volume et l'ensemble des solides présents en ce pP ; nil et -1 hors d'atteinte.
local function resoudre(t, pP)
  for S = 0, (1 << #t.precipites) - 1 do
    local V = volume_S(t, pP, S)
    -- au départ, V vaut zéro à l'arrondi près, et peut sortir à -1e-16 : c'est zéro
    if V < 0 and V > -1e-12 * t.V0 then V = 0 end
    if coherent(t, pP, V, S) then return V, S end
  end
  return nil, -1
end

--- Le volume versé qui amène l'ion commun à ce pP, ou nil.
function M.volume_precipitation(t, pP) return (resoudre(t, pP)) end

--- Les partenaires dont le solide est présent au point de ce pP.
function M.solides(t, pP)
  local _, S = resoudre(t, pP)
  local noms = {}
  for i, p in ipairs(t.precipites) do
    if S >= 0 and present(S, i) then noms[#noms + 1] = p.partenaire end
  end
  return noms
end

local function derivee_volume_S(t, pP, S)
  local D = residu_p(t, t.burette, pP, S)
  return -(t.V0 + volume_S(t, pP, S)) * derivee_residu_p(t, pP, S) / D
end

--- dV/dpP sur le morceau de courbe où tombe ce pP.
function M.derivee_volume_precipitation(t, pP)
  local _, S = resoudre(t, pP)
  if S < 0 then return 0 / 0 end
  return derivee_volume_S(t, pP, S)
end

function M.pente_precipitation(t, pP) return 1 / M.derivee_volume_precipitation(t, pP) end

local function haut_p(t)
  local h = 0
  for _, s in ipairs({ t.becher, t.burette }) do
    h = math.max(h, s.commun)
    for _, x in ipairs(s.partenaires) do h = math.max(h, x[2]) end
  end
  return math.log(h + 1e-300, 10) + 1
end

--- Le pP de l'ion commun après avoir versé V : a + sum r max(0, T_X - x_sat) - T_P croît avec a.
function M.pP(t, V)
  local TP = (t.V0 * t.becher.commun + V * t.burette.commun) / (t.V0 + V)
  local TX = {}
  for i, p in ipairs(t.precipites) do TX[i] = totale_p(t, p.partenaire, V) end
  return -dichotomie(function(u)
    local r = 10 ^ u - TP
    for i, p in ipairs(t.precipites) do r = r + p.r * math.max(0, TX[i] - x_sat(p, -u)) end
    return r
  end, -M.P_MAX, haut_p(t))
end

--- Le cologarithme d'un partenaire en solution, toutes formes protonées comprises.
function M.p_partenaire(t, V, nom)
  for _, p in ipairs(t.precipites) do
    if p.partenaire == nom then
      return -math.log(math.min(totale_p(t, nom, V), x_sat(p, M.pP(t, V))), 10)
    end
  end
  error(("titrages : %s n'est pas un partenaire du titrage"):format(nom))
end

local function bornes_p(t, Vmax, depart)
  local p0 = M.pP(t, 0)
  if p0 > M.P_MAX - 1 then p0 = M.pP(t, depart * Vmax) end
  return p0, M.pP(t, Vmax)
end

--- La courbe { {V, pP}, … } à longueur d'arc constante.
function M.courbe_precipitation(t, Vmax, n)
  n = n or 300
  local p0, p1 = bornes_p(t, Vmax, 1e-6)
  local m = 40 * n
  local g, V, s = {}, {}, { 0 }
  local dp = math.max(math.abs(p1 - p0), 1e-12)
  for i = 1, m do
    g[i] = p0 + (p1 - p0) * (i - 1) / (m - 1)
    V[i] = M.volume_precipitation(t, g[i]) or 0
  end
  for i = 2, m do
    local a, b = (V[i] - V[i - 1]) / Vmax, (g[i] - g[i - 1]) / dp
    s[i] = s[i - 1] + math.sqrt(a * a + b * b)
  end
  local points, k = {}, 2
  for i = 1, n do
    local cible = s[m] * (i - 1) / (n - 1)
    while k < m and s[k] < cible do k = k + 1 end
    local a = s[k] > s[k - 1] and (cible - s[k - 1]) / (s[k] - s[k - 1]) or 0
    local p = g[k - 1] + a * (g[k] - g[k - 1])
    points[i] = { M.volume_precipitation(t, p) or 0, p }
  end
  return points
end

--- Les apparitions et disparitions de solides : { {V, pP, partenaire, apparition}, … }.
function M.evenements(t, Vmax, n)
  n = n or 4000
  local p0, p1 = bornes_p(t, Vmax, 1e-9)
  local g, S = {}, {}
  for i = 1, n do
    g[i] = p0 + (p1 - p0) * (i - 1) / (n - 1)
    local _, s = resoudre(t, g[i])
    S[i] = s
  end
  local resultat = {}
  for k = 1, n - 1 do
    local avant, apres = S[k], S[k + 1]
    if avant ~= apres and avant >= 0 and apres >= 0 then
      for i, p in ipairs(t.precipites) do
        if present(avant ~ apres, i) then
          local function sursaturation(q)
            local V = resoudre(t, q)
            return totale_p(t, p.partenaire, V or 0 / 0) - x_sat(p, q)
          end
          local q = changement(sursaturation, g[k], g[k + 1])
          resultat[#resultat + 1] = { M.volume_precipitation(t, q), q, p.partenaire, present(apres, i) }
        end
      end
    end
  end
  return resultat
end

--- Les points d'inflexion { {V, pP, pente}, … }, à l'intérieur des morceaux lisses.
function M.sauts_precipitation(t, Vmax, n)
  n = n or 4000
  local p0, p1 = bornes_p(t, Vmax, 1e-6)
  local g, S, d = {}, {}, {}
  for i = 1, n do
    g[i] = p0 + (p1 - p0) * (i - 1) / (n - 1)
    local _, s = resoudre(t, g[i])
    S[i] = s
    d[i] = (s >= 0) and math.abs(derivee_volume_S(t, g[i], s)) or 0 / 0
  end
  local resultat = {}
  for i = 2, n - 1 do
    local s = S[i]
    if s >= 0 and S[i - 1] == s and S[i + 1] == s and d[i] < d[i - 1] and d[i] <= d[i + 1] then
      local q = minimum_lisse(function(x) return math.abs(derivee_volume_S(t, x, s)) end, g[i - 1], g[i + 1], 1e-6)
      local V = volume_S(t, q, s)
      if V == V and V > 0 and V < Vmax then
        resultat[#resultat + 1] = { V, q, 1 / derivee_volume_S(t, q, s) }
      end
    end
  end
  return resultat
end

--- Les équivalences : l'ion commun titre les partenaires du bécher, du premier qui
--- précipite au dernier ; ou un partenaire seul titre l'ion commun restant (Volhard).
function M.equivalences_precipitation(t)
  local b, u = t.becher, t.burette
  if u.commun > 0 and #u.partenaires == 0 and b.commun == 0 then
    local presents = {}
    for _, p in ipairs(t.precipites) do
      local C = concentration(b.partenaires, p.partenaire)
      if C > 0 then
        presents[#presents + 1] = { p = p, C = C, apparition = (p.mu * (math.log(C, 10) + p.la) + p.pKs) / p.nu }
      end
    end
    table.sort(presents, function(x, y) return x.apparition > y.apparition end)
    local volumes, cumul = {}, 0
    for _, x in ipairs(presents) do
      cumul = cumul + x.p.r * x.C * t.V0
      volumes[#volumes + 1] = cumul / u.commun
    end
    return volumes
  end
  if #u.partenaires == 1 and u.commun == 0 and b.commun > 0 then
    local nom, C = u.partenaires[1][1], u.partenaires[1][2]
    local reste, r = b.commun, nil
    for _, p in ipairs(t.precipites) do
      reste = reste - p.r * concentration(b.partenaires, p.partenaire)
      if p.partenaire == nom then r = p.r end
    end
    if reste > 0 then return { t.V0 * reste / (r * C) } end
  end
  error("titrages : équivalences définies pour l'ion commun qui titre des partenaires, ou un partenaire seul qui titre l'ion commun restant")
end

-- --- l'interface avec TeX -----------------------------------------------------
--
-- Le paquet .sty appelle ces fonctions. Elles ne touchent à tex que dans sortir(),
-- ce qui laisse le reste testable hors de TeX.

M.reglages = { chiffres = 2, pHmin = 0, pHmax = 14 }

local function nombre(x) return ("%.6g"):format(x) end

local function sortir(texte)
  if tex and tex.sprint then tex.sprint(texte) else io.write(texte, "\n") end
end
M.sortir = sortir

--- Enregistre un titrage décrit par le .sty.
function M.definir(nom, spec)
  local function sol(d)
    local especes = {}
    for _, e in ipairs(d.especes or {}) do
      local x = M.espece(e.pKa, e.C, e.q)
      x.z0, x.lamci = e.z0 or 0, e.lamci
      if e.lam and #e.lam > 0 then x.lam = e.lam end
      especes[#especes + 1] = x
    end
    local s = M.solution(especes, (d.base or 0) - (d.acide or 0), spec.pKe)
    s.acide, s.base, s.lamA, s.lamB = d.acide or 0, d.base or 0, d.lamA, d.lamB
    return s
  end
  local t = M.titrage(sol(spec.becher or {}), sol(spec.burette or {}), spec.V0 or 10, spec.dilution)
  t.points = spec.points or 300
  t.corrigee = spec.corrigee
  t.unite = spec.unite or "S/m"
  t.abscisse = spec.abscisse or "volume"
  t.echelle = spec.echelle or 1
  if spec.Vmax and spec.Vmax > 0 then
    t.Vmax = spec.Vmax
  else
    local ok, eq = pcall(M.equivalences, t)
    t.Vmax = (ok and #eq > 0) and (eq[#eq] + eq[1]) or 2 * t.V0
  end
  M.titrages[nom] = t
  return t
end

--- Enregistre un titrage potentiométrique.
function M.definir_redox(nom, spec)
  local s = spec.pente or M.NERNST
  local function sol(d)
    local especes = {}
    for _, e in ipairs(d.especes or {}) do
      local c = M.couple_redox(e.E0, e.n, e.H, e.st)
      local q = e.q
      if q == "ox" then q = c.N elseif q == "red" or q == nil then q = 0 end
      especes[#especes + 1] = M.espece_redox(c, e.C, q)
    end
    return M.solution_redox(especes, spec.pH or 0, s)
  end
  local t = M.titrage(sol(spec.becher or {}), sol(spec.burette or {}), spec.V0 or 10)
  t.redox = true
  t.ref = type(spec.reference) == "number" and spec.reference or M.REFERENCES[spec.reference or "ENH"]
  if not t.ref then error("titrages : référence inconnue « " .. tostring(spec.reference) .. " »") end
  t.points = spec.points or 300
  t.abscisse = spec.abscisse or "volume"
  t.echelle = spec.echelle or 1
  if spec.Vmax and spec.Vmax > 0 then
    t.Vmax = spec.Vmax
  else
    local ok, eq = pcall(M.equivalences_redox, t)
    t.Vmax = (ok and #eq > 0) and (eq[#eq] + eq[1]) or 2 * t.V0
  end
  M.titrages[nom] = t
  return t
end

--- Les réglages communs à tous les titrages : points, abscisse, volume maximal.
local function reglages(t, spec, equivalences)
  t.points = spec.points or 300
  t.abscisse = spec.abscisse or "volume"
  t.echelle = spec.echelle or 1
  if spec.Vmax and spec.Vmax > 0 then
    t.Vmax = spec.Vmax
  else
    local ok, eq = pcall(equivalences, t)
    t.Vmax = (ok and #eq > 0) and (eq[#eq] + eq[1]) or 2 * t.V0
  end
end

--- Enregistre un titrage complexométrique.
function M.definir_complexation(nom, spec)
  local function sol(d) return M.solution_complexation(d.metaux, d.ligand) end
  local t = M.titrage_complexation(sol(spec.becher or {}), sol(spec.burette or {}), spec.V0 or 10,
    spec.complexes or {}, spec.pKa, spec.pH, spec.suivi)
  t.suivi = t.metal
  reglages(t, spec, M.equivalences_complexation)
  M.titrages[nom] = t
  return t
end

--- Enregistre un titrage par précipitation.
function M.definir_precipitation(nom, spec)
  local function sol(d) return M.solution_precipitation(d.commun, d.partenaires) end
  local t = M.titrage_precipitation(sol(spec.becher or {}), sol(spec.burette or {}), spec.V0 or 10,
    spec.precipites or {}, spec.pH)
  t.commun = spec.commun or "Ag+"
  if spec.suivi and spec.suivi ~= t.commun then
    M.p_partenaire(t, 0, spec.suivi) -- vérifie que c'est un partenaire
    t.suivi = spec.suivi
  end
  reglages(t, spec, M.equivalences_precipitation)
  M.titrages[nom] = t
  return t
end

--- La grandeur mesurée après avoir versé V : pM du métal suivi, pP de l'ion commun, ou pX
--- d'un partenaire. espece remplace celle du titrage.
local function mesure(t, V, espece)
  if t.type == "complexation" then return M.pM(t, V, espece or t.suivi) end
  espece = espece or t.suivi
  if espece and espece ~= t.commun then return M.p_partenaire(t, V, espece) end
  return M.pP(t, V)
end

--- d pX/dV d'un partenaire : -nu/mu d pP/dV au contact de son solide, la dilution sinon.
local function pente_partenaire(t, V, pP, nom)
  for _, p in ipairs(t.precipites) do
    if p.partenaire == nom then
      local T = totale_p(t, nom, V)
      if T >= x_sat(p, pP) then
        return -(concentration(t.burette.partenaires, nom) - T) / ((t.V0 + V) * T * LN10)
      end
      return -p.nu / p.mu * M.pente_precipitation(t, pP)
    end
  end
end

local EQUIVALENCES = {
  complexation = M.equivalences_complexation,
  precipitation = M.equivalences_precipitation,
}

local function trouver(nom)
  local t = M.titrages[nom]
  if not t then error(("titrages : aucun titrage « %s » ; le définir par \\titragedefinir"):format(nom)) end
  return t
end

--- Le diviseur des volumes : la clé echelle, ou le premier volume équivalent (abscisse n).
local function echelle(t)
  if t.abscisse ~= "n" then return t.echelle or 1 end
  local eq = t.redox and M.equivalences_redox(t) or (EQUIVALENCES[t.type] or M.equivalences)(t)
  if #eq == 0 then error("titrages : abscisse n sans équivalence") end
  return eq[1]
end

local function coordonnees(liste)
  local morceaux = {}
  for i, c in ipairs(liste) do morceaux[i] = "(" .. nombre(c[1]) .. "," .. nombre(c[2]) .. ")" end
  return table.concat(morceaux, " ")
end

-- \addplot+ pour les courbes : le style du cycle de pgfplots s'applique, les options
-- s'y ajoutent ; \addplot seul pour les constructions, qui n'ont pas à le suivre.
local function addplot(options, liste, cycle)
  local commande = cycle and "\\addplot+[no markers," or "\\addplot["
  sortir(commande .. options .. "] coordinates {" .. coordonnees(liste) .. "};")
end

--- \addplottitrage
function M.tex_courbe(nom, options)
  local t = trouver(nom)
  local e, liste = echelle(t), {}
  if t.redox then
    for i, c in ipairs(M.courbe_redox(t, t.Vmax, t.points)) do liste[i] = { c[1] / e, c[2] - t.ref } end
  elseif t.type == "complexation" then
    for i, c in ipairs(M.courbe_complexation(t, t.Vmax, t.points, t.suivi)) do liste[i] = { c[1] / e, c[2] } end
  elseif t.type == "precipitation" then
    for i, c in ipairs(M.courbe_precipitation(t, t.Vmax, t.points)) do
      liste[i] = { c[1] / e, t.suivi and M.p_partenaire(t, c[1], t.suivi) or c[2] }
    end
  else
    for i, c in ipairs(M.courbe(t, t.Vmax, t.points)) do liste[i] = { c[1] / e, c[2] } end
  end
  addplot(options, liste, true)
end

local FACTEURS = { ["S/m"] = 1, ["mS/cm"] = 10, ["mS/m"] = 1000 }

--- \addplottitrageconductivite : sigma(V), à pas de volume constant.
function M.tex_conductivite(nom, options)
  local t = trouver(nom)
  if t.redox or t.type then error("titrages : la conductimétrie n'est prise en charge qu'en acide-base") end
  local e, liste, f = echelle(t), {}, FACTEURS[t.unite] or 1
  local n = math.max(t.points // 2, 50)
  for i = 1, n do
    local V = t.Vmax * (i - 1) / (n - 1)
    local sigma = M.conductivite(t, V)
    if t.corrigee then sigma = sigma * (t.V0 + V) / t.V0 end
    liste[i] = { V / e, f * sigma }
  end
  addplot(options, liste, true)
end

--- \addplottitragederivee
function M.tex_derivee(nom, options)
  local t = trouver(nom)
  local e, liste = echelle(t), {}
  if t.redox then
    for i, c in ipairs(M.courbe_redox(t, t.Vmax, t.points)) do liste[i] = { c[1] / e, M.pente_redox(t, c[1]) * e } end
    return addplot(options, liste, true)
  end
  if t.type == "complexation" then
    for i, c in ipairs(M.courbe_complexation(t, t.Vmax, t.points, t.suivi)) do
      liste[i] = { c[1] / e, M.pente_complexation(t, c[3], t.suivi) * e }
    end
    return addplot(options, liste, true)
  end
  if t.type == "precipitation" then
    for i, c in ipairs(M.courbe_precipitation(t, t.Vmax, t.points)) do
      local d = t.suivi and pente_partenaire(t, c[1], c[2], t.suivi) or M.pente_precipitation(t, c[2])
      liste[i] = { c[1] / e, d * e }
    end
    return addplot(options, liste, true)
  end
  for i, c in ipairs(M.courbe(t, t.Vmax, t.points)) do
    liste[i] = { c[1] / e, M.pente(t, c[2]) * e }
  end
  addplot(options, liste, true)
end

--- \addplottitragesauts
function M.tex_sauts(nom, options)
  local t = trouver(nom)
  local e, liste = echelle(t), {}
  if t.redox then
    for i, s in ipairs(M.sauts_redox(t, t.Vmax)) do liste[i] = { s[1] / e, s[2] - t.ref } end
  elseif t.type == "complexation" then
    for i, s in ipairs(M.sauts_complexation(t, t.Vmax, nil, t.suivi)) do liste[i] = { s[1] / e, s[2] } end
  elseif t.type == "precipitation" then
    for i, s in ipairs(M.sauts_precipitation(t, t.Vmax)) do liste[i] = { s[1] / e, mesure(t, s[1]) } end
  else
    for i, s in ipairs(M.sauts(t, t.Vmax)) do liste[i] = { s[1] / e, s[2] } end
  end
  addplot("only marks," .. options, liste)
end

--- \addplottitrageevenements : les points anguleux d'un titrage par précipitation.
function M.tex_evenements(nom, options)
  local t = trouver(nom)
  if t.type ~= "precipitation" then error("titrages : « " .. nom .. " » n'est pas un titrage par précipitation") end
  local e, liste = echelle(t), {}
  for i, v in ipairs(M.evenements(t, t.Vmax)) do liste[i] = { v[1] / e, mesure(t, v[1]) } end
  addplot("only marks," .. options, liste)
end

--- \addplottitragetangentes : en pH ou en E (par rapport à la référence).
function M.tex_tangentes(nom, V1, options)
  local t = trouver(nom)
  if t.type then error("titrages : la méthode des tangentes n'est prise en charge qu'en pH-métrie et en potentiométrie") end
  local e = echelle(t)
  local r, y1, y2, y
  if t.redox then
    r = M.tangentes_redox(t, V1 * e, t.Vmax)
    y1, y2, y = r.E1 - t.ref, r.E2 - t.ref, r.E - t.ref
  else
    r = M.tangentes(t, V1 * e, t.Vmax)
    y1, y2, y = r.pH1, r.pH2, r.pH
  end
  local marge = 0.15 * (r.V2 - r.V1)
  local a, b = r.V1 - marge, r.V2 + marge
  local Vm, ym = (r.V1 + r.V2) / 2, (y1 + y2) / 2
  for _, P in ipairs({ { r.V1, y1 }, { r.V2, y2 }, { Vm, ym } }) do
    addplot(options, {
      { a / e, P[2] + r.pente * (a - P[1]) },
      { b / e, P[2] + r.pente * (b - P[1]) },
    })
  end
  addplot("only marks," .. options, { { r.V1 / e, y1 }, { r.V2 / e, y2 }, { r.V / e, y } })
end

--- \addplotfractions : une courbe par forme, pH en abscisse.
function M.tex_fractions(pKa, options, n)
  n = n or 200
  local a, b = M.reglages.pHmin, M.reglages.pHmax
  local courbes = {}
  for i = 1, n do
    local p = a + (b - a) * (i - 1) / (n - 1)
    for j, x in ipairs(M.fractions(pKa, p)) do
      courbes[j] = courbes[j] or {}
      courbes[j][i] = { p, x }
    end
  end
  for _, liste in ipairs(courbes) do addplot(options, liste, true) end
end

local function arrondi(x)
  return ("%." .. M.reglages.chiffres .. "f"):format(x)
end

--- Les nombres, pour \titragepH, \titragevolume, \titrageequivalence, \titragesaut.
--- Le volume versé qui amène au potentiel E (V/ENH), par dichotomie en log V.
function M.volume_redox(t, E)
  local V = volume_brut(t, E)
  if math.abs(M.E_enh(t, V) - E) > 1e-6 then return nil end
  return V
end

--- Le volume versé qui amène la grandeur mesurée à y, par dichotomie sur [0, Vmax] :
--- pM et pP sont monotones le long d'un titrage ; nil hors de la courbe.
local function volume_mesure(t, y, espece)
  local function f(V) return mesure(t, V, espece) - y end
  local a, b = f(0), f(t.Vmax)
  if a ~= a or b ~= b or (a > 0) == (b > 0) then return nil end
  return changement(f, 0, t.Vmax)
end

local function nombre_cp(t, e, quoi, arg, champ)
  if quoi == "pM" then
    return mesure(t, arg * e, champ ~= "" and champ or nil)
  elseif quoi == "volume" then
    return (volume_mesure(t, arg) or 0 / 0) / e
  elseif quoi == "equivalence" then
    return (EQUIVALENCES[t.type](t)[math.floor(arg)] or 0 / 0) / e
  elseif quoi == "saut" then
    local liste = (t.type == "complexation") and M.sauts_complexation(t, t.Vmax, nil, t.suivi)
      or M.sauts_precipitation(t, t.Vmax)
    local s = liste[math.floor(arg)]
    if not s then error(("titrages : pas de saut n° %d"):format(arg)) end
    if champ == "pente" and t.type == "precipitation" and t.suivi then
      return pente_partenaire(t, s[1], s[2], t.suivi) * e
    end
    return ({ volume = s[1] / e, pM = mesure(t, s[1]), pente = s[3] * e })[champ]
  elseif quoi == "evenement" then
    if t.type ~= "precipitation" then error("titrages : les événements sont ceux d'un titrage par précipitation") end
    local v = M.evenements(t, t.Vmax)[math.floor(arg)]
    if not v then error(("titrages : pas d'événement n° %d"):format(arg)) end
    return ({ volume = v[1] / e, pM = mesure(t, v[1]) })[champ]
  end
  error("titrages : « " .. quoi .. " » n'a pas de sens pour ce titrage ; lire la courbe par \\titragepM")
end

function M.tex_nombre(nom, quoi, arg, champ)
  local t = trouver(nom)
  local e = echelle(t)
  local x
  if t.type then
    x = nombre_cp(t, e, quoi, arg, champ)
  elseif quoi == "pM" or quoi == "evenement" then
    error("titrages : « " .. nom .. " » se lit par \\titragepH ou \\titrageE")
  elseif t.redox then
    if quoi == "E" or quoi == "pH" then
      if quoi == "pH" then error("titrages : un titrage rédox se lit par \\titrageE") end
      x = M.E_enh(t, arg * e) - t.ref
    elseif quoi == "volume" then
      x = (M.volume_redox(t, arg + t.ref) or 0 / 0) / e
    elseif quoi == "equivalence" then
      x = (M.equivalences_redox(t)[math.floor(arg)] or 0 / 0) / e
    elseif quoi == "saut" then
      local s = M.sauts_redox(t, t.Vmax)[math.floor(arg)]
      if not s then error(("titrages : pas de saut n° %d"):format(arg)) end
      x = ({ volume = s[1] / e, E = s[2] - t.ref, pente = s[3] * e })[champ]
    elseif quoi == "tangentes" then
      local r = M.tangentes_redox(t, arg * e, t.Vmax)
      x = ({ volume = r.V / e, E = r.E - t.ref, pente = r.pente * e })[champ]
    else
      error("titrages : « " .. quoi .. " » n'a pas de sens pour un titrage rédox")
    end
  elseif quoi == "pH" then
    x = M.pH(t, arg * e)
  elseif quoi == "E" then
    error("titrages : un titrage acide-base se lit par \\titragepH")
  elseif quoi == "volume" then
    x = (M.volume(t, arg) or 0 / 0) / e
  elseif quoi == "equivalence" then
    x = (M.equivalences(t)[math.floor(arg)] or 0 / 0) / e
  elseif quoi == "saut" then
    local s = M.sauts(t, t.Vmax)[math.floor(arg)]
    if not s then error(("titrages : pas de saut n° %d"):format(arg)) end
    x = ({ volume = s[1] / e, pH = s[2], pente = s[3] * e })[champ]
  elseif quoi == "tangentes" then
    local r = M.tangentes(t, arg * e, t.Vmax)
    x = ({ volume = r.V / e, pH = r.pH, pente = r.pente * e })[champ]
  elseif quoi == "conductivite" then
    local V = arg * e
    x = M.conductivite(t, V) * (FACTEURS[t.unite] or 1)
    if t.corrigee then x = x * (t.V0 + V) / t.V0 end
  end
  if x == nil then error("titrages : champ inconnu « " .. tostring(champ) .. " »") end
  sortir(arrondi(x))
end

return M
