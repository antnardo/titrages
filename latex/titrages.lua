-- titrages.lua — le moteur de calcul du paquet LaTeX titrages.
--
-- Le même algorithme que le paquet Python titrages (https://pypi.org/project/titrages),
-- en Lua 5.3 pur : il tourne sous LuaTeX comme sous un lua ordinaire, ce qui permet de
-- le tester hors de TeX (tests/test_titrages.lua).
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
-- Copyright (c) 2026 antnardo — licence MIT (voir LICENSE).

local M = {}
M.version = "0.1.0"
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

local function trouver(nom)
  local t = M.titrages[nom]
  if not t then error(("titrages : aucun titrage « %s » ; le définir par \\titragedefinir"):format(nom)) end
  return t
end

--- Le diviseur des volumes : la clé echelle, ou le premier volume équivalent (abscisse n).
local function echelle(t)
  if t.abscisse ~= "n" then return t.echelle or 1 end
  local eq = t.redox and M.equivalences_redox(t) or M.equivalences(t)
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
  else
    for i, c in ipairs(M.courbe(t, t.Vmax, t.points)) do liste[i] = { c[1] / e, c[2] } end
  end
  addplot(options, liste, true)
end

local FACTEURS = { ["S/m"] = 1, ["mS/cm"] = 10, ["mS/m"] = 1000 }

--- \addplottitrageconductivite : sigma(V), à pas de volume constant.
function M.tex_conductivite(nom, options)
  local t = trouver(nom)
  if t.redox then error("titrages : la conductimétrie d'un titrage rédox n'est pas prise en charge") end
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
  else
    for i, s in ipairs(M.sauts(t, t.Vmax)) do liste[i] = { s[1] / e, s[2] } end
  end
  addplot("only marks," .. options, liste)
end

--- \addplottitragetangentes : en pH ou en E (par rapport à la référence).
function M.tex_tangentes(nom, V1, options)
  local t = trouver(nom)
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

function M.tex_nombre(nom, quoi, arg, champ)
  local t = trouver(nom)
  local e = echelle(t)
  local x
  if t.redox then
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
