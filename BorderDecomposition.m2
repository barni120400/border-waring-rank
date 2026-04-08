load "GenerateGorensteinAlgebra.m2"
load "ExtractBorderForms.m2"

-- ============================================
-- Step 1: Compute the target form F
-- ============================================
idealDegree = 5;
numVars = 2;
idealNumber = 2;
d = idealDegree - 2;  -- form degree = 3

pair = generateGorensteinAlgebra(idealDegree, numVars, idealNumber, deformed=>false);
algebra = pair#0;
variableWeights = pair#1;

F = extractBorderForms(algebra, variableWeights, d);
numXVars = degree algebra;  -- = 5

print "Target form F:";
print tex F;
print("Number of x-variables: " | toString(numXVars));
print("Form degree: " | toString(d));

-- Store F's coefficients for later matching
formRing = ring F;
(fMonMat, fCoeffMat) = coefficients(F);
fMons = flatten entries fMonMat;
fCoeffs = flatten entries fCoeffMat;
fTable = new HashTable from apply(#fMons, i -> fMons#i => fCoeffs#i);

-- ============================================
-- Step 2: Build ring tower with 72 unknowns
-- ============================================
r = 12;  -- number of linear forms
kk = ZZ/p;

-- Fix L_1 = x_0 (by GL action): b_1 = 0, c_(1,j) = 0
-- Unknowns: l_i (scalar coefficients), b_i for i>=2, c_(i,j) for i>=2
lVars = toList apply(1..r, i -> l_i);
bVars = toList apply(2..r, i -> b_i);
cVars = flatten toList apply(2..r, i -> toList apply(1..(numXVars-1), j -> c_(i,j)));
allUnknowns = lVars | bVars | cVars;

print("");
print("Number of unknowns: " | toString(#allUnknowns));

-- Tower: coeffR[eps][x_0..x_4]
coeffR = kk[allUnknowns];
epsR = coeffR[eps];
R = epsR[x_0..x_(numXVars-1)];

-- Promote eps into R
epsInR = sub(eps, R);

-- ============================================
-- Step 3: Build linear forms and compute sum of cubes
-- ============================================
print "Computing sum of weighted cubes...";

-- L_1 = x_0 (fixed), so l_1 * L_1^3 = l_1 * x_0^3
totalSum = sub(l_1, R) * x_0^3 + sum toList apply(2..r, i -> (
    Li := (1 + sub(b_i, R) * epsInR) * x_0
          + sum toList apply(1..(numXVars-1), j -> sub(c_(i,j), R) * epsInR * x_j);
    sub(l_i, R) * Li^3
));

print "Done computing.";

-- ============================================
-- Step 4: Extract equations
-- ============================================
print "Extracting equations...";

-- Stage 1: extract coefficients of x-monomials
xVarsList = toList apply(0..(numXVars-1), i -> x_i);
(monMat, coeffMat) = coefficients(totalSum, Variables => xVarsList);
xMonomials = flatten entries monMat;
epsPolys = flatten entries coeffMat;

-- Stage 2: for each x-monomial, extract eps^k coefficients
equations = {};
scan(#xMonomials, idx -> (
    xMon := xMonomials#idx;
    epsPoly := sub(epsPolys#idx, epsR);

    scan(toList(0..d), k -> (
        ck := coefficient(eps^k, epsPoly);

        if k < d then (
            -- eps^0, eps^1, eps^2: must be zero
            if ck != 0 then
                equations = append(equations, (xMon, k, ck));
        ) else (
            -- eps^3: must match F
            xMonInFormRing := sub(xMon, formRing);
            fVal := if fTable#?xMonInFormRing then sub(fTable#xMonInFormRing, coeffR) else 0_coeffR;
            residual := ck - fVal;
            if residual != 0 then
                equations = append(equations, (xMon, k, residual));
        );
    ));
));

print("Number of non-trivial equations: " | toString(#equations));

-- ============================================
-- Step 5: Export to msolve format
-- ============================================
print "Exporting equations to msolve format...";

eqPolys = apply(equations, eq -> eq#2);

-- Build variable name list matching coeffR's generators
varNames = apply(gens coeffR, v -> (
    s := toString v;
    -- msolve doesn't like subscript notation; replace _ with nothing or use flat names
    s
));

outFile = "border_waring.ms" << "";

-- Line 1: comma-separated variable names
outFile << concatenate between(",", varNames) << endl;

-- Line 2: the prime
outFile << toString(char kk) << endl;

-- Line 3+: the polynomials, comma-separated, ending with newline
scan(#eqPolys, i -> (
    outFile << toString eqPolys#i;
    if i < #eqPolys - 1 then outFile << "," << endl
    else outFile << endl;
));

outFile << close;

print("Exported " | toString(#eqPolys) | " equations in " | toString(#varNames) | " unknowns to border_waring.ms");
print("Prime: " | toString(char kk));

-- ============================================
-- Step 6: LaTeX output
-- ============================================
cleanTex = x -> (
    s := tex x;
    if #s >= 2 and s#0 == "$" and s#(#s-1) == "$" then substring(1, #s-2, s) else s
);

<< endl;
<< "% =============================================" << endl;
<< "% Border Waring Decomposition Equations" << endl;
<< "% =============================================" << endl;
<< endl;
<< "% Form: F = " << cleanTex F << endl;
<< "% Decomposition: sum_{i=1}^{12} l_i * L_i(eps)^3" << endl;
<< "% L_i(eps) = (1 + b_i eps) x_0 + c_{i,1} eps x_1 + ... + c_{i,4} eps x_4" << endl;
<< "% Unknowns: l_i, b_i, c_{i,j} (" << #allUnknowns << " total)" << endl;
<< endl;

scan(toList(0..d), k -> (
    eqs := select(equations, eq -> eq#1 == k);
    if #eqs > 0 then (
        << "% --- epsilon^" << k << " equations (" << #eqs << " equations) ---" << endl;
        << "\\subsection*{$\\varepsilon^{" << k << "}$ equations}" << endl;
        << "\\begin{align*}" << endl;
        scan(#eqs, idx -> (
            (myMon, myEpsPow, myPoly) := eqs#idx;
            << "  \\text{Coeff of } " << cleanTex myMon << ": \\quad & ";
            << cleanTex myPoly;
            if idx < #eqs - 1 then << " = 0 \\\\" << endl
            else << " = 0" << endl;
        ));
        << "\\end{align*}" << endl;
        << endl;
    );
));
