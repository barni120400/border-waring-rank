load "HelperFunctions.m2"
load "GenerateGorensteinAlgebra.m2"
load "ExtractBorderForms.m2"

-- ============================================================
-- Parameters (edit these)
-- ============================================================
idealDegree = 5;
numVars = 1;
idealNumber = 1;

-- ============================================================
-- Helper functions
-- ============================================================

-- Compute the ungraded socle degree of an Artinian algebra A.
-- This is the largest k such that m^k != 0, where m is the maximal ideal.
computeSocleDegree = A -> (
    m := ideal vars A;
    k := 1;
    while ideal(0_A) != m^k do k = k + 1;
    return k - 1;
);

-- Strip $...$ from M2's tex output
cleanTex = x -> (
    s := tex x;
    s = if #s >= 2 and s#0 == "$" and s#(#s-1) == "$" then substring(1, #s-2, s) else s;
    -- Replace M2's \mathit{alpha} with \alpha
    s = replace("\\\\mathit\\{alpha\\}", "\\\\alpha", s);
    s
);

-- ============================================================
-- Generic-d form formatter
-- ============================================================
-- Given a homogeneous form F of degree d in x_0,...,x_k,
-- produce LaTeX using \binom{d}{...} notation.
-- Each monomial c * x_0^a * x_{i1}^{e1} * ... has c = multinomial(a, e1, ...).
-- We write: \binom{d}{e1, ...} x_0^{d - (e1+...)} x_{i1}^{e1} ...

genericDForm = (F, d) -> (
    R := ring F;
    n := numgens R;
    termsList := terms F;
    -- Sort terms: first by x_0 power (descending = most x_0 first),
    -- then by variable index
    -- Sort terms: by x_0 power descending, then by highest non-x_0 variable descending
    termsList = sort(termsList, t -> (
        exps := (exponents t)_0;
        a := exps#0;
        -- Find the highest non-x_0 variable index with nonzero exponent
        maxVar := 0;
        for i from 1 to #exps - 1 do if exps#i > 0 then maxVar = i;
        (-a, -maxVar)
    ));
    result := "";
    for t in termsList do (
        exps := (exponents t)_0;
        coeff := leadCoefficient t;
        a := exps#0;  -- x_0 exponent
        -- Collect non-x_0 exponents with their variable indices
        nonZeroExps := {};
        for i from 1 to n-1 do (
            if exps#i > 0 then nonZeroExps = append(nonZeroExps, (i, exps#i));
        );
        -- Build the \binom{d}{...} part, accounting for non-multinomial coefficients
        j := d - a;  -- sum of non-x_0 exponents
        binomStr := "";
        if j == 0 then (
            -- Pure x_0^d term
            coeffStr := cleanTex coeff;
            if coeff != 1 then binomStr = coeffStr | " ";
        ) else (
            -- Compute the expected multinomial coefficient d!/(a! * e1! * ... * ek!)
            allExps := {a} | apply(nonZeroExps, pair -> pair#1);
            expectedMultinom := (d)! / product apply(allExps, e -> e!);
            -- Ratio of actual coefficient to expected multinomial
            ratio := coeff / expectedMultinom;
            nonZeroExpValues := apply(nonZeroExps, pair -> pair#1);
            binomPart := if #nonZeroExpValues == 1 then (
                "\\binom{d}{" | toString(nonZeroExpValues#0) | "}"
            ) else (
                expStrs := apply(nonZeroExpValues, e -> toString e);
                "\\binom{d}{" | concatenate between(", ", expStrs) | "}"
            );
            if ratio == 1 then (
                binomStr = binomPart;
            ) else if ratio == -1 then (
                binomStr = "-" | binomPart;
            ) else (
                ratioStr := cleanTex ratio;
                binomStr = ratioStr | " " | binomPart;
            );
        );
        -- Build the x_0 power part (always use generic d notation)
        x0Str := "";
        if j == 0 then (
            x0Str = "x_0^d";
        ) else if j == 1 then (
            x0Str = "x_0^{d-1}";
        ) else (
            x0Str = "x_0^{d - " | toString(j) | "}";
        );
        -- Build the non-x_0 variable parts with spaces
        varStr := "";
        for pair in nonZeroExps do (
            idx := pair#0;
            e := pair#1;
            if varStr != "" then varStr = varStr | " ";
            if e == 1 then
                varStr = varStr | "x_" | toString(idx)
            else
                varStr = varStr | "x_" | toString(idx) | "^{" | toString(e) | "}";
        );
        -- Combine with spaces between parts
        termStr := "";
        if binomStr != "" then termStr = termStr | binomStr | " ";
        if x0Str != "" then termStr = termStr | x0Str;
        if varStr != "" then (
            if x0Str != "" then termStr = termStr | " ";
            termStr = termStr | varStr;
        );
        -- Add to result with +/- separator
        if result == "" then
            result = termStr
        else if #termStr > 0 and termStr#0 == "-" then
            result = result | "\n- " | substring(1, termStr)
        else
            result = result | "\n+ " | termStr;
    );
    return result;
);

-- ============================================================
-- Generic-d apolar ideal formatter
-- ============================================================
-- For each generator of the apolar ideal, replace x_0^d with x_0^{d}
-- and x_0^{d-1} with x_0^{d - 1}.

genericDApolarGen = (g, d) -> (
    R := ring g;
    s := cleanTex g;
    -- Check if this generator involves x_0^d or x_0^{d-1}
    exps := (exponents leadMonomial g)_0;
    a := exps#0;  -- x_0 power of the lead monomial
    -- For generators with large x_0 powers, do symbolic replacement
    termList := terms g;
    if #termList == 1 and a == d then (
        -- Pure x_0^d
        return "x_0^{d}";
    );
    if #termList == 1 and a == d - 1 then (
        -- x_0^{d-1} * x_i
        result := "x_0^{d - 1}";
        for i from 1 to numgens R - 1 do (
            if exps#i > 0 then (
                if exps#i == 1 then
                    result = result | " x_" | toString(i)
                else
                    result = result | " x_" | toString(i) | "^" | toString(exps#i);
            );
        );
        return result;
    );
    -- For generators involving a mix (e.g., x_1*x_3 - x_0*x_4), use M2's tex
    -- but clean up formatting
    return s;
);

-- Format the full apolar ideal, with truncation generators grouped at end.
-- Truncation generators are single-term monomials: x_0^d or x_0^{d-1} * x_i.
genericDApolarIdeal = (I, d) -> (
    genList := flatten entries gens I;
    isTruncationGen := g -> (
        if #(terms g) > 1 then return false;
        exps := (exponents g)_0;
        exps#0 >= d - 1
    );
    lowGens := select(genList, g -> not isTruncationGen(g));
    highGens := select(genList, g -> isTruncationGen(g));
    allGens := lowGens | highGens;
    strs := apply(allGens, g -> genericDApolarGen(g, d));
    return concatenate between(",\\;\n", strs);
);

-- ============================================================
-- Multiplication operator LaTeX formatter
-- ============================================================

texMatrix = M -> (
    nR := numrows M;
    nC := numcols M;
    result := "\\begin{pmatrix}\n";
    for i from 0 to nR - 1 do (
        row := "";
        for j from 0 to nC - 1 do (
            entry := M_(i,j);
            entryStr := if entry == 0 then "0" else cleanTex entry;
            if j > 0 then row = row | "&";
            row = row | entryStr;
        );
        if i < nR - 1 then row = row | "\\\\";
        result = result | row | "\n";
    );
    result = result | "\\end{pmatrix}";
    return result;
);

-- ============================================================
-- Intrinsic connected sum detection
-- ============================================================
-- Checks if A is a connected sum, independent of the choice of generators.
--
-- Criterion: A is a connected sum (1+(n-1) split) iff there exists
-- v in m/m^2 such that:
--   (1) rank(L_v) <= 1, where L_v: m/m^2 -> A is multiplication by v, AND
--   (2) v^2 != 0 in A (so the kernel of L_v is a complement to v).
--
-- This detects 1+(n-1) splits only. For multiplicity <= 9 this is sufficient:
-- the smallest indecomposable algebra with >= 2 generators has dim 6 (A_{6,2,4}),
-- so a k+(n-k) split with k >= 2 where both sides are indecomposable requires
-- dim >= 6+6-2 = 10 > 9. Hence within the Casnati classification (dim <= 9),
-- every connected sum has at least one single-generator component.
--
-- Method: polynomial system in n+2 variables (alpha_i, dI1, dI2).
-- All 2x2 minors of the multiplication matrix M(alpha) = 0 (rank <= 1),
-- plus Rabinowitsch constraints for v != 0 and v^2 != 0.
-- Solved exactly by Groebner basis over ZZ/p.
isConnectedSum = A -> (
    varList := flatten entries vars A;
    basisA := flatten entries basis A;
    ambA := ambient A;
    n := #varList;
    if n <= 1 then return false;
    -- Polynomial ring: alpha_0..alpha_{n-1}, dI1 (v != 0), dI2 (v^2 != 0)
    -- Use the algebra's coefficient ring (may contain alpha parameter)
    baseField := coefficientRing A;
    varNames := apply(n, i -> (symbol al)_i) | {(symbol dI1), (symbol dI2)};
    S := baseField[varNames];
    aVars := apply(n, i -> S_i);
    dI1v := S_n;
    dI2v := S_(n+1);
    -- Build multiplication matrix M(alpha): M_{k,j} = sum_i alpha_i * coeff(basisA_k, e_i*e_j)
    M := matrix apply(#basisA, k -> apply(n, j -> (
        c := 0_S;
        for i from 0 to n-1 do (
            prod := varList#i * varList#j;
            coeff := coefficient(lift(basisA#k, ambA), lift(prod, ambA));
            c = c + aVars#i * sub(coeff, S);
        );
        c
    )));
    -- Rank <= 1: all 2x2 minors = 0
    minorIdeal := minors(2, M);
    minorEqns := select(flatten entries gens minorIdeal, eq -> eq != 0_S);
    -- v^2 = sum_{i,j} alpha_i * alpha_j * (e_i * e_j), extract each basis coordinate
    vSquaredCoeffs := apply(#basisA, k -> (
        c := 0_S;
        for i from 0 to n-1 do for j from 0 to n-1 do (
            prod := varList#i * varList#j;
            coeff := coefficient(lift(basisA#k, ambA), lift(prod, ambA));
            c = c + aVars#i * aVars#j * sub(coeff, S);
        );
        c
    ));
    nonzeroQs := select(vSquaredCoeffs, q -> q != 0_S);
    if #nonzeroQs == 0 then return false;
    -- Try each alpha_i (v != 0) and each Q_k (v^2 != 0)
    for i from 0 to n-1 do (
        for qk in nonzeroQs do (
            testEqns := minorEqns | {aVars#i * dI1v - 1_S, qk * dI2v - 1_S};
            I := ideal testEqns;
            if (1_S % I) != 0 then return true;
        );
    );
    return false;
);

-- Old naive detection (on given generators only) — kept for the decomposition output.
-- Returns list of components, each a list of 0-based variable indices.
detectConnectedSumNaive = A -> (
    varList := flatten entries vars A;
    n := #varList;
    if n <= 1 then return {toList(0..n-1)};
    adj := new MutableList from apply(n, i -> {});
    for i from 0 to n-1 do (
        for j from i+1 to n-1 do (
            if varList#i * varList#j != 0_A then (
                adj#i = append(adj#i, j);
                adj#j = append(adj#j, i);
            );
        );
    );
    visited := new MutableList from apply(n, i -> false);
    csComponents := {};
    for start from 0 to n-1 do (
        if not visited#start then (
            comp := {};
            queue := {start};
            while #queue > 0 do (
                v := queue#0;
                queue = drop(queue, 1);
                if not visited#v then (
                    visited#v = true;
                    comp = append(comp, v);
                    for u in adj#v do (
                        if not visited#u then queue = append(queue, u);
                    );
                );
            );
            csComponents = append(csComponents, comp);
        );
    );
    return csComponents;
);

-- Compute the dimension of the subalgebra of A generated by a set of variables.
-- Works by finding the span of all monomials in these variables within A.
subalgebraDim = (A, componentIndices) -> (
    varList := flatten entries vars A;
    compVars := apply(componentIndices, i -> varList#i);
    -- Start with {1}, iteratively add products with component variables
    basisSoFar := {1_A};
    frontier := {1_A};
    while #frontier > 0 do (
        newElems := {};
        for b in frontier do (
            for v in compVars do (
                prod := b * v;
                if prod != 0_A then (
                    -- Check linear independence with existing basis
                    testBasis := basisSoFar | {prod};
                    if rank matrix apply(testBasis, b ->
                        apply(flatten entries basis A, m -> coefficient(m, lift(b, ambient A)))
                    ) > #basisSoFar then (
                        basisSoFar = append(basisSoFar, prod);
                        newElems = append(newElems, prod);
                    );
                );
            );
        );
        frontier = newElems;
    );
    return #basisSoFar;
);

-- ============================================================
-- Identify the Casnati type of a multi-generator subalgebra component.
-- Given algebra A and component indices, compute the subalgebra's Hilbert
-- function and match against all valid Casnati types for that (dim, ngens).
-- Returns (idealNumber, label, name) or (0, label, name) if unknown.
-- ============================================================

-- Valid Casnati types for each (e, n) pair, from the Ideal*.m2 constraints
validTypes = (d, k) -> (
    result := {};
    -- Ideal 2: n >= 2, e >= n+2
    if k >= 2 and d >= k + 2 then result = append(result, 2);
    -- Ideal 3: n >= 2, e >= n+4
    if k >= 2 and d >= k + 4 then result = append(result, 3);
    -- Ideal 4: n >= 2, e = n+4
    if k >= 2 and d == k + 4 then result = append(result, 4);
    -- Ideal 5: n >= 2, e >= n+5 (e=n+4 is isomorphic to type 3, see rem:type-3-5-isomorphism)
    if k >= 2 and d >= k + 5 then result = append(result, 5);
    -- Ideal 6: n >= 2, e = n+6
    if k >= 2 and d == k + 6 then result = append(result, 6);
    -- Ideal 7: n >= 2, e = n+6
    if k >= 2 and d == k + 6 then result = append(result, 7);
    -- Ideal 8: n >= 2, e = n+6
    if k >= 2 and d == k + 6 then result = append(result, 8);
    -- Ideal 9: n >= 3, e = n+5 (needs alpha)
    if k >= 3 and d == k + 5 then result = append(result, 9);
    -- Ideals 10-14: n >= 3, e = n+5
    if k >= 3 and d == k + 5 then result = result | {10, 11, 12, 13, 14};
    -- Ideals 15-19: n = 2, e = 9
    if k == 2 and d == 9 then result = result | {15, 16, 17, 18, 19};
    -- Ideals 20-26: n = 3, e = 9
    if k == 3 and d == 9 then result = result | {20, 21, 22, 23, 24, 25, 26};
    return result;
);

-- Compute the Hilbert function of the subalgebra of A generated by component vars.
subalgebraHilbert = (A, componentIndices) -> (
    varList := flatten entries vars A;
    compVars := apply(componentIndices, i -> varList#i);
    socleDeg := computeSocleDegree(A);
    -- For each degree k, count the number of linearly independent monomials
    -- of degree k in the component variables
    ambBasis := flatten entries basis A;
    basisSoFar := {1_A};
    frontier := {1_A};
    while #frontier > 0 do (
        newElems := {};
        for b in frontier do (
            for v in compVars do (
                prod := b * v;
                if prod != 0_A then (
                    testBasis := basisSoFar | {prod};
                    if rank matrix apply(testBasis, b ->
                        apply(ambBasis, m -> coefficient(m, lift(b, ambient A)))
                    ) > #basisSoFar then (
                        basisSoFar = append(basisSoFar, prod);
                        newElems = append(newElems, prod);
                    );
                );
            );
        );
        frontier = newElems;
    );
    -- Compute degree of each basis element in terms of the component variables
    -- (degree = min number of component variable multiplications to reach it)
    -- Actually, just use the standard grading of the ambient ring
    hilbert := apply(0..socleDeg, deg -> (
        #select(basisSoFar, b -> (degree b)_0 == deg)
    ));
    -- Trim trailing zeros
    while #hilbert > 0 and last hilbert == 0 do hilbert = drop(hilbert, -1);
    return hilbert;
);

-- Check if k[z_1..z_k]/I is isomorphic to k[z_1..z_k]/J via a linear
-- change of variables. Uses a FLAT ring (all variables in one ring) so
-- that ideal membership and Groebner basis computations are correct.
-- Iff criterion: solutions exist over the algebraic closure iff 1 ∉ ideal.
areIsomorphicIdeals = (I, J) -> (
    R := ring I;
    k := numgens R;
    baseField := coefficientRing R;
    -- Create flat ring: baseField[dd, mm_{i,j}, zz_0..zz_{k-1}]
    -- dd is used for the Rabinowitsch trick to enforce det(M) != 0
    matVarNames := flatten apply(k, i -> apply(k, j -> (symbol mm)_(i,j)));
    zVarNames := apply(k, i -> (symbol zz)_i);
    flatR := baseField[{(symbol detInv)} | matVarNames | zVarNames, MonomialOrder => Eliminate (1 + k*k)];
    -- dd is at index 0, mm at indices 1..k^2, zz at indices k^2+1..k^2+k
    detInvVar := flatR_0;
    zVarsFlat := apply(k, i -> flatR_(1 + k*k + i));
    mmVarsFlat := apply(k, i -> apply(k, j -> flatR_(1 + i*k + j)));
    -- Build substitution: z_i -> sum_j mm_{i,j} * zz_j
    transformedZ := apply(k, i -> sum apply(k, j -> mmVarsFlat#i#j * zVarsFlat#j));
    -- Map generators of I to flatR, substitute, get transformed generators
    mapRtoFlat := map(flatR, R, zVarsFlat);
    gensI := flatten entries gens I;
    gensJ := flatten entries gens J;
    -- Transform each generator of I by the linear map
    transformedGensI := apply(gensI, g -> (
        gFlat := mapRtoFlat g;
        sub(gFlat, apply(k, i -> zVarsFlat#i => transformedZ#i))
    ));
    -- J in flat ring
    JFlat := ideal apply(gensJ, g -> mapRtoFlat g);
    -- For each transformed generator, reduce modulo J and collect equations
    -- The remainder must be zero: extract z-monomial coefficients and set to 0
    eqnList := flatten apply(transformedGensI, gT -> (
        remPoly := gT % JFlat;
        (monMat, coeffMat) := coefficients(remPoly, Variables => zVarsFlat);
        flatten entries coeffMat
    ));
    -- Add Rabinowitsch constraint: det(M) * dd - 1 = 0 to enforce M invertible
    M := matrix apply(k, i -> apply(k, j -> mmVarsFlat#i#j));
    detConstraint := det(M) * detInvVar - 1_(flatR);
    eqnList = eqnList | {detConstraint};
    if #eqnList == 0 then return true;
    eqnIdeal := ideal eqnList;
    -- 1 not in ideal iff solutions exist over algebraic closure
    result := (1_(flatR) % eqnIdeal) != 0;
    return result;
);

-- Identify the type of a subalgebra component.
-- Returns (dim, ngens, typeNumber) where typeNumber=0 means unknown.
identifyComponentType = (A, componentIndices) -> (
    d := subalgebraDim(A, componentIndices);
    k := #componentIndices;
    if k == 1 then return (d, 1, 1);
    -- Compute the subalgebra ideal via kernel of inclusion map
    varList := flatten entries vars A;
    compVars := apply(componentIndices, i -> varList#i);
    baseField := ZZ/p;
    subR := baseField[z_1..z_k, MonomialOrder=>GRevLex];
    phi := map(A, subR, compVars);
    subIdeal := trim kernel phi;
    -- Try all valid candidate types with matching Hilbert function first (cheap filter)
    subAlg := subR / subIdeal;
    subSocleDeg := computeSocleDegree(subAlg);
    subHilbert := apply(0..subSocleDeg, deg -> hilbertFunction(deg, subAlg));
    candidates := validTypes(d, k);
    -- First pass: check for EXACT ideal equality (preferred — gives natural type label)
    for t in candidates do (
        candidatePair := try generateGorensteinAlgebra(d, k, t, deformed=>false) else continue;
        candAlg := candidatePair#0;
        ambCand := ambient candAlg;
        mapToSubR := map(subR, ambCand, gens subR);
        candGens := flatten entries gens ideal candAlg;
        candIdealInSubR := ideal apply(candGens, g -> mapToSubR lift(g, ambCand));
        if subIdeal == candIdealInSubR then return (d, k, t);
    );
    -- Second pass: check GL_k orbit isomorphism (for when ideals differ but algebras are isomorphic)
    for t in candidates do (
        candidatePair := try generateGorensteinAlgebra(d, k, t, deformed=>false) else continue;
        candAlg := candidatePair#0;
        candSocleDeg := computeSocleDegree(candAlg);
        candHilbert := apply(0..candSocleDeg, deg -> hilbertFunction(deg, candAlg));
        if subHilbert != candHilbert then continue;
        ambCand := ambient candAlg;
        mapToSubR := map(subR, ambCand, gens subR);
        candGens := flatten entries gens ideal candAlg;
        candIdealInSubR := ideal apply(candGens, g -> mapToSubR lift(g, ambCand));
        if try areIsomorphicIdeals(subIdeal, candIdealInSubR) else false then (
            return (d, k, t);
        );
    );
    return (d, k, 0);
);

-- ============================================================
-- Main computation (skipped when TESTING flag is set)
-- ============================================================
if not (try TESTING else false) then (

pair = generateGorensteinAlgebra(idealDegree, numVars, idealNumber, deformed=>false);
algebra = pair#0;
variableWeights = pair#1;
formDegree = computeSocleDegree(algebra);
dimA = degree algebra;

print("Algebra dimension: " | toString(dimA));
print("Socle degree: " | toString(formDegree));

form = extractBorderForms(algebra, variableWeights, formDegree);
apolarIdeal = perpIdeal(form);

-- Get basis sorted by weight (matching the order used by extractBorderForms)
idealA = ideal algebra;
rawBasisA = flatten entries basis algebra;
if variableWeights != {} then (
    -- Use sortByWeight to replicate the exact ordering used in extractBorderForms
    -- sortByWeight returns monomials in a k[y_1..y_k] ring, so we substitute back
    baseCoeffRing := coefficientRing algebra;
    tempRing := baseCoeffRing[y_1..y_(length variableWeights), MonomialOrder=>GRevLex];
    monListTemp := apply(rawBasisA, m -> substitute(m, tempRing));
    sortedTemp := sortByWeight(monListTemp, variableWeights);
    sortedBasisA = apply(sortedTemp, m -> substitute(m, ring rawBasisA#0));
) else (
    sortedBasisA = sort rawBasisA;
);

-- Compute multiplication operators
R := ring form;
dehomForm = sub(form, {x_0 => 1});
dehomRing = ring dehomForm;
(H, ML) = hankelOperator(dehomForm, formDegree);
B = {1_dehomRing} | toList apply(1..dimA-1, i -> x_i_dehomRing);
varsList = toList apply(1..dimA-1, i -> x_i_dehomRing);
MTs = apply(varsList, v -> multiplicationOperatorTranspose(v, H, B, ML));

-- ============================================================
-- Output: example.tex (paper-facing, no multiplication operators)
-- ============================================================

autoFileName = "output/auto_e" | toString(idealDegree) | "_n" | toString(numVars) | "_i" | toString(idealNumber) | ".tex";
outFile = autoFileName << "";

-- Label and algebra name
algebraName = "[A_{" | toString(idealDegree) | "," | toString(numVars) | "," | toString(idealNumber) | "}]";
labelStr = "alg:e" | toString(idealDegree) | "_n" | toString(numVars) | "_i" | toString(idealNumber);

-- Format ideal generators inline
idealGens = flatten entries gens idealA;
idealStr = concatenate between(",\\; ", apply(idealGens, g -> cleanTex g));

-- Format variable list for the ring
varListStr = concatenate between(", ", apply(flatten entries vars(ambient algebra), v -> cleanTex v));

-- Algebra name as header
outFile << "\\item \\label{" << labelStr << "} $" << algebraName << "$." << endl;
outFile << endl;

-- Description list for structured fields
outFile << "\\begin{description}[leftmargin=1.5cm, style=sameline]" << endl;
outFile << "\\item[Ideal] $(" << idealStr << ") \\subset \\bK[" << varListStr << "]$" << endl;

-- Basis
basisStr = concatenate between(",\\; ", apply(sortedBasisA, b -> cleanTex b));
outFile << "\\item[Basis] $" << basisStr << "$" << endl;

-- Weights (or note about no grading)
if variableWeights != {} then (
    weightVals = apply(dimA, i -> toString(
        if i == 0 then 0
        else (
            exps := flatten exponents sortedBasisA#i;
            sum apply(length variableWeights, j -> exps#j * variableWeights#j)
        )
    ));
    outFile << "\\item[Weights] $\\vw = (" << concatenate between(", ", weightVals) << ")$" << endl;
) else (
    outFile << "\\item[Weights] This algebra does not admit a grading." << endl;
);

-- Connected sum detection (intrinsic, basis-independent)
isCS = isConnectedSum(algebra);
-- Use naive decomposition for displaying the components (works when given generators split)
csComponents = if isCS then detectConnectedSumNaive(algebra) else {toList(0..numgens(ambient algebra)-1)};
if isCS then (
    compDescs := apply(#csComponents, idx -> (
        (d, k, t) := identifyComponentType(algebra, csComponents#idx);
        if t > 0 then (
            summandLabel := "alg:e" | toString(d) | "_n" | toString(k) | "_i" | toString(t);
            summandName := "[A_{" | toString(d) | "," | toString(k) | "," | toString(t) | "}]";
            "\\hyperref[" | summandLabel | "]{" | summandName | "}"
        ) else (
            "[A_{" | toString(d) | "," | toString(k) | ",?}]"
        )
    ));
    outFile << "\\item[Decomposition] $" << algebraName << " = " << concatenate between(" \\# ", compDescs) << "$" << endl;
);

-- Canonical form
-- Canonical form (use display math for line wrapping)
outFile << "\\item[Canonical form]" << endl;
outFile << "\\[" << endl;
outFile << "F = " << genericDForm(form, formDegree) << endl;
outFile << "\\]" << endl;

-- Apolar ideal (use display math for line wrapping)
outFile << "\\item[Apolar ideal]" << endl;
outFile << "\\[" << endl;
outFile << "(" << endl;
outFile << genericDApolarIdeal(apolarIdeal, formDegree) << endl;
outFile << ")" << endl;
outFile << "\\]" << endl;

outFile << "\\end{description}" << endl;

outFile << close;
print("Wrote " | autoFileName);

-- ============================================================
-- Output: mult_*.tex (multiplication operators)
-- ============================================================

multFileName = "output/mult_e" | toString(idealDegree) | "_n" | toString(numVars) | "_i" | toString(idealNumber) | ".tex";
if MTs =!= null then (
    multFile = multFileName << "";
    multFile << "% Multiplication operators for e=" << idealDegree << ", n=" << numVars << ", idealNumber=" << idealNumber << endl;
    multFile << "Multiplication operators:" << endl;
    multFile << "\\[" << endl;
    multFile << "\\begin{aligned}" << endl;
    for i from 0 to #MTs - 1 do (
        mtPair = MTs#i;
        mtNum = mtPair#0;
        mtDen = mtPair#1;
        mtMatrix = (1/mtDen) * mtNum;
        varIdx = i + 1;
        sep = if i < #MTs - 1 and i % 2 == 0 then ","
              else if i < #MTs - 1 and i % 2 == 1 then ""
              else ".";
        lineBreak = if i < #MTs - 1 and i % 2 == 1 then "\n\\\\[1.2em]" else "";
        qquad = if i % 2 == 1 or i == #MTs - 1 then "" else "\n\\qquad";
        multFile << "M_{x_" << varIdx << "}^T &=" << endl;
        multFile << texMatrix(mtMatrix) << sep << qquad << lineBreak << endl;
    );
    multFile << "\\end{aligned}" << endl;
    multFile << "\\]" << endl;
    multFile << close;
    print("Wrote " | multFileName);
) else (
    print("Skipped multiplication operators (symbolic coefficients)");
);
); -- end of TESTING guard
