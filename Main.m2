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
    if #s >= 2 and s#0 == "$" and s#(#s-1) == "$" then substring(1, #s-2, s) else s
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
        -- Build the \binom{d}{...} part
        j := d - a;  -- sum of non-x_0 exponents
        binomStr := "";
        if j == 0 then (
            -- Pure x_0^d term — coefficient should be 1
            binomStr = "";
        ) else (
            nonZeroExpValues := apply(nonZeroExps, pair -> pair#1);
            if #nonZeroExpValues == 1 then (
                e1 := nonZeroExpValues#0;
                binomStr = "\\binom{d}{" | toString(e1) | "}";
            ) else (
                -- Multiple non-x_0 variables: \binom{d}{e1, e2, ...}
                expStrs := apply(nonZeroExpValues, e -> toString e);
                binomStr = "\\binom{d}{" | concatenate between(", ", expStrs) | "}";
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
        -- Add to result with + separator
        if result == "" then
            result = termStr
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
-- Main computation
-- ============================================================

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

outFile << "\\item Parameters from Casnati's paper:" << endl;
outFile << "\\[" << endl;
outFile << "e = " << idealDegree << ", \\quad n = " << numVars << ", \\quad \\text{idealNumber} = " << idealNumber << endl;
outFile << "\\]" << endl;
outFile << endl;

outFile << "Ideal: " << endl;
outFile << "\\[" << endl;
-- Format ideal generators
idealGens = flatten entries gens idealA;
idealStr = concatenate between(",\\; ", apply(idealGens, g -> cleanTex g));
outFile << "(" << idealStr << ")" << endl;
outFile << "\\]" << endl;
outFile << endl;

outFile << "Generators:" << endl;
outFile << "\\[" << endl;
basisStr = concatenate between(",\\; ", apply(sortedBasisA, b -> cleanTex b));
outFile << basisStr << endl;
outFile << "\\]" << endl;
outFile << endl;

outFile << "Form:" << endl;
outFile << "\\[" << endl;
outFile << genericDForm(form, formDegree) << endl;
outFile << "\\]" << endl;
outFile << endl;

if variableWeights != {} then (
    outFile << "Variable weights:" << endl;
    outFile << "\\[" << endl;
    weightStrs = apply(dimA, i -> "q_" | toString(i) | " = " | toString(
        if i == 0 then 0
        else (
            exps := flatten exponents sortedBasisA#i;
            sum apply(length variableWeights, j -> exps#j * variableWeights#j)
        )
    ));
    outFile << concatenate between(",\\; ", weightStrs) << endl;
    outFile << "\\]" << endl;
    outFile << endl;
) else (
    outFile << "This algebra does not admit a grading." << endl;
    outFile << endl;
);

outFile << "Apolar ideal:" << endl;
outFile << "\\[" << endl;
outFile << "(" << endl;
outFile << genericDApolarIdeal(apolarIdeal, formDegree) << endl;
outFile << ")" << endl;
outFile << "\\]" << endl;

outFile << close;
print("Wrote " | autoFileName);

-- ============================================================
-- Output: mult_*.tex (multiplication operators)
-- ============================================================

multFileName = "output/mult_e" | toString(idealDegree) | "_n" | toString(numVars) | "_i" | toString(idealNumber) | ".tex";
multFile = multFileName << "";

multFile << "% Multiplication operators for e=" << idealDegree << ", n=" << numVars << ", idealNumber=" << idealNumber << endl;
multFile << "Multiplication operators:" << endl;
multFile << "\\[" << endl;
multFile << "\\begin{aligned}" << endl;

for i from 0 to #MTs - 1 do (
    mtPair = MTs#i;
    mtNum = mtPair#0;
    mtDen = mtPair#1;
    -- Compute the actual matrix M^T = mtNum / mtDen
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
