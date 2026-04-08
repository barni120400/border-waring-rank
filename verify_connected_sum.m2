p = 61;
load "GenerateGorensteinAlgebra.m2"

logFile = "/tmp/intrinsic_cs_correct_results.txt" << "";
logFile << "CORRECT intrinsic connected sum detection (1+(n-1) split)" << endl;
logFile << "=========================================================" << endl;
logFile << "Criterion: exists v in m/m^2 with rank(L_v) <= 1 AND v^2 != 0" << endl;
logFile << "Variables: n+2 (alpha_i, dI1, dI2). Exact Groebner basis." << endl << endl << flush;

isConnectedSum = (A) -> (
    varList := flatten entries vars A;
    basisA := flatten entries basis A;
    ambA := ambient A;
    n := #varList;
    if n <= 1 then return false;
    -- Polynomial ring: alpha_0..alpha_{n-1}, dI1 (for alpha_i != 0), dI2 (for v^2 != 0)
    varNames := apply(n, i -> (symbol al)_i) | {(symbol dI1), (symbol dI2)};
    S := (ZZ/p)[varNames];
    aVars := apply(n, i -> S_i);
    dI1v := S_n;
    dI2v := S_(n+1);
    -- Build multiplication matrix M(alpha): M_{k,j} = sum_i alpha_i * coeff(basisA_k, e_i*e_j)
    M := matrix apply(#basisA, k -> apply(n, j -> (
        c := 0_S;
        for i from 0 to n-1 do (
            prod := varList#i * varList#j;
            coeff := coefficient(lift(basisA#k, ambA), lift(prod, ambA));
            c = c + aVars#i * sub(coeff, ZZ/p);
        );
        c
    )));
    -- Rank <= 1: all 2x2 minors = 0
    minorIdeal := minors(2, M);
    minorEqns := select(flatten entries gens minorIdeal, eq -> eq != 0_S);
    -- v^2 = sum_{i,j} alpha_i * alpha_j * (e_i * e_j) in A
    -- Extract coefficient of each basis element -> Q_k(alpha)
    vSquaredCoeffs := apply(#basisA, k -> (
        c := 0_S;
        for i from 0 to n-1 do for j from 0 to n-1 do (
            prod := varList#i * varList#j;
            coeff := coefficient(lift(basisA#k, ambA), lift(prod, ambA));
            c = c + aVars#i * aVars#j * sub(coeff, ZZ/p);
        );
        c
    ));
    -- Remove zero Q_k
    nonzeroQs := select(vSquaredCoeffs, q -> q != 0_S);
    if #nonzeroQs == 0 then return false; -- v^2 = 0 for all v -> no CS
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

-- ALL algebras up to dim 8 (excluding alpha-parameterized Ideal9)
allExamples = {
    (3,1,1), (4,1,1), (4,2,2),
    (5,1,1), (5,2,2), (5,3,2),
    (6,1,1), (6,2,2), (6,2,3), (6,2,4), (6,3,2), (6,4,2),
    (7,1,1), (7,2,2), (7,2,3), (7,2,5), (7,3,2), (7,3,3), (7,3,4), (7,4,2), (7,5,2),
    (8,1,1), (8,2,2), (8,2,3), (8,2,5), (8,2,6), (8,2,7), (8,2,8),
    (8,3,2), (8,3,3), (8,3,5), (8,3,10), (8,3,11), (8,3,12), (8,3,13), (8,3,14),
    (8,4,2), (8,4,3), (8,4,4), (8,5,2), (8,6,2)
};

for dim from 3 to 8 do (
    group := select(allExamples, t -> t#0 == dim);
    logFile << "=== Dimension " << dim << " ===" << endl << flush;
    for triple in group do (
        (e, n, t) = toSequence triple;
        logFile << "  A_{" << e << "," << n << "," << t << "} (n=" << n << "): " << flush;
        pair := try generateGorensteinAlgebra(e, n, t, deformed=>false) else (logFile << "FAILED" << endl << flush; continue);
        result := try isConnectedSum(pair#0) else "ERROR";
        logFile << toString result << endl << flush;
    );
    logFile << endl << flush;
);

logFile << "DONE." << endl;
logFile << close;
