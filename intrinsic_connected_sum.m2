-- Intrinsic connected sum detection for 2-generator algebras
--
-- For A = K[y1,y2]/I, check if there exist linearly independent
-- a = alpha*y1 + beta*y2, b = gamma*y1 + delta*y2 with a*b = 0 in A.
--
-- Method: expand a*b = alpha*gamma*y1^2 + (alpha*delta+beta*gamma)*y1*y2 + beta*delta*y2^2
-- Express in the basis of A and set each coordinate to zero.
-- Add Rabinowitsch constraint alpha*delta - beta*gamma != 0 for linear independence.
-- Check if the polynomial system has solutions (1 not in ideal).
--
-- This is exact: solutions exist over the algebraic closure iff 1 not in ideal.
-- Does NOT depend on the choice of generators (intrinsic property).

p = 61;

isConnectedSum2 = (A) -> (
    basisA := flatten entries basis A;
    ambA := ambient A;
    y1 := (vars A)_(0,0);
    y2 := (vars A)_(0,1);
    prods := {y1^2, y1*y2, y2^2};
    S := (ZZ/p)[(symbol al), (symbol be), (symbol ga), (symbol de), (symbol dI)];
    alv := S_0; bev := S_1; gav := S_2; dev := S_3; dIv := S_4;
    eqns := {};
    for m in basisA do (
        c1 := coefficient(lift(m, ambA), lift(y1^2, ambA));
        c2 := coefficient(lift(m, ambA), lift(y1*y2, ambA));
        c3 := coefficient(lift(m, ambA), lift(y2^2, ambA));
        eq := sub(c1, ZZ/p) * alv*gav + sub(c2, ZZ/p) * (alv*dev + bev*gav) + sub(c3, ZZ/p) * bev*dev;
        if eq != 0_S then eqns = append(eqns, eq);
    );
    eqns = append(eqns, (alv*dev - bev*gav)*dIv - 1_S);
    I := ideal eqns;
    return (1_S % I) != 0;
);

-- Run on all n=2 algebras up to dim 8
load "GenerateGorensteinAlgebra.m2"

logFile = "intrinsic_connected_sum_output.txt" << "";
logFile << "Intrinsic connected sum detection for all n=2 algebras" << endl;
logFile << "======================================================" << endl << endl;

for params in {(6,2,2),(6,2,3),(6,2,4),(6,2,5),(7,2,2),(7,2,3),(7,2,5),(8,2,2),(8,2,3),(8,2,5),(8,2,6),(8,2,7),(8,2,8)} do (
    (e, n, t) = toSequence params;
    logFile << "A_{" << e << "," << n << "," << t << "}: " << flush;
    pair := generateGorensteinAlgebra(e, n, t, deformed=>false);
    result := isConnectedSum2(pair#0);
    logFile << "connected sum = " << toString result << endl << flush;
);

logFile << endl << "DONE." << endl;
logFile << close;
