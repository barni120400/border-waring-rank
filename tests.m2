-- Exhaustive test suite for border-waring-rank computations
-- Run: M2 --script tests.m2

-- We need to source Main.m2's function definitions without running the computation.
-- Set a flag so Main.m2 can skip its computation section.
TESTING = true;

load "Main.m2"

-- Test infrastructure
passed = 0;
failed = 0;
assert' = (name, condition) -> (
    if condition then (
        stdio << "  PASS: " << name << endl;
        passed = passed + 1;
    ) else (
        stdio << "  *** FAIL: " << name << " ***" << endl;
        failed = failed + 1;
    );
);

-- ============================================================
stdio << "=== computeSocleDegree ===" << endl;
-- ============================================================

pair311 = generateGorensteinAlgebra(3, 1, 1, deformed=>false);
assert'("A_{3,1,1} socle degree = 2", computeSocleDegree(pair311#0) == 2);

pair511 = generateGorensteinAlgebra(5, 1, 1, deformed=>false);
assert'("A_{5,1,1} socle degree = 4", computeSocleDegree(pair511#0) == 4);

pair522 = generateGorensteinAlgebra(5, 2, 2, deformed=>false);
assert'("A_{5,2,2} socle degree = 3", computeSocleDegree(pair522#0) == 3);

pair532 = generateGorensteinAlgebra(5, 3, 2, deformed=>false);
assert'("A_{5,3,2} socle degree = 2", computeSocleDegree(pair532#0) == 2);

pair828 = generateGorensteinAlgebra(8, 2, 8, deformed=>false);
sd828 = computeSocleDegree(pair828#0);
assert'("A_{8,2,8} socle degree is positive integer", sd828 > 0);

-- ============================================================
stdio << endl << "=== cleanTex ===" << endl;
-- ============================================================

R = (ZZ/p)[x_1];
s1 = cleanTex(x_1_R);
assert'("cleanTex strips $", #s1 > 0 and s1#0 != "$");

-- Alpha replacement (needs alpha in scope from Ideal9 loading)
K = frac(ZZ/p[alpha]);
RK = K[x_1];
s2 = cleanTex(alpha * x_1_RK);
assert'("cleanTex replaces mathit{alpha} with \\alpha", match("\\\\alpha", s2));

-- ============================================================
stdio << endl << "=== genericDForm ===" << endl;
-- ============================================================

-- Standard form: all multinomial coefficients
R1 = (ZZ/p)[x_0, x_1, x_2, x_3, x_4];
F1 = 4*x_0^3*x_4 + 12*x_0^2*x_1*x_3 + 6*x_0^2*x_2^2 + 12*x_0*x_1^2*x_2 + x_1^4;
gd1 = genericDForm(F1, 4);
assert'("Standard form has \\binom{d}{1}", match("\\\\binom\\{d\\}\\{1\\}", gd1));
assert'("Standard form has no alpha", not match("alpha", gd1));
assert'("Standard form has no + -", not match("\\+ -", gd1));

-- Alpha form: coefficient differs from multinomial
RK2 = K[x_0, x_1, x_2];
F2 = alpha * x_2^3 + 3*x_0*x_1^2*x_2;  -- at d=3
gd2 = genericDForm(F2, 3);
assert'("Alpha form contains \\alpha", match("\\\\alpha", gd2));

-- Negative coefficient
R3 = (ZZ/p)[x_0, x_1, x_2];
F3 = -3*x_0*x_1*x_2 + x_1^3;  -- at d=3, multinomial(1,1,1)=6, ratio=-3/6=-1/2
gd3 = genericDForm(F3, 3);
assert'("Negative term uses - not + -", not match("\\+ -", gd3));

-- ============================================================
stdio << endl << "=== isConnectedSum (intrinsic) ===" << endl;
-- ============================================================

-- Positive cases (must return true)
A422 = (generateGorensteinAlgebra(4, 2, 2, deformed=>false))#0;
assert'("A_{4,2,2} is CS", isConnectedSum(A422));

A623 = (generateGorensteinAlgebra(6, 2, 3, deformed=>false))#0;
assert'("A_{6,2,3} is CS", isConnectedSum(A623));

-- A_{6,2,5} at e=n+4 is removed from the code, construct manually
A625manual = (ZZ/p)[y_1, y_2] / ideal(y_1^2*y_2, y_2^2 - y_1^2);
assert'("A_{6,2,5} manual is CS (hidden decomposition)", isConnectedSum(A625manual));

assert'("A_{5,3,2} is CS", isConnectedSum(pair532#0));

-- Negative cases (must return false)
assert'("A_{3,1,1} is NOT CS", not isConnectedSum(pair311#0));
assert'("A_{5,1,1} is NOT CS", not isConnectedSum(pair511#0));

A624 = (generateGorensteinAlgebra(6, 2, 4, deformed=>false))#0;
assert'("A_{6,2,4} is NOT CS", not isConnectedSum(A624));

A725 = (generateGorensteinAlgebra(7, 2, 5, deformed=>false))#0;
assert'("A_{7,2,5} is NOT CS", not isConnectedSum(A725));

assert'("A_{8,2,8} is NOT CS", not isConnectedSum(pair828#0));

A8310 = (generateGorensteinAlgebra(8, 3, 10, deformed=>false))#0;
assert'("A_{8,3,10} is NOT CS", not isConnectedSum(A8310));

A8314 = (generateGorensteinAlgebra(8, 3, 14, deformed=>false))#0;
assert'("A_{8,3,14} is NOT CS (v^2=0 trap)", not isConnectedSum(A8314));

A839 = (generateGorensteinAlgebra(8, 3, 9, deformed=>false))#0;
assert'("A_{8,3,9} (alpha param) is NOT CS", not isConnectedSum(A839));

-- ============================================================
stdio << endl << "=== detectConnectedSumNaive ===" << endl;
-- ============================================================

comps511 = detectConnectedSumNaive(pair511#0);
assert'("A_{5,1,1} naive: 1 component", #comps511 == 1);

comps422 = detectConnectedSumNaive(A422);
assert'("A_{4,2,2} naive: 2 components", #comps422 == 2);

comps532 = detectConnectedSumNaive(pair532#0);
assert'("A_{5,3,2} naive: 3 components", #comps532 == 3);

comps624 = detectConnectedSumNaive(A624);
assert'("A_{6,2,4} naive: 1 component (y1*y2 != 0)", #comps624 == 1);

comps625 = detectConnectedSumNaive(A625manual);
assert'("A_{6,2,5} naive: 1 component (hidden CS missed)", #comps625 == 1);

-- ============================================================
stdio << endl << "=== subalgebraDim ===" << endl;
-- ============================================================

A522 = pair522#0;
assert'("A_{5,2,2} subalg({y1}) dim = 4", subalgebraDim(A522, {0}) == 4);
assert'("A_{5,2,2} subalg({y2}) dim = 3", subalgebraDim(A522, {1}) == 3);

A734 = (generateGorensteinAlgebra(7, 3, 4, deformed=>false))#0;
assert'("A_{7,3,4} subalg({y1,y2}) dim = 6", subalgebraDim(A734, {0, 1}) == 6);
assert'("A_{7,3,4} subalg({y3}) dim = 3", subalgebraDim(A734, {2}) == 3);

-- ============================================================
stdio << endl << "=== validTypes ===" << endl;
-- ============================================================

vt62 = validTypes(6, 2);
assert'("validTypes(6,2) contains 2", member(2, vt62));
assert'("validTypes(6,2) contains 3", member(3, vt62));
assert'("validTypes(6,2) contains 4", member(4, vt62));
assert'("validTypes(6,2) does NOT contain 5 (removed at e=n+4)", not member(5, vt62));

vt72 = validTypes(7, 2);
assert'("validTypes(7,2) contains 5 (e > n+4)", member(5, vt72));

vt93 = validTypes(9, 3);
assert'("validTypes(9,3) contains 20", member(20, vt93));
assert'("validTypes(9,3) contains 26", member(26, vt93));

vt91 = validTypes(9, 1);
assert'("validTypes(9,1) is empty (n=1 handled separately)", vt91 == {});

-- ============================================================
stdio << endl << "=== identifyComponentType ===" << endl;
-- ============================================================

-- Single-generator components are always type 1
(d1, k1, t1) = identifyComponentType(A522, {0});
assert'("A_{5,2,2} comp {y1}: type = (4,1,1)", d1 == 4 and k1 == 1 and t1 == 1);

(d2, k2, t2) = identifyComponentType(A522, {1});
assert'("A_{5,2,2} comp {y2}: type = (3,1,1)", d2 == 3 and k2 == 1 and t2 == 1);

-- Multi-generator component
(d3, k3, t3) = identifyComponentType(A734, {0, 1});
assert'("A_{7,3,4} comp {y1,y2}: dim=6, ngens=2", d3 == 6 and k3 == 2);
assert'("A_{7,3,4} comp {y1,y2}: type=4", t3 == 4);

-- ============================================================
stdio << endl << "=== areIsomorphicIdeals ===" << endl;
stdio << "  REMOVED: function was broken (false positives for n>=3)." << endl;
stdio << "  Use standalone verify_iso_*.m2 scripts for isomorphism checks." << endl;
-- ============================================================

-- ============================================================
stdio << endl << "========================================" << endl;
stdio << "RESULTS: " << passed << " passed, " << failed << " failed" << endl;
stdio << "========================================" << endl;
if failed > 0 then stdio << "*** SOME TESTS FAILED ***" << endl
else stdio << "All tests passed." << endl;
