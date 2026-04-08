"""
Tests for the connected sum algorithm using the decomposition validator.
"""

import logging
from sympy import simplify, Rational
from graded_decomposition import GradedDecomposition, decomposition_for_Gn
from connected_sum import connected_sum
from validate_decomposition import validate_and_extract, polynomial_to_string

logging.basicConfig(level=logging.WARNING)

passed = 0
failed = 0


def check(name, condition):
    global passed, failed
    if condition:
        print(f"  PASS: {name}")
        passed += 1
    else:
        print(f"  *** FAIL: {name} ***")
        failed += 1


# ============================================================
print("=== Test 1: G_3 # G_3 = A_{4,2,2} ===")
# Expected: f = x_2 + ½x_1² + ½y_1², 4 summands, weights (0,1,2,1)
# ============================================================
D_A = decomposition_for_Gn(3)
D_B = decomposition_for_Gn(3)
result = connected_sum(D_A, D_B)

check("G3#G3: 4 summands", result.num_summands == 4)

valid, poly, errors = validate_and_extract(result, verbose=False)
check("G3#G3: valid decomposition (low-weight = 0)", valid)

# Check polynomial
expected_poly = {(2, 0, 0): 1, (0, 1, 0): 1, (0, 0, 2): 1}
check("G3#G3: polynomial matches f = x_2 + x_1² + y_1²",
      all(simplify(poly.get(k, 0) - v) == 0 for k, v in expected_poly.items())
      and all(k in expected_poly for k in poly))

print(f"  Polynomial: {polynomial_to_string(poly, result.get_var_names())}")

# ============================================================
print("\n=== Test 2: G_4 # G_3 = A_{5,2,2} ===")
# G_4 has dim 4, G_3 has dim 3. Result should have dim 5, summands 4+3-2=5
# Expected: f = f_{A_{5,2,2}} — the canonical form of the connected sum algebra
# ============================================================
D_A = decomposition_for_Gn(4)
D_B = decomposition_for_Gn(3)
result = connected_sum(D_A, D_B)

check("G4#G3: 5 summands", result.num_summands == 5)

valid, poly, errors = validate_and_extract(result, verbose=False)
check("G4#G3: valid decomposition (low-weight = 0)", valid)

if not valid:
    for e, m, w in errors:
        print(f"    Error: <c,v^{e}> = {m} (weight {w})")

print(f"  Polynomial: {polynomial_to_string(poly, result.get_var_names())}")
print(f"  Weights: {result.q}")
print(f"  Variables: {result.get_var_names()}")

# ============================================================
print("\n=== Test 3: G_3 # G_3 # G_3 (iterated) ===")
# G_3 # G_3 = A_{4,2,2} (4 summands)
# (G_3 # G_3) # G_3 should have 4+3-2 = 5 summands
# dim = 3+3+3 - 2*2 = 5
# ============================================================
D_AB = connected_sum(decomposition_for_Gn(3), decomposition_for_Gn(3))
D_C = decomposition_for_Gn(3)
result3 = connected_sum(D_AB, D_C)

check("G3#G3#G3: 5 summands", result3.num_summands == 5)

valid3, poly3, errors3 = validate_and_extract(result3, verbose=False)
check("G3#G3#G3: valid decomposition (low-weight = 0)", valid3)

if not valid3:
    for e, m, w in errors3:
        print(f"    Error: <c,v^{e}> = {m} (weight {w})")

print(f"  Polynomial: {polynomial_to_string(poly3, result3.get_var_names())}")
print(f"  Weights: {result3.q}")
print(f"  Variables: {result3.get_var_names()}")

# ============================================================
print("\n=== Test 4: G_5 # G_3 ===")
# G_5 has dim 5 (5 summands), G_3 has dim 3 (3 summands)
# Result: 5+3-2 = 6 summands, dim = 5+3-2 = 6
# ============================================================
D_A = decomposition_for_Gn(5)
D_B = decomposition_for_Gn(3)
result4 = connected_sum(D_A, D_B)

check("G5#G3: 6 summands", result4.num_summands == 6)

valid4, poly4, errors4 = validate_and_extract(result4, verbose=False)
check("G5#G3: valid decomposition (low-weight = 0)", valid4)

if not valid4:
    for e, m, w in errors4:
        print(f"    Error: <c,v^{e}> = {m} (weight {w})")

print(f"  Polynomial: {polynomial_to_string(poly4, result4.get_var_names())}")
print(f"  Weights: {result4.q}")

# ============================================================
print("\n=== Test 5: G_4 # G_4 ===")
# G_4 has dim 4 (4 summands). Result: 4+4-2 = 6 summands, dim = 6
# ============================================================
D_A = decomposition_for_Gn(4)
D_B = decomposition_for_Gn(4)
result5 = connected_sum(D_A, D_B)

check("G4#G4: 6 summands", result5.num_summands == 6)

valid5, poly5, errors5 = validate_and_extract(result5, verbose=False)
check("G4#G4: valid decomposition (low-weight = 0)", valid5)

if not valid5:
    for e, m, w in errors5:
        print(f"    Error: <c,v^{e}> = {m} (weight {w})")

print(f"  Polynomial: {polynomial_to_string(poly5, result5.get_var_names())}")
print(f"  Weights: {result5.q}")

# ============================================================
print(f"\n{'='*60}")
print(f"RESULTS: {passed} passed, {failed} failed")
print(f"{'='*60}")
