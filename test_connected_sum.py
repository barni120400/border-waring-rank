"""
Tests for the connected sum algorithm using the decomposition validator.
"""

import logging
from sympy import simplify
from graded_decomposition import decomposition_for_Gn
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
print("=== Test 1: G_3 # G_3 ===")
# f = x_3 + x_1² + x_2², 4 summands, weights (0,1,1,2)
# ============================================================
result = connected_sum(decomposition_for_Gn(3), decomposition_for_Gn(3))
check("G3#G3: 4 summands", result.num_summands == 4)
check("G3#G3: weights (0,1,1,2)", result.q == [0, 1, 1, 2])
valid, poly, _ = validate_and_extract(result, verbose=False)
check("G3#G3: valid (low-weight = 0)", valid)
expected = {(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 1): 1}
check("G3#G3: f = x_3 + x_1² + x_2²",
      all(simplify(poly.get(k, 0) - v) == 0 for k, v in expected.items())
      and all(k in expected for k in poly))
print(f"  Polynomial: {polynomial_to_string(poly, result.get_var_names())}")

# ============================================================
print("\n=== Test 2: G_4 # G_3 ===")
# 5 summands, weights should be sorted
# ============================================================
result2 = connected_sum(decomposition_for_Gn(4), decomposition_for_Gn(3))
check("G4#G3: 5 summands", result2.num_summands == 5)
valid2, poly2, errors2 = validate_and_extract(result2, verbose=False)
check("G4#G3: valid (low-weight = 0)", valid2)
if not valid2:
    for e, m, w in errors2[:3]:
        print(f"    Error: e={e} m={m} w={w}")
print(f"  Polynomial: {polynomial_to_string(poly2, result2.get_var_names())}")
print(f"  Weights: {result2.q}")

# ============================================================
print("\n=== Test 3: G_3 # G_3 # G_3 (iterated) ===")
# f = x_4 + x_1² + x_2² + x_3², 5 summands, weights (0,1,1,1,2)
# ============================================================
D33 = connected_sum(decomposition_for_Gn(3), decomposition_for_Gn(3))
result3 = connected_sum(D33, decomposition_for_Gn(3))
check("G3#G3#G3: 5 summands", result3.num_summands == 5)
check("G3#G3#G3: weights (0,1,1,1,2)", result3.q == [0, 1, 1, 1, 2])
valid3, poly3, errors3 = validate_and_extract(result3, verbose=False)
check("G3#G3#G3: valid (low-weight = 0)", valid3)
expected3 = {(2, 0, 0, 0): 1, (0, 2, 0, 0): 1, (0, 0, 2, 0): 1, (0, 0, 0, 1): 1}
check("G3#G3#G3: f = x_4 + x_1² + x_2² + x_3²",
      all(simplify(poly3.get(k, 0) - v) == 0 for k, v in expected3.items())
      and all(k in expected3 for k in poly3))
print(f"  Polynomial: {polynomial_to_string(poly3, result3.get_var_names())}")

# ============================================================
print("\n=== Test 4: G_5 # G_3 ===")
# 6 summands
# ============================================================
result4 = connected_sum(decomposition_for_Gn(5), decomposition_for_Gn(3))
check("G5#G3: 6 summands", result4.num_summands == 6)
valid4, poly4, errors4 = validate_and_extract(result4, verbose=False)
check("G5#G3: valid (low-weight = 0)", valid4)
if not valid4:
    for e, m, w in errors4[:3]:
        print(f"    Error: e={e} m={m} w={w}")
print(f"  Polynomial: {polynomial_to_string(poly4, result4.get_var_names())}")
print(f"  Weights: {result4.q}")

# ============================================================
print("\n=== Test 5: G_4 # G_4 ===")
# 6 summands
# ============================================================
result5 = connected_sum(decomposition_for_Gn(4), decomposition_for_Gn(4))
check("G4#G4: 6 summands", result5.num_summands == 6)
valid5, poly5, errors5 = validate_and_extract(result5, verbose=False)
check("G4#G4: valid (low-weight = 0)", valid5)
if not valid5:
    for e, m, w in errors5[:3]:
        print(f"    Error: e={e} m={m} w={w}")
print(f"  Polynomial: {polynomial_to_string(poly5, result5.get_var_names())}")
print(f"  Weights: {result5.q}")

# ============================================================
print("\n=== Test 6: G_3 # G_3 # G_3 # G_3 (4-fold) ===")
# 6 summands, dim = 4*3 - 3*2 = 6
# ============================================================
D33 = connected_sum(decomposition_for_Gn(3), decomposition_for_Gn(3))
D333 = connected_sum(D33, decomposition_for_Gn(3))
result6 = connected_sum(D333, decomposition_for_Gn(3))
check("G3^#4: 6 summands", result6.num_summands == 6)
valid6, poly6, errors6 = validate_and_extract(result6, verbose=False)
check("G3^#4: valid (low-weight = 0)", valid6)
expected6 = {(2,0,0,0,0): 1, (0,2,0,0,0): 1, (0,0,2,0,0): 1, (0,0,0,2,0): 1, (0,0,0,0,1): 1}
check("G3^#4: f = x_5 + x_1² + x_2² + x_3² + x_4²",
      all(simplify(poly6.get(k, 0) - v) == 0 for k, v in expected6.items())
      and all(k in expected6 for k in poly6))
print(f"  Polynomial: {polynomial_to_string(poly6, result6.get_var_names())}")

# ============================================================
print(f"\n{'='*60}")
print(f"RESULTS: {passed} passed, {failed} failed")
print(f"{'='*60}")
