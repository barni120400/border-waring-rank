"""
Diagonalization of graded border Waring decompositions.

Implements Algorithm 1 from "Fixed parameter debordering of Waring rank"
(arXiv:2502.03150). Transforms a normalized decomposition into diagonalized
form where each variable x_k "owns" one linear form L_k:

  L_0 = x_0
  L_k = x_0 + ε^{q_k} x_k + (lower ε-order terms in x_1,...,x_{k-1})

The algorithm applies a change of variables A_0 (invertible, no ε) to the
x-variables via column operations on V. The original decomposition is for f,
the diagonalized one is for f(A_0 x). We track A_0 so it can be reverted.
"""

from sympy import Matrix, zeros, eye, simplify
from graded_decomposition import GradedDecomposition


def is_diagonalized(decomp: GradedDecomposition) -> bool:
    """Check if a decomposition is in diagonalized form.

    Diagonalized means:
    - Row 0 of V: (1, 0, 0, ..., 0) — L_0 = x_0
    - For k = 1,...,n: V[k, k] = 1 (diagonal entry)
    - For k = 1,...,n and i > k: if q[i] >= q[k] then V[k, i] = 0
      (upper triangle is zero for variables at equal or higher weight)

    Note: V[k, i] may be nonzero when q[i] < q[k] (variable i has
    strictly lower weight than k). This is the lower-triangular
    structure in the WEIGHT ordering, not necessarily the index ordering.
    """
    V = decomp.V
    n = decomp.n_vars
    q = decomp.q

    # Check row 0
    if V[0, 0] != 1:
        return False
    for i in range(1, V.cols):
        if V[0, i] != 0:
            return False

    # Check rows 1..n
    for k in range(1, min(n + 1, V.rows)):
        if V[k, k] != 1:
            return False
        for i in range(k + 1, V.cols):
            if q[i] >= q[k] and V[k, i] != 0:
                return False

    return True


def diagonalize(decomp: GradedDecomposition):
    """Diagonalize a normalized decomposition.

    Returns (diag_decomp, A_0) where:
    - diag_decomp: diagonalized GradedDecomposition
    - A_0: the invertible change-of-variables matrix (constant part)
           The diagonalized decomposition is for f(A_0 x).

    Method: Gaussian elimination on V using COLUMN operations only
    (= changes of x-variables). Process columns left to right, grouped
    by weight. For each variable x_k:
    1. Find an unassigned row with nonzero V[j, k]
    2. Swap that row to position k
    3. Eliminate V[assigned_row, k] for all assigned rows above k
    4. Scale column k so V[k, k] = 1
    Track all column operations in A_0.
    """
    V = Matrix(decomp.V)  # ensure mutable copy
    c = list(decomp.c)
    q = list(decomp.q)
    n = decomp.n_vars

    # A_0 tracks cumulative column operations: V_new = V_orig * A_0
    A_0 = eye(n + 1)

    # Ensure row 0 is (1, 0, ..., 0) — the pure x_0 summand
    if not all(V[0, i] == 0 for i in range(1, n + 1)):
        # Find a row that is (1, 0, ..., 0)
        found = False
        for j in range(V.rows):
            if V[j, 0] == 1 and all(V[j, i] == 0 for i in range(1, n + 1)):
                if j != 0:
                    V.row_swap(0, j)
                    c[0], c[j] = c[j], c[0]
                found = True
                break
        if not found:
            raise ValueError("No pure x_0 summand found")

    # assigned_rows: set of row indices already assigned to a variable
    assigned_rows = {0}

    # Group variables by weight
    weight_groups = {}
    for i in range(1, n + 1):
        w = q[i]
        if w not in weight_groups:
            weight_groups[w] = []
        weight_groups[w].append(i)

    for w in sorted(weight_groups.keys()):
        for var_idx in weight_groups[w]:
            # Step 1: Find pivot row — an unassigned row with nonzero V[j, var_idx]
            pivot_row = None
            for j in range(V.rows):
                if j not in assigned_rows and V[j, var_idx] != 0:
                    pivot_row = j
                    break

            if pivot_row is None:
                raise ValueError(
                    f"No pivot for x_{var_idx} (weight {w})")

            # Step 2: Swap pivot row into position var_idx
            if pivot_row != var_idx:
                V.row_swap(var_idx, pivot_row)
                c[var_idx], c[pivot_row] = c[pivot_row], c[var_idx]
                if pivot_row in assigned_rows:
                    assigned_rows.discard(pivot_row)
                    assigned_rows.add(var_idx)

            # Step 3: Eliminate entries — zero out V[k, var_idx]
            # for all assigned rows k where q[k] <= q[var_idx]
            # (eliminate variables at same or lower weight)
            for k in sorted(assigned_rows):
                if k == var_idx:
                    continue
                if q[k] > q[var_idx]:
                    continue  # skip higher-weight variables
                entry = V[k, var_idx]
                if entry != 0:
                    for j in range(V.rows):
                        V[j, var_idx] = simplify(V[j, var_idx] - entry * V[j, k])
                    for i in range(n + 1):
                        A_0[i, var_idx] = simplify(A_0[i, var_idx] - entry * A_0[i, k])

            # Step 4: Scale column var_idx so V[var_idx, var_idx] = 1
            pivot_val = V[var_idx, var_idx]
            if pivot_val == 0:
                raise ValueError(
                    f"Zero pivot at ({var_idx}, {var_idx}) after elimination")
            if pivot_val != 1:
                scale = 1 / pivot_val
                for j in range(V.rows):
                    V[j, var_idx] = simplify(V[j, var_idx] * scale)
                for i in range(n + 1):
                    A_0[i, var_idx] = simplify(A_0[i, var_idx] * scale)

            assigned_rows.add(var_idx)

    return GradedDecomposition(c=c, V=V, q=q, d=decomp.d, var_names=decomp.var_names), A_0


if __name__ == "__main__":
    from graded_decomposition import decomposition_for_Gn

    for n in [3, 4, 5]:
        print(f"\n=== Diagonalizing G_{n} ===")
        D = decomposition_for_Gn(n)
        print(f"Original V:\n{D.V}")
        print(f"Is diagonalized: {is_diagonalized(D)}")

        D_diag, A0 = diagonalize(D)
        print(f"Diagonalized V:\n{D_diag.V}")
        print(f"Is diagonalized: {is_diagonalized(D_diag)}")
        print(f"A_0:\n{A0}")
        print(f"Coefficients: {D_diag.c}")

        # Verify moments preserved
        print(f"Moment <c, v^(1,0,...)>: {D_diag.moment([1] + [0]*(n-2))}")
        print(f"Moment <c, v^(0,...,0,1)>: {D_diag.moment([0]*(n-2) + [1])}")
