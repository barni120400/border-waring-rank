"""
Connected sum algorithm for graded border Waring decompositions.

Given two graded decompositions D_A and D_B, produces a graded
decomposition of A # B with r_A + r_B summands (saving 2 from naive).

Two-step process:
(a) Compute decomposition for f_A + f_B in DISJOINT variables
    (following the proof of thm:connected-sum-graded-border-wr)
(b) Apply change of variables to match the CANONICAL form f_{A#B}
    (merges socle variables, renames to standard ordering)
"""

import logging
from sympy import Matrix, zeros, eye, simplify, Rational, Symbol, sqrt
from graded_decomposition import GradedDecomposition
from diagonalize import diagonalize, is_diagonalized

logger = logging.getLogger(__name__)


def connected_sum_disjoint(D_A: GradedDecomposition, D_B: GradedDecomposition) -> GradedDecomposition:
    """Step (a): Compute decomposition of f_A + f_B in disjoint variables.

    D_A uses variables x_0, x_1, ..., x_{n_A}
    D_B uses variables x_0, y_1, ..., y_{n_B}
    (x_0 is shared, all others disjoint)

    Output has variables x_0, x_1, ..., x_{n_A}, y_1, ..., y_{n_B}
    and r_A + r_B summands.
    """
    assert D_A.is_normalized(), "D_A must be normalized"
    assert D_B.is_normalized(), "D_B must be normalized"

    n_A = D_A.n_vars
    n_B = D_B.n_vars
    r_A = D_A.r
    r_B = D_B.r

    logger.info(f"Connected sum: D_A({r_A+1} summands, {n_A} vars, q̂={D_A.weighted_socle_degree}) "
                f"# D_B({r_B+1} summands, {n_B} vars, q̂={D_B.weighted_socle_degree})")

    # Step 0a: Rescale weights to common socle degree q̂ = lcm(q̂_A, q̂_B)
    q_hat_A = D_A.weighted_socle_degree
    q_hat_B = D_B.weighted_socle_degree
    from math import gcd
    q_hat = (q_hat_A * q_hat_B) // gcd(q_hat_A, q_hat_B)
    scale_A = q_hat // q_hat_A
    scale_B = q_hat // q_hat_B

    # Rescale weight vectors
    D_A_scaled = GradedDecomposition(
        c=D_A.c, V=D_A.V,
        q=[w * scale_A for w in D_A.q],
        d=D_A.d, var_names=D_A.var_names
    )
    D_B_scaled = GradedDecomposition(
        c=D_B.c, V=D_B.V,
        q=[w * scale_B for w in D_B.q],
        d=D_B.d, var_names=D_B.var_names
    )
    logger.info(f"Rescaled weights to q̂={q_hat}: q_A={D_A_scaled.q}, q_B={D_B_scaled.q}")

    # Step 0b: Diagonalize
    D_A_diag, A0_A = diagonalize(D_A_scaled)
    D_B_diag, A0_B = diagonalize(D_B_scaled)
    assert is_diagonalized(D_A_diag)
    assert is_diagonalized(D_B_diag)
    logger.info(f"Diagonalized. V_A:\n{D_A_diag.V}\nV_B:\n{D_B_diag.V}")

    alpha = list(D_A_diag.c)
    beta = list(D_B_diag.c)
    V_A = D_A_diag.V
    V_B = D_B_diag.V

    # Step 1: Rescale D_A so α_1 = Σ β_j (j>=1)
    S = sum(beta[1:])
    mu = simplify(S / alpha[1])
    alpha = [simplify(mu * a) for a in alpha]
    logger.info(f"Rescaled D_A by μ={mu}. Now α_1={alpha[1]}, Σβ_j={S}")

    # Step 2: Merge x_0^d
    alpha_0_new = simplify(alpha[0] + beta[0])
    logger.info(f"Merged x_0^d: α_0={alpha_0_new}")

    # Step 3: Eliminate L_1, compensate
    alpha_0_new = simplify(alpha_0_new + alpha[1])
    logger.info(f"Eliminated L_1. α_0={alpha_0_new}")

    # Step 4: Assemble
    # Variables: x_0, x_1,...,x_{n_A}, y_1,...,y_{n_B} — ALL disjoint (no socle merge)
    n_total = n_A + n_B + 1  # total columns including x_0
    num_summands = r_A + r_B

    c_fin = [alpha_0_new] + alpha[2:] + beta[1:]

    V_fin = zeros(num_summands, n_total)

    # Summand 0: x_0 only
    V_fin[0, 0] = 1

    # Summands 1..r_A-1: from A's summands 2..r_A (A-variables in cols 0..n_A)
    for idx, a_idx in enumerate(range(2, r_A + 1)):
        row = idx + 1
        for col in range(n_A + 1):
            V_fin[row, col] = V_A[a_idx, col]

    # Summands r_A..r_A+r_B-1: from B's summands 1..r_B
    # B-variables y_1,...,y_{n_B} go to cols n_A+1,...,n_A+n_B
    # Plus x_1 shift from Step 3 (col 1)
    for idx, b_idx in enumerate(range(1, r_B + 1)):
        row = (r_A - 1) + idx + 1
        # x_0 from B
        V_fin[row, 0] = V_B[b_idx, 0]
        # x_1 shift (the compensation from eliminating L_1)
        V_fin[row, 1] = 1
        # y_1,...,y_{n_B} from B (B's cols 1..n_B → combined cols n_A+1..n_A+n_B)
        for b_col in range(1, n_B + 1):
            V_fin[row, n_A + b_col] = V_B[b_idx, b_col]

    # Weights: [0, q_A[1],...,q_A[n_A], q_B[1],...,q_B[n_B]]
    q_fin = D_A_diag.q + D_B_diag.q[1:]

    # Variable names
    var_names = [f"x_{i}" for i in range(n_A + 1)] + [f"y_{j}" for j in range(1, n_B + 1)]

    logger.info(f"Assembled: {num_summands} summands, {n_total} vars, q={q_fin}")
    logger.info(f"V:\n{V_fin}")

    # Step 5: Revert diagonalization and μ-rescaling
    # Build combined revert matrix acting on columns of V_fin
    q_hat_A = D_A_diag.weighted_socle_degree
    A0_A_inv = A0_A.inv()
    A0_B_inv = A0_B.inv()

    # μ-revert for A-variables: x_i → μ^{-q_i/q̂_A} x_i
    mu_revert_A = eye(n_A + 1)
    for i in range(1, n_A + 1):
        power = Rational(D_A_diag.q[i], q_hat_A)
        mu_revert_A[i, i] = simplify(mu ** (-power))

    revert_A = simplify(mu_revert_A * A0_A_inv)

    # Combined revert: block diagonal [revert_A | A0_B_inv_block]
    revert = eye(n_total)
    # A-block: cols 0..n_A
    for i in range(n_A + 1):
        for j in range(n_A + 1):
            revert[i, j] = revert_A[i, j]
    # B-block: cols n_A+1..n_A+n_B correspond to y_1..y_{n_B}
    # A0_B_inv acts on (x_0, y_1,...,y_{n_B}), we only need the y-part
    for i in range(1, n_B + 1):
        for j in range(1, n_B + 1):
            revert[n_A + i, n_A + j] = A0_B_inv[i, j]

    V_reverted = simplify(V_fin * revert)
    logger.info(f"Reverted V:\n{V_reverted}")

    return GradedDecomposition(c=c_fin, V=V_reverted, q=q_fin, d=D_A.d, var_names=var_names)


def drop_socle_B(decomp: GradedDecomposition, n_A: int, n_B: int) -> GradedDecomposition:
    """Step (b): Set y_{n_B} = 0 (drop B's socle variable).

    The socle is already represented by x_{n_A} in f_A.
    Setting y_{n_B} = 0 removes the duplicate socle contribution from f_B.
    This is just deleting the y_{n_B} column from V.
    """
    col_socle_B = n_A + n_B  # y_{n_B} column (last B-variable)

    V = decomp.V[:, :col_socle_B].row_join(decomp.V[:, col_socle_B + 1:])
    q_new = decomp.q[:col_socle_B] + decomp.q[col_socle_B + 1:]
    var_names_new = decomp.var_names[:col_socle_B] + decomp.var_names[col_socle_B + 1:]

    return GradedDecomposition(c=decomp.c, V=V, q=q_new, d=decomp.d, var_names=var_names_new)


def sort_variables_by_weight(decomp: GradedDecomposition) -> GradedDecomposition:
    """Sort variables by weight and rename to x_0, x_1, x_2, ...

    Permutes columns of V to put weights in non-decreasing order.
    This is a change of variables (column permutation) which is its
    own inverse — no separate revert needed since it's just renaming.

    For variables at the same weight, the original order is preserved
    (stable sort).
    """
    q = decomp.q
    n = decomp.n_vars

    # Build sort permutation: indices sorted by weight (stable)
    # Column 0 (x_0, weight 0) always stays first
    perm = [0] + sorted(range(1, n + 1), key=lambda i: q[i])

    # Permute columns of V
    V_sorted = decomp.V[:, perm]

    # Permute weights
    q_sorted = [q[i] for i in perm]

    # New variable names: x_0, x_1, ..., x_n
    var_names = [f"x_{i}" for i in range(n + 1)]

    return GradedDecomposition(c=decomp.c, V=V_sorted, q=q_sorted, d=decomp.d, var_names=var_names)


def connected_sum(D_A: GradedDecomposition, D_B: GradedDecomposition) -> GradedDecomposition:
    """Full connected sum: compute decomposition of canonical f_{A#B}.

    1. Compute decomposition for f_A + f_B (disjoint variables)
    2. Drop B's socle variable (y_{n_B} = 0)
    3. Sort variables by weight, rename to x_0, x_1, ..., x_n
    """
    n_A = D_A.n_vars
    n_B = D_B.n_vars

    # Step (a): disjoint variable decomposition
    D_disjoint = connected_sum_disjoint(D_A, D_B)

    # Step (b): drop B's socle variable (y_{n_B} = 0)
    D_dropped = drop_socle_B(D_disjoint, n_A, n_B)

    # Step (c): sort variables by weight, rename to x_0, x_1, ...
    D_sorted = sort_variables_by_weight(D_dropped)

    return D_sorted


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    from graded_decomposition import decomposition_for_Gn

    print("=" * 60)
    print("G_3 # G_3 = A_{4,2,2}")
    print("=" * 60)
    D_A = decomposition_for_Gn(3)
    D_B = decomposition_for_Gn(3)
    result = connected_sum(D_A, D_B)

    print(f"\nResult: {result}")
    print(f"Coefficients: {result.c}")
    print(f"V matrix:\n{result.V}")
    print(f"Weights: {result.q}")
    print(f"Variables: {result.get_var_names()}")
    print(f"Summands: {result.num_summands} (expected 4)")

    # Expected canonical form: f = x_2 + ½x_1² + ½y_1²
    # Variables: x_0(w=0), x_1(w=1), x_2(w=2, socle), y_1(w=1)
    # Only check moments at weight ≤ socle degree (= 2)
    print(f"\nMoment verification (f = x_2 + ½x_1² + ½y_1²):")
    print(f"Variables: {result.get_var_names()}, weights: {result.q}")
    for e, expected in [([1,0,0], 0), ([0,1,0], 1), ([0,0,1], 0),
                         ([2,0,0], 1), ([0,0,2], 1)]:
        m = simplify(result.moment(e))
        status = "✓" if m == expected else f"✗ (got {m})"
        print(f"  <c, v^{e}> = {m}  expected {expected}  {status}")
