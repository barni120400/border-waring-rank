"""
Connected sum algorithm for graded border Waring decompositions.

Given two graded decompositions D_A and D_B, produces a graded
decomposition of A # B with r_A + r_B summands (saving 2 from naive).

Follows the proof of thm:connected-sum-graded-border-wr in connected sum.tex.

The algorithm:
0. Diagonalize both inputs (track change of basis A_0)
1. Rescale D_A so α_1 = Σ β_j (multiply all c_A by μ)
2. Merge x_0^d: α_0 ← μα_0 + β_0, remove R_0
3. Eliminate L_1: remove, update α_0 ← α_0 + α_1, shift R_j ← R_j + ε^{w_1} x_1
4. Assemble with shared socle (x_{n_A} = y_{n_B}, delete y_{n_B})
5. Revert diagonalization and μ-rescaling

Output variables: x_0, x_1, ..., x_{n_A}, y_1, ..., y_{n_B-1}
"""

import logging
from sympy import Matrix, zeros, eye, simplify, Rational, Symbol
from graded_decomposition import GradedDecomposition
from diagonalize import diagonalize, is_diagonalized
from math import gcd

logger = logging.getLogger(__name__)


def connected_sum(D_A: GradedDecomposition, D_B: GradedDecomposition) -> GradedDecomposition:
    """Compute the connected sum decomposition D_A # D_B.

    Input: Two graded decompositions (must be normalized: v_0 = 1).
    Output: A GradedDecomposition for A # B with r_A + r_B summands.

    The output is in variables x_0, x_1,...,x_{n_A}, y_1,...,y_{n_B-1}
    where x_0 is shared and x_{n_A} = y_{n_B} (shared socle, y_{n_B} deleted).
    """
    assert D_A.is_normalized(), "D_A must be normalized (v_0 = 1)"
    assert D_B.is_normalized(), "D_B must be normalized (v_0 = 1)"

    n_A = D_A.n_vars  # number of A-variables beyond x_0
    n_B = D_B.n_vars  # number of B-variables beyond x_0
    r_A = D_A.r
    r_B = D_B.r
    d = D_A.d  # degree (should be same for both, or use Symbol)

    logger.info(f"Connected sum: D_A has {r_A+1} summands, {n_A} vars; "
                f"D_B has {r_B+1} summands, {n_B} vars")

    # ============================================================
    # Step 0: Diagonalize both
    # ============================================================
    logger.info("Step 0: Diagonalizing...")
    D_A_diag, A0_A = diagonalize(D_A)
    D_B_diag, A0_B = diagonalize(D_B)
    assert is_diagonalized(D_A_diag), "Diagonalization of D_A failed"
    assert is_diagonalized(D_B_diag), "Diagonalization of D_B failed"
    logger.info(f"  D_A diagonalized: V =\n{D_A_diag.V}")
    logger.info(f"  D_B diagonalized: V =\n{D_B_diag.V}")

    # ============================================================
    # Step 1: Rescale D_A so α_1 = Σ_{j=1}^{r_B} β_j
    # ============================================================
    alpha = list(D_A_diag.c)
    beta = list(D_B_diag.c)
    V_A = D_A_diag.V.copy()
    V_B = D_B_diag.V.copy()

    S = sum(beta[1:])  # Σ β_j for j >= 1
    mu = simplify(S / alpha[1])
    logger.info(f"Step 1: Rescaling D_A by μ = {mu}")
    logger.info(f"  α_1 = {alpha[1]}, Σβ_j = {S}")
    alpha = [simplify(mu * a) for a in alpha]
    # Now alpha[1] = S = Σ β_j
    # The decomposition now represents μ·f_A (not f_A)

    # ============================================================
    # Step 2: Merge x_0^d terms
    # ============================================================
    # Both have L_0 = R_0 = x_0. Combine α_0 + β_0, drop R_0.
    alpha_0_new = simplify(alpha[0] + beta[0])
    logger.info(f"Step 2: Merged x_0^d: α_0 = {alpha[0]} + β_0 = {beta[0]} → {alpha_0_new}")

    # ============================================================
    # Step 3: Eliminate L_1 and compensate
    # ============================================================
    # Remove α_1·L_1^d. Compensate:
    #   (a) x_0^d part: α_0 ← α_0 + α_1
    #   (b) x_1 part: for each R_j (j=1..r_B), add 1 to the x_1 coordinate
    alpha_0_new = simplify(alpha_0_new + alpha[1])
    logger.info(f"Step 3: Eliminated L_1. New α_0 = {alpha_0_new}")
    logger.info(f"  Shifting all R_j by x_1 (adding ε^{{w_1}} x_1 to each B-form)")

    # ============================================================
    # Step 4: Assemble the combined decomposition
    # ============================================================
    # Summands from A: index 0 (x_0, with new α_0), indices 2..r_A
    # Summands from B: indices 1..r_B (with x_1 shift)
    # Total: 1 + (r_A - 1) + r_B = r_A + r_B

    # Coefficients
    c_fin = [alpha_0_new] + alpha[2:] + beta[1:]
    logger.info(f"Step 4: Final coefficients: {c_fin}")
    logger.info(f"  {1 + (r_A - 1) + r_B} = {len(c_fin)} summands (should be {r_A + r_B})")

    # Variables: x_0, x_1, ..., x_{n_A}, y_1, ..., y_{n_B-1}
    # Note: x_{n_A} = y_{n_B} (shared socle), y_{n_B} is deleted
    # Total variables: 1 + n_A + (n_B - 1) = n_A + n_B
    n_combined = n_A + n_B  # number of variables (including x_0)
    num_summands = r_A + r_B

    # Weight vector
    q_A = D_A_diag.q
    q_B = D_B_diag.q
    # x_0 has weight 0, then A-weights, then B-weights (excluding y_{n_B})
    q_fin = q_A + q_B[1:-1]  # [0, q_A1, ..., q_A_{n_A}, q_B1, ..., q_B_{n_B-1}]
    logger.info(f"  Combined weights: {q_fin}")

    # Build the V matrix
    # Rows = summands, Cols = variables
    # Col indices: 0=x_0, 1..n_A = x_1..x_{n_A}, n_A+1..n_A+n_B-1 = y_1..y_{n_B-1}

    V_fin = zeros(num_summands, n_combined)

    # Summand 0: x_0 (the merged term) — just x_0
    V_fin[0, 0] = 1

    # Summands 1..(r_A-1): from A's summands 2..r_A
    for idx, a_idx in enumerate(range(2, r_A + 1)):
        row = idx + 1
        # Copy A's V entries for x_0, x_1, ..., x_{n_A}
        for col in range(n_A + 1):
            V_fin[row, col] = V_A[a_idx, col]

    # Summands r_A..(r_A+r_B-1): from B's summands 1..r_B
    for idx, b_idx in enumerate(range(1, r_B + 1)):
        row = (r_A - 1) + idx + 1  # offset by the A summands
        # x_0 entry from B
        V_fin[row, 0] = V_B[b_idx, 0]
        # x_1 entry: add 1 (the shift from Step 3)
        V_fin[row, 1] = 1
        # y_1, ..., y_{n_B-1} entries from B (columns 1..n_B-1 of V_B)
        for b_col in range(1, n_B):  # B-columns 1 to n_B-1 (excluding socle col n_B)
            combined_col = n_A + b_col  # position in combined V
            V_fin[row, combined_col] = V_B[b_idx, b_col]
        # Socle: y_{n_B} = x_{n_A}, so V_B[b_idx, n_B] goes to column n_A
        V_fin[row, n_A] = simplify(V_fin[row, n_A] + V_B[b_idx, n_B])

    logger.info(f"  Combined V:\n{V_fin}")

    # Variable names
    var_names = [f"x_{i}" for i in range(n_A + 1)] + [f"y_{j}" for j in range(1, n_B)]

    # ============================================================
    # Step 5: Revert changes of variables
    # ============================================================
    # The diagonalized decomposition is for f_A(A0_A x) and f_B(A0_B y).
    # After the μ-rescaling, it's for μ·f_A(A0_A x) + f_B(A0_B y).
    # To get f_A + f_B in original variables:
    #   (a) Revert A0_A on x-variables: x → A0_A^{-1} x
    #   (b) Revert A0_B on y-variables: y → A0_B^{-1} y
    #   (c) Revert μ-rescaling on A-variables:
    #       x_i → μ^{-q_i/q̂_A} x_i for i = 1,...,n_A
    #       (weight-preserving diagonal scaling)
    #
    # Build the combined revert matrix:
    # It acts on columns of V_fin.

    q_hat_A = D_A_diag.weighted_socle_degree

    # Revert diagonalization: A0_A^{-1} for A-cols, A0_B^{-1} for B-cols
    A0_A_inv = A0_A.inv()
    A0_B_inv = A0_B.inv()

    # Build the μ-revert diagonal for A-variables
    # x_i → μ^{-q_i/q̂_A} x_i
    # Since μ may not have clean roots, we absorb this into A0_A_inv
    mu_revert = eye(n_A + 1)
    for i in range(1, n_A + 1):
        power = Rational(q_A[i], q_hat_A)
        mu_revert[i, i] = simplify(mu ** (-power))
    logger.info(f"  μ-revert diagonal:\n{mu_revert}")

    # Combined revert for A-variables: mu_revert * A0_A_inv
    revert_A = simplify(mu_revert * A0_A_inv)

    # Build the combined revert matrix (block diagonal)
    # Cols 0..n_A use revert_A, cols n_A+1..end use A0_B_inv (excluding x_0 and socle)
    revert = eye(n_combined)
    # A-block: cols 0..n_A
    for i in range(n_A + 1):
        for j in range(n_A + 1):
            revert[i, j] = revert_A[i, j]
    # B-block: cols n_A+1..n_combined-1 correspond to y_1..y_{n_B-1}
    # A0_B_inv acts on (x_0, y_1, ..., y_{n_B}) but we only need y_1..y_{n_B-1}
    # The B-revert ignores x_0 (handled by A) and socle (identified with x_{n_A})
    for i in range(1, n_B):  # y_1 to y_{n_B-1}
        for j in range(1, n_B):  # same range
            revert[n_A + i, n_A + j] = A0_B_inv[i, j]

    # Apply revert: V_fin_reverted = V_fin * revert^T
    # (column operations = right-multiply V by revert)
    V_reverted = simplify(V_fin * revert)
    logger.info(f"Step 5: Reverted V:\n{V_reverted}")

    return GradedDecomposition(
        c=c_fin,
        V=V_reverted,
        q=q_fin,
        d=d,
        var_names=var_names
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    from graded_decomposition import decomposition_for_Gn

    print("=== G_3 # G_3 ===")
    D_A = decomposition_for_Gn(3)
    D_B = decomposition_for_Gn(3)
    result = connected_sum(D_A, D_B)
    print(f"\nResult: {result}")
    print(f"Coefficients: {result.c}")
    print(f"V matrix:\n{result.V}")
    print(f"Weights: {result.q}")
    print(f"Variable names: {result.get_var_names()}")
    print(f"Expected: 4 summands (3 + 3 - 2 = 4)")
    print(f"Actual: {result.num_summands} summands")

    # Verify moments
    print(f"\nMoments:")
    # For G_3 # G_3 = A_{4,2,2}: f = x_2 + (1/2)x_1^2 + y_1^2/2 + x_2*y_1^...
    # Actually f_A + f_B = (x_2 + (1/2)x_1^2) + (y_2 + (1/2)y_1^2)
    # But x_2 = y_2 (shared socle), so f = x_2 + (1/2)x_1^2 + (1/2)y_1^2
    # Moments: <c, v^(1,0,0)> = 0 (x_1 coeff), <c, v^(0,1,0)> = 1 (x_2 coeff)
    #          <c, v^(0,0,1)> = 0 (y_1 coeff), <c, v^(2,0,0)> = 1 (x_1^2 coeff)
    for e in [[1,0,0], [0,1,0], [0,0,1], [2,0,0], [0,0,2]]:
        print(f"  <c, v^{e}> = {result.moment(e)}")
