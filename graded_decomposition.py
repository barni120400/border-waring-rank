"""
Graded local border Waring decomposition data structure.

A graded decomposition (c, V, q)_d represents:
    ε^a · F + O(ε^{a+1}) = Σ_j c_j ℓ_j^d
where ℓ_j = Σ_i ε^{q_i} V[j,i] x_i.

Convention:
    - V has shape (r+1) x (n+1): rows = summands, columns = variables
    - q[0] = 0 always (x_0 has weight 0)
    - V[:,0] = [1,...,1] for normalized decompositions (v_0 = 1)
    - The moment condition: <c, v^e> = coeff_{x^e}(f) / binom(d, e)
      where v^e = Π_i V[:,i]^{e_i} (entrywise product, then sum with c)
"""

from sympy import (
    Matrix, Symbol, Rational, Integer, binomial, factorial,
    exp, I, pi, cos, sin, simplify, nsimplify, sqrt,
    Poly, degree, LC
)
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from itertools import product as iter_product
from math import gcd
from functools import reduce


def lcm(a, b):
    return abs(a * b) // gcd(a, b)


@dataclass
class GradedDecomposition:
    """A graded local border Waring decomposition (c, V, q)_d.

    c: list of coefficients [c_0, ..., c_r] (SymPy expressions)
    V: matrix where V[j, i] = coefficient of x_i in linear form j
       Shape: (r+1) x (n+1). Row j = summand j, column i = variable x_i.
    q: list of weights [q_0=0, q_1, ..., q_n] (non-negative integers)
    d: degree (SymPy Symbol 'd' for generic, or integer for concrete)
    var_names: optional variable names for display
    """
    c: list
    V: Matrix
    q: list
    d: object
    var_names: Optional[List[str]] = None

    def __post_init__(self):
        assert self.q[0] == 0, "q[0] must be 0 (x_0 has weight 0)"
        assert self.V.rows == len(self.c), f"V rows ({self.V.rows}) must equal len(c) ({len(self.c)})"
        assert self.V.cols == len(self.q), f"V cols ({self.V.cols}) must equal len(q) ({len(self.q)})"

    @property
    def r(self) -> int:
        """Number of summands minus 1 (so r+1 = total summands)."""
        return len(self.c) - 1

    @property
    def n_vars(self) -> int:
        """Number of variables beyond x_0."""
        return self.V.cols - 1

    @property
    def num_summands(self) -> int:
        return len(self.c)

    @property
    def weighted_socle_degree(self) -> int:
        return max(self.q)

    def moment(self, exponent_vec: list):
        """Compute <c, v^e> = Σ_j c_j Π_i V[j,i]^{e_i}.

        exponent_vec: [e_1, ..., e_n] (length n, not n+1 — no x_0 exponent)
        """
        assert len(exponent_vec) == self.n_vars, \
            f"exponent_vec length {len(exponent_vec)} != n_vars {self.n_vars}"
        result = 0
        for j in range(self.num_summands):
            term = self.c[j]
            for i, e in enumerate(exponent_vec):
                if e > 0:
                    term *= self.V[j, i + 1] ** e  # i+1 because col 0 is x_0
            result += term
        return simplify(result)

    def weighted_degree(self, exponent_vec: list) -> int:
        """Compute Σ q_i * e_i for the exponent vector."""
        return sum(self.q[i + 1] * e for i, e in enumerate(exponent_vec))

    def is_normalized(self) -> bool:
        """Check if v_0 = 1 (first column of V is all ones)."""
        return all(self.V[j, 0] == 1 for j in range(self.num_summands))

    def get_var_names(self) -> List[str]:
        if self.var_names:
            return self.var_names
        return [f"x_{i}" for i in range(self.n_vars + 1)]

    def __repr__(self):
        return (f"GradedDecomposition(r={self.r}, n_vars={self.n_vars}, "
                f"q={self.q}, d={self.d}, summands={self.num_summands})")


def root_of_unity(n, k=1):
    """Return ζ_n^k = e^{2πik/n}."""
    if n == 1:
        return 1
    if n == 2:
        return (-1) ** k
    if n == 4:
        return I ** k
    return exp(2 * pi * I * k / n)


def decomposition_for_Gn(n: int) -> GradedDecomposition:
    """Construct the canonical graded decomposition of G_n = K[t]/(t^{n-1}).

    G_n has dim = n, weights (0, 1, 2, ..., n-1), socle degree = n-1.
    The decomposition uses n-th roots of unity with n summands:

    F = -x_0^d + (1/(n-1)) Σ_{k=0}^{n-2} (x_0 + ζ_{n-1}^k ε x_1 + ζ_{n-1}^{2k} ε^2 x_2 + ... + ζ_{n-1}^{(n-2)k} ε^{n-2} x_{n-1})^d

    Wait — looking at the actual example for G_5 (decomp_e5_n1_i1.tex):
    F = -x_0^d + (1/4) Σ_{k=0}^3 (x_0 + ζ_4^k ε x_1 + ζ_4^{2k} ε^2 x_2 + ζ_4^{3k} ε^3 x_3 + ζ_4^{4k} ε^4 x_4)^d

    So for G_n (dim n, socle degree n-1):
    - n summands total: 1 (the -x_0^d term) + (n-1) (the root-of-unity terms)
    - Weights: q = (0, 1, 2, ..., n-1)
    - Coefficients: c_0 = -1, c_k = 1/(n-1) for k = 1, ..., n-1
    - V[0,:] = (1, 0, 0, ..., 0) — just x_0
    - V[k,:] = (1, ζ_{n-1}^{k-1}, ζ_{n-1}^{2(k-1)}, ..., ζ_{n-1}^{(n-1)(k-1)}) for k = 1,...,n-1
      Wait, that gives n-1 variable entries but there should be n variables (x_0 through x_{n-1}).
    """
    dim = n
    socle_deg = n - 1
    num_vars = n  # x_0, x_1, ..., x_{n-1}
    q = list(range(num_vars))  # [0, 1, 2, ..., n-1]

    # From decomp_e5_n1_i1.tex: G_5 has 5 summands
    # c_0 = -1, c_1 = c_2 = c_3 = c_4 = 1/4 (using 4th roots of unity)
    # General: G_n has n summands with (n-1)th roots of unity

    num_summands = n
    c = [Integer(-1)] + [Rational(1, n - 1)] * (n - 1)

    # V matrix: (n) x (n)
    # Row 0: (1, 0, 0, ..., 0) — the x_0^d term
    # Row k (1-indexed): (1, ζ^{k-1}, ζ^{2(k-1)}, ..., ζ^{(n-1)(k-1)})
    # where ζ = ζ_{n-1} = e^{2πi/(n-1)}

    V_data = []
    # Row 0: just x_0
    V_data.append([1] + [0] * (n - 1))
    # Rows 1 through n-1: roots of unity
    for k in range(1, n):
        row = [1]  # x_0 coefficient is always 1 (normalized)
        for i in range(1, n):
            # x_i coefficient in summand k: ζ_{n-1}^{i*(k-1)}
            zeta_power = i * (k - 1)
            row.append(root_of_unity(n - 1, zeta_power))
        V_data.append(row)

    V = Matrix(V_data)
    d = Symbol('d')
    var_names = [f"x_{i}" for i in range(n)]

    return GradedDecomposition(c=c, V=V, q=q, d=d, var_names=var_names)


if __name__ == "__main__":
    # Test G_3
    print("=== G_3 (dim 3) ===")
    D3 = decomposition_for_Gn(3)
    print(D3)
    print(f"Coefficients: {D3.c}")
    print(f"V matrix:\n{D3.V}")
    print(f"Weights: {D3.q}")
    print(f"Is normalized: {D3.is_normalized()}")

    # Test moments: for G_3, f = x_2 + (1/2)x_1^2
    # <c, v^{(0,0,1)}> should give the coefficient of x_2 in f / binom(d,1) = 1/binom(d,1)
    # Actually <c, v^e> = coeff_{x^e}(f) directly (without the binom)
    # For e = (1,0): <c, v_1> = Σ c_j V[j,1] = -1*0 + 1/2*1 + 1/2*(-1) = 0
    # For e = (0,1): <c, v_2> = Σ c_j V[j,2] = -1*0 + 1/2*1 + 1/2*1 = 1
    # For e = (2,0): <c, v_1^2> = Σ c_j V[j,1]^2 = -1*0 + 1/2*1 + 1/2*1 = 1
    m_01 = D3.moment([0, 1])
    m_10 = D3.moment([1, 0])
    m_20 = D3.moment([2, 0])
    print(f"\nMoments:")
    print(f"  <c, v^(1,0)> = {m_10} (expect 0)")
    print(f"  <c, v^(0,1)> = {m_01} (expect 1)")
    print(f"  <c, v^(2,0)> = {m_20} (expect 1)")

    # Test G_5
    print("\n=== G_5 (dim 5) ===")
    D5 = decomposition_for_Gn(5)
    print(D5)
    print(f"Coefficients: {D5.c}")
    print(f"Is normalized: {D5.is_normalized()}")
