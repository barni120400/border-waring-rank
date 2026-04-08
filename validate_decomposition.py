"""
Validate a graded border Waring decomposition and extract its polynomial.

For a graded decomposition (c, V, q)_d with socle degree q̂ = max(q):
- All moments <c, v^e> at weight w(e) < q̂ must be 0
- Moments at weight w(e) = q̂ give the coefficients of the polynomial f
- The polynomial f is: f = Σ_e <c, v^e> · x^e / e!  (where the sum is over
  exponent vectors e with w(e) = q̂)

Actually, from rem:coeff-formula:
  coeff_{x_0^{d-|e|} x^e}(F) = binom(d, e) · <c, v^e>
So f (dehomogenized, x_0=1) has:
  coeff_{x^e}(f) = binom(d, e) · <c, v^e> / binom(d, |e|) ...

Wait — more simply: <c, v^e> is the moment, and the polynomial is
  f = Σ_{w(e)=q̂} <c, v^e> · x^e
where x^e = x_1^{e_1} ... x_n^{e_n} (no factorials, no binomials).
The binomial factors appear in the HOMOGENIZED form F = H_d(f).
"""

from sympy import Symbol, symbols, simplify, Rational, Poly, nsimplify, N, Abs
from itertools import product as iter_product
from graded_decomposition import GradedDecomposition


def enumerate_exponents(n_vars: int, max_total: int):
    """Enumerate all exponent vectors (e_1,...,e_n) with e_i >= 0 and sum <= max_total."""
    if n_vars == 0:
        yield []
        return
    for e0 in range(max_total + 1):
        for rest in enumerate_exponents(n_vars - 1, max_total - e0):
            yield [e0] + rest


def validate_and_extract(decomp: GradedDecomposition, verbose: bool = True):
    """Validate a decomposition and extract the polynomial it represents.

    Returns (is_valid, polynomial_dict, errors) where:
    - is_valid: True if all low-weight moments are 0
    - polynomial_dict: {exponent_tuple: coefficient} for the polynomial f
    - errors: list of (exponent, moment, weight) for failed validations
    """
    n = decomp.n_vars
    q = decomp.q
    q_hat = decomp.weighted_socle_degree

    # We need to check all exponent vectors e with weight w(e) <= q_hat
    # and sum(e) reasonable (bounded by q_hat / min(q_i for q_i > 0))
    min_weight = min(w for w in q[1:] if w > 0)
    max_total_degree = q_hat // min_weight + 1

    errors = []
    polynomial = {}

    for e in enumerate_exponents(n, max_total_degree):
        if all(ei == 0 for ei in e):
            continue  # skip the zero exponent

        w = sum(q[i + 1] * e[i] for i in range(n))

        if w > q_hat:
            continue  # beyond socle degree, not checked

        moment_raw = decomp.moment(e)
        # Try symbolic simplification, then numerical if needed
        moment = simplify(moment_raw)
        if moment != 0:
            # SymPy sometimes can't simplify sums of roots of unity
            # Check numerically
            num_val = complex(N(moment))
            if abs(num_val) < 1e-10:
                moment = 0

        if w < q_hat:
            # Low-weight: must be 0
            if moment != 0:
                errors.append((tuple(e), moment, w))
                if verbose:
                    print(f"  FAIL: <c, v^{e}> = {moment} (weight {w} < {q_hat}, should be 0)")
        else:
            # At socle weight: coefficient of f
            if moment != 0:
                polynomial[tuple(e)] = moment
                if verbose:
                    print(f"  f coeff x^{e} = {moment} (weight {w})")

    is_valid = len(errors) == 0
    if verbose:
        if is_valid:
            print(f"  VALID: all {q_hat}-1 low-weight moments are 0")
        else:
            print(f"  INVALID: {len(errors)} low-weight moments are nonzero")

    return is_valid, polynomial, errors


def polynomial_to_string(poly_dict: dict, var_names: list) -> str:
    """Convert polynomial dict to readable string."""
    terms = []
    for e, coeff in sorted(poly_dict.items()):
        monomial_parts = []
        for i, ei in enumerate(e):
            if ei > 0:
                name = var_names[i + 1] if i + 1 < len(var_names) else f"x_{i+1}"
                if ei == 1:
                    monomial_parts.append(name)
                else:
                    monomial_parts.append(f"{name}^{ei}")
        monomial = " ".join(monomial_parts) if monomial_parts else "1"

        if coeff == 1:
            terms.append(monomial)
        elif coeff == Rational(1, 2):
            terms.append(f"(1/2){monomial}")
        else:
            terms.append(f"{coeff}·{monomial}")

    return " + ".join(terms) if terms else "0"


if __name__ == "__main__":
    from graded_decomposition import decomposition_for_Gn

    print("=== Validating G_3 decomposition ===")
    D3 = decomposition_for_Gn(3)
    valid, poly, errors = validate_and_extract(D3)
    print(f"Polynomial: {polynomial_to_string(poly, D3.get_var_names())}")
    print()

    print("=== Validating G_4 decomposition ===")
    D4 = decomposition_for_Gn(4)
    valid, poly, errors = validate_and_extract(D4)
    print(f"Polynomial: {polynomial_to_string(poly, D4.get_var_names())}")
    print()

    print("=== Validating G_5 decomposition ===")
    D5 = decomposition_for_Gn(5)
    valid, poly, errors = validate_and_extract(D5)
    print(f"Polynomial: {polynomial_to_string(poly, D5.get_var_names())}")
