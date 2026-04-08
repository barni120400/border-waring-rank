"""
Generate graded decomposition LaTeX files for ALL algebras.

1. A_{e,1,1} algebras: root-of-unity pattern (direct formula)
2. Connected sum algebras: use the connected_sum algorithm
"""

import logging
import os
from sympy import simplify, N, Abs, Rational
from graded_decomposition import decomposition_for_Gn, GradedDecomposition
from connected_sum import connected_sum
from validate_decomposition import validate_and_extract, polynomial_to_string

logging.basicConfig(level=logging.WARNING)

OUTPUT_DIR = os.path.expanduser("~/Projects/Rafael-Aviv/examples")


def latex_for_Gn(e: int) -> str:
    """Generate LaTeX for A_{e,1,1} = K[y]/(y^e) using root-of-unity pattern."""
    n = e - 1  # ζ_{n} roots
    lines = [r"\textbf{Graded decomposition:}"]
    lines.append(r"\[")
    lines.append(r"\begin{aligned}")
    lines.append(r"&-x_0^d \\")

    # The sum: (1/n) Σ_{k=0}^{n-1} (x_0 + ζ_n^k ε x_1 + ... + ζ_n^{(e-1)k} ε^{e-1} x_{e-1})^d
    terms = []
    for i in range(e):
        if i == 0:
            terms.append("x_0")
        else:
            terms.append(rf"\zeta_{{{n}}}^{{{i}k}}\,\epsilon^{{{i}}} x_{{{i}}}")

    inner = "\n".join([terms[0]] + [f"+ {t}" for t in terms[1:]])
    lines.append(rf"&+ \frac{{1}}{{{n}}}\sum_{{k=0}}^{{{n-1}}}\Bigl(")
    lines.append(inner)
    lines.append(r"\Bigr)^d")

    lines.append(r"\end{aligned}")
    lines.append(r"\]")
    return "\n".join(lines)


def latex_for_connected_sum(decomp: GradedDecomposition) -> str:
    """Generate LaTeX from a GradedDecomposition object."""
    lines = [r"\textbf{Graded decomposition:}"]
    lines.append(r"\[")
    lines.append(r"\begin{aligned}")

    var_names = decomp.get_var_names()
    q = decomp.q

    for j in range(decomp.num_summands):
        # Build the linear form string
        lf_parts = []
        for i in range(decomp.V.cols):
            v = decomp.V[j, i]
            if v == 0:
                continue
            # Variable with epsilon power
            if q[i] == 0:
                var_str = var_names[i]
            else:
                var_str = rf"\varepsilon^{{{q[i]}}}" + " " + var_names[i]

            # Coefficient
            v_simplified = simplify(v)
            # Try numerical simplification for complex roots
            v_num = complex(N(v_simplified))
            if abs(v_num.imag) < 1e-10:
                v_num = v_num.real
                if abs(v_num - round(v_num)) < 1e-10:
                    v_coeff = int(round(v_num))
                elif abs(v_num * 2 - round(v_num * 2)) < 1e-10:
                    v_coeff = Rational(round(v_num * 2), 2)
                else:
                    v_coeff = v_simplified
            else:
                v_coeff = v_simplified

            if v_coeff == 1:
                lf_parts.append(var_str)
            elif v_coeff == -1:
                lf_parts.append(f"-{var_str}")
            else:
                v_str = str(v_coeff)
                # Clean up SymPy output
                if isinstance(v_coeff, (int, Rational)):
                    v_str = str(v_coeff)
                else:
                    v_str = str(simplify(v_coeff))
                lf_parts.append(f"{v_str} {var_str}")

        linear_form = " + ".join(lf_parts).replace("+ -", "- ")

        # Coefficient
        c = decomp.c[j]
        c_num = complex(N(c))
        if abs(c_num.imag) < 1e-10:
            c_num = c_num.real
            if abs(c_num - round(c_num)) < 1e-10:
                c_val = int(round(c_num))
            elif abs(c_num * 2 - round(c_num * 2)) < 1e-10:
                c_val = Rational(round(c_num * 2), 2)
            elif abs(c_num * 3 - round(c_num * 3)) < 1e-10:
                c_val = Rational(round(c_num * 3), 3)
            elif abs(c_num * 4 - round(c_num * 4)) < 1e-10:
                c_val = Rational(round(c_num * 4), 4)
            else:
                c_val = c
        else:
            c_val = c

        # Format the summand
        if j == 0 and linear_form == var_names[0]:
            # Pure x_0^d term
            if c_val == -1:
                summand = f"-{var_names[0]}^d"
            elif c_val == 1:
                summand = f"{var_names[0]}^d"
            else:
                summand = f"{c_val}\\,{var_names[0]}^d"
        else:
            if c_val == 1:
                summand = f"({linear_form})^d"
            elif c_val == -1:
                summand = f"-({linear_form})^d"
            elif isinstance(c_val, Rational) and c_val.q != 1:
                summand = rf"\frac{{{c_val.p}}}{{{c_val.q}}}({linear_form})^d"
            else:
                summand = f"{c_val}\\,({linear_form})^d"

        if j == 0:
            lines.append(f"&{summand} \\\\")
        else:
            sign = "+" if str(summand)[0] != '-' else ""
            lines.append(f"&\\quad {sign} {summand} \\\\")

    # Remove trailing \\  from last line
    lines[-1] = lines[-1].rstrip(" \\\\")

    lines.append(r"\end{aligned}")
    lines.append(r"\]")
    return "\n".join(lines)


def write_decomp(name: str, latex: str):
    """Write a decomp file."""
    path = os.path.join(OUTPUT_DIR, f"decomp_{name}.tex")
    with open(path, 'w') as f:
        f.write(latex + "\n")
    print(f"  Wrote decomp_{name}.tex")


def main():
    # ============================================================
    # Part 1: A_{e,1,1} for e = 3,...,9 (root-of-unity pattern)
    # ============================================================
    print("=== Generating A_{e,1,1} decompositions ===")
    for e in range(3, 10):
        latex = latex_for_Gn(e)
        write_decomp(f"e{e}_n1_i1", latex)

    # ============================================================
    # Part 2: Connected sum decompositions
    # ============================================================
    print("\n=== Generating connected sum decompositions ===")

    # Build decomposition objects for the indecomposable single-generator algebras
    G = {}
    for e in range(3, 10):
        G[e] = decomposition_for_Gn(e)

    # All connected sums we can build from G_n # G_m # ...
    # From the classification, connected sum algebras at each dimension:

    # Dim 4: A_{4,2,2} = G_3 # G_3
    cs_examples = [
        ("e4_n2_i2", [3, 3], "G_3 \\# G_3"),
        # Dim 5
        ("e5_n2_i2", [4, 3], "G_4 \\# G_3"),
        ("e5_n3_i2", [3, 3, 3], "G_3^{\\#3}"),
        # Dim 6
        ("e6_n2_i2", [5, 3], "G_5 \\# G_3"),
        ("e6_n2_i3", [4, 4], "G_4 \\# G_4"),
        ("e6_n3_i2", [4, 3, 3], "G_4 \\# G_3^{\\#2}"),
        ("e6_n4_i2", [3, 3, 3, 3], "G_3^{\\#4}"),
        # Dim 7
        ("e7_n2_i2", [6, 3], "G_6 \\# G_3"),
        ("e7_n2_i3", [5, 4], "G_5 \\# G_4"),
        ("e7_n3_i2", [5, 3, 3], "G_5 \\# G_3^{\\#2}"),
        ("e7_n3_i3", [4, 4, 3], "G_4 \\# G_4 \\# G_3"),
        ("e7_n4_i2", [4, 3, 3, 3], "G_4 \\# G_3^{\\#3}"),
        ("e7_n5_i2", [3, 3, 3, 3, 3], "G_3^{\\#5}"),
        # Dim 8
        ("e8_n2_i2", [7, 3], "G_7 \\# G_3"),
        ("e8_n2_i3", [6, 4], "G_6 \\# G_4"),
        ("e8_n2_i6", [5, 5], "G_5 \\# G_5"),
        ("e8_n3_i2", [6, 3, 3], "G_6 \\# G_3^{\\#2}"),
        ("e8_n3_i3", [5, 4, 3], "G_5 \\# G_4 \\# G_3"),
        ("e8_n3_i5", [7, 3], None),  # Skip — type 5 component, need different input
        ("e8_n3_i12", [4, 4, 4], "G_4^{\\#3}"),
        ("e8_n4_i2", [5, 3, 3, 3], "G_5 \\# G_3^{\\#3}"),
        ("e8_n4_i3", [4, 4, 3, 3], "G_4^{\\#2} \\# G_3^{\\#2}"),
        ("e8_n5_i2", [4, 3, 3, 3, 3], "G_4 \\# G_3^{\\#4}"),
        ("e8_n6_i2", [3, 3, 3, 3, 3, 3], "G_3^{\\#6}"),
    ]

    for name, parts, label in cs_examples:
        if label is None:
            print(f"  SKIPPED {name} (non-G_n component)")
            continue

        print(f"  Computing {name} = {label}...")
        try:
            # Build iterated connected sum
            D = G[parts[0]]
            for p in parts[1:]:
                D = connected_sum(D, G[p])

            # Validate
            valid, poly, errors = validate_and_extract(D, verbose=False)
            if not valid:
                print(f"    WARNING: validation failed for {name}")
                for e, m, w in errors[:2]:
                    print(f"      e={e} m={m} w={w}")

            # Generate LaTeX
            latex = latex_for_connected_sum(D)
            write_decomp(name, latex)
            print(f"    {D.num_summands} summands, valid={valid}")
        except Exception as ex:
            print(f"    ERROR: {ex}")

    # Also handle connected sums with non-G_n indecomposable components
    # A_{7,3,4} = A_{6,2,4} # G_3 — but we don't have A_{6,2,4} decomposition yet
    # A_{8,3,13} = A_{6,2,4} # G_4 — same
    # A_{8,4,4} = A_{6,2,4} # G_3 # G_3 — same
    # These need A_{6,2,4}'s graded decomposition as input, which is handcrafted
    print("\n=== Connected sums with non-G_n components (need handcrafted input) ===")
    print("  e7_n3_i4 = A_{6,2,4} # G_3  — needs A_{6,2,4} decomp")
    print("  e8_n3_i13 = A_{6,2,4} # G_4  — needs A_{6,2,4} decomp")
    print("  e8_n3_i14 = ???  — needs investigation")
    print("  e8_n4_i4 = A_{6,2,4} # G_3^2  — needs A_{6,2,4} decomp")


if __name__ == "__main__":
    main()
