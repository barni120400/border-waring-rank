-- Ideal2.m2

load "../Parameters.m2"

generateIdeal2 = (e, n) -> (
    if not (4 <= n+2 and n+2 <= e and e <= 9) then (
        error "Ideal2 requires 4 ≤ n+2 ≤ e ≤ 9";
    );
    -- Create coefficient ring with epsilon as a parameter (transcendental extension)
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = flatten for j from 2 to n list (
        (for i from 1 to j-1 list (y_i - epsilon^2)*y_j) | { y_j^2 - y_1^(e-n) - 8*epsilon^6 }
    );
    I = ideal Glist;
    return {R, I};
)
