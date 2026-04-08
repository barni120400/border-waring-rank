-- Ideal1.m2

load "../Parameters.m2"

generateIdeal1 = (e, n) -> (
    if not (n == 1 and 1 <= e and e <= 9) then (
        error "Ideal1 requires n = 1 and 1 ≤ e ≤ 9";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1];
    I = ideal(y_1^e);
    return {R, I};
)
