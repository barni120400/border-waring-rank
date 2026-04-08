-- Ideal21.m2

load "../Parameters.m2"

generateIdeal21 = (e, n) -> (
    if not (9 == n+6 and n+6 == e) then (
        error "Ideal21 requires 9 = n+6 = e";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1*y_2 + y_3^2, y_1*y_3, y_2^2 - y_3^2 - y_1^3 };
    I = ideal Glist;
    return {R, I};
)
