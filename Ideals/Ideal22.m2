-- Ideal22.m2

load "../Parameters.m2"

generateIdeal22 = (e, n) -> (
    if not (9 == n+6 and n+6 == e) then (
        error "Ideal22 requires 9 = n+6 = e";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1*y_2, y_1*y_3, y_2*y_3,  y_2^3- y_1^4, y_3^3 - y_1^4};
    I = ideal Glist;
    return {R, I};
)
