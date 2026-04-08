-- Ideal18.m2

load "../Parameters.m2"

generateIdeal18 = (e, n) -> (
    if not (9 == n+7 and n+7 == e) then (
        error "Ideal18 requires 9 = n+7 = e";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1^3, y_2^3 - y_1^2 * y_2 };
    
    I = ideal Glist;
    return {R, I};
)
