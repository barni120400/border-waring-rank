-- Ideal15.m2

load "../Parameters.m2"

generateIdeal15 = (e, n) -> (
    if not (9 == n+7 and n+7 == e) then (
        error "Ideal15 requires 9 = n+7 = e";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1 * y_2, y_2^4 - y_1^5 };
    
    I = ideal Glist;
    return {R, I};
)
