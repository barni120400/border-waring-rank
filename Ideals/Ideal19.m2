-- Ideal19.m2

load "../Parameters.m2"

generateIdeal19 = (e, n) -> (
    if not (9 == n+7 and n+7 == e) then (
        error "Ideal19 requires 9 = n+7 = e";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1^3, y_2^3 };
    
    I = ideal Glist;
    return {R, I};
)
