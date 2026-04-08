-- Ideal8.m2

load "../Parameters.m2"

generateIdeal8 = (e, n) -> (
    if not (8 <= n+6 and n+6 == e and e <= 9) then (
        error "Ideal8 requires 8 ≤ n+6 = e ≤ 9";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1^3*y_2 - y_1^4, y_2^2 - y_1^3, y_1^5 };
    for j from 3 to n do (
        for i from 1 to j-1 do Glist = append(Glist, y_i*y_j);
        Glist = append(Glist, y_j^2 - y_1^4)
    );
    I = ideal Glist;
    return {R, I};
)
