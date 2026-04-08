-- Ideal14.m2

load "../Parameters.m2"

generateIdeal14 = (e, n) -> (
    if not (8 <= n+5 and n+5 == e and e <= 9) then (
        error "Ideal14 requires 8 <= n+5 = e <= 9";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1^2, y_1*y_2, 2*y_1*y_3 + y_2^2, y_3^3, y_2*y_3^2 };
    for j from 4 to n do (
        for i from 1 to j-1 do (
            Glist = append(Glist, y_i * y_j)
        );
        Glist = append(Glist, y_j^2 - y_1*y_3^2)
    );
    I = ideal Glist;
    return {R, I};
)
