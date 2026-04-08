-- Ideal5.m2

load "../Parameters.m2"

generateIdeal5 = (e, n) -> (
    if not (6 <= n+4 and n+4 <= e and e <= 9) then (
        error "Ideal5 requires 6 ≤ n+4 ≤ e ≤ 9";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    Glist = { y_1^2*y_2, y_2^2 - y_1^(e-n-2) };
    for j from 3 to n do (
        for i from 1 to j-1 do Glist = append(Glist, y_i*y_j);
        Glist = append(Glist, y_j^2 - y_1^(e-n-1))
    );
    I = ideal Glist;
    return {R, I};
)
