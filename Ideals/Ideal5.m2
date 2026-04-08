-- Ideal5.m2

load "../Parameters.m2"

generateIdeal5 = (e, n) -> (
    if not (7 <= n+5 and n+5 <= e and e <= 9) then (
        error "Ideal5 requires 7 ≤ n+5 ≤ e ≤ 9 (e=n+4 is isomorphic to type 3)";
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
