-- Ideal9.m2

load "../Parameters.m2"

generateIdeal9 = (e, n, alpha) -> (
    if not (8 <= n+5 and n + 5 == e and e <= 9) then (
        error "Ideal9 requires 8 <= n+5 = e <= 9";
    );
    K = frac(ZZ/p[epsilon]);
    R = K[y_1..y_n];
    -- Initial generators for the 10th form
    Glist = { y_1*y_2 + y_3^2, y_1*y_3, y_1^2 + y_2^2 - alpha*y_3^2 };
    -- If n > 3, add cross terms and cubics for higher variables
    for j from 4 to n do (
        for i from 1 to j-1 do (
            Glist = append(Glist, y_i*y_j)
        );
        Glist = append(Glist, y_j^2 - y_1^3)
    );
    I = ideal Glist;
    return {R, I};
)
