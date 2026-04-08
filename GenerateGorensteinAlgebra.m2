-- GenerateIdeal.m2

load "Ideals/Ideal1.m2";
load "Ideals/Ideal2.m2";
load "Ideals/Ideal3.m2";
load "Ideals/Ideal4.m2";
load "Ideals/Ideal5.m2";
load "Ideals/Ideal6.m2";
load "Ideals/Ideal7.m2";
load "Ideals/Ideal8.m2";
load "Ideals/Ideal9.m2";
load "Ideals/Ideal10.m2";
load "Ideals/Ideal11.m2";
load "Ideals/Ideal12.m2";
load "Ideals/Ideal13.m2";
load "Ideals/Ideal14.m2";
load "Ideals/Ideal15.m2";
load "Ideals/Ideal16.m2";
load "Ideals/Ideal17.m2";
load "Ideals/Ideal18.m2";
load "Ideals/Ideal19.m2";
load "Ideals/Ideal20.m2";
load "Ideals/Ideal21.m2";
load "Ideals/Ideal22.m2";
load "Ideals/Ideal23.m2";
load "Ideals/Ideal24.m2";
load "Ideals/Ideal25.m2";
load "Ideals/Ideal26.m2";

-- computeWeights.m2
-- Given a polynomial ring R and an ideal I,
-- compute positive integer weights w_i so that every *binomial* generator
-- of I is weighted-homogeneous.  Monomial-only gens are ignored.

computeWeights = (R, I) -> (
    -- Substitute epsilon=0 in the ideal to compute weights without the parameter
    coeffRing := coefficientRing R;
    try (
        I = substitute(I, {epsilon => 0});
    );
    
    gensList   := flatten entries gens I;
    varList    := vars R;

    -- Build the matrix M representing the equations required for homogenization
    equationsList = toList flatten apply(gensList, f -> (
        monomialsList = terms f;
        if #monomialsList == 1 then return {};

        toList apply((1..(#monomialsList - 1)), r -> (
            eL := exponents monomialsList#(r - 1);
            eR := exponents monomialsList#(r);
            equations = apply(length flatten entries varList, j -> (flatten eL)#j - (flatten eR)#j)
        ))        
    ));

    equationsList = select(equationsList, e -> e != {});

    if #equationsList == 0 then (
        print "NOTE: No equations for homogenizing ideal";
        return apply(length gensList, i -> 1);
    );

    M = matrix equationsList;

    -- Solve for the nullspace
    K := kernel M;

    if rank K == 0 then error("The ideal cannot being homogenized");

    if rank K == 1 then (
        return flatten entries (gens K)_{0};
    );

    -- Multiple weight assignments: pick a column with all positive entries
    print ("NOTE: Multiple weight assignments (kernel rank " | toString(rank K) | "), picking one");
    G := gens K;
    numCols := numColumns G;
    -- Try each column, also its negation
    for i from 0 to numCols - 1 do (
        col := flatten entries G_{i};
        if all(col, c -> c > 0) then return col;
        negCol := apply(col, c -> -c);
        if all(negCol, c -> c > 0) then return negCol;
    );
    -- Try linear combinations: sum, difference, and small multiples
    for a from -3 to 3 do (
        for b from -3 to 3 do (
            if a == 0 and b == 0 then continue;
            if numCols < 2 and b != 0 then continue;
            combo := apply(numRows G, i -> (
                val := a * (flatten entries G_{0})#i;
                if numCols >= 2 then val = val + b * (flatten entries G_{1})#i;
                val
            ));
            if all(combo, c -> c > 0) then return combo;
        );
    );
    -- Last resort: return first column (may have negatives)
    return flatten entries G_{0};
);

generateGorensteinAlgebra = method(Options => {deformed => false});
generateGorensteinAlgebra (ZZ, ZZ, ZZ) := opts -> (e, n, f) -> (
    -- call the right form-generator, which returns {R,I}
    pair := if     f == 1 then generateIdeal1(e, n)
            else if f == 2 then generateIdeal2(e, n)
            else if f == 3 then generateIdeal3(e, n)
            else if f == 4 then generateIdeal4(e, n)
            else if f == 5 then generateIdeal5(e, n)
            else if f == 6 then generateIdeal6(e, n)
            else if f == 7 then generateIdeal7(e, n)
            else if f == 8 then generateIdeal8(e, n)
            else if f == 9 then generateIdeal9(e, n, null) -- alpha kept symbolic
            else if f == 10 then generateIdeal10(e, n)
            else if f == 11 then generateIdeal11(e, n)
            else if f == 12 then generateIdeal12(e, n)
            else if f == 13 then generateIdeal13(e, n)
            else if f == 14 then generateIdeal14(e, n)
            else if f == 15 then generateIdeal15(e, n)
            else if f == 16 then generateIdeal16(e, n)
            else if f == 17 then generateIdeal17(e, n, null) -- alpha kept symbolic
            else if f == 18 then generateIdeal18(e, n)
            else if f == 19 then generateIdeal19(e, n)
            else if f == 20 then generateIdeal20(e, n)
            else if f == 21 then generateIdeal21(e, n)
            else if f == 22 then generateIdeal22(e, n)
            else if f == 23 then generateIdeal23(e, n)
            else if f == 24 then generateIdeal24(e, n)
            else if f == 25 then generateIdeal25(e, n)
            else if f == 26 then generateIdeal26(e, n)
            else error "Ideal number must be between 1 and 26";

    R := pair#0;       -- ambient ring
    I := pair#1;       -- defining ideal

    -- If deformed is false, substitute epsilon=0 in both ring and ideal
    if not opts.deformed then (
        try (
            -- Substitute epsilon=0 in the ideal
            I = substitute(I, {epsilon => 0});
            varList := gens R;
            -- Ideal9 and Ideal17 use alpha as a symbolic parameter
            newCoeffRing := if f == 9 or f == 17 then frac(ZZ/p[alpha]) else ZZ/p;
            newR := newCoeffRing[varList];
            -- Substitute the ideal into the new ring
            I = substitute(I, newR);
            R = newR;
        );
    );

    print "Ideal:";
    print tex I;

    variableWeights = try computeWeights(R, I) else {};

    return {R/I, variableWeights};        -- return the quotient and variable weights
)
