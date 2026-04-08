needsPackage "SpechtModule";
needsPackage "CorrespondenceScrolls";

-- The following function computes the annihilator ideal of a polynomial.
-- It has been taken from the PhD Thesis of Alessandro Oneto titled "Waring-type problems for polynomials: Algebra meets Geometry".

perpIdeal = method();
perpIdeal (RingElement) := F-> (
    if (not isHomogeneous(F)) then
      error "expected homogeneous polynomial";
    P := ring(F);
    D := first degree F;
    L := for j from 0 to D list (
        basis(j,P) * mingens kernel transpose diff(transpose basis(j,P), diff(basis(D-j,P),F))
        );
    return trim ideal L;
);

-- Compute the apolar inner product of two homogeneous polynomials f and g
homogeneousApolarInnerProduct = (f, g) -> (
    R := ring f;
    d := degree f;
    if degree g != d then error "Polynomials must have the same degree";
    Pmatrix := matrix{{f,g}};
    (M, C) := coefficients Pmatrix;
    monomialList := flatten entries M;
    innerProductResult := 0;
    for i from 0 to #monomialList - 1 do (
        Evector := (exponents monomialList_i)_0;
        monDenominator := multinomial Evector;
        innerProductResult = innerProductResult + C_0_i * C_1_i / monDenominator;
    );
    return innerProductResult;
)

-- Test cases for homogeneousApolarInnerProduct
-- R = QQ[x,y, z];
-- f1 = x^2*y + x*y^2 + z^3;
-- g1 = x^2*y + 2*x*y^2 + 3*z^3;
-- result1 = homogeneousApolarInnerProduct(f1, g1);
-- print "Apolar inner product of f1 and g1:";
-- print result1;
-- Expected output: 1/3 + 2/3 + 3 = 4


-- Apolar inner product (not necessarily homogeneous)
-- be careful with this function, because if we extend the Hankel then the first entries are no longer the apolar pairing, since the pairing changes based on degree
-- this is fine, since in BT'20 they only want to extend the Hankel operator of a polynomial, and that may not correspond to the hankel of another polynomial
basicApolarPairing = method();
basicApolarPairing (RingElement, RingElement) := (f, g) -> (
    R := ring f;
    -- we are considering the apolar pairing defined on the space of polynomials up to degree maxDeg
    maxDeg := max((degree f)_0, (degree g)_0);
    Pmatrix := matrix{{f,g}};
    (M, C) := coefficients Pmatrix;
    monomialList := flatten entries M;
    innerProductResult := 0;
    for i from 0 to #monomialList - 1 do (
        Evector := (exponents monomialList_i)_0 | {maxDeg - sum((exponents monomialList_i)_0)};
        monDenominator := multinomial Evector;
        innerProductResult = innerProductResult + C_0_i * C_1_i / monDenominator;
    );
    return innerProductResult;
)

-- Test cases for basicApolarPairing
-- R = QQ[x,y];
-- f1 = x^2*y + x*y + y^2;
-- g1 = x + 2*x*y + 3*y^3;
-- result1 = basicApolarPairing(f1, g1);
-- should be 1/3
-- f2 = x^3 + x^2*y + x + y^3;
-- g2 = 2*x^2*y + y^2 + 3;
-- result2 = basicApolarPairing(f2, g2);
-- should be 2/3

-- This is the code for the general apolar pairing where we specify the degree of the ring under which we are taking the pairing from
apolarPairing = method();
apolarPairing (RingElement, RingElement, ZZ) := (f, g, ringDegree) -> (
    if ringDegree < 0 then error "ringDegree must be non-negative";
    R := ring f;
    if (degree f)_0 > ringDegree then return 0_R;
    if (degree g)_0 > ringDegree then return 0_R;
    
    -- Optimized case: if g is a monomial, compute directly without loop
    if size g == 1 then (
        coeffG := leadCoefficient g;
        monG := leadMonomial g;
        coeffF := coefficient(monG, f);
        Evector := (exponents monG)_0 | {ringDegree - sum((exponents monG)_0)};
        monDenominator := multinomial Evector;
        return coeffF * coeffG / monDenominator;
    );
    
    -- General case: compute using the loop
    -- we are considering the apolar pairing defined on the space of polynomials up to degree ringDegree
    Pmatrix := matrix{{f,g}};
    (M, C) := coefficients Pmatrix;
    monomialList := flatten entries M;
    innerProductResult := 0;
    for i from 0 to #monomialList - 1 do (
        Evector := (exponents monomialList_i)_0 | {ringDegree - sum((exponents monomialList_i)_0)};
        monDenominator := multinomial Evector;
        innerProductResult = innerProductResult + C_0_i * C_1_i / monDenominator;
    );
    return innerProductResult;
)

-- Test cases for apolarPairing
-- R = QQ[x,y];
-- f1 = x^2*y + x*y + y^2;
-- g1 = x + 2*x*y + 3*y^3;
-- result1 = apolarPairing(f1, g1,3);
-- should be 1/3
-- f2 = x^3 + x^2*y + x + y^3;
-- g2 = 2*x^2*y + y^2 + 3;
-- result2 = apolarPairing(f2, g2, 4);
-- should be 1/6

-- Question: can we simplify the above even further by using either diff or contract?
-- see https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/Macaulay2Doc/html/_diff_spand_spcontract.html
-- seems a bit hard to do directly since we need to divide by the multinomial coefficients


-- Compute the Hankel operator for a homogeneous polynomial polyF up to degree operatorDegree
-- want to return the monomial basis as well, so we can interpret the rows and columns
-- The Hankel operator will always be based at the variable x_0 of the ring of polyF
-- Important notes: the ring of polyF must have variables named x_0, x_1, ..., x_n for some n
-- operatorDegree must be at least the degree of polyF, so we can see all non-zero entries of the operator
-- I am leaving the variable x_0 there because in the definition of apolar inner product we need to remember the degree of the polynomial
-- Important observation: we cannot use the method Hankel to construct the Hankel operator here, since the Hankel method constructs a Hankel matrix, and we do not necessarily want that
hankelOperator = (polyF, operatorDegree) -> (
    R := ring polyF;
    degF := (degree polyF)_0;
    if operatorDegree < degF then error "operatorDegree must be at least the degree of polyF";
    monomialList := flatten entries basis(operatorDegree, R);
    -- print monomialList;
    dehomogenizedMonomialList := apply(monomialList, P -> sub(P, {x_0 => 1}));
    -- print dehomogenizedMonomialList;
    -- dehomogenizedHankelEntriesDifferentials := matrix{dehomogenizedMonomialList} * transpose matrix{dehomogenizedMonomialList};
    -- hankelEntries := apply(monomialList, myEntry -> homogeneousApolarInnerProduct(polyF, myEntry));
    outputHankelOperator := matrix apply(dehomogenizedMonomialList, rowMon -> apply(dehomogenizedMonomialList, colMon -> apolarPairing(polyF, rowMon * colMon, operatorDegree)));
    -- print "Hankel operator:";
    -- print outputHankelOperator;
    return (outputHankelOperator, dehomogenizedMonomialList);
)

-- TODO: modify the above to allow for computation of rectangular Hankel operators (i.e., different row and column sizes
-- TODO: modify to allow for the Hankel operator to to show zero entries beyond degree of polyF - perhaps a contraction with respect to x_0 should be enough for every monomial in monomialList?
-- TODO: figure out what is the best monomial ordering for the Hankel operator, as x_0 is always the base variable (local decomposition)

monomialBasisIndices = (basisList, monList) -> (
    return positions(monList, i -> isMember(i, basisList));
)

multiplicationOperatorTranspose = (myVar, hankelOp, basisList, monList) -> (
    R := ring myVar;
    basisIndices := monomialBasisIndices(basisList, monList);
    -- print basisIndices;
    shiftedBasisIndices := monomialBasisIndices(apply(basisList, b -> b * myVar), monList);
    -- print shiftedBasisIndices;
    shiftedBasisOperator := hankelOp_shiftedBasisIndices^basisIndices;
    basisOp := hankelOp_basisIndices^basisIndices;
    -- print "Basis operator:";
    -- print basisOp;
    -- print "Shifted basis operator:";
    -- print shiftedBasisOperator;
    detBasis := det basisOp;
    if detBasis == 0 then (
        error "Basis operator has zero determinant";
    );
    -- Compute MT = shiftedBasisOperator * (basisOp)^(-1)
    -- Use classical adjoint formula: A^(-1) = adj(A) / det(A)
    -- We return (shiftedBasisOperator * adj(basisOp), detBasis) as the pair (numerator, denominator)
    n := numRows basisOp;
    adj := mutableMatrix(ring basisOp, n, n);
    for i from 0 to n-1 do (
        for j from 0 to n-1 do (
            rows := toList (0..j-1) | toList (j+1..n-1);
            cols := toList (0..i-1) | toList (i+1..n-1);
            minorMat := basisOp_cols^rows;
            cofactor := (-1)^(i+j) * det minorMat;
            adj_(i,j) = cofactor;
        );
    );

    numerator := shiftedBasisOperator * matrix adj;
    return (numerator, detBasis);
)

-- Test cases for hankelOperator
-- R = QQ[x_0, x_1];
-- F1 = (x_0 + x_1)^2*(x_0 - x_1);
-- (H1, ML1) = hankelOperator(F1, 3);
-- Code works on first example for sure 
-- TODO: check other examples
-- F2 = (x_0 + x_1)^3*(x_0 - x_1);
-- (H2, ML2) = hankelOperator(F2, 4);
-- F3 = (x_0 + x_1)^3*(x_0 - x_1)^2;
-- (H3, ML3) = hankelOperator(F3, 5);

-- Example from Bernardi-Taufer paper
-- F = x^5 + 32*x^4*y  - 36*x^4*z - 62*x^3*y^2 + 220*x^3*y*z - 154*x^3*z^2 + 172*x^2*y^3 + 1140*x^2*y*z^2 - 556*x^2*z^3 - 157*x*y^4 + 948*x*y^3*z - 2118*x*y^2*z^2 + 2132*x*y*z^3 - 744*x^2*y^2*z - 799*x*z^4 + 64*y^5 - 482*y^4*z + 1448*y^3*z^2 - 2172*y^2*z^3 + 1628*y*z^4 - 488*z^5;
-- Doing it here with x_0, x_1, x_2 instead of x, y, z
-- R = QQ[x_0, x_1, x_2];
-- F = sub(x_0^5 + 32*x_0^4*x_1- 36*x_0^4*x_2 - 62*x_0^3*x_1^2 + 220*x_0^3*x_1*x_2 - 154*x_0^3*x_2^2 + 172*x_0^2*x_1^3 + 1140*x_0^2*x_1*x_2^2 - 556*x_0^2*x_2^3 - 157*x_0*x_1^4 + 948*x_0*x_1^3*x_2 - 2118*x_0*x_1^2*x_2^2 + 2132*x_0*x_1*x_2^3 - 744*x_0^2*x_1^2*x_2 - 799*x_0*x_2^4 + 64*x_1^5 - 482*x_1^4*x_2 + 1448*x_1^3*x_2^2 - 2172*x_1^2*x_2^3 + 1628*x_1*x_2^4 - 488*x_2^5, {x_0 => 1});
-- (HF, MLF) = hankelOperator(F, 5);
-- B = {1_R} | gens R | {x_1^2, x_1 * x_2};
-- B = delete(x_0, B); -- dehomogenized basis
-- MT = apply({x_1, x_2}, myVar -> multiplicationOperatorTranspose (myVar, HF, B, MLF));
-- Combine the multiplication operators to get common eigenvalues and generalized eigenvectors
-- genCombMT = MT_0 + 2*MT_1;


-- Example from Alessandra's slide
-- I am not getting the same result as in the slides, need to check with her
-- R = QQ[x_0, x_1, x_2];
-- Fslides = sub(-2*x_0^7 - 4*x_0^6*x_1 + 92*x_0^6*x_2 + 15*x_0^5*x_1^2 - 675*x_0^5*x_2^2 - 20*x_0^4*x_1^3 + 2700*x_0^4*x_2^3 + 15*x_0^3*x_1^4 - 6075*x_0^3*x_2^4 - 6*x_0^2*x_1^5 + 7290*x_0^2*x_2^5 + x_0*x_1^6 - 3645*x_0*x_2^6, {x_0 => 1});
-- (Hslides, MLslides) = hankelOperator(Fslides, 7);
-- Bslides = delete(x_0, {1_R} | gens R | {x_1^2, x_2^2, x_1^3});
-- MTslides = apply({x_1, x_2}, myVar -> multiplicationOperatorTranspose (myVar, Hslides, Bslides, MLslides));
-- Combine the multiplication operators to get common eigenvalues and generalized eigenvectors
-- genCombMTslides = MTslides_0 + 2*MTslides_1;



-- Casnati cases with schemes of length 4
-- problematic because it has infinitely many decompositions
-- R = QQ[x_0, x_1, x_2, x_3];
-- F42 = sub(3*x_0*x_1^2+3*x_0*x_2^2+3*x_0^2*x_3, {x_0 => 1 + x_1 + x_2 + x_3});
-- (H42, ML42) = hankelOperator(F42, 3);
-- B = {1_R} | gens R;
-- MT1 = multiplicationOperatorTranspose (x_1, H42, B, ML42)
-- MT2 = multiplicationOperatorTranspose (x_2, H42, B, ML42)
-- MT3 = multiplicationOperatorTranspose (x_3, H42, B, ML42)
-- alternatively 
-- MT = apply({x_1, x_2, x_3}, myVar -> multiplicationOperatorTranspose (myVar, H42, B, ML42))
-- Combine the multiplication operators to get common eigenvalues and generalized eigenvectors
-- genCombMT = MT1 + 2*MT2 + 3*MT3;
-- Jordan decomposition of genCombMT
-- J = matrix {
--     {6, 1, 0, 0},
--     {0, 6, 1, 0},
--     {0, 0, 6, 0},
--     {0, 0, 0, 6}
-- };

-- P = matrix {
--     {-5, -6,  1,  0},
--     {-5, -7,  0, -2},
--     {-5, -8,  0,  1},
--     {-5, -6,  0,  0}
-- };
-- F42prime = sub((1/3) * diff(x_0, 3*x_0*x_1^2+3*x_0*x_2^2+3*x_0^2*x_3), {x_0 => 1, x_1 => x_1 + 1, x_2 => x_2 + 1});
-- (H42prime, ML42prime) = hankelOperator(F42prime, 2);
-- Bprime = {1_R} | gens R;
-- MTprime = apply({x_1, x_2, x_3}, myVar -> multiplicationOperatorTranspose (myVar, H42prime, Bprime, ML42prime));
-- Combine the multiplication operators to get common eigenvalues and generalized eigenvectors
-- genCombMTprime = MTprime_0 + 2*MTprime_1 + 3*MTprime_2;
-- Jordan decomposition of genCombMTprime
-- Jprime = matrix {
--     {0, 1, 0, 0},
--     {0, 0, 1, 0},
--     {0, 0, 0, 0},
--     {0, 0, 0, 0}
-- };
-- Pprime = matrix {
--     {5, 0,  0,  0},
--     {0, 1, -3, -2},
--     {0, 2,  0,  1},
--     {0, 0,  1,  0}
-- };


-- Casnati cases with schemes of length 5
-- R = QQ[x_0, x_1, x_2, x_3, x_4];
-- B = delete(x_0, {1_R} | gens R);
-- d51 = 4;
-- F51 = (binomial(d51,1) * x_0^(d51-1) * x_4 + 4*multinomial({d51-2,1,1})* x_0^(d51-2) * x_1 * x_3 + 6*binomial(d51,2)* x_0^(d51-2) * x_2^2 + 12*multinomial({d51-3,2,1})* x_0^(d51-3) * x_1^2 * x_2 + 24*binomial(d51,4)* x_0^(d51-4) * x_1^4)/d51;
-- (H51, ML51) = hankelOperator(sub(F51, {x_0 => 1}), d51);
-- MT51 = apply({x_1, x_2, x_3, x_4}, myVar -> multiplicationOperatorTranspose (myVar, H51, B, ML51));
-- Combine the multiplication operators to get common eigenvalues and generalized eigenvectors
-- genCombMT51 = MT51_0 + 2*MT51_1 + 3*MT51_2 + 4*MT51_3;
-- F52 = diff(x_0, 2*x_0*x_1^3+3*x_0^2*x_2^2+6*x_0^2*x_1*x_3+2*x_0^3*x_4)/6;
-- -- Example without derivative
-- F52 = sub(2*x_0*x_1^3+3*x_0^2*x_2^2+6*x_0^2*x_1*x_3+2*x_0^3*x_4, {x_0 => 1});
-- (H52, ML52) = hankelOperator(F52, 4);
-- MT52 = apply({x_1, x_2, x_3, x_4}, myVar -> multiplicationOperatorTranspose (myVar, H52, B, ML52));
-- Combine the multiplication operators to get common eigenvalues and generalized eigenvectors
-- genCombMT52 = MT52_0 + 2*MT52_1 + 3*MT52_2 + 4*MT52_3;
-- Jordan decomposition of genCombMT52 
-- TODO 
-- something else (this is the first polynomial in the Casnati examples with e=5, n=2, idealNumber=2)
-- F521 = sub(sub(2*x_0*x_1^3+3*x_0^2*x_2^2+6*x_0^2*x_1*x_3+2*x_0^3*x_4, {x_0 => x_0 + x_1 }), {x_0 => 1, x_1 => x_1 +1, x_2 => x_2 +1, x_3 => x_3 +1, x_4 => x_4 +1});
-- inpDeg = 4;
-- F521 = sub(binomial(inpDeg,1)*x_0^(inpDeg-1)*x_4 -2 *binomial(inpDeg,2)*x_0^(inpDeg-2)*x_2^2 + 2*multinomial({1,1,inpDeg-2})*x_0^(inpDeg-2)*x_1*x_3 + binomial(inpDeg,3)*x_0^(inpDeg-3)*x_1^3, {x_0 => x_0 + x_1});
-- (H521, ML521) = hankelOperator(F521, 4);
-- B = {1_R} | gens R;
-- MT521 = apply({x_1, x_2, x_3, x_4}, myVar -> multiplicationOperatorTranspose (myVar, H521, B, ML521));
-- Combine the multiplication operators to get common eigenvalues and generalized eigenvectors
-- genCombMT521 = MT521_0 + 2*MT521_1 + 3*MT521_2 + 4*MT521_3;
-- Jordan decomposition of genCombMT521 
-- multiplicationOperatorTranspose (x_1, H521, B, ML521)
-- TODO: think about this later - there seems to be a problem with the code, since we are finding a Waring decomposition for this polynomial of length 5


-- Function to display matrix with correct epsilon powers
-- This works around Macaulay2's display bug for polynomial ring towers
displayMatrix = method();
displayMatrix (Matrix, ZZ, ZZ) := (M, maxRows, maxCols) -> (
    numRows := min(maxRows, numrows M);
    numCols := min(maxCols, numcols M);
    print("Matrix (showing " | numRows | " x " | numCols | " block):");
    for i from 0 to numRows-1 do (
        rowStr := "| ";
        for j from 0 to numCols-1 do (
            entry := M_(i,j);
            entryStr := toString(entry);
            -- Pad to reasonable width
            while #entryStr < 20 do entryStr = entryStr | " ";
            rowStr = rowStr | entryStr | " ";
        );
        rowStr = rowStr | "|";
        print rowStr;
    );
);