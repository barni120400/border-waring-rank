-- Sort a list of monomials by weight given as a list of weights.

-- monList    : Sequence of monomials in some ring R
-- varList    : List or Sequence of exactly those variables you have weights for
-- weightList : List of integers, same length as varList
sortByWeight = (monList, weightList) -> (
    -- Decorate each monomial with its weighted degree
    decorated = toList apply(monList, m -> (
        exps = flatten exponents m;
        wt = sum apply(length weightList, i -> (
            exps#i * weightList#i
        ));
        (wt, m)
    ));
    
    -- Use coefficient ring from the first monomial's ring
    coeffRing = coefficientRing ring monList#0;
    R = coeffRing[y_1..y_(length weightList), MonomialOrder=>GRevLex];
    decorated = apply(decorated, pair -> (pair#0, substitute(pair#1, R)));
    
    -- -- Sort by weight, breaking ties by the monomial order
    sortedDec = sort(decorated, MonomialOrder=>GRevLex);

    print "Weights:";
    print tex apply(sortedDec, pair -> pair#0);
    

    -- -- Return only the monomials
    return apply(sortedDec, p -> p#1);
);


-- Function extract border forms
-- Input: a Gorenstein Artin Algebra A with degree myVars
-- Output: a list of forms in ZZ/p[x_0..x_(myVars-1)] with border rank myVars
extractBorderForms = method();
extractBorderForms (QuotientRing, List, ZZ) := (A, variableWeights, inputDegree) -> (
    if (dim A != 0) then
        error "expected Artinian Gorenstein algebra";
    if (inputDegree <= 0) then
        error "expected positive degree";
    myDegree = degree A;
    -- construct the extended ring which will give us the border forms
    -- Use the same coefficient ring as the input algebra
    baseCoeffRing = coefficientRing A;
    myBaseRing = baseCoeffRing[x_0..x_(myDegree-1), MonomialOrder=>GRevLex];
    fullRing = myBaseRing ** A;
    -- get basis of the algebra inside of fullRing
    tempVarsA = flatten entries vars A;
    varsA = toList apply(tempVarsA, a -> a_fullRing); 
    basisA = flatten entries basis(Variables => varsA, fullRing);

    -- basisA lives in fullRing, but you only want to weight by the original varsA
    if variableWeights != {} then {
        monList = apply(basisA, m -> substitute(m, baseCoeffRing[y_1..y_(length variableWeights), MonomialOrder=>GRevLex]));
        basisA = sortByWeight(monList, variableWeights);
        basisA = apply(basisA, m -> substitute(m, fullRing));
    } else {
        basisA = sort(basisA, MonomialOrder=>GRevLex);
    };

    print "Generators:";
    print tex basisA;

    print "Socle:";
    print tex annihilator ideal varsA;

    -- get the border forms
    ell = 0_fullRing; 
    apply(0..(myDegree-1), i -> (ell = ell + (basisA_i) * (x_i)_fullRing));
    fullRingFormsList = flatten entries (coefficients(ell^inputDegree, Variables => varsA))_1;
    -- convert the forms to the base ring
    tempBorderFormsList = apply(fullRingFormsList, f -> sub(f, myBaseRing));
    -- filter the forms to keep only those depending on all the variables of the base ring
    borderFormsList = select(tempBorderFormsList, f -> if (basis(1, perpIdeal(f)) == 0) then true else false);

    -- We expect to have only one form at this point
    if (length borderFormsList != 1) then
            error "expected one form";

    return borderFormsList#0;
);


-- Function to compute numerical invariants of a list of polynomials
apolarityInvariants = method();
apolarityInvariants (List) := (formList) -> (
    -- compute the annihilator ideals
    annihilatorList = apply(formList, f -> perpIdeal(f));
    listLength = length annihilatorList;
    -- compute the Hilbert series of the annihilator ideals
    hilbertSeriesList = apply(0..(listLength-1), i -> apply(0..((degree formList_i)_0), j -> hilbertFunction(j, annihilatorList_i)));
    return {formList, annihilatorList, hilbertSeriesList};
);


-- Example to apply the above functions
-- construct a Gorenstein Artin algebra from the Casnati paper, for instance the one given below
-- A = ZZ/p[y_1, y_2]/ideal(y_1*y_2, y_2^2-y_1^3);
-- given the degree of the forms that we want to extract, in this case I picked 5
-- extractBorderForms(ZZ/2311[y_1, y_2]/ideal(y_1*y_2, y_2^2-y_1^3) , 5)