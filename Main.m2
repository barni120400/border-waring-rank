load "HelperFunctions.m2"
load "GenerateGorensteinAlgebra.m2"
load "ExtractBorderForms.m2"

-- Compute the ungraded socle degree of an Artinian algebra A.
-- This is the largest k such that m^k != 0, where m is the maximal ideal.
computeSocleDegree = A -> (
    m := ideal vars A;
    k := 1;
    while ideal(0_A) != m^k do k = k + 1;
    return k - 1;
);

-- Parameters (edit these)
idealDegree = 5;
numVars = 1;
idealNumber = 1;

pair = generateGorensteinAlgebra(idealDegree, numVars, idealNumber, deformed=>false);
algebra = pair#0;
variableWeights = pair#1;
formDegree = computeSocleDegree(algebra);

print("Algebra dimension: " | toString(degree algebra));
print("Socle degree: " | toString(formDegree));
print("Variable weights: " | toString(variableWeights));

form = extractBorderForms(algebra, variableWeights, formDegree);
print("Form: ");
print tex form;

apolarIdeal = perpIdeal(form);
print("Apolar ideal: ");
print tex apolarIdeal;
