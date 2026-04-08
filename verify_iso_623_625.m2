-- VERIFICATION: A_{6,2,3} ≅ A_{6,2,5}
-- 
-- Type 3: K[x,y]/(xy, y^3 - x^3)
-- Type 5: K[x,y]/(x^2*y, y^2 - x^2)
--
-- Claim: phi: type5 -> type3 via x -> x+y, y -> -x+y is an isomorphism
-- Claim: psi: type3 -> type5 via x -> (x-y)/2, y -> (x+y)/2 is its inverse

logFile = "/tmp/verify_iso_623_625.txt" << "";
logFile << "VERIFICATION: A_{6,2,3} isomorphic to A_{6,2,5}" << endl;
logFile << "================================================" << endl << endl;

R = QQ[x, y];

-- Build the two algebras
i3 = ideal(x*y, y^3 - x^3);
i5 = ideal(x^2*y, y^2 - x^2);
t3 = R / i3;
use R;
t5 = R / i5;

logFile << "Type 3: QQ[x,y]/(xy, y^3-x^3)" << endl;
logFile << "Type 5: QQ[x,y]/(x^2*y, y^2-x^2)" << endl;
logFile << "dim(Type 3) = " << degree t3 << endl;
logFile << "dim(Type 5) = " << degree t5 << endl << endl;

-- CHECK 1: phi is well-defined
logFile << "CHECK 1: phi: Type5 -> Type3 via x->x+y, y->-x+y" << endl;
phi = map(t3, t5, {sub(x,t3)+sub(y,t3), -sub(x,t3)+sub(y,t3)});
wd1 = isWellDefined phi;
logFile << "  isWellDefined(phi) = " << toString wd1 << endl;

-- CHECK 2: phi has zero kernel
logFile << "CHECK 2: kernel(phi) = 0?" << endl;
ker1 = kernel phi;
zeroKer1 = (ker1 == ideal(0_t5));
logFile << "  kernel(phi) == 0: " << toString zeroKer1 << endl;

-- CHECK 3: psi is well-defined
logFile << "CHECK 3: psi: Type3 -> Type5 via x->(x-y)/2, y->(x+y)/2" << endl;
use R; t3b = R/i3; use R; t5b = R/i5;
psi = map(t5b, t3b, {(1/2)*(sub(x,t5b)-sub(y,t5b)), (1/2)*(sub(x,t5b)+sub(y,t5b))});
wd2 = isWellDefined psi;
logFile << "  isWellDefined(psi) = " << toString wd2 << endl;

-- CHECK 4: psi has zero kernel  
logFile << "CHECK 4: kernel(psi) = 0?" << endl;
ker2 = kernel psi;
zeroKer2 = (ker2 == ideal(0_t3b));
logFile << "  kernel(psi) == 0: " << toString zeroKer2 << endl;

-- CHECK 5: phi o psi = id on Type3
logFile << "CHECK 5: phi(psi(x)) = x and phi(psi(y)) = y in Type3?" << endl;
use R; t3c = R/i3; use R; t5c = R/i5;
phi2 = map(t3c, t5c, {sub(x,t3c)+sub(y,t3c), -sub(x,t3c)+sub(y,t3c)});
psi2 = map(t5c, t3c, {(1/2)*(sub(x,t5c)-sub(y,t5c)), (1/2)*(sub(x,t5c)+sub(y,t5c))});
comp1 = phi2 * psi2;  -- composition: Type3 -> Type5 -> Type3
checkX = comp1(sub(x,t3c)) == sub(x,t3c);
checkY = comp1(sub(y,t3c)) == sub(y,t3c);
logFile << "  phi(psi(x)) == x: " << toString checkX << endl;
logFile << "  phi(psi(y)) == y: " << toString checkY << endl;

-- CHECK 6: psi o phi = id on Type5
logFile << "CHECK 6: psi(phi(x)) = x and psi(phi(y)) = y in Type5?" << endl;
comp2 = psi2 * phi2;  -- composition: Type5 -> Type3 -> Type5
checkX2 = comp2(sub(x,t5c)) == sub(x,t5c);
checkY2 = comp2(sub(y,t5c)) == sub(y,t5c);
logFile << "  psi(phi(x)) == x: " << toString checkX2 << endl;
logFile << "  psi(phi(y)) == y: " << toString checkY2 << endl;

-- SUMMARY
logFile << endl << "SUMMARY" << endl;
logFile << "=======" << endl;
allPass = wd1 and zeroKer1 and wd2 and zeroKer2 and checkX and checkY and checkX2 and checkY2;
logFile << "All 8 checks pass: " << toString allPass << endl;
if allPass then logFile << "CONCLUSION: A_{6,2,3} and A_{6,2,5} are isomorphic over QQ." << endl
else logFile << "CONCLUSION: VERIFICATION FAILED." << endl;

logFile << close;
