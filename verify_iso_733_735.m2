-- VERIFICATION: A_{7,3,3} ≅ A_{7,3,5} over an algebraically closed field
--
-- Type 3 at (7,3): K[x,y,z]/(xy, y^3-x^3, xz, yz, z^2-x^3)
-- Type 5 at (7,3): K[x,y,z]/(x^2*y, y^2-x^2, xz, yz, z^2-x^3)
--
-- The isomorphism requires sqrt(2) for the z-scaling:
--   phi: x -> x+y, y -> -x+y, z -> sqrt(2)*z
--   psi: x -> (x-y)/2, y -> (x+y)/2, z -> z/sqrt(2)
-- Works over any algebraically closed field of char != 2.

logFile = "verify_iso_733_735_output.txt" << "";
logFile << "VERIFICATION: A_{7,3,3} isomorphic to A_{7,3,5}" << endl;
logFile << "================================================" << endl;
logFile << "Working over QQ(sqrt(2)) since the isomorphism requires sqrt(2)." << endl << endl;

K = toField(QQ[s]/(s^2-2));
R = K[x, y, z];

i3 = ideal(x*y, y^3-x^3, x*z, y*z, z^2-x^3);
use R;
i5 = ideal(x^2*y, y^2-x^2, x*z, y*z, z^2-x^3);
t3 = R / i3;
use R;
t5 = R / i5;

sval = sub(s, K);
sinv = sub(s/2, K);  -- 1/sqrt(2) = sqrt(2)/2

logFile << "Type 3: K[x,y,z]/(xy, y^3-x^3, xz, yz, z^2-x^3)" << endl;
logFile << "Type 5: K[x,y,z]/(x^2*y, y^2-x^2, xz, yz, z^2-x^3)" << endl;
logFile << "K = QQ(sqrt(2)), sqrt(2) = s" << endl;
logFile << "dim(Type 3) = " << degree t3 << endl;
logFile << "dim(Type 5) = " << degree t5 << endl << endl;

-- CHECK 1: phi well-defined
logFile << "CHECK 1: phi: Type5 -> Type3 via x->x+y, y->-x+y, z->s*z" << endl;
phi = map(t3, t5, {sub(x,t3)+sub(y,t3), -sub(x,t3)+sub(y,t3), sval*sub(z,t3)});
wd1 = isWellDefined phi;
logFile << "  isWellDefined(phi) = " << toString wd1 << endl;

-- CHECK 2: kernel(phi) = 0
logFile << "CHECK 2: kernel(phi) = 0?" << endl;
ker1 = kernel phi;
zeroKer1 = (ker1 == ideal(0_t5));
logFile << "  kernel(phi) == 0: " << toString zeroKer1 << endl;

-- CHECK 3: psi well-defined
logFile << "CHECK 3: psi: Type3 -> Type5 via x->(x-y)/2, y->(x+y)/2, z->(s/2)*z" << endl;
use R; t3b = R/i3; use R; t5b = R/i5;
psi = map(t5b, t3b, {(1/2)*(sub(x,t5b)-sub(y,t5b)), (1/2)*(sub(x,t5b)+sub(y,t5b)), sinv*sub(z,t5b)});
wd2 = isWellDefined psi;
logFile << "  isWellDefined(psi) = " << toString wd2 << endl;

-- CHECK 4: kernel(psi) = 0
logFile << "CHECK 4: kernel(psi) = 0?" << endl;
ker2 = kernel psi;
zeroKer2 = (ker2 == ideal(0_t3b));
logFile << "  kernel(psi) == 0: " << toString zeroKer2 << endl;

-- CHECK 5: phi o psi = id
logFile << "CHECK 5: phi o psi = id on Type3?" << endl;
use R; t3c = R/i3; use R; t5c = R/i5;
phi2 = map(t3c, t5c, {sub(x,t3c)+sub(y,t3c), -sub(x,t3c)+sub(y,t3c), sval*sub(z,t3c)});
psi2 = map(t5c, t3c, {(1/2)*(sub(x,t5c)-sub(y,t5c)), (1/2)*(sub(x,t5c)+sub(y,t5c)), sinv*sub(z,t5c)});
comp1 = phi2 * psi2;
chkX = comp1(sub(x,t3c)) == sub(x,t3c);
chkY = comp1(sub(y,t3c)) == sub(y,t3c);
chkZ = comp1(sub(z,t3c)) == sub(z,t3c);
logFile << "  phi(psi(x)) == x: " << toString chkX << endl;
logFile << "  phi(psi(y)) == y: " << toString chkY << endl;
logFile << "  phi(psi(z)) == z: " << toString chkZ << endl;

-- CHECK 6: psi o phi = id
logFile << "CHECK 6: psi o phi = id on Type5?" << endl;
comp2 = psi2 * phi2;
chkX2 = comp2(sub(x,t5c)) == sub(x,t5c);
chkY2 = comp2(sub(y,t5c)) == sub(y,t5c);
chkZ2 = comp2(sub(z,t5c)) == sub(z,t5c);
logFile << "  psi(phi(x)) == x: " << toString chkX2 << endl;
logFile << "  psi(phi(y)) == y: " << toString chkY2 << endl;
logFile << "  psi(phi(z)) == z: " << toString chkZ2 << endl;

-- SUMMARY
logFile << endl << "SUMMARY" << endl << "=======" << endl;
allPass = wd1 and zeroKer1 and wd2 and zeroKer2 and chkX and chkY and chkZ and chkX2 and chkY2 and chkZ2;
logFile << "All 10 checks pass: " << toString allPass << endl;
if allPass then logFile << "CONCLUSION: A_{7,3,3} and A_{7,3,5} are isomorphic over any algebraically closed field of char != 2." << endl
else logFile << "CONCLUSION: VERIFICATION FAILED." << endl;

logFile << close;
