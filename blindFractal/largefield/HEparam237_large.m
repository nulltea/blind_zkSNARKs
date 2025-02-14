//// Setting HE 237 parameters

//m := 2^15;
m := 2^14;
n := EulerPhi(m);

RNSmodbits := 60;

Zx<x> := PolynomialRing(Integers());
Phi := Zx ! CyclotomicPolynomial(m);

// packs 512 elements of Fp for 269bit p
b := 332;
p := b^(2^5) + 1;
//t := x^(2^9) - b;
t := x^(2^8) - b;

R := Resultant(Phi, t); // number of elements in pt space
len_CLPX_ptspace := Valuation(R,p); // number of elements of size p in pt space

// now we check how these are divided over the slots
Fp := GF(p);
Fpz<z> := PolynomialRing(Fp);
T := GCD(Fpz ! t, Fpz ! Phi); // Bootland trick
facs := Factorisation(T);
CLPXpack := 0;
for fi in facs do
    degi := Degree(fi[1]);
    //print degi;
    CLPXpack +:= 1;
end for;

print "CLPXpack =", CLPXpack;

//// Get specs for HE scheme

load "../HEspecs.m";

