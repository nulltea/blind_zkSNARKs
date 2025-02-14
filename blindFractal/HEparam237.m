//// Setting HE 237 parameters

p := 2^64 - 2^32 + 1;
i := 2; // (min,max) = (0,5)
j := 11; // (min,max) = (6,16)
m := (2^j)*3*7;
n := EulerPhi(m);

RNSmodbits := 60;
if vcoed then
    qbits := 375;
else
    qbits := 382;
end if;
qbits_rs := 96;

Zm := IntegerRing(m);
d := Order(Zm ! p);
BFVpack := n/d;

print "BFVpack =", BFVpack;

Zx<x> := PolynomialRing(Integers());
Phi := Zx ! CyclotomicPolynomial(m);
t := x^(7*(2^(i+j-6))) - 2^(2^i);

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
    assert degi eq 2;
    CLPXpack +:= 1;
end for;
assert len_CLPX_ptspace/2 eq CLPXpack;

print "CLPXpack =", CLPXpack;

//// Get specs for HE scheme

load "HEspecs.m";

N0, Nadd, Naut, Nptct, Nctct, Nms := getHEnoise(t);

