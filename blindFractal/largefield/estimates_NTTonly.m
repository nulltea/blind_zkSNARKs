clear;

vcoed := true; // if false: zkdel

par_factor := 32; // parallelization factor
print_nice := false;
print_info := true;

//// Load HE specifications
load "HEparam237_large.m";

qbits := 120;

N0 := 7.8;
Nms := 3.4;
Nadd := 0.5;
Naut := 0;
Nptct := 13.6;
Nctct := 17.8;

//// Define FRI specs
len_z := 2^18; // pol degree
rho := 1/2; // rate of RS codes
len_L := Integers() ! (len_z/rho); // size of RS domain

//// pack: whatever packing size of z is at start of computation
pack := 2^Floor(Log(2,CLPXpack));
print "Actual pack =", pack;

if pack eq 2^9 then
    b1 := 2^3;
    w1 := 2^0;
elif pack eq 2^8 then
    b1 := 2^2;
    w1 := 2^0;
end if;

load "../utils.m";

//// Start estimating

// set counters to zero
Cadd := 0;
Cptct := 0;
Caut := 0;
Cctct := 0;
T := 0;
MSeffect := 0;

////
// we assume transcript arrives in normal in ROW-MAJOR packing order
////

//// Computing ct[f_Mz] for M = (_,A,B,C)

// inverse NTTs in base b1 with w1 ciphertexts in row
C_NTT(~Cadd, ~Cptct, ~Caut, len_z, b1 : w := w1);
N := N0 + N_NTT(len_z, b1 : w := w1);

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing $\\ct[\Coeff{f_z}]$/$\\ct[\Coeff{f_{Mz}}]$", ~T, ~qbits, ~MSeffect, ~N);


print "Total execution time (s):", RealField(10) ! T;
print "Total parallelized execution time (s):", RealField(10) ! T/par_factor;
print "Final remaining noise (bits):", RealField(10) ! N, "in", qbits, "bits modulus";

Ctotal := Cadd + Cptct + Caut + Cctct;
print "Total number of operations:", Round(Ctotal);

