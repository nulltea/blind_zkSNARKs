clear;

vcoed := true; // if false: zkdel

par_factor := 32; // parallelization factor
print_nice := false;
print_info := true;

//// Load HE specifications
load "HEparam237_large.m";

N0 := 7.8;
Nms := 3.4;
Nadd := 0.5;
Naut := 0;
Nptct := 13.6;
Nctct := 17.8;

//// Define FRI specs
len_z := 2^15; // pol degree
//len_z := 2^20; // pol degree
rho := 1/2; // rate of RS codes
len_L := Integers() ! (len_z/rho); // size of RS domain

//// pack: whatever packing size of z is at start of computation
pack := 2^Floor(Log(2,CLPXpack));
print "Actual pack =", pack;

if len_z eq 2^15 then
    if pack eq 2^9 then
        b1 := 2^2;
        w1 := 2^0;
        b2 := 2^2;
        w2 := 2^1;
        qbits := 180;
     elif pack eq 2^8 then
        b1 := 2^2;
        w1 := 2^1;
        b2 := 2^2;
        w2 := 2^2;
        qbits := 180;
     end if;
elif len_z eq 2^20 then
    if pack eq 2^9 then
        b1 := 2^2;
        w1 := 2^1;
        b2 := 2^2;
        w2 := 2^2;
        qbits := 360;
     elif pack eq 2^8 then
        b1 := 2^2;
        w1 := 2^0;
        b2 := 2^2;
        w2 := 2^1;
        qbits := 300;
    end if;
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

// NTTs in base b2 with w2 ciphertexts in row
assert w2 gt w1;
assert w2/w1 eq len_L/len_z;

C_NTT(~Cadd, ~Cptct, ~Caut, len_L, b2 : w := w2, halfzero := true);
N := N + N_NTT(len_L, b2 : w := w2);

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing $\\ct[\vec{f_z}]$/$\\ct[\vec{f_{Mz}}]$", ~T, ~qbits, ~MSeffect, ~N);

//// Computing ct[f_FRI/2^k]

// computing ct[f_FRI/2^k]
deg_FRI := len_z;
depth_FRI := Min(Log(2,deg_FRI), Log(2,len_L/pack)); // stops at size pack

// full composition (compute every layer from top)
for i in [1 .. depth_FRI] do
    tmp := (len_L/2^i)/pack;
    Cadd +:= tmp*(2^i-1);
    Cptct +:= tmp*(2^i);
    if i eq depth_FRI then
        N +:= Nptct + i*Nadd;
    end if;
end for;

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing FRI", ~T, ~qbits, ~MSeffect, ~N);

print "Total execution time (s):", RealField(10) ! T;
print "Total parallelized execution time (s):", RealField(10) ! T/par_factor;
print "Final remaining noise (bits):", RealField(10) ! N, "in", qbits, "bits modulus";

Ctotal := Cadd + Cptct + Caut + Cctct;
print "Total number of operations:", Round(Ctotal);

