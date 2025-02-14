clear;

vcoed := true; // if false: zkdel

par_factor := 32; // parallelization factor
RS := true;
print_nice := true;
print_info := false;

//// Load HE specifications
load "HEparam237.m";

//// Define Fractal specifications
len_z := 2^20; // circuit size (+ Fractal b)
rho := 1/2; // rate of RS codes
len_L := Integers() ! (len_z/rho); // size of RS domain
avghw := 3; // avg hamming weight of a R1CS matrix row

//// pack: whatever packing size of z is at start of computation
pack := 2^Floor(Log(2,CLPXpack));
print "Actual pack =", pack;
assert BFVpack/CLPXpack in Integers();
actualBFVpack := pack * BFVpack/CLPXpack;
print "Actual BFVpack =", actualBFVpack;

// Load blind Fractal parameters
load "BFparam.m";

load "utils.m";

//// Start estimating

BFVcts := Ceiling(len_z/actualBFVpack);
print "Encryption time (s):", RealField(10) ! BFVcts*Tenc;

comsize := BFVcts*2*n*qbits;
print "full encrypted transcript size in BFV (MB):", RealField(10) ! comsize/(8*10^6);
print "";

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

// unpacking
//
Caut +:= BFVcts*Ceiling(actualBFVpack/pack);
Cptct +:= BFVcts*Ceiling(actualBFVpack/pack);
Nz := N0 + Naut + Nptct;

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Unpacking", ~T, ~qbits, ~MSeffect, ~Nz);

//// Computing ct[Mz] for M = (A,B,C)

C_MV(~Cadd, ~Cptct, ~Caut, len_z : avghw := avghw, rep := 3);
N_Mz := Nz + N_MV(len_z: avghw := avghw);

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing $\\ct[Mz]$", ~T, ~qbits, ~MSeffect, ~N_Mz);
modswitch(~Nz, MSeffect);

//// Computing ct[f_Mz] for M = (_,A,B,C)

// inverse NTTs in base b1 with w1 ciphertexts in row
for i in [1..4] do
    C_NTT(~Cadd, ~Cptct, ~Caut, len_z, b1 : w := w1);
end for;
N_coeffMz := N_Mz + N_NTT(len_z, b1 : w := w1);
N_coeffz := Nz + N_NTT(len_z, b1 : w := w1);

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing $\\ct[\Coeff{f_z}]$/$\\ct[\Coeff{f_{Mz}}]$", ~T, ~qbits, ~MSeffect, ~N_coeffMz);
modswitch(~N_coeffz, MSeffect);

// NTTs in base b2 with w2 ciphertexts in row
assert w2 gt w1;
assert w2/w1 eq len_L/len_z;

for i in [1..4] do
    C_NTT(~Cadd, ~Cptct, ~Caut, len_L, b2 : w := w2, halfzero := true);
end for;
N_fMz := N_coeffMz + N_NTT(len_L, b2 : w := w2);
N_fz := N_coeffz + N_NTT(len_L, b2 : w := w2);

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing $\\ct[\vec{f_z}]$/$\\ct[\vec{f_{Mz}}]$", ~T, ~qbits, ~MSeffect, ~N_fMz);
modswitch(~N_fz, MSeffect);


//// Computing ct[f_sc]

// Computing ct[g]
degf := 2*len_z;

// assumes we have Fractal t(x) in NTT form
// computing ct[t*f_z]
Cptct +:= degf/pack;
N_tfz := N_fz + Nptct; // can not be squashed in, ct[fz] used later

// computing ct[f_sc]
Cptct +:= 3*(len_L/pack);
Cadd +:= 3*(len_L/pack);
N_fMz_sc := N_fMz + Nptct; // can not be squashed in, ct[fMz] used later

// computing coef(f_sc)
C_NTT(~Cadd, ~Cptct, ~Caut, degf, b3 : w := w3);
N_fsc := Max([N_tfz, N_fMz_sc]);
N_fsc +:= N_NTT(degf, b3 : w := w3);

// computing coef(Xg)
assert degf mod len_z eq 0;
assert degf/len_z eq 2;

assert w3 mod 2 eq 0 or w3 eq 1;
if w3 eq 1 then
    Caut +:= degf/pack;
    Cadd +:= len_z/pack;
    Cptct +:= degf/pack;
    N_g := N_fsc + Naut + Nadd + Nptct;
else
    N_g := N_fsc;
end if;
Cadd +:= len_z/pack;
N_g +:= Nadd;

// computing ct[g]
C_NTT(~Cadd, ~Cptct, ~Caut, len_L, b4 : w := w4, halfzero := true);
N_g +:= N_NTT(len_L, b4 : w := w4); // ptct for div by X is squashed into NTT

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing $\\ct[\vec{g}]$", ~T, ~qbits, ~MSeffect, ~N_g);
modswitch(~N_fsc, MSeffect);
modswitch(~N_fz, MSeffect);
modswitch(~N_fMz, MSeffect);

//// Computing ct[f_FRI/2^k]

// computing ct[f_AzXfBz]
Cctct +:= (len_L/pack);
N_fAzXfBz := N_fMz + Nctct;

// computing ct[f_FRI]
Cptct +:= 8*(len_L/pack); // (one is plaintext)
Cadd +:= 8*(len_L/pack); // (one is plaintext)
// ptct can be squashed into previous for some
N_fFRI := Max([N_fAzXfBz, N_g, N_fsc, N_fz + Nptct, N_fMz + Nptct]);

cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
"Computing $\\ct[\vec{\fFRI}]$", ~T, ~qbits, ~MSeffect, ~N_fFRI);

// computing ct[f_FRI/2^k]
deg_FRI := len_z;
depth_FRI := Min(Log(2,deg_FRI), Log(2,len_L/pack)); // stops at size pack

N := N_fFRI;
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

//// Ring switching to committed ciphertexts

// First we need to modswitch for HE security!
assert qbits gt qbits_rs;
MSeffect := qbits - qbits_rs;
qbits := qbits_rs;
modswitch(~N, MSeffect);
pN := N;

FRI_open := 5; // number of Fp2 elts to open when opening top FRI layer (in Fractal)

if RS then
    for i in [0..depth_FRI] do
        tmp := (len_L/pack)/(2^i);
        if i eq 0 then
            tmp *:= FRI_open;
        end if;
        Cptct +:= 4*tmp;
        Caut +:= 12*tmp;
    end for;
    N +:= Nptct + Naut;
    cumulate(Cadd, Cptct, Caut, Cctct, ~pN, ~pCadd, ~pCptct, ~pCaut, ~pCctct,\
    "Ringswitching", ~T, ~qbits, ~MSeffect, ~N);
    print "";
else
    print "WARNING: not taking into account ring switching time";
    print "";
end if;

//// How many ciphertexts to decrypt?

l := 101; // number of repetitions of FRI needed for 100bit sec.

FRI_hits := function(len : batch := 1)
    batch := Integers() ! batch;
    assert len mod batch eq 0;
    Slen := len/batch;
    num_ct_reveals := 0;
    for i in [1..Log(2,Slen)] do
        temp := (Slen/(2^i))*(1-(1-(2^i)/Slen)^l);
        if i eq 1 then temp *:= FRI_open; end if;
        num_ct_reveals +:= 2*temp;
    end for;
    num_ct_reveals +:= 1;
    return num_ct_reveals;
end function;

print "Total amount of Fp2 elements to reveal:", Round(FRI_hits(len_L : batch := 1));
print "Expected total amount of ciphertext to reveal:", Round(FRI_hits(len_L : batch := pack/(2^2)));
print "";
print "Total execution time (h):", RealField(10) ! T/3600;
print "Total parallelized execution time (min):", RealField(10) ! T/(60*par_factor);
print "Final remaining noise (bits):", RealField(10) ! N, "in", qbits, "bits modulus";

sccomsize:= FRI_hits(len_L : batch := pack/(2^2))*2*3072*qbits;
print "server-to-client communication in GBFV (MB):", RealField(10) ! sccomsize/(8*10^6);

