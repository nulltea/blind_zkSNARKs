from math import sqrt

# Create a header file with proof system parameters for
# proving knowledge of a witness s in Rp^n (Rp = Zp[X]/(X^d + 1))
# such that
#
#   1. s satisfies a linear relation over Rp: As + t = 0
#   2. each element in a partition of s either ..
#      2.1 has binary coefficients only
#      2.2 satisfies an l2-norm bound

vname = "param"               # variable name

ct_count = 2                  # number of ciphertexts
deg   = 2048                  # ring Rp degree d
# mod   = 1099511629429                 # CT modulus (41 bit)
# mod   = 9223372036854776261           # CT modulus (64 bit)
mod   = 582331390981917839              # CT modulus (60 bit)
mod_t = 5 # Plaintext modulus (prime)
dim   = (2,3)                 # dimensions of A in Rp^(m,n)

v_e = 1

# B_e
# sigma_err = 1.55*2 # this set in grandom_static since we used log2o=1
# B_e = sqrt(2*deg)*sigma_err
# print("B_e = ", B_e)

# B_v
# delta_phi_m = deg # to-do check!!?
# inf_norm_t = mod_t
# B_v = mod / (2*delta_phi_m*inf_norm_t) - 0.5
# print("B_v = ", B_v)

#       [ s,    e,  v_inhv]
wpart = [ [0], [1], [2] ]   # partition of s
wl2   = [ 0 , sqrt(2*deg)*(1.55*2), mod / (2*deg*mod_t) - 0.5  ]  # l2-norm bounds: l2(s) <= sqrt(2048)
wbin  = [ 1, 0, 0       ]  # binary coeffs
wrej  = [ 0             ]  # rejection sampling


# Optional: some linf-norm bound on s.
# Tighter bounds result in smaller proofs.
# If not specified, the default is the naive bound max(1,floor(max(wl2))).
# wlinf = 1
