name = "params1"                                # param variable name

# log2q should be greater than log2(fhe_ctxt_mod) + 2 + log2((d+3)/2) 
# log2((d+3)/2) is 6(if d=64) or 7(if d=128)
log2q = 67                                       # ring modulus bits
d = 64                                           # ring degree

m1 = 64                                          # length of bounded message s1
#m1 = 10
alpha = sqrt(d*m1 * 3 ^ 2)                       # bound on s1: l2(s1) <= alpha
l = 32                                           # length of unbounded message m
#l = 2

nbin = 0                                         # set length of vector with binary coefficients

n = []                                 # set lengths of vectors bounded in l2 norm
B = [1]    # set l2 norm bounds
#n = [2, 1]                  # set lengths of vectors bounded in l2 norm
#B = [sqrt(128), sqrt(64)]   # set l2 norm bounds


nprime = 1                                       # set length of vector bounded in linf norm
Bprime = 4                                       # set linf norm bound

# gamma1 = 14
# gamma2 = 1
# gamma3 = 5
# gamma4 = 5
