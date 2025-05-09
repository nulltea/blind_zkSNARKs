// auto-generated by lnp-quad-eval-codegen.sage from ../tests/lnp-quad-eval-params1.sage.
//
// protocol is statistically complete with correctness error >= 1 - 2^(-6)
// protocol is simulatable under MLWE(26,23,[-1,1])
// protocol is knowledge-sound with knowledge error <= 2^(-127.0)
//
// Ring
// degree d = 64
// modulus q = 1099511627917, log(q) ~ 40.0
// factors q = q1
//
// Compression
// D = 8
// gamma = 130404, log(gamma) ~ 16.992629
// m = (q-1)/gamma = 8431579, log(m) ~ 23.007371
//
// Dimensions of secrets
// s1: m1 = 10
// m: l = 3
// s2: m2 = 49
//
// Size of secrets
// l2(s1) <= alpha = 75.894664
// m unbounded
// s2 uniform in [-nu,nu] = [-1,1]
//
// Challenge space
// c uniform in [-omega,omega] = [-8,8], o(c)=c, sqrt(l1(o(c)*c)) <= eta = 140
//
// Standard deviations
// stdev1 = 101580.8, log(stdev1/1.55) = 16.0
// stdev2 = 50790.4, log(stdev2/1.55) = 15.0
//
// Repetition rate
// M1 = 4.0724839
// M2 = 7.9736184
// total = 32.472433
//
// Security
// MSIS dimension: 17
// MSIS root hermite factor: 1.004346
// MLWE dimension: 26
// MLWE root hermite factor: 1.0042737
//
// 23 bit moduli for degree 64: [8386817, 8386177, 8385281, 8384641, 8383489]
// bit length of products: [22, 45, 68, 91, 114]
// inverses: [1, 1664132, -20283, -3834820, -819883]

#include "lazer.h"
static const limb_t params1_q_limbs[] = {1099511627917UL};
static const int_t params1_q = {{(limb_t *)params1_q_limbs, 1, 0}};
static const limb_t params1_qminus1_limbs[] = {1099511627916UL};
static const int_t params1_qminus1 = {{(limb_t *)params1_qminus1_limbs, 1, 0}};
static const limb_t params1_m_limbs[] = {8431579UL};
static const int_t params1_m = {{(limb_t *)params1_m_limbs, 1, 0}};
static const limb_t params1_mby2_limbs[] = {0};
static const int_t params1_mby2 = {{(limb_t *)params1_mby2_limbs, 1, 0}};
static const limb_t params1_gamma_limbs[] = {130404UL};
static const int_t params1_gamma = {{(limb_t *)params1_gamma_limbs, 1, 0}};
static const limb_t params1_gammaby2_limbs[] = {65202UL};
static const int_t params1_gammaby2 = {{(limb_t *)params1_gammaby2_limbs, 1, 0}};
static const limb_t params1_pow2D_limbs[] = {256UL};
static const int_t params1_pow2D = {{(limb_t *)params1_pow2D_limbs, 1, 0}};
static const limb_t params1_pow2Dby2_limbs[] = {128UL};
static const int_t params1_pow2Dby2 = {{(limb_t *)params1_pow2Dby2_limbs, 1, 0}};
static const limb_t params1_Bsq_limbs[] = {45753870618406UL, 0UL};
static const int_t params1_Bsq = {{(limb_t *)params1_Bsq_limbs, 2, 0}};
static const limb_t params1_scM1_limbs[] = {13768581241400741304UL, 1337092322823884803UL, 4UL};
static const int_t params1_scM1 = {{(limb_t *)params1_scM1_limbs, 3, 0}};
static const limb_t params1_scM2_limbs[] = {9341264381171379285UL, 17960088946314419876UL, 7UL};
static const int_t params1_scM2 = {{(limb_t *)params1_scM2_limbs, 3, 0}};
static const limb_t params1_stdev1sq_limbs[] = {10318658929UL, 0UL};
static const int_t params1_stdev1sq = {{(limb_t *)params1_stdev1sq_limbs, 2, 0}};
static const limb_t params1_stdev2sq_limbs[] = {2579664732UL, 0UL};
static const int_t params1_stdev2sq = {{(limb_t *)params1_stdev2sq_limbs, 2, 0}};
static const limb_t params1_inv2_limbs[] = {549755813958UL};
static const int_t params1_inv2 = {{(limb_t *)params1_inv2_limbs, 1, 1}};
static const polyring_t params1_ring = {{params1_q, 64, 41, 6, moduli_d64, 5, params1_inv2}};
static const dcompress_params_t params1_dcomp = {{ params1_q, params1_qminus1, params1_m, params1_mby2, params1_gamma, params1_gammaby2, params1_pow2D, params1_pow2Dby2, 8, 1, 24 }};
static const abdlop_params_t params1_quad_eval = {{ params1_ring, params1_dcomp, 10, 49, 3, 3, 17, params1_Bsq, 1, 8, 5, 140, 1, 16, params1_scM1, params1_stdev1sq, 1, 15, params1_scM2, params1_stdev2sq}};
static const abdlop_params_t params1_quad_many = {{ params1_ring, params1_dcomp, 10, 49, 5, 1, 17, params1_Bsq, 1, 8, 5, 140, 1, 16, params1_scM1, params1_stdev1sq, 1, 15, params1_scM2, params1_stdev2sq}};
static const lnp_quad_eval_params_t params1 = {{ params1_quad_eval, params1_quad_many, 4}};
