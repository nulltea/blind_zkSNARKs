#include <stdio.h>
#include "lazer.h"
#include "../src/brandom.h"
#include "../src/memory.h"
// #include "brandom.h"
// #include "vdec_params_tbox.h"
#include "vdec_params.h"
//#include "lnp-quad-eval-params1.h"
// #include "vdec_ct_newq.h"
//#include "vdec_ct.h"
// #include "vdec_ct_62bits.h"
//#include "vdec_ct_60bits.h"
#include "vdec_ct_gbfv_60bits.h"
#include <mpfr.h>
#include <sys/time.h>

#define N 1 /* number of quadratic equations */
#define M 1 /* number of quadratic eval equations */
#define CT_COUNT 64 /* number of ciphertexts */
#define GBFV 1 /* if using BFV instead, set to 0 (changes rotation function) */
#define DEGREE 12288 /* fhe degree */
//#define DEGREE 2048

/* Number of elements in an n x n (upper) diagonal matrix. */
#define NELEMS_DIAG(n) (((n) * (n) - (n)) / 2 + (n))

static void vdec_lnp_tbox (uint8_t seed[32], const lnp_quad_eval_params_t params, 
                           polyvec_t sk, int8_t sk_sign[], polyvec_t ct0, polyvec_t ct1, 
                           polyvec_t m_delta, unsigned int fhe_degree);

static inline void _expand_R_i2 (int8_t *Ri, unsigned int ncols, unsigned int i,
                                const uint8_t cseed[32]);

static void __shuffleauto2x2submatssparse (spolymat_t a);
static void print_exec_time(struct timeval start, struct timeval end, const char *str);
static void __shuffleautovecsparse (spolyvec_t r);
static void __schwartz_zippel_accumulate (
    spolymat_ptr R2i, spolyvec_ptr r1i, poly_ptr r0i, spolymat_ptr Rprime2i[],
    spolyvec_ptr rprime1i[], poly_ptr rprime0i[], unsigned int M_alt,
    const intvec_t v, const lnp_quad_eval_params_t params);
static void __schwartz_zippel_auto (spolymat_ptr R2i, spolyvec_ptr r1i,
                                    poly_ptr r0i, spolymat_ptr R2i2,
                                    spolyvec_ptr r1i2, poly_ptr r0i2,
                                    const lnp_quad_eval_params_t params);
static void __schwartz_zippel_accumulate2 (
    spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[],
    spolymat_ptr R2i2[], spolyvec_ptr r1i2[], poly_ptr r0i2[],
    spolymat_ptr R2primei[], spolyvec_ptr r1primei[], poly_ptr r0primei[],
    unsigned int M_alt, const uint8_t seed[32], uint32_t dom,
    const lnp_quad_eval_params_t params);
static void __schwartz_zippel_accumulate_beta (
    spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[],
    spolymat_ptr R2i2[], spolyvec_ptr r1i2[], poly_ptr r0i2[],
    spolymat_ptr R2t, spolyvec_ptr r1t, poly_ptr r0t, const uint8_t seed[32],
    uint32_t dom, const lnp_quad_eval_params_t params, const unsigned int nbounds);
static void __schwartz_zippel_accumulate_z (
    spolymat_ptr R2i[], spolyvec_ptr r1i[], poly_ptr r0i[], 
    spolymat_ptr R2i2[], spolyvec_ptr r1i2[], poly_ptr r0i2[],
    spolymat_ptr R2t, spolyvec_ptr r1t, poly_ptr r0t, 
    intvec_t u_, polyvec_t z4, intvec_t ct1, const uint8_t seed[32],
    uint32_t dom, const lnp_quad_eval_params_t params,
    const unsigned int nprime);

static inline void
__evaleq (poly_ptr res, spolymat_ptr Rprime2, spolyvec_ptr rprime1,
          poly_ptr rprime0, polyvec_ptr s)
{
  polyring_srcptr Rq = rprime0->ring;
  polyvec_t tmp;

  ASSERT_ERR (res != rprime0);

  polyvec_alloc (tmp, Rq, spolymat_get_nrows (Rprime2));

  if (rprime0 != NULL)
    poly_set (res, rprime0);
  else
    poly_set_zero (res);

  if (rprime1 != NULL)
    poly_adddot2 (res, rprime1, s, 0);

  if (Rprime2 != NULL)
    {
      polyvec_mulsparse (tmp, Rprime2, s);
      polyvec_fromcrt (tmp);
      poly_adddot (res, s, tmp, 0);
    }
  poly_fromcrt (res);
  poly_mod (res, res);
  poly_redc (res, res);

  polyvec_free (tmp);
}




int main(void)
{
    lazer_init();

    /* importing ciphertext materials */
    /* init Rq */
    abdlop_params_srcptr abdlop = params1->quad_eval;
    polyring_srcptr Rq = abdlop->ring;
    const unsigned int proof_degree = Rq->d;
    // printf("\nring modulus bits = %d\n", Rq->d);
    
    /* fhe parameters */
    const unsigned int fhe_degree = DEGREE;
    //const int_t fhe_modulus;
    //int_alloc(fhe_modulus, Rq->q->nlimbs);
    //int_set_i64(fhe_modulus, 2^54+1);
    printf("\nproof modulus: \n");

    /* INIT sk */
    POLYVEC_T(sk_vec_polys, Rq, fhe_degree/proof_degree);
    printf("fhe_degree: %d\n", fhe_degree);
    printf("proof_degree: %d\n", proof_degree);
    printf("fhe_degree/proof_degree: %d\n", fhe_degree/proof_degree);
    poly_ptr poly;
    intvec_ptr coeffs;
    for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
        poly = polyvec_get_elem(sk_vec_polys, i);
        coeffs = poly_get_coeffvec (poly);
        for (size_t j=0; j<proof_degree ; j++) {
            intvec_set_elem_i64(coeffs,j,static_sk[j+i*proof_degree]);
            //printf("sk elem (%d): %d\n", j+i*proof_degree, static_sk[j+i*proof_degree]);
        }
    }
    //print_polyvec_element("First element of sk", sk_vec_polys, 0, 64);

    /* INIT ct0 */
    // in this section, we flatten the ct0[][] 2D array to a large 1D array of polynomials
    // Calculate total polynomials needed
    size_t polys_per_ct = fhe_degree / proof_degree;
    size_t total_polys = CT_COUNT * polys_per_ct;
    polyvec_t ct0_vec_polys;
    polyvec_alloc(ct0_vec_polys, Rq, total_polys);
    // POLYVEC_T(ct0_vec_polys, Rq, total_polys);
    
    // You can treat static_ct0 as a one-dimensional array in memory 
    // since C stores arrays in row-major order. Cast it to a one-dimensional array as follows:
    const int64_t *flat_static_ct0 = (const int64_t *)static_ct0;

    size_t total_elements = CT_COUNT * fhe_degree;
    size_t element_index = 0;

    // flatten
    // for (size_t i = 0; i < total_polys; i++) {
    //   poly = polyvec_get_elem(ct0_vec_polys, i);
    //   coeffs = poly_get_coeffvec(poly);
    //   for (size_t j = 0; j < proof_degree; j++) {
    //       if (element_index < total_elements) {
    //           int64_t value = flat_static_ct0[element_index];
    //           intvec_set_elem_i64(coeffs, j, value);
    //           element_index++;
    //       } else {
    //           // Handle error: element_index exceeds total_elements
    //           fprintf(stderr, "Error: element_index exceeds total_elements\n");
    //           exit(1);
    //       }
    //   }
    // }

    for (size_t k = 0; k < CT_COUNT; k++) {
      for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
          poly = polyvec_get_elem(ct0_vec_polys, k*polys_per_ct + i);
          coeffs = poly_get_coeffvec (poly);
          for (size_t j=0; j<proof_degree ; j++) {
              intvec_set_elem_i64(coeffs,j,static_ct0[j+i*proof_degree]);
              // printf("%d, %d, %d\n", k,i,j);
          }
      }
    }
    // print_polyvec_element("First element of ct0", ct0_vec_polys, 0, 64);
    // some manual tests:
    // print_polyvec_element("First element of ct0", ct0_vec_polys, 0, 1);
    // print_polyvec_element("Last element of ct0", ct0_vec_polys, ct0_vec_polys->nelems-1, 64);
    // printf("\nNumber of polynomials in ct0_vec_polys is %d", ct0_vec_polys->nelems);
    // printf("\nNumber of polynomials in ct0_vec_polys should be %d\n", total_polys);



    /* INIT ct1 */
    // This is the same as ct0, so we need to flatten it
    // size calculations
    polys_per_ct = fhe_degree / proof_degree;
    total_polys = CT_COUNT * polys_per_ct;
    polyvec_t ct1_vec_polys;
    polyvec_alloc(ct1_vec_polys, Rq, total_polys);
    // POLYVEC_T(ct1_vec_polys, Rq, total_polys);

    // one-dimensional array in memory of C
    const int64_t *flat_static_ct1 = (const int64_t *)static_ct1;

    total_elements = CT_COUNT * fhe_degree;
    element_index = 0;

    // flatten
    // for (size_t i = 0; i < total_polys; i++) {
    //   poly = polyvec_get_elem(ct1_vec_polys, i);
    //   coeffs = poly_get_coeffvec(poly);
    //   for (size_t j = 0; j < proof_degree; j++) {
    //       if (element_index < total_elements) {
    //           int64_t value = flat_static_ct1[element_index];
    //           intvec_set_elem_i64(coeffs, j, value);
    //           element_index++;
    //       } else {
    //           // Handle error: element_index exceeds total_elements
    //           fprintf(stderr, "Error: element_index exceeds total_elements\n");
    //           exit(1);
    //       }
    //   }
    // }
    for (size_t k = 0; k < CT_COUNT; k++) {
      for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
          poly = polyvec_get_elem(ct1_vec_polys, k*polys_per_ct + i);
          coeffs = poly_get_coeffvec (poly);
          for (size_t j=0; j<proof_degree ; j++) {
              intvec_set_elem_i64(coeffs,j,static_ct1[j+i*proof_degree]);
              // printf("sk elem: %d\n", static_ct0[j+i*proof_degree]);
          }
      }
    }
    // print_polyvec_element("First element of ct1", ct1_vec_polys, 0, 64);

    // some manual tests:
    // print_polyvec_element("First element of ct1", ct1_vec_polys, 0, 1);
    // print_polyvec_element("First element of ct1", ct1_vec_polys, 32, 1);
    // //print_polyvec_element("Last element of ct1", ct1_vec_polys, ct1_vec_polys->nelems-1, 64);
    // printf("\nNumber of polynomials in ct1_vec_polys is %d", ct1_vec_polys->nelems);
    // printf("\nNumber of polynomials in ct1_vec_polys should be %d\n", total_polys);

    /* INIT m_delta */
    // This is the same as ct0, so we need to flatten it
    // size calculations
    size_t polys_per_mdelta = fhe_degree / proof_degree;
    total_polys = CT_COUNT * polys_per_mdelta;
    polyvec_t mdelta_vec_polys;
    polyvec_alloc(mdelta_vec_polys, Rq, total_polys);
    // POLYVEC_T(mdelta_vec_polys, Rq, total_polys);

    // one-dimensional array in memory of C
    const int64_t *flat_static_mdelta = (const int64_t *)static_m_delta;

    total_elements = CT_COUNT * fhe_degree;
    element_index = 0;

    // flatten
    // for (size_t i = 0; i < total_polys; i++) {
    //   poly = polyvec_get_elem(mdelta_vec_polys, i);
    //   coeffs = poly_get_coeffvec(poly);
    //   for (size_t j = 0; j < proof_degree; j++) {
    //       if (element_index < total_elements) {
    //           int64_t value = flat_static_mdelta[element_index];
    //           intvec_set_elem_i64(coeffs, j, value);
    //           element_index++;
    //       } else {
    //           // Handle error: element_index exceeds total_elements
    //           fprintf(stderr, "Error: element_index exceeds total_elements\n");
    //           exit(1);
    //       }
    //   }
    // }

    for (size_t k = 0; k < CT_COUNT; k++) {
      for (size_t i=0; i<(fhe_degree/proof_degree) ; i++) {
          poly = polyvec_get_elem(mdelta_vec_polys, k*polys_per_ct + i);
          coeffs = poly_get_coeffvec (poly);
          for (size_t j=0; j<proof_degree ; j++) {
              intvec_set_elem_i64(coeffs,j,static_m_delta[j+i*proof_degree]);
              // printf("sk elem: %d\n", static_ct0[j+i*proof_degree]);
          }
      }
    }
    // some manual tests:
    // print_polyvec_element("First element of m_delta", mdelta_vec_polys, 0, 1);
    // print_polyvec_element("Last element of m_delta", mdelta_vec_polys, mdelta_vec_polys->nelems-1, 64);
    // printf("\nNumber of polynomials in mdelta_vec_polys is %d", mdelta_vec_polys->nelems);
    // printf("\nNumber of polynomials in mdelta_vec_polys should be %d\n", total_polys);




    uint8_t seed[32] = { 0 };
    seed[0] = 2;

    vdec_lnp_tbox (seed, params1, sk_vec_polys, static_sk, ct0_vec_polys, ct1_vec_polys, 
                   mdelta_vec_polys, fhe_degree);

    mpfr_free_cache();
    printf("Finished.\n");
}


/* R2 != R2_ */
static void
_scatter_smat(spolymat_ptr R2, spolymat_ptr R2_, unsigned int m1,
              unsigned int Z, unsigned int l)
{
    const unsigned int nelems = R2_->nelems;
    unsigned int i, row, col;
    poly_ptr poly, poly2;

    //   ASSERT_ERR (R2->nelems_max >= R2_->nelems_max);
    //   ASSERT_ERR (spolymat_is_upperdiag (R2_));

    (void)l; /* unused */

    for (i = 0; i < nelems; i++)
    {
        poly = spolymat_get_elem(R2_, i);
        row = spolymat_get_row(R2_, i);
        col = spolymat_get_col(R2_, i);

        //   ASSERT_ERR (row < 2 * (m1 + l));
        //   ASSERT_ERR (col < 2 * (m1 + l));
        //   ASSERT_ERR (col >= row);

        if (col >= 2 * m1)
            col += 2 * Z;
        if (row >= 2 * m1)
            row += 2 * Z;

        poly2 = spolymat_insert_elem(R2, row, col);
        poly_set(poly2, poly);
    }
    R2->sorted = 0;
    spolymat_sort(R2);
    //   ASSERT_ERR (spolymat_is_upperdiag (R2));
}

/* r1, r1_ may not overlap */
static void
_scatter_vec(spolyvec_ptr r1, spolyvec_ptr r1_, unsigned int m1,
             unsigned int Z)
{
    const unsigned int nelems = r1_->nelems;
    unsigned int i, elem;
    poly_ptr poly, poly2;

    //   ASSERT_ERR (r1->nelems_max >= r1_->nelems_max);

    for (i = 0; i < nelems; i++)
    {
        poly = spolyvec_get_elem(r1_, i);
        elem = spolyvec_get_elem_(r1_, i);

        if (elem >= 2 * m1)
            elem += 2 * Z;

        poly2 = spolyvec_insert_elem(r1, elem);
        poly_set(poly2, poly);
    }
    r1->sorted = 1;
}

static void vdec_lnp_tbox(uint8_t seed[32], const lnp_quad_eval_params_t params,  
                          polyvec_t sk, int8_t sk_sign[], polyvec_t ct0, polyvec_t ct1, 
                          polyvec_t m_delta, unsigned int fhe_degree)
{
  struct timeval start_proof, end_proof, start_rot, end_rot, start_debug, end_debug;
  struct timeval start_bound, end_bound, start_eval, end_eval, start_quadmany, end_quadmany;
  long seconds_proof, useconds_proof, seconds_rot, useconds_rot,seconds_debug, useconds_debug;
  long seconds_bound, useconds_bound, seconds_eval, useconds_eval,seconds_quadmany, useconds_quadmany;
  double wall_time_proof, wall_time_rot, wall_time_debug;
  double wall_time_bound, wall_time_eval, wall_time_quadmany;
  gettimeofday(&start_proof, NULL);  // Start timing
  gettimeofday(&start_bound, NULL);  // Start timing
    /************************************************************************/
    /*                                                                      */
    /*    OUR CUSTOM PROOF: committing to witness + computing u vectors     */
    /*                                                                      */
    /************************************************************************/
    abdlop_params_srcptr abdlop = params->quad_eval;
    uint8_t hashp[32] = { 0 };
    uint8_t hashv[32] = { 0 };
    polyring_srcptr Rq = abdlop->ring;
    const unsigned int lambda = params->lambda;
    // const unsigned int N_ = lambda / 2;
    INT_T (lo, Rq->q->nlimbs);
    INT_T (hi, Rq->q->nlimbs);
    int b;
    // uint8_t buf[2];
    uint32_t dom;
    unsigned int i, j, k;

    polyvec_t subv;
    int_ptr coeff;
    polymat_t A1, A2prime, Bprime;
    polyvec_t s1, s2, m, tA1, tA2, tB, z1, z21, hint, h, s, tmp;
    poly_t c;
    poly_ptr poly;
    // const unsigned int n = 2 * (abdlop->m1 + abdlop->l) + params->lambda;
    // const unsigned int np = 2 * (abdlop->m1 + abdlop->l);
    int b1 = 1, b2 = 1;

    dom = 0;

    poly_alloc (c, Rq);
    polyvec_alloc (s1, Rq, abdlop->m1);
    polyvec_alloc (s2, Rq, abdlop->m2);
    polyvec_alloc (m, Rq, abdlop->l + params->lambda / 2 + 1);
    polyvec_alloc (tA1, Rq, abdlop->kmsis);
    polyvec_alloc (tA2, Rq, abdlop->kmsis);
    polyvec_alloc (tB, Rq, abdlop->l + abdlop->lext);
    polyvec_alloc (z1, Rq, abdlop->m1);
    polyvec_alloc (z21, Rq, abdlop->m2 - abdlop->kmsis);
    polyvec_alloc (hint, Rq, abdlop->kmsis);
    polyvec_alloc (h, Rq, params->lambda / 2);
    polyvec_alloc (s, Rq, 2 * (abdlop->m1 + abdlop->l));
    polyvec_alloc (tmp, Rq, 2 * (abdlop->m1 + abdlop->l));
   
    polymat_alloc (A1, Rq, abdlop->kmsis, abdlop->m1);
    polymat_alloc (A2prime, Rq, abdlop->kmsis, abdlop->m2 - abdlop->kmsis);
    polymat_alloc (Bprime, Rq, abdlop->l + abdlop->lext,
                    abdlop->m2 - abdlop->kmsis);


    polyvec_t h_our;
    polyvec_alloc (h_our, Rq, params->lambda / 2);

    // size of original bdlop message - without y's and beta's
    const unsigned int short_l = 0;

    const unsigned int d = polyring_get_deg (Rq);
    const unsigned int m1 = abdlop->m1;
    const unsigned int l = abdlop->l;
    const unsigned int nbounds = 1; // TODO: number of u vectors we want to proof are small - will change to 1
    //const unsigned int nprime = ct0->nelems * CT_COUNT;
    const unsigned int nprime = fhe_degree/d * CT_COUNT;

    int_set_i64 (lo, -1);
    int_set_i64 (hi, 1);
    polyvec_urandom_bnd (s2, lo, hi, seed, dom++);

    printf("ajtai size: %d, bdlop size: %d, lext:%d, lambda:%d\n", m1, l, abdlop->lext, lambda);
    printf("quad-many l: %d, quad-many lext:%d\n\n", params->quad_many->l, params->quad_many->lext);

    // #region Committing to witness


    // build witness and commit
    polyvec_t tobe_sk; // vectors to be committed using abdlop
    polyvec_get_subvec(tobe_sk, s1, 0, sk->nelems, 1);
    polyvec_set(tobe_sk, sk);


    // generate abdlop keys and commit to sk (ajtai part) and the vinh (bdlop part)
    abdlop_keygen (A1, A2prime, Bprime, seed, abdlop);
    abdlop_commit (tA1, tA2, tB, s1, m, s2, A1, A2prime, Bprime, abdlop);
    // printf("tA1 size:%d, tA2 size:%d, tB size:%d\n", tA1->nelems, tA2->nelems, tB->nelems);
    // printf("Bprime: %d rows, %d cols\n", Bprime->nrows, Bprime->ncols);

    // #endregion

    // #region build u vectors


    // build u vector - u_v (and temporary u_s = sk in coefficient form)
    intvec_t u_v_vec, u_s_vec;
    intvec_alloc(u_v_vec, d * ct0->nelems, Rq->q->nlimbs);
    // intvec_alloc(u_s_vec, d * sk->nelems, Rq->q->nlimbs);
    //INTVEC_T(u_v_vec, d * ct0->nelems, Rq->q->nlimbs);
    //INTVEC_T(u_s_vec, d * sk->nelems, Rq->q->nlimbs);
    intvec_ptr u_v = &u_v_vec;
    intvec_ptr u_s = &u_s_vec;

    poly_ptr poly_tmp;
    intvec_ptr coeffs;

    // u_s
    // for (i=0; i<sk->nelems; i++) {
    //     poly_tmp = polyvec_get_elem(sk, i);     
    //     coeffs = poly_get_coeffvec (poly_tmp);
    //     for (j=0; j<d; j++) {
    //         intvec_set_elem(u_s_vec, i*d+j, intvec_get_elem(coeffs, j));
    //     }
    // }

    // u_v  
    // printf("start u_v build\n");
    polyvec_t c0_m;
    polyvec_alloc(c0_m, Rq, ct0->nelems);
    printf("\nct0->nelems: %d", ct0->nelems);
    polyvec_sub(c0_m, ct0, m_delta, 0);
    printf("\nc0_m->nelems: %d", c0_m->nelems);

    // generate intvec with coeffs of ct0 - delta_m
    intvec_t sum_tmp_vec;
    intvec_alloc(sum_tmp_vec, d * c0_m->nelems, Rq->q->nlimbs);
    // INTVEC_T(sum_tmp_vec, d * c0_m->nelems, Rq->q->nlimbs);
    printf("\nsum_tmp_vec->nelems: %d\n", sum_tmp_vec->nelems);
    intvec_ptr sum_tmp = &sum_tmp_vec;
    for (i=0; i<c0_m->nelems; i++) {
        poly_tmp = polyvec_get_elem(c0_m, i); 
        coeffs = poly_get_coeffvec (poly_tmp);
        for (j=0; j<d; j++) {
            intvec_set_elem(sum_tmp, i*d+j, intvec_get_elem(coeffs, j));
        }
    }
    // intvec_dump(sum_tmp_vec);
    //print_polyvec_element("element of vinh", c0_m_v, 31, 64);
    //printf("sum_tmp nelems = %d\n", sum_tmp->nelems);
    // print_polyvec_element("element of ct0", ct0, 60, 3);
    // print_polyvec_element("element of m_delta", m_delta, 60, 3);
    // printf("\nsum_tmp:\n");
    // for (i=0; i<1; i++) {
    //     printf("%lld ", intvec_get_elem_i64(sum_tmp, 60*d+2));
    // }
    // printf("\n\n");

    // For each ct in ct1:
    // generate intvec with coeffs of ct1, do rotations and
    // dot product with u_s
    // Calculate sizes (based on the paper):
    polyvec_ptr ct1_ptr = &ct1;
    poly_ptr first_ct1_ptr = polyvec_get_elem(ct1, 0);
    size_t n = fhe_degree/d; //first_ct1_ptr->coeffs->nelems; // fhe_degree/proof_degree
    size_t r = CT_COUNT;
    printf("\nn: %d", n);
    printf("\nr (CT_COUNT): %d\n", CT_COUNT);

    // polymat_t Ds;
    // polymat_alloc (Ds, Rq, CT_COUNT*n*d, m1);

    intvec_t w_sk, rot_s_vec;
    intvec_alloc(w_sk, CT_COUNT*d*n, Rq->q->nlimbs);
    intvec_alloc(rot_s_vec, d * n, Rq->q->nlimbs);
    // INTVEC_T(w_sk, CT_COUNT*d*n, Rq->q->nlimbs);
    // INTVEC_T(rot_s_vec, d * n, Rq->q->nlimbs);
    intvec_ptr rot_s = &rot_s_vec;
    INT_T (new, 2 * Rq->q->nlimbs);

    intvec_t ct1_allcoeffs;
    intvec_alloc(ct1_allcoeffs, d * ct1->nelems, Rq->q->nlimbs);

    // gettimeofday(&start_rot, NULL);  // Start timing
    // #pragma omp parallel for private(k,i, subv) shared(ct1, ct1_allcoeffs)
    for (k=0; k<CT_COUNT; k++) {
      // getting k-th ct1 coeffs
      INTVEC_T(ct1_coeffs_vec, d * n, Rq->q->nlimbs);
      intvec_ptr ct1_coeffs = &ct1_coeffs_vec;
      for (i=0; i<n; i++) {
          poly_tmp = polyvec_get_elem(ct1, k*n + i);     
          coeffs = poly_get_coeffvec (poly_tmp);
          for (j=0; j<d; j++) {
              intvec_set_elem(ct1_coeffs, i*d+j, intvec_get_elem(coeffs, j));
          }
      }
      intvec_get_subvec(subv, ct1_allcoeffs, k*n*d, n*d, 1);
      intvec_set(subv, ct1_coeffs);
      // INTVEC_T(row, d * m1, Rq->q->nlimbs);
      // gbfv_rot_row_fast(row, ct1_coeffs, 2000, Rq);
      // for (j=0; j<10; j++)
      //     printf("%d ", intvec_get_elem_i64(row, j));
      // printf("\n");

      //test_gbfv_rot(Rq);
      if (GBFV == 1) {
        // printf("start u_v build with GBFV\n");
        // intmat_t Ds_coeffs_vec;
        // intmat_alloc(Ds_coeffs_vec, fhe_degree, m1*d, Rq->q->nlimbs);
        // intmat_ptr Ds_coeffs = &Ds_coeffs_vec;
        // INTVEC_T(ct1_coeffs_vec2, d * n, Rq->q->nlimbs);
        // intvec_ptr ct1_coeffs2 = &ct1_coeffs_vec2;
        // INTVEC_T(ct1_coeffs_vec3, d * n, Rq->q->nlimbs);
        // intvec_ptr ct1_coeffs3 = &ct1_coeffs_vec3;

        // gbfv_rot_col(ct1_coeffs2, ct1_coeffs, 2000, Rq);
        // for (j=0; j<10; j++)
        //     printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, j));
        // printf("\n");

        // gbfv_rot(Ds_coeffs, ct1_coeffs, Rq);
        // intmat_get_col(ct1_coeffs3, Ds_coeffs, 2000);
        // for (j=0; j<10; j++)
        //     printf("%lld ", intvec_get_elem_i64(ct1_coeffs3, j));
        // printf("\n");

        // printf("col is correct : %d\n", intvec_eq(ct1_coeffs3, ct1_coeffs2));


        // intvec_t rot_coeffvec;
        // poly_ptr Ds_elem;
        // rotating coeffs of k-th ct1 and multiplying with u_s
        INT_T (new, 2 * Rq->q->nlimbs);
        for (i=0; i<(d * n); i++) {
            int_set_zero(new);
            gbfv_rot_sum(new, ct1_coeffs, i, sk_sign, Rq);
            
            // intmat_get_row(ct1_coeffs2, Ds_coeffs, i);
            // // // gbfv_rot_row(ct1_coeffs2, ct1_coeffs, i, Rq); substitute previous line with this once function works.
            // for (j=0; j<(ct1_coeffs2->nelems)/d; j++) {
            //     intvec_get_subvec (rot_coeffvec, ct1_coeffs2, 0+j*d, d, 1); // XXX correct
            //     Ds_elem = polymat_get_elem(Ds, k*(n*d)+i, j);
            //     poly_set_coeffvec(Ds_elem, rot_coeffvec);
            // }

            // intvec_dot(new, ct1_coeffs2, u_s); // TO UNCOMMENT
            int_mod(new, new, Rq->q);
            int_redc(new, new, Rq->q);
            intvec_set_elem(rot_s, i, new);
        }
      }
      // else if (GBFV == 0) {
      //   printf("start u_v build with BFV\n");
      //   intvec_t rot_coeffvec;
      //   poly_ptr Ds_elem;
      //   // INTVEC_T(ct1_coeffs_vec2, d * n, Rq->q->nlimbs);
      //   intvec_t ct1_coeffs_vec2; 
      //   intvec_alloc(ct1_coeffs_vec2, d * n, Rq->q->nlimbs);
      //   intvec_ptr ct1_coeffs2 = &ct1_coeffs_vec2;
      //   // rotating coeffs of k-th ct1 and multiplying with u_s
      //   intvec_reverse(ct1_coeffs, ct1_coeffs);
      //   for (i=0; i<(d * n); i++) {
      //       intvec_lrot(ct1_coeffs2, ct1_coeffs, i+1);
      //       intvec_neg_self(ct1_coeffs2);

      //       for (j=0; j<(ct1_coeffs2->nelems)/d; j++) {
      //           intvec_get_subvec (rot_coeffvec, ct1_coeffs2, 0+j*d, d, 1); // XXX correct
      //           Ds_elem = polymat_get_elem(Ds, k*(n*d)+i, j);
      //           poly_set_coeffvec(Ds_elem, rot_coeffvec);
      //       }
      //       //printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, i));
      //       // INTVEC_T(fake_u_s, d * sk->nelems, Rq->q->nlimbs);
      //       // intvec_set_ones(fake_u_s);
      //       // intvec_dot(new, ct1_coeffs2, fake_u_s);
      //       intvec_dot(new, ct1_coeffs2, u_s); // TO UNCOMMENT
      //       // do we need to do mod and redc?
      //       //printf("new1: %lld\n", int_get_i64(new));
      //       int_mod(new, new, Rq->q);
      //       //printf("new2: %lld\n", int_get_i64(new));
      //       int_redc(new, new, Rq->q);
      //       //printf("new3: %lld\n", int_get_i64(new));
      //       intvec_set_elem(rot_s, i, new);

      //       //printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, i));
      //   }
      // }

      for (i=0; i<(n*d); i++) {
        intvec_set_elem(w_sk, k*(d*n) + i, intvec_get_elem(rot_s, i));
      }
    }
    // gettimeofday(&end_rot, NULL);  // End timing
    // // Compute the time difference in microseconds
    // seconds_rot  = end_rot.tv_sec  - start_rot.tv_sec;
    // useconds_rot = end_rot.tv_usec - start_rot.tv_usec;
    // wall_time_rot = seconds_rot + useconds_rot/1e6;  // Convert to seconds
    // printf(" ---------> (Rotations and w_sk): execution time: %f seconds\n", wall_time_rot);


    intvec_add(u_v, w_sk, sum_tmp);
    // printf("dumping w_sk\n");
    // coeff = intvec_get_elem(w_sk, 0);
    // int_dump(coeff);
    // coeff = intvec_get_elem(w_sk, n*d+0);
    // int_dump(coeff);

    // printf("dumping u_v\n");
    // for (i=0; i<10; i++)
    //     printf("%d ", intvec_get_elem_i64(u_v, i));
    // printf("\n");

    printf("finished u_v build\n");

    // #endregion

    /************************************************************************/
    /*                                                                      */
    /*       PROOF OF L2 NORM BOUND: computing z's: z_s, z_l and z_v        */
    /*                                                                      */
    /************************************************************************/
    
    // prepare randomness sent to lnp_tbox_prove from lnp-tbox-test
    memset (hashp, 0xff, 32); 
    // seed is also sent, but it is already declared before
    
    // things from lnp_tbox_prove
    shake128_state_t hstate;
    coder_state_t cstate;
    uint8_t hash0[32];
    uint8_t expseed[3 * 32];
    // const uint8_t *seed_rej34 = expseed;
    const uint8_t *seed_cont = expseed + 32;
    const uint8_t *seed_cont2 = expseed + 64;
    /* buff for encoding of tg */
    const unsigned int log2q = polyring_get_log2q (Rq);
    uint8_t out[CEIL(log2q * d * lambda / 2, 8) + 1];

    // this is done before calling compute_z34
    /*
     * Expand input seed into two seeds: one for rejection sampling on z3, z4
     * and one for continuing the protocol.
    */
    shake128_init (hstate);
    shake128_absorb (hstate, seed, 32);
    shake128_squeeze (hstate, expseed, sizeof (expseed));
    shake128_clear (hstate);
    // compute_z34 is then called. hash and seed_rej34 are sent to compute_z34
    // hstate and cstate are declared again inside compute_z34

    
    // things from compute_z34
    //const unsigned int log2q = polyring_get_log2q (Rq);  
    const unsigned int kmsis = abdlop->kmsis;
    const unsigned int m2 = abdlop->m2;
    // todo: why not instantiating s3coeffs?
    INTVEC_T (yv_coeffs, 256, int_get_nlimbs (Rq->q));
    INTVEC_T (zv_coeffs, 256, int_get_nlimbs (Rq->q));
    polyvec_t s21, yv_, tyv, tbeta, beta, yv,// s1_, m_ are probably not needed
              zv, zv_;
    // intvec_ptr coeffs; // already declared
    polymat_t Byv, Bbeta;
    shake128_state_t hstate_z34;
    coder_state_t cstate_z34;
    rng_state_t rstate_signs;
    rng_state_t rstate_rej;
    //uint32_t dom = 0; // already declared
    uint8_t rbits;
    unsigned int nrbits, outlen, loff;
    uint8_t out_z34[CEIL (256 * 2 * log2q + d * log2q, 8) + 1];
    uint8_t cseed[32]; /* seed for challenge */
    //poly_ptr poly; // already declared
    int beta_v;
    int rej;

    // what is this for??
    memset (out_z34, 0, CEIL (256 * 2 * log2q + d * log2q, 8) + 1); // XXX
    
    polyvec_alloc (yv, Rq, 256 / d);
    polyvec_alloc (zv, Rq, 256 / d);
    polyvec_alloc (zv_, Rq, 256 / d);

    // printf("allocated space for building zv\n");

    /* s1 = s1_,upsilon, m = m_,y3_,y4_,beta */
    // from lnp-tbox: s1_ m_ probably not needed
    //polyvec_get_subvec (s1_, s1, 0, m1, 1);
    //polyvec_get_subvec (m_, m, 0, l, 1);
    //printf("m length %d, l=%d, beta pos %d\n", m->nelems, l, short_l + (256 / d) * nbounds);
    polyvec_get_subvec (beta, m, short_l + (256 / d) * nbounds, 1, 1);
    polyvec_set_zero (beta);
    polyvec_get_subvec (s21, s2, 0, m2 - kmsis, 1);


    /* tB = tB_,ty,tbeta */
    loff = 0;
    polyvec_get_subvec (yv_, m, short_l + loff, 256 / d, 1);
    polyvec_get_subvec (tyv, tB, short_l + loff, 256 / d, 1);
    polymat_get_submat (Byv, Bprime, short_l + loff, 0, 256 / d, m2 - kmsis, 1, 1);
    polyvec_set_coeffvec2 (yv, yv_coeffs);
    polyvec_set_coeffvec2 (zv_, zv_coeffs);
    loff += 256 / d;

    polyvec_get_subvec (tbeta, tB, short_l + loff, 1, 1);
    polymat_get_submat (Bbeta, Bprime, short_l + loff, 0, 1, m2 - kmsis, 1, 1);
    // printf("tbeta in pos: %d out of %d\n", short_l+loff, tB->nelems);

    // printf("extracted subvecs for yv and beta commitments\n");


    // what is this for? -> rejection sampling state
    nrbits = 0;
    rng_init (rstate_rej, seed, dom++); // seed was seed_tbox in compute_z34
    rng_init (rstate_signs, seed, dom++);


    // #region Rejection 
    // --------------------------------------------------------
    // REJECTION SAMPLING BLOCK (big while loop in compute_z34)
    // --------------------------------------------------------
    while (1) // only breaks loop once rejection sampling succeeds
    {
        /* sample signs */
        if (nrbits == 0)
        {
            rng_urandom (rstate_signs, &rbits, 1);
            nrbits = 8;
        }
        // printf("sampled signs\n");

        /* yv, append to m  */
        polyvec_grandom (yv, params->log2stdev4, seed, dom++); // stdev4 or 3??
        //polyvec_set_zero(yv); // DEBUG to remove
        polyvec_set (yv_, yv);
        /* tyv */
        polyvec_set (tyv, yv);
        polyvec_addmul (tyv, Byv, s21, 0);
        polyvec_mod (tyv, tyv);
        polyvec_redp (tyv, tyv);
        /* beta_v  */
        beta_v = (rbits & (1 << (8 - nrbits + 1))) >> (8 - nrbits + 1);
        beta_v = 1 - 2 * beta_v; /* {0,1} -> {1,-1} */
        // printf ("beta4 %d\n", beta4);
        nrbits -= 1;

         

        /* tbeta */ 
        poly = polyvec_get_elem (beta, 0);
        coeffs = poly_get_coeffvec (poly);
        intvec_set_elem_i64 (coeffs, 0, beta_v);
        polyvec_set (tbeta, beta);
        polyvec_addmul (tbeta, Bbeta, s21, 0);
        polyvec_mod (tbeta, tbeta);
        polyvec_redp (tbeta, tbeta);


        /* encode ty, tbeta, hash of encoding is seed for challenges */
        coder_enc_begin (cstate_z34, out_z34);
        coder_enc_urandom3 (cstate_z34, tyv, Rq->q, log2q);
        coder_enc_urandom3 (cstate_z34, tbeta, Rq->q, log2q);
        coder_enc_end (cstate_z34);

        outlen = coder_get_offset (cstate_z34);
        ASSERT_ERR (outlen % 8 == 0);
        ASSERT_ERR (outlen / 8 <= CEIL (256 * 2 * log2q + d * log2q, 8) + 1);
        outlen >>= 3; /* nbits to nbytes */

        shake128_init (hstate_z34);
        shake128_absorb (hstate_z34, hashp, 32);
        shake128_absorb (hstate_z34, out_z34, outlen);
        shake128_squeeze (hstate_z34, cseed, 32);

        // printf("created tbeta\n");

        // calculate zv
        INT_T (beta_v_Rij_uv_j, int_get_nlimbs (Rq->q));
        int8_t Ri_v[u_v->nelems];
        int_ptr uv_coeff, R_uv_coeff;

        //polyvec_fromcrt (s3);
        polyvec_fromcrt (yv);

        polyvec_set (zv_, yv);
        intvec_set_zero (yv_coeffs);


        for (i = 0; i < 256; i++)
        {
            R_uv_coeff = intvec_get_elem (yv_coeffs, i);

            // should I use this or _expand_Rprime_i?

            //printf("Expanding...\n");
            _expand_R_i2 (Ri_v, u_v->nelems, i, cseed);
            // if (i == 2) {
            //   printf(" - - Ri 1:  ");
            //   for (j = 0; j < u_v->nelems; j++)
            //     printf("%d", Ri_v[j]);
            //   printf("\n");
            // }
            //printf("Expanded\n");

            for (j = 0; j < u_v->nelems; j++)
            {
                if (Ri_v[j] == 0)
                {
                }
                else
                {
                    ASSERT_ERR (Ri_v[j] == 1 || Ri_v[j] == -1);

                    uv_coeff = intvec_get_elem (u_v_vec, j);
                    
                    int_set (beta_v_Rij_uv_j, uv_coeff);
                    int_mul_sgn_self (beta_v_Rij_uv_j, Ri_v[j]);
                    int_add (R_uv_coeff, R_uv_coeff, beta_v_Rij_uv_j);
                }
            }
        }
        intvec_mul_sgn_self (yv_coeffs, beta_v);
        intvec_add (zv_coeffs, zv_coeffs, yv_coeffs);
        // printf("created z_v\n");



        /* rejection sampling */

        intvec_mul_sgn_self (yv_coeffs, beta_v); /* revert mul by beta3 */
        rej = rej_bimodal (rstate_rej, zv_coeffs, yv_coeffs, params1_scM4, params1_stdev4sq);
        if (rej) {
            DEBUG_PRINTF (DEBUG_PRINT_REJ, "%s", "reject u_v");
            continue;
        }
        printf("did rejection sampling for z_v\n");
        

        break;
    }

    // #endregion








    /* update fiat-shamir hash */
    memcpy (hashp, cseed, 32); // hashp is called hash in compute_z34

    /* output proof (h,c,z1,z21,hint,z3,z4) */
    polyvec_set (zv, zv_);

    /* cleanup */
    // printf("will cleanup after calculating z's\n");
    rng_clear (rstate_signs);
    rng_clear (rstate_rej);
    polyvec_free (yv);
    polyvec_free (zv_);
    // printf("finished cleaning up after z's\n");

    gettimeofday(&end_bound, NULL);  // End timing
    // Compute the time difference in microseconds
    seconds_bound  = end_bound.tv_sec  - start_bound.tv_sec;
    useconds_bound = end_bound.tv_usec - start_bound.tv_usec;
    wall_time_bound = seconds_bound + useconds_bound/1e6;  // Convert to seconds
    printf(" ---------> (Bound proof): execution time: %f seconds\n", wall_time_bound);

    /************************************************************************/
    /*                                                                      */
    /*             BUILDING STATEMENTS FOR NEXT PARTS OF PROOF              */
    /*                                                                      */
    /************************************************************************/
    gettimeofday(&start_eval, NULL);  // Start timing


    //const uint8_t *seed_cont = expseed + 32;
    polyvec_t subv2, subv_auto, tg, s2_;
    polymat_t Bextprime;

    printf("start building statements for next parts of the proof\n");
    polyvec_free (s); // freeing this because this is used in the original proof from quad_eval_test. Later we can remove this
    polyvec_alloc (s, Rq, 2 * (m1 + params->quad_many->l)); // double check this l (from quad-many and not quad)


    // stuff from lnp-tbox (with adaptation)
    memcpy (hash0, hashp, 32); // save this level of FS hash for later

    /* tB = (tB_,tg,t) */
    polyvec_get_subvec (tg, tB, l, lambda / 2, 1); // is l the correct row -> it should be
    /* Bprime = (Bprime_,Bext,bext) */
    polymat_get_submat (Bextprime, Bprime, l, 0, lambda / 2, abdlop->m2 - abdlop->kmsis, 1, 1);


    // printf("building s - calculating automorphisms\n");
    /* BUILDING: s = (<s1>,<m>,<y_v>,<beta_v>) */ // should this block be here? calculating automorphisms
    // automorphism of ajtai part (s1)
    polyvec_get_subvec (subv, s, 0, m1, 2);
    polyvec_get_subvec (subv_auto, s, 1, m1, 2);
    polyvec_set (subv, s1);
    polyvec_auto (subv_auto, s1);

    // automorphism of original bdlop part (m)
    if (short_l > 0)
    {
        polyvec_get_subvec (subv, s, (m1) * 2, short_l, 2);
        polyvec_get_subvec (subv_auto, s, (m1) * 2 + 1, short_l, 2);
        polyvec_get_subvec (subv2, m, 0, short_l, 1);
        polyvec_set (subv, subv2);
        polyvec_auto (subv_auto, subv2);
    }

    // automorphism of extended bdlop part (y_v, beta_v)
    polyvec_get_subvec (subv, s, (m1 + short_l) * 2, loff + 1, 2);
    polyvec_get_subvec (subv_auto, s, (m1 + short_l) * 2 + 1, loff + 1, 2);
    polyvec_get_subvec (subv2, m, short_l, loff + 1, 1);

    polyvec_set (subv, subv2);
    polyvec_auto (subv_auto, subv2);
    // end of block for calculating automorphisms


    // printf("generate random h with coeffs 0 and d/s == 0\n");
    /* generate uniformly random h=g with coeffs 0 and d/2 == 0 */
    for (i = 0; i < lambda / 2; i++)
    {
        poly = polyvec_get_elem (h_our, i);
        coeffs = poly_get_coeffvec (poly);

        intvec_urandom (coeffs, Rq->q, log2q, seed_cont, i);
        intvec_set_elem_i64 (coeffs, 0, 0);
        intvec_set_elem_i64 (coeffs, d / 2, 0);
    }


    // printf("append g to bdlop part and commit\n");
    /* append g to message m */
    polyvec_get_subvec (subv, m, l, lambda / 2, 1);
    polyvec_set (subv, h_our);

    /* tg = Bexptprime*s2 + g */
    polyvec_set (tg, h_our);
    polyvec_get_subvec (s2_, s2, 0, abdlop->m2 - abdlop->kmsis, 1);
    polyvec_addmul (tg, Bextprime, s2_, 0);

    /* encode and hash tg */
    polyvec_mod (tg, tg);
    polyvec_redp (tg, tg);

    coder_enc_begin (cstate, out);
    coder_enc_urandom3 (cstate, tg, Rq->q, log2q);
    coder_enc_end (cstate);

    outlen = coder_get_offset (cstate);
    outlen >>= 3; /* nbits to nbytes */

    shake128_init (hstate);
    shake128_absorb (hstate, hashp, 32);
    shake128_absorb (hstate, out, outlen);
    shake128_squeeze (hstate, hashp, 32);
    shake128_clear (hstate);


    // printf("going into quad and quad_eval eqs setup\n");
    // instantiating QUAD + QUAD_EVAL eqs
    spolymat_t R2t;
    spolyvec_t r1t;
    poly_t r0t;
    const unsigned int n_ = 2 * (m1 + l);
    const unsigned int np2 = 2 * (m1 + params->quad_many->l);
    spolymat_ptr R2prime_sz[lambda / 2 + 1], R2primei; // double check: the +1 should be for beta
    spolyvec_ptr r1prime_sz[lambda / 2 + 1], r1primei; 
    poly_ptr r0prime_sz[lambda / 2 + 1], r0primei; 
    spolymat_ptr R2prime_sz2[lambda / 2];
    spolyvec_ptr r1prime_sz2[lambda / 2];
    poly_ptr r0prime_sz2[lambda / 2];
    const unsigned int ibeta = (m1 + short_l + loff) * 2;

    // printf("allocate space for 1 quad eq - R2t, r1t, r0t\n");
    /* allocate tmp space for 1 quadrativ eq */
    spolymat_alloc (R2t, Rq, n_, n_, NELEMS_DIAG (n_));
    spolymat_set_empty (R2t);

    spolyvec_alloc (r1t, Rq, n_, n_);
    spolyvec_set_empty (r1t);

    poly_alloc (r0t, Rq);
    poly_set_zero (r0t);

    // printf("Allocate lambda/2 eqs for sz accumulators\n");
    /* allocate lambda/2 eqs (schwarz-zippel accumulators) */
    for (i = 0; i < lambda / 2; i++) {
      R2primei = alloc_wrapper(sizeof (spolymat_t));
      // printf("Allocating\n");
      spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
      R2prime_sz[i] = R2primei;
      spolymat_set_empty (R2prime_sz[i]);

      R2prime_sz[i]->nrows = n_;
      R2prime_sz[i]->ncols = n_;
      R2prime_sz[i]->nelems_max = NELEMS_DIAG (n_);

      r1primei = alloc_wrapper (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np2, np2);
      r1prime_sz[i] = r1primei;
      spolyvec_set_empty (r1prime_sz[i]);

      r1prime_sz[i]->nelems_max = n_;

      r0primei = alloc_wrapper (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz[i] = r0primei;
      poly_set_zero (r0prime_sz[i]);
    }
    // printf("Allocate lambda/2 eqs for sz accumulators - part2\n");
    for (i = 0; i < lambda / 2; i++) {
      R2primei = alloc_wrapper (sizeof (spolymat_t));
      spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
      R2prime_sz2[i] = R2primei;
      spolymat_set_empty (R2prime_sz2[i]);

      R2prime_sz2[i]->nrows = n_;
      R2prime_sz2[i]->ncols = n_;
      R2prime_sz2[i]->nelems_max = NELEMS_DIAG (n_);

      r1primei = alloc_wrapper (sizeof (spolyvec_t));
      spolyvec_alloc (r1primei, Rq, np2, np2);
      r1prime_sz2[i] = r1primei;
      spolyvec_set_empty (r1prime_sz2[i]);

      r1prime_sz[i]->nelems_max = n_;

      r0primei = alloc_wrapper (sizeof (poly_t));
      poly_alloc (r0primei, Rq);
      r0prime_sz2[i] = r0primei;
      poly_set_zero (r0prime_sz2[i]);
    }

    /* set up quad eqs for lower level protocol */
    // allocate 2 eqs in beta,o(beta): (actually 1 in our case - only beta4)
    // prove beta3^2-1=0 over Rq -> (i2*beta+i2*o(beta))^2 - 1 = i4*beta^2 +
    // i2*beta*o(beta) + i4*o(beta)^2 - 1 == 0 terms: R2: 3, r1: 0, r0: 1 | * 1
    // prove beta4^2-1=0 over Rq -> (-i2*x^(d/2)*beta+i2*x^(d/2)*o(beta))^2 =
    // -i4*beta^2 + i2*beta*o(beta) -i4*o(beta)^2 - 1 == 0 terms: R2: 3, r1: 0,
    // r0: 1 | * 1

    // printf("Setting up quad eqs for beta_v\n");
    i = lambda / 2;

    R2primei = alloc_wrapper (sizeof (spolymat_t));
    spolymat_alloc (R2primei, Rq, np2, np2, NELEMS_DIAG (np2));
    R2prime_sz[i] = R2primei;
    spolymat_set_empty (R2prime_sz[i]);
    poly = spolymat_insert_elem (R2prime_sz[i], ibeta, ibeta);
    poly_set_zero (poly);
    coeff = poly_get_coeff (poly, 0);
    int_set (coeff, params->inv4); // double check inv4 which was added to quad_eval_params
    poly = spolymat_insert_elem (R2prime_sz[i], ibeta, ibeta + 1);
    poly_set_zero (poly);
    coeff = poly_get_coeff (poly, 0);
    int_set (coeff, Rq->inv2);
    poly = spolymat_insert_elem (R2prime_sz[i], ibeta + 1, ibeta + 1);
    poly_set_zero (poly);
    coeff = poly_get_coeff (poly, 0);
    int_set (coeff, params->inv4);
    R2prime_sz[i]->sorted = 1;

    r1prime_sz[i] = NULL;

    r0primei = alloc_wrapper (sizeof (poly_t));
    poly_alloc (r0primei, Rq);
    r0prime_sz[i] = r0primei;
    poly_set_zero (r0prime_sz[i]);
    coeff = poly_get_coeff (r0prime_sz[i], 0);
    int_set_i64 (coeff, -1);



    /* creating Ds, Dm, u .. */ 
    // (now done inside sz accumulate function)
    // polymat_t Dm;
    // polymat_t oDs;
    // polymat_t oDm;

    // if (Ds != NULL)
    // {
    //   polymat_alloc (oDs, Rq, nprime, m1);
    //   polymat_auto (oDs, Ds);
    // }
    // if (short_l > 0 && Dm != NULL)
    // {
    //   polymat_alloc (oDm, Rq, nprime, l);
    //   polymat_auto (oDm, Dm);
    // }



    /* accumulate schwarz-zippel .. */
    // #region sz accumulate

    printf("accumulating beta...\n");
    __schwartz_zippel_accumulate_beta ( // what should be here instead of params?
          R2prime_sz, r1prime_sz, r0prime_sz, R2prime_sz2, r1prime_sz2,
          r0prime_sz2, R2t, r1t, r0t, hashp, 0, params, nprime);



    printf("accumulating z4...\n");
    __schwartz_zippel_accumulate_z (R2prime_sz, r1prime_sz, r0prime_sz,
                                     R2prime_sz2, r1prime_sz2, r0prime_sz2,
                                     R2t, r1t, r0t, sum_tmp, zv,
                                     ct1_allcoeffs, hash0, d - 1, params, nprime);

    printf("schwartz zippel auto...\n");
    for (i = 0; i < lambda / 2; i++) {
      __schwartz_zippel_auto (R2prime_sz[i], r1prime_sz[i], r0prime_sz[i],
                              R2prime_sz2[i], r1prime_sz2[i], r0prime_sz2[i],
                              params);
    }


    POLY_T (tmp1, Rq);
    /* compute/output hi and set up quadeqs for lower level protocol */
    for (i = 0; i < lambda / 2; i++) {
        polyvec_get_subvec (subv, s, 0, n_, 1);

        __evaleq (tmp1, R2prime_sz[i], r1prime_sz[i], r0prime_sz[i], subv);
        poly = polyvec_get_elem (h_our, i); /* gi */
        poly_add (poly, poly, tmp1, 0);  /* hi = gi + schwarz zippel */

        /* build quadeqs */
        DEBUG_PRINTF (DEBUG_LEVEL >= 2, "set up quadeq %u", i);

        /* r0 */
        poly_sub (r0prime_sz[i], r0prime_sz[i], poly, 0); /* r0i -= -hi */

        /* r1 */
        r1prime_sz[i]->nelems_max = np2;
        poly = spolyvec_insert_elem (r1prime_sz[i],
                                     2 * (abdlop->m1 + abdlop->l + i));
        poly_set_one (poly);
        r1prime_sz[i]->sorted = 1; /* above appends */

        /* R2 only grows by lambda/2 zero rows/cols */
        R2prime_sz[i]->nrows = np2;
        R2prime_sz[i]->ncols = np2;
        R2prime_sz[i]->nelems_max = NELEMS_DIAG (np2);
    }
    // # endregion
    gettimeofday(&end_eval, NULL);  // End timing
    // Compute the time difference in microseconds
    seconds_eval  = end_eval.tv_sec  - start_eval.tv_sec;
    useconds_eval = end_eval.tv_usec - start_eval.tv_usec;
    wall_time_eval = seconds_eval + useconds_eval/1e6;  // Convert to seconds
    printf(" ---------> (Eval proof): execution time: %f seconds\n", wall_time_eval);


    gettimeofday(&start_quadmany, NULL);  // Start timing
    printf("quad many prove...\n");
    memcpy(hashv, hashp, 32);
    lnp_quad_many_prove (hashp, tB, c, z1, z21, hint, s1, m, s2, tA2, A1, A2prime,
                        Bprime, R2prime_sz, r1prime_sz, lambda / 2 + 1,
                        seed_cont2, params->quad_many);
    printf("finished proof generation\n\n");

    gettimeofday(&end_quadmany, NULL);  // End timing
    // Compute the time difference in microseconds
    seconds_quadmany  = end_quadmany.tv_sec  - start_quadmany.tv_sec;
    useconds_quadmany = end_quadmany.tv_usec - start_quadmany.tv_usec;
    wall_time_quadmany = seconds_quadmany + useconds_quadmany/1e6;  // Convert to seconds
    printf(" ---------> (Quadmany proof): execution time: %f seconds\n", wall_time_quadmany);


    gettimeofday(&end_proof, NULL);  // End timing
    // Compute the time difference in microseconds
    seconds_proof  = end_proof.tv_sec  - start_proof.tv_sec;
    useconds_proof = end_proof.tv_usec - start_proof.tv_usec;
    wall_time_proof = seconds_proof + useconds_proof/1e6;  // Convert to seconds
    printf(" --------->  Proof generation:  %f seconds\n", wall_time_proof);


    /************************************************************************/
    /*                                                                      */
    /*                            OUR VERIFICATION                          */
    /*                                                                      */
    /************************************************************************/

    // INT_T (linf, int_get_nlimbs (Rq->q));
    // polyvec_fromcrt (zv);
    // polyvec_linf (linf, zv);
    // b = (int_le (linf, params->Bz4));
    // // int_dump(linf);
    // // int_dump(params->Bz4);
    // printf("--> zv bound verification result: %d\n", b);

    INT_T (l2, 2*int_get_nlimbs (Rq->q));
    polyvec_fromcrt (zv);
    polyvec_l2sqr (l2, zv);
    b = (int_le (l2, params->Bz4));
    // int_dump(l2);
    // int_dump(params->Bz4);
    printf("--> zv bound verification result: %d\n", b);

    // int b1 = 1, b2 = 1;
    b1 = 1; b2 = 1;
    for (i = 0; i < lambda / 2; i++)
    {
      poly = polyvec_get_elem (h_our, i);
      coeff = poly_get_coeff (poly, 0);
      if (int_eqzero (coeff) != 1) {
        b1 = 0;
        printf("coeff 0 is ");
        int_dump(coeff);
      }
      coeff = poly_get_coeff (poly, d / 2);
      if (int_eqzero (coeff) != 1) {
        b2 = 0;
        printf("coeff d/2 is ");
        int_dump(coeff);
      }
    }
    printf("--> h_our coeff verification result: %d, %d\n", b1, b2);


    /* expect successful verification */
    //memset (hashv, 0xff, 32);
    //printf("verifying our quad_many proof\n");
    b = lnp_quad_many_verify (hashv, c, z1, z21, hint, tA1, tB, A1, A2prime,
                              Bprime, R2prime_sz, r1prime_sz, r0prime_sz,
                              lambda / 2 + 1, params->quad_many);

    printf("--> quad_many verification result: %d\n", b);

    /************************************************************************/
    /*                                                                      */
    /*                        END OF OUR CUSTOM PROOF                       */
    /*                                                                      */
    /************************************************************************/

    poly_free (c);
    polyvec_free (s1);
    polyvec_free (s2);
    polyvec_free (m);
    polyvec_free (tA1);
    polyvec_free (tA2);
    polyvec_free (tB);
    polyvec_free (z1);
    polyvec_free (z21);
    polyvec_free (hint);
    polyvec_free (h);
    polyvec_free (s);
    polyvec_free (tmp);
    polymat_free (A1);
    polymat_free (A2prime);
    polymat_free (Bprime);

}


// Function to print an array of uint8_t values with a description
void print_uint8_array(const char *description, const uint8_t *array, size_t length)
{
    // Print the description
    printf("\n%s = ", description);

    // Loop through the array and print each value
    for (size_t i = 0; i < length; i++)
    {
        printf("%u ", array[i]);
    }

    // Print a new line at the end
    printf("\n");
}

// Function to print an array of int64_t values with a description
void print_int64_array(const char *description, const int64_t *array, size_t length)
{
    // Print the description
    printf("\n%s: ", description);

    // Loop through the array and print each value
    for (size_t i = 0; i < length; i++)
    {
        printf("%lld ", (long long)array[i]);
    }

    // Print a new line at the end
    printf("\n");
}

// Function to print an intvec_t values with a description
void print_polyvec_element(const char *description, const polyvec_t vec, size_t pos, size_t length)
{
    // Print the description
    printf("\n%s: ", description);

    // Loop through the array and print each value
    poly_ptr poly;
    intvec_ptr coeffs;
    poly = polyvec_get_elem(vec, pos);
    coeffs = poly_get_coeffvec (poly);
    for (size_t i = 0; i < length; i++)
    {
        printf("%lld ", (long long)intvec_get_elem_i64(coeffs, i));
    }

    // Print a new line at the end
    printf("\n");
}

void
intvec_lrot_pos (intvec_t r, const intvec_t a, unsigned int n)
{
  const unsigned int nelems = r->nelems;
  unsigned int i;
  int_ptr t;
  INTVEC_T (tmp, r->nelems, r->nlimbs);

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);
  ASSERT_ERR (n < nelems);

  for (i = 1; i <= n; i++)
    {
      t = intvec_get_elem (tmp, n - i);
      int_set (t, intvec_get_elem_src (a, nelems - i));
    }
  for (i = n; i < nelems; i++)
    {
      intvec_set_elem (tmp, i, intvec_get_elem_src (a, i - n));
    }

  intvec_set (r, tmp);
}

void
intvec_reverse (intvec_t r, const intvec_t a)
{
  const unsigned int nelems = r->nelems;
  unsigned int i;
  //int_ptr t;
  INTVEC_T (tmp, r->nelems, r->nlimbs);

  ASSERT_ERR (r->nelems == a->nelems);
  ASSERT_ERR (r->nlimbs == a->nlimbs);

  for (i = 0; i < nelems; i++)
    {
      intvec_set_elem (tmp, i, intvec_get_elem_src (a, nelems-1-i));
    }

  intvec_set (r, tmp);
}

/* expand i-th row of R from cseed and i */
static inline void
_expand_R_i2 (int8_t *Ri, unsigned int ncols, unsigned int i,
             const uint8_t cseed[32])
{
//   _brandom (Ri, ncols, 1, cseed, i);
    brandom_wrapper(Ri, ncols, 1, cseed, i);
}

void test_lrot (const polyring_t ring) {
    printf("\n- Testig rotation ...\n");
    poly_t poly, poly2, poly3;
    poly_alloc(poly, ring);
    poly_alloc(poly2, ring);
    poly_alloc(poly3, ring);

    // Create a vector with 5 elements
    INTVEC_T(my_vector, ring->d, ring->q->nlimbs);
    // Set elements to 1, 2, 3, 4, 5
    for (size_t i = 0; i < ring->d; i++)
    {
        intvec_set_elem_i64(my_vector, i, i+1);
    }
    intvec_set_elem_i64(my_vector, 0, 1);
    printf("my_vector1: " );
    intvec_dump(my_vector); // prints vector

    poly_set_coeffvec(poly, my_vector);
    for (size_t i = 0; i < poly->coeffs->nelems; i++)
        printf("%lld ", (long long)intvec_get_elem_i64(poly->coeffs, i));
    printf("\n");

    INTVEC_T(my_vector2, ring->d, ring->q->nlimbs);
    // Set elements to 1, 2, 3, 4, 5
    for (size_t i = 0; i < ring->d; i++)
    {
        intvec_set_elem_i64(my_vector2, i, i+1);
    }
    intvec_set_elem_i64(my_vector2, 0, 1);
    printf("my_vector2: " );
    intvec_dump(my_vector2); // prints vector

    poly_set_coeffvec(poly2, my_vector2);
    for (size_t i = 0; i < poly2->coeffs->nelems; i++)
        printf("%lld ", (long long)intvec_get_elem_i64(poly2->coeffs, i));
    printf("\n");

    INTVEC_T(my_vector_rotated, ring->d, ring->q->nlimbs);
    intvec_ptr my_vector_rotated_ptr = &my_vector_rotated;


    intvec_reverse(my_vector, my_vector);

    // Rot(my_vector) * my_vector2
    INTVEC_T(rot_s_vec, ring->d, ring->q->nlimbs);
    INT_T (new, 2*ring->q->nlimbs);
    for (int i=0; i<my_vector_rotated_ptr->nelems; i++) {
        intvec_lrot(my_vector_rotated_ptr, my_vector, i+1);
        intvec_neg_self(my_vector_rotated_ptr);
        //intvec_dump(my_vector_rotated);
        intvec_dot(new, my_vector_rotated_ptr, my_vector2);

        // do we need to do mod and redc?
        //printf("new1: %lld\n", int_get_i64(new));
        //int_mod(new, new, ring->q);
        //printf("new2: %lld\n", int_get_i64(new));
        //int_redp(new, new, ring->q);
        //printf("new3: %lld\n", int_get_i64(new));
        intvec_set_elem(&rot_s_vec, i, new);

        //printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, i));
    }

    printf("Rot(my_vector) * my_vector2: " );
    intvec_dump(rot_s_vec); // prints vector
    printf("\n");

    poly_mul(poly3, poly, poly2);
    poly_redc(poly3, poly3);
    intvec_dump(poly3->coeffs);

    printf("- Testig lrot ended.\n");
}

void test_gbfv_rot (const polyring_t ring) {
    printf("\n- Testig rotation ...\n");
    poly_t poly, poly2, poly3;
    poly_alloc(poly, ring);
    poly_alloc(poly2, ring);
    poly_alloc(poly3, ring);

    // prepare materials from polynomial modulus
    // INTVEC_T(signs, 7, ring->q->nlimbs);
    // intvec_ptr signs_ptr = &signs;
    // INTVEC_T(signs_scaled, 7, 2*ring->q->nlimbs);
    // intvec_ptr signs_scaled_ptr = &signs_scaled;
    int offsets[] = {3, 5, 7, 9, 10, 15, 62};
    int signs[] = {1, -1, 1, -1, 1, -1, 1};
    // intvec_set_elem_i64(signs, 0, 1);
    // intvec_set_elem_i64(signs, 1, -1);
    // intvec_set_elem_i64(signs, 2, 1);
    // intvec_set_elem_i64(signs, 3, -1);
    // intvec_set_elem_i64(signs, 4, 1);
    // intvec_set_elem_i64(signs, 5, -1);
    // intvec_set_elem_i64(signs, 6, 1);

    //intvec_dump(signs);


    // Create a vector with 5 elements
    INTVEC_T(my_vector, ring->d, ring->q->nlimbs);
    // Set elements to 1, 2, 3, 4, 5
    for (size_t i = 0; i < ring->d; i++)
    {
        intvec_set_elem_i64(my_vector, i, i+1);
    }
    intvec_set_elem_i64(my_vector, 0, 1);
    printf("my_vector1: " );
    intvec_dump(my_vector); // prints vector

    poly_set_coeffvec(poly, my_vector);
    for (size_t i = 0; i < poly->coeffs->nelems; i++)
        printf("%lld ", (long long)intvec_get_elem_i64(poly->coeffs, i));
    printf("\n");

    INTVEC_T(my_vector2, ring->d, ring->q->nlimbs);
    // Set elements to 1, 2, 3, 4, 5
    for (size_t i = 0; i < ring->d; i++)
    {
        intvec_set_elem_i64(my_vector2, i, i+1);
    }
    intvec_set_elem_i64(my_vector2, 0, 1);
    printf("my_vector2: " );
    intvec_dump(my_vector2); // prints vector

    poly_set_coeffvec(poly2, my_vector2);
    for (size_t i = 0; i < poly2->coeffs->nelems; i++)
        printf("%lld ", (long long)intvec_get_elem_i64(poly2->coeffs, i));
    printf("\n");

    INTVEC_T(my_vector_rotated, ring->d, ring->q->nlimbs);
    intvec_ptr my_vector_rotated_ptr = &my_vector_rotated;


    //intvec_reverse(my_vector, my_vector);
    int_ptr multiplier, coeff1, coeff2;
    // Rot(my_vector) * my_vector2
    INTVEC_T(rot_s_vec, ring->d, ring->q->nlimbs);
    INT_T (new, 2*ring->q->nlimbs);
    for (int i=0; i<my_vector_rotated_ptr->nelems; i++) {
        multiplier = intvec_get_elem(my_vector_rotated, 0);
        int_dump(multiplier);
        intvec_lrot(my_vector_rotated_ptr, my_vector, i+1);
        int_dump(multiplier);

        //intvec_scale(signs_scaled_ptr, multiplier, signs_ptr);
        //intvec_dump(signs_scaled);
        printf("len signs = %d\n",sizeof(signs)/sizeof(signs[0]));

        for (int j=0; j<sizeof(signs)/sizeof(signs[0]); j++) {
          //printf("i,j = %d, %d\n", i, j);
          coeff1 = intvec_get_elem(my_vector_rotated_ptr, offsets[j]);
          //int_dump(coeff1);
          int_set(new, multiplier);
          int_mul_sgn_self(new, signs[j]);
          //coeff2 = intvec_get_elem(signs_scaled_ptr, j);
          //int_dump(coeff2);
          // printf("hi5");
          int_add(coeff1, coeff1, new);
          //intvec_set_elem(my_vector_rotated_ptr, offsets+1, new);
        }

        //intvec_neg_self(my_vector_rotated_ptr);
        intvec_dump(my_vector_rotated);
        //intvec_dot(new, my_vector_rotated_ptr, my_vector2);
        
        // do we need to do mod and redc?
        //printf("new1: %lld\n", int_get_i64(new));
        //int_mod(new, new, ring->q);
        //printf("new2: %lld\n", int_get_i64(new));
        //int_redp(new, new, ring->q);
        //printf("new3: %lld\n", int_get_i64(new));
        //intvec_set_elem(&rot_s_vec, i, new);

        //printf("%lld ", intvec_get_elem_i64(ct1_coeffs2, i));
    }

    // printf("Rot(my_vector) * my_vector2: " );
    // intvec_dump(rot_s_vec); // prints vector
    // printf("\n");

    // poly_mul(poly3, poly, poly2);
    // poly_redc(poly3, poly3);
    // intvec_dump(poly3->coeffs);

    printf("- Testig lrot ended.\n");
}

void gbfv_rot (intmat_t r, const intvec_t c1, const polyring_t ring) {
    printf("\n- Building GBFV rotation ...\n");

    intmat_t tmp_mat;
    intmat_alloc(tmp_mat, r->nrows, r->ncols, ring->q->nlimbs);

    intvec_t colvec;
    intvec_alloc(colvec, r->nrows, ring->q->nlimbs);
    //INTVEC_T(colvec, r->nrows, ring->q->nlimbs);
    intvec_ptr col = &colvec;

    intmat_get_col(col, tmp_mat, 0);
    intvec_set(col, c1);

    // prepare materials from polynomial modulus
    int offsets[] = {1024, 3072, 4096, 6144, 8192, 9216, 11264};
    int signs[] = {1, -1, -1, 1, -1, -1, 1};
    // int offsets[] = {0, 0, 0, 0, 0, 0, 0};
    // int signs[] = {1, 0, 0, 0, 0, 0, 0};


    // prepare materials from polynomial modulus
    intvec_t my_vector_rotated;
    intvec_alloc(my_vector_rotated, r->nrows, ring->q->nlimbs);
    //INTVEC_T(my_vector_rotated, r->nrows, ring->q->nlimbs);
    intvec_ptr my_vector_rotated_ptr = &my_vector_rotated;

    intvec_t previous_rotated;
    intvec_alloc(previous_rotated, r->nrows, ring->q->nlimbs);
    //INTVEC_T(previous_rotated, r->nrows, ring->q->nlimbs);
    intvec_ptr previous_rotated_ptr = &previous_rotated;

    //intvec_reverse(my_vector, my_vector);
    int_ptr multiplier, coeff1, coeff2;
    INT_T (new, 2*ring->q->nlimbs);
    for (int i=0; i<my_vector_rotated_ptr->nelems-1; i++) {
        //printf("i = %d\n", i);
        multiplier = intvec_get_elem(my_vector_rotated, 0);
        //int_dump(multiplier);
        intmat_get_col(previous_rotated, tmp_mat, i);
        intvec_lrot(my_vector_rotated_ptr, previous_rotated, 1);
        //int_dump(multiplier);

        for (int j=0; j<sizeof(signs)/sizeof(signs[0]); j++) {
          //if (signs[j] == 1 || signs[j] == -1) {
          // printf("i,j = %d, %d\n", i, j);
          coeff1 = intvec_get_elem(my_vector_rotated_ptr, offsets[j]);
          //int_dump(coeff1);
          int_set(new, multiplier);
          int_mul_sgn_self(new, signs[j]);
          //int_dump(coeff2);
          int_add(coeff1, coeff1, new);
          // int_mod(coeff1, coeff1, ring->q);
          // int_redc(coeff1, coeff1, ring->q);
          //}
        }
        //intvec_dump(my_vector_rotated);

        intmat_get_col(col, tmp_mat, i+1);
        intvec_set(col, my_vector_rotated);
    }

    // intmat_get_col(col, tmp_mat, 0);
    // for (int i=0; i<10; i++)
    //   printf("%d ", intvec_get_elem_i64(col, i));
    // printf("\n\n");

    intmat_set(r, tmp_mat);
    printf("- Building rot ended.\n");
}

void gbfv_rot_col (intvec_t row, const intvec_t c1, unsigned int l, const polyring_t ring) {
    // printf("\n- Building GBFV row (2) %d ...\n", l);

    // prepare materials from polynomial modulus
    int delta = 1024;
    int ind = 0;
    int n = c1->nelems;

    int matrix[12][12] = {{-1, 1, -1, 0, 0, 0, 0, -1, 1, -1, 0, 0}, 
                          {-1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0}, 
                          {0, -1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1}, 
                          {1, -1, 0, 0, 0, -1, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, 
                          {0, -1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0}, 
                          {1, -1, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 1, -1, 1, 0, 0, 0}};

    // result vector
    intvec_t tmp_vec;
    intvec_alloc(tmp_vec, c1->nelems, ring->q->nlimbs);

    INT_T(tmp, ring->q->nlimbs);
    INT_T(tmp2, ring->q->nlimbs);
    int_ptr coeff, coeff2;

    // int rk = l / delta;
    int startj, endj;

    for (int i=0; i<n; i++) {
        // printf("i = %d\n", i);
        int_set_zero(tmp);
        int rk = i / delta;
        // if (i<=20)
        //   printf("rk: %d\n", rk);

        // get term in x^k from x^i*c1(x)
        if (i >= l) {
          coeff2 = intvec_get_elem(c1, i-l);
          int_set(tmp, coeff2);
        }

        startj = ceil((float)(n-i) / (float)delta);
        endj = floor((float)(n-i+l-1) / (float)delta);
        // if (i<=20)
        //   printf("limits %d, %d\n", startj, endj);

        for (int j=startj; j<=endj; j++) {
          if (matrix[rk][j-startj] != 0) {
            // printf("getting element %d\n", l+j*delta-i);
            coeff2 = intvec_get_elem(c1, i+j*delta-l);
            int_set(tmp2, coeff2);
            // printf("getting sign %d, %d\n", rk, j-startj);
            int_mul_sgn_self(tmp2, matrix[rk][j-startj]);
            int_add(tmp, tmp, tmp2);
            // if (i<=10)
            //   int_dump(tmp);
          }
        }

        coeff2 = intvec_get_elem(tmp_vec, i);
        int_set(coeff2, tmp);
    }
    
    intvec_set(row, tmp_vec);
    intvec_free(tmp_vec);
    // printf("- Building rot row ended.\n");
}


void gbfv_rot_col_degree2 (intvec_t row, const intvec_t c1, unsigned int l, const polyring_t ring) {
    // printf("\n- Building GBFV row (2) %d ...\n", l);

    // prepare materials from polynomial modulus
    int delta = 256;
    int ind = 0;
    int n = c1->nelems;

    int matrix[12][12] = {{-1, 1, -1, 0, 0, 0, 0, -1, 1, -1, 0, 0}, 
                          {-1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0}, 
                          {0, -1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1}, 
                          {1, -1, 0, 0, 0, -1, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, 
                          {0, -1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0}, 
                          {1, -1, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 1, -1, 1, 0, 0, 0}};

    // result vector
    intvec_t tmp_vec;
    intvec_alloc(tmp_vec, c1->nelems, ring->q->nlimbs);

    INT_T(tmp, ring->q->nlimbs);
    INT_T(tmp2, ring->q->nlimbs);
    int_ptr coeff, coeff2;

    // int rk = l / delta;
    int startj, endj;

    for (int i=0; i<n; i++) {
        // printf("i = %d\n", i);
        int_set_zero(tmp);
        int rk = i / delta;
        // if (i<=20)
        //   printf("rk: %d\n", rk);

        // get term in x^k from x^i*c1(x)
        if (i >= l) {
          coeff2 = intvec_get_elem(c1, i-l);
          int_set(tmp, coeff2);
        }

        startj = ceil((float)(n-i) / (float)delta);
        endj = floor((float)(n-i+l-1) / (float)delta);
        // if (i<=20)
        //   printf("limits %d, %d\n", startj, endj);

        for (int j=startj; j<=endj; j++) {
          if (matrix[rk][j-startj] != 0) {
            // printf("getting element %d\n", l+j*delta-i);
            coeff2 = intvec_get_elem(c1, i+j*delta-l);
            int_set(tmp2, coeff2);
            // printf("getting sign %d, %d\n", rk, j-startj);
            int_mul_sgn_self(tmp2, matrix[rk][j-startj]);
            int_add(tmp, tmp, tmp2);
            // if (i<=10)
            //   int_dump(tmp);
          }
        }

        coeff2 = intvec_get_elem(tmp_vec, i);
        int_set(coeff2, tmp);
    }
    
    intvec_set(row, tmp_vec);
    intvec_free(tmp_vec);
    // printf("- Building rot row ended.\n");
}


void gbfv_rot_row (intvec_t row, const intvec_t c1, unsigned int n, const polyring_t ring) {
    printf("\n- Building GBFV row ...\n");

    intvec_t tmp_mat;
    intvec_alloc(tmp_mat, row->nelems, ring->q->nlimbs);

    intvec_t colvec;
    intvec_alloc(colvec, row->nelems, ring->q->nlimbs);
    //INTVEC_T(colvec, r->nrows, ring->q->nlimbs);
    intvec_ptr col = &colvec;

    int_ptr coeff, coeff2;
    coeff = intvec_get_elem(row, 0);
    coeff2 = intvec_get_elem(c1, n);
    int_set(coeff, coeff2);

    // prepare materials from polynomial modulus
    int offsets[] = {1024, 3072, 4096, 6144, 8192, 9216, 11264};
    int signs[] = {1, -1, -1, 1, -1, -1, 1};
    // int offsets[] = {0, 0, 0, 0, 0, 0, 0};
    // int signs[] = {1, 0, 0, 0, 0, 0, 0};

    printf("prepare materials\n");
    // prepare materials from polynomial modulus
    intvec_t my_vector_rotated;
    intvec_alloc(my_vector_rotated, row->nelems, ring->q->nlimbs);
    //INTVEC_T(my_vector_rotated, r->nrows, ring->q->nlimbs);
    intvec_ptr my_vector_rotated_ptr = &my_vector_rotated;

    intvec_t previous_rotated;
    intvec_alloc(previous_rotated, row->nelems, ring->q->nlimbs);
    //INTVEC_T(previous_rotated, r->nrows, ring->q->nlimbs);
    intvec_ptr previous_rotated_ptr = &previous_rotated;

    //intvec_reverse(my_vector, my_vector);
    int_ptr multiplier, coeff1;
    INT_T (new, 2*ring->q->nlimbs);
    for (int i=0; i<my_vector_rotated_ptr->nelems-1; i++) {
        //printf("i = %d\n", i);
        multiplier = intvec_get_elem(my_vector_rotated, 0);
        //int_dump(multiplier);
        //intmat_get_col(previous_rotated, tmp_mat, i);

        if (i == 0) {
          intvec_set(previous_rotated, c1);
        }
        else
          intvec_set(previous_rotated, col);

        // if (i < 4) {
        // for (int j=0; j<10; j++)
        //   printf("%lld ", intvec_get_elem_i64(col, j));
        // printf("\n");
        // }
        // if (i < 4) {
        // for (int j=0; j<10; j++)
        //   printf("%lld ", intvec_get_elem_i64(previous_rotated, j));
        // printf("\n\n");
        // }

        intvec_lrot(my_vector_rotated_ptr, previous_rotated, 1);
        //int_dump(multiplier);

        for (int j=0; j<sizeof(signs)/sizeof(signs[0]); j++) {
          //if (signs[j] == 1 || signs[j] == -1) {
          // printf("i,j = %d, %d\n", i, j);
          coeff1 = intvec_get_elem(my_vector_rotated_ptr, offsets[j]);
          //int_dump(coeff1);
          int_set(new, multiplier);
          int_mul_sgn_self(new, signs[j]);
          //int_dump(coeff2);
          int_add(coeff1, coeff1, new);
          // int_mod(coeff1, coeff1, ring->q);
          // int_redc(coeff1, coeff1, ring->q);
          //}
        }
        //intvec_dump(my_vector_rotated);

        coeff = intvec_get_elem(tmp_mat, i+1);
        coeff2 = intvec_get_elem(my_vector_rotated, i+1);
        int_set(coeff, coeff2);
        intvec_set(col, my_vector_rotated);
    }

    // intmat_get_col(col, tmp_mat, 0);
    // for (int i=0; i<10; i++)
    //   printf("%d ", intvec_get_elem_i64(col, i));
    // printf("\n\n");

    intvec_set(row, tmp_mat);
    printf("- Building rot row ended.\n");
}

void gbfv_rot_row_fast (intvec_t row, const intvec_t c1, unsigned int l, const polyring_t ring) {
    // printf("\n- Building GBFV row (2) %d ...\n", l);

    // prepare materials from polynomial modulus
    int delta = 1024;
    int ind = 0;
    int n = c1->nelems;

    int matrix[12][12] = {{-1, 1, -1, 0, 0, 0, 0, -1, 1, -1, 0, 0}, 
                          {-1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0}, 
                          {0, -1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1}, 
                          {1, -1, 0, 0, 0, -1, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, 
                          {0, -1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0}, 
                          {1, -1, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 1, -1, 1, 0, 0, 0}};

    // result vector
    intvec_t tmp_vec;
    intvec_alloc(tmp_vec, c1->nelems, ring->q->nlimbs);

    INT_T(tmp, ring->q->nlimbs);
    INT_T(tmp2, ring->q->nlimbs);
    int_ptr coeff, coeff2;

    int rk = l / delta;
    int startj, endj;

    for (int i=0; i<n; i++) {
        // printf("i = %d\n", i);
        int_set_zero(tmp);

        // get term in x^k from x^i*c1(x)
        if (l >= i) {
          coeff2 = intvec_get_elem(c1, l-i);
          int_set(tmp, coeff2);
        }

        startj = ceil((float)(n-l) / (float)delta);
        endj = floor((float)(n-l+i-1) / (float)delta);
        // if (i<=10)
        //   printf("limits %d, %d\n", startj, endj);

        for (int j=startj; j<=endj; j++) {
          if (matrix[rk][j-startj] != 0) {
            // printf("getting element %d\n", l+j*delta-i);
            coeff2 = intvec_get_elem(c1, l+j*delta-i);
            int_set(tmp2, coeff2);
            // printf("getting sign %d, %d\n", rk, j-startj);
            int_mul_sgn_self(tmp2, matrix[rk][j-startj]);
            int_add(tmp, tmp, tmp2);
            // if (i<=10)
            //   int_dump(tmp);
          }
        }

        coeff2 = intvec_get_elem(tmp_vec, i);
        int_set(coeff2, tmp);
    }
    
    intvec_set(row, tmp_vec);
    intvec_free(tmp_vec);
    // printf("- Building rot row ended.\n");
}

void gbfv_rot_sum (int_t sum, const intvec_t c1, unsigned int l, int8_t sk[], const polyring_t ring) {
    // printf("\n- Building GBFV row (2) %d ...\n", l);

    // prepare materials from polynomial modulus
    int delta = 1024;
    int ind = 0;
    int n = c1->nelems;

    int matrix[12][12] = {{-1, 1, -1, 0, 0, 0, 0, -1, 1, -1, 0, 0}, 
                          {-1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0}, 
                          {0, -1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1}, 
                          {1, -1, 0, 0, 0, -1, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, 
                          {0, -1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0}, 
                          {1, -1, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 1, -1, 1, 0, 0, 0}};

    // result vector
    INT_T(tmp, 2*ring->q->nlimbs);
    INT_T(tmp2, 2*ring->q->nlimbs);
    int_ptr coeff, coeff2;

    int rk = l / delta;
    int startj, endj;

    for (int i=0; i<n; i++) {
      if (sk[i] != 0) {
        // printf("i = %d\n", i);
        int_set_zero(tmp);

        // get term in x^k from x^i*c1(x)
        if (l >= i) {
          coeff2 = intvec_get_elem(c1, l-i);
          int_set(tmp, coeff2);
        }

        startj = ceil((float)(n-l) / (float)delta);
        endj = floor((float)(n-l+i-1) / (float)delta);
        // printf("limits %d, %d\n", startj, endj);

        for (int j=startj; j<=endj; j++) {
          if (matrix[rk][j-startj] != 0) {
            // printf("getting element %d\n", l+j*delta-i);
            coeff2 = intvec_get_elem(c1, l+j*delta-i);
            int_set(tmp2, coeff2);
            // printf("getting sign %d, %d\n", rk, j-startj);
            int_mul_sgn_self(tmp2, matrix[rk][j-startj]);
            int_add(tmp, tmp, tmp2);
          }
        }

        int_mul_sgn_self(tmp, sk[i]);
        int_add(sum, sum, tmp);

        // int_mod(sum, sum, ring->q);
        // int_redc(sum, sum, ring->q);
      }
    }

    // printf("- Building rot row ended.\n");
}


void gbfv_rot_sum_degree2 (int_t sum, const intvec_t c1, unsigned int l, int8_t sk[], const polyring_t ring) {
    // printf("\n- Building GBFV row (2) %d ...\n", l);

    // prepare materials from polynomial modulus
    int delta = 256;
    int ind = 0;
    int n = c1->nelems;

    int matrix[12][12] = {{-1, 1, -1, 0, 0, 0, 0, -1, 1, -1, 0, 0}, 
                          {-1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0}, 
                          {0, -1, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1}, 
                          {1, -1, 0, 0, 0, -1, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, 
                          {0, -1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0}, 
                          {1, -1, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0}, 
                          {1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}, 
                          {0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 
                          {-1, 1, 0, 0, 0, 0, 1, -1, 1, 0, 0, 0}};

    // result vector
    INT_T(tmp, 2*ring->q->nlimbs);
    INT_T(tmp2, 2*ring->q->nlimbs);
    int_ptr coeff, coeff2;

    int rk = l / delta;
    int startj, endj;

    for (int i=0; i<n; i++) {
      if (sk[i] != 0) {
        // printf("i = %d\n", i);
        int_set_zero(tmp);

        // get term in x^k from x^i*c1(x)
        if (l >= i) {
          coeff2 = intvec_get_elem(c1, l-i);
          int_set(tmp, coeff2);
        }

        startj = ceil((float)(n-l) / (float)delta);
        endj = floor((float)(n-l+i-1) / (float)delta);
        // printf("limits %d, %d\n", startj, endj);

        for (int j=startj; j<=endj; j++) {
          if (matrix[rk][j-startj] != 0) {
            // printf("getting element %d\n", l+j*delta-i);
            coeff2 = intvec_get_elem(c1, l+j*delta-i);
            int_set(tmp2, coeff2);
            // printf("getting sign %d, %d\n", rk, j-startj);
            int_mul_sgn_self(tmp2, matrix[rk][j-startj]);
            int_add(tmp, tmp, tmp2);
          }
        }

        int_mul_sgn_self(tmp, sk[i]);
        int_add(sum, sum, tmp);

        // int_mod(sum, sum, ring->q);
        // int_redc(sum, sum, ring->q);
      }
    }

    // printf("- Building rot row ended.\n");
}


/* expand i-th row of Rprime from cseed and 256 + i */
// static inline void
// _expand_Rprime_i (int8_t *Rprimei, unsigned int ncols, unsigned int i,
//                   const uint8_t cseed[32])
// {
//   _brandom (Rprimei, ncols, 1, cseed, 256 + i);
// }







/* swap row and col iff row > col */
static inline void
_diag (unsigned int *row, unsigned int *col, unsigned int r, unsigned int c)
{
  if (r > c)
    {
      *row = c;
      *col = r;
    }
  else
    {
      *row = r;
      *col = c;
    }
}

#define MAX(x, y) ((x) >= (y) ? (x) : (y))

/*
 * r = U^T*auto(a) = U*auto(a)
 * for each dim 2 subvec:
 * (a,b) -> auto((b,a))
 */
static void
__shuffleautovecsparse (spolyvec_t r)
{
  poly_ptr rp;
  unsigned int i, elem;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  _SVEC_FOREACH_ELEM (r, i)
  {
    rp = spolyvec_get_elem (r, i);
    elem = spolyvec_get_elem_ (r, i);

    poly_auto_self (rp);
    spolyvec_set_elem_ (r, i, elem % 2 == 0 ? elem + 1 : elem - 1);
  }

  r->sorted = 0; // XXX simpler sort possible
  spolyvec_sort (r);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}

/*
 *
 * r = U^T*auto(a)*U = U*auto(a)*U
 * for each 2x2 submat on or above the main diagonal:
 * [[a,b],[c,d]] -> auto([[d,c],[b,a]])
 * r != a
 */
static void
__shuffleauto2x2submatssparse (spolymat_t a)
{
  poly_ptr ap;
  unsigned int i, arow, acol;

  ASSERT_ERR (spolymat_get_nrows (a) % 2 == 0);
  ASSERT_ERR (spolymat_get_ncols (a) % 2 == 0);
  ASSERT_ERR (spolymat_is_upperdiag (a));

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  _SMAT_FOREACH_ELEM (a, i)
  {
    ap = spolymat_get_elem (a, i);
    arow = spolymat_get_row (a, i);
    acol = spolymat_get_col (a, i);

    if (arow % 2 == 0 && acol % 2 == 0)
      {
        spolymat_set_row (a, i, arow + 1);
        spolymat_set_col (a, i, acol + 1);
        poly_auto_self (ap);
      }
    else if (arow % 2 == 1 && acol % 2 == 1)
      {
        spolymat_set_row (a, i, arow - 1);
        spolymat_set_col (a, i, acol - 1);
        poly_auto_self (ap);
      }
    else if (arow % 2 == 1 && acol % 2 == 0)
      {
        spolymat_set_row (a, i, arow - 1);
        spolymat_set_col (a, i, acol + 1);
        poly_auto_self (ap);
      }
    else
      {
        /*
         * arow % 2 == 0 && acol % 2 == 1
         * This element's automorphism may land in the subdiagonal -1
         * if the 2x2 submat is on the main diagonal.
         * Check for this case and keep the matrix upper diagonal.
         */
        if (arow + 1 > acol - 1)
          {
            poly_auto_self (ap);
          }
        else
          {
            spolymat_set_row (a, i, arow + 1);
            spolymat_set_col (a, i, acol - 1);
            poly_auto_self (ap);
          }
      }
  }
  a->sorted = 0;
  spolymat_sort (a);
  ASSERT_ERR (spolymat_is_upperdiag (a));

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


static void
__schwartz_zippel_accumulate (spolymat_ptr R2i, spolyvec_ptr r1i, poly_ptr r0i,
                              spolymat_ptr Rprime2i[], spolyvec_ptr rprime1i[],
                              poly_ptr rprime0i[], unsigned int M_alt,
                              const intvec_t v,
                              const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  polyring_srcptr Rq = quad_eval->ring;
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  const unsigned int n = 2 * (m1 + l);
  spolyvec_t u0, u1, u2;
  spolymat_t t0, t1, t2;
  unsigned int j;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  spolyvec_alloc (u0, Rq, n, n);
  spolyvec_alloc (u1, Rq, n, n);
  spolyvec_alloc (u2, Rq, n, n);
  spolymat_alloc (t0, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t1, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t2, Rq, n, n, (n * n - n) / 2 + n);

  /* R2i */

  spolymat_set (t0, R2i);
  for (j = 0; j < M_alt; j++)
    {
      if (Rprime2i[j] == NULL)
        continue;

      spolymat_fromcrt (Rprime2i[j]);

      spolymat_scale (t1, intvec_get_elem (v, j), Rprime2i[j]);
      spolymat_add (t2, t0, t1, 0);
      spolymat_set (t0, t2);
    }
  spolymat_mod (R2i, t0);

  /* r1i */

  spolyvec_set (u0, r1i);
  for (j = 0; j < M_alt; j++)
    {
      if (rprime1i[j] == NULL)
        continue;

      spolyvec_fromcrt (rprime1i[j]);

      spolyvec_scale (u1, intvec_get_elem (v, j), rprime1i[j]);
      spolyvec_add (u2, u0, u1, 0);
      spolyvec_set (u0, u2);
    }
  spolyvec_mod (r1i, u0);

  if (r0i != NULL)
    {
      for (j = 0; j < M; j++)
        {
          if (rprime0i[j] == NULL)
            continue;

          poly_fromcrt (rprime0i[j]);
          poly_addscale (r0i, intvec_get_elem (v, j), rprime0i[j], 0);
        }
      poly_mod (r0i, r0i);
    }

  spolyvec_free (u0);
  spolyvec_free (u1);
  spolyvec_free (u2);
  spolymat_free (t0);
  spolymat_free (t1);
  spolymat_free (t2);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


/* just add equations that have already been multiplied by a challenge. */
static void
__schwartz_zippel_accumulate_ (spolymat_ptr R2i, spolyvec_ptr r1i,
                               poly_ptr r0i, spolymat_ptr Rprime2i[],
                               spolyvec_ptr rprime1i[], poly_ptr rprime0i[],
                               unsigned int M_alt,
                               const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  polyring_srcptr Rq = quad_eval->ring;
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  const unsigned int n = 2 * (m1 + l);
  spolyvec_t u0, u2;
  spolymat_t t0, t2;
  unsigned int j;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);
  // printf("  -   - inside accumulate function\n");
  // printf("  -   - m1=%d, l=%d, n=%d\n", m1, l, n);

  spolyvec_alloc (u0, Rq, n, n);
  spolyvec_alloc (u2, Rq, n, n);
  spolymat_alloc (t0, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t2, Rq, n, n, (n * n - n) / 2 + n);

  /* R2i */
  // printf("  -   - into R2i, M_alt = %d\n", M_alt);
  spolymat_set (t0, R2i);
  for (j = 0; j < M_alt; j++)
    {
      if (Rprime2i[j] == NULL)
        continue;

      spolymat_fromcrt (Rprime2i[j]);

      spolymat_add (t2, t0, Rprime2i[j], 0);
      spolymat_set (t0, t2);
    }
  spolymat_mod (R2i, t0);

  /* r1i */
  // printf("  -   - into r1i\n");
  spolyvec_set (u0, r1i);
  for (j = 0; j < M_alt; j++)
    {
      if (rprime1i[j] == NULL)
        continue;

      spolyvec_fromcrt (rprime1i[j]);

      spolyvec_add (u2, u0, rprime1i[j], 0);
      spolyvec_set (u0, u2);
    }
  spolyvec_mod (r1i, u0);

  // printf("  -   - into r0i\n");
  if (r0i != NULL)
    {
      for (j = 0; j < M; j++)
        {
          if (rprime0i[j] == NULL)
            continue;

          poly_fromcrt (rprime0i[j]);
          poly_add (r0i, r0i, rprime0i[j], 0);
        }
      poly_mod (r0i, r0i);
    }

  spolyvec_free (u0);
  spolyvec_free (u2);
  spolymat_free (t0);
  spolymat_free (t2);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


static void
__schwartz_zippel_auto (spolymat_ptr R2i, spolyvec_ptr r1i, poly_ptr r0i,
                        spolymat_ptr R2i2, spolyvec_ptr r1i2, poly_ptr r0i2,
                        const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  polyring_srcptr Rq = quad_eval->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int m1 = quad_eval->m1;
  const unsigned int l = quad_eval->l;
  const unsigned int n = 2 * (m1 + l);
  spolyvec_t u0, u1, u2;
  spolymat_t t0, t1, t2;
  poly_t tpoly;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  poly_alloc (tpoly, Rq);

  spolyvec_alloc (u0, Rq, n, n);
  spolyvec_alloc (u1, Rq, n, n);
  spolyvec_alloc (u2, Rq, n, n);
  spolymat_alloc (t0, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t1, Rq, n, n, (n * n - n) / 2 + n);
  spolymat_alloc (t2, Rq, n, n, (n * n - n) / 2 + n);

  /* R2i */

  spolymat_fromcrt (R2i);
  spolymat_fromcrt (R2i2);

  spolymat_set (t0, R2i);
  __shuffleauto2x2submatssparse (t0);

  spolymat_add (t1, R2i, t0, 0); // t1 = R2i + Uo(R2i)U

  spolymat_set (t0, R2i2);
  spolymat_lrot (t2, t0, d / 2);

  spolymat_add (R2i, t1, t2, 0); // R2i = R2i + Uo(R2i)U + R2i2X^(d/2)

  spolymat_set (t0, R2i2);
  __shuffleauto2x2submatssparse (t0);
  spolymat_lrot (t1, t0, d / 2);
  spolymat_add (t0, R2i, t1,
                0); // t0 = R2i + Uo(R2i)U + R2i2X^(d/2) +  Uo(R2i2)UX^(d/2)

  spolymat_scale (R2i, Rq->inv2, t0);

  /* r1i */

  spolyvec_fromcrt (r1i);
  spolyvec_fromcrt (r1i2);

  spolyvec_set (u0, r1i);
  __shuffleautovecsparse (u0);

  spolyvec_add (u1, r1i, u0, 0); // u1 = r1i + Uo(r1i)U

  spolyvec_set (u0, r1i2);
  spolyvec_lrot (u2, u0, d / 2);

  spolyvec_add (r1i, u1, u2, 0); // r1i = r1i + Uo(r1i)U + r1i2X^(d/2)

  spolyvec_set (u0, r1i2);
  __shuffleautovecsparse (u0);
  spolyvec_lrot (u1, u0, d / 2);
  spolyvec_add (u0, r1i, u1,
                0); // t0 = r1i + Uo(r1i)U + r1i2X^(d/2) +  Uo(r1i2)UX^(d/2)

  spolyvec_scale (r1i, Rq->inv2, u0);

  /* r0i */
  if (r0i != NULL)
    {
      poly_fromcrt (r0i);
      poly_fromcrt (r0i2);

      poly_auto (tpoly, r0i);
      poly_add (r0i, r0i, tpoly, 0); // r0i = r0i + o(r0i)

      poly_lrot (tpoly, r0i2, d / 2);
      poly_add (r0i, r0i, tpoly, 0); // r0i = r0i + o(r0i) + r0i2X^(d/2)

      poly_auto (tpoly, r0i2);
      poly_lrot (tpoly, tpoly, d / 2);
      poly_add (r0i, r0i, tpoly,
                0); // r0i = r0i + o(r0i) + r0i2X^(d/2) + o(r0i2)X^(d/2)

      poly_scale (r0i, Rq->inv2, r0i);
    }

  spolyvec_free (u0);
  spolyvec_free (u1);
  spolyvec_free (u2);
  spolymat_free (t0);
  spolymat_free (t1);
  spolymat_free (t2);
  poly_free (tpoly);

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


/*
 * R2i, r1i, r0i: first accumulator (lambda/2 eqs)
 * R2i2, r1i2, r0i2: second accumulator (lambda/2 eqs)
 * R2primei r1primei, r0primei: input eqs (M eqs)
 * Result is in first accumulator.
 */
static void
__schwartz_zippel_accumulate2 (spolymat_ptr R2i[], spolyvec_ptr r1i[],
                               poly_ptr r0i[], spolymat_ptr R2i2[],
                               spolyvec_ptr r1i2[], poly_ptr r0i2[],
                               spolymat_ptr R2primei[],
                               spolyvec_ptr r1primei[], poly_ptr r0primei[],
                               unsigned int M_alt, const uint8_t seed[32],
                               uint32_t dom,
                               const lnp_quad_eval_params_t params)
{
  abdlop_params_srcptr quad_eval = params->quad_eval;
  const unsigned int lambda = params->lambda;
  polyring_srcptr Rq = quad_eval->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  INTVEC_T (V, 2 * M_alt, Rq->q->nlimbs);
  intvec_t subv1, subv2;
  unsigned int i;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s begin", __func__);

  intvec_get_subvec (subv1, V, 0, M_alt, 1);
  intvec_get_subvec (subv2, V, M_alt, M_alt, 1);
  intvec_urandom (V, q, log2q, seed, dom);

  for (i = 0; i < lambda / 2; i++)
    {
      __schwartz_zippel_accumulate (R2i[i], r1i[i], r0i[i], R2primei, r1primei,
                                    r0primei, M_alt, subv1, params);
      __schwartz_zippel_accumulate (R2i2[i], r1i2[i], r0i2[i], R2primei,
                                    r1primei, r0primei, M_alt, subv2, params);
    }

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s end", __func__);
}


static void
__schwartz_zippel_accumulate_beta (spolymat_ptr R2i[], spolyvec_ptr r1i[],
                                    poly_ptr r0i[], spolymat_ptr R2i2[],
                                    spolyvec_ptr r1i2[], poly_ptr r0i2[],
                                    UNUSED spolymat_ptr R2t, spolyvec_ptr r1t,
                                    UNUSED poly_ptr r0t,
                                    const uint8_t seed[32], uint32_t dom,
                                    const lnp_quad_eval_params_t params,
                                    const unsigned int nprime)
{
  abdlop_params_srcptr tbox = params->quad_eval;
  const unsigned int Z = 0;
  polyring_srcptr Rq = tbox->ring;
  const unsigned int d = polyring_get_deg (Rq);
  const unsigned int m1 = tbox->m1 - Z;
  const unsigned int nex = 0;
  //const unsigned int l = tbox->l;
  const unsigned int l = 0;
  //const unsigned int nprime = params->nprime;
  const unsigned int loff
      = (nprime > 0 ? 256 / d : 0) + (nex > 0 ? 256 / d : 0);
  const unsigned int ibeta = (m1 + Z + l + loff) * 2;
  poly_ptr poly;
  int_ptr coeff;
  unsigned int i;
  spolymat_ptr R2tptr[1];
  spolyvec_ptr r1tptr[1];
  poly_ptr r0tptr[1];

  // printf("checking r1i[0] spolyvec inside acc beta 1. nelems_max = %d\n", r1i[0]->nelems_max);
  // unsigned int elem, j;
  // _SVEC_FOREACH_ELEM (r1i[0], j)
  // {
  //   elem = spolyvec_get_elem_ (r1i[0], j);
  //   if (elem >= r1i[0]->nelems_max)
  //     printf("current elem is %d / %d\n", elem, r1i[0]->nelems_max);
  // }
  // printf("ibeta is %d, loff is %d, l is %d, m1 is %d\n", ibeta, loff, l, m1);

  // d-1 eval eqs in beta,o(beta), for i=1,...,d-1:
  // prove const coeff of X^i * beta4 = 0 -> -i2*x^i*x^(d/2)*beta +
  // i2*x^i*x^(d/2)*o(beta) = 0 terms: R2: 0, r1: 2, r0: 0 | * (d-1)
  for (i = 1; i < d; i++)
    {
      R2tptr[0] = NULL;

      spolyvec_set_empty (r1t);
      poly = spolyvec_insert_elem (r1t, ibeta);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, i);
      int_set (coeff, Rq->inv2);
      poly = spolyvec_insert_elem (r1t, ibeta + 1);
      poly_set_zero (poly);
      coeff = poly_get_coeff (poly, i);
      int_set (coeff, Rq->inv2);
      // spolyvec_set_empty (r1t);
      // poly = spolyvec_insert_elem (r1t, ibeta);
      // poly_set_zero (poly);
      // if (i < d / 2)
      //   {
      //     coeff = poly_get_coeff (poly, i + d / 2);
      //     int_neg (coeff, Rq->inv2);
      //   }
      // else
      //   {
      //     coeff = poly_get_coeff (poly, i + d / 2 - d);
      //     int_set (coeff, Rq->inv2);
      //   }
      // poly = spolyvec_insert_elem (r1t, ibeta + 1);
      // poly_set_zero (poly);
      // if (i < d / 2)
      //   {
      //     coeff = poly_get_coeff (poly, i + d / 2);
      //     int_set (coeff, Rq->inv2);
      //   }
      // else
      //   {
      //     coeff = poly_get_coeff (poly, i + d / 2 - d);
      //     int_neg (coeff, Rq->inv2);
      //   }
      r1t->sorted = 1;
      r1tptr[0] = r1t;

      r0tptr[0] = NULL;

      __schwartz_zippel_accumulate2 (R2i, r1i, r0i, R2i2, r1i2, r0i2, R2tptr,
                                     r1tptr, r0tptr, 1, seed, dom + i,
                                     params);
    }
}


static void
__schwartz_zippel_accumulate_z (spolymat_ptr R2i[], spolyvec_ptr r1i[],
                                 poly_ptr r0i[], spolymat_ptr R2i2[],
                                 spolyvec_ptr r1i2[], poly_ptr r0i2[],
                                 spolymat_ptr R2t, spolyvec_ptr r1t,
                                 poly_ptr r0t, 
                                 intvec_t u_, 
                                 polyvec_t z4, intvec_t ct1, const uint8_t seed[32],
                                 uint32_t dom, const lnp_quad_eval_params_t params,
                                 const unsigned int nprime)
{
  abdlop_params_srcptr tbox = params->quad_eval;
  polyring_srcptr Rq = tbox->ring;
  int_srcptr q = Rq->q;
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int d = Rq->d;
  const unsigned int Z = 0;
  const unsigned int m1 = tbox->m1 - Z;
  //const unsigned int l = tbox->l;
  const unsigned int l = 0;
  const unsigned int nex = 0;
  //const unsigned int nprime = params->nprime;
  const unsigned int loff3 = (nex > 0 ? 256 / d : 0);
  const unsigned int loff4 = (nprime > 0 ? 256 / d : 0);
  const unsigned int loff = loff3 + loff4;
  const unsigned int ibeta = (m1 + Z + l + loff) * 2; // XXX correct
  const unsigned int is1 = 0;
  const unsigned int im = (m1 + Z) * 2;
  const unsigned int iy4 = (m1 + Z + l + loff3) * 2; // XXX correct
  unsigned int i, j, k;
  int8_t Rprimei[nprime * d];
  const unsigned int lambda = params->lambda;
  spolymat_ptr R2tptr[1];
  spolyvec_ptr r1tptr[1];
  poly_ptr r0tptr[1];
  int_ptr chal, acc;
  polymat_t mat, vRDs, vRDm; // vRpol;
  // polyvec_t subv1, subv2;
  int_ptr coeff1, coeff2;
  poly_ptr poly, poly2, poly3;
  int_srcptr inv2 = Rq->inv2;
  intvec_t row1;

  INT_T (tmp, 2 * Rq->q->nlimbs);
  //INTVEC_T (u_, nprime * d, Rq->q->nlimbs);
  INTVEC_T (z4_, 256, Rq->q->nlimbs);
  INTMAT_T (V, lambda, 256, Rq->q->nlimbs);
  // INTMAT_T (vR_, lambda, nprime * d, Rq->q->nlimbs);
  // INTMAT_T (vR, lambda, nprime * d, 2 * Rq->q->nlimbs);
  INTVEC_T (vRu, lambda, 2 * Rq->q->nlimbs);
  intmat_urandom (V, q, log2q, seed, dom);

  intmat_t vR_, vR;
  intmat_alloc(vR_, lambda, nprime * d, Rq->q->nlimbs);
  intmat_alloc(vR, lambda, nprime * d, 2*Rq->q->nlimbs);

  // for (i=0; i<V->nrows; i++) {
  //   for (j=0; j<V->ncols; j++) {
  //     coeff1 = intmat_get_elem(V, i, j);
  //     int_set_one(coeff1);
  //   }
  // }
  // intmat_dump(V);

  // polymat_alloc (mat, Rq, nprime, MAX (m1, l));
  //polymat_alloc (vRpol, Rq, lambda, nprime);
  polymat_alloc (vRDs, Rq, lambda, m1);
  if (l > 0)
    polymat_alloc (vRDm, Rq, lambda, l);

  // eval eqs in s1,o(s1),m,o(m),y4,o(y4),beta,o(beta) (approx. range proof
  // inf) prove z4 = y4 + beta4*Rprime*s4 over int vec of dim 256
  //       y4 - z4 + beta4*Rprime*s4 = 0
  //       y4 - z4 + beta4*Rprime*(Ds*s1+Dm*m+u) = 0  # view Ds resp Dm as int
  //       rotation matrices in Z^(d*nprimexd*m1) resp Z^(d*nprimexd*l)
  //   (Ds*s1+Dm*m+u is small, but Ds*s1, Dm*m, u do not need to be ..., so the
  //   below holds mod q)
  //       y4 - z4 + beta4*(Rprime*Ds)*s1 + beta4*(Rprime*Dm)*m +
  //       beta4*(Rprime*u) = 0
  // for i in 0,...,255
  //       y4i - z4i + beta4*(Rprime*Ds)i*s1 + beta4*(Rprime*Dm)i*m +
  //       beta4*(Rprime*u)i = 0
  // terms: R2: 2m1+2l, r1: 3, r0: 1 | * 256

  // compute vR=v*Rprime
  // then vR*Ds, vR*Dm, vR*u
  // printf("  - m1=%d, l=%d\n", m1, l);
  // printf("  - is1=%d, im=%d, iy4=%d, ibeta=%d\n", is1, im, iy4, ibeta);

  // printf("start accumulating z_v\n");

  // printf("  - building z4_\n");
  // instantiates z4_ intvec of coefficients of z4
  for (i = 0; i < loff4; i++)
    {
      poly = polyvec_get_elem (z4, i);
      for (j = 0; j < d; j++)
        {
          coeff1 = poly_get_coeff (poly, j);
          coeff2 = intvec_get_elem (z4_, i * d + j);
          int_set (coeff2, coeff1);
        }
    }
  //printf("\n");


  intmat_set_zero (vR);
  // vR is lambda * nprime*d matrix where challenge k out of lambda multiplies line k of matrix R
  // printf("  - building vR\n");
  // #region building vR + vR_
  for (i = 0; i < 256; i++)
    {
      //printf("row %d\n", i);
      _expand_R_i2 (Rprimei, nprime * d, i, seed);

      for (k = 0; k < lambda; k++)
        {
          chal = intmat_get_elem (V, k, i);

          for (j = 0; j < nprime * d; j++)
            {
              if (Rprimei[j] == 0)
                {
                }
              else
                {
                  ASSERT_ERR (Rprimei[j] == 1 || Rprimei[j] == -1);

                  acc = intmat_get_elem (vR, k, j);

                  int_set (tmp, chal);
                  int_mul_sgn_self (tmp, Rprimei[j]);
                  int_add (acc, acc, tmp);
                }
            }
        }
    }
  // printf("finished building vR\n");



  printf("  - computing vR_\n");
  // vR_ is same as vR but with correct number of limbs (after mod q)
  _MAT_FOREACH_ELEM (vR, i, j)
  {
    coeff1 = intmat_get_elem (vR, i, j);
    coeff2 = intmat_get_elem (vR_, i, j); // XXX correct
    int_mod (coeff2, coeff1, q);
  }
  // #endregion

  // #region vRu
  printf("  - computing vRu, vR_cols=%d, u_rows=%d\n", vR_->ncols, u_->nelems);
  // generates u_, intvec with coefficients of elements in u.
  // vRu is intvec of lambda entries, vRu = vR_ * u_
  if (u_ != NULL)
    {
      for (k = 0; k < lambda; k++)
        {
          intmat_get_row (row1, vR_, k);
          coeff1 = intvec_get_elem (vRu, k); // XXX correct
          intvec_dot2 (coeff1, row1, u_, Rq);
        }
    }
#if 0
  int32_t RPRIME[nprime * d * 256];
  INTMAT_T (Rprime, 256, nprime * d, 1);
  for (i = 0; i < 256; i++)
    {
      _expand_Rprime_i (Rprimei, nprime * d, i, seed);
      for (j = 0; j < nprime * d; j++)
        RPRIME[i * (nprime * d) + j] = Rprimei[j];
    }
  intmat_set_i32 (Rprime, RPRIME);
  intmat_dump (Rprime);
  //intmat_dump (V);
  //intmat_dump (vR);
#endif
  // #endregion


  // #region vRDs - current version
  // building vRDs and vRDm 
  // consists of the ring elements from multiplying the integers in vR_ and Ds/Dm
  printf("  - building vRDs\n");
  // printf("  - m1*d %d\n", m1*d);
  polymat_t ovRDs;
  polymat_alloc(ovRDs, Rq, vRDs->nrows, vRDs->ncols);
  // polymat_auto (ovRDs, vRDs);

  intmat_t vRDs_coeffs;
  intmat_alloc(vRDs_coeffs, vRDs->nrows, d * vRDs->ncols, Rq->q->nlimbs);
  
  struct timeval start, end;
  long seconds, useconds;
  double wall_time;
  gettimeofday(&start, NULL);  // Start timing


  // for (k = 0; k < lambda; k++) {
  //     printf("lambda = %d\n", k);
  //     intvec_t row11;
  //     intmat_get_row(row11, vR_, k);

  //     #pragma omp parallel for private(j) shared(vRDs, Ds, row11)
  //     for (i = 0; i < Ds->ncols; i++) {
  //       // printf("Hello world!! from: %d, value:%d\n", omp_get_thread_num(), k);
      
  //       poly_ptr polyy = polymat_get_elem (vRDs, k, i);
  //       poly_set_zero(polyy);
  //       for (j = 0; j < Ds->nrows; j++) {
  //         poly_ptr poly22 = polymat_get_elem (Ds, j, i);
  //         int_ptr coeff11 = intvec_get_elem(row11, j);
  //         poly_addscale(polyy, coeff11, poly22, 0);
  //       }
  //     }
  // }
  // poly = polymat_get_elem(vRDs, 0, 0);
  // poly_dump(poly);

  

  int_ptr coeff_vRDs;
  
  // #pragma omp parallel for private(j,k, i, coeff_vRDs) shared(vRDs_coeffs, vR_)
  for (i=0; i< m1*d; i++) {

    for (k=0; k<CT_COUNT; k++) {
      intvec_t subvec, ct_vec;
      intvec_get_subvec(subvec, ct1, k*m1*d, m1*d, 1);
      intvec_alloc(ct_vec, d * m1, Rq->q->nlimbs);
      intvec_ptr ct = &ct_vec; 
      gbfv_rot_col(ct_vec, subvec, i, Rq);

      for (j=0; j<lambda; j++) {
        INTVEC_T(vR_row_vec, d * m1, Rq->q->nlimbs);
        intvec_ptr vR_row = &vR_row_vec;
        intmat_get_row(vR_row, vR_, j);
        intvec_get_subvec(subvec, vR_row, k*m1*d, m1*d, 1);

        int_t new, new_;
        int_alloc(new, 2 * Rq->q->nlimbs);
        int_alloc(new_, Rq->q->nlimbs);
        // INT_T (new, 2 * Rq->q->nlimbs);
        intvec_dot2(new, subvec, ct_vec, Rq);
        int_mod(new_, new, Rq->q);

        coeff_vRDs = intmat_get_elem(vRDs_coeffs, j, i);
        int_add(coeff_vRDs, coeff_vRDs, new_);
        // int_mod(coeff_vRDs, coeff_vRDs, Rq->q);

        int_free(new);
      }
      intvec_free(ct_vec);
    }
  }
  
  // for (k=0; k<CT_COUNT; k++) {
  //   printf("vRDs current ciphertext %d\n", k);
  //   // getting k-th ct1 coeffs
  //   INTVEC_T(ct1_coeffs_vec, d * m1, Rq->q->nlimbs);
  //   intvec_ptr ct1_coeffs = &ct1_coeffs_vec;
  //   for (i=0; i< m1; i++) {
  //       // printf("%d\n", i);
  //       poly = polyvec_get_elem(ct1, k*m1 + i);     
  //       coeffs = poly_get_coeffvec (poly);
  //       for (j=0; j<d; j++) {
  //           intvec_set_elem(ct1_coeffs, i*d+j, intvec_get_elem(coeffs, j));
  //       }
  //   }
  //   // printf("going to multiply\n");
  //   // #pragma omp parallel for private(i,j,row1,coeff1) shared(vRDs_coeffs, vR_, ct1_coeffs)
  //   for (i=0; i<m1*d; i++) {
  //     // printf("%d\n", i);
  //     intvec_t row11;
  //     int_ptr coeff11;
  //     INTVEC_T(row, d * m1, Rq->q->nlimbs);
  //     INTVEC_T(scaled_row, d * m1, 2 * Rq->q->nlimbs);
  //     INTVEC_T(scaled_row2, d * m1, Rq->q->nlimbs);
  //     gbfv_rot_row_fast(row, ct1_coeffs, i, Rq);
  //     for (j=0; j<lambda; j++) {
  //       // printf("lambda %d\n", j);
  //       intmat_get_row(row11, vRDs_coeffs, j);
  //       coeff11 = intmat_get_elem(vR_, j, i + k*(d*m1));
  //       intvec_scale(scaled_row, coeff11, row);
  //       intvec_mod(scaled_row2, scaled_row, q);

  //       //#pragma omp critical
  //       intvec_add(row11, row11, scaled_row2);
        
  //       if (i%10 == 0)
  //         intvec_mod(row11, row11, q);
  //       // intvec_redc(row11, row11, q);
  //       // if (i == 0) {
  //       //   for (int j2=0; j2<10; j2++)
  //       //     printf("%lld ", intvec_get_elem_i64(row11, j2));
  //       //   printf("\n\n\n\n");
  //       // }
  //     }
  //   }
  // }
  intvec_t subvec;
  printf("putting coeffs into vRDs\n");
  for (i=0; i<vRDs_coeffs->nrows; i++) {
    intmat_get_row(row1, vRDs_coeffs, i);
    for (j=0; j<(row1->nelems)/d; j++) {
        intvec_get_subvec (subvec, row1, 0+j*d, d, 1); // XXX correct
        poly = polymat_get_elem(vRDs, i, j);
        poly_set_coeffvec(poly, subvec);
    }
  }
  
  printf("  - computing o(vRDs)\n");
  polymat_auto (ovRDs, vRDs);
  // polymat_lrot (vRDs, vRDs, d / 2); // * X^(d/2)  XXX correct


  gettimeofday(&end, NULL);  // End timing
  // Compute the time difference in microseconds
  seconds  = end.tv_sec  - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;
  wall_time = seconds + useconds/1e6;  // Convert to seconds
  printf(" ---------> (lambda & ovRDs): execution time: %f seconds\n", wall_time);


  // printf("  - building vRDm (old version)\n");
  // if (l > 0 && Dm != NULL) {
  //   for (k = 0; k < lambda; k++) {
  //       intmat_get_row(row1, vR_, k);

  //       for (i = 0; i < Dm->ncols; i++) {
  //           poly = polymat_get_elem (vRDm, k, i);
  //           poly_set_zero(poly);

  //           for (j = 0; j < Dm->nrows; j++) {
  //             poly2 = polymat_get_elem (Dm, j, i);
  //             coeff1 = intmat_get_elem(vR_, k, j);
  //             poly_addscale(poly, coeff1, poly2, 0);
  //           }
  //       }
  //   }
  //   printf("  - computing o(vRDm)\n");
  //   polymat_auto (vRDm, vRDm);
  //   // polymat_lrot (vRDm, vRDm, d / 2); // * X^(d/2)  XXX correct
  // }
  // #endregion


  // use previously built matrices to compute R2t, r1t and r0t
  for (k = 0; k < lambda; k++)
    {

      spolymat_set_empty (R2t);
      spolyvec_set_empty (r1t);
      poly_set_zero (r0t);

      R2tptr[0] = R2t;

      // get elements that multiply with s1 and with beta and o(beta)
      // set (s1,ibeta) element to be 1/2 * vRDs
      // set (s1,ibeta+1) element to be 1/2 * vRDs
      for (i = 0; i < m1; i++)
        {
          poly = spolymat_insert_elem (R2t, is1 + 2 * i, ibeta);
          poly2 = spolymat_insert_elem (R2t, is1 + 2 * i, ibeta + 1);

          poly3 = polymat_get_elem (ovRDs, k, i);
          poly_set (poly2, poly3);

          poly_scale (poly2, inv2, poly2);
          poly_set (poly, poly2);
        }


      R2t->sorted = 1;

      r1tptr[0] = r1t;

      // go to the elements of r1t that will multiply o(y4) and multiply the coefficients with challenges v
      for (i = 0; i < loff4; i++)
        {
          poly = spolyvec_insert_elem (r1t, iy4 + 1 + 2 * i);
          for (j = 0; j < d; j++)
            {
              coeff1 = poly_get_coeff (poly, j);
              coeff2 = intmat_get_elem (V, k, i * d + j);
              int_set (coeff1, coeff2);
              int_redc (coeff1, coeff1, q);
            }
        }

      // set ibeta  coefficient d/2 to be -1/2 * vRu
      // set ibeta+1 coefficient d/2 to be 1/2 * vRu 
      if (u_ != NULL)
        {
          poly = spolyvec_insert_elem (r1t, ibeta);
          poly2 = spolyvec_insert_elem (r1t, ibeta + 1);

          poly_set_zero (poly2);
          coeff1 = poly_get_coeff (poly2, 0);
          coeff2 = intvec_get_elem (vRu, k);
          int_mod (coeff1, coeff2, q);
          int_mul (tmp, inv2, coeff1);
          int_mod (coeff1, tmp, q);
          int_redc (coeff1, coeff1, q);

          poly_set_zero (poly);
          coeff2 = poly_get_coeff (poly, 0);
          //int_neg (coeff2, coeff1);
          int_set (coeff2, coeff1); // instead of line above
          int_redc (coeff2, coeff2, q);
        }

      r1t->sorted = 1;

      r0tptr[0] = r0t;

      poly_set_zero (r0t); // correct

      intmat_get_row (row1, V, k);
      intvec_dot2 (tmp, z4_, row1, Rq);
      // printf("dump tmp ");
      // int_mod(tmp, tmp, q);
      // int_redc(tmp, tmp, q);
      // int_dump(tmp);
      coeff1 = poly_get_coeff (r0t, 0);
      int_mod (coeff1, tmp, q);
      int_neg_self (coeff1);
      int_redc (coeff1, coeff1, q);

      // printf("  - calling accumulate functions\n");
      if (k % 2 == 0)
        __schwartz_zippel_accumulate_ (R2i[k / 2], r1i[k / 2], r0i[k / 2],
                                       R2tptr, r1tptr, r0tptr, 1,
                                       params);
      else
        __schwartz_zippel_accumulate_ (R2i2[k / 2], r1i2[k / 2], r0i2[k / 2],
                                       R2tptr, r1tptr, r0tptr, 1,
                                       params);
    }



  polymat_free (vRDs);
  if (l > 0)
    polymat_free (vRDm);
  //polymat_free (vRpol);
  // polymat_free (mat);
}


static void print_exec_time(struct timeval start, struct timeval end, const char *str) {
  long seconds  = end.tv_sec  - start.tv_sec;
  long useconds = end.tv_usec - start.tv_usec;
  double wall_time = seconds + useconds/1e6;  // Convert to seconds
  printf(" ---------> (%s): execution time: %f seconds\n", str, wall_time);
}