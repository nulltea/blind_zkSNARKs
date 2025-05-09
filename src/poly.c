#include "poly.h"
#include "brandom.h"
#include "grandom.h"
#include "intvec.h"
#include "lazer.h"
#include "memory.h"
// #include "ntt.h"
#include "urandom.h"
#include <math.h>

void
poly_alloc (poly_ptr r, const polyring_t ring)
{
  void *mem;

  mem = _alloc (_sizeof_poly_data (ring));

  _poly_init (r, ring, mem);
  r->mem = mem;
}

void
poly_free (poly_ptr r)
{
  if (r == NULL)
    return;

  if (r->crt)
    _free (r->crtrep, _sizeof_crtrep_data (r->ring));
  _free (r->mem, _sizeof_poly_data (r->ring));
}

void
poly_set (poly_t r, const poly_t a)
{
  polyring_srcptr ring = a->ring;
  intvec_srcptr acoeffs;

  // XXXASSERT_ERR (r->ring == a->ring);

  r->crt = a->crt;
  if (r->crt)
    {
      if (r->crtrep == NULL)
        r->crtrep = _alloc (_sizeof_crtrep_data (ring));

      memcpy (r->crtrep, a->crtrep,
              sizeof (crtcoeff_t) * ring->d * ring->nmoduli);
    }
  else
    {
      acoeffs = _get_coeffvec_src (a);
      poly_set_coeffvec (r, acoeffs);
    }
}

void
poly_add (poly_t r, poly_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  const unsigned int deg = polyring_get_deg (ring);
  crtcoeff_t *ca, *cb, *cr;
  unsigned int i;

  ASSERT_ERR (ring == poly_get_ring (a));
  ASSERT_ERR (ring == poly_get_ring (b));

  r->crt = _to_same_dom (a, b, crt);
  if (r->crt)
    {
      if (r->crtrep == NULL)
        r->crtrep = _alloc (_sizeof_crtrep_data (ring));

      _POLYRING_FOREACH_P (ring, i)
      {
#ifdef XXX
        crtcoeff_t tmp;
        unsigned int j;

        for (j = 0; j < ring->d; j++)
          {
            ca = _get_crtcoeff (a->crtrep, i, j, deg);
            cb = _get_crtcoeff (b->crtrep, i, j, deg);
            cr = _get_crtcoeff (r->crtrep, i, j, deg);

            tmp = (crtcoeff_t)(*ca) + (*cb);
            *cr = _mont_rsum (tmp, ring->moduli[i]->p);
          }
#else
        ca = _get_crtcoeff (a->crtrep, i, 0, deg);
        cb = _get_crtcoeff (b->crtrep, i, 0, deg);
        cr = _get_crtcoeff (r->crtrep, i, 0, deg);
        hexl_ntt_add (cr, ca, cb, deg, ring->moduli[i]->p);
#endif
      }
    }
  else
    {
      intvec_add (r->coeffs, a->coeffs, b->coeffs);
      intvec_redc (r->coeffs, r->coeffs, ring->q);
    }
}

void
poly_sub (poly_t r, poly_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  const unsigned int deg = polyring_get_deg (ring);
  crtcoeff_t *ca, *cb, *cr;
  unsigned int i;

  ASSERT_ERR (ring == poly_get_ring (a));
  ASSERT_ERR (ring == poly_get_ring (b));

  r->crt = _to_same_dom (a, b, crt);
  if (r->crt)
    {
      if (r->crtrep == NULL)
        r->crtrep = _alloc (_sizeof_crtrep_data (ring));

      _POLYRING_FOREACH_P (ring, i)
      {
#ifdef XXX
        crtcoeff_t tmp;
        unsigned int j;

        for (j = 0; j < ring->d; j++)
          {
            ca = _get_crtcoeff (a->crtrep, i, j, deg);
            cb = _get_crtcoeff (b->crtrep, i, j, deg);
            cr = _get_crtcoeff (r->crtrep, i, j, deg);

            tmp = (crtcoeff_t)(*ca) - (*cb);
            *cr = _mont_rsum (tmp, ring->moduli[i]->p);
          }
#else
        ca = _get_crtcoeff (a->crtrep, i, 0, deg);
        cb = _get_crtcoeff (b->crtrep, i, 0, deg);
        cr = _get_crtcoeff (r->crtrep, i, 0, deg);
        hexl_ntt_sub (cr, ca, cb, deg, ring->moduli[i]->p);
#endif
      }
    }
  else
    {
      intvec_sub (r->coeffs, a->coeffs, b->coeffs);
      intvec_redc (r->coeffs, r->coeffs, ring->q);
    }
}

void
poly_neg (poly_t r, poly_t b)
{
#if ASSERT == ASSERT_ENABLED
  polyring_srcptr ring = poly_get_ring (r);
#endif
  ASSERT_ERR (ring == poly_get_ring (b));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (b->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt neg");
#endif
  _fromcrt (b); // XXX can also do this in crt dom

  r->crt = 0;
  intvec_neg (r->coeffs, b->coeffs);
}

void
poly_scale (poly_t r, const int_t a, poly_t b)
{
  unsigned int i;
  polyring_srcptr ring = poly_get_ring (r);
  crtcoeff_t *cb, *cr;
  crtcoeff_t _a;

  ASSERT_ERR (ring == poly_get_ring (b));
  ASSERT_ERR (a->nlimbs == b->coeffs->nlimbs);

#if 0
  if (b->crt == 1)
    {
      const unsigned int deg = polyring_get_deg (ring);

      if (r->crtrep == NULL)
        {
          ASSERT_ERR (r->crt == 0);
          r->crtrep = _alloc (_sizeof_crtrep_data (ring));
        }
      r->crt = 1;

      _POLYRING_FOREACH_P (ring, i)
      {
        _a = int_mod_XXX_hexl (a, ring->moduli[i]->p);

        cb = _get_crtcoeff (b->crtrep, i, 0, deg);
        cr = _get_crtcoeff (r->crtrep, i, 0, deg);
        hexl_ntt_scale (cr, _a, cb, deg, ring->moduli[i]->p, 1);
      }
    }
  else
    {
#endif
  INTVEC_T (tmp, polyring_get_deg (ring), 2 * a->nlimbs);

  _fromcrt (b); // XXX can also do this in crt dom

  r->crt = 0;
  intvec_scale (tmp, a, _get_coeffvec_src (b));
  intvec_mod (r->coeffs, tmp, ring->q);
  intvec_redc (r->coeffs, r->coeffs, ring->q);
  //}
}

void
poly_rshift (poly_t r, poly_t a, unsigned int n)
{
  ASSERT_ERR (poly_get_ring (r) == poly_get_ring (a));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt rshift");
#endif
  _fromcrt (a);

  r->crt = 0;
  intvec_rshift (r->coeffs, a->coeffs, n);
}

void
poly_lshift (poly_t r, poly_t a, unsigned int n)
{
  ASSERT_ERR (poly_get_ring (r) == poly_get_ring (a));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt lshift");
#endif
  _fromcrt (a);

  r->crt = 0;
  intvec_lshift (r->coeffs, a->coeffs, n);
}

void
poly_rrot (poly_t r, poly_t a, unsigned int n)
{
  ASSERT_ERR (poly_get_ring (r) == poly_get_ring (a));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt rrot");
#endif
  _fromcrt (a);

  r->crt = 0;
  intvec_rrot (r->coeffs, a->coeffs, n);
}

void
poly_lrot (poly_t r, poly_t a, unsigned int n)
{
  ASSERT_ERR (poly_get_ring (r) == poly_get_ring (a));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt lrot");
#endif
  _fromcrt (a);

  r->crt = 0;
  intvec_lrot (r->coeffs, a->coeffs, n);
}

void
poly_mod (poly_t r, poly_t a)
{
  polyring_srcptr ring = poly_get_ring (r);

  ASSERT_ERR (poly_get_ring (r) == poly_get_ring (a));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt mod");
#endif
  _fromcrt (a);

  r->crt = 0;
  intvec_mod (r->coeffs, a->coeffs, ring->q);
}

void
poly_redc (poly_t r, poly_t a)
{
  polyring_srcptr ring = poly_get_ring (r);

  ASSERT_ERR (poly_get_ring (r) == poly_get_ring (a));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt redc");
#endif
  _fromcrt (a);

  r->crt = 0;
  intvec_redc (r->coeffs, a->coeffs, ring->q);
}

void
poly_redp (poly_t r, poly_t a)
{
  polyring_srcptr ring = poly_get_ring (r);

  ASSERT_ERR (poly_get_ring (r) == poly_get_ring (a));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt redp");
#endif
  _fromcrt (a);

  r->crt = 0;
  intvec_redp (r->coeffs, a->coeffs, ring->q);
}

int
poly_eq (poly_t a, poly_t b)
{
  ASSERT_ERR (poly_get_ring (a) == poly_get_ring (b));

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt eq");
#endif
  _fromcrt (a);
#if DEBUGINFO == DEBUGINFO_ENABLED
  if (b->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt eq");
#endif
  _fromcrt (b);

  return intvec_eq (_get_coeffvec_src (a), _get_coeffvec_src (b));
}

void
poly_toisoring (polyvec_t vec, poly_t a)
{
  polyring_srcptr Rprime = poly_get_ring (a);
  polyring_srcptr R = polyvec_get_ring (vec);
  const unsigned int dprime = polyring_get_deg (Rprime);
  const unsigned int d = polyring_get_deg (R);
  const unsigned int k = dprime / d;
  unsigned int i, j, l;
  poly_ptr poly;
  int_ptr coeffprime, coeff;

  ASSERT_ERR (!a->crt);
  ASSERT_ERR (d * k == dprime);
  ASSERT_ERR (polyvec_get_nelems (vec) == k);
  // XXXASSERT_ERR (int_eq (polyring_get_mod (Rprime), polyring_get_mod (R)) ==
  // 1);

  for (i = 0; i < k; i++)
    {
      poly = polyvec_get_elem (vec, i);
      poly->crt = 0;
      for (j = 0; j < d; j++)
        {
          l = k * j + i;

          coeffprime = poly_get_coeff (a, l);
          coeff = poly_get_coeff (poly, j);
          int_set (coeff, coeffprime);
        }
    }
}

void
poly_fromisoring (poly_t a, polyvec_t vec)
{
  polyring_srcptr Rprime = poly_get_ring (a);
  polyring_srcptr R = polyvec_get_ring (vec);
  const unsigned int dprime = polyring_get_deg (Rprime);
  const unsigned int d = polyring_get_deg (R);
  const unsigned int k = dprime / d;
  unsigned int i, j, l;
  poly_ptr poly;
  int_ptr coeffprime, coeff;

  ASSERT_ERR (d * k == dprime);
  ASSERT_ERR (polyvec_get_nelems (vec) == k);
  // XXXASSERT_ERR (int_eq (polyring_get_mod (Rprime), polyring_get_mod (R)) ==
  // 1);

  polyvec_fromcrt (vec);

  for (i = 0; i < k; i++)
    {
      poly = polyvec_get_elem (vec, i);
      for (j = 0; j < d; j++)
        {
          l = k * j + i;

          coeff = poly_get_coeff (poly, j);
          coeffprime = poly_get_coeff (a, l);
          int_set (coeffprime, coeff);
        }
    }
  a->crt = 0;
}

void
poly_mul (poly_t r, poly_t a, poly_t b)
{
  polyring_srcptr ring = r->ring;
  const unsigned int deg = polyring_get_deg (ring);
  unsigned int i;
  crtcoeff_t *ca, *cb, *cr;

  ASSERT_ERR (r->ring == a->ring);
  ASSERT_ERR (r->ring == b->ring);

  _tocrt (a);
  _tocrt (b);

  if (r->crtrep == NULL)
    {
      ASSERT_ERR (r->crt == 0);
      r->crtrep = _alloc (_sizeof_crtrep_data (ring));
    }
  r->crt = 1;

  _POLYRING_FOREACH_P (ring, i)
  {
#ifdef XXX
    crtcoeff_t p = ring->moduli[i]->p;
    crtcoeff_t pinv = ring->moduli[i]->mont_pinv;
    unsigned int j;
    crtcoeff_dbl_t tmp;

    for (j = 0; j < deg; j++)
      {
        ca = _get_crtcoeff (a->crtrep, i, j, deg);
        cb = _get_crtcoeff (b->crtrep, i, j, deg);
        cr = _get_crtcoeff (r->crtrep, i, j, deg);

        tmp = (crtcoeff_dbl_t)(*ca) * (*cb);
        *cr = _mont_rprod (tmp, p, pinv);
      }
#else
    ca = _get_crtcoeff (a->crtrep, i, 0, deg);
    cb = _get_crtcoeff (b->crtrep, i, 0, deg);
    cr = _get_crtcoeff (r->crtrep, i, 0, deg);
    hexl_ntt_mul (cr, ca, cb, 1, deg, ring->moduli[i]->p);
#endif
  }
}

void
poly_addmul (poly_t r, poly_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  poly_mul (tmp, a, b);
  poly_add (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_submul (poly_t r, poly_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  poly_mul (tmp, a, b);
  poly_sub (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_addmul2 (poly_t r, polymat_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, 1);

  ASSERT_ERR (polymat_get_nrows (a) == 1);

  polyvec_mul (tmp, a, b);
  poly_add (r, r, polyvec_get_elem (tmp, 0), crt);

  polyvec_free (tmp);
}

void
poly_submul2 (poly_t r, polymat_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  polyvec_t tmp;

  polyvec_alloc (tmp, ring, 1);

  ASSERT_ERR (polymat_get_nrows (a) == 1);

  polyvec_mul (tmp, a, b);
  poly_sub (r, r, polyvec_get_elem (tmp, 0), crt);

  polyvec_free (tmp);
}

void
poly_adddot (poly_t r, polyvec_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  polyvec_dot (tmp, a, b);
  poly_add (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_adddot2 (poly_t r, spolyvec_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  polyvec_dot2 (tmp, a, b);
  poly_add (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_subdot (poly_t r, polyvec_t a, polyvec_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  polyvec_dot (tmp, a, b);
  poly_sub (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_addscale (poly_t r, int_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  poly_scale (tmp, a, b);
  poly_add (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_subscale (poly_t r, int_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  poly_scale (tmp, a, b);
  poly_sub (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_tocrt (poly_t r)
{
  if (r->crt == 1)
    return;

  polyring_srcptr ring = poly_get_ring (r);
  const unsigned int deg = polyring_get_deg (ring);
  intvec_srcptr coeffs = r->coeffs;
  crtcoeff_t *crtrep;
  unsigned int i;
  void **hexl_ntts;
  // intvec_t tmp1;
  // INT_T (p, 1);
  crtcoeff_t *crtcoeff;

  DEBUG_PRINTF (DEBUG_PRINT_TOCRT, "crt %p", (void *)r);

  // intvec_alloc (tmp1, deg, 1); // XXX costs more than NTTs ...

  if (r->crtrep == NULL)
    r->crtrep = _alloc (_sizeof_crtrep_data (ring));
  r->crt = 1;

  hexl_ntts = (deg == 64 ? hexl_ntt_d64 : hexl_ntt_d128);

  crtrep = r->crtrep;
  _POLYRING_FOREACH_P (ring, i)
  {
#if 0
    int_set_i64 (p, ring->moduli[i]->p);

    intvec_mod (tmp1, coeffs, p); /* reduce to (-p,p) */
    // not needed, ntt input can be in (-p,p)
    // intvec_redc (tmp1, tmp1, p); /* reduce to (-(p-1)/2,(p-1)/2) */

    crtcoeff = _get_crtcoeff (crtrep, i, 0, deg);
#if CRTCOEFF_NBITS == 32
    intvec_get_i32 (crtcoeff, tmp1);
#else
#error "not implemented"
#endif
#else
    // XXX
    unsigned int j;
    for (j = 0; j < deg; j++)
      {
        crtcoeff = _get_crtcoeff (crtrep, i, j, deg);
        int_ptr coeff = intvec_get_elem (coeffs, j);
        *crtcoeff = int_mod_XXX_hexl (coeff, ring->moduli[i]->p);
      }
    crtcoeff = _get_crtcoeff (crtrep, i, 0, deg);
    // XXX
    //_ntt (crtcoeff, ring->moduli[i], deg);

    void *hexl_ntt = hexl_ntts[i];

    crtcoeff = _get_crtcoeff (crtrep, i, 0, deg);
    hexl_ntt_fwd (hexl_ntt, crtcoeff, 1, crtcoeff, 1);
#endif
  }
  // XXXintvec_free (tmp1);
}

void
poly_fromcrt (poly_t r)
{
  if (r->crt == 0)
    return;

  polyring_srcptr ring = poly_get_ring (r);
  const unsigned int deg = polyring_get_deg (ring);
  int_srcptr mod = polyring_get_mod (ring);
  intvec_ptr coeffs = _get_coeffvec (r);
  const unsigned int s = ring->nmoduli;
  crtcoeff_t *crtrep = r->crtrep;
  modulus_srcptr crtmod = ring->moduli[s - 1];
  crtcoeff_t *crtcoeff, x;
  crtcoeff_t p;
  int_ptr coeff;
  void **hexl_ntts;
  unsigned int i, j;
  float z;
  long round_z;
  INT_T (Pz, crtmod->P->nlimbs * 2);
  INT_T (tmp2, crtmod->P->nlimbs * 2);
  INT_T (tmp1, crtmod->P->nlimbs);
  INT_T (tmp3, mod->nlimbs);

  DEBUG_PRINTF (DEBUG_PRINT_FROMCRT, "icrt %p", (void *)r);

  crtmod = ring->moduli[s - 1];
  hexl_ntts = (deg == 64 ? hexl_ntt_d64 : hexl_ntt_d128);

  _POLYRING_FOREACH_P (ring, i)
  {
    crtcoeff = _get_crtcoeff (crtrep, i, 0, deg);
    // XXX_intt (crtcoeff, ring->moduli[i], deg);

    void *hexl_ntt = hexl_ntts[i];

    crtcoeff = _get_crtcoeff (crtrep, i, 0, deg);
    hexl_ntt_inv (hexl_ntt, crtcoeff, 1, crtcoeff, 1);
    for (j = 0; j < deg; j++)
      {
        crtcoeff = _get_crtcoeff (crtrep, i, j, deg);
        if (*crtcoeff > (ring->moduli[i]->p - 1) / 2)
          *crtcoeff = *crtcoeff + ring->moduli[i]->p;
      }
  }

  _VEC_FOREACH_ELEM (coeffs, j)
  {
    coeff = intvec_get_elem (coeffs, j);
    int_set_i64 (Pz, 0);
    z = 0;

    _POLYRING_FOREACH_P (ring, i)
    {
      p = ring->moduli[i]->p;

      crtcoeff = _get_crtcoeff (crtrep, i, j, deg);
      // if (*crtcoeff > (p - 1) / 2)
      //   *crtcoeff -= p;
      // else if (*crtcoeff < -(p - 1) / 2)
      //   *crtcoeff += p;
      //  XXX

      x = ((crtcoeff_dbl_t)crtmod->k[i] * (*crtcoeff)) % p; // XXX in (-p,p)
      // XXX
      // if (x > (p - 1) / 2)
      //  x -= p;
      // else if (x < -(p - 1) / 2)
      //  x += p;
      // XXX

      int_set_i64 (tmp1, x);
      int_mul (tmp2, crtmod->Pp[i], tmp1);
      int_add (Pz, Pz, tmp2);

      z += (float)x / p;
    }

    round_z = lroundf (z);

    int_set_i64 (tmp1, round_z);
    int_mul (tmp2, crtmod->P, tmp1);

    int_sub (Pz, Pz, tmp2);

    int_mod (tmp3, Pz, mod);
    int_redc (tmp3, tmp3, mod);

    int_set (coeff, tmp3);
  }

  _free (r->crtrep, _sizeof_crtrep_data (ring));
  r->crtrep = NULL;
  r->crt = 0;
}

void
poly_grandom (poly_t r, unsigned int log2o, const uint8_t seed[32],
              uint32_t dom)
{
  _poly_grandom (r, log2o, seed, dom);
}

void
poly_urandom (poly_t r, const int_t mod, unsigned int log2mod,
              const uint8_t seed[32], uint32_t dom)
{
  _poly_urandom (r, mod, log2mod, seed, dom);
}

void
poly_urandom_autostable (poly_t r, int64_t bnd, unsigned int log2,
                         const uint8_t seed[32], uint32_t dom)
{
  _poly_urandom_autostable (r, bnd, log2, seed, dom);
}

void
poly_urandom_bnd (poly_t r, const int_t lo, const int_t hi,
                  const uint8_t seed[32], uint32_t dom)
{
  _poly_urandom_bnd (r, lo, hi, seed, dom);
}

void
poly_brandom (poly_t r, unsigned int k, const uint8_t seed[32], uint32_t dom)
{
  _poly_brandom (r, k, seed, dom);
}

/* output coeffs in (-q,q) */
void
poly_tracemap (poly_t r, poly_t a)
{
  INT_T (prod, r->ring->q->nlimbs * 2);
  intvec_ptr rcoeffs = r->coeffs;
  intvec_ptr acoeffs = a->coeffs;
  const unsigned int nelems = rcoeffs->nelems / 2;
  unsigned int i;
  int_ptr ri, ai;
  int_srcptr ris, ais;

  ASSERT_ERR (r->ring == a->ring);

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt tracemap");
#endif
  _fromcrt (a);

  intvec_set_elem (rcoeffs, 0, intvec_get_elem_src (acoeffs, 0));

  for (i = 1; i < nelems; i++)
    {
      ri = intvec_get_elem (rcoeffs, i);
      ai = intvec_get_elem (acoeffs, i);
      ais = intvec_get_elem_src (acoeffs, acoeffs->nelems - i);
      int_sub (ri, ai, ais);
      int_mul (prod, ri, r->ring->inv2);
      int_mod (ri, prod, r->ring->q);
      int_redc (ri, ri, r->ring->q);
    }
  intvec_set_elem_i64 (rcoeffs, nelems, 0);
  for (i = nelems + 1; i < rcoeffs->nelems; i++)
    {
      ri = intvec_get_elem (rcoeffs, i);
      ris = intvec_get_elem_src (rcoeffs, 2 * nelems - i);
      int_neg (ri, ris);
    }

  r->crt = 0;
}

void
poly_auto (poly_t r, poly_t a)
{
  ASSERT_ERR (r != a);

#ifdef XXX
  /* o-1 auto in crt form reverses crt coeffs */
  if (a->crt == 1)
    {
      crtcoeff_t *ca, *cr;
      polyring_srcptr ring = poly_get_ring (a);
      const unsigned int deg = polyring_get_deg (ring);
      unsigned int i, j;

      if (r->crtrep == NULL)
        r->crtrep = _alloc (_sizeof_crtrep_data (ring));
      r->crt = 1;

      _POLYRING_FOREACH_P (ring, i)
      {
        for (j = 0; j < deg; j++)
          {
            cr = _get_crtcoeff (r->crtrep, i, j, deg);
            ca = _get_crtcoeff (a->crtrep, i, deg - j - 1, deg);

            *cr = *ca;
          }
      }
    }
  else
    {
      intvec_ptr rcoeffvec = _get_coeffvec (r);
      intvec_ptr acoeffvec = _get_coeffvec (a);

      intvec_auto (rcoeffvec, acoeffvec);
    }

  /* o-1 auto in crt form reverses crt coeffs */
  if (a->crt == 1)
    {
      crtcoeff_t *ca, *cb, tmp;
      polyring_srcptr ring = poly_get_ring (a);
      const unsigned int deg = polyring_get_deg (ring);
      unsigned int i, j;

      _POLYRING_FOREACH_P (ring, i)
      {
        for (j = 0; j < deg / 2; j++)
          {
            ca = _get_crtcoeff (a->crtrep, i, j, deg);
            cb = _get_crtcoeff (a->crtrep, i, deg - j - 1, deg);

            tmp = *cb;
            *cb = *ca;
            *ca = tmp;
          }
      }
    }
  else
    {
      intvec_ptr rcoeffvec = _get_coeffvec (r);
      intvec_ptr acoeffvec = _get_coeffvec (a);

      intvec_auto (rcoeffvec, acoeffvec);
    }
#endif

  intvec_ptr rcoeffvec = _get_coeffvec (r);
  intvec_ptr acoeffvec = _get_coeffvec (a);

  ASSERT_ERR (r->ring == a->ring);

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt auto");
#endif
  _fromcrt (a);

  intvec_auto (rcoeffvec, acoeffvec);
  r->crt = 0;
}

void
poly_auto_self (poly_t r)
{
  /* o-1 auto in crt form reverses crt coeffs */
  if (r->crt == 1)
    {
      crtcoeff_t *crb, *cre, tmp;
      polyring_srcptr ring = poly_get_ring (r);
      const unsigned int deg = polyring_get_deg (ring);
      unsigned int i, j;

      _POLYRING_FOREACH_P (ring, i)
      {
        for (j = 0; j < deg / 2; j++)
          {
            crb = _get_crtcoeff (r->crtrep, i, j, deg);
            cre = _get_crtcoeff (r->crtrep, i, deg - j - 1, deg);

            tmp = *crb;
            *crb = *cre;
            *cre = tmp;
          }
      }
    }
  else
    {
      intvec_ptr rcoeffvec = _get_coeffvec (r);

      intvec_auto_self (rcoeffvec);
    }
}

void
poly_dcompress_power2round (poly_t r, poly_t a,
                            const dcompress_params_t params)
{
  intvec_ptr rcoeffvec = _get_coeffvec (r);
  intvec_ptr acoeffvec = _get_coeffvec (a);

  ASSERT_ERR (r->ring == a->ring);

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (a->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s",
                  "implicit icrt dcompress power2round");
#endif
  _fromcrt (a);

  dcompress_power2round (rcoeffvec, acoeffvec, params);
  r->crt = 0;
}

void
poly_dcompress_decompose (poly_t r1, poly_t r0, poly_t r,
                          const dcompress_params_t params)
{
  intvec_ptr r0coeffvec = _get_coeffvec (r0);
  intvec_ptr r1coeffvec = _get_coeffvec (r1);
  intvec_ptr rcoeffvec = _get_coeffvec (r);

  ASSERT_ERR (r->ring == r0->ring);
  ASSERT_ERR (r->ring == r1->ring);

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (r->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt dcompress decompose");
#endif
  _fromcrt (r);

  dcompress_decompose (r1coeffvec, r0coeffvec, rcoeffvec, params);
  r1->crt = 0;
  r0->crt = 0;
}

void
poly_dcompress_use_ghint (poly_t ret, poly_t y, poly_t r,
                          const dcompress_params_t params)
{
  intvec_ptr retcoeffvec = _get_coeffvec (ret);
  intvec_ptr ycoeffvec = _get_coeffvec (y);
  intvec_ptr rcoeffvec = _get_coeffvec (r);

  ASSERT_ERR (ret->ring == y->ring);
  ASSERT_ERR (ret->ring == r->ring);

#if DEBUGINFO == DEBUGINFO_ENABLED
  if (y->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt dcompress use ghint");
#endif
  _fromcrt (y);
#if DEBUGINFO == DEBUGINFO_ENABLED
  if (r->crt == 1)
    DEBUG_PRINTF (DEBUG_LEVEL >= 2, "%s", "implicit icrt dcompress use ghint");
#endif
  _fromcrt (r);

  dcompress_use_ghint (retcoeffvec, ycoeffvec, rcoeffvec, params);
  ret->crt = 0;
}

void
poly_dcompress_make_ghint (poly_t ret, poly_t z, poly_t r,
                           const dcompress_params_t params)
{
  intvec_ptr retcoeffvec = _get_coeffvec (ret);
  intvec_ptr zcoeffvec = _get_coeffvec (z);
  intvec_ptr rcoeffvec = _get_coeffvec (r);

  ASSERT_ERR (ret->ring == z->ring);
  ASSERT_ERR (ret->ring == r->ring);

  _fromcrt (z);
  _fromcrt (r);

  dcompress_make_ghint (retcoeffvec, zcoeffvec, rcoeffvec, params);
  ret->crt = 0;
}

void
poly_l2sqr (int_t r, poly_t a)
{
  _fromcrt (a);

  intvec_l2sqr (r, poly_get_coeffvec (a));
}

void
poly_linf (int_t r, poly_t a)
{
  _fromcrt (a);

  intvec_linf (r, poly_get_coeffvec (a));
}

size_t
poly_out_str (FILE *stream, int base, poly_t a)
{
  _fromcrt (a);
  return intvec_out_str (stream, base, a->coeffs);
}

void
poly_dump (poly_t a)
{
  poly_out_str (stdout, 10, a);
  fprintf (stdout, "\n");
  fflush (stdout);
}

static void
_poly_urandom (poly_t r, const int_t mod, unsigned int log2mod,
               const uint8_t seed[32], uint64_t dom)
{
  intvec_ptr coeffvec = _get_coeffvec (r);

  r->crt = 0;
  _intvec_urandom (coeffvec, mod, log2mod, seed, dom);
}

static void
_poly_brandom (poly_t r, unsigned int k, const uint8_t seed[32], uint64_t dom)
{
  intvec_ptr coeffvec = _get_coeffvec (r);

  r->crt = 0;
  _intvec_brandom (coeffvec, k, seed, dom);
}

static void
_poly_grandom (poly_t r, unsigned int log2o, const uint8_t seed[32],
               uint64_t dom)
{
  intvec_ptr coeffvec = _get_coeffvec (r);

  r->crt = 0;
  _intvec_grandom (coeffvec, log2o, seed, dom);
}

static void
_poly_urandom_autostable (poly_t r, int64_t bnd, unsigned int log2,
                          const uint8_t seed[32], uint64_t dom)
{
  _intvec_urandom_autostable (r->coeffs, bnd, log2, seed, dom);
  r->crt = 0;
}

static void
_poly_urandom_bnd (poly_t r, const int_t lo, const int_t hi,
                   const uint8_t seed[32], uint64_t dom)
{
  intvec_ptr coeffvec = _get_coeffvec (r);

  r->crt = 0;
  _urandom_bnd (coeffvec, lo, hi, seed, dom);
}
