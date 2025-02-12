// #define LAZER_BENCH_SAMPLE 1
#include "lazer-bench.h"

// XXX#include "bench-params.h"
#include "../python/demo_params.h" //XXX
#include "lazer.h"
#include <errno.h>
#include <getopt.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>

#define HELP                                                                  \
  "Usage: bench [-hvs] [-o <file>] [-p <pre>]\n"                              \
  "Options:\n"                                                                \
  "  -h, --help           Display this information.\n"                        \
  "  -v, --verbose        Verbose output.\n"                                  \
  "  -s, --stack          Use a stack allocator.\n"                           \
  "  -o, --output <file>  Write benchmark results to <file>.\n"               \
  "  -p, --prefix <pre>   Run only benchmarks with prefix <pre>.\n"

#define STACK_SIZE 1000000000
#define SP_ALIGN 8
#define ALIGN(sp)                                                             \
  (uint8_t *)((((size_t)(sp) + SP_ALIGN - 1) / SP_ALIGN) * SP_ALIGN)

/* automake test exit status */
#define RC_BENCH_PASS 0
#define RC_BENCH_SKIP 77
#define RC_BENCH_FAIL 1
#define RC_BENCH_ERROR 99

#define BENCH_SUCC() exit (RC_BENCH_PASS)
#define BENCH_FAIL() exit (RC_BENCH_FAIL)
#define BENCH_ERR(fmt, ...)                                                   \
  do                                                                          \
    {                                                                         \
      fprintf (stderr, "%s: " fmt "\n", prg, __VA_ARGS__);                    \
      exit (RC_BENCH_ERROR);                                                  \
    }                                                                         \
  while (0)

struct benchmark
{
  const char *name;
  void (*func) (const char *, unsigned long);
  unsigned long miniter;
  const char *description;
};

static void clear (void);

static void *stack_alloc (size_t);
static void *stack_realloc (void *, size_t, size_t);
static void stack_free (void *, size_t);

static void bench_shake128_absorb512 (const char *, unsigned long);
static void bench_shake128_squeeze8 (const char *, unsigned long);
static void bench_shake128_squeeze16 (const char *, unsigned long);
static void bench_shake128_squeeze32 (const char *, unsigned long);
static void bench_shake128_squeeze64 (const char *, unsigned long);
static void bench_shake128_squeeze128 (const char *, unsigned long);
static void bench_shake128_squeeze256 (const char *, unsigned long);
static void bench_shake128_squeeze512 (const char *, unsigned long);
static void bench_rng8 (const char *, unsigned long);
static void bench_rng16 (const char *, unsigned long);
static void bench_rng32 (const char *, unsigned long);
static void bench_rng64 (const char *, unsigned long);
static void bench_rng128 (const char *, unsigned long);
static void bench_rng256 (const char *, unsigned long);
static void bench_rng512 (const char *, unsigned long);
static void bench_poly_urandom (const char *, unsigned long);
static void bench_poly_brandom (const char *, unsigned long);
static void bench_poly_grandom (const char *, unsigned long);
static void bench_poly_add (const char *, unsigned long);
static void bench_poly_mul (const char *, unsigned long);

static int verbose_p;
static const char *out_file, *pre_str;
static FILE *out_fp;
static const char *prg;

static const struct option opts[] = { { "help", no_argument, NULL, 'h' },
                                      { "verbose", no_argument, NULL, 'v' },
                                      { "stack", no_argument, NULL, 's' },
                                      { "out", required_argument, NULL, 'o' },
                                      { "pre", required_argument, NULL, 'n' },
                                      { NULL, 0, NULL, 0 } };

static struct benchmark benchmarks[] = {
  { "shake128-a512", bench_shake128_absorb512, 100000, "absorb 512 bytes" },
  { "shake128-s8", bench_shake128_squeeze8, 100000, "squeeze 8 bytes" },
  { "shake128-s16", bench_shake128_squeeze16, 100000, "squeeze 16 bytes" },
  { "shake128-s32", bench_shake128_squeeze32, 100000, "squeeze 32 bytes" },
  { "shake128-s64", bench_shake128_squeeze64, 100000, "squeeze 64 bytes" },
  { "shake128-s128", bench_shake128_squeeze128, 100000, "squeeze 128 bytes" },
  { "shake128-s256", bench_shake128_squeeze256, 100000, "squeeze 256 bytes" },
  { "shake128-s512", bench_shake128_squeeze512, 100000, "squeeze 512 bytes" },
  { "rng_urandom8", bench_rng8, 100000, "8 random bytes" },
  { "rng_urandom16", bench_rng16, 100000, "16 random bytes" },
  { "rng_urandom32", bench_rng32, 100000, "32 random bytes" },
  { "rng_urandom64", bench_rng64, 100000, "64 random bytes" },
  { "rng_urandom128", bench_rng128, 100000, "128 random bytes" },
  { "rng_urandom256", bench_rng256, 100000, "256 random bytes" },
  { "rng_urandom512", bench_rng512, 100000, "512 random bytes" },
  { "poly_urandom", bench_poly_urandom, 10000, "urandom polynomial" },
  { "poly_brandom", bench_poly_brandom, 10000, "brandom polynomial" },
  { "poly_grandom", bench_poly_grandom, 10000, "grandom polynomial" },
  { "poly_add", bench_poly_add, 10000, "polynomial addition" },
  { "poly_mul", bench_poly_mul, 10000, "polynomial multiplication" },
};
static size_t nbenchmarks = sizeof (benchmarks) / sizeof (benchmarks[0]);

int
main (int argc, char *argv[])
{
  prg = argv[0];
  int c, skip_p;
  size_t i;

  atexit (clear);

  while (1)
    {
      c = getopt_long (argc, argv, "hvso:p:", opts, NULL);
      if (c == -1)
        break;

      switch (c)
        {
        case 'h':
          printf (HELP);
          BENCH_SUCC ();
          break;
        case 'v':
          verbose_p = 1;
          break;
        case 's':
          mp_set_memory_functions (stack_alloc, stack_realloc, stack_free);
          lazer_set_memory_functions (stack_alloc, stack_realloc, stack_free);
          break;
        case 'o':
          out_file = optarg;
          break;
        case 'p':
          pre_str = optarg;
          break;
        default:
          BENCH_FAIL ();
        }
    }

  if (out_file != NULL)
    {
      out_fp = fopen (out_file, "w");
      if (out_fp == NULL)
        BENCH_ERR ("could not open %s: %s", out_file, strerror (errno));
    }

  for (i = 0; i < nbenchmarks; i++)
    {
      if (pre_str != NULL
          && strncmp (benchmarks[i].name, pre_str, strlen (pre_str)) != 0)
        skip_p = 1;
      else
        skip_p = 0;

      if (verbose_p)
        {
          printf ("%s benchmark (%lu/%lu): %s: %s\n",
                  skip_p ? "Skipping" : "Running", i + 1, nbenchmarks,
                  benchmarks[i].name, benchmarks[i].description);
        }

      if (skip_p)
        continue;

      benchmarks[i].func (benchmarks[i].name, benchmarks[i].miniter);
    }
  if (verbose_p)
    printf ("\n");

  if (out_fp != NULL)
    bench_fprint (out_fp);
  else
    bench_fprint (stdout);
  BENCH_SUCC ();
}

static void
clear (void)
{
  if (out_fp != NULL)
    fclose (out_fp);

  // XXX mpfr_free_cache ();
}

static void *
stack_alloc (size_t len)
{
  static uint8_t stack[STACK_SIZE];
  static uint8_t *sp = stack;
  uint8_t *sp_aligned = ALIGN (sp);

  sp = sp_aligned + len;
  if (sp > stack + STACK_SIZE)
    return NULL;
  return sp_aligned;
}

static void *
stack_realloc (void *omem, size_t olen, size_t nlen)
{
  void *nmem;

  nmem = stack_alloc (nlen);
  memcpy (nmem, omem, nlen < olen ? nlen : olen);
  return nmem;
}

static void
stack_free (void *mem, size_t len)
{
  (void)mem; /* unused */
  (void)len; /* unused */
}

static void
bench_shake128_absorb512 (const char *name, unsigned long miniter)
{
  shake128_state_t state;
  uint8_t buf[512];

  shake128_init (state);

  BENCH_BEGIN (name, miniter, 1);

  BENCH_CLOCK_START ();
  shake128_absorb (state, buf, sizeof (buf));
  BENCH_CLOCK_STOP ();
  BENCH_END ();

  shake128_clear (state);
}

#define BENCH_SHAKE128_SQUEEZE(x)                                             \
  static void bench_shake128_squeeze##x (const char *name,                    \
                                         unsigned long miniter)               \
  {                                                                           \
    shake128_state_t state;                                                   \
    uint8_t buf[x];                                                           \
                                                                              \
    shake128_init (state);                                                    \
    shake128_squeeze (state, buf, sizeof (buf)); /* finalize */               \
                                                                              \
    BENCH_BEGIN (name, miniter, 1);                                           \
                                                                              \
    BENCH_CLOCK_START ();                                                     \
    shake128_squeeze (state, buf, sizeof (buf));                              \
    BENCH_CLOCK_STOP ();                                                      \
    BENCH_END ();                                                             \
                                                                              \
    shake128_clear (state);                                                   \
  }

BENCH_SHAKE128_SQUEEZE (8)
BENCH_SHAKE128_SQUEEZE (16)
BENCH_SHAKE128_SQUEEZE (32)
BENCH_SHAKE128_SQUEEZE (64)
BENCH_SHAKE128_SQUEEZE (128)
BENCH_SHAKE128_SQUEEZE (256)
BENCH_SHAKE128_SQUEEZE (512)
#undef BENCH_SHAKE128_SQUEEZE

#define BENCH_RNG(x)                                                          \
  static void bench_rng##x (const char *name, unsigned long miniter)          \
  {                                                                           \
    rng_state_t state;                                                        \
    uint8_t seed[32] = { 0 }, buf[x] = { 0 };                                 \
    uint64_t dom = 0;                                                         \
                                                                              \
    rng_init (state, seed, dom);                                              \
                                                                              \
    BENCH_BEGIN (name, miniter, 1);                                           \
                                                                              \
    BENCH_CLOCK_START ();                                                     \
    rng_urandom (state, buf, sizeof (buf));                                   \
    BENCH_CLOCK_STOP ();                                                      \
    BENCH_END ();                                                             \
                                                                              \
    rng_clear (state);                                                        \
  }

BENCH_RNG (8)
BENCH_RNG (16)
BENCH_RNG (32)
BENCH_RNG (64)
BENCH_RNG (128)
BENCH_RNG (256)
BENCH_RNG (512)

#undef BENCH_RNG

static void
bench_poly_urandom (const char *name, unsigned long miniter)
{
  polyring_srcptr Rq = _p1_params->quad_eval->quad_many->ring;
  int_srcptr mod = polyring_get_mod (Rq);
  uint32_t dom = 0;
  uint8_t seed[32] = { 0 };
  POLY_T (r, Rq);

  bytes_urandom (seed, sizeof (seed));

  BENCH_BEGIN (name, miniter, 1);
  BENCH_CLOCK_START ();
  dom++;
  poly_urandom (r, mod, Rq->log2q, seed, dom);
  BENCH_CLOCK_STOP ();
  BENCH_END ();
}

static void
bench_poly_brandom (const char *name, unsigned long miniter)
{
  polyring_srcptr Rq = _p1_params->quad_eval->quad_many->ring;
  uint32_t dom = 0;
  uint8_t seed[32] = { 0 };
  POLY_T (r, Rq);

  bytes_urandom (seed, sizeof (seed));

  BENCH_BEGIN (name, miniter, 1);
  BENCH_CLOCK_START ();
  dom++;
  poly_brandom (r, 1, seed, dom);
  BENCH_CLOCK_STOP ();
  BENCH_END ();
}

static void
bench_poly_grandom (const char *name, unsigned long miniter)
{
  polyring_srcptr Rq = _p1_params->quad_eval->quad_many->ring;
  int_srcptr mod = polyring_get_mod (Rq);
  uint32_t dom = 0;
  uint8_t seed[32] = { 0 };
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  INT_T (one, mod->nlimbs);
  POLY_T (r, Rq);

  int_set_i64 (one, 1);

  int_set (hi, Rq->q);
  int_sub (hi, hi, one);
  int_rshift (hi, hi, 1);
  int_neg (lo, hi);

  bytes_urandom (seed, sizeof (seed));

  BENCH_BEGIN (name, miniter, 1);
  BENCH_CLOCK_START ();
  dom++;
  poly_grandom (r, 20, seed, dom);
  BENCH_CLOCK_STOP ();
  BENCH_END ();
}

static void
bench_poly_add (const char *name, unsigned long miniter)
{
  polyring_srcptr Rq = _p1_params->quad_eval->quad_many->ring;
  int_srcptr mod = polyring_get_mod (Rq);
  uint32_t dom;
  uint8_t seed[32] = { 0 };
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  POLY_T (r, Rq);

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 1));

  bytes_urandom (seed, sizeof (seed));

  dom = 1;
  poly_urandom_bnd (a, lo, hi, seed, dom);

  dom = 2;
  poly_urandom_bnd (b, lo, hi, seed, dom);

  BENCH_BEGIN (name, miniter, 1);
  BENCH_CLOCK_START ();
  poly_add (r, a, b, 0);
  BENCH_CLOCK_STOP ();
  BENCH_END ();
}

static void
bench_poly_mul (const char *name, unsigned long miniter)
{
  polyring_srcptr Rq = _p1_params->quad_eval->quad_many->ring;
  int_srcptr mod = polyring_get_mod (Rq);
  uint32_t dom;
  uint8_t seed[32] = { 0 };
  INT_T (lo, mod->nlimbs);
  INT_T (hi, mod->nlimbs);
  POLY_T (a, Rq);
  POLY_T (b, Rq);
  POLY_T (_a, Rq);
  POLY_T (_b, Rq);
  POLY_T (_r, Rq);

  int_set_i64 (lo, -((int64_t)1 << (Rq->log2q - 1)));
  int_set_i64 (hi, (int64_t)1 << (Rq->log2q - 1));

  bytes_urandom (seed, sizeof (seed));

  dom = 1;
  poly_urandom_bnd (a, lo, hi, seed, dom);

  dom = 2;
  poly_urandom_bnd (b, lo, hi, seed, dom);

  BENCH_BEGIN (name, miniter, 1);

  //poly_set (_a, a);
  //poly_set (_b, b);

  BENCH_CLOCK_START ();
  //poly_tocrt (_a);
  //poly_tocrt (_b);
  poly_mul (_r, _a, _b);
  //poly_fromcrt (_r);
  BENCH_CLOCK_STOP ();

  BENCH_END ();
}
