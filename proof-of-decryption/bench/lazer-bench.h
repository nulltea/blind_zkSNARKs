#ifndef LAZER_BENCH_H
#define LAZER_BENCH_H
#include <errno.h>
#include "lazer.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Define for sampling */
// #define LAZER_BENCH_SAMPLE

#define CONSTRUCTOR __attribute__ ((constructor))
#define UNUSED __attribute__ ((unused))

#define BENCH_MAX_SAMPLES 100000 /* must be even */
#define BENCH_NAN (0.0 / 0.0)
#define BENCH_FILE                                                            \
  (strrchr (__FILE__, '/') ? strrchr (__FILE__, '/') + 1 : __FILE__)

CONSTRUCTOR static void bench_init (void);
UNUSED static int bench_compare (const void *op1, const void *op2);
static long double bench_get_cpu_freq (void);
UNUSED static unsigned long long bench_diff (const struct timespec *,
                                             const struct timespec *);
UNUSED static long double bench_median (const unsigned long long *, size_t);
UNUSED static long double bench_average (const unsigned long long *, size_t);
UNUSED static void bench_fprint (FILE *stream);

UNUSED static unsigned long long bench_samples[BENCH_MAX_SAMPLES];
UNUSED struct timespec bench_start, bench_stop;
UNUSED static long double bench_q1, bench_q2, bench_q3, bench_avg;
UNUSED static const char *bench_name;
UNUSED unsigned long bench_max;
UNUSED unsigned long bench_i;
static char bench_buf[1000000];
static FILE *bench_fp;
static long double bench_mhz;

#define BENCH_BEGIN(name, iter, run_p)                                        \
  bench_name = (name);                                                        \
  if (!run_p)                                                                 \
    {                                                                         \
      bench_max = 1;                                                          \
    }                                                                         \
  else                                                                        \
    {                                                                         \
      bench_max = ((iter) + 1) / 2 * 2;                                       \
      if (bench_max < 2)                                                      \
        bench_max = 2;                                                        \
      if (bench_max > BENCH_MAX_SAMPLES)                                      \
        bench_max = BENCH_MAX_SAMPLES;                                        \
    }                                                                         \
  for (bench_i = 0; bench_i < bench_max; bench_i++)                           \
    {                                                                         \
      (void)0

#ifndef LAZER_BENCH_SAMPLE
#define BENCH_CLOCK_START() clock_gettime (CLOCK_MONOTONIC_RAW, &bench_start)
#else
#define BENCH_CLOCK_START() (void)0
#endif

#ifndef LAZER_BENCH_SAMPLE
#define BENCH_CLOCK_STOP()                                                    \
  clock_gettime (CLOCK_MONOTONIC_RAW, &bench_stop);                           \
  bench_samples[bench_i] = bench_diff (&bench_start, &bench_stop)
#else
#define BENCH_CLOCK_STOP() (void)0
#endif

#ifndef LAZER_BENCH_SAMPLE
#define BENCH_END()                                                           \
  }                                                                           \
  if (bench_max > 1)                                                          \
    {                                                                         \
      qsort (bench_samples, bench_max, sizeof (bench_samples[0]),             \
             bench_compare);                                                  \
      bench_q1 = bench_median (&bench_samples[0], bench_max / 2);             \
      bench_q2 = bench_median (&bench_samples[0], bench_max);                 \
      bench_q3 = bench_median (&bench_samples[bench_max / 2], bench_max / 2); \
      bench_avg = bench_average (&bench_samples[0], bench_max);               \
      bench_q1 /= (long double)1000;                                          \
      bench_q2 /= (long double)1000;                                          \
      bench_q3 /= (long double)1000;                                          \
      bench_avg /= (long double)1000;                                         \
      fprintf (                                                               \
          bench_fp,                                                           \
          "%-14s, %14.2Lf, %14.2Lf, %14.2Lf, %14.2Lf, %14.2Lf, %14.2Lf, "     \
          "%14.2Lf, %14.2Lf\n",                                               \
          bench_name, bench_avg *bench_mhz, bench_q1 *bench_mhz,              \
          bench_q2 *bench_mhz, bench_q3 *bench_mhz, bench_avg, bench_q1,      \
          bench_q2, bench_q3);                                                \
    }                                                                         \
  (void)0
#else
#define BENCH_END()                                                           \
  }                                                                           \
  (void)0
#endif

static void
bench_fprint (FILE *stream)
{
  if (bench_fp != NULL)
    fclose (bench_fp);

  /* no benchmarks have been run. */
  if (strcmp (bench_buf, "") == 0)
    return;

  fprintf (stream,
           "%-14s, %-14s, %-14s, %-14s, %-14s, %-14s, %-14s, %-14s, %-14s\n",
           "benchmark", "avg [cycles]", "Q1 [cycles]", "Q2 [cycles]",
           "Q3 [cycles]", "avg [us]", "Q1 [us]", "Q2 [us]", "Q3 [us]");
  fprintf (stream, "%s", bench_buf);
}

static void
bench_init (void)
{
  bench_mhz = bench_get_cpu_freq ();
  bench_fp = fmemopen (bench_buf, sizeof (bench_buf), "w");
  if (bench_fp == NULL)
    {
      fprintf (stderr, "%s:%d: fmemopen failed: %s\n", BENCH_FILE, __LINE__,
               strerror (errno));
      abort ();
    }
}

static unsigned long long
bench_diff (const struct timespec *start, const struct timespec *stop)
{
  return ((unsigned long long)stop->tv_sec - start->tv_sec) * 1000000000ULL
         + ((unsigned long long)stop->tv_nsec - start->tv_nsec);
}

static long double
bench_median (const unsigned long long *list, size_t len)
{
  long double median;

  if (len % 2 == 0)
    median = ((long double)list[len / 2 - 1] + list[len / 2]) / 2;
  else
    median = (long double)list[len / 2];
  return median;
}

static long double
bench_average (const unsigned long long *list, size_t len)
{
  unsigned long long sum = 0;
  size_t i;

  for (i = 0; i < len; i++)
    sum += (unsigned long long)list[i];
  return (long double)sum / len;
}

static int
bench_compare (const void *op1, const void *op2)
{
  const unsigned long long a = *((unsigned long long *)op1);
  const unsigned long long b = *((unsigned long long *)op2);

  if (a < b)
    return -1;
  else if (a > b)
    return 1;
  else
    return 0;
}

static long double
bench_get_cpu_freq (void)
{
  double mhz;
  FILE *fh;
  int rv;

  mhz = -1;

  fh = fopen ("/proc/cpuinfo", "r");
  if (fh == NULL)
    goto ret;

  while (1)
    {
      rv = fscanf (fh, "cpu MHz : %lf", &mhz);
      if (rv > 0)
        goto ret;
      if (rv == 0)
        rv = fscanf (fh, "%*[^\n]\n");
      if (rv < 0)
        goto ret;
    }
ret:
  if (fh != NULL)
    fclose (fh);
  if (mhz <= 0)
    {
      fprintf (stderr, "bench: failed to retrieve cpu MHz\n");
      mhz = BENCH_NAN;
    }
  return (long double)mhz;
}

#endif
