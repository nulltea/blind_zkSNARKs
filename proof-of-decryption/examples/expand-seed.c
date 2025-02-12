/**
 * Securely generate, expand and erase a seed.
 */
#include "lazer.h"
#include <stdio.h>

#define SEEDLEN 32
#define OUTLEN 32

int
main (void)
{
  unsigned char seed[SEEDLEN], out[OUTLEN];
  shake128_state_t state;
  int rc = 1, i;
  size_t len;

  /* Generate a 32 byte uniform seed. */
  bytes_urandom (seed, SEEDLEN);

  printf ("seed:\n");

  /* Print seed in hex. */
  len = bytes_out_str (stdout, seed, SEEDLEN);
  if (len != SEEDLEN * 2) /* expect two hex digits per byte */
    goto ret;

  /* Use shake128 to expand the seed. */
  shake128_init (state);

  /* Absorb seed in two steps (just to show that one can do
   * absorb multiple times). */
  shake128_absorb (state, seed, SEEDLEN / 2);
  shake128_absorb (state, seed + SEEDLEN / 2, SEEDLEN / 2);

  printf ("\n\nuniform bytes:\n");

  /* Squeeze four times 32 uniformly (pseudo) random bytes. */
  for (i = 0; i < 4; i++)
    {
      shake128_squeeze (state, out, OUTLEN);
      len = bytes_out_str (stdout, out, OUTLEN);
      if (len != SEEDLEN * 2) /* expect two hex digits per byte */
        goto ret;
      printf ("\n");
    }

  rc = 0;
ret:
  /* Erase seed and hash function state securely.
   * A call to memset would probably be optimized out here. */
  shake128_clear (state);
  bytes_clear (seed, SEEDLEN);
  if (rc != 0)
    printf ("\n\nwriting to stdout failed.\n");
  return rc;
}
