/**
 * Working with integers.
 */
#include "lazer.h"
#include <stdio.h>

int
main (void)
{
  int_t a, b, c, d;

  int_init (a, 128);
  int_init (b, 128);
  int_init (c, 128);
  int_init (d, 256);

  int_set_i64 (a, -1);
  int_set_i64 (b, 0);
  int_set_i64 (c, 1);

  int_set_i64 (a, 7);
  int_set_i64 (b, 3);
  int_set_i64 (c, 11);
  int_mul (d, a, b);
  int_mod (a, d, c);

  int_set_i64 (a, 4);
  int_set_i64 (b, -5);
  int_set_i64 (c, 11);
  int_mul (d, a, b);
  int_mod (a, d, c);

  int_set_i64 (a, -4);
  int_mul (d, a, a);

  int_set_i64 (a, -4);
  int_set_i64 (b, 5);
  int_mul (d, a, b);

  int_set_i64 (a, -3);
  int_set_i64 (b, 2);
  int_mul (d, a, a);

  int_set_i64 (a, 0);
  int_set_i64 (b, -1);
  int_mul (d, a, b); //XXX -0

  int_set_i64 (b, 7);
  int_set_i64 (c, 17);
  int_add (a, b, c);

  int_set_i64 (b, -7);
  int_set_i64 (c, 17);
  int_add (b, b, c);

  int_set_i64 (b, -17);
  int_set_i64 (c, 7);
  int_add (c, b, c);

  int_set_i64 (b, 17);
  int_set_i64 (c, -7);
  int_add (b, b, c);

  int_set_i64 (b, 7);
  int_set_i64 (c, -17);
  int_add (c, b, c);

  int_set_i64 (c, -17);
  int_add (c, c, c);

  int_set_i64 (c, 17);
  int_add (c, c, c);

  int_clear (a);
  int_clear (b);
  int_clear (c);
  return 0;
}
