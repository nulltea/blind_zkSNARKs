/**
 * Query lazer version information.
 */
#include "lazer.h"
#include <stdio.h>

int
main (void)
{
  /* print header version */
  printf ("lazer.h version: %s\n", LAZER_VERSION);

  /* print library version */
  printf ("liblazer.a version: %s\n", lazer_get_version ());

  /* check if versions match */
  if (LAZER_VERSION_MAJOR != lazer_get_version_major ()
      || LAZER_VERSION_MINOR != lazer_get_version_minor ()
      || LAZER_VERSION_PATCH != lazer_get_version_patch ())
    {
      printf ("header and library versions do not match.\n");
      return 1;
    }
  return 0;
}
