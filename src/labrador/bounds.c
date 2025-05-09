#include <stdint.h>
#include <stdio.h>
#include "data.h"

int main(void) {
  int32_t i;
  int16_t tx,minx,maxx;
  int32_t tz,minz,maxz;
  const pdata *prime = &primes[5];
  const int32_t bound = 3*(1 << 29) - 1;

  minx = INT16_MAX;
  maxx = INT16_MIN;
  minz = INT32_MAX;
  maxz = INT32_MIN;
  for(i=-bound;i<=bound;i++) {
    tx = i*prime->pinv;
    tx = (i - (int32_t)tx*prime->p) >> 16;
    if(tx < minx) minx = tx;
    if(tx > maxx) maxx = tx;

    tz = i + (1 << 29);
    tz = tz >> 30;
    tz = i - Q*tz;
    if(tz < minz) minz = tz;
    if(tz > maxz) maxz = tz;
  }

  printf("x: %d..%d\n",minx,maxx);
  printf("z: %d..%d\n",minz,maxz);

  return 0;
}
