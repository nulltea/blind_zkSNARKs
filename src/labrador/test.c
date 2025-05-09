#include <stdint.h>
#include <stdio.h>

const uint16_t p = 16001;
const uint16_t pinv = 4*385;

int main(void) {
  uint16_t i,j;
  uint16_t m;
  int16_t r,min,max;

  min = INT16_MAX;
  max = INT16_MIN;
  for(i=0;i<(1 << 14);i++) {
    for(j=0;j<(1 << 14);j++) {
      m = i*pinv;  // 4m
      r = j - ((uint32_t)m*p >> 16);
      if(r < min) min = r;
      if(r > max) max = r;

      if((((int32_t)i << 0) + ((int32_t)j << 14) - ((int32_t)r << 14)) % p) printf("FAILURE\n");
    }
  }

  printf("%d..%d\n",min,max);
  return 0;
}
