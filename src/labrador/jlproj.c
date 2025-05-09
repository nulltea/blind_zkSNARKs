#include <stdint.h>
#include <immintrin.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"
#include "polz.h"
#include "jlproj.h"

void poly_jlproj_add(int32_t r[256], const poly *p, const uint8_t mat[256*N/8]) {
  int i;
  __m512i a,b,c,d;
  __m512i e,f,g,h;
  __m512i k,l,m,n;
  __m512i s,t,u,v;
  __m512i w,x,y,z;
  const __m512i zero = _mm512_setzero_si512();

  a = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&p->vec->c[32]));
  b = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&p->vec->c[48]));

  e = _mm512_add_epi32(b,a);     // ++ ( b+a)
  f = _mm512_sub_epi32(b,a);     // +- ( b-a)
  g = _mm512_sub_epi32(a,b);     // -+ (-b+a)
  h = _mm512_sub_epi32(zero,e);  // -- (-b-a)

  a = _mm512_shuffle_i64x2(e,f,0x44);  // ++,+-
  b = _mm512_shuffle_i64x2(e,f,0xEE);  // ++,+-
  c = _mm512_shuffle_i64x2(g,h,0x44);  // -+,--
  d = _mm512_shuffle_i64x2(g,h,0xEE);  // -+,--

  e = _mm512_shuffle_i64x2(a,c,0x88);  // ++,+-,-+,--
  f = _mm512_shuffle_i64x2(a,c,0xDD);  // ++,+-,-+,--
  g = _mm512_shuffle_i64x2(b,d,0x88);  // ++,+-,-+,--
  h = _mm512_shuffle_i64x2(b,d,0xDD);  // ++,+-,-+,--

  a = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&p->vec->c[ 0]));
  b = _mm512_cvtepi16_epi32(_mm256_load_si256((__m256i*)&p->vec->c[16]));

  c = _mm512_add_epi32(b,a);  // ++ (b+a)
  d = _mm512_sub_epi32(b,a);  // +- (b-a)

  k = _mm512_shuffle_i64x2(c,c,0x00);  // ++
  l = _mm512_shuffle_i64x2(d,d,0x00);  // +-
  m = _mm512_shuffle_i64x2(c,c,0x55);
  n = _mm512_shuffle_i64x2(d,d,0x55);
  s = _mm512_shuffle_i64x2(c,c,0xAA);
  t = _mm512_shuffle_i64x2(d,d,0xAA);
  u = _mm512_shuffle_i64x2(c,c,0xFF);
  v = _mm512_shuffle_i64x2(d,d,0xFF);

  w = _mm512_add_epi32(e,k);  // ++++,+-++,-+++,--++
  x = _mm512_sub_epi32(e,k);  // ++--,+---,-+--,----
  y = _mm512_add_epi32(e,l);  // +++-,+-+-,-++-,--+-
  z = _mm512_sub_epi32(e,l);  // ++-+,+--+,-+-+,---+
  a = _mm512_add_epi32(f,m);  // ++++,+-++,-+++,--++
  b = _mm512_sub_epi32(f,m);  // ++--,+---,-+--,----
  c = _mm512_add_epi32(f,n);  // +++-,+-+-,-++-,--+-
  d = _mm512_sub_epi32(f,n);  // ++-+,+--+,-+-+,---+
  k = _mm512_add_epi32(g,s);  // ++++,+-++,-+++,--++
  l = _mm512_sub_epi32(g,s);  // ++--,+---,-+--,----
  m = _mm512_add_epi32(g,t);  // +++-,+-+-,-++-,--+-
  n = _mm512_sub_epi32(g,t);  // ++-+,+--+,-+-+,---+
  e = _mm512_add_epi32(h,u);  // ++++,+-++,-+++,--++
  f = _mm512_sub_epi32(h,u);  // ++--,+---,-+--,----
  g = _mm512_add_epi32(h,v);  // +++-,+-+-,-++-,--+-
  h = _mm512_sub_epi32(h,v);  // ++-+,+--+,-+-+,---+

  s = _mm512_unpacklo_epi32(w,y);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  t = _mm512_unpackhi_epi32(w,y);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  u = _mm512_unpacklo_epi32(z,x);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  v = _mm512_unpackhi_epi32(z,x);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  w = _mm512_unpacklo_epi32(a,c);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  x = _mm512_unpackhi_epi32(a,c);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  y = _mm512_unpacklo_epi32(d,b);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  z = _mm512_unpackhi_epi32(d,b);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  a = _mm512_unpacklo_epi32(k,m);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  b = _mm512_unpackhi_epi32(k,m);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  c = _mm512_unpacklo_epi32(n,l);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  d = _mm512_unpackhi_epi32(n,l);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  k = _mm512_unpacklo_epi32(e,g);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  l = _mm512_unpackhi_epi32(e,g);  // ++++,+++-,+-++,+-+-,-+++,-++-,--++,--+-
  m = _mm512_unpacklo_epi32(h,f);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----
  n = _mm512_unpackhi_epi32(h,f);  // ++-+,++--,+--+,+---,-+-+,-+--,---+,----

  e = _mm512_unpacklo_epi64(s,u);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  f = _mm512_unpackhi_epi64(s,u);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  g = _mm512_unpacklo_epi64(t,v);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  h = _mm512_unpackhi_epi64(t,v);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  s = _mm512_unpacklo_epi64(w,y);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  t = _mm512_unpackhi_epi64(w,y);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  u = _mm512_unpacklo_epi64(x,z);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  v = _mm512_unpackhi_epi64(x,z);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  w = _mm512_unpacklo_epi64(a,c);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  x = _mm512_unpackhi_epi64(a,c);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  y = _mm512_unpacklo_epi64(b,d);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  z = _mm512_unpackhi_epi64(b,d);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  a = _mm512_unpacklo_epi64(k,m);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  b = _mm512_unpackhi_epi64(k,m);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  c = _mm512_unpacklo_epi64(l,n);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----
  d = _mm512_unpackhi_epi64(l,n);  // ++++,+++-,++-+,++--,+-++,+-+-,+--+,+---,-+++,-++-,-+-+,-+--,--++,--+-,---+,----

  for(i=0;i<256/16;i++) {
    k = _mm512_loadu_si512((__m512i*)&mat[128*i+ 0]);
    l = _mm512_loadu_si512((__m512i*)&r[16*i]);
    m = _mm512_permutexvar_epi32(k,e);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,f);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,g);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,h);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,s);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,t);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,u);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,v);
    l = _mm512_add_epi32(l,m);
    k = _mm512_loadu_si512((__m512i*)&mat[128*i+64]);
    m = _mm512_permutexvar_epi32(k,w);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,x);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,y);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,z);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,a);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,b);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,c);
    l = _mm512_add_epi32(l,m);
    k = _mm512_srli_epi32(k,4);
    m = _mm512_permutexvar_epi32(k,d);
    l = _mm512_add_epi32(l,m);
    _mm512_storeu_si512((__m512i*)&r[16*i],l);
  }
}

void polyvec_jlproj_add(int32_t r[256], const poly *p, size_t len, const uint8_t *mat) {
  size_t i;

  for(i=0;i<len;i++)
    poly_jlproj_add(r,&p[i],&mat[256*N/8*i]);
}

static inline __m512i _mm512_hsum_epi32(__m512i a) {
  __m256i t;
  __m128i u;

  t = _mm256_add_epi32(_mm512_castsi512_si256(a),_mm512_extracti64x4_epi64(a,1));
  u = _mm_add_epi32(_mm256_castsi256_si128(t),_mm256_extracti64x2_epi64(t,1));
  u = _mm_add_epi32(u,_mm_unpackhi_epi64(u,u));
  u = _mm_add_epi32(u,_mm_srli_epi64(u,32));
  return _mm512_broadcastd_epi32(u);
}

void polz_jlproj_collapsmat24(polz *r, const uint8_t mat[256*N/8], const uint8_t buf[3*256]) {
  int i,j,k;
  __m512i f[16],g,sum;
  __m512i s,t,u,v;
  __mmask16 mask0,mask1;
  const __mmask32 alternate = _cvtu32_mask32(0x55555555);
  const __m512i mask14 = _mm512_set1_epi32(0x3FFF);
  const __m512i vpermbidx = _mm512_set_epi8(-1,47,46,45,-1,44,43,42,-1,41,40,39,-1,38,37,36,
                                            -1,35,34,33,-1,32,31,30,-1,29,28,27,-1,26,25,24,
                                            -1,23,22,21,-1,20,19,18,-1,17,16,15,-1,14,13,12,
                                            -1,11,10, 9,-1, 8, 7, 6,-1, 5, 4, 3,-1, 2, 1, 0);
  const __mmask64 vpermbmask = _knot_mask64(_mm512_movepi8_mask(vpermbidx));

  s = _mm512_loadu_si512(&buf[0]);
  s = _mm512_maskz_permutexvar_epi8(vpermbmask,vpermbidx,s);
  for(i=1;i<256/16;i++) {
    t = _mm512_loadu_si512(&buf[48*i]);
    t = _mm512_maskz_permutexvar_epi8(vpermbmask,vpermbidx,t);
    s = _mm512_add_epi32(s,t);
  }
  sum = _mm512_hsum_epi32(s);

  for(j=0;j<N/16;j++) {
    for(k=0;k<16;k++)
      f[k] = _mm512_setzero_si512();

    for(i=0;i<256/16;i++) {
      g = _mm512_loadu_si512(&buf[48*i]);
      g = _mm512_maskz_permutexvar_epi8(vpermbmask,vpermbidx,g);

      s = _mm512_loadu_si512((__m512i*)&mat[128*i+ 0]);
      t = _mm512_loadu_si512((__m512i*)&mat[128*i+64]);
      for(k=0;k<8;k++) {
        u = _mm512_slli_epi32(s,31-(4*k+j));
        v = _mm512_slli_epi32(t,31-(4*k+j));
        mask0 = _mm512_movepi32_mask(u);
        mask1 = _mm512_movepi32_mask(v);
        f[k] = _mm512_mask_add_epi32(f[k],mask0,f[k],g);
        f[8+k] = _mm512_mask_add_epi32(f[8+k],mask1,f[8+k],g);
      }
    }

    for(k=0;k<8;k++) {
      u = _mm512_unpacklo_epi32(f[2*k+0],f[2*k+1]);
      v = _mm512_unpackhi_epi32(f[2*k+0],f[2*k+1]);
      f[k] = _mm512_add_epi32(u,v);
    }
    for(k=0;k<4;k++) {
      u = _mm512_unpacklo_epi64(f[2*k+0],f[2*k+1]);
      v = _mm512_unpackhi_epi64(f[2*k+0],f[2*k+1]);
      f[k] = _mm512_add_epi32(u,v);
    }
    for(k=0;k<2;k++) {
      u = _mm512_shuffle_i64x2(f[2*k+0],f[2*k+1],0x44);
      v = _mm512_shuffle_i64x2(f[2*k+0],f[2*k+1],0xEE);
      f[k] = _mm512_add_epi32(u,v);
    }
    u = _mm512_shuffle_i64x2(f[0],f[1],0x88);
    v = _mm512_shuffle_i64x2(f[0],f[1],0xDD);
    u = _mm512_add_epi32(u,v);
    v = _mm512_sub_epi32(sum,u);
    f[0] = _mm512_sub_epi32(v,u);

    u = _mm512_and_si512(f[0],mask14);
    _mm512_mask_compressstoreu_epi16(&r->limbs[0].c[16*j],alternate,u);
    f[0] = _mm512_srai_epi32(f[0],14);
    u = _mm512_and_si512(f[0],mask14);
    _mm512_mask_compressstoreu_epi16(&r->limbs[1].c[16*j],alternate,u);
    u = _mm512_srai_epi32(f[0],14);
    _mm512_mask_compressstoreu_epi16(&r->limbs[2].c[16*j],alternate,u);
  }
}

void polxvec_jlproj_collapsmat24(polx *r, const uint8_t *mat, size_t len, const uint8_t buf[3*256]) {
  size_t i;
  polz t[16];

  while(len >= 16) {
    for(i=0;i<16;i++)
      polz_jlproj_collapsmat24(&t[i],&mat[i*256*N/8],buf);
    polzvec_sigmam1(t,t,16);
    polzvec_topolxvec(r,t,16);
    r += 16;
    mat += 16*256*N/8;
    len -= 16;
  }

  for(i=0;i<len;i++)
    polz_jlproj_collapsmat24(&t[i],&mat[i*256*N/8],buf);
  polzvec_sigmam1(t,t,len);
  polzvec_topolxvec(r,t,len);
}

void jlproj_collapsproj24(int64_t *r,const int32_t p[256],const uint8_t buf[3*256]) {
  size_t i;
  int64_t u,v;

  // TODO: Vectorize
  u = 0;
  for(i=0;i<256;i++) {
    v  = buf[3*i+0];
    v |= (int64_t)buf[3*i+1] <<  8;
    v |= (int64_t)buf[3*i+2] << 16;
    u += v*p[i];
  }

  *r = u;
}

uint64_t jlproj_normsq(const int32_t p[256]) {
  size_t i;
  __m512i f,g,s;

  s = _mm512_setzero_si512();
  for(i=0;i<256/16;i++) {
    f = _mm512_loadu_si512(&p[16*i]);
    g = _mm512_mul_epi32(f,f);
    f = _mm512_srli_epi64(f,32);
    f = _mm512_mul_epi32(f,f);
    f = _mm512_add_epi64(f,g);
    s = _mm512_add_epi64(f,s);
  }

  __m256i t;
  __m128i u;
  t = _mm256_add_epi64(_mm512_castsi512_si256(s),_mm512_extracti64x4_epi64(s,1));
  u = _mm_add_epi64(_mm256_castsi256_si128(t),_mm256_extracti64x2_epi64(t,1));
  u = _mm_add_epi64(u,_mm_unpackhi_epi64(u,u));
  return _mm_extract_epi64(u,1);
}
