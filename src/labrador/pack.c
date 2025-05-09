#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "malloc.h"
#include "labrador.h"
#include "chihuahua.h"
#include "dachshund.h"
#include "greyhound.h"
#include "pack.h"

void free_composite(composite *p) {
  size_t i;

  free_witness(&p->owt);
  for(i=0;i<p->l;i++) {
    free_proof(p->pi[i]);
    free(p->pi[i]);
  }
}

static void composite_prove(composite *p, statement *tst, witness *twt, double *twtsize) {
  size_t i = 0;
  double pisize;

  while(p->l < 16) {
    p->pi[p->l] = _malloc(sizeof(proof));
    prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],0);
    pisize = print_proof_pp(p->pi[p->l]);
    print_statement_pp(&tst[i^1]);
    twtsize[i^1] = print_witness_pp(&twt[i^1]);
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      free_proof(p->pi[p->l]);
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
      break;
    }

    free_statement(&tst[i]);
    free_witness(&twt[i]);
    p->size += pisize;
    p->l += 1;
    i ^= 1;
  }

  if(p->l < 16) {
    prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],1);
    pisize = print_proof_pp(p->pi[p->l]);
    print_statement_pp(&tst[i^1]);
    twtsize[i^1] = print_witness_pp(&twt[i^1]);
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      free_proof(p->pi[p->l]);
      free(p->pi[p->l]);
    }
    else {
      p->size += pisize;
      p->l += 1;
      i ^= 1;
    }
  }

  free_statement(&tst[i^1]);
  free_witness(&twt[i^1]);
  free_statement(&tst[i]);
  p->owt = twt[i];
  p->size += twtsize[i];
}

double composite_prove_principle(composite *p, const prncplstmnt *st, const witness *wt) {
  statement tst[2];
  witness twt[2];
  double twtsize[2];
  clock_t t;

  t = clock();
  p->pi[0] = _malloc(sizeof(proof));
  principle_prove(tst,twt,p->pi[0],st,wt,0);
  p->size = print_proof_pp(p->pi[0]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  p->l = 1;

  composite_prove(p,tst,twt,twtsize);
  t = clock() - t;
  printf("Commitment key length: %zu\n",comkey_len);
  printf("Chihuahua Pack members: %zu\n",p->l);
  printf("Chihuahua Pack size: %.2f KB\n",p->size);
  printf("Chihuahua Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return p->size;
}

double composite_prove_simple(composite *p, commitment *com, const smplstmnt *st, const witness *wt) {
  statement tst[2];
  witness twt[2];
  double twtsize[2];
  clock_t t;

  t = clock();
  p->pi[0] = _malloc(sizeof(proof));
  simple_prove(tst,twt,p->pi[0],com,st,wt,0);
  p->size = print_proof_pp(p->pi[0]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  p->l = 1;

  composite_prove(p,tst,twt,twtsize);
  t = clock() - t;
  printf("Commitment key length: %zu\n",comkey_len);
  printf("Dachshund Pack members: %zu\n",p->l);
  printf("Dachshund Pack size: %.2f KB\n",p->size);
  printf("Dachshund Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return p->size;
}

double composite_prove_polcom(composite *p, polcomprf *ppi, polcomctx *ctx, uint32_t x, uint32_t y) {
  prncplstmnt tst0[1];
  statement tst[2];
  witness twt[2];
  double twtsize[2];
  clock_t t;

  t = clock();
  polcom_eval(&twt[1],ppi,ctx,x,y);
  polcom_reduce(tst0,ppi);

  p->pi[0] = _malloc(sizeof(proof));
  principle_prove(tst,twt,p->pi[0],tst0,&twt[1],0);
  p->size = print_proof_pp(p->pi[0]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  free_prncplstmnt(tst0);
  free_witness(&twt[1]);
  p->l = 1;

  composite_prove(p,tst,twt,twtsize);
  t = clock() - t;
  printf("Commitment key length: %zu\n",comkey_len);
  printf("Greyhound Pack members: %zu\n",p->l);
  printf("Greyhound Pack size: %.2f KB\n",p->size);
  printf("Greyhound Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return p->size;
}

static int composite_verify(const composite *p, statement *tst) {
  size_t i,j;
  int ret = 0;

  i = 0;
  for(j=1;j<p->l;j++) {
    ret = reduce(&tst[i^1],p->pi[j],&tst[i]);
    free_statement(&tst[i]);
    i ^= 1;
    if(ret) {
      printf("ERROR: reduce to member %zu failed: %d\n",j,ret);
      goto end;
    }
  }

  ret = verify(&tst[i],&p->owt);
  if(ret) {
    printf("ERROR: verification of final statement failed: %d\n",ret);
    goto end;
  }

end:
  free_statement(&tst[i]);
  return ret;
}

int composite_verify_principle(const composite *p, const prncplstmnt *st) {
  int ret;
  statement tst[2];
  clock_t t;

  t = clock();
  ret = principle_reduce(tst,p->pi[0],st);
  if(ret) {
    printf("ERROR: principle_reduce failed: %d\n",ret);
    free_statement(tst);
    return ret;
  }

  ret = composite_verify(p,tst);
  t = clock() - t;
  printf("Chihuahua Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return ret;
}

int composite_verify_simple(const composite *p, const commitment *com, const smplstmnt *st) {
  int ret;
  statement tst[2];
  clock_t t;

  t = clock();
  ret = simple_reduce(tst,p->pi[0],com,st);
  if(ret) {
    printf("ERROR: simple_reduce failed: %d\n",ret);
    free_statement(tst);
    return ret;
  }

  ret = composite_verify(p,tst);
  t = clock() - t;
  printf("Dachshund Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return ret;
}

int composite_verify_polcom(const composite *p, const polcomprf *ppi) {
  int ret;
  prncplstmnt tst0[1];
  statement tst[2];
  clock_t t;

  t = clock();
  polcom_reduce(tst0,ppi);

  ret = principle_reduce(tst,p->pi[0],tst0);
  free_prncplstmnt(tst0);
  if(ret) {
    printf("ERROR: principle_reduce failed: %d\n",ret);
    free_statement(tst);
    return ret;
  }

  ret = composite_verify(p,tst);
  t = clock() - t;
  printf("Greyhound Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return ret;
}
