//<require "prep.pl"//>//
//<=$warn_dont_modify //>

//__INSERT_LICENSE__
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <mkl_cblas.h>

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fastlib2.h>

int FASTMAT2_USE_DGEMM=1;
int FASTMAT2_PROD_WAS_SUPERLINEAR=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class prod_subcache : public FastMatSubCache {
public:
  int superlinear, lda, ldb, ldc, nra, nca, ncb;
  CBLAS_TRANSPOSE transa, transb;
  double *paa0, *pbb0, *pcc0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Helper class, masks the relation betewen index
// and start of `LineCache' or `target' for the function
// `is_superl()' below. 
class superl_helper_t {
public:
  enum  mode_t { a,b,c,none };
  FastMatCache *cache;
  LineCache *line_cache_start;
  superl_helper_t(FastMatCache *cache_a)
    : cache(cache_a), 
      line_cache_start(cache->line_cache_start) { }
  double *address(int j, mode_t mode) {
    if (mode==a) 
      return *line_cache_start[j].starta;
    else if (mode==b) 
      return *line_cache_start[j].startb;
    else if (mode==c) 
      return line_cache_start[j].target;
    else assert(0);
    return NULL;
  }
  int has_rmo(int &m,int &n,int &incrow, int &inccol);
#define MODE_SHORT 0
#define MODE_LONG 1
  int has_rmo2(mode_t mode,
               int &nrow,int &ncol,
               int &incrow, int &inccol,
               int &byrows, int &trvmode);
  int has_rmo3(int &nrow,int &ncol);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int superl_helper_t
::has_rmo(int &m,int &n,int &incrow, int &inccol) {

  LineCache *line_cache_start=NULL;
  int N = cache->nlines;
  if (N<1) return 0;
  
  mode_t mode=c;
  double *aa00 = address(0,mode);
  inccol = 1;
  if (N>1) inccol = address(1,mode) - aa00; 
  int dp, j;
  for (j=2; j<N; j++) {
    dp = address(j,mode) - aa00;
    if (dp != j*inccol) break;
  }
  n=j;

  m = N/n;
  incrow = dp;
  if (m==1) incrow = dp+1;
  int l = 0;
  for (int j=0; j<m; j++) {
    for (int k=0; k<n; k++) {
      dp = address(l,mode) - aa00; 
      if (dp != incrow*j + inccol*k)
        return 0;
      l++;
    }
  }

  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Determines whether a given set of addresses &a(j), with j
// in [0,N) are in RMO (Rectangular Matrix Order). i.e. there
// exists m,n such that if j = k*n+l and
//
//   &aa(k,l) = &aa(0,0) + inccol*k + rowcol*l
int superl_helper_t
::has_rmo2(mode_t mode,int &nrow,int &ncol,
           int &incrow, int &inccol,
           int &byrows, int &trvmode) {

  int inc1=0, inc2=0, size1, size2;
  int N = cache->nlines;

  double *aa00 = address(0,mode);
  inc1 = 1;
  if (N>1) inc1 = address(1,mode) - aa00; 
  inccol = (mode==a? cache->line_cache_start->inca :
            mode==b? cache->line_cache_start->incb : -1);
  assert(inccol!=-1);

  ncol = cache->line_size;

  int dp, j;
  for (j=2; j<N; j++) {
    dp = address(j,mode) - aa00;
    if (dp != j*inc1) break;
  }

  size2=j;
  if (N%size2 != 0) return 0;
  size1 = N/size2;
  inc2 = dp;
  if (size1==1) inc2 = !dp;
  // if (size1==1 && inc1==0) inc2 = size2*;
  int l = 0;
  for (int j=0; j<size1; j++) {
    for (int k=0; k<size2; k++) {
      dp = address(l,mode) - aa00; 
      if (dp != inc2*j + inc1*k)
        return 0;
      l++;
    }
  }
  
  // It is a matrix
  assert(inc1==0 || inc2==0);
  
  byrows = 1;
  if (inc1==0) {
    trvmode = MODE_SHORT;
    incrow = inc2;
    nrow = size1;
  } else {
    trvmode = MODE_LONG;
    incrow = inc1;
    nrow = size2;
  }

  if (inccol!=1) {
    if (incrow!=1) return 0;
    else {
      // swap rows with cols
      swap(inccol,incrow);
      swap(ncol,nrow);
      byrows = !byrows;
    }
  }

#if 0  
  printf("nrow %d, ncol %d, incrow %d, inccol %d, trvmode %s, "
         "byrows %d\n",nrow,ncol,incrow,inccol,
         (trvmode==MODE_SHORT? "short" : "long"),byrows);
#endif

  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int superl_helper_t
::has_rmo3(int &nrow,int &ncol) {

  int inc1=0, inc2=0, size1, size2;
  int N = cache->nlines;

  if (N<=2) return 0;
  double *aa00 = address(0,superl_helper_t::a);
  double *bb00 = address(0,superl_helper_t::b);
  double *cc00 = address(0,superl_helper_t::c);
  
  int inc1a = address(1,superl_helper_t::a) - aa00; 
  int inc1b = address(1,superl_helper_t::b) - bb00; 
  int inc1c = address(1,superl_helper_t::c) - cc00; 

  int j, inc2a, inc2b, inc2c;
  int nrowopa = -1, ncolopb;
  for (j=1; j<N; j++) 
    if (address(j,superl_helper_t::a) - aa00) break;
  ncolopb=j;
  if (N % ncolopb != 0) return 0;
  nrowopa = N/ncolopb;

  LineCache *lc0 = cache->line_cache_start;
  int 
    inca = lc0->inca,
    lda = address(ncolopb,superl_helper_t::a) - aa00,
    incb = address(1,superl_helper_t::b) - bb00,
    ldb = lc0->incb,
    incc = address(1,superl_helper_t::c) - cc00,
    ldc = address(ncolopb,superl_helper_t::c) - cc00;
  int dp, l=0;
  for (j=0; j<nrowopa; j++) {
    for (int k=0; k<ncolopb; k++) {
      dp = address(l,superl_helper_t::a) - aa00; 
      if (dp != lda*k) return 0;

      dp = address(l,superl_helper_t::b) - bb00; 
      if (dp != incb*j) return 0;

      dp = address(l,superl_helper_t::c) - cc00; 
      if (dp != incc*j + ldc*k) return 0;

      l++;
    }
  }

  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//<$prod = <<'//EOF';
//EOF
_//>

//< print template_subst($prod); //>//

//<=$warn_dont_modify //>
