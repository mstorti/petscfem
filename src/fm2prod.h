// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id config-1.0.7-74-gc7a455f Mon Oct 29 23:28:05 2007 -0300$
#ifndef PETSCFEM_FM2PROD_H
#define PETSCFEM_FM2PROD_H

extern FastMat2Stats glob_fm2stats;
extern int FASTMAT2_USE_DGEMM;
#define USE_MPROD_FOR_2MATS

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class prod_subcache_t : public FastMatSubCache {
public:
  int superlinear, lda, ldb, ldc, nra, nca, ncb,
    not_superlinear_ok_flag;
#ifdef USE_MKL
  CBLAS_TRANSPOSE transa, transb;
#endif
  double *paa0, *pbb0, *pcc0;
  enum  mode_t { a,b,c,none };
  FastMatCache *cache;
  LineCache *line_cache_start;
  prod_subcache_t(FastMatCache *cache_a)
    : superlinear(0),
      cache(cache_a), 
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
  void ident();
  void dgemm();
  void print();
  void ok() { not_superlinear_ok_flag=1; }
  int not_superlinear_ok();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
struct fastmat_stats_t {
  int ncall, ncall_not_cached;
  double tcall, tcall_not_cached;
  void reset() {
    ncall=0;
    ncall_not_cached=0;
    tcall=0.0;
    tcall_not_cached=0.0;
  }
  void print();
  fastmat_stats_t() { reset(); }
};

extern fastmat_stats_t fastmat_stats;

#endif
