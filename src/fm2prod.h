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
  int print_prod_order;
  int ncall, ncall_not_cached;
  double tcall, tcall_not_cached;
  double tpart;
  void reset() {
    ncall=0;
    ncall_not_cached=0;
    tcall=0.0;
    tcall_not_cached=0.0;
    tpart=0.0;
  }
  void print();
  fastmat_stats_t() : print_prod_order(0) { reset(); }
};

extern fastmat_stats_t fastmat_stats;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// We keep a set of this pairs for the active matrices.
// Each pair contains the key (an index in `mat_info_cont'
// and the position in the active matrix set).
// This last one should be irrelevant and is kept only
// in order to make the computation more readable.
// Also could have an impact in performance because
// the product may be done with `dgemm' or not. 
class active_mat_info_t {
public:
  int key,position;
  // Comparison function
  bool operator<(const active_mat_info_t y) const {
    if (position != y.position) 
      return position < y.position;
    else 
      return key < y.key;
  }
  active_mat_info_t(int k,int p) 
  : key(k), position(p) { }
};

typedef std::set<active_mat_info_t> active_mat_info_cont_t;

#endif
