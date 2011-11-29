// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id config-1.0.7-74-gc7a455f Mon Oct 29 23:28:05 2007 -0300$
#ifndef PETSCFEM_FM2PROD_H
#define PETSCFEM_FM2PROD_H

// #define CALL_BLAS_LAPACK_FROM_PETSC
#ifdef CALL_BLAS_LAPACK_FROM_PETSC
#include <petsc.h>
#include <petscblaslapack.h>
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                      AtlasConj=114};
#else
#ifdef USE_MKL
#include <mkl_cblas.h>
#else
extern "C" {
#include <cblas.h>
}
#endif
#endif

extern FastMat2Stats glob_fm2stats;
extern int FASTMAT2_USE_DGEMM;
extern int FASTMAT2_USE_PROD2;
extern const int nmax;
// If not set: for product of 2 matrices a special code is used
// avoiding the overhead of the multiprod algorithm
// Currently it can't be deactivated
#define USE_MPROD_FOR_2MATS

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// We store a vector of these structures in the cache, and
// then `make_prod()' has the stored the order
// in which the products must be done, and which matrices
// are involved.
// The `mat_info's are stored for both the
// original matrices, and the temporaries that
// are created. So if we have `n' matrices we make
// `n-1' products, and we need a mat_info for the result
// so we have `2*n-1' mat_infos. 
// From the `n-1' results, the first `n-2' are temporaries,
// and the last one goes to the final output result. 
struct mat_info {
  // Pointers to old matrices should be `const'
  FastMat2 *Ap;
  // The vector that indicates the contractions
  // to be performed (the args to the low-level
  // prod())
  vector<int> contract;
  // The dims of the involved matrices
  vector<int> dims;
  // type: may be OLD, TMP or UNKNOWN
  // is_active: when we make a product the two involved
  //            matrices are marked as INACTIVE and the new
  //            inserted to ACTIVE
  // position: stores the position in the matrix list.
  //           We try to preserve the position so that
  //           the process is more clear. 
  int type, is_active, position;
  mat_info();
};

typedef vector<mat_info> mat_info_cont_t;

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
// This is the cache for the product operation
class multiprod_subcache_t : public FastMatSubCache {
public:
  enum { OLD=0, TMP=1, UNKNOWN=2 };
  enum { INACTIVE=0, ACTIVE=1, UNDEF=-1 };
  // Number of matrices involved in this product.
  // Must be >=2
  int nmat;
  // A vector of structures containing information
  // for each involved matrix (including temporaries)
  mat_info_cont_t mat_info_cont;
  // A table that stores in which orders must peformed
  // the products
  vector<int> order;
  multiprod_subcache_t(FastMatCache *cache_a) { }
  ~multiprod_subcache_t();
  // This makes the product when cached
  void make_prod();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This class stores the information for the wrapper
// function prod().
class mprodwrp_subcache_t : public FastMatSubCache {
public:
  // The list of matrices involved
  vector<const FastMat2 *> mat_list;
  // The vector of contraction indices passed
  Indx indx;
  mprodwrp_subcache_t() { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
struct fastmat_stats_t {
  int print_prod_order;
  int ncall, ncall_not_cached;
  double tcall, tcall_not_cached;
  double tpart;
  const char **labels;
  int nmat;
  void reset() {
    ncall=0;
    ncall_not_cached=0;
    tcall=0.0;
    tcall_not_cached=0.0;
    tpart=0.0;
  }
  void print();
  fastmat_stats_t() 
  : print_prod_order(0), 
    labels(NULL),
    nmat(-1)
  { reset(); }
};

extern fastmat_stats_t fastmat_stats;

#if 0
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
typedef std::set<active_mat_info_t> active_mat_info_cont_t;
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int compute_opcount(const mat_info &qmi,
                    const mat_info &rmi,
                    mat_info &smi,
                    int &qfree,int &rfree,int &qctr);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
intmax_t 
compute_optimal_order(const mat_info_cont_t &mat_info_cont,
                      vector<int> &order);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
intmax_t 
compute_heuristic_order(const mat_info_cont_t &mat_info_cont,
                        vector<int> &order);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
intmax_t 
compute_natural_order(const mat_info_cont_t &mat_info_cont,
                      vector<int> &order,int reverse=0);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class prod2_subcache_t : public FastMatSubCache {
public:

  // ----------- FMGEMM STUFF ----------------
  // If set: do some statistics about how many computations
  // call the FMGEMM functions
  // #define DO_SIZE_STATS
  // Utility macro for simplifiying the declaration of the FMGEMM functions
  // in `mygmdefs.h'
#define DECLFUN(fun)                                                    \
  static void fun(double *__restrict__ a,double * __restrict__ b,double * __restrict__ c)
  // This header contains the definitions of all the FMGEMM functions 
  // (there are NMAX*NMAX*NMAX*4 such functions
#include "./mygmdefs.h"  
  // nmax is a variable that can be changed at runtime by the user, so
  // regulating the size of the matrices for which the FMGEMM functions
  // are called
  // NMAX stores the maximum order for the FMGEMM matrices that have
  // been compiled, so nmax must be nmax<=NMAX
  // NMAX is set by the $nmax variable in the `mkmydgcode.pl' Perl script
  static int nmax, NMAX;
  // This is the typical FMGEMM function that is called
  // pointers to this functions are saved in a vector.
  // The pointer to the function is saved in the cache
  typedef void (*gemm_fun_t)(double *a,double *b,double *c);
  // Pointers to the functions are stored in this table.
  // Later this table is used as a multiarray vector of shape
  // [NMAX,NMAX,NMAX,2,2]
  static vector<gemm_fun_t> gemm_fun_table;
  // This is a static flag that indicates whether this table
  // has been initialiazed or not
  static int gemm_fun_table_was_initialized;
  // If set to 0 then FMGEMM are not called. This can be changed
  // at runtime
  static int FASTMAT2_USE_FMGEMM;

#ifdef DO_SIZE_STATS
  struct mat_sz_t { 
    int m,n,p; 
    bool operator<(const mat_sz_t &rhs) const;
  };
  struct stats_t { 
    int count; 
    double time;
    stats_t() : count(0), time(0.0) { }
  };
  static map<mat_sz_t,stats_t> stat_table;
  static void report_stats();
  typedef pair< mat_sz_t,stats_t> spair_t;
  static bool comp(const spair_t &a,const spair_t &b);
  static int do_size_stats, // choose at runtime whether to do stats or not
    total_calls,            // total FastMat2 prod calls
    fmgemm_calls,           // FastMat2 prod calls that use FMGEMM
    last_call_used_fmgemm;  // flags where the last call used FMGEMM or not
#endif

  static int gemm_fun_table_indx(int n,int m,int p,int jat,int jbt);
  static gemm_fun_t get_fun(int n,int m,int p,int jat,int jbt);
  static void gemm_fun_table_load(int n,int m,int p,int jat,int jbt,gemm_fun_t f);
  static void init_funs();
  static void load_funs();

  vector<double *> ap,bp,cp;  
  double *Ap,*Bp,*Cp;
  vector<double> a,b,c;
  int nA,nB,nC, 
    nrowa,ncola,
    nrowb,ncolb,
    nrowc,ncolc;
  int asl_ok,lda,
    bsl_ok,ldb,
    csl_ok,ldc;
  CBLAS_TRANSPOSE transa,transb,transc;
  int use_fmgemm;
  gemm_fun_t gfun;
  void init(const FastMat2 &A,const FastMat2 &B,
            FastMat2 &C,
            vector<int> &ixa, 
            vector<int> &ixb);
  void make_prod();
  prod2_subcache_t() : gfun(NULL) {}
};

#endif
