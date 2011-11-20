//__INSERT_LICENSE__
#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fm2stats.h>
#include <src/fastlib2.h>
#include <src/fm2prod.h>

int prod2_subcache_t::NMAX=PF_MYDGEMM_NMAX;
int prod2_subcache_t::nmax=PF_MYDGEMM_NMAX;
int prod2_subcache_t::gemm_fun_table_was_initialized=0;
int prod2_subcache_t::FASTMAT2_USE_MYDGEMM=1;
vector<prod2_subcache_t::gemm_fun_t> prod2_subcache_t::gemm_fun_table;

void prod2_subcache_t::load_funs() {
  gemm_fun_table.resize(NMAX*NMAX*NMAX);
#define LOADFUN(n,m,p) gemm_fun_table_load(n,m,p,&prod2_subcache_t::p_##n##_##m##_##p)
#include "./mygmload.h"
}

#define DEFFUN2(fun) \
  void prod2_subcache_t::fun(double *__restrict__ a,double  *__restrict__ b,double *__restrict__ c)
#include "./mygmcode.h"
