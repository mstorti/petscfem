//__INSERT_LICENSE__
// #include <src/fem.h>
// #include <src/fastmat2.h>
// #include <src/fm2stats.h>
// #include <src/fastlib2.h>
// #include <src/fm2prod.h>

#include <vector>
using namespace std;

// #define DEFFUN2(fun) 
//   void prod2_subcache_t::fun(double *__restrict__ a,double  *__restrict__ b,double *__restrict__ c)

typedef void (*gemm_fun_t)(double *a,double *b,double *c);

#define DEFFUN2(fun) \
  static void prod2_subcache_##fun(double *__restrict__ a,double  *__restrict__ b,double *__restrict__ c)
#include "./fmgemmcode.h"

void prod2_subcache_load_funs(vector<gemm_fun_t> &table) {
#define LOADFUN(n,m,p,jat,jbt,tindx,fun) table[tindx] = prod2_subcache_##fun
#include "./fmgemmload.h"
}
