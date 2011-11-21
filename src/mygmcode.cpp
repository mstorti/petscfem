//__INSERT_LICENSE__
#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fm2stats.h>
#include <src/fastlib2.h>
#include <src/fm2prod.h>

#define DEFFUN2(fun) \
  void prod2_subcache_t::fun(double *__restrict__ a,double  *__restrict__ b,double *__restrict__ c)
#include "./mygmcode.h"
