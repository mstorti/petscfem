//__INSERT_LICENSE__
// $Id: getsurf2.cpp,v 1.1 2005/01/14 21:01:21 mstorti Exp $

#include <cstdio>
#include <libguile.h>

typedef SCM(*scm_fun)();

extern "C"
void getsurf(SCM fun) {
  double xx = 2.3;
  SCM x = scm_make_real(xx);
  SCM ret = scm_apply(fun,x,SCM_EOL);
  printf("hi in getsurf, ret = %f\n",scm_num2dbl(ret,"ja ja"));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
init_femref(void) {
  scm_c_define_gsubr("getsurf",1,0,0,scm_fun(getsurf));
}
