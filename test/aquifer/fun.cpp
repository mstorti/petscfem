//__INSERT_LICENSE__
//$Id: fun.cpp,v 1.2 2002/02/10 02:32:58 mstorti Exp $

//#define USE_ADVANCED_DL_EFN

#include <math.h>

#ifdef USE_ADVANCED_DL_EFN
#include <src/fem.h>
#include <src/getprop.h>
#include <src/ampli.h>
#endif

double f0, f1, t0, t1, slope;

#ifdef USE_ADVANCED_DL_EFN
extern "C" void init_fun(void *t) {
  TextHashTable *thash = (TextHashTable *) t;
  int ierr;
  TGETOPTDEF_ND(thash,double,t0,0.);
  TGETOPTDEF_ND(thash,double,t1,1.);
  TGETOPTDEF_ND(thash,double,f0,0.);
  TGETOPTDEF_ND(thash,double,f1,1.);
#else
extern "C" void init_fun(void *t) {
  t0=0.;
  t1=1.;
  f0=0.;
  f1=1.;
#endif
  slope = (f1-f0)/(t1-t0);
}

extern "C" double eval_fun(double t) {
  if (t<t0) return f0;
  else if (t>t1) return f1;
  else return f0 + slope *(t-t0);
}
