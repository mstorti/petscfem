//__INSERT_LICENSE__
//$Id: fun.cpp,v 1.3 2002/02/10 02:49:57 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/ampli.h>

double f0, f1, t0, t1, slope;

extern "C" void init_fun(TextHashTable *thash) {
  int ierr;
  TGETOPTDEF_ND(thash,double,t0,0.);
  TGETOPTDEF_ND(thash,double,t1,1.);
  TGETOPTDEF_ND(thash,double,f0,0.);
  TGETOPTDEF_ND(thash,double,f1,1.);
  slope = (f1-f0)/(t1-t0);
}

extern "C" double eval_fun(double t) {
  if (t<t0) return f0;
  else if (t>t1) return f1;
  else return f0 + slope *(t-t0);
}
