//__INSERT_LICENSE__
//$Id: fun.cpp,v 1.1 2002/02/10 00:18:15 mstorti Exp $

#include <math.h>
#include <dlfcn.h>

extern "C" double eval_fun(double t) {
  double f0, f1, t0, t1, slope;
  f0=1.;
  f1=0.3;
  t0=0.;
  t1=1.;
  slope = (f1-f0)/(t1-t0);
  if (t<t0) return f0;
  else if (t>t1) return f1;
  else return f0 + slope *(t-t0);
}
