//__INSERT_LICENSE__
//$Id: fun.cpp,v 1.4 2002/02/10 03:17:28 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/ampli.h>

struct  MyFunData {
  double f0, f1, t0, t1, slope;
};

extern "C" void init_fun(TextHashTable *thash) {
  int ierr;
  TGETOPTDEF(thash,double,t0,0.);
  TGETOPTDEF(thash,double,t1,1.);
  TGETOPTDEF(thash,double,f0,0.);
  TGETOPTDEF(thash,double,f1,1.);

  d = new MyFunData;
  d->f0 = f0;
  d->f1 = f1;
  d->t0 = t0;
  d->t1 = t1;
  d->slope = (f1-f0)/(t1-t0);
}

extern "C" double eval_fun(double t) {
  if (t < d->t0) return d->f0;
  else if (t > d->t1) return d->f1;
  else return d->f0 + d->slope *(t - d->t0);
}

extern "C" void (clear) {
  delete d;
}
