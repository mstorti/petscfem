//__INSERT_LICENSE__
//$Id: fun.cpp,v 1.8 2002/02/10 23:03:52 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/ampli.h>

struct  MyFunData {
  double f0, f1, t0, t1, slope;
};

INIT_FUN {
  int ierr;
  TGETOPTDEF(thash,double,t0,0.);
  TGETOPTDEF(thash,double,t1,1.);
  TGETOPTDEF(thash,double,f0,0.);
  TGETOPTDEF(thash,double,f1,1.);

  MyFunData *d = new MyFunData;
  fun_data = d;
  d->f0 = f0;
  d->f1 = f1;
  d->t0 = t0;
  d->t1 = t1;
  d->slope = (f1-f0)/(t1-t0);
}

EVAL_FUN {
  MyFunData *d = (MyFunData *) fun_data;
  if (t < d->t0) return d->f0;
  else if (t > d->t1) return d->f1;
  else return d->f0 + d->slope *(t - d->t0);
}

CLEAR_FUN {
  MyFunData *d = (MyFunData *) fun_data;
  delete d;
  fun_data=NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
INIT_FUN1(smramp) {
  int ierr;
  TGETOPTDEF(thash,double,t0,0.);
  TGETOPTDEF(thash,double,t1,1.);
  TGETOPTDEF(thash,double,f0,0.);
  TGETOPTDEF(thash,double,f1,1.);

  MyFunData *d = new MyFunData;
  fun_data = d;
  d->f0 = f0;
  d->f1 = f1;
  d->t0 = t0;
  d->t1 = t1;
  d->slope = (f1-f0)/(t1-t0);
}

EVAL_FUN1(smramp) {
  MyFunData *d = (MyFunData *) fun_data;
  if (t < d->t0) return d->f0;
  else if (t > d->t1) return d->f1;
  else {
    double th = (d->t1 + d->t0)/2;
    double xi = 2.*(t - th)/(d->t1 - d->t0);
    double poli =0.5*(3.*xi-cube(xi));
    double dh = (d->f1 + d->f0) / 2;
    return dh + (d->f1 - d->f0)/2. * poli;
  }
}

CLEAR_FUN1(smramp) {
  MyFunData *d = (MyFunData *) fun_data;
  delete d;
  fun_data=NULL;
}
