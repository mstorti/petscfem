//__INSERT_LICENSE__
//$Id: fun2.cpp,v 1.2 2002/02/10 20:11:48 mstorti Exp $

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
  printf("in init: t0 %f,t1 %f,f0 %f,f1 %f\n",t0,t1,f0,f1);

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
  // printf("in eval_fun: t0 %f,t1 %f,f0 %f,f1 %f\n",d->t0,d->t1,d->f0,d->f1);
  if (t < d->t0) return d->f0;
  else if (t > d->t1) return d->f1;
  else return d->f0 + d->slope *(t - d->t0);
}

CLEAR_FUN {
  MyFunData *d = (MyFunData *) fun_data;
  delete d;
  fun_data=NULL;
}
