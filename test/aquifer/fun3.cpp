//__INSERT_LICENSE__
//$Id: fun3.cpp,v 1.2 2002/02/10 23:03:52 mstorti Exp $

#include <src/ampli.h>

class  smramp {
public:
  double f0, f1, t0, t1, slope;
  void init(TextHashTable *thash);
  double eval_fun(double);
};

void smramp::init(TextHashTable *thash) {
  int ierr;
  TGETOPTDEF_ND(thash,double,t0,0.);
  TGETOPTDEF_ND(thash,double,t1,1.);
  TGETOPTDEF_ND(thash,double,f0,0.);
  TGETOPTDEF_ND(thash,double,f1,1.);
  slope = (f1-f0)/(t1-t0);
}

double  smramp::eval_fun(double t) {
  if (t < t0) return f0;
  else if (t > t1) return f1;
  else return f0 + slope *(t - t0);
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(smramp);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class tanh_ramp {
public:
  double A, t0, delta, base;
  void init(TextHashTable *thash);
  double eval_fun(double);
};

void tanh_ramp::init(TextHashTable *thash) {
  int ierr;
  TGETOPTDEF_ND(thash,double,base,0.);
  TGETOPTDEF_ND(thash,double,A,1.);
  TGETOPTDEF_ND(thash,double,t0,0.);
  TGETOPTDEF_ND(thash,double,delta,1.);
  assert(delta>0.);
}

double  tanh_ramp::eval_fun(double t) {
  if (t < t0) return base;
  else return base + A * tanh((t-t0)/delta);
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(tanh_ramp);
