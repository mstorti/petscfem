//__INSERT_LICENSE__
//$Id: fun3.cpp,v 1.1 2002/02/10 22:33:28 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
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

//#define fun_obj_class_name smramp

INIT_FUN1(smramp) {
  smramp *fun_obj_ptr = new smramp;
  fun_data = fun_obj_ptr;
  fun_obj_ptr->init(thash);
}

EVAL_FUN1(smramp) {
  smramp *fun_obj_ptr = (smramp *) fun_data;
  return fun_obj_ptr->eval_fun(t);
}

CLEAR_FUN1(smramp) {
  smramp *fun_obj_ptr = (smramp *) fun_data;
  delete fun_obj_ptr;
}
