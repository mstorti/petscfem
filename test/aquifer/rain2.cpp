//__INSERT_LICENSE__
//$Id: rain2.cpp,v 1.1 2002/09/08 20:37:01 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/elemset.h>

class rain {
public:
  void init(TextHashTable *thash) {}
  double eval(double t,double val) { 
    return (t<0.01? 2. : 0.);
  }
};

PROP_INIT_FUN(rain) {
  rain *fun_obj = new rain;
  fun_data = fun_obj;
  fun_obj->init(thash);
}

PROP_EVAL_FUN(rain) {
  rain *fun_obj = (rain *)fun_data;
  return fun_obj->eval(t,val);
}

PROP_CLEAR_FUN(rain) {
  delete (rain *)fun_data;
}
