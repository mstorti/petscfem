//__INSERT_LICENSE__
//$Id: rain.cpp,v 1.2 2002/09/08 21:59:01 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/elemset.h>

PROP_INIT_FUN(rain) {
}

PROP_EVAL_FUN(rain) {
  return (t<0.01? 2. : 0.);
}

PROP_CLEAR_FUN(rain) {
}
