//__INSERT_LICENSE__
//$Id: rain.cpp,v 1.1 2002/09/08 16:28:11 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/elemset.h>

PROP_INIT_FUN(rain) {
}

PROP_EVAL_FUN(rain) {
  return 2.;
}

PROP_CLEAR_FUN(rain) {
}
