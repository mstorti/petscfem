//__INSERT_LICENSE__
//$Id: rain2.cpp,v 1.3 2002/09/08 23:29:16 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/elemset.h>

class rain {
private:
  double rain0,T0;
public:
  void init(TextHashTable *thash) {
    int ierr;
    TGETOPTDEF_ND(thash,double,rain0,0.);
    TGETOPTDEF_ND(thash,double,T0,0.);
  }
  double eval(double t,double val) { 
    return (t<T0? rain0 : 0.);
  }
};

PROPERTY_TEMP_FUN(rain);
