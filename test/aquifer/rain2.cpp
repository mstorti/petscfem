//__INSERT_LICENSE__
//$Id: rain2.cpp,v 1.2 2002/09/08 21:59:01 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/elemset.h>

class rain {
private:
  double rain0;
public:
  void init(TextHashTable *thash) {
    int ierr;
    TGETOPTDEF_ND(thash,double,rain0,0.);
  }
  double eval(double t,double val) { 
    return (t<0.01? rain0 : 0.);
  }
};

PROPERTY_TEMP_FUN(rain);
