//__INSERT_LICENSE__
//$Id: gbody.cpp,v 1.2 2005/01/29 16:51:34 mstorti Exp $

#include <math.h>

#include <src/fem.h>
#include <src/getprop.h>
#include <src/elemset.h>

class gbody_scale {
private:
  double T,G0,omega;
public:
  void init(TextHashTable *thash) {
    T = 4;
    G0 = 3;
    omega = 2.0*M_PI/T;
  }
  double eval(double t,double val) { 
    return G0*sin(omega*t);
  }
};

PROPERTY_TEMP_FUN(gbody_scale);
