// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: randomg.h,v 1.1 2001/09/22 04:29:36 mstorti Exp $
#ifndef RANDOMG_H
#define RANDOMG_H

#include <cstdlib>
#include <cmath>

namespace Random {

  inline double drand() {  
    return ((double)(rand()))/((double)(RAND_MAX));
  }

  inline int irand(int imin,int imax) {
    return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
  }

  class Generator {
  public:
    double get() {return map(drand());};
    virtual double map(double x) {return x;};
  };

  extern Generator uniform;

}

#endif
