// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: randomg.h,v 1.2 2001/09/22 20:40:43 mstorti Exp $
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

  class LogGen : public Generator {
    double map(double x) {return exp(-20.*x);};
  };

  class IntGen : public Generator {
    int start,end;
    double map(double x) {return double(floor((end-start+1)*x)+start);}
  public:
    IntGen(int s=1,int e=100) : start(s), end(e) {};
  };

  extern Generator uniform;
  extern LogGen log_gen;
  extern IntGen int_gen;

}

#endif
