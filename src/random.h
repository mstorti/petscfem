// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: random.h,v 1.1 2001/09/21 23:14:48 mstorti Exp $
#ifndef RANDOM_H
#define RANDOM_H

namespace Random {

  double drand() {  
    return ((double)(rand()))/((double)(RAND_MAX));
  }

  int irand(int imin,int imax) {
    return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
  }

  class Generator {
  public:
    double get() {return map(drand());};
    virtual double map(double x) {return x;};
  }
 
  Generator uniform;

}

#endif
