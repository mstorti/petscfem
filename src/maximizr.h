// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: maximizr.h,v 1.1 2001/08/01 22:43:26 mstorti Exp $
#ifndef MAXIMIZR_H
#define MAXIMIZR_H

template<class T> 
class Maximizer {
private:
  int ini;
  T t_;
  double max_;
public:
  void reset() {ini=0;};
  const T &t;
  const double &max;
  Maximizer() : ini(0), max_(0), max(max_), t(t_) {};
  void scan(T &tt,double val);
};

template<class T> 
void Maximizer<T>::scan(T &tt,double val) {
  if (!ini || val>max_) {
    t_=tt; 
    max_=val;
    ini=1;
  }
};

#endif
