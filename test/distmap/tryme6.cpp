//__INSERT_LICENSE__
//$Id: tryme6.cpp,v 1.3 2003/01/08 15:49:04 mstorti Exp $
#include <cstdio>
#include <vector>

template < typename Array, typename Scalar >
class SumArray : public Array {
public:
  Scalar sum() {
    int k,n = size();
    Scalar sum_c = 0;
    for (k=0; k<n; k++) 
      sum_c += val(k);
    return sum_c;
  }
};

template < typename Scalar >
class Array : public vector<Scalar> {
public:
  Scalar val(int j) const { return (*this)[j]; }
};

template < typename Scalar>
class Sample {
public:
  class Fun {
  public:
    virtual Scalar f(Scalar x)=0;
  } *fun;
  int nstep;
  Scalar start, step;
  int size() { return nstep; }
  Scalar val(int j) const { return fun->f(start+j*step); }
};

class Sqr : public Sample<double>::Fun {
public:
  double f(double x) { return x*x; }
} sqr;

class Cub : public Sample<double>::Fun {
public:
  double f(double x) { return x*x*x; }
} cub;

int main() {
  SumArray<Array<int>,int> a;
  int n=5;
  a.resize(n);
  for (int j=0; j<n; j++)
    a[j] = j*j;
  printf("suma: %d\n",a.sum());

  SumArray < Sample <double>, double > s;
  s.start = 0.;
  s.step  = 0.001;
  s.nstep = 1001;
  s.fun = &sqr;
  printf("suma de sqr * h (should be 0.3333) = %f\n",s.sum()*s.step);
  s.fun = &cub;
  printf("suma de cub * h (should be 0.25)   = %f\n",s.sum()*s.step);
}
