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

template < typename Scalar, typename Fun >
class Sample {
public:
  Fun *fun;
  int nstep;
  Scalar start, step;
  int size() { return nstep; }
  Scalar val(int j) const { return fun->f(start+j*step); }
};

class DoubleFun {
public:
  virtual double f(double x)=0;
};

class Sqr : public DoubleFun {
public:
  double f(double x) { return x*x; }
} sqr;

int main() {
  SumArray<Array<int>,int> a;
  int n=5;
  a.resize(n);
  for (int j=0; j<n; j++)
    a[j] = j*j;
  printf("suma: %d\n",a.sum());

  SumArray < Sample <double, DoubleFun>, double > s;
  s.start = 0.;
  s.step  = 0.001;
  s.nstep = 1001;
  s.fun = &sqr;
  printf("suma de sqr * h (should be 0.3333) = %f\n",s.sum()*s.step);

}
