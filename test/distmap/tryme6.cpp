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

template < typename Scalar >
class Sample {
public:
  int nstep;
  Scalar start,step;
  int size() { return nstep; }
  Scalar val(int j) const { return fun(start+j*step); }
};

class DoubleFun : public SumArray< Sample<double>, double > {
public:
  virtual double fun(double x)=0;
};

class Sqr : public DoubleFun {
public:
  double fun(double x) { return x*x; }
} sqr;

int main() {
  SumArray<Array<int>,int> a;
  int n=5;
  a.resize(n);
  for (int j=0; j<n; j++)
    a[j] = j*j;
  printf("suma: %d\n",a.sum());

  sqr.start = 0.;
  sqr.step  = 0.01;
  sqr.nstep = 100;
  printf("suma de sqr= %f\n",sqr.sum());

}
