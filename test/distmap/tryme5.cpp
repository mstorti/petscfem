#include <cstdio>
#include <vector>

template < typename Scalar >
class Array : public vector<Scalar> {
public:
  Scalar val(int j) const { return (*this)[j]; }
};

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

int main() {
  SumArray<Array<int>,int> a;
  int n=5;
  a.resize(n);
  for (int j=0; j<n; j++)
    a[j] = j*j;
  printf("suma: %d\n",a.sum());
}
