#include <iostream>

class Transformer {
  virtual void transform(int &b) {};
}

template<class B>
class A {
  Transformer *t;
public:
  int transform(B &b);
  int foo(B b) {
    transform(b);
    cout << b << endl;
    return 0;
  }
};
  
int D::transform(int &b) { b=2*b; return 0;}

int main() {
  D a;
  a.foo(2);
}
