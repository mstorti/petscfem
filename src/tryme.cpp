#include <cstdio>
#include "secant.h"

class myfun : public Secant {
public:
  myfun(double x0) : Secant(x0) {};
  double residual(double x);
};

double myfun::residual(double x) {
  return .7-tanh(x);
}

int main() {
  Secant *s = new myfun(0.9);
  s->omega=0.5;
  s->x0 = 0.9;
  s->tol=1e-12;
  double x = s->sol();
  printf("sol: %12.7f, err: %e, its: %d\n",x,s->f,s->its);
}
