/*__INSERT_LICENSE__*/
// $Id: sp.cpp,v 1.16 2001/09/29 00:57:41 mstorti Exp $

#include <cmath>
#include <vector>

#include <sparse.h>

using namespace Sparse;

class Poly : public ScalarFunObj {
public:
  ~Poly() {};
  double a,b,c;
  double fun(double v) const {return (a*v+b)*v+c;}
} poly;

void my_set(Vec &v,int k=1) {
  v.clear();
  int l = v.length();
  for (int j=0; j<l; j+= k) 
    v.set(j,double(j));
}

double power_nth(double v,void *u) {
  return pow(v,*(double *)u);
}

#define M 1000
int main() {
  int j,k,m,N;

  Vec v(5),w(4),u,res_v;
  GenVec &res = res_v;
  FullVec uu_f;
  GenVec &uu = uu_f;
  Mat a,b,c;
  Indx I(14,17),J(8,14,2),K(0,4,2);
  double d1,d2,tol=1e-11,p,err;

  // Setting individual values
  v.set(2,1.).set(3,2.).print("set values ...  ");
  v.set(10,10.).set(11,11.).print("set at 10 11 ...  ");
  v.grow(0).set(5,5.).print("set at 5 ...  ");
  // v.set(20,20.).print(); // this should give an error
  v.resize(10).print("resized to 10... ");
  v.clear().print("cleared ...  ");
  // Setting values in w
  w.set(0,14.).set(1,15.).set(2,16.).set(3,17.).print("w: ");
  // Inject values of w in v
  v.resize(20).set(I,w).print("setting v[I] = w  ");

  u.resize(20).set(J,v,I).print("u[I] = v[J]");
  u.scale(2.).print("u*2 = ");

  u.print("u:");
  v.print("v:");
  u.axpy(3.,v).print("u = u+3*v ... ");

  printf("u empty OK ? %d\n",!u.empty());
  u.clear();
  printf("u empty after clear() OK ? %d\n",u.empty());

  m=5;
  a.resize(m,m);
  for (j=0; j<m; j++) 
    for (k=0; k<m; k++) 
      if ((j+k) % 2) a.set(j,k,double(j*10+k));
  a.print_f("checkerboard mat");
  
  u.resize(5);
  my_set(u);
  u.print("after my_set; ");
  for (j=0; j<2; j++) 
    a.setr(j,u);
  a.print_f("rows 0 1  set to u: ");

  for (j=3; j<5; j++) 
    a.setc(j,u);
  a.print_f("cols 3 4 set to u: ");

  u.set(a,3).print("u = 3rd row of a\n");
  u.setc(a,3).print("u = 3rd col of a\n");

  b.setr(a,K).print_f("b set to rows 0 2 4 of a");
  printf("size of b: %d\n",b.size());
  b.setc(a,K).print_f("b set to cols 0 2 4 of a");
  printf("size of b: %d\n",b.size());
  
  printf("b not empty OK ? %d\n",!b.empty());
  b.clear();
  printf("b empty OK ? %d\n",b.empty());

  u.clear().resize(100).random_fill(.1,Random::log_gen)
    .print("random fill 100 elements 10\% fill");
  u.purge(1e-2).print("purged to 1e-2");

  a.resize(4,4).id(2.3).print("a set to 2.3*Id(4,4): ");
  a.resize(3,3).id().print("a set to Id(3,3): ");

  u.clear().resize(10).random_fill(.5)
    .print("----\nu set to random_fill(10): ");
  a.diag(u).print("a = diag(u): ");

  u.resize(100).random_fill(.3);
  w.resize(100).random_fill(.2);
  d1 = w.dot(u);
  d2 = u.dot(w);
  printf("u.w: %f, u2.w: %f, err = %g, err<tol ? %d\n",
	 d1,d2,fabs(d1-d2),fabs(d1-d2)<tol);

  a.resize(5,5).random_fill(.2,int_gen).print_f("a with random entries: ");
  u.resize(5).random_fill(.2,int_gen).print_f("u with random entries: ");
  w.resize(5).prod(a,u).print_f("w = a * u :");

  b.resize(5,5).random_fill(.3,int_gen).print_f("b with random entries: ");
  c.resize(5,5).prod(a,b).print_f("c = a * b :");

  poly.a=1.; poly.b=2.; poly.c=3.;
  u.print("u: ");
  u.apply(poly).print("u.^2+2*u+3: ");
  u.apply(sqrt).print("sqrt(u): ");

  p = 1.5;
  u.apply(power_nth,&p).print("u = u^(3/2): ");
  p = 1./1.5;
  u.apply(power_nth,&p).print("u = u^(2/3): ");

  printf("sum(u): %f\n",u.sum());
  printf("max(u): %f\n",u.max());
  printf("min(u): %f\n",u.min()); // should return 0.
  printf("sum_abs(u): %f\n",u.sum_abs());
  printf("sum_sq(u): %f\n",u.sum_sq());
  printf("sum_pow(u,3.): %f\n",u.sum_pow(3.));
  printf("sum_pow(u,5.): %f\n",u.sum_pow(5.));
  printf("max_abs(u): %f\n",u.max_abs());

  u.scale(-1.);
  printf("\n\n---u = -u\n"
	 "max(u): %f\n",u.max()); // should return 0.
  printf("min(u): %f\n",u.min());
  printf("sum_abs(u): %f\n",u.sum_abs());
  printf("sum_sq(u): %f\n",u.sum_sq());
  printf("sum_pow(u,3.): %f\n",u.sum_pow(3.));
  printf("sum_pow(u,5.): %f\n",u.sum_pow(5.));
  printf("max_abs(u): %f\n",u.max_abs());
  
  a.print_f("\n\n----------\na: ");
  a.apply(poly).print_f("a.^2 + 2*a + 3: ");
  a.apply(sqrt).print_f("sqrt(a): ");

  p = 1.5;
  a.apply(power_nth,&p).print_f("a = a^(3/2): ");
  p = 1./1.5;
  a.apply(power_nth,&p).print_f("a = a^(2/3): ");

  printf("sum(a): %f\n",a.sum());
  printf("max(a): %f\n",a.max());
  printf("min(a): %f\n",a.min()); // should return 0.
  printf("sum_abs(a): %f\n",a.sum_abs());
  printf("sum_sq(a): %f\n",a.sum_sq());
  printf("sum_pow(a,3.): %f\n",a.sum_pow(3.));
  printf("sum_pow(a,5.): %f\n",a.sum_pow(5.));
  printf("max_abs(a): %f\n",a.max_abs());

  a.scale(-1.).print_f("\n\n------------\na = -a: ");
  printf("max(a): %f\n",a.max()); // should return 0.
  printf("min(a): %f\n",a.min());
  printf("sum_abs(a): %f\n",a.sum_abs());
  printf("sum_sq(a): %f\n",a.sum_sq());
  printf("sum_pow(a,3.): %f\n",a.sum_pow(3.));
  printf("sum_pow(a,5.): %f\n",a.sum_pow(5.));
  printf("max_abs(a): %f\n",a.max_abs());

  N=5;
  b.clear().resize(N,N).id(1.3);
  a.clear().random_fill(.6).scale(.1).print("0.1*rand(5) (fill .6):");
  a.axpy(1.,b).print_f("Id + 0.1 * rand: ");

  res_v.resize(N);

  u.resize(N).random_fill(.5).print("b:");
  uu.set(u);
  a.solve(uu_f);
  uu.print("sol:");
  res_v.prod(a,uu);
  res_v.axpy(-1.,u);
  err = res_v.sum_abs();
  printf("err < tol ? %d, error: %g \n",err<tol,err);

  u.clear().random_fill(.5).print("b:");
  uu.set(u);
  a.solve(uu_f);
  uu.print("sol:");
  res_v.prod(a,uu);
  res_v.axpy(-1.,u);
  err = res_v.sum_abs();
  printf("err < tol ? %d, error: %g \n",err<tol,err);

#if 1
  N = 500;
  while (1) {
    printf("New matrix...\n");
    b.clear().resize(N,N).id(1.3);
    a.clear().resize(N,N).random_fill(.1).scale(.1).axpy(1.,b);
    b.clear();
    for (int l=0; l<10; l++) {
      u.clear().resize(N).random_fill(.5);
      uu.set(u);
      a.solve(uu_f);
      res_v.resize(N).clear().prod(a,uu);
      res_v.axpy(-1.,u);
      err = res_v.sum_abs();
      printf("size %d\n",N);
      printf("err < tol ? %d, error: %g \n",err<tol,err);
    }
  }
#endif
}
