//__INSERT_LICENSE__
//$Id: fm2funm.cpp,v 1.1.4.1 2004/05/21 23:46:28 mstorti Exp $

#include <src/fastmat2.h>
#include "./fm2funm.h"

void FastMat2_fund::init(FastMat2 &A) {
  assert(A.n()==2);
  assert(A.dim(1)==A.dim(2));
  m = A.dim(1);
  V.resize(2,m,m);
  Vi.resize(2,m,m);
  D.resize(2,m,m).set(0.);
// corregido por Beto solo para salir del paso , luego que Mario vea como generalizarlo
// mergearlo con la otra version. En el caso no simetrico "_ns" se requiere que lambda
// y flambda sea de 2 x m (autovalores complejos)
  lambda.resize(2,2,m);
  flambda.resize(1,m);
  tmp.resize(2,m,m);
}

void FastMat2_fund::apply(const FastMat2 &A,FastMat2 &fA) {
  lambda.seig(A,V);
  f(lambda,flambda);
  D.d(1,2).set(flambda).rs();
  tmp.prod(D,V,1,-1,2,-1);
  fA.prod(V,tmp,1,-1,-1,2);
}

void FastMat2_fund::apply_ns(const FastMat2 &A,FastMat2 &fA) {
  double tol=1.e-6;
  lambda.eig(A,&V,NULL,1,0);
  lambda.ir(1,2); 
  double lambda_imag_mod = lambda.sum_square_all();
  if(lambda_imag_mod>tol) lambda.print(" Imag(lambda)= ");
  assert(lambda_imag_mod<tol);
  lambda.ir(1,1);
  f(lambda,flambda);
  lambda.rs();
  D.d(1,2).set(flambda).rs();
  Vi.inv(V);
  tmp.prod(D,Vi,1,-1,-1,2);
  fA.prod(V,tmp,1,-1,-1,2);
}

double FastMat2_funm_ff(double x,void *a) {
  FastMat2_funm *fff = (FastMat2_funm *)a;
  return fff->f(x);
}

void FastMat2_funm::f(const FastMat2 &D,FastMat2 &fD) {
  fD.set(D).fun(FastMat2_funm_ff,this);
}

