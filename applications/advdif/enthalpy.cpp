//__INSERT_LICENSE__
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"
#include "../../src/fastmat2.h"

#include "nwadvdif.h"

void GlobalScalarEF::init(int ndof,int ndim,int nel,double Cp_=1.) {
  Cp=Cp_;
  eye_ndof.resize(2,ndof,ndof).set(0.).eye(Cp);
  htmp1.resize(1,nel);
  htmp2.resize(2,nel,nel);
}

void GlobalScalarEF::enthalpy(FastMat2 &H, FastMat2 &U) {
  H.set(U).scale(Cp);
}

void GlobalScalarEF::comp_W_Cp_N(FastMat2 &W_Cp_N,
				  FastMat2 &W,FastMat2 &N,double w) {
  htmp1.set(N).scale(w);
  htmp2.prod(W,htmp1,1,2); // tmp12 = SHAPE' * SHAPE
  W_Cp_N.prod(htmp2,eye_ndof,1,3,2,4); // tmp13 = SHAPE' * SHAPE * I
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ScalarperFieldEF::init(int ndof,int ndim,int nel) {
  Cp.resize(1,ndof);
}

void ScalarperFieldEF::update(const double *ejac) {
  Cp.set(ejac);
}

void ScalarperFieldEF::enthalpy(FastMat2 &H, FastMat2 &U) {
  H.set(U).mult(Cp);
}

void ScalarperFieldEF::comp_W_Cp_N(FastMat2 &W_Cp_N,
				  FastMat2 &W,FastMat2 &N,double w) {
  htmp1.resize(1,ndof);
  htmp2.resize(2,nel,nel);

  htmp1.set(N).scale(w);
  htmp2.prod(W,N,1,2);
  W_Cp_N.set(0.).d(4,2);
  W_Cp_N.prod(htmp2,Cp,1,3,2).rs();
}
