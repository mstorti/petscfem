//__INSERT_LICENSE__
#include <src/fem.h>
#include <src/utils.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "nwadvdif.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void GlobalScalarEF::init"
void GlobalScalarEF::init(int ndof,int ndim,int nel,double Cp_=1.) {
  Cp=Cp_;
  eye_ndof.resize(2,ndof,ndof).set(0.).eye(1.);
  htmp1.resize(1,nel);
  htmp2.resize(2,nel,nel);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void GlobalScalarEF::enthalpy"
void GlobalScalarEF::enthalpy(FastMat2 &H) {
  H.set(UU).scale(Cp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void GlobalScalarEF::comp_W_Cp_N"
void GlobalScalarEF::comp_W_Cp_N(FastMat2 &W_Cp_N,
				 const FastMat2 &W,const FastMat2 &N,double w) {
  htmp1.set(N).scale(w*Cp);
  htmp2.prod(W,htmp1,1,2); // tmp12 = SHAPE' * SHAPE
  W_Cp_N.prod(htmp2,eye_ndof,1,3,2,4); // tmp13 = SHAPE' * SHAPE * I
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void GlobalScalarEF::comp_P_Cp"
void GlobalScalarEF::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.set(P_supg).scale(Cp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void ScalarPerFieldEF::init"
void ScalarPerFieldEF::init(int ndof_,int ndim,int nel) {
  ndof=ndof_;
  Cp.resize(1,ndof);
  htmp1.resize(1,nel);
  htmp2.resize(2,nel,nel);
}

#undef __FUNC__
#define __FUNC__ "void ScalarPerFieldEF::update"
void ScalarPerFieldEF::update(const double *ejac) {
  Cp.set(ejac);
}

#undef __FUNC__
#define __FUNC__ "void ScalarPerFieldEF::enthalpy"
void ScalarPerFieldEF::enthalpy(FastMat2 &H) {
  H.set(UU).mult(Cp);
}

#undef __FUNC__
#define __FUNC__ "void ScalarPerFieldEF::comp_W_Cp_N"
void ScalarPerFieldEF::comp_W_Cp_N(FastMat2 &W_Cp_N,
				   const FastMat2 &W,const FastMat2 &N,double w) {
  htmp1.set(N).scale(w);
  htmp2.prod(W,htmp1,1,2);
  W_Cp_N.set(0.).d(2,4);
  W_Cp_N.prod(htmp2,Cp,1,3,2).rs();
}

#undef __FUNC__
#define __FUNC__ "void ScalarPerFieldEF::comp_P_Cp"
void ScalarPerFieldEF::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.set(P_supg);
  for (int k=1; k<=ndof; k++) {
    P_Cp.ir(2,k).scale(Cp.get(k));
  }
  P_Cp.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void FullEF::init"
void FullEF::init(int ndof,int ndim,int nel) {
  Cp.resize(2,ndof,ndof);
  htmp1.resize(1,nel);
  htmp2.resize(2,nel,nel);
}

#undef __FUNC__
#define __FUNC__ "void FullEF::enthalpy"
void FullEF::enthalpy(FastMat2 &H) {
  H.prod(Cp,UU,1,-1,-1);
}

#undef __FUNC__
#define __FUNC__ "void FullEF::update"
void FullEF::update(const double *ejac) {
  Cp.set(ejac);
}

#undef __FUNC__
#define __FUNC__ "void FullEF::comp_W_Cp_N"
void FullEF::comp_W_Cp_N(FastMat2 &W_Cp_N,
			 const FastMat2 &W,const FastMat2 &N,double w) {
  htmp1.set(N).scale(w);
  htmp2.prod(W,htmp1,1,2);
  W_Cp_N.prod(htmp2,Cp,1,3,2,4);
}

#undef __FUNC__
#define __FUNC__ "void FullEF::comp_P_Cp"
void FullEF::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}
