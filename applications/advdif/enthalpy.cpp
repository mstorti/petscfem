//__INSERT_LICENSE__
//$Id: enthalpy.cpp,v 1.15 2003/07/03 04:32:11 mstorti Exp $
#include <src/fem.h>
#include <src/utils.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "nwadvdif.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void GlobalScalarEF::init"
void GlobalScalarEF::init(int ndof,int ndim,int nel,double Cp_) {
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void user_def_ef_t::init(int ndof,int ndim,int nel,
                         const NewElemset *elemset_) {
  elemset = elemset_;
  int ierr;
  EGETOPTDEF_ND(elemset,double,Cp1,NAN);
  EGETOPTDEF_ND(elemset,double,Cp2,NAN);
  EGETOPTDEF_ND(elemset,double,L,NAN);
  EGETOPTDEF_ND(elemset,double,Tf,NAN);
  EGETOPTDEF_ND(elemset,double,delta,NAN);
  EGETOPTDEF_ND(elemset,double,Cp3,NAN);

  static int flag=0;
  if (!flag) {
    flag = 1;
    SHV(Cp1);
    SHV(Cp2);
    SHV(L);
    SHV(Tf);
    SHV(delta);
    SHV(Cp3);
  }

  eye_ndof.resize(2,ndof,ndof).set(0.).eye(1.);
  htmp1.resize(1,nel);
  htmp2.resize(2,nel,nel);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void user_def_ef_t::update(const double *ejac) {
  mat_indx = int(*ejac);
}

// //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// void user_def_ef_t::set_state(const FastMat2 &U) {
//   UU.set(U);
// }
double user_def_ef_t::hfun(double T) {
  // printf("mat_indx %d\n",mat_indx);
  if (mat_indx==0) {
    double
      hreg = T<Tf? Cp1*T : Cp1*Tf+Cp2*(T-Tf),
      hlat = L*pf_regheavis(T-Tf,-0.5*delta,0.5*delta);
    return hreg + hlat;
  } else {
    return Cp3*T;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
double user_def_ef_t::Cpfun(double T) {
  double epsln=1e-5;
  return (hfun(T+epsln)-hfun(T-epsln))/(2.0*epsln);
 // return (hfun(T+epsln)-hfun(T-epsln))/(2.0*epsln);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void user_def_ef_t::enthalpy(FastMat2 &H) {
  double T=*UU.data();
  H.set(hfun(T));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void user_def_ef_t::comp_W_Cp_N(FastMat2 &W_Cp_N,
                                const FastMat2 &W,
                                const FastMat2 &N,
                                double w) {
  double T=*UU.data();
  double Cp=Cpfun(T);
  htmp1.set(N).scale(w*Cp);
  htmp2.prod(W,htmp1,1,2); // tmp12 = SHAPE' * SHAPE
  W_Cp_N.prod(htmp2,eye_ndof,1,3,2,4); // tmp13 = SHAPE' * SHAPE * I
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void user_def_ef_t::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  double T=*UU.data();
  double Cp=Cpfun(T);
  P_Cp.set(P_supg).scale(Cp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void user_def_ef_t::get_Cp(FastMat2 &Cp_a) {
  double T=*UU.data();
  double Cp=Cpfun(T);
  Cp_a.eye(Cp);
}
