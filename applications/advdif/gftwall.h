// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: gftwall.h,v 1.1 2005/01/07 22:24:44 mstorti Exp $
#ifndef PETSCFEM_ADVDIF_LAGMUL_H
#define PETSCFEM_ADVDIF_LAGMUL_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class gas_flow_twall : public AdvdifLagrangeMult {
private:
  int ndim;
public:
  ~gas_flow_twall() {} 
  int nres() { return 1; }
  void lag_mul_dof(int jr,int &node,int &dof) {
    assert(jr==0);
    node = 2; dof=1;
  }
  void lm_initialize() { }
  void init() {
    EGETOPTDEF_ND(elemset,int,ndim,0);
    //o Constant of a particular gas for ideal gas
    //  law (state equation for the gas)
    EGETOPTDEF_ND(elemset,double,Rgas,287.);
  }
  void res(int k,FastMat2 &U,FastMat2 & r,
	   FastMat2 & w,FastMat2 & jac) {
    double rho = U.getel(1,1);
    double p = U.getel(1,ndim+2);
    double T = U.getel(2,2);
    assert(r.n()==2 && r.dim(1)==1 && r.dim(2)==1);
    r.setel(p-rho*Rgas*T,1,1);
    w.set(0.);
    w.setel(1.0,1,ndim+2,1);
    jac.setel(Rgas*T,1,1,1);
    jac.setel(1.0,1,1,ndim+2);
  }
  void close() {}
};

#endif
