// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: smoke.h,v 1.1 2003/01/15 23:50:33 mstorti Exp $
#ifndef SMOKE_H
#define SMOKE_H

#include "advective.h"

#define GETOPT_PROP(type,name,default) elemset->get_prop(name##_prop,#name) //nd

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The flux function for flow in a channel with arbitrary shape and
    using the Kinematic Wave Model.
*/ 
class smoke_ff : public AdvDifFFWEnth {
private:
  const int ndim;
  Property u_prop;
  const NewAdvDif *new_adv_dif_elemset;
  double phi;
  FastMat2 u;
public:
  smoke_ff(const NewAdvDif *e) : AdvDifFFWEnth(e), ndim(2) {}
  ~smoke_ff();
  
  void start_chunk(int &ret_options) {
    new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);     
    elemset->get_prop(u_prop,"u");
    assert(u_prop.length==ndim);
    u.resize(1,ndim);
    // A.resize(3,ndim,1,1);
  }
  
  void element_hook(ElementIterator &element) {
    u.set(new_adv_dif_elemset->prop_array(element_m,u_prop));
    // A.ir(2,1).ir(3,1).set(u).rs();
  }

  // void set_state(const FastMat2 &U,const FastMat2 &grad_U);

  void set_state(const FastMat2 &U) { phi = U.get(1); }

  void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
    // A_grad_N.prod(A,grad_N,-1,2,3,-1,1);
  }

  void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
    assert(0);			// fixme:= Not implemented yet, used for
				// absorbing boundary conditions
  }

  void compute_flux(COMPUTE_FLUX_ARGS) {
  }

  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & dshapex,double w);

  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);

  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w);

  void enthalpy(FastMat2 &H);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);

  int dim() const { return 1; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class smoke : public NewAdvDif {
public:
  smoke() :  NewAdvDif(new smoke_ff(this)) {};
};

#endif
