// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: smoke.h,v 1.3 2003/05/26 03:08:06 mstorti Exp $
#ifndef SMOKE_H
#define SMOKE_H

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include "advective.h"
#include "stream.h"

#define GETOPT_PROP(type,name,default) elemset->get_prop(name##_prop,#name) //nd

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The flux function for flow in a channel with arbitrary shape and
    using the Kinematic Wave Model.
*/ 
class smoke_ff : public AdvDifFFWEnth {
private:
  const int ndim;
  Property u_prop;
  double phi;
  FastMat2 u, U, Cp, W_N, A, Uintri;
  int nel,ndof,nelprops;
  ElementIterator element_m;
public:
  smoke_ff(const NewAdvDif *e) : AdvDifFFWEnth(e), ndim(2) {}
  ~smoke_ff();
  
  void start_chunk(int &ret_options);
  
  void element_hook(ElementIterator &element);

  void set_state(const FastMat2 &U);

  void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N);

  void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal);

  void compute_flux(COMPUTE_FLUX_ARGS);

  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & dshapex,double w);

  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);

  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w);

  void enthalpy(FastMat2 &H);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);

  // int dim() const { return 1; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class smoke : public NewAdvDif {
public:
  smoke() :  NewAdvDif(new smoke_ff(this)) {};
};

#endif
