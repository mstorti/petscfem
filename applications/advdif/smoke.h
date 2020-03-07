// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: smoke.h,v 1.7 2005/02/23 01:40:34 mstorti Exp $
#ifndef PETSCFEM_SMOKE_H
#define PETSCFEM_SMOKE_H

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/fm2temp.h>
#include "advective.h"
#include "stream.h"

#define GETOPT_PROP(type,name,default) elemset->get_prop(name##_prop,#name) //nd

class AJacsm {
public:
  virtual ~AJacsm() { }
public:
  virtual void comp_A_grad_N(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_A_grad_U(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal)=0;
  virtual void comp_Uintri(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_flux(FastMat2 & A,FastMat2 & B) =0 ;
  virtual void comp_vel_per_field(FastMat2 &vel_per_field)=0;
  virtual void comp_vel_vec_per_field(FastMat2 &vel_vec_per_field) {
    assert("Not defined comp_vel_vec_per_field() fun.");
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The flux function for flow in a channel with arbitrary shape and
    using the Kinematic Wave Model.
*/ 
class smoke_ff : public AdvDifFFWEnth {
private:
  int ndim, use_nodal_vel, nodal_vel_indx, fs_advdif_supg;
  Property u_prop, G_prop;
  double phi, omega, drdphi, Cr, phieq, diff_max, diff, 
    tau_fac, diffusivity0, Dt;
  FastMat2 u, s, U, Cp, W_N, A, Uintri, tmp2, tmp0, tmp3;
  int nel,ndof,nelprops;
  ElementIterator element;
  const double *advjac;
  Property advective_jacobians_prop;
  FastMat2Tmp tmp;

  //  AJacsm *a_jac;

public:
  smoke_ff(const NewElemset *e) : AdvDifFFWEnth(e) {}
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class smoke_bcconv : public NewBcconv {
public:
  smoke_bcconv() : NewBcconv(new smoke_ff(this)) {};
  // smoke_bcconv() : NewBcconv(NULL) {};
};

#endif
