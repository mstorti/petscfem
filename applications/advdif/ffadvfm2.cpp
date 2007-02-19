//__INSERT_LICENSE__
//$Id: ffadvfm2.cpp,v 1.54.8.1 2007/02/19 20:23:56 mstorti Exp $

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>
//#include <string>

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "nwadvdif.h"

NewAdvDifFF::~NewAdvDifFF() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"newadvecfm2_ff_t::newadvecfm2_ff_t(NewElemset *elemset_)"
newadvecfm2_ff_t::newadvecfm2_ff_t(NewElemset *elemset_) 
  : NewAdvDifFF(elemset_), 
    null_source_term(*this), gscalar_source_term(*this),
    full_source_term(*this),
    full_c_jac(*this), scalar_c_jac(*this),
    null_c_jac(*this),
    scalar_per_field_c_jac(*this),
    null_a_jac(*this),
    u_global(*this), u_per_field(*this),
    full_adv_jac(*this), full_dif_jac(*this),
    global_dif_tensor(*this), per_field_dif_tensor(*this),
    null_d_jac(*this), global_scalar_djac(*this),
    scalar_dif_per_field(*this)
{}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__				\
"void newadvecfm2_ff_t::NullSourceTerm::	\
add_source_term(FastMat2 &G_source)"
void newadvecfm2_ff_t::NullSourceTerm::
add_source_term(FastMat2 &G_source) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__				\
"void newadvecfm2_ff_t::GScalarSourceTerm::	\
add_source_term(FastMat2 &G_source)"
void newadvecfm2_ff_t::GScalarSourceTerm::
add_source_term(FastMat2 &G_source) {
  G_source.add(*ff.s_body);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullSourceTerm::\
add_source_term(FastMat2 &G_source)"
void newadvecfm2_ff_t::FullSourceTerm::
add_source_term(FastMat2 &G_source) {
  G_source.add(ff.S_body);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__					\
"void newadvecfm2_ff_t::NullCJac"			\
"::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w)"
void newadvecfm2_ff_t::NullCJac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  N_N_C.set(0.);
}

#undef __FUNC__
#define __FUNC__					\
"void newadvecfm2_ff_t::NullCJac::"			\
"comp_G_source(FastMat2 &G_source, FastMat2 &U)"
void newadvecfm2_ff_t::NullCJac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.set(0.);
}


#undef __FUNC__
#define __FUNC__				\
"void newadvecfm2_ff_t::NullCJac::"		\
"comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,"	\
 "FastMat2 &N,double w)"
void newadvecfm2_ff_t::NullCJac::
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	   FastMat2 &N,double w) {
  N_P_C.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullCJac\
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w)"
void newadvecfm2_ff_t::FullCJac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.set(N).scale(w);
  tmp.prod(tmp2,N,1,2);
  N_N_C.prod(tmp,ff.C_jac,1,3,2,4);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullCJac::\
comp_G_source(FastMat2 &G_source, FastMat2 &U)"
void newadvecfm2_ff_t::FullCJac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.prod(ff.C_jac,U,1,-1,-1);
}


#undef __FUNC__
#define __FUNC__ \
"newadvecfm2_ff_t::FullCJac::comp_N_P_C()"
void newadvecfm2_ff_t::FullCJac::
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	   FastMat2 &N,double w) {
  tmp26.set(P_supg).scale(w);
  tmp27.prod(N,ff.C_jac,1,2,3); // tmp27 = N * C_jac 
  N_P_C.prod(tmp26,tmp27,1,-1,2,-1,3); // tmp28 = P_supg * C_jac * N
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarCJac\
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w)"
void newadvecfm2_ff_t::ScalarCJac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.set(N).scale(w*(*ff.reacjac));
  tmp.prod(tmp2,N,1,2);
  N_N_C.prod(tmp,ff.eye_ndof,1,3,2,4);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarCJac::\
comp_G_source(FastMat2 &G_source, FastMat2 &U)"
void newadvecfm2_ff_t::ScalarCJac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.set(U).scale(*ff.reacjac);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarCJac::\
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,\
	   FastMat2 &N,double w)"
void newadvecfm2_ff_t::ScalarCJac::
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	   FastMat2 &N,double w) {
  tmp26.set(P_supg).scale(w*(*ff.reacjac));
  // perhaps this may be optimized by computing the following
  // operations in reverse order
  tmp27.prod(N,ff.eye_ndof,1,2,3);        // tmp27 = N * C_jac 
  N_P_C.prod(tmp26,tmp27,1,-1,2,-1,3); // tmp28 = P_supg * C_jac * N
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarPerFieldCjac\
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w)"
void newadvecfm2_ff_t::ScalarPerFieldCjac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.set(N).scale(w);
  tmp.prod(tmp2,N,1,2);
  N_N_C.set(0.).d(4,2);
  N_N_C.prod(tmp,ff.C_jac,1,2,3);
  N_N_C.rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarPerFieldCjac::\
comp_G_source(FastMat2 &G_source, FastMat2 &U)"
void newadvecfm2_ff_t::ScalarPerFieldCjac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.set(U).mult(ff.C_jac);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarPerFieldCjac::\
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,\
	   FastMat2 &N,double w)"
void newadvecfm2_ff_t::ScalarPerFieldCjac::
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	   FastMat2 &N,double w) {
  tmp26.set(P_supg).scale(w);
  ff.N_C.set(0.).d(3,2).prod(N,ff.C_jac,1,2).rs();
  N_P_C.prod(tmp26,ff.N_C,1,-1,2,-1,3); // tmp28 = P_supg * C_jac * N
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullDjac\
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U)"
void newadvecfm2_ff_t::NullDjac
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.set(0.);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullDjac\
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,\
		       FastMat2 & dshapex,double w)"
void newadvecfm2_ff_t::NullDjac
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.set(0.);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullDjac\
::comp_dif_per_field(FastMat2 &dif_per_field)"
void newadvecfm2_ff_t::NullDjac
::comp_dif_per_field(FastMat2 &dif_per_field) {
  dif_per_field.set(0.);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::GlobalScalar\
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U)"
void newadvecfm2_ff_t::GlobalScalar
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.t().set(grad_U).scale(*(ff.difjac)).rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::GlobalScalar\
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,\
		       FastMat2 & dshapex,double w)"
void newadvecfm2_ff_t::GlobalScalar
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  tmp.prod(dshapex,dshapex,-1,1,-1,2).scale(w*(*(ff.difjac)));
  grad_N_D_grad_N.prod(tmp,ff.eye_ndof,1,3,2,4);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::GlobalScalar\
::comp_dif_per_field(FastMat2 &dif_per_field)"
void newadvecfm2_ff_t::GlobalScalar
::comp_dif_per_field(FastMat2 &dif_per_field) {
  dif_per_field.set(*(ff.difjac));
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::GlobalDifTensor\
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U)"
void newadvecfm2_ff_t::GlobalDifTensor
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.prod(ff.D_jac,grad_U,2,-1,-1,1);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::GlobalDifTensor\
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,\
		       FastMat2 & dshapex,double w)"
void newadvecfm2_ff_t::GlobalDifTensor
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  tmp1.prod(ff.D_jac,dshapex,1,-1,-1,2).scale(w);
  tmp2.prod(dshapex,tmp1,-1,1,-1,2);
  grad_N_D_grad_N.prod(tmp2,ff.eye_ndof,1,3,2,4);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::GlobalDifTensor\
::comp_dif_per_field(FastMat2 &dif_per_field)"
void newadvecfm2_ff_t::GlobalDifTensor
::comp_dif_per_field(FastMat2 &dif_per_field) {
  double dd = ff.D_jac.d(2,1).sum_all()/double(ff.ndim);
  ff.D_jac.rs();
  dif_per_field.set(dd);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::PerFieldDifTensor\
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U)"
void newadvecfm2_ff_t::PerFieldDifTensor
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  for (int k=1; k<=ff.ndof; k++) {
    ff.D_jac.ir(1,k);
    grad_U.ir(2,k);
    fluxd.ir(1,k).prod(ff.D_jac,grad_U,1,-1,-1);
  }
  fluxd.rs();
  grad_U.rs();
  ff.D_jac.rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::PerFieldDifTensor\
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,\
		       FastMat2 & dshapex,double w)"
void newadvecfm2_ff_t::PerFieldDifTensor
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.set(0.);
  for (int k=1; k<=ff.ndof; k++) {
    ff.D_jac.ir(1,k);
    tmp.prod(ff.D_jac,dshapex,1,-1,-1,2).scale(w);
    grad_N_D_grad_N.ir(2,k).ir(4,k).prod(dshapex,tmp,-1,1,-1,2);
  }
  ff.D_jac.rs();
  grad_N_D_grad_N.rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::PerFieldDifTensor\
::comp_dif_per_field(FastMat2 &dif_per_field)"
void newadvecfm2_ff_t::PerFieldDifTensor
::comp_dif_per_field(FastMat2 &dif_per_field) {
  ff.D_jac.d(3,2);
  dif_per_field.sum(ff.D_jac,1,-1).scale(1./double(ff.ndim));
  ff.D_jac.rs();
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarDifPerField\
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U)"
void newadvecfm2_ff_t::ScalarDifPerField
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.t().set(grad_U).rs();
  for (int j=1; j<=ff.ndof; j++) 
    fluxd.ir(1,j).scale(ff.D_jac.get(j));
  fluxd.rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarDifPerField\
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,\
		       FastMat2 & dshapex,double w)"
void newadvecfm2_ff_t::ScalarDifPerField
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.set(0.);
  tmp.set(ff.D_jac).scale(w);
  grad_N_grad_N.prod(dshapex,dshapex,-1,1,-1,2);
  grad_N_D_grad_N.d(2,4).prod(grad_N_grad_N,tmp,1,3,2).rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::ScalarDifPerField\
::comp_dif_per_field(FastMat2 &dif_per_field)"
void newadvecfm2_ff_t::ScalarDifPerField
::comp_dif_per_field(FastMat2 &dif_per_field) {
  //  dif_per_field.set(*(ff.difjac));
  dif_per_field.set((ff.difjac));
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullDifJac\
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U)"
void newadvecfm2_ff_t::FullDifJac
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.prod(ff.D_jac,grad_U,1,-2,2,-1,-1,-2);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullDifJac\
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,\
		       FastMat2 & dshapex,double w)"
void newadvecfm2_ff_t::FullDifJac
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  D_grad_N.prod(ff.D_jac,dshapex,2,3,1,-1,-1,4)
    .scale(w);
  grad_N_D_grad_N.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullDifJac\
::comp_dif_per_field(FastMat2 &dif_per_field)"
void newadvecfm2_ff_t::FullDifJac
::comp_dif_per_field(FastMat2 &dif_per_field) {
  ff.D_jac.d(2,1).d(4,3);
  dif_per_field.sum_abs(ff.D_jac,1,-1).scale(1./ff.ndim);
  ff.D_jac.rs();
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullAdvJac::\
comp_flux(FastMat2 &flux,FastMat2 &U)"
void newadvecfm2_ff_t::FullAdvJac::
comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.prod(ff.u,U,2,1,-1,-1);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullAdvJac::\
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal)"
void newadvecfm2_ff_t::FullAdvJac::
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.prod(ff.u,normal,-1,1,2,-1);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullAdvJac::\
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U)"
void newadvecfm2_ff_t::FullAdvJac::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.prod(ff.u,grad_U,-1,1,-2,-1,-2);
}
  
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullAdvJac::\
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex)"
void newadvecfm2_ff_t::FullAdvJac::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  A_grad_N.prod(dshapex,ff.u,-1,1,-1,2,3);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullAdvJac::\
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco)"
void newadvecfm2_ff_t::FullAdvJac::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  // Here we take the diagonal of each Jacobian as the component
  // of a per field vector velocity
  ff.u.d(3,2);
  Uintri.prod(iJaco,ff.u,2,-1,-1,1);
  ff.u.rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::FullAdvJac::\
comp_vel_per_field(FastMat2 &vel_per_field)"
void newadvecfm2_ff_t::FullAdvJac::
comp_vel_per_field(FastMat2 &vel_per_field) {
  // This is approximate. We tak as velocity for a field
  // the sum over spatial dimensions of the digonal
  // elements of the jacobians:
  ff.u.d(3,2);
  vel_per_field.sum_square(ff.u,-1,1).fun(sqrt);
  ff.u.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullAJac::comp_flux(FastMat2 &flux,FastMat2 &U)"
void newadvecfm2_ff_t::NullAJac::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.set(0.);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullAJac::\
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal)"
void newadvecfm2_ff_t::NullAJac::
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.set(0.);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullAJac::\
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U)"
void newadvecfm2_ff_t::NullAJac::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.set(0.);
}
  
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullAJac::\
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex)"
void newadvecfm2_ff_t::NullAJac::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  A_grad_N.set(0.);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullAJac::\
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco)"
void newadvecfm2_ff_t::NullAJac::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  Uintri.set(0.);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::NullAJac::\
comp_vel_per_field(FastMat2 &vel_per_field)"
void newadvecfm2_ff_t::NullAJac::
comp_vel_per_field(FastMat2 &vel_per_field) {
  vel_per_field.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UPerField::\
comp_flux(FastMat2 &flux,FastMat2 &U)"
void newadvecfm2_ff_t::UPerField::
comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.set(ff.u);
  for (int j=1; j<=ff.ndof; j++) {
    flux.ir(1,j).scale(U.get(j));
  }
  flux.rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UPerField::\
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal)"
void newadvecfm2_ff_t::UPerField::
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.set(0.).d(2,1).prod(ff.u,normal,1,-1,-1);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UPerField::\
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U)"
void newadvecfm2_ff_t::UPerField::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  for (int j=1; j<=ff.ndof; j++) {
    grad_U.ir(2,j);
    ff.u.ir(1,j);
    tmp.prod(grad_U,ff.u,-1,-1);
    A_grad_U.setel(tmp.get(),j);
  }
  grad_U.rs();
  ff.u.rs();
  A_grad_U.rs();
}
  
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UPerField::\
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex)"
void newadvecfm2_ff_t::UPerField::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  A_grad_N.set(0.).d(3,2).prod(dshapex,ff.u,-1,1,2,-1).rs();
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UPerField::\
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco)"
void newadvecfm2_ff_t::UPerField::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  Uintri.prod(iJaco,ff.u,2,-1,1,-1);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UPerField::\
comp_vel_per_field(FastMat2 &vel_per_field)"
void newadvecfm2_ff_t::UPerField::
comp_vel_per_field(FastMat2 &vel_per_field) {
  vel_per_field.sum_square(ff.u,1,-1).fun(sqrt);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UGlobal::\
comp_flux(FastMat2 &flux,FastMat2 &U)"
void newadvecfm2_ff_t::UGlobal::
comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.prod(U,ff.u,1,2);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UGlobal::\
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal)"
void newadvecfm2_ff_t::UGlobal::
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  tmp5.prod(ff.u,normal,-1,-1);
  double un = tmp5.get();
  A_jac_n.eye(un);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UGlobal::\
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U)"
void newadvecfm2_ff_t::UGlobal::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.prod(ff.u,grad_U,-1,-1,1);
}
  
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UGlobal::\
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex)"
void newadvecfm2_ff_t::UGlobal::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  tmp.prod(ff.u,dshapex,-1,-1,1);
  for (int j=1; j<=ff.nel; j++) {
    A_grad_N.ir(1,j).eye(tmp.get(j));
  }
  A_grad_N.rs();
}
  
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UGlobal::\
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco)"
void newadvecfm2_ff_t::UGlobal::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  tmp3.prod(iJaco,ff.u,1,-1,-1);
  Uintri.prod(ff.tmp2,tmp3,1,2);
}

#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::UGlobal::\
comp_vel_per_field(FastMat2 &vel_per_field)"
void newadvecfm2_ff_t::UGlobal::
comp_vel_per_field(FastMat2 &vel_per_field) {
  double vv= ff.u.sum_square_all();
  vel_per_field.set(sqrt(vv));
}

#undef __FUNC__
#define __FUNC__ "comp_vel_vec_per_field()"
void newadvecfm2_ff_t::UGlobal::
comp_vel_vec_per_field(FastMat2 &vel_vec_per_field) {
  for (int k=1; k<=ff.ndof; k++) {
    vel_vec_per_field.ir(1,k).set(ff.u);
  }
  vel_vec_per_field.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void advecfm2_ff_t::element_hook(ElementIterator &)"
// This is to pass to the advective function the element
#undef __FUNC__
#define __FUNC__ \
"void newadvecfm2_ff_t::element_hook(ElementIterator &element_)"
void newadvecfm2_ff_t::element_hook(ElementIterator &element_) {
  element = element_;

  advjac = elemset->prop_array(element,advective_jacobians_prop);
  u.set(advjac);

  difjac_mol = elemset->prop_array(element,diffusive_jacobians_mol_prop);
  difjac = elemset->prop_array(element,diffusive_jacobians_prop);
  //  difjac = difjac_mol;
  d_jac->update(difjac);

  reacjac = elemset->prop_array(element,reactive_jacobians_prop);
  C_jac.set(reacjac);

  s_body = elemset->prop_array(element,source_term_prop);
  S_body.set(s_body);

  e_jac = elemset->prop_array(element,enthalpy_jacobians_prop);
  enthalpy_fun->update(e_jac);
}  

//  enum advective_jacobian_type {
//    undefined, global_vector, vector_per_field, full};

//  enum diffusive_jacobians_type {
//    undefined, global_scalar, scalar_per_field, full};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void advecfm2_ff_t::start_chunk(int ret_options)"
void newadvecfm2_ff_t::start_chunk(int &ret_options) {
  int ierr;
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset); 
  EGETOPTDEF_ND(elemset,int,ndim,0); //nd

  //o Add LES for this particular elemset.
  EGETOPTDEF_ND(elemset,int,LES,0);
  //o van Driest constant for the damping law.
  EGETOPTDEF_ND(elemset,double,A_van_Driest,0);
  assert(A_van_Driest==0.);
  //o Smagorinsky constant.
  EGETOPTDEF_ND(elemset,double,C_smag,0.18);
  //o turbulent Prandtl
  EGETOPTDEF_ND(elemset,double,Pr_t,1.0);


  //o Scale the SUPG upwind term. Set to 0 in order to
  //  not to include the upwind term. 
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);

  elemset->elem_params(nel,ndof,nelprops);
  Uintri.resize(2,ndof,ndim);
  tmp2.resize(1,ndof).set(1.);
  tmp3.set(tmp2);
  dif_per_field.resize(1,ndof);
  vel_per_field.resize(1,ndof);
  vel_vec_per_field.resize(2,ndof,ndim);
  eye_ndof.resize(2,ndof,ndof).eye();

  strain_rate.resize(2,ndim,ndim);

  ret_options &= !SCALAR_TAU; // tell the advective element routine

  // Read enthalpy jacobians
  //o _T: double[var_len]
  //  _N: enthalpy_jacobians _D: no default  _DOC: 
  //i_tex ../../doc/advdifop.tex enthalpy_jacobians
  //  _END
  elemset->get_prop(enthalpy_jacobians_prop,"enthalpy_jacobians");

  //o Set enthalpy jacobian to the desired type. May be one of 
  // `` #ident# '', `` #global_scalar# '',
  // `` #scalar_per_field# '' or `` #full# ''. 
  // See documentation for the  #enthalpy_jacobians#  option. 
  EGETOPTDEF(elemset,string,enthalpy_jacobians_type,string("undefined"));
  string enthalpy_jacobians_type_s = enthalpy_jacobians_type;

  if (enthalpy_jacobians_type==string("undefined")) {
    if (enthalpy_jacobians_prop.length == 0) {
      enthalpy_jacobians_type=string("ident");
    } else if (enthalpy_jacobians_prop.length == 1) {
      enthalpy_jacobians_type=string("global_scalar");
    } else if (enthalpy_jacobians_prop.length == ndof) {
      enthalpy_jacobians_type=string("scalar_per_field");
    } else if (enthalpy_jacobians_prop.length == ndof*ndof) {
      enthalpy_jacobians_type=string("full");
    }
  }
  if (enthalpy_jacobians_type==string("ident") &&
      enthalpy_jacobians_prop.length == 0) {
    enthalpy_fun = &identity_ef;
    identity_ef.init(ndof,ndim,nel);
  } else if (enthalpy_jacobians_type==string("global_scalar") &&
      enthalpy_jacobians_prop.length == 1) {
    enthalpy_fun = &global_scalar_ef;
    global_scalar_ef.init(ndof,ndim,nel);
  } else if (enthalpy_jacobians_type==string("scalar_per_field") &&
	     enthalpy_jacobians_prop.length == ndof) {
    enthalpy_fun = &scalar_per_field_ef;
    scalar_per_field_ef.init(ndof,ndim,nel);
  } else if (enthalpy_jacobians_type==string("full") &&
	     enthalpy_jacobians_prop.length == ndof*ndof) {
    enthalpy_fun = &full_ef;
    full_ef.init(ndof,ndim,nel);
  } else 
    PETSCFEM_ERROR("Not implemented yet\n"
		   "Enthalpy Jacobian does not fit in \n"
		   "a predefined name/length combination\n"
		   "name entered: \"%s\"\n"
		   "length entered: %d\n",
		   enthalpy_jacobians_type_s.c_str(),
		   enthalpy_jacobians_prop.length);

  // Read advective jacobians
  //o _T: double[var_len]
  //  _N: advective_jacobians _D: no default  _DOC: 
  //i_tex ../../doc/advdifop.tex advective_jacobians
  //  _END
  elemset->get_prop(advective_jacobians_prop,"advective_jacobians");

  //o Set advective jacobian to the desired type. May be one of 
  // `` #null# '', `` #global_vector# '',
  // `` #vector_per_field# '' or `` #full# ''. 
  // See documentation for the  #advective_jacobians#  option. 
  EGETOPTDEF(elemset,string,advective_jacobians_type,string("undefined"));
  string advective_jacobians_type_s = advective_jacobians_type;

  if (advective_jacobians_type==string("undefined")) {
    if (advective_jacobians_prop.length == 0) {
      advective_jacobians_type=string("null");
    } else if (advective_jacobians_prop.length == ndim) {
      advective_jacobians_type=string("global_vector");
    } else if (advective_jacobians_prop.length == ndim*ndof) {
      advective_jacobians_type=string("vector_per_field");
    } else if (advective_jacobians_prop.length == ndim*ndof*ndof) {
      advective_jacobians_type=string("full");
    } 
  }

  if (advective_jacobians_type==string("null") &&
      advective_jacobians_prop.length == 0) {
    a_jac =  &null_a_jac;
  } else if (advective_jacobians_type==string("global_vector") &&
      advective_jacobians_prop.length == ndim) {
    u.resize(1,ndim);
    a_jac =  &u_global;
  } else if (advective_jacobians_type==string("vector_per_field") &&
	     advective_jacobians_prop.length == ndim*ndof) {
    u.resize(2,ndof,ndim);
    a_jac =  &u_per_field;
  } else if (advective_jacobians_type==string("full") &&
	     advective_jacobians_prop.length == ndim*ndof*ndof) {
    u.resize(3,ndim,ndof,ndof);
    a_jac =  &full_adv_jac;
  } else {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Advective Jacobian does not fit in \n"
		"a predefined name/length combination\n"
		"name entered: \"%s\"\n"
		"length entered: %d\n",
		advective_jacobians_type_s.c_str(),
		advective_jacobians_prop.length);
    abort(); // Not defined type of Jacobian
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Read diffusive jacobians (diffusivity matrices)
  //o _T: double[var_len] 
  //  _N: diffusive_jacobians _D: no default  _DOC: 
  //i_tex ../../doc/advdifop.tex diffusive_jacobians
  //  _END
  //
  // molecular diffusive jacobians needed to add the turbulence correction afterwards
  elemset->get_prop(diffusive_jacobians_mol_prop,"diffusive_jacobians_mol");
  //
  elemset->get_prop(diffusive_jacobians_prop,"diffusive_jacobians");

  assert(diffusive_jacobians_prop.length==diffusive_jacobians_mol_prop.length);

  //o Set diffusive jacobian to the desired type
  //  May be one of 
  // `` #null# '', 
  // `` #global_scalar# '', `` #global_tensor# '', 
  // `` #tensor_per_field# '', `` #scalar_per_field# ''
  //   or `` #null# ''
  // See documentation for the  #diffusive_jacobians#  option. 
  EGETOPTDEF_ND(elemset,string,diffusive_jacobians_type,string("undefined"));
  string diffusive_jacobians_type_s=diffusive_jacobians_type;

  if (diffusive_jacobians_type==string("undefined")) {
    if (diffusive_jacobians_prop.length == 0) {
      diffusive_jacobians_type=string("null");
    } else if (diffusive_jacobians_prop.length == 1) {
      diffusive_jacobians_type=string("global_scalar");
    } else if (diffusive_jacobians_prop.length == ndim*ndim) {
      diffusive_jacobians_type=string("global_tensor");
    } else if (diffusive_jacobians_prop.length == ndim*ndim*ndof) {
      diffusive_jacobians_type=string("tensor_per_field");
    } else if (diffusive_jacobians_prop.length == ndof) {
      diffusive_jacobians_type=string("scalar_per_field");
    } else if (diffusive_jacobians_prop.length == ndim*ndim*ndof*ndof) {
      diffusive_jacobians_type=string("full");
    }
  }
  if (diffusive_jacobians_type==string("null") &&
      diffusive_jacobians_prop.length == 0) {
    d_jac =  &null_d_jac;
  } else if (diffusive_jacobians_type==string("global_scalar") &&
      diffusive_jacobians_prop.length == 1) {
    d_jac =  &global_scalar_djac;
  } else if (diffusive_jacobians_type==string("global_tensor") &&
	     diffusive_jacobians_prop.length == ndim*ndim) {
    D_jac.resize(2,ndim,ndim);
    d_jac =  &global_dif_tensor;
  } else if (diffusive_jacobians_type==string("tensor_per_field") &&
	     diffusive_jacobians_prop.length == ndim*ndim*ndof) {
    D_jac.resize(3,ndof,ndim,ndim);
    d_jac =  &per_field_dif_tensor;
  } else if (diffusive_jacobians_type==string("scalar_per_field") &&
	     diffusive_jacobians_prop.length == ndof) {
    D_jac.resize(1,ndof);
    d_jac =  &scalar_dif_per_field;
  } else if (diffusive_jacobians_type==string("full") &&
	     diffusive_jacobians_prop.length == ndim*ndim*ndof*ndof) {
    D_jac.resize(4,ndof,ndof,ndim,ndim);
    d_jac =  &full_dif_jac;
  } else {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Diffusive Jacobian does not fit in \n"
		"a predefined name/length combination\n"
		"name entered: \"%s\"\n"
		"length entered: %d\n",
		diffusive_jacobians_type_s.c_str(),
		diffusive_jacobians_prop.length);
    assert(0);
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Read reactive jacobians (reactive  matrix)
  //o _T: double[var_len] 
  //  _N: reactive_jacobians _D: all zero  _DOC: 
  //i_tex ../../doc/advdifop.tex reactive_jacobians
  //  _END
  elemset->get_prop(reactive_jacobians_prop,"reactive_jacobians");

  //o Set reactive jacobian to the desired type. 
  //  May be one of 
  // `` #null# '', 
  // `` #global_scalar# '', 
  // `` #scalar_per_field# ''
  // or `` #null# ''
  // See documentation for the  #reactive_jacobians#  option. 
  EGETOPTDEF(elemset,string,reactive_jacobians_type,string("undefined"));
  string reactive_jacobians_type_s=reactive_jacobians_type;

  if (reactive_jacobians_type==string("undefined")) {
    if (reactive_jacobians_prop.length == 0) {
      reactive_jacobians_type=string("null");
    } else if (reactive_jacobians_prop.length == 1) {
      reactive_jacobians_type=string("global_scalar");
    } else if (reactive_jacobians_prop.length == ndof) {
      reactive_jacobians_type=string("scalar_per_field");
    } else if (reactive_jacobians_prop.length == ndof*ndof) {
      reactive_jacobians_type=string("full");
    }
  }
  if (reactive_jacobians_type==string("null") &&
      reactive_jacobians_prop.length == 0) {
    c_jac =  &null_c_jac;
  } else if (reactive_jacobians_type==string("global_scalar") &&
      reactive_jacobians_prop.length == 1) {
    c_jac =  &scalar_c_jac;
  } else if (reactive_jacobians_type==string("scalar_per_field") &&
	     reactive_jacobians_prop.length == ndof) {
    N_C.resize(3,nel,ndof,ndof);
    C_jac.resize(1,ndof);
    c_jac =  &scalar_per_field_c_jac;
  } else if (reactive_jacobians_type==string("full") &&
	     reactive_jacobians_prop.length == ndof*ndof) {
    C_jac.resize(2,ndof,ndof);
    c_jac =  &full_c_jac;
  } else {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Reactive Jacobian does not fit in \n"
		"a predefined name/length combination\n"
		"name entered: \"%s\"\n"
		"length entered: %d\n",
		reactive_jacobians_type_s.c_str(),
		reactive_jacobians_prop.length);
    assert(0);
  }

  // Source term
  //o _T: double[var_len] 
  //  _N: source term  _D: all zero  _DOC: 
  //i_tex ../../doc/advdifop.tex source_term
  //  _END
  elemset->get_prop(source_term_prop,"source_term");

  //o Set source term to the desired type
  //  May be one of 
  // `` #null# '', 
  // `` #global_scalar# '', 
  // or `` #null# ''
  // See documentation for the `` #source_term# '' option. 
  EGETOPTDEF(elemset,string,source_term_type,string("undefined"));
  string source_term_type_s=source_term_type;

  if (source_term_type==string("undefined")) {
    if (source_term_prop.length == 0) {
      source_term_type=string("null");
    } else if (source_term_prop.length == 1) {
      source_term_type=string("global_scalar");
    } else if (source_term_prop.length == ndof) {
      source_term_type=string("full");
    }
  }
  if (source_term_type==string("null") &&
      source_term_prop.length == 0) {
    source_term= &null_source_term;
  } else if (source_term_type==string("global_scalar") &&
      source_term_prop.length == 1) {
    source_term= &gscalar_source_term;
  } else if (source_term_type==string("full") &&
	     source_term_prop.length == ndof) {
    S_body.resize(1,ndof);
    source_term =  &full_source_term;
  } else {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Source term does not fit in \n"
		"a predefined name/length combination\n"
		"name entered: \"%s\"\n"
		"length entered: %d\n",
		source_term_type_s.c_str(),
		source_term_prop.length);
    assert(0);
  }

  if(LES) {      
    assert(diffusive_jacobians_type==string("scalar_per_field")||
	   diffusive_jacobians_type==string("global_scalar"));
  }

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::compute_flux(COMPUTE_FLUX_ARGS) {
  int ierr;
  double tau_a, tau_delta, gU, A01v[9];

  //  FastMat2 strain_rate(2,ndim,ndim);

  // Unfortunately we have to use copies of U and
  // iJaco due to const'ness restrictions.

  Ucpy.set(U);
  a_jac->comp_flux(flux,Ucpy);
  
  
  if (options & COMP_UPWIND) {
    
    iJaco_cpy.set(iJaco);

    // A_grad_U es ndof x 1
    // A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);
    a_jac->comp_A_grad_U(A_grad_U,grad_U);

    lam_max = 0.;
    tau_supg.set(0.);

    a_jac->comp_Uintri(Uintri,iJaco_cpy);
    a_jac->comp_vel_per_field(vel_per_field);
    if (new_adv_dif_elemset) {
      const NewAdvDif *e = new_adv_dif_elemset;
      if (e->use_ALE()) {
        a_jac->comp_vel_vec_per_field(vel_vec_per_field);
        for (int k=1; k<=ndof; k++) {
          vel_vec_per_field.ir(1,k);
          tmp16.set(vel_vec_per_field);
          tmp16.rest(e->v_mesh);
          vel_per_field.setel(tmp16.norm_p_all(),k);
          Uintri.ir(1,k).prod(iJaco_cpy,tmp16,1,-1,-1);
        }
        vel_vec_per_field.rs();
        Uintri.rs();
      }
    }

    double visco_t = 0.0;
    if (LES){

      advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
      assert(advdf_e);
      
      double Volume = advdf_e->volume();
      double Delta;
      if (ndim==2) Delta = sqrt(Volume);
      if (ndim==3) Delta = cbrt(Volume);

      strain_rate.set(grad_H);
      grad_H.t();
      strain_rate.add(grad_H).scale(0.5);
      grad_H.rs();

      double tr = (double) tmp15.prod(strain_rate,strain_rate,-1,-2,-1,-2);
      double van_D = 1.0;  // fixme using Van Driest law
      visco_t = SQ(C_smag*Delta*van_D)*sqrt(2*tr);

      if(diffusive_jacobians_type==string("global_scalar")){
	*(double *)difjac = *difjac_mol +  visco_t/Pr_t;
      } else {
	double *difjac_aux = (double *)difjac;
	for (int k=1; k<=ndof; k++) 
	  difjac_aux[k-1] = difjac_mol[k-1] + visco_t/Pr_t;	  
      }

      d_jac->update(difjac);
    }
    d_jac->comp_fluxd(fluxd,grad_U);
    d_jac->comp_dif_per_field(dif_per_field);
        
    for (int k=1; k<=ndof; k++) {
      
      Uintri.ir(1,k);
      double alpha=dif_per_field.get(k);
      double vel= vel_per_field.get(k);

      // double h_supg = 2.*vel/sqrt(Uintri.sum_square_all());
      double Uh = sqrt(Uintri.sum_square_all()); // this is
				// approx. 2*U/h
      double tau;
      FastMat2::branch();
      if (vel*vel > 20*Uh*alpha) { // remove singularity when D=0
	FastMat2::choose(0);	// magic=1
	tau = tau_fac/Uh;	// intrinsic time
      } else if (vel*vel > 1e-5*Uh*alpha) {		
        // remove singularity when v=0
	FastMat2::choose(1);
	double Pe  = vel*vel/(Uh*alpha);	// Peclet number
	// magic function
	double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
	tau = tau_fac/Uh*magic; // intrinsic time
        // printf("Pe %f, tau %f\n",Pe,tau);
      } else {
	FastMat2::choose(2);
	double h = 2./sqrt(tmp0.sum_square(iJaco,1,-1).max_all());
	tau = tau_fac*h*h/(12.*alpha);
      }
      FastMat2::leave();
      tau_supg.setel(tau,k,k); // Setting in the intrinsic time
				// matrix
      if (vel>lam_max) lam_max = vel;
    }
    Uintri.rs();
    delta_sc = 0.;
  }
  
  if (options & COMP_SOURCE) {
    c_jac->comp_G_source(G_source,Ucpy);
    G_source.scale(-1.);
    source_term->add_source_term(G_source);
  }
}
