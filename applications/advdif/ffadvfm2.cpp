/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>
//#include <string>

#include "../../src/fem.h"
#include "../../src/texthash.h"
#include "../../src/getprop.h"
#include "../../src/utils.h"
#include "../../src/util2.h"
#include "../../src/fastmat2.h"

#include "nwadvdif.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
newadvecfm2_ff_t::newadvecfm2_ff_t(NewAdvDif *elemset_) 
  : NewAdvDifFF(elemset_), u_per_field(*this), u_global(*this), 
  full_adv_jac(*this), full_dif_jac(*this),
  scalar_dif_per_field(*this), global_scalar_djac(*this),
  global_dif_tensor(*this), per_field_dif_tensor(*this),
  full_c_jac(*this), scalar_c_jac(*this),
  scalar_per_field_c_jac(*this),
  null_c_jac(*this)
{};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::NullCJac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  N_N_C.set(0.);
}

void newadvecfm2_ff_t::NullCJac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.set(0.);
}


void newadvecfm2_ff_t::NullCJac::
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	   FastMat2 &N,double w) {
  N_P_C.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::FullCJac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.set(N).scale(w);
  tmp.prod(tmp2,N,1,2);
  N_N_C.prod(tmp,ff.C_jac,1,3,2,4);
}

void newadvecfm2_ff_t::FullCJac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.prod(ff.C_jac,U,1,-1,-1);
}


void newadvecfm2_ff_t::FullCJac::
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	   FastMat2 &N,double w) {
  tmp26.set(P_supg).scale(w);
  tmp27.prod(N,ff.C_jac,1,2,3); // tmp27 = N * C_jac 
  N_P_C.prod(tmp26,tmp27,1,-1,2,-1,3); // tmp28 = P_supg * C_jac * N
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::ScalarCJac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.set(N).scale(w*(*ff.reacjac));
  tmp.prod(tmp2,N,1,2);
  N_N_C.prod(tmp,ff.eye_ndof,1,3,2,4);
}

void newadvecfm2_ff_t::ScalarCJac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.set(U).scale(*ff.reacjac);
}

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
void newadvecfm2_ff_t::ScalarPerFieldCjac
::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.set(N).scale(w);
  tmp.prod(tmp2,N,1,2);
  N_N_C.set(0.).d(4,2);
  N_N_C.prod(tmp,ff.C_jac,1,2,3);
  N_N_C.rs();
}

void newadvecfm2_ff_t::ScalarPerFieldCjac::
comp_G_source(FastMat2 &G_source, FastMat2 &U) {
  G_source.set(U).mult(ff.C_jac);
}

void newadvecfm2_ff_t::ScalarPerFieldCjac::
comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	   FastMat2 &N,double w) {
  tmp26.set(P_supg).scale(w);
  ff.N_C.set(0.).d(3,2).prod(N,ff.C_jac,1,2).rs();
  N_P_C.prod(tmp26,ff.N_C,1,-1,2,-1,3); // tmp28 = P_supg * C_jac * N
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::GlobalScalar
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.t().set(grad_U).scale(*(ff.difjac)).rs();
}

void newadvecfm2_ff_t::GlobalScalar
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  tmp.prod(dshapex,dshapex,-1,1,-1,2).scale(w*(*(ff.difjac)));
  grad_N_D_grad_N.prod(tmp,ff.eye_ndof,1,3,2,4);
}

void newadvecfm2_ff_t::GlobalScalar
::comp_dif_per_field(FastMat2 &dif_per_field) {
  dif_per_field.set(*(ff.difjac));
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::GlobalDifTensor
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.prod(ff.D_jac,grad_U,2,-1,-1,1);
}

void newadvecfm2_ff_t::GlobalDifTensor
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  tmp1.prod(ff.D_jac,dshapex,1,-1,-1,2).scale(w);
  tmp2.prod(dshapex,tmp1,-1,1,-1,2);
  grad_N_D_grad_N.prod(tmp2,ff.eye_ndof,1,3,2,4);
}

void newadvecfm2_ff_t::GlobalDifTensor
::comp_dif_per_field(FastMat2 &dif_per_field) {
  double dd = ff.D_jac.d(2,1).sum_all()/double(ff.elemset->ndim);
  dif_per_field.set(dd);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::PerFieldDifTensor
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  for (int k=1; k<=ff.elemset->ndof; k++) {
    ff.D_jac.ir(1,k);
    grad_U.ir(2,k);
    fluxd.ir(1,k).prod(ff.D_jac,grad_U,1,-1,-1);
  }
  fluxd.rs();
  grad_U.rs();
  ff.D_jac.rs();
}

void newadvecfm2_ff_t::PerFieldDifTensor
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.set(0.);
  for (int k=1; k<=ff.elemset->ndof; k++) {
    ff.D_jac.ir(1,k);
    tmp.prod(ff.D_jac,dshapex,1,-1,-1,2).scale(w);
    grad_N_D_grad_N.ir(2,k).ir(4,k).prod(dshapex,tmp,-1,1,-1,2);
  }
  ff.D_jac.rs();
  grad_N_D_grad_N.rs();
}

void newadvecfm2_ff_t::PerFieldDifTensor
::comp_dif_per_field(FastMat2 &dif_per_field) {
  ff.D_jac.d(3,2);
  dif_per_field.sum(ff.D_jac,1,-1).scale(1./double(ff.elemset->ndim));
  ff.D_jac.rs();
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::ScalarDifPerField
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.t().set(grad_U).rs();
  for (int j=1; j<=ff.elemset->ndof; j++) 
    fluxd.ir(1,j).scale(ff.D_jac.get(j));
  fluxd.rs();
}

void newadvecfm2_ff_t::ScalarDifPerField
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.set(0.);
  tmp.set(ff.D_jac).scale(w);
  grad_N_grad_N.prod(dshapex,dshapex,-1,1,-1,2);
  grad_N_D_grad_N.d(2,4).prod(grad_N_grad_N,tmp,1,3,2).rs();
}

void newadvecfm2_ff_t::ScalarDifPerField
::comp_dif_per_field(FastMat2 &dif_per_field) {
  dif_per_field.set(*(ff.difjac));
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::FullDifJac
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.prod(ff.D_jac,grad_U,1,-2,2,-1,-1,-2);
}

void newadvecfm2_ff_t::FullDifJac
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  D_grad_N.prod(ff.D_jac,dshapex,2,3,1,-1,-1,4)
    .scale(w);
  grad_N_D_grad_N.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
}

void newadvecfm2_ff_t::FullDifJac
::comp_dif_per_field(FastMat2 &dif_per_field) {
  ff.D_jac.d(2,1).d(4,3);
  dif_per_field.sum_abs(ff.D_jac,1,-1).scale(1./ff.elemset->ndim);
  ff.D_jac.rs();
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::FullAdvJac::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.prod(ff.u,U,2,1,-1,-1);
}

void newadvecfm2_ff_t::FullAdvJac::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.prod(ff.u,grad_U,-1,1,-2,-1,-2);
}
  
void newadvecfm2_ff_t::FullAdvJac::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  A_grad_N.prod(dshapex,ff.u,-1,1,-1,2,3);
}

void newadvecfm2_ff_t::FullAdvJac::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  // Here we take the diagonal of each Jacobian as the component
  // of a per field vector velocity
  ff.u.d(3,2);
  Uintri.prod(iJaco,ff.u,2,-1,-1,1);
  ff.u.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::UPerField::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.set(ff.u);
  for (int j=1; j<=ff.elemset->ndof; j++) {
    flux.ir(1,j).scale(U.get(j));
  }
  flux.rs();
}

void newadvecfm2_ff_t::UPerField::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  for (int j=1; j<=ff.elemset->ndof; j++) {
    grad_U.ir(2,j);
    ff.u.ir(1,j);
    tmp.prod(grad_U,ff.u,-1,-1);
    A_grad_U.setel(tmp.get(),j);
  }
  grad_U.rs();
  ff.u.rs();
  A_grad_U.rs();
}
  
void newadvecfm2_ff_t::UPerField::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  A_grad_N.set(0.).d(3,2).prod(dshapex,ff.u,-1,1,2,-1).rs();
}

void newadvecfm2_ff_t::UPerField::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  Uintri.prod(iJaco,ff.u,2,-1,1,-1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
void newadvecfm2_ff_t::UGlobal::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.prod(U,ff.u,1,2);
}

void newadvecfm2_ff_t::UGlobal::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.prod(ff.u,grad_U,-1,-1,1);
}
  
void newadvecfm2_ff_t::UGlobal::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  tmp.prod(ff.u,dshapex,-1,-1,1);
  for (int j=1; j<=ff.elemset->nel; j++) {
    A_grad_N.ir(1,j).eye(tmp.get(j));
  }
  A_grad_N.rs();
}
  
void newadvecfm2_ff_t::UGlobal::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  tmp3.prod(iJaco,ff.u,1,-1,-1);
  Uintri.prod(ff.tmp2,tmp3,1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void advecfm2_ff_t::element_hook(ElementIterator &)"
// This is to pass to the advective function the element
void newadvecfm2_ff_t::element_hook(ElementIterator &element_) {
  element = element_;

  advjac = elemset->prop_array(element,advective_jacobians_prop);
  u.set(advjac);

  difjac = elemset->prop_array(element,diffusive_jacobians_prop);
  d_jac->update(difjac);

  reacjac = elemset->prop_array(element,reactive_jacobians_prop);
  C_jac.set(reacjac);
}  

//  enum advective_jacobian_type {
//    undefined, global_vector, vector_per_field, full};

//  enum diffusive_jacobians_type {
//    undefined, global_scalar, scalar_per_field, full};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void advecfm2_ff_t::start_chunk(int ret_options)"
void newadvecfm2_ff_t::start_chunk(int ret_options) {
  int ierr;

  // FastMat2Shell *A_jac_l = elemset->A_jac;
  ndim = elemset->ndim;
  ndof = elemset->ndof;
  Uintri.resize(2,ndof,ndim);
  tmp2.resize(1,ndof).set(1.);
  tmp3.set(tmp2);
  dif_per_field.resize(1,ndof);
  eye_ndof.resize(2,ndof,ndof).eye();

  ret_options &= !SCALAR_TAU; // tell the advective element routine

  // Read advective jacobians
  //o _T: double[ndim]/double[ndim*ndof]/double[ndim*ndof*ndof] 
  //  _N: advective_jacobians _D: no default  _DOC: 
  //i_tex ../../doc/advdifop.tex advective_jacobians
  //  _END
  elemset->get_prop(advective_jacobians_prop,"advective_jacobians");

  //o Set diffusive jacobian to the desired type
  EGETOPTDEF(elemset,string,advective_jacobians_type,string("undefined"));
  string advective_jacobians_type_s = advective_jacobians_type;

  if (advective_jacobians_type==string("undefined")) {
    if (advective_jacobians_prop.length == ndim) {
      advective_jacobians_type=string("global_vector");
    } else if (advective_jacobians_prop.length == ndim*ndof) {
      advective_jacobians_type=string("vector_per_field");
    } else if (advective_jacobians_prop.length == ndim*ndof*ndof) {
      advective_jacobians_type=string("full");
    } 
  }

  if (advective_jacobians_type==string("global_vector") &&
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
    PetscPrintf(PETSC_COMM_WORLD,
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
  elemset->get_prop(diffusive_jacobians_prop,"diffusive_jacobians");

  //o Set diffusive jacobian to the desired type
  EGETOPTDEF(elemset,string,diffusive_jacobians_type,string("undefined"));
  string diffusive_jacobians_type_s=diffusive_jacobians_type;

  if (diffusive_jacobians_type==string("undefined")) {
    if (diffusive_jacobians_prop.length == 1) {
      diffusive_jacobians_type=string("global_scalar");
    } else if (diffusive_jacobians_prop.length == ndim*ndim) {
      diffusive_jacobians_type=string("global_tensor");
    } else if (diffusive_jacobians_prop.length == ndim*ndim*ndof) {
      diffusive_jacobians_type=string("per_field_tensor");
    } else if (diffusive_jacobians_prop.length == ndof) {
      diffusive_jacobians_type=string("scalar_per_field");
    } else if (diffusive_jacobians_prop.length == ndim*ndim*ndof*ndof) {
      diffusive_jacobians_type=string("full");
    }
  }
  if (diffusive_jacobians_type==string("global_scalar") &&
      diffusive_jacobians_prop.length == 1) {
    d_jac =  &global_scalar_djac;
  } else if (diffusive_jacobians_type==string("global_tensor") &&
	     diffusive_jacobians_prop.length == ndim*ndim) {
    D_jac.resize(2,ndim,ndim);
    d_jac =  &global_dif_tensor;
  } else if (diffusive_jacobians_type==string("per_field_tensor") &&
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
    PetscPrintf(PETSC_COMM_WORLD,
		"Diffusive Jacobian does not fit in \n"
		"a predefined name/length combination\n"
		"name entered: \"%s\"\n"
		"length entered: %d\n",
		diffusive_jacobians_type_s.c_str(),
		diffusive_jacobians_prop.length);
    assert(0);
  }

  // Read reactive jacobians (reactive  matrix)
  //o _T: double[var_len] 
  //  _N: reactive_jacobians _D: all zero  _DOC: 
  //  FIXME:= TO BE DOCUMENTED LATER
  //  _END
  elemset->get_prop(reactive_jacobians_prop,"reactive_jacobians");

  //o Set reactive jacobian to the desired type
  EGETOPTDEF(elemset,string,reactive_jacobians_type,string("undefined"));
  string reactive_jacobians_type_s=reactive_jacobians_type;

  if (reactive_jacobians_type==string("undefined")) {
    if (reactive_jacobians_prop.length == 1) {
      reactive_jacobians_type=string("global_scalar");
    } else if (reactive_jacobians_prop.length == ndof) {
      reactive_jacobians_type=string("scalar_per_field");
    } else if (reactive_jacobians_prop.length == ndof*ndof) {
      reactive_jacobians_type=string("full");
    } else if (reactive_jacobians_prop.length == 0) {
      reactive_jacobians_type=string("null");
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
    N_C.resize(3,elemset->nel,ndof,ndof);
    C_jac.resize(1,ndof);
    c_jac =  &scalar_per_field_c_jac;
  } else if (reactive_jacobians_type==string("full") &&
	     reactive_jacobians_prop.length == ndof*ndof) {
    C_jac.resize(2,ndof,ndof);
    c_jac =  &full_c_jac;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
		"Reactive Jacobian does not fit in \n"
		"a predefined name/length combination\n"
		"name entered: \"%s\"\n"
		"length entered: %d\n",
		reactive_jacobians_type_s.c_str(),
		reactive_jacobians_prop.length);
    assert(0);
  }

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::compute_flux(COMPUTE_FLUX_ARGS) {
  int ierr;
  double tau_a, tau_delta, gU, A01v[9];

  // Unfortunately we have to use copies of U and
  // iJaco due to const'ness restrictions.

  Ucpy.set(U);
  iJaco_cpy.set(iJaco);
  a_jac->comp_flux(flux,Ucpy);

  if (options & COMP_UPWIND) {

    d_jac->comp_fluxd(fluxd,grad_U);

    // A_grad_U es ndof x 1
    // A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);
    a_jac->comp_A_grad_U(A_grad_U,grad_U);

    lam_max = 0.;
    tau_supg.set(0.);

    a_jac->comp_Uintri(Uintri,iJaco_cpy);
    d_jac->comp_dif_per_field(dif_per_field);

    for (int k=1; k<=ndof; k++) {
      
      Uintri.ir(1,k);
      double alpha=dif_per_field.get(k);
      double vel = sqrt(u.sum_square_all());

      // double h_supg = 2.*vel/sqrt(Uintri.sum_square_all());
      double Uh = sqrt(Uintri.sum_square_all()); // this is
				// approx. 2*U/h
      double tau;
      FastMat2::branch();
      if (vel*vel > 1e-5*Uh*alpha) {		// remove singularity when v=0
	FastMat2::choose(0);
	double Pe  = vel*vel/(Uh*alpha);	// Peclet number
	// magic function
	double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
	tau = tau_fac/Uh*magic; // intrinsic time
      } else {
	FastMat2::choose(1);
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
#if 0
    G_source.set(0.);		// Only null source term allowed
				// right now!!

    C_jac.set(C_jac_l);
    G_source.prod(C_jac,U,1,-1,-1).scale(-1.);
#endif
    c_jac->comp_G_source(G_source,Ucpy);
  }
}
