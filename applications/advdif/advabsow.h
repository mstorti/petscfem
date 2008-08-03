// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: advabso.h,v 1.11 2006/03/29 14:25:35 mstorti Exp $
#ifndef PETSCFEM_ADVABSOW_H
#define PETSCFEM_ADVABSOW_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>
#include "./advective.h"
#include "./lagmul.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class AdvectiveAbsoWall : public AdvdifLagrangeMult {
private:
  int ndim,nel,ndof;
  // Pointer to adv-diff flux fun
  NewAdvDifFF *adv_diff_ff;
  FastMat2 dummy,flux,fluxd,A_grad_U,Uref_glob,
    grad_U,Uold,normal,vmesh,vnor,A_jac,S,invS,c,
    Pi_m,Pi_p,Uref,tmp1,Ulambda,Uo,Ufluid,
    dU,Cp,invCp,unor,mask,rlam;
  int use_old_state_as_ref, use_uref_glob;
  int switch_to_ref_on_incoming;
  int ALE_flag;			
  // for ndimel~=ndim elemsets and flux functions
  int ndimel;
  // Per node normal
  Property normal_prop;
  // Velocity of mesh in ALE
  Property vmesh_prop;
  // Turn to wall boundary condition
  int turn_wall, vel_indx;
public:
  AdvectiveAbsoWall(NewAdvDifFF *ff) 
    : adv_diff_ff(ff),
      use_old_state_as_ref(0),
      turn_wall(0) {} 
  ~AdvectiveAbsoWall() { delete adv_diff_ff; } 
  int nres() { return ndof; }
  void lag_mul_dof(int jr,int &node,int &dof) {
    node = 2; dof=jr;
  }
  //element init
  virtual void 
  element_hook(ElementIterator &element);
  void lm_initialize();
  void init();
  void res(int k,FastMat2 &U,FastMat2 &r,
	   FastMat2 &w,FastMat2 &jac);
  void close() {}
};

#endif
