// -*- mode: C++ -*-

#ifndef PETSCFEM_GODUNOV_H  
#define PETSCFEM_GODUNOV_H

#include "../../src/fem.h"
#include "advective.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/* this element solves the Riemann problem with the godunov scheme in
   a structured layer of ncells attached to a ficticious boundary.
*/
// classe para fijar caracteriticas de referencia para el prob de abs en adv-diff
class AdvDiff_Godunov : public NewElemset {
  //o Index (field numbers) of $u$ and $h$, number of dimensions
  int ndim,ndof,nu,nh,nel;
  //o Number of properties and restrictions
  int nprops;
  //o cells in godunov element
  int ncells;
  //o Gravity
  double gravity;
  //o Flux_Star tolerance
  double fs_tol;
  //o node normal 
  Property normal_prop;
  FastMat2 normal;
  //o Elementset iterator 
  ElementIterator element_m;
  //o Pointer to adv-diff flux function
  NewAdvDifFF *adv_diff_ff;
public:
  ASK_FUNCTION;
  NewAssembleFunction new_assemble;
  void init();
  //o element init
  void element_hook(ElementIterator &element);
  //o computes the residual and jacobian of the function to be imposed
  //  void comp_flux_star(ElementIterator &element,FastMat2 &U,FastMat2 &f_s);
  //o Constructor from the pointer to the flux function
  AdvDiff_Godunov(NewAdvDifFF *adv_diff_ff_=NULL) : 
    adv_diff_ff(adv_diff_ff_) {};
  //o Destructor
  ~AdvDiff_Godunov() {delete adv_diff_ff;};
};

#endif

