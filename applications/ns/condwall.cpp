//__INSERT_LICENSE__
// $Id: condwall.cpp,v 1.3 2005/03/28 16:42:53 mstorti Exp $

#include "./condwall.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int cond_wall::nres() { return ndof*2; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void cond_wall::
lag_mul_dof(int jr,int &node,int &dof) {
  if (jr<=ndof) {
    node = 3; dof=jr;
  } else {
    node = 4; dof=jr-ndof;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void cond_wall::
lm_initialize() { 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void cond_wall::
init() { 
  int nelprops;
  elem_params(nel,ndof,nelprops);
  int ierr;
  //o Resistance of the membrane (fixme:=
  //  so far may be only ON(R>0) / OFF(R==0) 
  NSGETOPTDEF(double,resistance,0.);
  //o Dimension of the problem
  NSGETOPTDEF_ND(int,ndim,0);
  assert(ndof==ndim+1); // Only NS incompressible so far
  R = resistance;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void cond_wall::
res(int k,FastMat2 &U,FastMat2 &r,
    FastMat2 &w,FastMat2 &jac) { 
  U.ir(1,1);
  U1.set(U);
  U.ir(1,2);
  U2.set(U);
  U.rs();
  if (R>0) {
    // Closed
    assert(0); // not implemented yet
  } else {
    // Open
    r.set(0.).is(1,1,ndof).set(U1).rest(U2).rs();
    w.set(0.).is(3,1,ndof)
      .ir(1,1).eye()
      .ir(1,2).eye(-1.0).rs();
    jac.set(0.)
      .is(1,1,ndof).ir(2,1).eye()
      .ir(2,2).eye(-1.)
      .rs();
    // printf("hi hi \n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void cond_wall::
element_hook(ElementIterator &element) { }
