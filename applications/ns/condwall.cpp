//__INSERT_LICENSE__
// $Id: condwall.cpp,v 1.5 2005/03/29 04:01:50 mstorti Exp $

#include "./condwall.h"

dvector<double> cond_wall_resistance;

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
lm_initialize() { }

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
  //o Flags taking a per element resistance from
  //  global #cond_wall_resistance# vector. 
  NSGETOPTDEF_ND(int,use_vector_resistance,0);
  double r = 0.;
  if (use_vector_resistance)
    cond_wall_resistance.resize(size(),r);
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
  if (use_vector_resistance) 
    R = cond_wall_resistance.ref(k);
  if (R>0) {
    // Closed
    U1.is(1,1,ndim);
    r.is(1,1,ndim).set(U1).rs();
    U1.rs();
    U2.is(1,1,ndim);
    r.rs().is(1,ndof+1,ndof+ndim).set(U2).rs();
    U2.rs();

    w.set(0.);
    w.is(2,1,ndim).ir(1,1).is(3,1,ndim).eye().rs();
    w.is(2,1,ndim).ir(1,2).is(3,ndof+1,ndof+ndim).eye().rs();

    jac.set(0.).is(1,1,ndim).ir(2,1)
      .is(3,1,ndim).eye().rs();
    jac.is(1,ndof+1,ndof+ndim).ir(2,2)
      .is(3,1,ndim).eye().rs();
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
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void cond_wall::
element_hook(ElementIterator &element) { }