//__INSERT_LICENSE__
// $Id: condwall.cpp,v 1.8 2005/04/01 02:39:16 mstorti Exp $

#include "./condwall.h"

cond_wall_data_map_t cond_wall_data_map;

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

  string ename = name();
  if(cond_wall_data_map.find(ename)
     != cond_wall_data_map.end()) {
    data_p = &cond_wall_data_map[ename];
    // printf("in cond_wall::init(), data_p %p\n",data_p);
  }
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
  if (data_p && data_p->Rv.size()>0) {
    if (k==0) assert(data_p->Rv.size()==size());
    R = data_p->Rv.ref(k);
    // printf("k %d, R %f\n",k,R);
  }
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
