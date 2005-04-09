//__INSERT_LICENSE__
// $Id: condwallpen.cpp,v 1.2 2005/04/09 11:36:52 mstorti Exp $

#include "./condwallpen.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int CondWallRestriction::
init(int nel_a,int ndof_a,TextHashTable *thash) { 
  nel = nel_a;
  ndof = ndof_a;
  assert(nel==2);
  return ndof;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void CondWallRestriction::
res(int k,FastMat2 &U,FastMat2 & r,
    FastMat2 & w,FastMat2 & jac) {
  U.ir(1,1);
  r.set(0.).is(1,1,ndof).set(U);
  U.ir(1,2);
  r.rest(U).rs();
  U.rs();
  w.set(0.).is(3,1,ndof)
    .ir(1,1).eye()
    .ir(1,2).eye(-1.0).rs();
  jac.set(0.)
    .is(1,1,ndof).ir(2,1).eye()
    .ir(2,2).eye(-1.).rs();
  
}
