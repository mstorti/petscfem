//__INSERT_LICENSE__
// $Id: condwallpen.cpp,v 1.12.10.1 2007/02/19 20:23:56 mstorti Exp $

#include "./condwallpen.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int CondWallRestriction::
init(int nel_a,int ndof_a,
     TextHashTable *thash,const char *name) { 
  nel = nel_a;
  ndof = ndof_a;
  assert(nel==2 || nel==4);

  int ierr;
  //o Resistance of the membrane (fixme:=
  //  so far may be only ON(R>0) / OFF(R==0) 
  TGETOPTDEF(thash,double,resistance,0.);
  //o Dimension of the problem
  TGETOPTDEF_ND(thash,int,ndim,0);
  //o Use a continuous value for #R, i.e.
  //  for #0<R<infty# we have a partial resistance
  //  to the flow. 
  TGETOPTDEF_ND(thash,int,use_partial_resistance,0);
  assert(ndof==ndim+1); // Only NS incompressible so far
  R = resistance;
  u1.resize(1,ndim);
  u2.resize(1,ndim);
  U1.resize(1,ndof);
  U2.resize(1,ndof);

  string ename = name;
  if(cond_wall_data_map.find(ename)
     != cond_wall_data_map.end()) {
    data_p = &cond_wall_data_map[ename];
    if (!MY_RANK)
      printf("in cond_wall::init(), data_p %p\n",(void*)data_p);
  }

  return 2*ndof;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void CondWallRestriction::
res(int k,FastMat2 &U,FastMat2 & r,
    FastMat2 & w,FastMat2 & jac) {
  if (use_partial_resistance)
    res_new(k,U,r,w,jac);
  else res_old(k,U,r,w,jac);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void CondWallRestriction::
res_old(int k,FastMat2 &U,FastMat2 & r,
    FastMat2 & w,FastMat2 & jac) {
  U.ir(1,1);
  U1.set(U);
  U.ir(1,2);
  U2.set(U);
  U.rs();
  if (data_p && data_p->Rv.size()>0) {
    assert(k<data_p->Rv.size());
    R = data_p->Rv.ref(k);
  }
  // printf("k %d, R %f\n",k,R);
  if (R>0) {
    u1.set(0.0);
    u2.set(0.0);
    if (data_p && data_p->u1.size()>0) {
      assert(k*ndim<data_p->u1.size());
      assert(k*ndim<data_p->u2.size());
      u1.set(&data_p->u1.e(k,0));
      u2.set(&data_p->u2.e(k,0));
    }
    // Closed
    r.set(0.0);
    U1.is(1,1,ndim);
    r.is(1,1,ndim).set(U1).minus(u1).rs();
    U1.rs();
    U2.is(1,1,ndim);
    r.rs().is(1,ndof+1,ndof+ndim)
      .set(U2).minus(u2).rs();
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
    r.set(0.).is(1,1,ndof).set(U1).minus(U2).rs();
    w.set(0.).is(3,1,ndof)
      .ir(1,1).eye()
      .ir(1,2).eye(-1.0).rs();
    jac.set(0.)
      .is(1,1,ndof).ir(2,1).eye()
      .ir(2,2).eye(-1.).rs();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void CondWallRestriction::
set_ldf(FastMat2 &ldf_user,
	vector<double> &ldf) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void CondWallRestriction::
lag_mul_dof(int jr,int &node,int &dof) {
  assert(nel==4);
  if (jr<=ndof) {
    node = 3; dof=jr;
  } else {
    node = 4; dof=jr-ndof;
  }
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void CondWallRestriction::
res(int k,FastMat2 &U,FastMat2 & r,
    FastMat2 & w,FastMat2 & jac) {
  U.ir(1,1);
  r.set(0.).is(1,1,ndof).set(U);
  U.ir(1,2);
  r.minus(U).rs();
  U.rs();
  w.set(0.).is(3,1,ndof)
    .ir(1,1).eye()
    .ir(1,2).eye(-1.0).rs();
  jac.set(0.)
    .is(1,1,ndof).ir(2,1).eye()
    .ir(2,2).eye(-1.).rs();
  
}
#endif

