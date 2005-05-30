//__INSERT_LICENSE__
// $Id: condwall-new.cpp,v 1.3 2005/05/30 02:01:25 mstorti Exp $

#include "./condwall.h"
#include "./condwallpen.h"
extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void CondWallRestriction::
res_new(int k,FastMat2 &U,FastMat2 & r,
    FastMat2 & w,FastMat2 & jac) {
  if (data_p && data_p->Rv.size()>0) {
    assert(k<data_p->Rv.size());
    R = data_p->Rv.ref(k);
#if 1
    // Prints `R' table 
    static int flag=0;
    dvector<double> &Rv = data_p->Rv;
    if (!flag) {
      for (int j=0; j<Rv.size(); j++) 
	printf("j %d, R %g\n",j,Rv.ref(j));
      flag=1;
    }
#endif
  }

  U.ir(1,1);
  U1.set(U);
  U.ir(1,2);
  U2.set(U);
  U.rs();

  // The axis which is normal
  int axi=1;

  r.set(0.);
  jac.set(0.);
  w.set(0.);

  U1.is(1,1,ndim);
  U2.is(1,1,ndim);
  r.set(0.).is(1,1,ndim).set(U2)
    .rest(U1).rs();
  U1.rs();
  U2.rs();
  jac.is(1,1,ndim).ir(2,1).is(3,1,ndim).eye(-1).rs();
  jac.is(1,1,ndim).ir(2,2).is(3,1,ndim).eye(+1).rs();

  w.is(3,1,ndim).ir(1,1).is(2,1,ndim).eye(-1).rs();
  w.is(3,1,ndim).ir(1,2).is(2,1,ndim).eye(+1).rs();

  double 
    p1 = U1.get(ndof),
    p2 = U2.get(ndof),
    u1n = U1.get(axi),
    u2n = U2.get(axi);
  double res_darcy = p1-p2-R*(u1n+u2n)/2;
  r.setel(res_darcy,ndof);
  jac.ir(1,ndof)
    .setel(-R/2,1,axi).setel(+1,1,ndof)
    .setel(-R/2,2,axi).setel(-1,2,ndof).rs();
  w.setel(+1,1,ndof,ndof)
    .setel(-1,2,ndof,ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void cond_wall::
res_new(int k,FastMat2 &U,FastMat2 &r,
	FastMat2 &w,FastMat2 &jac) { 
  assert(0);
}

