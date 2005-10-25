//__INSERT_LICENSE__
//$Id: mmoveopt.cpp,v 1.1 2005/10/25 12:43:53 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/texthf.h>
#ifdef USE_NEWMAT
#include <newmatap.h>
#endif

#include "nsi_tet.h"
#include "adaptor.h"
#include "mmoveopt.h"

extern GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_opt::init() {

  int ierr;
  
  assert(ndof==ndim);
  assert(ndim==2 || ndim==3);
  assert(nel==ndim+1);

  dVdW.resize(2,ndim,ndim).set(0.);
  dSldW.resize(2,ndim,ndim).set(0.);
  dWdu.resize(4,ndim,ndim,nel,ndim).set(0.);
  d2VdW2.resize(4,ndim,ndim,ndim,ndim).set(0.);
  d2SldW2.resize(4,ndim,ndim,ndim,ndim).set(0.);
  w.resize(2,ndim,ndim).set(0.);
  w0.resize(2,ndim,ndim).set(0.);

  vaux1.resize(1,ndim).set(0.);
  vaux2.resize(2,ndim,3).set(0.);

  epsilon_LC.eps_LC();

  for(int i=1;i <= ndim;i++){
    for(int j=1;j <= ndim;j++) {
      dWdu.setel(-1.,i,j,1,j);
      for(int k=2;k <= nel;k++) {
	if (i+1 == k) {
	  dWdu.setel(1.,i,j,k,j);
	}
      }
    }
  }

  //o The functional to be minimized is $\Phi = \sum_{e=1,...,Nel} \q_e^r$,
  // where $\q_e = Vol/\sum_{i} l_i^{n_d}$,
  // and $r={\tt distor\_exp}$.
  TGETOPTDEF_ND(thash,double,distor_exp,1.);

  TGETOPTDEF_ND(thash,double,volume_exp,2.);
  //o Adds a term $\propto {\tt c\_volume}\,{\rm volume}$ to the functional. 
  TGETOPTDEF_ND(thash,double,c_volume,0.);
  //o Scales distortion function
  TGETOPTDEF_ND(thash,double,c_distor,1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_opt::
element_connector(const FastMat2 &xloc,
		  const FastMat2 &state_old,
		  const FastMat2 &state_new,
		  FastMat2 &res,FastMat2 &mat) {

  double C,V,Sl,Q,Vref;

  x.set(xloc).add(state_new);
  x0.set(xloc).add(state_old);
  for (int i=1;i<=ndim;i++){
    for (int j=1;j<=ndim;j++){
      w.setel(x.get(i+1,j)-x.get(1,j),i,j);
      w0.setel(x0.get(i+1,j)-x0.get(1,j),i,j);
    }
  }
    
  if (ndim == 2){
    C = 4.*sqrt(3.);

    V = w.det();
    V *= 0.5;

    Vref = w0.det();
    Vref = 0.5*abs(Vref);

    vaux.norm_p(w.ir(1,1),2);
    Sl  = pow(double(vaux),ndim);
    vaux.norm_p(w.ir(1,2),2);
    Sl += pow(double(vaux),ndim);
    w.rs();
    for (int j=1;j<=ndim;j++){
      vaux1.setel(w.get(2,j)-w.get(1,j),j);
    }
    vaux.norm_p(vaux1,2);
    Sl += pow(double(vaux),ndim);

    dVdW.setel(w.get(2,2),1,1);
    dVdW.setel(-1.*w.get(1,2),2,1);
    dVdW.setel(w.get(1,1),2,2);
    dVdW.setel(-1.*w.get(2,1),1,2);
    dVdW.scale(0.5);

    dSldW.setel(w.get(1,1)-vaux1.get(1),1,1);
    dSldW.setel(w.get(1,2)-vaux1.get(2),1,2);
    dSldW.setel(w.get(2,1)+vaux1.get(1),2,1);
    dSldW.setel(w.get(2,2)+vaux1.get(2),2,2);
    dSldW.scale(2.);

    d2VdW2.setel(0.5,1,1,2,2).setel(-0.5,1,2,2,1).setel(-0.5,2,1,1,2).setel(0.5,2,2,1,1);
    d2SldW2.setel(4.,1,1,1,1).setel(4.,1,2,1,2).setel(-2.,1,1,2,1).setel(-2.,1,2,2,2);
    d2SldW2.setel(-2.,2,1,1,1).setel(-2.,2,2,1,2).setel(4.,2,1,2,1).setel(4.,2,2,2,2);
  } else if (ndim == 3) {

    int ind[5]={3,1,2,3,1};

    C = 36.*sqrt(2.);

    V = w.det();
    V *= 1./6.;

    Vref = w0.det();
    Vref = 1./6.*abs(Vref);

    Sl = 0.;
    for (int i=1;i<=ndim;i++) {
      vaux.norm_p(w.ir(1,i),2);
      Sl += pow(double(vaux),ndim);
      w.rs();
      for (int j=1;j<=ndim;j++){
	vaux1.setel(w.get(ind[i+1],j)-w.get(ind[i],j),j);
      }
      vaux.norm_p(vaux1,2);
      Sl += pow(double(vaux),ndim);
    }

    dVdW.set(0.);
    d2VdW2.set(0.);
    for (int p=1;p<=ndim;p++) {
      for (int q=1;q<=ndim;q++) {
	for (int r=1;r<=ndim;r++) {
	  dVdW.ir(1,1).ir(2,p).set(dVdW.get(1,p)+epsilon_LC.get(p,q,r)*w.get(2,q)*w.get(3,r)).rs();
	  dVdW.ir(1,2).ir(2,q).set(dVdW.get(2,q)+epsilon_LC.get(p,q,r)*w.get(1,p)*w.get(3,r)).rs();
	  dVdW.ir(1,3).ir(2,r).set(dVdW.get(3,r)+epsilon_LC.get(p,q,r)*w.get(1,p)*w.get(2,q)).rs();

	  d2VdW2.ir(1,2).ir(2,q).ir(3,1).ir(4,p).set(d2VdW2.get(2,q,1,p)+epsilon_LC.get(p,q,r)*w.get(3,r)).rs();
	  d2VdW2.ir(1,3).ir(2,r).ir(3,1).ir(4,p).set(d2VdW2.get(3,r,1,p)+epsilon_LC.get(p,q,r)*w.get(2,q)).rs();
	  d2VdW2.ir(1,1).ir(2,p).ir(3,2).ir(4,q).set(d2VdW2.get(1,p,2,q)+epsilon_LC.get(p,q,r)*w.get(3,r)).rs();
	  d2VdW2.ir(1,3).ir(2,r).ir(3,2).ir(4,q).set(d2VdW2.get(3,r,2,q)+epsilon_LC.get(p,q,r)*w.get(1,p)).rs();
	  d2VdW2.ir(1,1).ir(2,p).ir(3,3).ir(4,r).set(d2VdW2.get(1,p,3,r)+epsilon_LC.get(p,q,r)*w.get(2,q)).rs();
	  d2VdW2.ir(1,2).ir(2,q).ir(3,3).ir(4,r).set(d2VdW2.get(2,q,3,r)+epsilon_LC.get(p,q,r)*w.get(1,p)).rs();
	}
      }
    }
    dVdW.scale(1./6.);
    d2VdW2.scale(1./6.);

    dSldW.set(0.);
    d2SldW2.set(0.);
    for (int i=1;i<=ndim;i++) {
      vaux2.ir(2,1).set(w.ir(1,i)).rs();
      vaux1.set(w.ir(1,ind[i+1]));
      vaux2.ir(2,2).set(vaux1.rest(w.ir(1,ind[i]))).rs();
      vaux1.set(w.ir(1,ind[i]));
      vaux2.ir(2,3).set(vaux1.rest(w.ir(1,ind[i-1]))).rs();
      for (int l=1;l<=ndim;l++) {
	vaux1.ir(1,l).norm_p(vaux2.ir(2,l),2).rs();
      }
      w.rs();
      vaux2.rs();
      for (int j=1;j<=ndim;j++) {
	dSldW.setel(vaux1.get(1)*w.get(i,j)-vaux1.get(2)*vaux2.get(j,2)+vaux1.get(3)*vaux2.get(j,3),i,j);

	for (int k=1;k<=ndim;k++) {
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)+w.get(i,k)*w.get(i,j)/vaux1.get(1),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)+vaux2.get(k,2)*vaux2.get(j,2)/vaux1.get(2),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)+vaux2.get(k,3)*vaux2.get(j,3)/vaux1.get(3),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,ind[i+1],k)-vaux2.get(k,2)*vaux2.get(j,2)/vaux1.get(2),i,j,ind[i+1],k);
	  d2SldW2.setel(d2SldW2.get(i,j,ind[i-1],k)-vaux2.get(k,3)*vaux2.get(j,3)/vaux1.get(3),i,j,ind[i-1],k);
	  if (k==j) {
	    d2SldW2.setel(d2SldW2.get(i,j,i,k)+vaux1.get(1)+vaux1.get(2)+vaux1.get(3),i,j,i,k);
	    d2SldW2.setel(d2SldW2.get(i,j,ind[i+1],k)-vaux1.get(2),i,j,ind[i+1],k);
	    d2SldW2.setel(d2SldW2.get(i,j,ind[i-1],k)-vaux1.get(3),i,j,ind[i-1],k);
	  }
	}
      }
    }
    dSldW.scale(3.);
    d2SldW2.scale(3.);

  }

  Q = C*V/Sl;

  dVdu.prod(dVdW,dWdu,-1,-2,-1,-2,1,2);
  dSldu.prod(dSldW,dWdu,-1,-2,-1,-2,1,2);

  tmp.prod(d2VdW2,dWdu,1,2,-1,-2,-1,-2,3,4);
  d2Vdu2.prod(tmp,dWdu,-1,-2,3,4,-1,-2,1,2);
  tmp.prod(d2SldW2,dWdu,1,2,-1,-2,-1,-2,3,4);
  d2Sldu2.prod(tmp,dWdu,-1,-2,3,4,-1,-2,1,2);

  dQ.set(dVdu).scale(Sl).rest(dSldu.scale(V)).scale(C/pow(Sl,2));

  dSldu.scale(1./V);

  res.set(dQ).scale(distor_exp*c_distor*pow(Q,distor_exp-1.));
  res.axpy(dVdu,volume_exp*c_volume/Vref*pow(V/Vref-1.,volume_exp-1.));

  d2Q.set(d2Vdu2).scale(Sl);
  mat1.prod(dVdu,dSldu,1,2,3,4);
  d2Q.axpy(mat1,1.);
  mat1.prod(dSldu,dVdu,1,2,3,4);
  d2Q.axpy(mat1,-1.);
  d2Q.axpy(d2Sldu2,-V).scale(C/pow(Sl,2));
  mat1.prod(dQ,dSldu,1,2,3,4);
  d2Q.axpy(mat1,-2./Sl);

  mat.prod(dQ,dQ,1,2,3,4).scale((distor_exp-1)/Q);
  mat.axpy(d2Q,1.).scale(distor_exp*c_distor*pow(Q,distor_exp-1.));

  mat1.prod(dVdu,dVdu,1,2,3,4);
  mat.axpy(mat1,volume_exp*c_volume/pow(Vref,volume_exp)*(volume_exp-1.)*pow(V-Vref,volume_exp-2.));
  mat.axpy(d2Vdu2,volume_exp*c_volume/pow(Vref,volume_exp)*pow(V-Vref,volume_exp-1.)).scale(-1.);

}
