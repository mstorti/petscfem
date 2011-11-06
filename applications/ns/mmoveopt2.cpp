//__INSERT_LICENSE__
//$Id$

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
#include "mmoveopt2.h"

extern int MY_RANK, SIZE;

#if 0
#undef NDEBUG
static double mypow(double x,double y) {
  assert(x>0.0);
  return exp(y*log(x));
}
#define pow mypow
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_opt2::init() {

  int ierr;
  
  PETSCFEM_ASSERT(ndof==ndim,"ndof %d, ndim %d",ndof,ndim);  
  PETSCFEM_ASSERT(ndim==2 || ndim==3,"ndim %d",ndim);  
  PETSCFEM_ASSERT(nel==ndim+1,"nel %d, ndim %d",nel,ndim);  

  dVdW.resize(2,ndim,ndim).set(0.);
  dSldW.resize(2,ndim,ndim).set(0.);
  dWdu.resize(4,ndim,ndim,nel,ndim).set(0.);
  d2VdW2.resize(4,ndim,ndim,ndim,ndim).set(0.);
  d2SldW2.resize(4,ndim,ndim,ndim,ndim).set(0.);
  w.resize(2,ndim,ndim).set(0.);
  w0.resize(2,ndim,ndim).set(0.);

  QQ.resize(2,ndim,ndim);
  D.resize(1,ndim);
  VV.resize(2,ndim,ndim);
  iTalpha.resize(2,ndim,ndim).eye();
  T1.resize(2,ndim,ndim);
  iVV.resize(2,ndim,ndim);
  tmp5.resize(2,ndim,ndim);
  tmp6.resize(2,ndim,ndim);

  tmp2.resize(2,ndim+1,ndim+1);
  xreg.resize(2,ndim+1,ndim+1);
  double xreg_v_tri[] = {0.,0.,1.0,1.0,0.,1.0,0.5,sqrt(3.0)/2.0,1.0};
  double xreg_v_tetra[] = {0.,0.,0.,1.0,
			   1.,0.,0.,1.0,
			   0.5,sqrt(3.0)/2.0,0.,1.0,
			   0.5,1.0/(sqrt(3.0)*2.0),sqrt(2.0/3.0),1.0};
  xreg.t().set(ndim==2? xreg_v_tri : xreg_v_tetra).rs();

#if 0
  xreg.print("xreg: ");
  FastMat2 a(1,ndim);
  xreg.is(1,1,ndim);
  for (int j=1; j<4; j++) {
    for (int k=1; k<4; k++) {
      xreg.ir(2,j);
      xreg.print("");
      a.set(0.0).add(xreg);
      xreg.ir(2,k);
      xreg.print("");
      a.set(0.0).minus(xreg);
      printf("length of edge: %f\n",sqrt(a.sum_square_all()));
    }
  }
  PetscFinalize();
  exit(0);
#endif

  tmp3.inv(xreg);
  tmp4.resize(2,ndim,ndim+1);

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
  TGETOPTDEF_ND(thash,double,distor_exp,-1.);

  TGETOPTDEF_ND(thash,double,volume_exp,2.);
  //o Adds a term $\propto {\tt c\_volume}\,{\rm volume}$ to the functional. 
  TGETOPTDEF_ND(thash,double,c_volume,0.);
  //o Scales distortion function
  TGETOPTDEF_ND(thash,double,c_distor,1.);
  //o Relaxation factor
  TGETOPTDEF_ND(thash,double,relax_factor,1.);
  //o If true, then the reference mesh is used as the optimal mesh. 
  TGETOPTDEF_ND(thash,double,use_ref_mesh,1.0);
  PETSCFEM_ASSERT(use_ref_mesh>=0 && use_ref_mesh<=1.0,
                  "use_ref_mesh should be in range [0,1]. use_ref_mesh: %f",
                  use_ref_mesh);  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_opt2::
element_connector(const FastMat2 &xloc,
		  const FastMat2 &state_old,
		  const FastMat2 &state_new,
		  FastMat2 &res,FastMat2 &mat) {

  double C=NAN,V=NAN,Sl=NAN,Q=NAN,Vref=NAN;
  double relax_factor_now = relax_factor;
  if (glob_param->inwt>0) 
    relax_factor_now = 1.0;

  // `y' coordinates are real
  // `x' coordinates are in the metric where the reference
  // element is `regular'
  xreg.is(1,1,ndim);
#if 0
  xref.ctr(xreg,2,1).scale(1.0-use_ref_mesh)
    .axpy(xloc,use_ref_mesh);
#else
  xref.set(xloc);
#endif
  xreg.rs();

  if (use_ref_mesh > 0.0) {
    tmp4.prod(xref,tmp3,-1,1,-1,2);
    tmp4.is(2,1,ndim);
    T1.set(tmp4);
    tmp4.rs();

    // We want a smoothly blended reference element
    // xref(alpha) that goes, from a pure regular element
    // for alpha=0 to the ref1 element for alpha=1.
    
    // T1 is the transformation from the regular to the
    // reference element, i.e. xref1 = T1 xreg + x0. 

    // Then T1 is decomposed as T1 = Q*O with
    // `O' orthogonal and `Q' spd.  Then we take the
    // ref(alpha) element as xref = T(alpha) xreg
    // with Talpha = Q^alpha * O.

    if (use_ref_mesh< 1.0) {
      // QQ = T1*T1'
      QQ.prod(T1,T1,1,-1,2,-1);
      // [V,D] = eig(QQ)
      D.seig(QQ,VV);
      iVV.inv(VV);
      // VV = VV * D^{-(1-alpha)/2}
      for (int j=1; j<=ndim; j++) {
        double dd = D.get(j);
        dd = pow(dd,-(1.0-use_ref_mesh)/2.0);
        VV.ir(2,j).scale(dd);
      }
      VV.rs();
      // tmp5 = Q^{-(1-alpha)}
      tmp5.prod(VV,iVV,1,-1,-1,2);
      // tmp6 = T(alpha) = Q^{-(1-alpha)} * T1
      tmp6.prod(tmp5,T1,1,-1,-1,2);
      // iTalpha = inv(T(alpha))
      iTalpha.inv(tmp6);
      
    } else iTalpha.inv(T1);

    y.set(xloc).add(state_new);
    y0.set(xloc).add(state_old);
    
    x.prod(iTalpha,y,2,-1,1,-1);
    x0.prod(iTalpha,y0,2,-1,1,-1);
  } else {
    x.set(xloc).add(state_new);
    x0.set(xloc).add(state_old);
  }

  if (relax_factor_now!=1.0) {
    dx.set(x).minus(x0);
    x.set(x0).axpy(dx,relax_factor_now);
  }

  for (int i=1;i<=ndim;i++){
    for (int j=1;j<=ndim;j++){
      w.setel(x.get(i+1,j)-x.get(1,j),i,j);
      w0.setel(x0.get(i+1,j)-x0.get(1,j),i,j);
    }
  }
    
  if (ndim == 2) {
    C = 4.*sqrt(3.);

    V = w.det();
    V *= 0.5;

    Vref = w0.det();
    Vref = 0.5*abs(Vref);

    vaux.norm_2(w.ir(1,1));
    Sl  = pow(double(vaux),ndim);
    vaux.norm_2(w.ir(1,2));
    Sl += pow(double(vaux),ndim);
    w.rs();
    for (int j=1;j<=ndim;j++){
      vaux1.setel(w.get(2,j)-w.get(1,j),j);
    }
    vaux.norm_2(vaux1);
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

    d2VdW2
      .setel(0.5,1,1,2,2)
      .setel(-0.5,1,2,2,1)
      .setel(-0.5,2,1,1,2)
      .setel(0.5,2,2,1,1);
    d2SldW2
      .setel(4.,1,1,1,1)
      .setel(4.,1,2,1,2)
      .setel(-2.,1,1,2,1)
      .setel(-2.,1,2,2,2);
    d2SldW2
      .setel(-2.,2,1,1,1)
      .setel(-2.,2,2,1,2)
      .setel(4.,2,1,2,1)
      .setel(4.,2,2,2,2);
  } else if (ndim == 3) {

    int ind[5]={3,1,2,3,1};

    C = 36.*sqrt(2.);

    V = w.det();
    V *= 1./6.;

    Vref = w0.det();
    Vref = 1./6.*abs(Vref);

    Sl = 0.;
    for (int i=1;i<=ndim;i++) {
      vaux.norm_2(w.ir(1,i));
      Sl += pow(double(vaux),ndim);
      w.rs();
      for (int j=1;j<=ndim;j++){
	vaux1.setel(w.get(ind[i+1],j)-w.get(ind[i],j),j);
      }
      vaux.norm_2(vaux1);
      Sl += pow(double(vaux),ndim);
    }

    dVdW.set(0.);
    d2VdW2.set(0.);

#define USE_VECMACROS
#ifdef USE_VECMACROS
    double *wp = w.storage_begin();
    double *dVdWp = dVdW.storage_begin();
    double *d2VdW2p = d2VdW2.storage_begin();
    double *dSldWp = dSldW.storage_begin();
    double *d2SldW2p = d2SldW2.storage_begin();
    double *epslcp = epsilon_LC.storage_begin();
    double *vaux1p = vaux1.storage_begin();
    double *vaux2p = vaux2.storage_begin();

#define W(j,k) VEC2(wp,j,k,ndim)
#define DVDW(j,k) VEC2(dVdWp,j,k,ndim)
#define D2VDW2(j,k,l,m) VEC4(d2VdW2p,j,k,ndim,l,ndim,m,ndim)
#define DSLDW(j,k) VEC2(dSldWp,j,k,ndim)
#define D2SLDW2(j,k,l,m) VEC4(d2SldW2p,j,k,ndim,l,ndim,m,ndim)
#define EPSLC(j,k,l) VEC3(epslcp,j,k,ndim,l,ndim)

#define VAUX1(j) (vaux1p[j])
#define VAUX2(j,k) VEC2(vaux2p,j,k,3)

    for (int p=0;p<ndim;p++) {
      for (int q=0;q<ndim;q++) {
	for (int r=0;r<ndim;r++) {
          DVDW(0,p) += EPSLC(p,q,r)*W(1,q)*W(2,r);
          DVDW(1,q) += EPSLC(p,q,r)*W(0,p)*W(2,r);
          DVDW(2,r) += EPSLC(p,q,r)*W(0,p)*W(1,q);
	  
          D2VDW2(1,q,0,p) += EPSLC(p,q,r)*W(2,r);
          D2VDW2(2,r,0,p) += EPSLC(p,q,r)*W(1,q);
          D2VDW2(0,p,1,q) += EPSLC(p,q,r)*W(2,r);
          D2VDW2(2,r,1,q) += EPSLC(p,q,r)*W(0,p);
          D2VDW2(0,p,2,r) += EPSLC(p,q,r)*W(1,q);
          D2VDW2(1,q,2,r) += EPSLC(p,q,r)*W(0,p);
	}
      }
    }
#else
    for (int p=1;p<=ndim;p++) {
      for (int q=1;q<=ndim;q++) {
	for (int r=1;r<=ndim;r++) {
	  dVdW.ir(1,1).ir(2,p)
            .set(dVdW.get(1,p)+epsilon_LC.get(p,q,r)
                 *w.get(2,q)*w.get(3,r)).rs();
	  dVdW.ir(1,2).ir(2,q)
            .set(dVdW.get(2,q)+epsilon_LC.get(p,q,r)
                 *w.get(1,p)*w.get(3,r)).rs();
	  dVdW.ir(1,3).ir(2,r)
            .set(dVdW.get(3,r)+epsilon_LC.get(p,q,r)
                 *w.get(1,p)*w.get(2,q)).rs();

	  d2VdW2.ir(1,2).ir(2,q).ir(3,1).ir(4,p)
            .set(d2VdW2.get(2,q,1,p)+epsilon_LC.get(p,q,r)*w.get(3,r)).rs();
	  d2VdW2.ir(1,3).ir(2,r).ir(3,1).ir(4,p)
            .set(d2VdW2.get(3,r,1,p)+epsilon_LC.get(p,q,r)*w.get(2,q)).rs();
	  d2VdW2.ir(1,1).ir(2,p).ir(3,2).ir(4,q)
            .set(d2VdW2.get(1,p,2,q)+epsilon_LC.get(p,q,r)*w.get(3,r)).rs();
	  d2VdW2.ir(1,3).ir(2,r).ir(3,2).ir(4,q)
            .set(d2VdW2.get(3,r,2,q)+epsilon_LC.get(p,q,r)*w.get(1,p)).rs();
	  d2VdW2.ir(1,1).ir(2,p).ir(3,3).ir(4,r)
            .set(d2VdW2.get(1,p,3,r)+epsilon_LC.get(p,q,r)*w.get(2,q)).rs();
	  d2VdW2.ir(1,2).ir(2,q).ir(3,3).ir(4,r)
            .set(d2VdW2.get(2,q,3,r)+epsilon_LC.get(p,q,r)*w.get(1,p)).rs();
	}
      }
    }
#endif
 
    dVdW.scale(1./6.);
    d2VdW2.scale(1./6.);

    dSldW.set(0.);
    d2SldW2.set(0.);
#ifdef USE_VECMACROS
    int ii,jj,kk,ll;
    for (int i=1;i<=ndim;i++) {
      ii = i-1;
      vaux2.ir(2,1).set(w.ir(1,i)).rs();
      vaux1.set(w.ir(1,ind[i+1]));
      vaux2.ir(2,2).set(vaux1.minus(w.ir(1,ind[i]))).rs();
      vaux1.set(w.ir(1,ind[i]));
      vaux2.ir(2,3).set(vaux1.minus(w.ir(1,ind[i-1]))).rs();
      for (int l=1;l<=ndim;l++) {
        ll = l-1;
        vaux2.ir(2,l);
        VAUX1(ll) = vaux2.norm_2_all();
      }
      vaux2.rs();
      w.rs();
      vaux2.rs();
      double val;
      for (int j=1;j<=ndim;j++) {
        jj = j-1;
        val = vaux1.get(1)*w.get(i,j)-vaux1.get(2)*vaux2.get(j,2)
          +vaux1.get(3)*vaux2.get(j,3);
	dSldW.setel(val,i,j);

	for (int k=1;k<=ndim;k++) {
          kk = k-1;
	  D2SLDW2(ii,jj,ii,kk) += W(ii,kk)*W(ii,jj)/VAUX1(0);
	  D2SLDW2(ii,jj,ii,kk) += VAUX2(kk,1)*VAUX2(jj,1)/VAUX1(1);
	  D2SLDW2(ii,jj,ii,kk) += VAUX2(kk,2)*VAUX2(jj,2)/VAUX1(2);
	  D2SLDW2(ii,jj,ind[i+1]-1,kk) -= VAUX2(kk,1)*VAUX2(jj,1)/VAUX1(1);
	  D2SLDW2(ii,jj,ind[i-1]-1,kk) -= VAUX2(kk,2)*VAUX2(jj,2)/VAUX1(2);

	  if (k==j) {
	    D2SLDW2(ii,jj,ii,kk) += VAUX1(0)+VAUX1(1)+VAUX1(2);
	    D2SLDW2(ii,jj,ind[i+1]-1,kk) -= VAUX1(1);
	    D2SLDW2(ii,jj,ind[i-1]-1,kk) -= VAUX1(2);
	  }
	}
      }
    }
#else
    for (int i=1;i<=ndim;i++) {
      vaux2.ir(2,1).set(w.ir(1,i)).rs();
      vaux1.set(w.ir(1,ind[i+1]));
      vaux2.ir(2,2).set(vaux1.minus(w.ir(1,ind[i]))).rs();
      vaux1.set(w.ir(1,ind[i]));
      vaux2.ir(2,3).set(vaux1.minus(w.ir(1,ind[i-1]))).rs();
      for (int l=1;l<=ndim;l++) {
	vaux1.ir(1,l).norm_2(vaux2.ir(2,l)).rs();
      }
      w.rs();
      vaux2.rs();
      for (int j=1;j<=ndim;j++) {
	dSldW.setel(vaux1.get(1)*w.get(i,j)-vaux1.get(2)
                    *vaux2.get(j,2)+vaux1.get(3)*vaux2.get(j,3),i,j);

	for (int k=1;k<=ndim;k++) {
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)
                        +w.get(i,k)*w.get(i,j)/vaux1.get(1),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)
                        +vaux2.get(k,2)*vaux2.get(j,2)
                        /vaux1.get(2),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)
                        +vaux2.get(k,3)*vaux2.get(j,3)
                        /vaux1.get(3),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,ind[i+1],k)
                        -vaux2.get(k,2)*vaux2.get(j,2)
                        /vaux1.get(2),i,j,ind[i+1],k);
	  d2SldW2.setel(d2SldW2.get(i,j,ind[i-1],k)
                        -vaux2.get(k,3)*vaux2.get(j,3)
                        /vaux1.get(3),i,j,ind[i-1],k);
	  if (k==j) {
	    d2SldW2.setel(d2SldW2.get(i,j,i,k)+vaux1.get(1)
                          +vaux1.get(2)+vaux1.get(3),i,j,i,k);
	    d2SldW2.setel(d2SldW2.get(i,j,ind[i+1],k)
                          -vaux1.get(2),i,j,ind[i+1],k);
	    d2SldW2.setel(d2SldW2.get(i,j,ind[i-1],k)
                          -vaux1.get(3),i,j,ind[i-1],k);
	  }
	}
      }
    }
#endif
    dSldW.scale(3.);
    d2SldW2.scale(3.);
  }
  if (V<=0.0) {
    printf("[%d] elem %d, vol %g\n",MY_RANK,elem,V);
    set_error(1);
  }

  Q = C*V/Sl;

  dVdu.prod(dVdW,dWdu,-1,-2,-1,-2,1,2);
  dSldu.prod(dSldW,dWdu,-1,-2,-1,-2,1,2);

  tmp.prod(d2VdW2,dWdu,1,2,-1,-2,-1,-2,3,4);
  d2Vdu2.prod(tmp,dWdu,-1,-2,3,4,-1,-2,1,2);
  tmp.prod(d2SldW2,dWdu,1,2,-1,-2,-1,-2,3,4);
  d2Sldu2.prod(tmp,dWdu,-1,-2,3,4,-1,-2,1,2);

  dQ.set(dVdu).scale(Sl).minus(dSldu.scale(V)).scale(C/pow(Sl,2));

  dSldu.scale(1./V);

  res.set(dQ).scale(distor_exp*c_distor*pow(Q,distor_exp-1.));
  res.axpy(dVdu,volume_exp*c_volume/Vref*pow(V/Vref-1.,volume_exp-1.));
  res.scale(1.0/relax_factor_now);

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
  double coef = volume_exp*c_volume/pow(Vref,volume_exp)
    *(volume_exp-1.)*pow(V-Vref,volume_exp-2.);
  mat.axpy(mat1,coef);
  coef = volume_exp*c_volume/pow(Vref,volume_exp)
    *pow(V-Vref,volume_exp-1.);
  mat.axpy(d2Vdu2,coef).scale(-1.);

  if (use_ref_mesh > 0.0) {
    res2.prod(res,iTalpha,1,-1,-1,2);
    res.set(res2);
    mat2.prod(mat,iTalpha,1,-1,3,4,-1,2);
    mat.prod(mat2,iTalpha,1,2,3,-1,-1,4);
  }
}
