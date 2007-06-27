//__INSERT_LICENSE__
//$Id: mmoveopt3.cpp,v 1.10.4.1 2007/02/27 01:15:22 mstorti Exp $

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
#include "mmoveopt3.h"

extern GlobParam *GLOB_PARAM;
extern double mmv_delta,mmv_d2fd,mmv_dfd;
extern double min_quality,min_volume;
extern double mmv_functional;
extern int    tarea,tangled_mesh;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_opt3::before_chunk(const char *jobinfo) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);

  PETSCFEM_ASSERT(comp_mat || comp_mat_res,
                  "Only jobinfo=\"comp_mat_res\" processed. \n"
                  "Received unrecognized jobinfo %s\n",jobinfo);  

  int ierr;

  PETSCFEM_ASSERT(ndof==ndim,
                  "Number of dofs must be equal to \n"
                  "number of dimensions for this elemset.\n"
                  "ndim %d, ndof %d\n",ndim,ndof);  
  PETSCFEM_ASSERT(ndim==2 || ndim==3,
                  "Only available for 2D/3D. ndim %d\n",
                  ndim);  
  PETSCFEM_ASSERT(nel==ndim+1,
                  "Only available for simplices \n"
                  "(triangles in 2D, tetras in 3D), ndim %d, nel %d.\n",
                  ndim,nel);  

  if (comp_mat_res) {
    res_h = get_arg_handle("res",
                           "No handle for `res'\n");
    res_delta_h = get_arg_handle("res_delta",
                                 "No handle for `res_delta'\n");
    mat_h = get_arg_handle("A","No handle for `A'\n");

    dVdW.resize(2,ndim,ndim).set(0.);
    dSldW.resize(2,ndim,ndim).set(0.);
    dWdu.resize(4,ndim,ndim,nel,ndim).set(0.);
    d2VdW2.resize(4,ndim,ndim,ndim,ndim).set(0.);
    d2SldW2.resize(4,ndim,ndim,ndim,ndim).set(0.);
    w.resize(2,ndim,ndim).set(0.);
    w0.resize(2,ndim,ndim).set(0.);

//     d2Vdu2.resize(4,nel,ndim,nel,ndim).set(0.);
//     d2Sldu2.resize(4,nel,ndim,nel,ndim).set(0.);

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
        a.set(0.0).rest(xreg);
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
  TGETOPTDEF_ND(thash,int,use_ref_mesh,1);
  //o Relaxation factor for matrix nondiagonal terms in the untangling stage
  TGETOPTDEF_ND(thash,double,relax_matrix_factor,1.0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_opt3::
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
  xref.set(xloc);
 
  if (use_ref_mesh) {
    tmp4.prod(xref,tmp3,-1,1,-1,2);
    tmp4.is(2,1,ndim);
    T0.set(tmp4);
    // T0.print("T0:");
    tmp4.rs();
    iT0.inv(T0);

    y.set(xloc).add(state_new);
    y0.set(xloc).add(state_old);
    
    x.prod(iT0,y,2,-1,1,-1);
    x0.prod(iT0,y0,2,-1,1,-1);
  } else {
    x.set(xloc).add(state_new);
    x0.set(xloc).add(state_old);
  }

  if (relax_factor_now!=1.0) {
    dx.set(x).rest(x0);
    x.set(x0).axpy(dx,relax_factor_now);
  }

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

    d2VdW2.setel(0.5,1,1,2,2).setel(-0.5,1,2,2,1)
      .setel(-0.5,2,1,1,2).setel(0.5,2,2,1,1);
    d2SldW2.setel(4.,1,1,1,1).setel(4.,1,2,1,2)
      .setel(-2.,1,1,2,1).setel(-2.,1,2,2,2);
    d2SldW2.setel(-2.,2,1,1,1).setel(-2.,2,2,1,2)
      .setel(4.,2,1,2,1).setel(4.,2,2,2,2);
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
	  dVdW.ir(1,1).ir(2,p)
            .set(dVdW.get(1,p) +epsilon_LC.get(p,q,r)
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
	dSldW.setel(vaux1.get(1)*w.get(i,j)
		    -vaux1.get(2)*vaux2.get(j,2)+vaux1.get(3)*vaux2.get(j,3),i,j);

	for (int k=1;k<=ndim;k++) {
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)
			+w.get(i,k)*w.get(i,j)/vaux1.get(1),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)
			+vaux2.get(k,2)*vaux2.get(j,2)/vaux1.get(2),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,i,k)
                        +vaux2.get(k,3)*vaux2.get(j,3)/vaux1.get(3),i,j,i,k);
	  d2SldW2.setel(d2SldW2.get(i,j,ind[i+1],k)
                        -vaux2.get(k,2)*vaux2.get(j,2)/vaux1.get(2),
                        i,j,ind[i+1],k);
	  d2SldW2.setel(d2SldW2.get(i,j,ind[i-1],k)-vaux2.get(k,3)
                        *vaux2.get(j,3)/vaux1.get(3),i,j,ind[i-1],k);
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
    dSldW.scale(3.);
    d2SldW2.scale(3.);

  }
  double el_quality = C*V/Sl;
  min_quality = (min_quality > el_quality ? el_quality : min_quality);
  min_volume  = (min_volume > V ? V : min_volume);

  // h = h(V,\delta)
  double h = 0.5*(V+pow(pow(V,2.)+4.*pow(mmv_delta,2.),0.5));
  double dh1 = 0.5*(1.+V/pow(pow(V,2.)+4.*pow(mmv_delta,2.),0.5));
  double dh2 = 2.*mmv_delta/pow(pow(V,2.)+4.*pow(mmv_delta,2.),0.5);
  double d2h11 = 2.*pow(mmv_delta,2.)/pow(pow(V,2.)+4.*pow(mmv_delta,2.),1.5);
  double d2h22 = 2.*pow(V,2.)/pow(pow(V,2.)+4.*pow(mmv_delta,2.),1.5);

  Q = C*h/Sl;

  mmv_functional += pow(Q,distor_exp);

  if (tarea > 0){

    dVdu.prod(dVdW,dWdu,-1,-2,-1,-2,1,2);
    dSldu.prod(dSldW,dWdu,-1,-2,-1,-2,1,2);

    tmp.prod(d2VdW2,dWdu,1,2,-1,-2,-1,-2,3,4);
    d2Vdu2.prod(tmp,dWdu,-1,-2,3,4,-1,-2,1,2);
    tmp.prod(d2SldW2,dWdu,1,2,-1,-2,-1,-2,3,4);
    d2Sldu2.prod(tmp,dWdu,-1,-2,3,4,-1,-2,1,2);
    
    dQ.set(dVdu).scale(dh1*Sl).rest(dSldu.scale(h)).scale(C/pow(Sl,2));

    dSldu.scale(1./h);
  
    d2Q.set(d2Vdu2).scale(Sl);
    mat1.prod(dVdu,dSldu,1,2,3,4);
    d2Q.axpy(mat1,1.);
    mat1.prod(dSldu,dVdu,1,2,3,4);
    d2Q.axpy(mat1,-1.).scale(dh1);
    d2Q.axpy(d2Sldu2,-h);
    mat1.prod(dVdu,dVdu,1,2,3,4);
    d2Q.axpy(mat1,d2h11*Sl).scale(C/pow(Sl,2));
    mat1.prod(dQ,dSldu,1,2,3,4);
    d2Q.axpy(mat1,-2./Sl);

    res.set(dQ).scale(distor_exp*c_distor*pow(Q,distor_exp-1.));
    res.axpy(dVdu,volume_exp*c_volume/Vref*pow(V/Vref-1.,volume_exp-1.));
    res.scale(1.0/relax_factor_now);
    
    mmv_dfd += distor_exp*c_distor*pow(Q,distor_exp-1.)*C/Sl*dh2;
    
    mat.prod(dQ,dQ,1,2,3,4).scale((distor_exp-1)/Q);
    mat.axpy(d2Q,1.).scale(distor_exp
			   *c_distor*pow(Q,distor_exp-1.));
    
    mat1.prod(dVdu,dVdu,1,2,3,4);
    mat.axpy(mat1,volume_exp*c_volume/pow(Vref,volume_exp)
	     *(volume_exp-1.)*pow(V-Vref,volume_exp-2.));
    mat.axpy(d2Vdu2,volume_exp*c_volume/pow(Vref,volume_exp)
	     *pow(V-Vref,volume_exp-1.)).scale(-1.);
    
    if (use_ref_mesh) {
      res2.prod(res,iT0,1,-1,-1,2);
      res.set(res2);
      mat2.prod(mat,iT0,1,-1,3,4,-1,2);
      mat.prod(mat2,iT0,1,2,3,-1,-1,4);
    }

    // Ver si esta va aca!!!
    // Relajo los terminos de acople de la matriz
    //  if (tangled_mesh){
    for (int i=1;i<=ndim;i++){
      for (int j=1;j<=ndim;j++){
	if (i!=j){
	  mat.ir(4,j).ir(2,i).scale(relax_matrix_factor);
	  mat.rs();
	}
      }
    }
    //  }

    //  res_delta.set(res).scale(2.0);
    mat2.set(dVdu).scale(-2.*mmv_delta*V*Sl/pow(pow(V,2.)+4.*pow(mmv_delta,2.),1.5));
    mat2.axpy(dSldu,dh2).scale(C/pow(Sl,2.));
    res_delta.set(dQ).scale((distor_exp-1)/Q*C/Sl*dh2);
    res_delta.axpy(mat2,-1.).scale(distor_exp*c_distor*pow(Q,distor_exp-1.)).scale(-1.);
    
    mmv_d2fd += -distor_exp*c_distor*pow(Q,distor_exp-1.)*((distor_exp-1.)/Q*pow(C/Sl*dh2,2.)+C/Sl*d2h22);

    int nen = nel*ndof;
    export_vals(res_h,res);
    export_vals(res_delta_h,res_delta);
    export_vals(mat_h,mat);
  }
}
