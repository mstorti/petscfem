//__INSERT_LICENSE__
/* $Id: lagmul.cpp,v 1.3 2005/01/26 18:41:16 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>
#define LagrangeMult GLagrangeMult

extern TextHashTable *GLOBAL_OPTIONS;

LagrangeMult::~LagrangeMult() {};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "LagrangeMult::assemble"
int LagrangeMult::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			   Dofmap *dofmap,const char *jobinfo,int myrank,
			   int el_start,int el_last,int iter_mode,
			   const TimeData *time_) {

  int comp_mat, comp_mat_res;
  get_comp_flags(jobinfo,comp_mat,comp_mat_res);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

  // nr:= number of restrictions
  int nr,ierr=0,jr,jfic,kfic,dofic;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsikeps\n");

  double *locst=NULL,*locst2=NULL,
    *retval=NULL,*retvalmat=NULL,lambda,rr;

  if (comp_mat_res) 
    get_data(arg_data_v,locst,retval,retvalmat);
  else if (comp_mat) 
    get_data(arg_data_v,retvalmat);

  //o Using Lagrange multipliers leads to diagonal null terms, which can
  // cause zero pivots when using direct methods. With this option
  // a small term is added to the diagonal in order to fix this. The
  // term is added only in the Jacobian or also in the residual (which
  // results would be non-consistent). See option
  //  #lagrange_residual_factor# . 
  TGETOPTDEF(thash,double,lagrange_diagonal_factor,0.);
  //o The diagonal term proportional to   #lagrange_diagonal_factor#  
  // may be also entered in the residual. If this is so
  // ( #lagrange_residual_factor=1# , then the
  // method is ``non-consistent'', i.e. the restriction is not exactly
  // satisfied by the non-linear scheme is exactly Newton-Raphson. If
  // not ( #lagrange_residual_factor=0# ) then the restriction is
  // consistently satisfied but with a non exact Newton-Raphson. 
  TGETOPTDEF(thash,double,lagrange_residual_factor,0.);
  //o Using Lagrange multipliers can lead to bad conditioning, which
  // causes poor convergence with iterative methods or amplification
  // of rounding errors. This factor scales the columns in the matrix
  // that correspond to the lagrange multipliers and can help in
  // better conditioning the system. 
  TGETOPTDEF(thash,double,lagrange_scale_factor,1.);
  //o Using Lagrange multipliers can lead to bad conditioning, which
  // causes poor convergence with iterative methods or amplification
  // of rounding errors. This factor scales the row in the matrix
  // that correspond to the new equation. 
  TGETOPTDEF(thash,double,lagrange_row_scale_factor,1.);
  // Dimension of the embedding space (position vector of nodes)
  TGETOPTDEF(thash,int,ndim,0); //nd

  // Call callback function defined by user initializing the elemset
  init();
  nr = nres(); // Get dimensions of problem

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,nel,ndof),R(2,nel,ndof);
  xloc_m.resize(2,nel,ndim);
  if (comp_mat) matloc_prof.set(1.);

  FastMat2 r(1,nr),w(3,nel,ndof,nr),jac(3,nr,nel,ndof);
  jac.set(0.);

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  int ielh=-1;
  xnod = nodedata->nodedata;
  nu=nodedata->nu;
  
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    elem = k;

    if(comp_mat) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      continue;
    }      

    U.set(&(LOCST(ielh,0,0)));
    matloc.set(0.);
    R.set(0.);

    if (comp_mat_res) {
      res(k,U,r,w,jac);
      for (jr=1; jr<=nr; jr++) {
	// get node/field of the Lag.mul.
	lag_mul_dof(jr,jfic,dofic);

	lambda = U.get(jfic,dofic);
	w.rs().ir(3,jr);
	R.axpy(w,lambda*lagrange_scale_factor);
	rr = r.get(jr) - lagrange_diagonal_factor
	  *lagrange_residual_factor*U.get(jfic,dofic);
	R.addel(rr*lagrange_row_scale_factor,jfic,dofic);

	matloc.rs().ir(3,jfic).
	  ir(4,dofic).axpy(w,-lagrange_scale_factor);
	matloc.rs().addel(lagrange_diagonal_factor
			  *lagrange_row_scale_factor,
			  jfic,dofic,jfic,dofic);
	jac.rs().ir(1,jr);
	matloc.ir(1,jfic).ir(2,dofic)
	  .axpy(jac,-lagrange_row_scale_factor);
      }
      
      R.rs().export_vals(&(RETVAL(ielh,0,0)));
      // matloc.set(0.);
      matloc.rs().export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  close();
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define NODEDATA(j,k) VEC2(xnod,j,k,nu)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Load local node coordinates in local vector
const FastMat2 &LagrangeMult::xloc() {
  for (int kloc=0; kloc<nel; kloc++) {
    int node = ICONE(elem,kloc);
    xloc_m.ir(1,kloc+1).set(&NODEDATA(node-1,0));
  }
  xloc_m.rs();
  return xloc_m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LagrangeMult::initialize() {
  lm_initialize();
}

