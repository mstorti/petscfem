//__INSERT_LICENSE__
/* $Id: lagmul.cpp,v 1.8 2002/09/16 00:16:25 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"

extern TextHashTable *GLOBAL_OPTIONS;

LagrangeMult::~LagrangeMult() {};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int LagrangeMult::ask(char *,int &)"
int LagrangeMult::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "LagrangeMult::assemble"
int LagrangeMult::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			   Dofmap *dofmap,const char *jobinfo,int myrank,
			   int el_start,int el_last,int iter_mode,
			   const TimeData *time_) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_ke);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_mat_res_ke);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

  // nr:= number of restrictions
  int nr,ierr=0,jr,jfic,kfic,dofic;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsikeps\n");

  double *locst,*locst2,*retval,*retvalmat,lambda,rr;
  GlobParam *glob_param;
  double *hmin,Dt,rec_Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res || comp_mat_res_ke) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    retvalmat = arg_data_v[ja++].retval;
    hmin = (arg_data_v[ja++].vector_assoc)->begin();
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  } else if (comp_mat || comp_mat_ke) {
    retvalmat = arg_data_v[0].retval;
  }

  //o Using Lagrange multipliers leads to diagonal null terms, which can
  // cause zero pivots when using direct methods. With this option
  // a small term is added to the diagonal in order to fix this. The
  // term is added only in the Jacobian or also in the residual (which
  // results would be non-consistent). See option
  // \verb+lagrange_residual_factor+. 
  TGETOPTDEF(thash,double,lagrange_diagonal_factor,0.);
  //o The diagonal term proportional to  \verb+lagrange_diagonal_factor+ 
  // may be also entered in the residual. If this is so
  // (\verb+lagrange_residual_factor=1+, then the
  // method is ``non-consistent'', i.e. the restriction is not exactly
  // satisfied by the non-linear scheme is exactly Newton-Raphson. If
  // not (\verb+lagrange_residual_factor=0+) then the restriction is
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

  init();
  nr = nres(); // Get dimensions of problem

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,nel,ndof),R(2,nel,ndof);
  if (comp_mat) matloc_prof.set(1.);

  FastMat2 r(1,nr),w(3,nel,ndof,nr),jac(3,nr,nel,ndof);
  jac.set(0.);

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  int elem;
  int ielh=-1;
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
  return 0;
}
