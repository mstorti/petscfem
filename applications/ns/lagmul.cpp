//__INSERT_LICENSE__
/* $Id: lagmul.cpp,v 1.3 2001/10/05 12:29:13 mstorti Exp $ */

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

  int nelr,nfic; // number of real/fictitious nodes

  // Verify that the number of nodes is even
  // (a fictitious node for each real node)
  PETSCFEM_ASSERT0(nel % 2 == 0,"");
  nel2 = nel/2;
  
  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_ke);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_mat_res_ke);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

  int ierr=0;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsikeps\n");

  double *locst,*locst2,*retval,*retvalmat;
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
  // cuases poor convergence with iterative methods or amplification
  // of rounding errors. This factor scales the columns in the matrix
  // that correspond to the lagrange multipliers and can help in
  // better conditioning the system. 
  TGETOPTDEF(thash,double,lagrange_scale_factor,1.);

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,nel,ndof),R(2,nel,ndof);
  if (comp_mat) matloc_prof.set(1.);

  init();
  int nr = nres();
  FastMat2 r(1,nr),lambda(3,nel2,ndof,nr),jac(3,nr,nel2,ndof);
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
      res(k,U,r,lambda,jac);
      U.is(1,nel2+1,nel).is(2,1,nr);
      R.is(1,1,nel2).prod(lambda,U,1,2,-1,-2,-1,-2).scale(lagrange_scale_factor);
      
      R.rs().is(1,nel2+1,nel).is(2,1,nr).set(r)
	.axpy(U,-lagrange_diagonal_factor*lagrange_residual_factor).rs();
      
      matloc.is(1,1,nel2).is(3,nel2+1,nel).is(4,1,nr).set(lambda)
	.scale(-lagrange_scale_factor).rs();
      matloc.is(1,nel2+1,nel).is(3,1,nel2).is(2,1,nr).set(jac).scale(-1.).rs();
      matloc.is(1,nel2+1,nel).is(3,nel2+1,nel).d(2,4).is(2,1,nr)
	.set(lagrange_diagonal_factor).rs();

      R.export_vals(&(RETVAL(ielh,0,0)));
      // matloc.set(0.);
      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}
