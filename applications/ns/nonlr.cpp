//__INSERT_LICENSE__
/* $Id: nonlr.cpp,v 1.13 2001/06/20 02:14:53 mstorti Exp $ */

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/util2.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/fastmat2.h"

#include "nsi_tet.h"

extern TextHashTable *GLOBAL_OPTIONS;

NonLinearRes::~NonLinearRes() {};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int ns_volume_element::ask(char *,int &)"
int NonLinearRes::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "nsi_tet_les_fm2::assemble"
int NonLinearRes::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			   Dofmap *dofmap,const char *jobinfo,int myrank,
			   int el_start,int el_last,int iter_mode,
			   const TimeData *time_) {

  PETSCFEM_ASSERT0(nel==2,"");

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

  //o Using Lagrange multipliers leads to diagonal null terms This can
  // lead, to zero pitvots when using direct methods. With this option
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

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,2,ndof),R(2,2,ndof);
  if (comp_mat) matloc_prof.set(1.);

  init();
  int nr = nres();
  FastMat2 r(1,nr),lambda(2,ndof,nr),jac(2,nr,ndof),
    pp(2,nel*ndof,nel*ndof);
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
      res(U,r,lambda,jac);
      U.ir(1,2).is(2,1,nr);
      R.ir(1,1).prod(lambda,U,1,-1,-1);
      
      R.rs().ir(1,2).is(2,1,nr).set(r)
	.axpy(U,-lagrange_diagonal_factor*lagrange_residual_factor).rs();
      
      matloc.ir(1,1).ir(3,2).is(4,1,nr).set(lambda).scale(-1.).rs();
      matloc.ir(1,2).ir(3,1).is(2,1,nr).set(jac).scale(-1.).rs();
      matloc.ir(1,2).ir(3,2).d(2,4).is(2,1,nr).
	set(lagrange_diagonal_factor).rs();

      R.export_vals(&(RETVAL(ielh,0,0)));
      // matloc.set(0.);
      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      // pp.set(matloc.storage_begin());
      // pp.print("matloc: ");
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

void wall_law_res::init() {
  int ierr;
  TGETOPTDEF_ND(thash,int,ndim,0); //nd
  PETSCFEM_ASSERT0(ndim>0,"ndim should be >0\n");
  nk = ndim+2;
  ne = ndim+3;
  wf->init();

  //o The $y^+$ coordinate of the computational boundary
  TGETOPTDEF_ND(thash,double,y_wall_plus,30.);
  wf->w(y_wall_plus,fwall,fprime);
  //o C_mu (turbulence constant)
  TGETOPTDEF_ND(thash,double,C_mu,0.09);
  //o viscosity of the fluid
  TGETOPTDEF_ND(thash,double,viscosity,1.);
  //o von Karman constant (law of the wall - turbulence model)
  TGETOPTDEF_ND(thash,double,von_Karman_cnst,0.4);
  //o Do not impose the relation between k,epsilon at the wall. 
  TGETOPTDEF_ND(thash,double,turbulence_coef,1.);

  coef_k = -2./(fwall*fwall*sqrt(C_mu));
  coef_e = 1./(int_pow(fwall,4)*von_Karman_cnst*y_wall_plus*viscosity);
};

void wall_law_res::res(FastMat2 & U,FastMat2 & r,
		       FastMat2 & lambda,FastMat2 & jac) {
  U.ir(1,1);
  double u = sqrt(U.is(1,1,ndim).sum_square_all());
  U.is(1);
  if (turbulence_coef != 0.) {
    // turbulent case
    double ustar = u/fwall;
    double kw = int_pow(ustar,2)/sqrt(C_mu);
    double epsw = int_pow(ustar,4)/(von_Karman_cnst*y_wall_plus*viscosity);
  
    double k = U.get(nk);
    double eps = U.get(ne);

    // Discard k and eps eqs.
    lambda.set(0.).setel(1.,nk,1).setel(1.,ne,2);
    
    // set U to the velocity part
    U.is(2,1,ndim);

    // k eq. on the wall: k = 
    r.setel(k-kw,1);
    jac.ir(1,1).is(2,1,ndim).set(U).scale(-coef_k)
      .is(2).setel(1.,nk);

    // eps eq. on the wall
    r.setel(eps-epsw,2);
    double cc = -4*u*u/(int_pow(fwall,4)*von_Karman_cnst*y_wall_plus*viscosity);
    jac.rs().ir(1,2).is(2,1,ndim).set(U).scale(cc).is(2).setel(1.,ne);

    jac.rs();
    U.rs();

  } else {

    // laminar case
    double k = U.get(nk);
    double eps = U.get(ne);

    lambda.set(0.).setel(1.,nk,1).setel(1.,ne,2);
    r.setel(k-1.,1);
    jac.ir(1,1).is(2).setel(1.,nk);

    // eps eq. on the wall
    r.setel(eps-1.,2);
    jac.rs().ir(1,2).setel(1.,ne);

    jac.rs();
    U.rs();
  }
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
