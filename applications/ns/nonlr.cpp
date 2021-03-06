//__INSERT_LICENSE__
/* $Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"

extern TextHashTable *GLOBAL_OPTIONS;

NonLinearRes::~NonLinearRes() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int ns_volume_element::ask(char *,int &)"
int NonLinearRes::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_prof);
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

  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_prof_ke);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_mat_res_ke);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

  int ierr=0;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsikeps\n");

  double *locst=NULL,*retval=NULL,*retvalmat=NULL;
  double PFUNUSED *locst2=NULL;
  GlobParam *glob_param=NULL;
  double PFUNUSED *hmin=NULL,Dt,rec_Dt;
  int PFUNUSED ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res || comp_mat_res_ke) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    retvalmat = arg_data_v[ja++].retval;
    hmin = &*(arg_data_v[ja++].vector_assoc)->begin();
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  } else if (comp_prof || comp_prof_ke) {
    retvalmat = arg_data_v[0].retval;
  }

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
  // cuases poor convergence with iterative methods or amplification
  // of rounding errors. This factor scales the columns in the matrix
  // that correspond to the lagrange multipliers and can help in
  // better conditioning the system. 
  TGETOPTDEF(thash,double,lagrange_scale_factor,1.);

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,2,ndof),R(2,2,ndof);
  if (comp_prof) matloc_prof.set(1.);

  init();
  int nr = nres();
  FastMat2 r(1,nr),lambda(2,ndof,nr),jac(2,nr,ndof),
    pp(2,nel*ndof,nel*ndof);
  jac.set(0.);

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  int PFUNUSED elem;
  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    elem = k;

    if(comp_prof) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      continue;
    }      

    U.set(&(LOCST(ielh,0,0)));
    matloc.set(0.);
    R.set(0.);

    if (comp_mat_res) {
      res(k,U,r,lambda,jac);
      U.ir(1,2).is(2,1,nr);
      R.ir(1,1).prod(lambda,U,1,-1,-1).scale(lagrange_scale_factor);
      
      R.rs().ir(1,2).is(2,1,nr).set(r)
	.axpy(U,-lagrange_diagonal_factor*lagrange_residual_factor).rs();
      
      matloc.ir(1,1).ir(3,2).is(4,1,nr).set(lambda)
	.scale(-lagrange_scale_factor).rs();
      matloc.ir(1,2).ir(3,1).is(2,1,nr).set(jac).scale(-1.).rs();
      matloc.ir(1,2).ir(3,2).d(2,4).is(2,1,nr)
	.set(lagrange_diagonal_factor).rs();

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
  u_wall.resize(1,ndim);
  du_wall.resize(1,ndim);

  // Physical properties
  int iprop=0; 
  //  DEFPROPN(u_wall,ndim);
  u_wall_indx = iprop; 
  ierr = get_prop(iprop,elem_prop_names,thash,elprpsindx,propel, "u_wall",(ndim));
  nprops=iprop;

  //o The $y^+$ coordinate of the computational boundary
  TGETOPTDEF_ND(thash,double,y_wall_plus,30.);
  //o The $y$ (normal) coordinate of the computational boundary. 
  // Only for laminar computations.
  TGETOPTDEF_ND(thash,double,y_wall,0.);

  //o C_mu (turbulence constant)
  TGETOPTDEF_ND(thash,double,C_mu,0.09);
  //o viscosity of the fluid
  TGETOPTDEF_ND(thash,double,viscosity,1.);
  //o Density
  SGETOPTDEF(double,rho,1.);
  //o von Karman constant (law of the wall - turbulence model)
  TGETOPTDEF_ND(thash,double,von_Karman_cnst,0.4);
  //o Do not impose the relation between k,epsilon at the wall. 
  TGETOPTDEF_ND(thash,double,turbulence_coef,1.);

  if (y_wall>0.) {
    wfs->nu = viscosity;
    wfs->y_wall = y_wall;
    wfs->rho = rho;
  } else {
    wf->w(y_wall_plus,fwall,fprime);
  }      

}

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)

void wall_law_res::res(int k,FastMat2 & U,FastMat2 & r,
		       FastMat2 & lambda,FastMat2 & jac) {

  load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
  u_wall.set(propel+u_wall_indx);
  U.ir(1,1).is(2,1,ndim);
  du_wall.set(U).minus(u_wall);
  double u = sqrt(du_wall.sum_square_all());
  double tau_w, ustar,yplus,dustar_du;
  U.is(2);
  if (turbulence_coef != 0.) {
    // turbulent case

    if (y_wall>0) {
      wfs->solve(u,ustar,tau_w,yplus,fwall,fprime,dustar_du);
    } else {
      ustar = u/fwall;
      dustar_du = 1./fwall;
    }

    double kw = int_pow(ustar,2)/sqrt(C_mu);
    double dkw_du = 2.*ustar/sqrt(C_mu)*dustar_du;

    double epsw = int_pow(ustar,4)/(von_Karman_cnst*y_wall_plus*viscosity);
    double depsw_du;
    if (y_wall>0) {
      depsw_du = 3*int_pow(ustar,2)/von_Karman_cnst/y_wall*dustar_du;
    } else {
      depsw_du = 4*epsw/ustar*dustar_du;
    }
  
    double k = U.get(nk);
    double eps = U.get(ne);

    // Discard k and eps eqs.
    lambda.set(0.).setel(1.,nk,1).setel(1.,ne,2);
    
    // set U to the velocity part
    U.is(2,1,ndim);

    // k eq. on the wall: k = 
    r.setel(k-kw,1);
    jac.ir(1,1).is(2,1,ndim).set(du_wall).scale(-dkw_du/u)
      .is(2).setel(1.,nk);

    // eps eq. on the wall
    r.setel(eps-epsw,2);
    jac.rs().ir(1,2).is(2,1,ndim).set(du_wall).scale(-depsw_du/u)
      .is(2).setel(1.,ne);

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

