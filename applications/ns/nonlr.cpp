//__INSERT_LICENSE__
/* $Id: nonlr.cpp,v 1.6 2001/06/01 03:30:50 mstorti Exp $ */

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/fastmat2.h"

#include "nsi_tet.h"

extern TextHashTable *GLOBAL_OPTIONS;

NonLinearRes::~NonLinearRes() {};

double int_pow(double base,int exp) {
  double r=1.;
  for (int j=0; j<exp; j++)
    r *= base;
  return r;
}

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
      
      R.rs().ir(1,2).is(2,1,nr).set(r).rs();
      
      matloc.ir(1,1).ir(3,2).is(4,1,nr).set(lambda).rs();
      matloc.ir(1,2).ir(3,1).is(2,1,nr).set(jac).rs();

      R.export_vals(&(RETVAL(ielh,0,0)));
      // matloc.set(0.);
      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      pp.set(matloc.storage_begin());
      pp.print("matloc: ");
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

void wall_law_res::init() {
  int ierr;
  TGETOPTNDEF(thash,int,ndim,none); //nd
  nk = ndim+2;
  ne = ndim+3;
};

void wall_law_res::res(FastMat2 & U,FastMat2 & r,
		       FastMat2 & lambda,FastMat2 & jac) {
  double k = U.get(1,nk);
  double eps = U.get(1,ne);
  double rr = eps - k;
  r.setel(rr,1);
  lambda.set(0.).setel(1.,ne,1);
  jac.setel(1.,1,nk);
  jac.setel(-1.,1,ne);
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
