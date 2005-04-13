//__INSERT_LICENSE__ $Id: fstepfm2.cpp,v 1.3 2002/07/26 00:57:31
//mstorti Exp $
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "fracstep.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int fracstep_fm2_cw::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat_prof);
  DONT_SKIP_JOBINFO(comp_res_mom);
  DONT_SKIP_JOBINFO(comp_mat_poi);
  DONT_SKIP_JOBINFO(comp_res_poi);
  DONT_SKIP_JOBINFO(comp_mat_prj);
  DONT_SKIP_JOBINFO(comp_res_prj);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "fracstep::assemble"
int fracstep::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		       Dofmap *dofmap,const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_) {

  assert(fractional_step);
  int ierr=0, axi;

  GET_JOBINFO_FLAG(comp_mat_prof);
  GET_JOBINFO_FLAG(comp_res_mom);
  GET_JOBINFO_FLAG(comp_mat_poi);
  GET_JOBINFO_FLAG(comp_res_poi);
  GET_JOBINFO_FLAG(comp_mat_prj);
  GET_JOBINFO_FLAG(comp_res_prj);
  GET_JOBINFO_FLAG(get_nearest_wall_element);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele) VEC2(retval,iele,0,nel*ndof)
#define RETVALMAT(iele)     VEC2(retvalmat,    iele,0,nen*nen)
#define RETVALMAT_MOM(iele) VEC2(retvalmat_mom,iele,0,nen*nen)
#define RETVALMAT_POI(iele) VEC2(retvalmat_poi,iele,0,nen*nen)
#define RETVALMAT_PRJ(iele) VEC2(retvalmat_prj,iele,0,nen*nen)

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)
#define NN_IDX(j) ELEMIPROPS_ADD(j,0)
#define IDENT(j,k) (ident[ndof*(j)+(k)]) 
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)
  
  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim,couple_velocity=0;
  ierr = get_int(thash,"couple_velocity",&couple_velocity,1); CHKERRA(ierr);
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  assert(nel==2);
  int nen = nel*ndof;

  // Unpack nodedata
  int nu=nodedata->nu;
  int nnod = dofmap->nnod;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  double *locst,*locst2,*retval,*retvalmat,*retvalmat_mom,*retvalmat_poi,
    *retvalmat_prj;
  WallData *wall_data;

  // rec_Dt is the reciprocal of Dt (i.e. 1/Dt)
  // for steady solutions it is set to 0. (Dt=inf)
  GlobParam *glob_param=NULL;
  double Dt;
  arg_data *A_mom_arg,*A_poi_arg,*A_prj_arg;
  if (comp_mat_prof) {
    int ja=0;
    retvalmat_mom = arg_data_v[ja++].retval;
    retvalmat_poi = arg_data_v[ja++].retval;
    retvalmat_prj = arg_data_v[ja++].retval;
  } else if (comp_res_mom) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    A_mom_arg = &arg_data_v[ja];
    retvalmat = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else if (comp_mat_poi) {
    int ja=0;
    A_poi_arg = &arg_data_v[ja];
    retvalmat_poi = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else if (comp_res_poi) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else if (comp_mat_prj) {
    int ja=0;
    A_prj_arg = &arg_data_v[ja];
    retvalmat_prj = arg_data_v[ja].retval;
  } else if (comp_res_prj) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    if (!reuse_mat) {
      A_prj_arg = &arg_data_v[ja];
      retvalmat_prj = arg_data_v[ja].retval;
      ja++;
    }
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else assert(0); // Not implemented yet!!

  FastMat2 veccontr(2,nel,ndof),
    locstate(2,nel,ndof),locstate2(2,nel,ndof);

  if (ndof != ndim+1) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FastMat2 matloc(4,nel,ndof,nel,ndof),
    mom_profile(4,nel,ndof,nel,ndof), du(1,ndim);

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  double penalization_factor = 1000.0;

  int ielh=-1;
  int SHV_debug=0;
#undef SHV
#define SHV(pp) { if (SHV_debug) pp.print(#pp); }
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    //if (epart[k] != myrank+1) continue;
    ielh++;

    if (comp_mat_prof) {
      matloc.set(1.);
      matloc.export_vals(&(RETVALMAT_MOM(ielh)));
      matloc.export_vals(&(RETVALMAT_PRJ(ielh)));
      matloc.export_vals(&(RETVALMAT_POI(ielh)));
      continue;
    } 

    if (comp_res_mom || comp_res_poi || comp_res_prj) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }
    veccontr.set(0.);

    if (comp_res_mom || comp_res_prj) {
      locstate.ir(1,1).is(2,1,ndim);
      du.set(locstate);
      locstate.ir(1,2);
      du.rest(locstate);
      locstate.rs();
      veccontr.rs().ir(1,1).is(2,1,ndim)
	.axpy(du,-penalization_factor);
      veccontr.rs().ir(1,2).is(2,1,ndim)
	.axpy(du,+penalization_factor);
      veccontr.export_vals(&(RETVAL(ielh)));

      matloc.set(0.).is(2,1,ndim).is(4,1,ndim)
	.ir(1,1).ir(3,1).eye()
	.ir(1,1).ir(3,2).eye(-1.)
	.ir(1,2).ir(3,1).eye(-1.)
	.ir(1,2).ir(3,2).eye(1.)
	.rs();
    } else if (comp_mat_poi) {
      matloc.set(0.).ir(2,ndof).ir(4,ndof)
	.setel(+1.0,1,1)
	.setel(-1.0,1,2)
	.setel(-1.0,2,1)
	.setel(+1.0,2,2);
      matloc.export_vals(&(RETVALMAT_POI(ielh)));
    } else if (comp_res_poi) {
      veccontr.rs().export_vals(&(RETVAL(ielh)));
    } else if (comp_mat_prj) {
      matloc.export_vals(&(RETVALMAT_PRJ(ielh)));
    } else assert(0);
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
