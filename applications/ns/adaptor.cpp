//__INSERT_LICENSE__
//$Id: adaptor.cpp,v 1.7 2003/09/16 21:17:29 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"

extern TextHashTable *GLOBAL_OPTIONS;
   
#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
adaptor::adaptor() : elem_init_flag(0) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "nsi_tet_les_fm2::assemble"
int adaptor::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time_) {

  int kloc,node;

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(get_nearest_wall_element);

  assert(!comp_res);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsilesther\n");

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)

  //o Number of Gauss points.
  TGETOPTDEF_ND(thash,int,npg,0);
  // ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  TGETOPTDEF_ND(thash,int,ndim,0); //nd
  assert(npg>=0);  
  assert(ndim>0);
  TGETOPTDEF_ND(thash,int,ndimel,ndim);
  assert(ndimel>=0 && ndimel<=ndim);
  int nen = nel*ndof;

  // Unpack nodedata
  int nu=nodedata->nu;

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;

  if (comp_mat) retvalmat = arg_data_v[0].retval;

  double *hmin,Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    retvalmat = arg_data_v[ja++].retval;
    ja++;
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  } 

  // allocate local vecs
  nen = nel*ndof;
  FastMat2 veccontr(2,nel,ndof),xloc(2,nel,ndim),locstate(2,nel,ndof), 
    locstate2(2,nel,ndof), matlocf(4,nel,ndof,nel,ndof),
    matloc_prof(4,nel,ndof,nel,ndof);

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];
  int nprops=iprop;
  
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  
  //GPdata gp_data(geom,ndim,nel,npg);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

  // Memory allocation and initialization 
  shape.resize(2,nel,npg);
  dshapexi.resize(3,ndimel,nel,npg);
  wpg.resize(1,npg);

  for (int ipg=0; ipg<npg; ipg++) {
    wpg.setel(WPG,ipg+1);
    shape.ir(2,ipg+1).set(SHAPE);
    dshapexi.ir(3,ipg+1).set(DSHAPEXI);
  }
  shape.rs();
  dshapexi.rs();

  matloc_prof.set(1.);

  if (comp_mat_res && !elem_init_flag) {
    elem_init_flag = 1;
    for (elem=el_start; elem<=el_last; elem++) {
      if (!compute_this_elem(elem,this,myrank,iter_mode)) continue;
      element_init();
    }
  }
  if (error_code) return error_code;

  // Users may use `init()' in order to perform calculations
  // *outside* the element loop
  if (comp_mat_res) init();

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    elem=k;
    // load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    if (comp_mat_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matlocf.set(0.);
    veccontr.set(0.);

    if(comp_mat)
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));

    if (comp_mat_res) {
      // Users have to implement this function with the physics of
      // the problem.
      element_connector(xloc,locstate2,locstate,veccontr,matlocf);
      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();

  if (comp_mat_res) clean();
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
