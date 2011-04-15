//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "nsid.h"

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define ICONE(j,k) VEC2(icone,j,k,nel)
int ns_id::real_nodes(int iele,const int *&node) {
  node = &ICONE(iele,0);
  return (part_include_fic ? nel : real_nodes_m);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ns_id::initialize() {
  int ierr;
  //o Do partitioning including fictitious nodes
  TGETOPTDEF_ND(thash,int,part_include_fic,1);
  //o Pass here the number of real nodes to be reported. 
  TGETOPTDEF(thash,int,real_nodes,0);
  real_nodes_m = real_nodes;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "ns_id::assemble"
int ns_id::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		    Dofmap *dofmap,const char *jobinfo,int myrank,
		    int el_start,int el_last,int iter_mode,
		    const TimeData *time_) {

  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(get_nearest_wall_element);

  // Essentially treat comp_res as comp_mat_res but
  // with the side effect of update_jacobian=1
  int update_jacobian=1;	
  if (comp_res) {
    comp_mat_res=1;
    update_jacobian=0;
  }

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;
  // PetscPrintf(PETSCFEM_COMM_WORLD,"entrando a nsi_tet\n");

  // TGETOPTDEF(thash,int,ndim,0); //nd
  int nen = nel*ndof;

  // Unpack Dofmap
  int *ident,neq,nnod;
  neq = dofmap->neq;
  nnod = dofmap->nnod;

  // Get arguments from arg_list
  double *locst=NULL,*locst2=NULL,*retval=NULL,*retvalmat=NULL;

  if (comp_mat_res) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    if (update_jacobian) retvalmat = arg_data_v[ja++].retval;
  } 

  //__USE_PIPE_FOR_VERB_CHAR__
  //o Residual value is
  // #ns_id_fac*(ns_id_cn *(x^n-x_ref) + ns_id_cn1 * (x^{n+1}-x_ref))#
  SGETOPTDEF(double,ns_id_fac,1.);
  //o see doc for #ns_id_fac#
  SGETOPTDEF(double,ns_id_cn,0.);
  //o see doc for #ns_id_fac#
  SGETOPTDEF(double,ns_id_cn1,1.);

  //o see doc for #ns_id_fac#
  SGETOPTDEF(double,ns_id_lumped_cn,0.);
  //o see doc for #ns_id_fac#
  SGETOPTDEF(double,ns_id_lumped_cn1,0.);

  FastMat2 x_ref(2,nel,ndof),x_ref1(1,ndof);
  //o _T: double[ndof] _N: state_ref _D: null vector 
  // _DOC: Reference state value. _END
  x_ref1.set(0.);
  ierr = get_double(thash,"state_ref",x_ref1.storage_begin(),1,ndof);
  for (int j=1; j<=nel; j++) x_ref.ir(1,j).set(x_ref1);
  x_ref.rs();

  // allocate local vecs
  FastMat2 veccontr(2,nel,ndof), locstate2(2,nel,ndof),
    locstate(2,nel,ndof),matloc_prof(nen,nen),
    matlocf(2,nen,nen);
  double lumped_vc; // lumped mass matrix contribution

  if (comp_prof) {
    matloc_prof.set(1.);
  }

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;

    if (comp_mat_res || comp_res) {
      locstate.set(&(LOCST(ielh,0,0))).minus(x_ref);
      locstate2.set(&(LOCST2(ielh,0,0))).minus(x_ref);
      lumped_vc = 0;
      if (ns_id_lumped_cn) 
	lumped_vc += ns_id_lumped_cn * locstate2.sum_all();
      if (ns_id_lumped_cn1) 
	lumped_vc += ns_id_lumped_cn1 * locstate.sum_all();

      veccontr.set(locstate2).scale(-ns_id_cn).axpy(locstate,-ns_id_cn1)
	.add(lumped_vc).scale(ns_id_fac).export_vals(&(RETVAL(ielh,0,0)));
      matlocf.eye(ns_id_cn1).add(ns_id_lumped_cn1).scale(ns_id_fac);
      if (update_jacobian) 
	matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}
