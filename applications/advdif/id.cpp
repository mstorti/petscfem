//$Id: id.cpp,v 1.1.2.1 2003/10/16 19:07:15 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "advective.h"
#include "id.h"

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "id::assemble"
int id::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		 Dofmap *dofmap,const char *jobinfo,int myrank,
		 int el_start,int el_last,int iter_mode,
		 const TimeData *time_) {

  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_res);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;
  int nen = nel*ndof;

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;

  if (comp_res) {
    int ja=0;
    locst2 = arg_data_v[ja++].locst;
    locst  = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    ja++;
    retvalmat = arg_data_v[ja++].retval;
  } 

  //o Residual value is #id_fac*(id_cn * x^n + id_cn1 * x^{n+1})#
  SGETOPTDEF(double,id_fac,1.);
  //o see doc for #id_fac#
  SGETOPTDEF(double,id_cn,0.);
  //o see doc for #id_fac#
  SGETOPTDEF(double,id_cn1,1.);

  //o see doc for #id_fac#
  SGETOPTDEF(double,id_lumped_cn,0.);
  //o see doc for #id_fac#
  SGETOPTDEF(double,id_lumped_cn1,0.);

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

    if (comp_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
      lumped_vc = id_lumped_cn * locstate2.sum_all()
	+ id_lumped_cn1 * locstate.sum_all();

      veccontr.set(locstate2).scale(-id_cn).axpy(locstate,-id_cn1)
	.add(lumped_vc).scale(id_fac).export_vals(&(RETVAL(ielh,0,0)));
      matlocf.eye(id_cn1*id_fac);
      matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}
