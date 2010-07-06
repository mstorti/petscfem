//$Id: id.cpp,v 1.3.86.1 2007/02/23 19:18:07 rodrigop Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "advective.h"
#include "volabso.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int volabso::ask"
int volabso::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "volabso::assemble"
int volabso::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
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
  double *locst=NULL,*locst2=NULL,*retval=NULL,*retvalmat=NULL;

  if (comp_res) {
    int ja=0;
    locst2 = arg_data_v[ja++].locst;
    locst  = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    ja++;
    retvalmat = arg_data_v[ja++].retval;
  } 

  //o Factor affecting the whole abso term
  SGETOPTDEF(double,abso_fac,1.0);
  //o Gravity
  SGETOPTDEF(double,gravity,NAN);
  PETSCFEM_ASSERT0(!isnan(gravity),"gravity is required");  
  PETSCFEM_ASSERT0(gravity>=0.0,"gravity must be nonnegative");  
  //o Water depth reference value
  SGETOPTDEF(double,h0,NAN);
  PETSCFEM_ASSERT0(!isnan(h0),"h0 is required");  
  PETSCFEM_ASSERT0(h0>=0.0,"h0 must be nonnegative");  

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

      // veccontr.set(locstate2).scale(-id_cn).axpy(locstate,-id_cn1)
      //   .add(lumped_vc).scale(id_fac).export_vals(&(RETVAL(ielh,0,0)));
      // matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}
