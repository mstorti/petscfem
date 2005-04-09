//__INSERT_LICENSE__
/* $Id: penalize.cpp,v 1.4 2005/04/09 11:36:55 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/penalize.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "Penalize::new_assemble"
void Penalize::
new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
	     const Dofmap *dofmap,const char *jobinfo,
	     const ElementList &elemlist,
	     const TimeData *time_data) try {

  nodedata_m = nodedata_m;
  int comp_mat, comp_mat_res;
  get_comp_flags(jobinfo,comp_mat,comp_mat_res);

  if (comp_mat_res) 
    get_data(arg_data_v,stateo,staten,retval,retvalmat);
  else if (comp_mat) 
    get_data(arg_data_v,retvalmat);

  int ierr;
  //o The residual computed by the #res()# function is scaled by this. 
  NSGETOPTDEF(double,penalization_factor,0.);
  double K = penalization_factor;
  // Dimension of the embedding space (position vector of nodes)
  NSGETOPTDEF(int,ndim,0); //nd
  // Use or not caches for the FastMat2 libray
  NSGETOPTDEF(int,use_fastmat2_cache,1);

  int nelprops,ndof;
  elem_params(nel,ndof,nelprops);
  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    matloc(4,nel,ndof,nel,ndof), U(2,nel,ndof),R(2,nel,ndof);
  // xloc_m.resize(2,nel,ndim);
  if (comp_mat) matloc_prof.set(1.);

  // Call callback function defined by user initializing the elemset
  int nr = restr->init(nel,ndof,option_table());

  int nu = nodedata->nu;
  int nH = nu-ndim;

  FastMat2 r(1,nr),w(3,nel,ndof,nr),jac(3,nr,nel,ndof);
  jac.set(0.);

  //#define COMPUTE_FD_RES_JACOBIAN
#ifdef COMPUTE_FD_RES_JACOBIAN
  FastMat2 res_fd_jac(3,nr,nel,ndof),
    d_res_fd_jac(3,nr,nel,ndof),
    res_pert(1,nr),U_pert(2,nel,ndof),
    lambda_pert(3,nel,ndof,nr),fd_jac(3,nr,nel,ndof);
  res_fd_jac.set(0.0);
  res_pert.set(0.0);
  U_pert.set(0.0);
  lambda_pert.set(0.0);
  fd_jac.set(0.0);
#endif
 
  FastMatCacheList cache_list;
  if (use_fastmat2_cache)
    FastMat2::activate_cache(&cache_list);
  int ielh=-1;
  nu=nodedata->nu;

  int k_chunk;
  for (element = elemlist.begin();
       element!=elemlist.end(); element++) try {

    element.position(elem,k_chunk);
    FastMat2::reset_cache();
    ielh++;

    element_hook(element);
    if(comp_mat) {
      matloc_prof.export_vals(retvalmat->profile);
      continue;
    }      

    U.set(element.vector_values(*staten));
    matloc.set(0.);
    R.set(0.);

    if (comp_mat_res) {
      restr->res(elem,U,r,w,jac);
      R.prod(w,r,1,2,-1,-1).scale(K);
      matloc.prod(w,jac,1,2,-1,-1,3,4).scale(K);
      R.print("R:");
      matloc.print("matloc:");
      R.rs()
	.export_vals(element
			 .ret_vector_values(*retval));
      matloc.rs()
	.export_vals(element
		     .ret_mat_values(*retvalmat));
      jac.rs();
      w.rs();
    }

#ifdef COMPUTE_FD_RES_JACOBIAN
    double eps_fd=1e-4;
    for (int jele=1; jele<=nel; jele++) {
      for (int jdof=1; jdof<=ndof; jdof++) {
	U_pert.set(U);      	
	U_pert.addel(eps_fd,jele,jdof);
	res(elem,U_pert,res_pert,lambda_pert,fd_jac);
	res_pert.rest(r).scale(1./eps_fd);
	res_fd_jac.ir(2,jele).ir(3,jdof)
	  .set(res_pert).rs();
      }
    }
    d_res_fd_jac
      .set(jac).rest(res_fd_jac);
    double erro = d_res_fd_jac.sum_abs_all();
    // printf("error %g\n",erro);
#endif	

  } catch (GenericError e) {
    set_error(1);
    return;
  } 
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  close();
} catch (GenericError e) {
  set_error(1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Restriction::~Restriction() { }
