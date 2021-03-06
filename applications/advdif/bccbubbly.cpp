//__INSERT_LICENSE__
//$Id: bccbubbly.cpp,v 1.3 2003/09/14 00:23:19 mstorti Exp $

extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
  
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include "nwadvdif.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int bubbly_bcconv::ask(char *,int &)"
int NewBcconv::ask(const char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   int ierr;
   NSGETOPTDEF(int,weak_form,1);
   if (!weak_form) return 0;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bubbly_bcconv::assemble"
void NewBcconv::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			     const Dofmap *dofmap,const char *jobinfo,
			     const ElementList &elemlist,
			     const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  assert(npg>0);
  assert(ndim>0);
  NSGETOPTDEF(int,weak_form,1);

  // Initialize flux functions
  int ret_options=0;
  adv_diff_ff->start_chunk(ret_options); 

  int ndimel = ndim-1;
  int nel,ndof,nelprops;
  elem_params(nel,ndof,nelprops);
  int nen = nel*ndof;

  // Unpack Dofmap
  int neq,nnod;
  neq = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),H(nH),grad_H;

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  int jdtmin;
  arg_data *staten,*stateo,*retval,*fdj_jac,*jac_prof,*Ajac;
  GlobParam *glob_param;
  FastMat2 matlocf(4,nel,ndof,nel,ndof);
  //  if (comp_mat_mass || comp_diag_mat_mass) {
  //	printf("BC_CONV is not designed for matrices assembly \n");
  //	exit(1);
  //  }

  if (comp_res) {
    int j=-1;
    stateo = &arg_data_v[++j];
    staten = &arg_data_v[++j];
    retval  = &arg_data_v[++j];
    jdtmin = ++j;
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
    if (ADVDIF_CHECK_JAC)
      fdj_jac = &arg_data_v[++j];
  }

  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matlocf.set(1.);
  }
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),
    locstate(nel,ndof),locstaten(nel,ndof),
    locstateo(nel,ndof);

  nen = nel*ndof;

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco, delta_sc;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FMatrix Jaco(ndimel,ndim),flux(ndof,ndim),fluxd(ndof,ndim),grad_U(ndim,ndof),
    A_grad_U(ndof),G_source(ndof),tau_supg(ndof,ndof),    
    fluxn(ndof),iJaco,normal(ndim),nor,lambda,Vr,Vr_inv,U(ndof);

  FastMat2 A_jac(3,ndim,ndof,ndof),D_jac(4,ndim,ndof,ndof),
    A_jac_n(2,ndof,ndof),C_jac(2,ndof,ndof);
  FastMat2 tmp1,tmp2,tmp3;

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  int start_chunk=1;
  for (ElementIterator element = elemlist.begin(); 
       element!=elemlist.end(); element++) {
    FastMat2::reset_cache();

    int k,ielh;
    element.position(k,ielh);
    // Initialize element
    adv_diff_ff->element_hook(element); 
    // Load local node coordinates in local vector
    element.node_data(nodedata,xloc.storage_begin(),
		       Hloc.storage_begin());
    
    if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
      continue;
    }

    locstateo.set(element.vector_values(*stateo)); // State at time t_n    
    locstaten.set(element.vector_values(*staten)); // State at time t_{n+1}

    // State at time t_{n+\alpha}
    locstate.set(0.).axpy(locstaten,ALPHA).axpy(locstateo,(1-ALPHA));

    matlocf.set(0.);
    veccontr.set(0.);

    // DUDA: esto no se puede sacar fuera del lazo de los elementos o es lo mismo ???
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      // normal:= normal vector times the surface of the element
      detJaco = mydetsur(Jaco,normal);
      normal.scale(-1.); // fixme:= This is to compensate a bug in mydetsur

      if (detJaco<=0.) {
	detj_error(detJaco,k);
	set_error(1);
      }

      // This is because I don't know how to use 0
      // dimension matrices. 
      if (nH>0) H.prod(SHAPE,Hloc,-1,-1,1);
      // grad_H = 0; // it is not used in this calling to flux_fun

      // state variables and gradient
      U.prod(SHAPE,locstate,-1,-1,1);
      grad_U.set(0.);  // it is not used in this calling to flux_fun

      delta_sc=0;
      double lambda_max_pg;
#if 1
      adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux,fluxd,
				A_grad_U,grad_U,G_source,
				tau_supg,delta_sc,
				lambda_max_pg, nor,lambda,Vr,Vr_inv,
				DEFAULT);
#else
      ierr =  (*adv_diff_ff)(this,U,ndim,iJaco,H,grad_H,flux,fluxd,A_jac,
			     A_grad_U,grad_U,G_source,D_jac,C_jac,tau_supg,delta_sc,
			     lambda_max_pg,nor,lambda,Vr,Vr_inv,
			     element.props(),NULL,DEFAULT,
			     start_chunk,ret_options);
#endif

      // normal = pvec(Jaco.SubMatrix(1,1,1,ndim).t(),
      // Jaco.SubMatrix(2,2,1,ndim).t());
      fluxn.prod(flux,normal,1,-1,-1);
      tmp1.prod(SHAPE,fluxn,1,2);
      veccontr.axpy(tmp1,WPG);

      adv_diff_ff->comp_A_jac_n(A_jac_n,normal);
      // A_jac_n.prod(A_jac,normal,-1,1,2,-1); // A_jac_n = (A_jac)_j n_j
      tmp2.prod(SHAPE,A_jac_n,1,2,3); //
      tmp3.prod(tmp2,SHAPE,1,2,4,3);
      matlocf.axpy(tmp3,-ALPHA*WPG);

    }

#if 0
    if (!weak_form) {
      veccontr.set(0.);
      matlocf.set(0.);
    }
#endif
    if (comp_res) {
      veccontr.export_vals(element.ret_vector_values(*retval));
      if (ADVDIF_CHECK_JAC)
        veccontr.export_vals(element.ret_fdj_values(*retval));
      if (comp_mat_each_time_step_g) 
	matlocf.export_vals(element.ret_mat_values(*Ajac));
    } else if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
