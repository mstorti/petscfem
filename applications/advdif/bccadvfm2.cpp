//__INSERT_LICENSE__
//$Id: bccadvfm2-ale.cpp,v 1.30.10.1 2007/02/23 19:18:07 rodrigop Exp $

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
#define __FUNC__ "int BcconvAdv::ask(char *,int &)"
int NewBcconv::ask(const char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   int ierr;
   NSGETOPTDEF(int,weak_form,1);
   if (!weak_form) return 0;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "NewBcconv::new_assemble"
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
  adv_diff_ff->start_chunk(ret_options); //por el termino de frontera (weak form)

  // int ndimel = ndim-1; ndimelf:ndimel de la frontera
  // es necesario ya que no se puede definir ndimelf en funcion de ndim en caso de probs 1D
  // entonces cdo tengamos bcconv's que sean puntos hay que poner el flag en elemset ndimel=0
  // en otro caso ndimelf=ndim-1=ndimel
  // Es necesario tambien conocer la dimension de elemento del cual es frontera
  // entonces debo pasarlo como param
  NSGETOPTDEF(int,ndimel,ndim);
  //o Dimension of bcconv element
  NSGETOPTDEF(int,ndimelf,ndimel-1);
  //o Reverse sense of normal. Useful if, for instance,
  //  you generated, by mistake, the wrong sense of
  //  numeration of nodes in the connectivities. 
  NSGETOPTDEF(int,reverse_normal,0);
  //o Flag to turn on ALE (Arbitrary Lagrangian-Eulerian) formulation. 
  NSGETOPTDEF(int,use_ALE_form,0);
  //o Pointer to old coordinates in
  //  #nodedata# array excluding the first "ndim" values
  NSGETOPTDEF(int,indx_ALE_xold,1);

  int  nel,ndof,nelprops;
  elem_params(nel,ndof,nelprops);
  int PFUNUSED nen = nel*ndof;
  FastMat2 bcconv_factor;

  adv_diff_ff->get_bcconv_factor(bcconv_factor);

  // Unpack Dofmap
  int nnod;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  int nH = nu-ndim;
  FastMat2  Hloc(2,nel,nH),H(1,nH),grad_H(2,ndimel,nH),
    vloc_mesh(2,nel,ndim), v_mesh(1,ndim),   
    Cp(2,ndof,ndof);

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  int PFUNUSED jdtmin;
  arg_data *staten=NULL,*stateo=NULL,*retval=NULL,
    *fdj_jac=NULL,*jac_prof=NULL,*Ajac=NULL;
  GlobParam *glob_param=NULL;
  FastMat2 matlocf(4,nel,ndof,nel,ndof);
  //  if (comp_mat_mass || comp_diag_mat_mass) {
  //	printf("BC_CONV is not designed for matrices assembly \n");
  //	exit(1);
  //  }

  double rec_Dt = NAN;
  if (comp_res) {
    int j=-1;
    stateo = &arg_data_v[++j];
    staten = &arg_data_v[++j];
    retval  = &arg_data_v[++j];
    jdtmin = ++j;
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
    rec_Dt = 1./DT;
    if (glob_param->steady) rec_Dt = 0.;
    if (ADVDIF_CHECK_JAC)
      fdj_jac = &arg_data_v[++j];
  }

  FastMat2 prof_nodes(2,nel,nel),prof_fields(2,ndof,ndof),matlocf_fix(4,nel,ndof,nel,ndof);
  FastMat2 Id_ndf(2,ndof,ndof),Id_nel(2,nel,nel),prof_fields_diag_fixed(2,ndof,ndof);

  adv_diff_ff->set_profile(prof_fields); // profile by equations (seteada en 1.)
  prof_nodes.set(1.);

  Id_ndf.eye();
  Id_nel.eye();
  
  prof_fields.d(1,2);
  prof_fields_diag_fixed.set(0.);
  prof_fields_diag_fixed.d(1,2);
  prof_fields_diag_fixed.set(prof_fields);
  prof_fields.rs();
  prof_fields_diag_fixed.rs(); //termina siendo deltaij
  prof_fields_diag_fixed.scale(-1.).add(Id_ndf);

  matlocf_fix.prod(prof_fields_diag_fixed,Id_nel,2,4,1,3);

  matlocf.prod(prof_fields,prof_nodes,2,4,1,3);
  
  if (0) {
    matlocf_fix.add(matlocf);
    if (comp_res) 
      matlocf_fix.export_vals(Ajac->profile);
    if (comp_prof) {
      jac_prof = &arg_data_v[0];
      matlocf_fix.export_vals(jac_prof->profile);
    }
  } else {
    matlocf.add(matlocf_fix);
    if (comp_res) 
      matlocf.export_vals(Ajac->profile);
    if (comp_prof) {
      jac_prof = &arg_data_v[0];
      matlocf.export_vals(jac_prof->profile);
    }
  }


  //  if (comp_prof) {
  //   jac_prof = &arg_data_v[0];
    //    matlocf.set(1.);
  // }


  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),
    xloc_new(nel,ndim), xloc_old(nel,ndim),
    locstate(nel,ndof),locstaten(nel,ndof),
    locstateo(nel,ndof);

  nen = nel*ndof;

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  // Gauss Point data
  //o Normal (only makes sense in 1D). 
  NSGETOPTDEF(double,normal_1d,0.);
  // hay que tener un elemset para los flujos entrantes y otro para los salientes
  // y se le pasa la normal por elemset
  if (ndimelf==0) { assert(normal_1d==1. || normal_1d==-1.); }
  
  //o Normal (only makes sense in 1D). 
  NSGETOPTDEF(int,use_fastmat2_cache,1);

  GPdata gp_data(geometry.c_str(),ndimelf,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double PFUNUSED detJaco, detJaco_new, detJaco_old, detJaco_mid, delta_sc;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FastMat2 Jaco(2,ndimelf,ndim),flux(2,ndof,ndimel),Jaco_new(2,ndimelf,ndim),
    fluxd(2,ndof,ndimel),grad_U(2,ndim,ndof),Jaco_old(2,ndimelf,ndim),
    A_grad_U(1,ndof),G_source(1,ndof),tau_supg(2,ndof,ndof),Jaco_mid(2,ndimelf,ndim),    
    fluxn(1,ndof),iJaco,normal(1,ndim),nor,lambda,Vr,Vr_inv,U(1,ndof),normal_mid(1,ndim),
    Halpha(1,ndof),ALE_flux(2,ndof,ndim),normal_new(1,ndim),normal_old(1,ndim);

  FastMat2 A_jac(3,ndim,ndof,ndof),D_jac(4,ndim,ndof,ndof),
    A_jac_n(2,ndof,ndof),C_jac(2,ndof,ndof);
  FastMat2 tmp1,tmp2,tmp3,tmp4;

  FastMatCacheList cache_list;
  if (use_fastmat2_cache) FastMat2::activate_cache(&cache_list);

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
    xloc_new.set(xloc);

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

    // nodal computation of mesh velocity
    if (use_ALE_form) {
      // Compute mesh velocity (at each element node)
      PETSCFEM_ASSERT0(nH >= ndim,
                       "Not enough columns for retrieving old coordinates");  
      PETSCFEM_ASSERT0(nH-indx_ALE_xold+1 >= ndim,
                       "Not enough columns for retrieving old coordinates");  
      Hloc.is(2,indx_ALE_xold,indx_ALE_xold+ndim-1);
      xloc_old.set(Hloc);
      Hloc.rs();
      xloc.scale(ALPHA).axpy(xloc_old,1-ALPHA);
      vloc_mesh.set(xloc_new).minus(xloc_old).scale(rec_Dt).rs();
    }

    // DUDA: esto no se puede sacar fuera del lazo
    // de los elementos o es lo mismo ???
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      if (ndimelf>0) {
	Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);  // xloc is at t_{n+alpha}
	detJaco = Jaco.detsur(&normal);

        if (use_ALE_form) {
          Jaco_new.prod(DSHAPEXI,xloc_new,1,-1,-1,2);
          Jaco_old.prod(DSHAPEXI,xloc_old,1,-1,-1,2);
          // normal:= normal vector times the surface of the element
          if (ndim == 2){ // we need only two points in 2D to integ temporal average
            detJaco_new = Jaco_new.detsur(&normal_new);
            detJaco_old = Jaco_old.detsur(&normal_old);
            normal.set(normal_new).add(normal_old).scale(0.5);
          } else if (ndim == 3) { // we need 3 points in 3D to integ
            // temporal average with
            // Gauss-Lobatto
            detJaco_new = Jaco_new.detsur(&normal_new);
            detJaco_old = Jaco_old.detsur(&normal_old);
            Jaco_mid.set(Jaco_new).add(Jaco_old).scale(0.5);
            detJaco_mid = Jaco_mid.detsur(&normal_mid);
            normal.set(0.).axpy(normal_new,1./6.).axpy(normal_old,1./6.)
              .axpy(normal_mid,4./6.);
          }
        }
	normal.scale(-1.); // fixme:= This is to compensate a bug in mydetsur
	if (detJaco<=0.) {
	  detj_error(detJaco,k);
	  set_error(1);
	}
      } else normal.resize(1,ndimel).set(normal_1d);
      if (reverse_normal) normal.scale(-1.);
      
      // This is because I don't know how to use 0
      // dimension matrices. 
      if (nH>0) {
	H.prod(SHAPE,Hloc,-1,-1,1);
	grad_H.set(0.0);  // no aporta al termino de contorno
      }
      
      // state variables and gradient
      U.prod(SHAPE,locstate,-1,-1,1); // U^{tn+alpha}
      grad_U.set(0.);  // it is not used in this calling to flux_fun

      delta_sc=0;
      double lambda_max_pg;
      adv_diff_ff->set_state(U,grad_U);
      adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux,fluxd,
				A_grad_U,grad_U,G_source,
				tau_supg,delta_sc,
				lambda_max_pg, nor,lambda,Vr,Vr_inv,
				DEFAULT);

      if (use_ALE_form) {
        // If ALE, then the flux must be corrected
        // by the mesh convected flux:
        // F = F - v_mesh * H
        // where H are the conservative variables
        // (generalized enthalpy)
        // Compute conservative variables
	adv_diff_ff->enthalpy_fun->enthalpy(Halpha,U);
        // Localize mesh velocity at Gauss point
        v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
        // Compute ALE flux and correct advective flux
        ALE_flux.prod(Halpha,v_mesh,1,2);
        flux.minus(ALE_flux);
      }

      // normal = pvec(Jaco.SubMatrix(1,1,1,ndim).t(),
      // Jaco.SubMatrix(2,2,1,ndim).t());
      fluxn.prod(flux,normal,1,-1,-1);
      tmp1.prod(SHAPE,fluxn,1,2);
      veccontr.axpy(tmp1,WPG);

      adv_diff_ff->comp_A_jac_n(A_jac_n,normal);

      if (use_ALE_form) {
        // If ALE, then we must modify the flux Jacobian
        // by the ALE term 
        adv_diff_ff->get_Cp(Cp);// Cp^{tn+alpha}
        tmp4.prod(v_mesh,normal,-1,-1);
        double vmeshnor = tmp4;
        A_jac_n.axpy(Cp,-vmeshnor);
      }
      
      // A_jac_n.prod(A_jac,normal,-1,1,2,-1); // A_jac_n = (A_jac)_j n_j
      tmp2.prod(SHAPE,A_jac_n,1,2,3); //
      tmp3.prod(tmp2,SHAPE,1,2,4,3);
      //      matlocf.axpy(tmp3,-ALPHA*WPG);
      matlocf.axpy(tmp3,-WPG);

    }

    // apply bcconv_factor to mask residual and matrix contributions
    for (int j=1; j<=ndof; j++) {      
      veccontr.ir(2,j).scale(bcconv_factor.get(j)).rs();
      matlocf.ir(2,j).scale(bcconv_factor.get(j)).rs();
    }

    if (comp_res) {
      veccontr.export_vals(element.ret_vector_values(*retval));
      if (ADVDIF_CHECK_JAC)
        veccontr.export_vals(element.ret_fdj_values(*fdj_jac));
      if (comp_mat_each_time_step_g) {
	matlocf.add(matlocf_fix);
	matlocf.export_vals(element.ret_mat_values(*Ajac));
      }
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
