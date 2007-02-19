//$Id: godunov.cpp,v 1.3.56.1 2007/02/19 20:23:56 mstorti Exp $
#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
//#include <src/cloud.h>
#include <src/elemset.h>

#include "advective.h"
#include "godunov.h"
#include "nwadvdif.h"

//extern Mesh *GLOBAL_MESH;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int NewAdvDif::ask(char *,int &)"
int AdvDiff_Godunov::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_prof);
  DONT_SKIP_JOBINFO(comp_res);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void NewAdvDif::assemble"
void AdvDiff_Godunov::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
				    const Dofmap *dofmap,const char *jobinfo,
				    const ElementList &elemlist,
				    const TimeData *time_data) {
  int nelprops,node,dof,ierr=0;
  elem_params(nel,ndof,nelprops);
  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_res);
  NSGETOPTDEF(int,ndim,0);
  assert(ndim>0);
  int nu=nodedata->nu;
  int nH = nu-ndim;
  int ret_options=0;
  assert(nel==2);
  //  adv_diff_ff->start_chunk(ret_options); //ff ini
  
  //aqui en advdif esta primero en el arg_data_v el estado en t_n y despues en t_n+1.
  arg_data *locstold,*locst,*retval,*fdj_jac,*jac_prof,*Ajac;
  GlobParam *glob_param;
  double *hmin,Dt,rec_Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_res) {
    int ja=-1;
    locstold = &arg_data_v[++ja];
    locst = &arg_data_v[++ja];
    retval = &arg_data_v[++ja];
    ++ja; // Not used. Should be used for computation of automatic Dt
    Ajac = &arg_data_v[++ja];
    glob_param = (GlobParam *)(arg_data_v[++ja].user_data);
#ifdef CHECK_JAC
    fdj_jac = &arg_data_v[++ja];
#endif
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.0;
  } else if (comp_prof) {
    jac_prof = &arg_data_v[0];
  }
  
  FastMat2 matlocf(4,nel,ndof,nel,ndof),matlocf_mass(4,nel,ndof,nel,ndof);
  FastMat2 prof_nodes(2,nel,nel),prof_fields(2,ndof,ndof),matlocf_fix(4,nel,ndof,nel,ndof);
  FastMat2 Id_ndf(2,ndof,ndof),Id_nel(2,nel,nel),prof_fields_diag_fixed(2,ndof,ndof);
  
  //o Use the weak form for the Galerkin part of the advective term.
  NSGETOPTDEF(int,weak_form,1);
  //o Weights the temporal term with $N+\beta P$, i.e.
  //o $\beta=0$ is equivalent to waight the temporal term a la
  //o Galerkin and $\beta=1$ is equivalent to do the consistent SUPG weighting.
  NSGETOPTDEF(double,beta_supg,1.);
  //o Use lumped mass (used mainly to avoid oscillations for small time steps).
  NSGETOPTDEF(int,lumped_mass,0);
  //o Add a shock capturing term
  NSGETOPTDEF(double,shocap,0.0);
  
  // Initialize flux functions
  ff_options=0;
  adv_diff_ff->start_chunk(ff_options);
  int ndimel = 1;//adv_diff_ff->dim();
  FMatrix grad_H(ndimel,nH);
  
  adv_diff_ff->set_profile(prof_fields); // profile by equations
  prof_nodes.set(1.);
  
  Id_ndf.eye();
  Id_nel.eye();
  
  prof_fields.d(1,2);
  prof_fields_diag_fixed.set(0.);
  prof_fields_diag_fixed.d(1,2);
  prof_fields_diag_fixed.set(prof_fields);
  prof_fields.rs();
  prof_fields_diag_fixed.rs();
  prof_fields_diag_fixed.scale(-1.).add(Id_ndf);
  
  matlocf_fix.prod(prof_fields_diag_fixed,Id_nel,2,4,1,3);
  matlocf.prod(prof_fields,prof_nodes,2,4,1,3);
  
  matlocf.add(matlocf_fix);
  if (comp_res)
    matlocf.export_vals(Ajac->profile);
  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matlocf.export_vals(jac_prof->profile);
  }
  int nlog_vars;
  const int *log_vars;
  adv_diff_ff->get_log_vars(nlog_vars,log_vars);
  //o Use log-vars for $k$ and $\epsilon$
  NSGETOPTDEF(int,use_log_vars,0);
  if (!use_log_vars) nlog_vars=0;
  
  // Not implemented yet:= not lumped_mass + log-vars
  assert(!use_log_vars || lumped_mass);
  // lumped_mass:= If this options is activated then all the inertia
  // term matrix comtributions are added to 'matlocf_mass' and the
  // vector contribution terms are discarded. Then at the last moment
  // matlocf_mass*(Un-Uo)/Dt is added.
  
  // Allocate local vecs
  FMatrix veccontr(nel,ndof),veccontr_mass(nel,ndof),
    xloc(nel,ndim),lstate(nel,ndof),
    lstateo(nel,ndof),lstaten(nel,ndof),dUloc_c(nel,ndof),
    dUloc(nel,ndof),matloc;
  FastMat2 true_lstate(2,nel,ndof),
    true_lstateo(2,nel,ndof),true_lstaten(2,nel,ndof);
  
  nen = nel*ndof;
  
  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);
  
  double detJaco, wpgdet, delta_sc;
  int elem, ipg,node, jdim, kloc,lloc,ldof;
  
  dshapex.resize(2,ndimel,nel);
  FMatrix Jaco(ndimel,ndim),Jaco_av(ndimel,ndim),
    iJaco(ndimel,ndimel),
    flux(ndof,ndimel),fluxd(ndof,ndimel),mass(nel,nel),
    grad_U(ndimel,ndof), A_grad_U(ndof),
    G_source(ndof), dUdt(ndof), Un(ndof),
    Ho(ndof),Hn(ndof);
  // These are declared but not used
  FMatrix nor,lambda,Vr,Vr_inv,U(ndof),Ualpha(ndof),
    lmass(nel),Id_ndof(ndof,ndof),
    tmp1,tmp2,tmp3,tmp4,tmp5,hvec(ndimel),tmp6,tmp7,
    tmp8,tmp9,tmp10,tmp11(ndof,ndimel),tmp12,tmp14,
    tmp15,tmp17,tmp19,tmp20,tmp21,tmp22,tmp23,
    tmp24,tmp_sc,tmp_sc_v;
  FastMat2 A_grad_N(3,nel,ndof,ndof),
    grad_N_D_grad_N(4,nel,ndof,nel,ndof),N_N_C(4,nel,ndof,nel,ndof),
    N_P_C(3,ndof,nel,ndof),N_Cp_N(4,nel,ndof,nel,ndof),
    P_Cp(2,ndof,ndof);
  Ao_grad_N.resize(3,nel,ndof,ndof);
  tau_supg.resize(2,ndof,ndof);
  P_supg.resize(3,nel,ndof,ndof);
  
  FMatrix Jaco_axi(2,2);
  int ind_axi_1, ind_axi_2;
  double detJaco_axi;
  
  //Elementser initialization
  init();
  
#define COMPUTE_FD_ADV_JACOBIAN
#ifdef COMPUTE_FD_ADV_JACOBIAN
  FastMat2 A_fd_jac(3,ndimel,ndof,ndof),U_pert(1,ndof),flux_pert(2,ndof,ndimel);
#endif
  FastMat2 U0(1,ndof);
  double u0[] = {1.,1.,1.,2.,0.,0.,0.1,0.1};
  U0.set(u0);
  
  Uo.resize(1,ndof);
  Id_ndof.set(0.);
  for (int j=1; j<=ndof; j++) Id_ndof.setel(1.,j,j);
  
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
  
  for (ElementIterator element = elemlist.begin();
       element!=elemlist.end(); element++) try {
	 
	 FastMat2::reset_cache();
	 
	 // Initialize element
	 adv_diff_ff->element_hook(element);
	 // Get nodedata info (coords. etc...)
	 element.node_data(nodedata,xloc.storage_begin(),
			   Hloc.storage_begin());
	 //Element subdiscretization
	 double e_lenght=sqrt(SQ(xloc.get(1,1)-xloc.get(2,1))+ \
			      SQ(xloc.get(1,2)-xloc.get(2,2)));
	 double cell_lenght=e_lenght/ncells;
	 
	 if (comp_prof) {
	   matlocf.export_vals(element.ret_mat_values(*jac_prof));
	   continue;
	 }
	 
	 if (comp_res) {
	   lambda_max=0;
	   lstateo.set(element.vector_values(*stateo));
	   lstaten.set(element.vector_values(*staten));
	 }
	 
	 // State at time t_{n+\alpha}
	 lstate.set(0.).axpy(lstaten,ALPHA).axpy(lstateo,(1-ALPHA));
	 log_transf(true_lstate ,lstate ,nlog_vars,log_vars);
	 
	 veccontr.set(0.);
	 mass.set(0.);
	 lmass.set(0.);
	 matlocf.set(0.);
	 if (lumped_mass) matlocf_mass.set(0.);
	 
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])
	 // loop over Gauss points
	 Jaco_av.set(0.);
	 for (ipg=0; ipg<npg; ipg++) {
	   Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
	   Jaco_av.add(Jaco);
	   
	   if (ndim==ndimel) {
	     iJaco.inv(Jaco);
	     detJaco = Jaco.det();
	   } else if (ndimel==1) {
	     // This allows to solve problems on streams like rivers or
	     // ducts or advective problems on plane surfaces (not
	     // implemented yet). Also, it could be used also for advective
	     // problems on arbitrary surfaces (ndim=3 and ndimel=2) but I
	     // don't know how to do that yet. (tensorial calculus...)
	     detJaco = Jaco.norm_p_all(2);
	     iJaco.setel(1./detJaco,1,1);
	   }
	   if (detJaco<=0.) {
	     int k,ielh;
	     element.position(k,ielh);
	     detj_error(detJaco,k);
	     set_error(1);
	   }
	   wpgdet = detJaco*WPG;
	   
	   // Set volume of element. This is somewhat incorrect
	   // because we should loop over Gauss points.
	   // Correct for simplices (tris and tetras) and not
	   // deformed quads and hexas.
	   Volume = double(npg)*wpgdet;
	   volume_flag = 1;
	   
	   // This is incorrect. Master elment volume is included in the
	   // Gauss point weight.
	   // Volume = double(npg)*wpgdet/gp_data.master_volume;
	   dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);
	   
	   if (nH>0) {
	     H.prod(SHAPE,Hloc,-1,-1,1);
	     grad_H.prod(dshapex,Hloc,1,-1,-1,2);
	   }
	   
	   if (comp_res) {
	     // state variables and gradient
	     Un.prod(SHAPE,lstaten,-1,-1,1);
	     adv_diff_ff->enthalpy_fun->enthalpy(Hn,Un);
	     Uo.prod(SHAPE,lstateo,-1,-1,1);
	     adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
	     Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
	     dUdt.set(Hn).rest(Ho).scale(rec_Dt_m);
	     for (int k=0; k<nlog_vars; k++) {
	       int jdof=log_vars[k];
	       double UU=exp(Ualpha.get(jdof));
	       dUdt.ir(1,jdof).scale(UU);
	     }
	     dUdt.rs();
	     
	     // Pass to the flux function the true positive values
	     U.prod(SHAPE,true_lstate,-1,-1,1);
	     grad_U.prod(dshapex,true_lstate,1,-1,-1,2);
	     
	     delta_sc=0;
	     double lambda_max_pg;
	     
	     // Compute A_grad_U in the `old' state
	     adv_diff_ff->set_state(Uo,grad_U); // fixme:= ojo que le pasamos
	     // grad_U (y no grad_Uold) ya que
	     // no nos interesa la parte difusiva
	     adv_diff_ff->compute_flux(Uo,iJaco,H,grad_H,flux,fluxd,
				       A_grad_U,grad_U,G_source,
				       tau_supg,delta_sc,
				       lambda_max_pg, nor,lambda,Vr,Vr_inv,
				       COMP_SOURCE | COMP_UPWIND);
	     adv_diff_ff->comp_A_grad_N(Ao_grad_N,dshapex);
#define USE_OLD_STATE_FOR_P_SUPG
#ifdef USE_OLD_STATE_FOR_P_SUPG
	     // This computes either the standard `P_supg' perturbation
	     // function or other written by the user in the
	     // flux-function.
	     adv_diff_ff->comp_P_supg(P_supg);
#endif
	     // Set the state of the fluid so that it can be used to
	     // compute matrix products
	     adv_diff_ff->set_state(U,grad_U);
	     adv_diff_ff->enthalpy_fun->set_state(U);
	     adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux,fluxd,
				       A_grad_U,grad_U,G_source,
				       tau_supg,delta_sc,
				       lambda_max_pg, nor,lambda,Vr,Vr_inv,
				       COMP_SOURCE | COMP_UPWIND);
	     
#ifdef COMPUTE_FD_ADV_JACOBIAN
	     double eps_fd=1e-4;
	     for (int jdof=1; jdof<=ndof; jdof++) {
	       U_pert.set(U).is(1,jdof).add(eps_fd).rs();
	       adv_diff_ff->set_state(U_pert,grad_U);
	       adv_diff_ff->compute_flux(U_pert,iJaco,H,grad_H,flux_pert,fluxd,
					 A_grad_U,grad_U,G_source,
					 tau_supg,delta_sc,
					 lambda_max_pg, nor,lambda,Vr,Vr_inv,0);
	       flux_pert.rest(flux).scale(1./eps_fd);
	       flux_pert.t();
	       A_fd_jac.ir(3,jdof).set(flux_pert).rs();
	       flux_pert.rs();
	     }
#endif
	     
	     if (lambda_max_pg>lambda_max) lambda_max=lambda_max_pg;
	     tmp10.set(G_source);	// tmp10 = G - dUdt
	     if (!lumped_mass) tmp10.rest(dUdt);
	     if (beta_supg==1.) {
	       tmp1.rs().set(tmp10).rest(A_grad_U); //tmp1= G - dUdt - A_grad_U
	     } else {
	       tmp1.set(dUdt).scale(-beta_supg).add(G_source);
	     }
	     // tmp15.set(SHAPE).scale(wpgdet*rec_Dt_m);
	     // tmp12.prod(SHAPE,tmp15,1,2); // tmp12 = SHAPE' * SHAPE
	     // tmp13.prod(tmp12,eye_ndof,1,3,2,4); // tmp13 = SHAPE' * SHAPE * I
	     adv_diff_ff->enthalpy_fun->comp_W_Cp_N(N_Cp_N,SHAPE,SHAPE,wpgdet*rec_Dt_m);
	     if (lumped_mass) {
	       matlocf_mass.add(N_Cp_N);
	     } else {
	       matlocf.add(N_Cp_N);
	     }
	     adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);
	     // Termino Galerkin
	     if (weak_form) {
	       assert(!lumped_mass && beta_supg==1.); // Not implemented yet!!
	       // weak version
	       tmp11.set(flux).rest(fluxd); // tmp11 = flux_c - flux_d
	       tmp23.set(SHAPE).scale(-wpgdet*ALPHA);
	       tmp14.prod(A_grad_N,tmp23,1,2,4,3);
	       matlocf.add(tmp14);
	     } else {
	       // tmp2.prod(SHAPE,tmp1,1,2); // tmp2= SHAPE' * (G - dUdt - A_grad_U)
	       tmp2.prod(SHAPE,A_grad_U,1,2); // tmp2= SHAPE' * A_grad_U
	       veccontr.axpy(tmp2,-wpgdet);
	       tmp11.set(0.).rest(fluxd); // tmp11 = - flux_d
	       tmp23.set(SHAPE).scale(wpgdet*ALPHA);
	       tmp14.prod(A_grad_N,tmp23,3,2,4,1);
	       matlocf.add(tmp14);
	     }
	     // tmp8= DSHAPEX * (w*flux_c - flux_d)
	     //            w = weak_form
	     tmp8.prod(dshapex,tmp11,-1,1,2,-1);
	     tmp9.prod(SHAPE,tmp10,1,2); // tmp9 = SHAPE' * (G - dUdt)
	     tmp8.add(tmp9);		// tmp8 = DSHAPEX * tmp11
	     veccontr.axpy(tmp8,wpgdet);
	     
	     // Diffusive term in matrix
	     adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
					       dshapex,wpgdet*ALPHA);
	     matlocf.add(grad_N_D_grad_N);
	     
	     // Reactive term in matrix (Galerkin part)
	     adv_diff_ff->comp_N_N_C(N_N_C,SHAPE,wpgdet*ALPHA);
	     matlocf.add(N_N_C);
	     
#ifndef USE_OLD_STATE_FOR_P_SUPG
	     // This computes either the standard `P_supg' perturbation
	     // function or other written by the user in the
	     // flux-function.
	     adv_diff_ff->comp_P_supg(P_supg);
#endif
	     // adding shock-capturing term
	     if (shocap>0.) {
	       tmp_sc.prod(dshapex,dshapex,-1,1,-1,2).scale(shocap*delta_sc*ALPHA*wpgdet);
	       for (int jdf=1; jdf<=ndof; jdf++) {
		 matlocf.ir(2,jdf).ir(4,jdf).add(tmp_sc).rs();
	       }
	       tmp_sc_v.prod(dshapex,grad_U,-1,1,-1,2).scale(-shocap*delta_sc*wpgdet);
	       veccontr.add(tmp_sc_v);
	       //	SHV(delta_sc);
	     }
	     for (int jel=1; jel<=nel; jel++) {
	       P_supg.ir(1,jel);
	       tmp4.prod(tmp1,P_supg,-1,1,-1);
	       //	  veccontr.ir(1,jel).axpy(tmp4,wpgdet).ir(1);
	       veccontr.ir(1,jel).axpy(tmp4,wpgdet);
	       
	       matlocf.ir(1,jel);
	       
	       tmp19.set(P_supg).scale(ALPHA*wpgdet);
	       tmp20.prod(tmp19,A_grad_N,1,-1,2,-1,3);
	       matlocf.add(tmp20);
	       
	       // Reactive term in matrix (SUPG term)
	       adv_diff_ff->comp_N_P_C(N_P_C,P_supg,SHAPE,wpgdet*ALPHA);
	       matlocf.add(N_P_C);
	       
	       tmp21.set(SHAPE).scale(beta_supg*wpgdet*rec_Dt_m);
	       adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
	       tmp22.prod(P_Cp,tmp21,1,3,2);
	       if (lumped_mass) {
		 // I think that if 'lumped_mass' is used then the
		 // contribution to the mass matrix from the SUPG
		 // perturbation term is null, but I include it.
		 matlocf_mass.ir(1,jel).add(tmp22).rs();
	       } else {
		 matlocf.add(tmp22);
	       }
	       matlocf.rs();
	       veccontr.rs();
	     }
	     P_supg.rs();
	   } else {
	     printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	     exit(1);
	   }
	   volume_flag=0;
	 }
	 if (comp_res) {
#if 0
	   // Compute the local critical time step. This is something
	   // like min (h/la_max) where h is the local mesh size, and
	   // la_max the local maximum eigenvalue. We do this for each
	   // Gauss point, and then take the minimum. The local h is
	   // estimated as twice the minimum length of the columns of
	   // Jaco. (The length of each column of Jaco is roughly half
	   // the corresponding h.) The maximum eigenvalue is estimated
	   // as the max over all the jacobian matrices of the norm1 od
	   // the matrix.
	   double hhh,hloc,dtloc;
	   Jaco_av.scale(1./double(npg));
	   hvec.sum_square(Jaco_av,-1,1);
	   hloc = 2.*sqrt(hvec.min_all());
	   
	   dtloc = hloc/lambda_max;
	   // PetscPrintf(PETSCFEM_COMM_WORLD,
	   // "On element %d, hloc %f, lambda_max %f, dtloc %f\n",
	   // k,hloc,lambda_max,dtloc);
	   
	   if (dtloc<DTMIN || !WAS_SET) {
	     DTMIN = dtloc;
	     WAS_SET = 1;
	     //   	PetscPrintf(PETSCFEM_COMM_WORLD,
	     //   		    "setting dtmin: %f, hloc %f, lambda_max: %f\n",
	     //   		    DTMIN,hloc,lambda_max);
	   }
	   
	   if (local_time_step_g) {
	     mass.scale(1./dtloc);
	     matloc.kron(mass,Id_ndof);
	   }
#endif
	   if (lumped_mass) {
	     // lump mass matrix
#if 1   // With this commented should be equivalent to no mass lumping
	     matlocf_mass.reshape(2,nen,nen);
	     for (int j=1; j<=nen; j++) {
	       double m = matlocf_mass.ir(1,j).sum_all();
	       matlocf_mass.set(0.).setel(m,j);
	     }
	     matlocf_mass.rs().reshape(4,nel,ndof,nel,ndof);
#endif
	     // Compute derivative of local element state
	     // The Dt is already included in the mass matrix
	     dUloc.set(lstaten).rest(lstateo);
	     dUloc_c.set(dUloc);
	     // correct the logarithmic variables for a factor $U^{n+\alpha}$
	     for (int k=0; k<nlog_vars; k++) {
	       int jdof=log_vars[k];
	       for (int j=1; j<nel; j++) {
		 true_lstate.ir(2,jdof);
		 dUloc_c.ir(2,jdof).mult(true_lstate);
	       }
	     }
	     dUloc_c.rs();
	     true_lstate.rs();
	     
	     // Compute inertia term with lumped mass
	     veccontr_mass.prod(matlocf_mass,dUloc_c,1,2,-1,-2,-1,-2);
	     // Add (rest) to vector contribution to be returned
	     veccontr.rest(veccontr_mass);
	     // Scale mass matrix by $U^{n+\alpha}/Dt*(1+alpha\DU)$
	     for (int k=0; k<nlog_vars; k++) {
	       int jdof=log_vars[k];
	       for (int j=1; j<=nel; j++) {
		 // Raw mass matrix value
		 double m = matlocf_mass.get(j,jdof,j,jdof);
		 // True (positive variable) value
		 double UU = true_lstate.get(j,jdof);
		 // Difference of log variable
		 double DU = dUloc.get(j,jdof);
		 // Correction factor
		 double f=UU*(1+ALPHA*DU);
		 matlocf_mass.setel(m*f,j,jdof,j,jdof);
		 // Correct the spatial jacobian
		 matlocf.ir(3,j).ir(4,jdof).scale(UU);
	       }
	     }
	     matlocf.rs();
	     // Add to matrix contribution to be returned
	     matlocf.add(matlocf_mass);
	   }
	   veccontr.export_vals(element.ret_vector_values(*retval));
#ifdef CHECK_JAC
	   veccontr.export_vals(element.ret_fdj_values(*fdj_jac));
#endif
	   if (comp_mat_each_time_step_g) {
	     matlocf.add(matlocf_fix);
	     matlocf.export_vals(element.ret_mat_values(*Ajac));
	   }
	 } else if (comp_prof) {
	   matlocf.export_vals(element.ret_mat_values(*jac_prof));
	 }
	 
       } catch (GenericError e) {
	 set_error(1);
	 return;
       }
  
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
} catch (GenericError e) {
  set_error(1);
}

void AdvDiff_Godunov::init() {
  get_prop(normal_prop,"normal");
  int ret_options=0;
  NSGETOPTDEF(int,ncells,0);
  assert(ncells);
  NSGETOPTDEF(double,fs_tol,1.0e-3);
}

void AdvDiff_Godunov::element_hook(ElementIterator &element){
  element_m = element;
  int pp=adv_diff_ff->dim();
  assert(normal_prop.length == pp);
  normal.set(prop_array(element_m,normal_prop));
}
