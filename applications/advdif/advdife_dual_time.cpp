//__INSERT_LICENSE__
//$Id: advdife.cpp,v 1.105 2005/05/26 22:07:26 mstorti Exp $
extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
extern int MY_RANK,SIZE;

#include <vector>
#include <string>

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>
#include <src/generror.h>

#include "nwadvdif.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void NewAdvDifFF::get_preco(FastMat2 &preco) {
  PetscPrintf(PETSC_COMM_WORLD,
	      "Using dual-time formulation needs definition for get_preco() virtual function\n"
	      "in the flux function object.\n");
  assert(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "NewAdvDif::assemble"
void NewAdvDif
::new_assemble_preco(arg_data_list &arg_data_v,const Nodedata *nodedata,
                     const Dofmap *dofmap,const char *jobinfo,
                     const ElementList &elemlist,
                     const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  int locdof,kldof,lldof;

  int BOOL_COMP_RES=1;
  int BOOL_COMP_JAC=1;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  assert(npg>0);
  assert(ndim>0);

  int nelprops,nel,ndof;
  elem_params(nel,ndof,nelprops);
  int PFUNUSED nen = nel*ndof;

  // Unpack Dofmap
  int nnod;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  // H is a generalized local property passed per node with the nodal
  // coordinates. In shallow water nH =1 and H is the depth. It is
  // needed in order to compute the source term. In 1D Euler it may be
  // the area section of the tube. Its gradient is needed for the
  // source term in the momentum eqs.
  int nH = nu-ndim;
  // FMatrix  Hloc(nel,nH),H(nH),vloc_mesh(nel,ndim),v_mesh(ndim);
  FMatrix  Hloc(nel,nH),H(nH),vloc_mesh(nel,ndim);
  //  FastMat2 v_mesh;

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  double *retvalt;
  time_m = double(* (const Time *) time_data);

  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step.
  vector<double> *dtmin;
  double lambda_max=NAN;
  int PFUNUSED jdtmin;
  GlobParam *glob_param=NULL;
  // The trapezoidal rule integration parameter
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)
#define DPT (glob_param->Dpt)
  arg_data *staten=NULL,*stateo=NULL,*statepo=NULL,*retval=NULL,
    *fdj_jac=NULL,*jac_prof=NULL,*Ajac=NULL;
  if (comp_res) {
    int j=-1;
    stateo  = &arg_data_v[++j]; //[0]
    statepo = &arg_data_v[++j]; //[1]
    staten  = &arg_data_v[++j]; //[2]
    retval  = &arg_data_v[++j]; //[3]
    jdtmin = ++j;//[4]
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];//[5]
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
    rec_Dt_m = 1./DT;
    if (glob_param->steady) rec_Dt_m = 0.;

// apply ALPHA to time step in order to take into account higher order temporal integration
    rec_Dt_m = rec_Dt_m/ALPHA;

    rec_Dpt = 1./DPT;
    if (glob_param->psteady) rec_Dpt = 0.;

    if (ADVDIF_CHECK_JAC)
      fdj_jac = &arg_data_v[++j];
  }

  FastMat2 matlocf(4,nel,ndof,nel,ndof),
    matlocf_mass(4,nel,ndof,nel,ndof);
  FastMat2 prof_nodes(2,nel,nel), prof_fields(2,ndof,ndof),
    matlocf_fix(4,nel,ndof,nel,ndof);
  FastMat2 Id_ndf(2,ndof,ndof),Id_nel(2,nel,nel),
    prof_fields_diag_fixed(2,ndof,ndof);

  //o Use the weak form for the Galerkin part of the advective term.
  NSGETOPTDEF(int,weak_form,1);
  // o Weights the temporal term with $N+\beta P$, i.e.
  // $\beta=0$ is equivalent to weight the temporal term a la
  // Galerkin and $\beta=1$ is equivalent to do the consistent SUPG weighting.
  //  NSGET OPTDEF(double,beta_supg,1.);
  //o Use lumped mass (used mainly to avoid oscillations for small time steps).
  NSGETOPTDEF(int,lumped_mass,0);
  //o Add a shock capturing term
  NSGETOPTDEF(double,shocap,0.0);
  //o Add an anisotropic shock capturing term
  NSGETOPTDEF(double,shocap_aniso,0.0);
  //o Use the advective Jacobian in the previous time step
  //  for the SUPG stabilization term. This accelerates
  //  convergence of the Newton iteration.
  NSGETOPTDEF(int,use_Ajac_old,0);
  //o Report jacobians on random elements (should be in range 0-1).
  NSGETOPTDEF(double,compute_fd_adv_jacobian_random,1.0);
  //o ALE_flag : flag to ON ALE computation
  NSGETOPTDEF(int,ALE_flag,0);
  //o indx_ALE_xold : pointer to old coordinates in
  //  NODEDATA array excluding the first "ndim" values
  NSGETOPTDEF(int,indx_ALE_xold,1);
  //o Compute the non-symmetric term
  //  correspoding to nonlinearities in the
  //  diffusive Jacobian.
  NSGETOPTDEF(int,compute_dDdU_term,1);

  //o Use 1 Gauss point for some matrix terms.
  // In this application for diffusive terms grad_N_D_grad_N
  // and for stabilized advection terms P_supg_A_grad_N.
  NSGETOPTDEF(int,use_low_gpdata,0);

  //o key for computing reactive terms or not
  NSGETOPTDEF(int,compute_reactive_terms,1);

  assert(compute_fd_adv_jacobian_random>0.
	 && compute_fd_adv_jacobian_random <=1.);

  //o Add axisymmetric version for this particular elemset.
  NSGETOPTDEF(string,axisymmetric,"none");
  assert(axisymmetric.length()>0);
  if (axisymmetric=="none") axi=0;
  else if (axisymmetric=="x") axi=1;
  else if (axisymmetric=="y") axi=2;
  else if (axisymmetric=="z") axi=3;
  else {
    PetscPrintf(PETSC_COMM_WORLD,
		"Invalid value for \"axisymmetric\" option\n"
		"axisymmetric=\"%s\"\n",axisymmetric.c_str());
    PetscFinalize();
    exit(0);
  }

  //o Uses operations caches for computations with the FastMat2
  //  library for internal computations. Note that this affects also the
  //  use of caches in routines like fluxes, etc...
  NSGETOPTDEF(int,use_fastmat2_cache,1);

  // Initialize flux functions
  ff_options=0;
  adv_diff_ff->start_chunk(ff_options);
  int ndimel = adv_diff_ff->dim();
  if (ndimel<0) ndimel = ndim;
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

#if 0
  if (use_log_vars) {
    // Bes sure that we are in shallow water.
    // This would be returned by the flux function
    if (ndof==5) {   // turbulent shallow water
      nlog_vars=2;
      log_vars = log_vars_swt; // Return dofs for k, epsilon
    } else if (ndof==1) { // thermal problem
      nlog_vars=1;
      log_vars = &log_vars_sc;
    } else { assert(0);}
  } else {
    nlog_vars=0;
    log_vars=NULL;
  }
#endif

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
    dUloc(nel,ndof),matloc,
    lstatepo(nel,ndof);
  FastMat2 true_lstate(2,nel,ndof),
    true_lstateo(2,nel,ndof),true_lstaten(2,nel,ndof),
    true_lstatepo(2,nel,ndof);

  FastMat2 true_lstate_abs(2,nel,ndof);

  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);
  GPdata gp_data_low(geometry.c_str(),ndimel,nel,1,GP_FASTMAT2);

  double detJaco, wpgdet, delta_sc, delta_sc_old;
  int elem, ipg,node, jdim, kloc,lloc,ldof;
  double lambda_max_pg;

  dshapex.resize(2,ndimel,nel);
  FMatrix Jaco(ndimel,ndim),Jaco_av(ndimel,ndim),
    iJaco(ndimel,ndimel),
    flux(ndof,ndimel),fluxd(ndof,ndimel),mass(nel,nel),
    grad_U(ndimel,ndof), A_grad_U(ndof),Ao_grad_U(ndof),
    G_source(ndof), dUdt(ndof), Un(ndof),
    Ho(ndof),Hn(ndof);

  // Agregado para computar los terminos pseudo-temporales
  FMatrix dUdpt(ndof),Upo(ndof),tmp100,tmp210,tmp220;
  FastMat2 preco(2,ndof,ndof),N_Ga_N(4,nel,ndof,nel,ndof),
    P_Ga(2,ndof,ndof);

  // These are declared but not used
  FMatrix nor,lambda,Vr,Vr_inv,U(ndof),Ualpha(ndof),
    lmass(nel),Id_ndof(ndof,ndof),
    tmp1,tmp2,tmp3,tmp4,tmp5,hvec(ndimel),tmp6,tmp7,
    tmp8,tmp9,tmp10,tmp11(ndof,ndimel),tmp12,tmp14,
    tmp15,tmp17,tmp19,tmp20,tmp21,tmp22,tmp23,tmp1_old,
    tmp24,tmp_sc,tmp_sc_v,tmp_shc_grad_U,
    tmp_j_grad_U(ndof),tmp_j_gradN,
    tmp_sc_aniso,tmp_matloc_aniso,
    tmp_sc_v_aniso;
  FMatrix tmp_ALE_01,tmp_ALE_02,
    tmp_ALE_03,tmp_ALE_04,tmp_ALE_05,
   tmp_ALE_06,tmp_ALE_07;
  FastMat2 A_grad_N(3,nel,ndof,ndof),
    grad_N_D_grad_N(4,nel,ndof,nel,ndof),N_N_C(4,nel,ndof,nel,ndof),
    N_P_C(3,ndof,nel,ndof),N_Cp_N(4,nel,ndof,nel,ndof),
    P_Cp(2,ndof,ndof),grad_N_dDdU_N(4,nel,ndof,nel,ndof);
  Ao_grad_N.resize(3,nel,ndof,ndof);
  tau_supg.resize(2,ndof,ndof);
  P_supg.resize(3,nel,ndof,ndof);
  grad_U_norm.resize(2,ndimel,ndof);

  FMatrix Jaco_axi(2,2);
  int ind_axi_1, ind_axi_2;
  double detJaco_axi;

  FastMat2 Cr(2,ndof,ndof);
  FastMat2 Cp_bis(2,ndof,ndof),
    Cp_bis_old(2,ndof,ndof),Ao(3,ndim,ndof,ndof);
  FastMat2 delta_sc_v(1,ndof);

  FMatrix Jaco_low(ndimel,ndim),iJaco_low(ndimel,ndimel);
  double detJaco_low=NAN,wpgdet_low;
  // dshapex_low.resize(2,ndimel,nel);

  if (axi) assert(ndim==3);

  // For the computation of the jacobian with
  // finite differences
  FastMat2 A_fd_jac(3,ndimel,ndof,ndof),U_pert(1,ndof),
    flux_pert(2,ndof,ndimel),A_jac_err, A_jac(3,ndimel,ndof,ndof),
    Id_ndim(2,ndim,ndim),jvec(1,ndim),jvec_old(1,ndim);
  Id_ndim.eye();

  double delta_aniso,delta_aniso_old;

  // Position of current element in elemset and in chunk
  int k_elem, k_chunk;

  v_mesh.resize(1,ndim);

  Uo.resize(1,ndof);
  Id_ndof.set(0.);
  for (int j=1; j<=ndof; j++) Id_ndof.setel(1.,j,j);

  FastMatCacheList cache_list;
  if (use_fastmat2_cache) FastMat2::activate_cache(&cache_list);

  // printf("[%d] %s start: %d last: %d\n",MY_RANK,jobinfo,el_start,el_last);
  for (ElementIterator element = elemlist.begin();
       element!=elemlist.end(); element++) try {

    element.position(k_elem,k_chunk);
    FastMat2::reset_cache();

    // Initialize element
    adv_diff_ff->element_hook(element);
    // Get nodedata info (coords. etc...)
    element.node_data(nodedata,xloc.storage_begin(),
		      Hloc.storage_begin());

    if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
      continue;
    }

    if (comp_res) {
      lambda_max=0;
      lstateo.set(element.vector_values(*stateo));
      lstaten.set(element.vector_values(*staten));
      lstatepo.set(element.vector_values(*statepo));
    }

    // State at time t_{n+\alpha}
    lstate.set(0.).axpy(lstaten,ALPHA).axpy(lstateo,(1-ALPHA));
    log_transf(true_lstate ,lstate ,nlog_vars,log_vars);
    log_transf(true_lstateo,lstateo,nlog_vars,log_vars);
    log_transf(true_lstatepo,lstatepo,nlog_vars,log_vars);

    veccontr.set(0.);
    mass.set(0.);
    lmass.set(0.);
    matlocf.set(0.);
    if (lumped_mass) matlocf_mass.set(0.);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE	 (*gp_data.FM2_shape[ipg])
#define WPG	 (gp_data.wpg[ipg])

#define DSHAPEXI_LOW (*gp_data_low.FM2_dshapexi[0])
#define SHAPE_LOW    (*gp_data_low.FM2_shape[0])
#define WPG_LOW	     (gp_data_low.wpg[0])

#undef SHV
#define SHV(mess,v)                                     \
     printf("[%d] elem %d, %s " #v " %f\n",             \
	    MY_RANK,k_elem,mess,v.sum_square_all(),     \
	    v.sum_square_all());
    
    // assert(ndim==ndimel);
    // vaux.set(DSHAPEXI_LOW);
    Jaco_low.prod(DSHAPEXI_LOW,xloc,1,-1,-1,2);
    if (ndim==ndimel) {
      detJaco_low = Jaco_low.det();
      iJaco_low.inv(Jaco_low);
    } else if (ndimel==1) {
      detJaco_low = Jaco_low.detsur();
      iJaco_low.setel(1./detJaco_low,1,1);
    }
    wpgdet_low = detJaco_low*WPG_LOW;
    dshapex.prod(iJaco_low,DSHAPEXI_LOW,1,-1,-1,2);
    
    if (axi >0){
      ind_axi_1 = (	axi   % 3)+1;
      ind_axi_2 = ((axi+1) % 3)+1;
      Jaco_axi.setel(Jaco_low.get(ind_axi_1,ind_axi_1),1,1);
      Jaco_axi.setel(Jaco_low.get(ind_axi_1,ind_axi_2),1,2);
      Jaco_axi.setel(Jaco_low.get(ind_axi_2,ind_axi_1),2,1);
      Jaco_axi.setel(Jaco_low.get(ind_axi_2,ind_axi_2),2,2);
      detJaco_axi = Jaco_axi.det();
      double wpgdet_axi = detJaco_axi*WPG_LOW;
      Volume = 0.5*fabs(wpgdet_axi);
      
    } else {
      Volume = wpgdet_low;
    }
    volume_flag = 1;

    // nodal computation of mesh velocity
    if (ALE_flag) {
      assert(nH >= ndim);
      assert(indx_ALE_xold <= nH+1-ndim);
      Hloc.is(2,indx_ALE_xold,indx_ALE_xold+ndim-1);
      vloc_mesh.set(xloc).minus(Hloc).scale(rec_Dt_m).rs();
      Hloc.rs();
    }
    
    if(use_low_gpdata) {
      
      if (comp_res) {
	
	// with old state
	Uo.prod(SHAPE_LOW,true_lstateo,-1,-1,1);
	adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
	grad_Uo.prod(dshapex,true_lstateo,1,-1,-1,2);

	true_lstate_abs.set(true_lstateo.fun(abs));
	grad_U_norm.prod(dshapex,true_lstate_abs,1,-1,-1,2);

	adv_diff_ff->set_state(Uo,grad_Uo);
	adv_diff_ff->compute_flux(Uo,iJaco_low,H,grad_H,flux,fluxd,
				  A_grad_U,grad_Uo,G_source,
				  tau_supg,delta_sc_old,
				  lambda_max_pg, nor,lambda,Vr,Vr_inv,
				  COMP_SOURCE | COMP_UPWIND);
	if (BOOL_COMP_JAC)
	  adv_diff_ff->comp_A_grad_N(Ao_grad_N,dshapex);
	
	adv_diff_ff->comp_P_supg(P_supg);

	if (use_Ajac_old) adv_diff_ff->get_Ajac(Ao);

	// with pseudo-old state
	Upo.prod(SHAPE_LOW,true_lstatepo,-1,-1,1);

	// with current state
	U.prod(SHAPE_LOW,true_lstate,-1,-1,1);
	adv_diff_ff->enthalpy_fun->enthalpy(Hn,U);
	grad_U.prod(dshapex,true_lstate,1,-1,-1,2);

	adv_diff_ff->set_state(U,grad_U);
	adv_diff_ff->compute_flux(U,iJaco_low,H,grad_H,flux,fluxd,
				  A_grad_U,grad_U,G_source,
				  tau_supg,delta_sc_old,
				  lambda_max_pg, nor,lambda,Vr,Vr_inv,
				  COMP_SOURCE | COMP_UPWIND);
	if (BOOL_COMP_JAC)
	  adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);
	
	adv_diff_ff->get_preco(preco);

	if (ALE_flag) {
	  v_mesh.prod(SHAPE_LOW,vloc_mesh,-1,-1,1);
	  adv_diff_ff->get_Cp(Cp_bis);
	  if (BOOL_COMP_JAC){
	    tmp_ALE_01.prod(v_mesh,dshapex,-1,-1,1);
	  }
	  tmp_ALE_02.prod(v_mesh,grad_U,-1,-1,1);
	}

	tmp11.set(0.).minus(fluxd);
	tmp8.prod(dshapex,tmp11,-1,1,2,-1);
	veccontr.axpy(tmp8,wpgdet_low);

	if (BOOL_COMP_JAC){
	  adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
					    dshapex,wpgdet_low);
	  matlocf.add(grad_N_D_grad_N);
	}
	
	// Psuedo-temporal derivate
	dUdpt.set(Un).minus(Upo).scale(rec_Dpt);
	tmp100.prod(preco,dUdpt,1,-1,-1);

	dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);
	
	tmp10.set(G_source);	// tmp10 = G - dUdt
	if (!lumped_mass) tmp10.minus(dUdt);
	tmp10.minus(tmp100);	// tmp10 = G - dUdt - Gamma*dUdtau
	
	tmp1.rs().set(tmp10).minus(A_grad_U); //tmp1= G - dUdt - Gamma*dUdtau - A_grad_U
	
	if (use_Ajac_old) {
	  Ao_grad_U.prod(Ao,grad_U,-1,1,-2,-1,-2);
	  tmp1_old.rs().set(tmp10).minus(Ao_grad_U); //tmp1= G - dUdt - Gamma*dUdtau - A_grad_U
	}
	
	for (int jel=1; jel<=nel; jel++) {
	  P_supg.ir(1,jel);
	  matlocf.ir(1,jel);
	  
	  veccontr.ir(1,jel);
	  
	  if (use_Ajac_old) {
	    tmp4.prod(tmp1_old,P_supg,-1,1,-1);
	  } else {
	    tmp4.prod(tmp1,P_supg,-1,1,-1);
	  }
	  
	  veccontr.axpy(tmp4,wpgdet_low);
	  
	  if (BOOL_COMP_JAC){
	    tmp19.set(P_supg).scale(wpgdet_low);
	    if (use_Ajac_old)
	      tmp20.prod(tmp19,Ao_grad_N,1,-1,2,-1,3);
	    else tmp20.prod(tmp19,A_grad_N,1,-1,2,-1,3);
	    matlocf.add(tmp20);
	  }
	  
	  if(!lumped_mass) {

	    adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
	   
	    if (BOOL_COMP_JAC){
	      if(compute_reactive_terms){
		// Reactive term in matrix (SUPG term)
		adv_diff_ff->comp_N_P_C(N_P_C,P_supg,SHAPE_LOW,wpgdet_low);
		matlocf.add(N_P_C);
	      }
	    
	      tmp21.set(SHAPE_LOW).scale(wpgdet_low*rec_Dt_m);
	      tmp22.prod(P_Cp,tmp21,1,3,2);
	      matlocf.add(tmp22);
	    
	      tmp210.set(SHAPE_LOW).scale(wpgdet_low*rec_Dpt);
	      adv_diff_ff->enthalpy_fun->comp_P_Gamma(P_Ga,P_supg);
	      tmp220.prod(P_Ga,tmp210,1,3,2);
	      matlocf.add(tmp220);
	    }

	    if (ALE_flag) {
	      if (BOOL_COMP_JAC){
		tmp_ALE_07.prod(P_Cp,tmp_ALE_01,1,3,2);
		matlocf.axpy(tmp_ALE_07,-wpgdet_low);
	      }
	      tmp_ALE_06.prod(P_Cp,tmp_ALE_02,1,-1,-1);
	      veccontr.axpy(tmp_ALE_06,wpgdet_low);
	    }
	    
	  }
	}
	
	matlocf.rs();
	P_supg.rs();
	veccontr.rs();
	
      }	      
    }

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

      dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);

      if (nH>0) {
	H.prod(SHAPE,Hloc,-1,-1,1);
	grad_H.prod(dshapex,Hloc,1,-1,-1,2);
      }

      if (comp_res) {

	// state variables and gradient
	Un.prod(SHAPE,lstate,-1,-1,1);
	adv_diff_ff->enthalpy_fun->enthalpy(Hn,Un);
	Uo.prod(SHAPE,lstateo,-1,-1,1);
	adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
	Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
	dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);

	Upo.prod(SHAPE,lstatepo,-1,-1,1);
	dUdpt.set(Un).minus(Upo).scale(rec_Dpt);

	for (int k=0; k<nlog_vars; k++) {
	  int jdof=log_vars[k];
	  double UU=exp(Ualpha.get(jdof));
	  dUdt.ir(1,jdof).scale(UU);
	}
	dUdt.rs();

	// Pass to the flux function the true positive values
	U.prod(SHAPE,true_lstate,-1,-1,1);
	grad_U.prod(dshapex,true_lstate,1,-1,-1,2);
	grad_Uo.prod(dshapex,true_lstateo,1,-1,-1,2);

	delta_sc=0;
	delta_sc_old=0;

	true_lstate_abs.set(true_lstateo.fun(abs));
	grad_U_norm.prod(dshapex,true_lstate_abs,1,-1,-1,2);

	// Compute A_grad_U in the `old' state
	adv_diff_ff->set_state(Uo,grad_Uo); // fixme:= ojo que le pasamos
					   // grad_U (y no grad_Uold) ya que
					   // no nos interesa la parte difusiva
	adv_diff_ff->compute_flux(Uo,iJaco,H,grad_H,flux,fluxd,
				  A_grad_U,grad_Uo,G_source,
				  tau_supg,delta_sc_old,
				  lambda_max_pg, nor,lambda,Vr,Vr_inv,
				  COMP_SOURCE | COMP_UPWIND);
	if (BOOL_COMP_JAC)	
	  adv_diff_ff->comp_A_grad_N(Ao_grad_N,dshapex);
	Ao_grad_U.set(A_grad_U);
	if (shocap>0. || shocap_aniso>0.)
	  adv_diff_ff->get_Cp(Cp_bis_old);
	if (use_Ajac_old) adv_diff_ff->get_Ajac(Ao);
	  adv_diff_ff
	    ->compute_shock_cap_aniso(delta_aniso_old,jvec_old);

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

	if (!lumped_mass) {
	  adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux,fluxd,
				    A_grad_U,grad_U,G_source,
				    tau_supg,delta_sc,
				    lambda_max_pg, nor,lambda,Vr,Vr_inv,
				    COMP_SOURCE | COMP_UPWIND);
	} else {
	  adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux,fluxd,
				    A_grad_U,grad_U,G_source,
				    tau_supg,delta_sc,
				    lambda_max_pg, nor,lambda,Vr,Vr_inv,
				    COMP_SOURCE_NOLUMPED | COMP_UPWIND);
	}

	if (compute_fd_adv_jacobian) {
	  int check_this = drand()<compute_fd_adv_jacobian_random;
	  comp_total++;
	  FastMat2::branch();
	  if (check_this) {
	    FastMat2::choose(0);
	    comp_checked++;
	    double &eps_fd = compute_fd_adv_jacobian_eps;
	    for (int jdof=1; jdof<=ndof; jdof++) {
	      U_pert.set(U).is(1,jdof).add(eps_fd).rs();
	      adv_diff_ff->set_state(U_pert,grad_U);
	      adv_diff_ff->compute_flux(U_pert,iJaco,H,grad_H,flux_pert,fluxd,
					A_grad_U,grad_U,G_source,
					tau_supg,delta_sc,
					lambda_max_pg, nor,lambda,Vr,Vr_inv,0);
	      flux_pert.minus(flux).scale(1./eps_fd);
	      flux_pert.t();
	      A_fd_jac.ir(3,jdof).set(flux_pert).rs();
	      flux_pert.rs();
	    }
	    for (int j=1; j<=ndim; j++) {
	      Id_ndim.ir(2,j);
	      A_jac.ir(1,j);
	      adv_diff_ff->comp_A_jac_n(A_jac,Id_ndim);
	    }
	    Id_ndim.rs();
	    A_jac.rs();
	    A_jac_err.set(A_jac).minus(A_fd_jac);
#define FM2_NORM sum_abs_all
	    double A_jac_norm = A_jac.FM2_NORM();
	    double A_jac_err_norm = A_jac_err.FM2_NORM();
	    double A_fd_jac_norm = A_fd_jac.FM2_NORM();
	    double A_rel_err = A_jac_err_norm/A_fd_jac_norm;

#define CHK_MAX(v) if(v>v##_max) v##_max=v
#define CHK_MIN(v) if(v<v##_min) v##_min=v
#define CHK(v) CHK_MAX(v); CHK_MIN(v)
	    CHK(A_jac_norm);
	    CHK(A_jac_err_norm);
	    CHK(A_fd_jac_norm);
	    CHK(A_rel_err);
#undef CHK_MAX
#undef CHK_MIN
#undef CHK

	    int print_this =
	      A_rel_err >= compute_fd_adv_jacobian_rel_err_threshold;
	    if (compute_fd_adv_jacobian>=2 && print_this) {
	      printf("elem %d, |A_a|=%g, |A_n|=%g, |A_a-A_n|=%g, (rel.err %g)\n",
		     k_elem,A_jac_norm,A_fd_jac_norm,A_jac_err_norm,
		     A_rel_err);
	    }
	    if (compute_fd_adv_jacobian>=3  && print_this) {
	      A_jac.print("A_a: ");
	      A_fd_jac.print("A_n: ");
	    }
	    // Reset state in flux function to state U
	    adv_diff_ff->set_state(U,grad_U);
	    adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux_pert,fluxd,
				      A_grad_U,grad_U,G_source,
				      tau_supg,delta_sc,
				      lambda_max_pg, nor,lambda,Vr,Vr_inv,0);
	  }
	  FastMat2::leave();
	}

	if (lambda_max_pg>lambda_max) lambda_max=lambda_max_pg;
	
	adv_diff_ff->get_preco(preco);
	tmp100.prod(preco,dUdpt,1,-1,-1);

	tmp10.set(G_source);	// tmp10 = G - dUdt
	if (!lumped_mass) tmp10.minus(dUdt);
	tmp10.minus(tmp100);	// tmp10 = G - dUdt - Gamma*dUdtau

	tmp1.rs().set(tmp10).minus(A_grad_U); //tmp1= G - dUdt - Gamma*dUdtau - A_grad_U

	if (use_Ajac_old) {
	  Ao_grad_U.prod(Ao,grad_U,-1,1,-2,-1,-2);
	  tmp1_old.rs().set(tmp10).minus(Ao_grad_U); //tmp1= G - dUdt - Gamma*dUdtau - A_grad_U
	}
	
	if (BOOL_COMP_JAC){
	  if (!lumped_mass) {
	    adv_diff_ff->enthalpy_fun
	      ->comp_W_Cp_N(N_Cp_N,SHAPE,SHAPE,
			    wpgdet*rec_Dt_m);
	    matlocf.add(N_Cp_N);
	  }
	  adv_diff_ff->enthalpy_fun->comp_W_Gamma_N(N_Ga_N,SHAPE,SHAPE,
						    wpgdet*rec_Dpt);
	  matlocf.add(N_Ga_N);	

	  // A_grad_N.prod(dshapex,A_jac,-1,1,-1,2,3);
	  adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);
	}

	// add ALE Galerkin terms
	if (ALE_flag) {
	  v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	  adv_diff_ff->get_Cp(Cp_bis);
	  tmp_ALE_02.prod(v_mesh,grad_U,-1,-1,1);
	  tmp_ALE_03.prod(SHAPE,Cp_bis,1,2,3);

	  tmp_ALE_04.prod(tmp_ALE_03,tmp_ALE_02,1,2,-1,-1);
	  veccontr.axpy(tmp_ALE_04,wpgdet);

	  if (BOOL_COMP_JAC){
	    tmp_ALE_01.prod(v_mesh,dshapex,-1,-1,1);
	    tmp_ALE_05.prod(tmp_ALE_03,tmp_ALE_01,1,2,4,3);
	    matlocf.axpy(tmp_ALE_05,-wpgdet);
	  }

	}

	// Termino Galerkin
	if (weak_form) {
	  if (use_low_gpdata){
	    tmp11.set(flux); // tmp11 = flux_c (viscous part is integrated in 1 PG)
	  } else {
	    tmp11.set(flux).minus(fluxd); // tmp11 = flux_c - flux_d
	  }

	  if (BOOL_COMP_JAC){
	    tmp23.set(SHAPE).scale(-wpgdet);
	    tmp14.prod(A_grad_N,tmp23,1,2,4,3);
	    matlocf.add(tmp14);
	  }
	} else {
	  // tmp2.prod(SHAPE,tmp1,1,2); // tmp2= SHAPE' * (G - dUdt - A_grad_U)
	  tmp2.prod(SHAPE,A_grad_U,1,2); // tmp2= SHAPE' * A_grad_U
	  veccontr.axpy(tmp2,-wpgdet);

	  //	  tmp11.set(0.).minus(fluxd); // tmp11 = - flux_d
	  if (use_low_gpdata){
	    tmp11.set(0.); // tmp11 = 0 (viscous part is integrated in 1 PG)
	  } else {
	    tmp11.set(0.).minus(fluxd); // tmp11 = - flux_d
	  }

	  if (BOOL_COMP_JAC){
	    tmp23.set(SHAPE).scale(wpgdet);
	    tmp14.prod(A_grad_N,tmp23,3,2,4,1);
	    matlocf.add(tmp14);
	  }
	}
	// tmp8= DSHAPEX * (w*flux_c - flux_d)
	//	      w = weak_form
	tmp8.prod(dshapex,tmp11,-1,1,2,-1);
	tmp9.prod(SHAPE,tmp10,1,2); // tmp9 = SHAPE' * (G - dUdt - Gamma*dUdtau)
	tmp8.add(tmp9);		// tmp8 = DSHAPEX * tmp11
	veccontr.axpy(tmp8,wpgdet);

	// Diffusive term in matrix
#if 0
	tmp17.set(dshapex).scale(wpgdet);
	adv_diff_ff->comp_D_grad_N(D_grad_N,tmp17);
	tmp18.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
#endif
	if (BOOL_COMP_JAC){
	  if(use_low_gpdata==0){
	    adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
					      dshapex,wpgdet);
	    matlocf.add(grad_N_D_grad_N);
	  }
	}

	// add non linear contributions of
	// diffusive matrix, Djac(U) ==> d Djac /dU
	if (BOOL_COMP_JAC){
	  if (compute_dDdU_term) {
	    int flag = adv_diff_ff
	      ->comp_grad_N_dDdU_N(grad_N_dDdU_N,grad_U,
				   dshapex,SHAPE,wpgdet);
	    if (flag) matlocf.add(grad_N_dDdU_N);
	  }
	}

	// Reactive term in matrix (Galerkin part)
#if 0
	tmp23.set(SHAPE).scale(wpgdet);
	tmp24.prod(SHAPE,tmp23,1,2); // tmp24 = SHAPE' * SHAPE
	tmp25.prod(tmp24,C_jac,1,3,2,4); // tmp25 = SHAPE' * SHAPE * C_jac
#endif

	if (BOOL_COMP_JAC){
	  if(compute_reactive_terms){
	    if (!lumped_mass) {
	      adv_diff_ff->comp_N_N_C(N_N_C,SHAPE,wpgdet);
	      matlocf.add(N_N_C);
	    }
	  }
	}

	// adding shock-capturing term
	delta_sc_v.set(0.0);
	if (shocap>0. ) {
	  adv_diff_ff->compute_delta_sc_v(delta_sc_v);
	  
	  tmp_shc_grad_U.prod(Cp_bis_old,grad_U,2,-1,1,-1);
	  for (int jdf=1; jdf<=ndof; jdf++) {
	    delta_sc_v.addel(delta_sc_old,jdf);
	  }

	  if (BOOL_COMP_JAC)
	    tmp_sc.prod(dshapex,dshapex,-1,1,-1,2).scale(shocap*wpgdet);
	  tmp_sc_v.prod(dshapex,tmp_shc_grad_U,-1,1,-1,2);
	  for (int jdf=1; jdf<=ndof; jdf++) {
	    double delta = (double)delta_sc_v.get(jdf);
	    if (BOOL_COMP_JAC){
	      for (int kdf=1; kdf<=ndof; kdf++) {
		double tmp_shc_1=Cp_bis_old.get(jdf,kdf);
		matlocf.ir(2,jdf).ir(4,kdf).axpy(tmp_sc,delta*tmp_shc_1).rs();
	      }
	    }
	    tmp_sc_v.ir(2,jdf)
	      .scale(-shocap*delta*wpgdet).rs();
	    delta_sc_v.rs();
	  }

	  veccontr.add(tmp_sc_v);
	}

	// adding ANISOTROPIC shock-capturing term
	// Falta usar el Cp_bis_old aca tambien
	if (shocap_aniso>0.) {
	  delta_aniso = 0.0;
	  adv_diff_ff
	    ->compute_shock_cap_aniso(delta_aniso,jvec);

	  tmp_shc_grad_U.prod(Cp_bis_old,grad_U,2,-1,1,-1);
	  tmp_j_grad_U.prod(jvec_old,tmp_shc_grad_U,-1,-1,1);
	  tmp_j_gradN.prod(jvec_old,dshapex,-1,-1,1);

	  tmp_sc_v_aniso.prod(tmp_j_gradN,tmp_j_grad_U,1,2);

	  if (BOOL_COMP_JAC){
	    tmp_sc_aniso.prod(tmp_j_gradN,tmp_j_gradN,1,2)
	      .scale(shocap_aniso*wpgdet);
	    tmp_matloc_aniso.prod(tmp_sc_aniso,Cp_bis_old,1,3,2,4);
	    matlocf.axpy(tmp_matloc_aniso,delta_aniso_old);
	  }
	  veccontr.axpy(tmp_sc_v_aniso,
			-shocap_aniso*delta_aniso_old*wpgdet);
	}

#ifndef USE_OLD_STATE_FOR_P_SUPG
	// This computes either the standard `P_supg' perturbation
	// function or other written by the user in the
	// flux-function.
	adv_diff_ff->comp_P_supg(P_supg);
#endif
	  
	  
	for (int jel=1; jel<=nel; jel++) {
	  P_supg.ir(1,jel);
	  
	  veccontr.ir(1,jel);
	  matlocf.ir(1,jel);
	  
	  if (use_Ajac_old) {
	    tmp4.prod(tmp1_old,P_supg,-1,1,-1);
	  } else {
	    tmp4.prod(tmp1,P_supg,-1,1,-1);
	  }
	  
	  veccontr.axpy(tmp4,wpgdet);
	  
	  if (BOOL_COMP_JAC){
	    tmp19.set(P_supg).scale(wpgdet);
	    if (use_Ajac_old)
	      tmp20.prod(tmp19,Ao_grad_N,1,-1,2,-1,3);
	    else tmp20.prod(tmp19,A_grad_N,1,-1,2,-1,3);
	    matlocf.add(tmp20);
	  }
	  
	  if(!lumped_mass) {

	    adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
	   
	    if (BOOL_COMP_JAC){   
	      if(compute_reactive_terms){
		// Reactive term in matrix (SUPG term)
		adv_diff_ff->comp_N_P_C(N_P_C,P_supg,SHAPE,wpgdet);
		matlocf.add(N_P_C);
	      }
	      
	      tmp21.set(SHAPE).scale(wpgdet*rec_Dt_m);
	      tmp22.prod(P_Cp,tmp21,1,3,2);
	      matlocf.add(tmp22);
	    
	      tmp210.set(SHAPE).scale(wpgdet*rec_Dpt);
	      adv_diff_ff->enthalpy_fun->comp_P_Gamma(P_Ga,P_supg);
	      tmp220.prod(P_Ga,tmp210,1,3,2);
	      matlocf.add(tmp220);
	    }

	    if (ALE_flag) {
	      if (BOOL_COMP_JAC){
		tmp_ALE_07.prod(P_Cp,tmp_ALE_01,1,3,2);
		matlocf.axpy(tmp_ALE_07,-wpgdet);
	      }
	      tmp_ALE_06.prod(P_Cp,tmp_ALE_02,1,-1,-1);
	      veccontr.axpy(tmp_ALE_06,wpgdet);
	    }
	  }
	}
	matlocf.rs();
	veccontr.rs();
	P_supg.rs();
	
      } else {
	
	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);
	
      }
      
      lmass.axpy(SHAPE,wpgdet);
      
    }  // ipg loop
      
    volume_flag=0;

    if (comp_res) {
      
      if (lumped_mass) {
	  
	// lumped terms treatment
	// temporal and source terms are lumped out of the gauss points loop
	for (int j=1; j<=nel; j++) {
	    
	  lstate.ir(1,j);
	  Un.set(lstate);
	  lstate.rs();
	    
	  lstateo.ir(1,j);
	  Uo.set(lstateo);
	  lstateo.rs();
	  //	adv_diff_ff->set_state(Uo,grad_U);
	  adv_diff_ff->set_state(Uo);
	    
	  // adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
	  adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
	    
	  adv_diff_ff->set_state(Un,grad_U);
	  adv_diff_ff->compute_flux(Un,iJaco,H,grad_H,flux,fluxd,
				    A_grad_U,grad_U,G_source,
				    tau_supg,delta_sc,
				    lambda_max_pg, nor,lambda,Vr,Vr_inv,
				    COMP_SOURCE_LUMPED);
	  adv_diff_ff->enthalpy_fun->enthalpy(Hn,Un);
	    
	  Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
	    
	  //	adv_diff_ff->enthalpy_fun->comp_P_Cp(Cp,Id_ndf);
	  adv_diff_ff->get_Cp(Cp_bis);
	  Cp_bis.scale(rec_Dt_m);
	    
	  dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);
	    
	  for (int k=0; k<nlog_vars; k++) {
	    int jdof=log_vars[k];
	    double UU=exp(Ualpha.get(jdof));
	    dUdt.ir(1,jdof).scale(UU);
	  }
	  dUdt.rs();
	    
	  tmp10.set(G_source).minus(dUdt);	// tmp10 = G - dUdt
	    
	  if (BOOL_COMP_JAC){
	    adv_diff_ff->get_C(Cr);
	    tmp12.set(Cr).add(Cp_bis);
	  }

	  double mass_node;
	  mass_node = double(lmass.get(j));
	  veccontr.ir(1,j).axpy(tmp10,mass_node).rs();
	  if (BOOL_COMP_JAC)
	    matlocf.ir(1,j).ir(3,j).axpy(tmp12,mass_node).rs();

	}
      }

      veccontr.export_vals(element.ret_vector_values(*retval));
      if (ADVDIF_CHECK_JAC)
        veccontr.export_vals(element.ret_fdj_values(*fdj_jac));
      if (BOOL_COMP_JAC){
	if (comp_mat_each_time_step_g) {
	  matlocf.add(matlocf_fix);
	  matlocf.export_vals(element.ret_mat_values(*Ajac));
	}
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
} 
