//__INSERT_LICENSE__
extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;

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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "NewAdvDif::assemble"
void NewAdvDif::
new_assemble_ALE_formulation(arg_data_list &arg_data_v,const Nodedata *nodedata,
			     const Dofmap *dofmap,const char *jobinfo,
			     const ElementList &elemlist,
			     const TimeData *time_data) try {
  static int use_ALE_formulation_reported = 0;
  if (!use_ALE_formulation_reported) {
    if (!MY_RANK) 
      printf("=======================================\n"
             "=======================================\n"
             "====   USING ALE+GCL Formulation   ====\n"
             "=======================================\n"
             "=======================================\n");
    use_ALE_formulation_reported = 1;
  }
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr = 0;
  int locdof, kldof, lldof;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  assert(npg > 0);
  assert(ndim > 0);

  int nelprops,nel,ndof;
  elem_params(nel, ndof, nelprops);
  int nen = nel*ndof;

  // Unpack Dofmap
  int neq, nnod;
  neq  = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu = nodedata->nu;
  // H is a generalized local property passed per node with the nodal
  // coordinates. In shallow water nH =1 and H is the depth. It is
  // needed in order to compute the source term. In 1D Euler it may be
  // the area section of the tube. Its gradient is needed for the
  // source term in the momentum eqs.
  int nH = nu-ndim;
  //  FMatrix  Hloc(nel,nH),H(nH),vloc_mesh(nel,ndim),v_mesh(ndim);
  FMatrix Hloc(nel,nH),H(nH),vloc_mesh(nel,ndim);
  //  FastMat2 v_mesh;

  if(nnod != nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  double *retvalt;
  time_m = double(* (const Time *) time_data);

  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step.
  vector<double> *dtmin;
  double         lambda_max  = NAN;
  int            jdtmin;
  GlobParam      *glob_param = NULL;

  // The trapezoidal rule integration parameter
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)

  arg_data *staten = NULL, *stateo = NULL, *retval = NULL,
    *fdj_jac = NULL, *jac_prof = NULL, *Ajac = NULL;

  if (comp_res) {
    int j  = -1;
    stateo = &arg_data_v[++j]; //[0]
    staten = &arg_data_v[++j]; //[1]
    retval = &arg_data_v[++j]; //[2]
    jdtmin = ++j;              //[3]

#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set

    Ajac       = &arg_data_v[++j]; //[4]
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
    rec_Dt_m   = 1./DT;
    if (glob_param->steady) rec_Dt_m = 0.;
    assert(ALPHA > 0.); // FIX THIS IN ORDER TO USE ALPHA=0 !!!!

#ifdef CHECK_JAC
    fdj_jac = &arg_data_v[++j];
#endif
  }

  FastMat2 matlocf(4,nel,ndof,nel,ndof),
    matlocf_mass(4,nel,ndof,nel,ndof);
  FastMat2 prof_nodes(2,nel,nel), prof_fields(2,ndof,ndof),
    matlocf_fix(4,nel,ndof,nel,ndof);
  FastMat2 Id_ndof(2,ndof,ndof),Id_nel(2,nel,nel),
    prof_fields_diag_fixed(2,ndof,ndof);

  //o Use the weak form for the Galerkin part of the advective term.
  NSGETOPTDEF(int,weak_form,1);
  // o Weights the temporal term with $N+\beta P$, i.e.
  // $\beta=0$ is equivalent to weight the temporal term a la
  // Galerkin and $\beta=1$ is equivalent to do the consistent SUPG weighting.
  //  NSGET OPTDEF(double,beta_supg,1.);
  //o Use lumped mass (used mainly to avoid oscillations for small time steps).
  NSGETOPTDEF(int,lumped_mass,0);
  //o Add a shock capturing term (delta_sc = shocap * delta_sc_ff + shocap_const)
  NSGETOPTDEF(double,shocap,0.0);
  //o Add a constant shock capturing term
  NSGETOPTDEF(double,shocap_const,0.0);
  //o Add an anisotropic shock capturing term 
  //o (delta_sc = shocap * delta_sc_ff + shocap_const)
  NSGETOPTDEF(double,shocap_aniso,0.0); 
  //o Add a constant anisotropic shock capturing term
  //o (delta_sc_aniso = shocap_aniso * delta_sc_aniso_ff + shocap_aniso_const)
  NSGETOPTDEF(double,shocap_aniso_const,0.0);
  //o Use the advective Jacobian in the previous time step
  //  for the SUPG stabilization term. This accelerates
  //  convergence of the Newton iteration.
  NSGETOPTDEF(int,use_Ajac_old,0);
  //o Report jacobians on random elements (should be in range 0-1).
  NSGETOPTDEF(double,compute_fd_adv_jacobian_random,1.0);
  //o Pointer to old coordinates in
  //  #nodedata# array excluding the first "ndim" values
  NSGETOPTDEF(int,indx_ALE_xold,1);
  //o Compute the term non-symmetric term
  //  correspoding to nonlinearities in the
  //  diffusive Jacobian.
  NSGETOPTDEF(int,compute_dDdU_term,1);
  //o Don't use the averaged Jacobian ALE fix for the DGCL. This is only
  //  for debugging purposes, you should use the fix!! 
  NSGETOPTDEF(int,dont_use_average_jaco_fix,0);

  //o key for computing reactive terms or not
  NSGETOPTDEF(int,compute_reactive_terms,1);

  static int ale_with_no_weak_form_warning_given = 0;
  if (!ale_with_no_weak_form_warning_given
      && !weak_form && use_ALE_form) {
    printf("=============================================\n"
           "=============================================\n"
           "WARNING: WEAK_FORM=0 AND ALE IS A BAD CHOICE \n"
           "=============================================\n"
           "=============================================\n");
    ale_with_no_weak_form_warning_given = 1;
  }
  assert(compute_fd_adv_jacobian_random>0.
	 && compute_fd_adv_jacobian_random <=1.);

  //o Add axisymmetric version for this particular elemset.
  NSGETOPTDEF(string,axisymmetric,"none");
  assert(axisymmetric.length() > 0);
  if (axisymmetric == "none")   axi = 0;
  else if (axisymmetric == "x") axi = 1;
  else if (axisymmetric == "y") axi = 2;
  else if (axisymmetric == "z") axi = 3;
  else {
    PetscPrintf(PETSCFEM_COMM_WORLD,
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
  ff_options = 0;
  adv_diff_ff->start_chunk(ff_options);
  int ndimel = adv_diff_ff->dim();
  if (ndimel < 0) ndimel = ndim;
  FMatrix grad_H(ndimel,nH);

  // Don't know how to handle the following combinations
  assert(!(use_ALE_form && ndim!=ndimel));
  assert(!(use_ALE_form && lumped_mass));

  adv_diff_ff->set_profile(prof_fields); // profile by equations
  prof_nodes.set(1.);

  Id_ndof.eye();
  Id_nel.eye();

  prof_fields.d(1,2);
  prof_fields_diag_fixed.set(0.);
  prof_fields_diag_fixed.d(1,2);
  prof_fields_diag_fixed.set(prof_fields);
  prof_fields.rs();
  prof_fields_diag_fixed.rs();
  prof_fields_diag_fixed.scale(-1.).add(Id_ndof);

  matlocf_fix.prod(prof_fields_diag_fixed,Id_nel,2,4,1,3);
  matlocf.prod(prof_fields,prof_nodes,2,4,1,3);

  matlocf.add(matlocf_fix);
  if (comp_res)
    matlocf.export_vals(Ajac->profile);
  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matlocf.export_vals(jac_prof->profile);
  }

  // lumped_mass:= If this options is activated then all the inertia
  // term matrix comtributions are added to 'matlocf_mass' and the
  // vector contribution terms are discarded. Then at the last moment
  // matlocf_mass*(Un-Uo)/Dt is added.

  // Allocate local vecs
  FMatrix veccontr(nel,ndof), veccontr_mass(nel,ndof),
    xloc(nel,ndim), xloc_old(nel,ndim), xloc_new(nel,ndim), lstate(nel,ndof),
    lstateo(nel,ndof), lstaten(nel,ndof), dUloc_c(nel,ndof),
    dUloc(nel,ndof), matloc, xloc_mid(nel,ndim);

  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);

  double detJaco, detJaco_new = NAN, detJaco_old = NAN, detJaco_mid = NAN,
    wpgdet, delta_sc, delta_sc_old, lambda_max_pg;
  int elem, ipg,node, jdim, kloc, lloc, ldof;

  dshapex.resize(2,ndimel,nel);
  dshapex_gcl.resize(2,ndimel,nel);
  FMatrix Jaco(ndimel,ndim), Jaco_av(ndimel,ndim),
    Jaco_new(ndimel,ndim), Jaco_old(ndimel,ndim),
    iJaco(ndimel,ndimel), iJaco_old(ndimel,ndimel),
    Q(ndimel,ndimel), Jaco_mid(ndimel,ndim),iJaco_mid(ndimel,ndim),
    iJaco_new(ndimel,ndimel), flux(ndof,ndimel),
    fluxd(ndof,ndimel), mass(nel,nel),
    grad_U(ndimel,ndof), A_grad_U(ndof),
    Ao_grad_U(ndof), G_source(ndof), 
    dHdt(ndof), dHdt2(ndof), Un(ndof),
    Ho(ndof),Hn(ndof),Halpha(ndof);
  // When ALE: dHdt is in fact d(JH)/dt, and dHdt2 is d(H)/dt
  // In the residual 
  // FMatrix grad_U_norm(ndimel,ndof);
  // These are declared but not used
  FMatrix nor, lambda, Vr, Vr_inv, U(ndof), Ualpha(ndof),
    lmass(nel), tmp1, tmp2, tmp3, tmp4, tmp5, hvec(ndimel),
    tmp6, tmp7, tmp8,tmp9, tmp10, tmp10j, tmp11(ndof,ndimel),
    tmp12, tmp14, tmp14b, tmp15, tmp17, tmp19, tmp20, tmp21, 
    tmp22, tmp23, tmp1_old, tmp24,tmp_sc, tmp_sc_v,tmp_shc_grad_U,
    tmp_j_grad_U(ndof), tmp_j_gradN, tmp_sc_aniso, tmp_matloc_aniso,
    tmp_sc_v_aniso, tmp_ALE_flux, tmp_ALE_jac, v_mesh_grad_N(nel),
    tmp_ALE_01,tmp_ALE_02, tmp_ALE_03, tmp_ALE_04, tmp_ALE_05,
    tmp_ALE_06,tmp_ALE_07;
  FastMat2 A_grad_N(3,nel,ndof,ndof),
    grad_N_D_grad_N(4,nel,ndof,nel,ndof),N_N_C(4,nel,ndof,nel,ndof),
    N_P_C(3,ndof,nel,ndof),N_Cp_N(4,nel,ndof,nel,ndof),
    P_Cp(2,ndof,ndof),grad_N_dDdU_N(4,nel,ndof,nel,ndof);

  Ao_grad_N.resize(3,nel,ndof,ndof);
  tau_supg.resize(2,ndof,ndof);
  P_supg.resize(3,nel,ndof,ndof);
  Cp.resize(2,ndof,ndof);
  Cp_old.resize(2,ndof,ndof);
  Uo.resize(1,ndof);

  FMatrix Jaco_axi(2,2);
  int ind_axi_1, ind_axi_2;
  double detJaco_axi;

  FastMat2 Cr(2,ndof,ndof), Ao(3,ndim,ndof,ndof), delta_sc_v(1,ndof);

  if (axi) assert(ndim == 3);

  // For the computation of the jacobian with
  // finite differences
  FastMat2 A_fd_jac(3,ndimel,ndof,ndof), U_pert(1,ndof),
    flux_pert(2,ndof,ndimel), A_jac_err, A_jac(3,ndimel,ndof,ndof),
    Id_ndim(2,ndim,ndim), jvec(1,ndim), jvec_old(1,ndim);
  double delta_aniso, delta_aniso_old;
  int k_elem, k_chunk;

  Id_ndim.eye();

  v_mesh.resize(1,ndim);
  PETSCFEM_ASSERT(nH >= ndim,"This element requires the old mesh "
                  "position to be passed as an H field. nH %d, ndim %d",
                  nH,ndim);
  PETSCFEM_ASSERT(indx_ALE_xold >= nH+1-ndim,
                  "bad indx_ALE_xold, not remaining enough columns. "
                  "indx_ALE_xold %d, nH %d, ndim %d",
                  indx_ALE_xold,nH,ndim);  

  FastMatCacheList cache_list;
  if (use_fastmat2_cache) FastMat2::activate_cache(&cache_list);

  // printf("[%d] %s start: %d last: %d\n",MY_RANK,jobinfo,el_start,el_last);
  for (ElementIterator element = elemlist.begin();
       element != elemlist.end(); element++) try {

      element.position(k_elem,k_chunk);
      FastMat2::reset_cache();

      // Initialize element
      adv_diff_ff->element_hook(element);
      // Get nodedata info (coords. etc...)
      element.node_data(nodedata,xloc.storage_begin(),
			Hloc.storage_begin());
      xloc_new.set(xloc);

      if (comp_prof) {
	matlocf.export_vals(element.ret_mat_values(*jac_prof));
	continue;
      }

      if (comp_res) {
	lambda_max = 0;
	lstateo.set(element.vector_values(*stateo));
	lstaten.set(element.vector_values(*staten));
      }

      // State at time t_{n+\alpha}
      lstate.set(0.).axpy(lstaten,ALPHA).axpy(lstateo,(1-ALPHA));

      veccontr.set(0.);
      mass.set(0.);
      lmass.set(0.);
      matlocf.set(0.);

      v_mesh.set(0.);

      if (lumped_mass) matlocf_mass.set(0.);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE	 (*gp_data.FM2_shape[ipg])
#define WPG	 (gp_data.wpg[ipg])

#undef SHV
#define SHV(mess,v)                                     \
      printf("[%d] elem %d, %s " #v " %f\n",		\
	     MY_RANK,k_elem,mess,v.sum_square_all(),	\
	     v.sum_square_all());
    
      volume_flag = 1;

      if (0){
	// DEBUG
	int kk,ielhh;
	element.position(kk,ielhh);
	printf("Element %d \n",kk);
	lstate.print("Estado :");
	// END DEBUG
      }
    
      // nodal computation of mesh velocity
      Hloc.is(2, indx_ALE_xold, indx_ALE_xold+ndim-1);
      xloc_old.set(Hloc);
      Hloc.rs();
      xloc.scale(ALPHA).axpy(xloc_old,1-ALPHA);
      vloc_mesh.set(xloc_new).rest(xloc_old).scale(rec_Dt_m).rs();
      if (0){
	int kk,ielhh;
	element.position(kk, ielhh);
	if (rand() % 50 == 0){
	  printf("Element %d \n",kk);
	  xloc.print(" xloc^(n+1) :");
	  Hloc.print(" Hloc^(n+1) :");
	  vloc_mesh.print(" vloc_mesh^(n+1) :");
	}
      }

      // loop over Gauss points
      Jaco_av.set(0.);
      for (ipg = 0; ipg < npg; ipg++) {
	//      Matrix xpg = SHAPE * xloc;
	Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2); // xloc is at t_{n+alpha}
	Jaco_av.add(Jaco);
	Jaco_new.prod(DSHAPEXI,xloc_new,1,-1,-1,2);
	Jaco_old.prod(DSHAPEXI,xloc_old,1,-1,-1,2);
      
	if (ndim == ndimel) {
	  iJaco.inv(Jaco);
	  detJaco     = Jaco.det();
	  detJaco_new = Jaco_new.det();
	  detJaco_old = Jaco_old.det();
	  //         printf("detjaco (start,new,old): %g %g %g\n",detJaco, 
	  //                detJaco_new, detJaco_old);
	  
	  if (ndim == 2){ // we need only two points in 2D to integ temporal average
	    iJaco_new.inv(Jaco_new);
	    Q.set(0.0).axpy(iJaco_new,detJaco_new*0.5);
	    iJaco_old.inv(Jaco_old);
	    Q.axpy(iJaco_old,detJaco_old*0.5);
	    Q.scale(1.0/detJaco);
	    dshapex_gcl.prod(Q,DSHAPEXI,1,-1,-1,2);
	  } else if (ndim == 3) { // we need 3 points in 3D to integ
				  // temporal average with
				  // Gauss-Lobatto
	    iJaco_new.inv(Jaco_new);
	    Q.set(0.0).axpy(iJaco_new,detJaco_new/6.0);
	    iJaco_old.inv(Jaco_old);
	    Q.axpy(iJaco_old,detJaco_old/6.0);
	    
	    // for 3D holy grial GCL mid point integ
	    xloc_mid.set(xloc_new).scale(0.5).axpy(xloc_old,0.5).rs(); 
	    Jaco_mid.prod(DSHAPEXI,xloc_mid,1,-1,-1,2);
	    detJaco_mid = Jaco_mid.det();
	    iJaco_mid.inv(Jaco_mid);
	    Q.axpy(iJaco_mid,detJaco_mid*2.0/3.0);
	    Q.scale(1.0/detJaco);
	    dshapex_gcl.prod(Q,DSHAPEXI,1,-1,-1,2);
	  }
          if (dont_use_average_jaco_fix) 
            dshapex_gcl.prod(iJaco,DSHAPEXI,1,-1,-1,2);
	} else if (ndimel == 1) {
	  // This allows to solve problems on streams like rivers or
	  // ducts or advective problems on plane surfaces (not
	  // implemented yet). Also, it could be used also for advective
	  // problems on arbitrary surfaces (ndim=3 and ndimel=2) but I
	  // don't know how to do that yet. (tensorial calculus...)
	  detJaco = Jaco.norm_p_all(2);
	  iJaco.setel(1./detJaco,1,1);
	}

	if (detJaco <= 0.) {
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
	  Un.prod(SHAPE,lstaten,-1,-1,1);//lstaten is at t_{n+1}
	  adv_diff_ff->enthalpy_fun->enthalpy(Hn,Un);
	  Uo.prod(SHAPE,lstateo,-1,-1,1);
	  adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
	  U.prod(SHAPE,lstate,-1,-1,1); // state at t_{n+alpha}
	  adv_diff_ff->enthalpy_fun->enthalpy(Halpha,U);
	  // This is d(JH)/dt (J is the jacobian of the
	  // ALE transformation
	  // Hn represents H in time t^{n+1} so that it must be affected
	  // by detJaco_new
	  dHdt.set(Hn).scale(detJaco_new/detJaco).axpy(Ho,-detJaco_old/detJaco)
	    .scale(rec_Dt_m);
	  // This is plain dH/dt
	  dHdt2.set(Hn).rest(Ho).scale(rec_Dt_m);

	  dHdt.rs();

	  // Pass to the flux function the true positive values
	  grad_U.prod(dshapex,lstate,1,-1,-1,2);// grad U at t_{n+alpha}
	  grad_Uo.prod(dshapex,lstateo,1,-1,-1,2);

	  delta_sc     = 0;
	  delta_sc_old = 0;

	  // Compute A_grad_U in the `old' state
	  adv_diff_ff->set_state(Uo,grad_Uo); // fixme:= ojo que le pasamos
	  // grad_U (y no grad_Uold) ya que no nos interesa la parte difusiva

	  v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	
	  adv_diff_ff->compute_flux(Uo,iJaco,H,grad_H,flux,fluxd,
				    A_grad_U,grad_Uo,G_source,
				    tau_supg,delta_sc_old,
				    lambda_max_pg, nor,lambda,Vr,Vr_inv,
				    COMP_SOURCE | COMP_UPWIND);
	  adv_diff_ff->comp_A_grad_N(Ao_grad_N,dshapex);
	  Ao_grad_U.set(A_grad_U);
	  
	  if (shocap > 0. || shocap_aniso > 0. || 
	      shocap_const > 0. || shocap_aniso_const > 0.)
	    adv_diff_ff->get_Cp(Cp_old);
	  
	  if (use_Ajac_old) adv_diff_ff->get_Ajac(Ao);
	  jvec_old.set(0.);
	  if (shocap_aniso > 0. || shocap_aniso_const > 0.){
	    jvec_old.set(NAN);
	    adv_diff_ff->compute_shock_cap_aniso(delta_aniso_old,jvec_old);
	    if (isnan(jvec_old.get(1))) {
	      if (shocap_aniso_const > 0.0) {
		jvec_old.set(0.0);
	      } else {
		jvec_old.set(0.0);
#if 0
		printf("Anisotropic Shock-capturing term NOT implemented for"
		       " this flux function\n");
		throw GenericError("not implemented shocap_aniso error");
#endif
	      }
	      
	    }
	  }
	  // Set the state of the fluid so that it can be used to
	  // compute matrix products
	  adv_diff_ff->set_state(U,grad_U); // U at t_{n+alpha}
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
	      for (int jdof = 1; jdof <= ndof; jdof++) {
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
	      for (int j = 1; j <= ndim; j++) {
		Id_ndim.ir(2,j);
		A_jac.ir(1,j);
		adv_diff_ff->comp_A_jac_n(A_jac,Id_ndim);
	      }
	      Id_ndim.rs();
	      A_jac.rs();
	      A_jac_err.set(A_jac).rest(A_fd_jac);
#define FM2_NORM sum_abs_all
	      double A_jac_norm     = A_jac.FM2_NORM();
	      double A_jac_err_norm = A_jac_err.FM2_NORM();
	      double A_fd_jac_norm  = A_fd_jac.FM2_NORM();
	      double A_rel_err      = A_jac_err_norm/A_fd_jac_norm;

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
	      if (compute_fd_adv_jacobian >= 2 && print_this) {
		printf("elem %d, |A_a|=%g, |A_n|=%g, |A_a-A_n|=%g, (rel.err %g)\n",
		       k_elem,A_jac_norm,A_fd_jac_norm,A_jac_err_norm,
		       A_rel_err);
	      }
	      if (compute_fd_adv_jacobian >= 3  && print_this) {
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

	  if (lambda_max_pg > lambda_max) lambda_max = lambda_max_pg;

	  tmp10.set(G_source);	// tmp10 = G - dHdt
	  if (!lumped_mass) tmp10.rest(dHdt);
	  tmp10j.set(G_source);	// tmp10 = G - dHdt
	  if (!lumped_mass) tmp10j.rest(dHdt2);
	
	  // For the stabilization term use dHdt2
	  // (i.e. not including dJ/dt term)
	  // tmp1= G - dHdt2 - A_grad_U
	  tmp1.rs().set(tmp10j).rest(A_grad_U); 
	
	  if (use_Ajac_old) {
	    Ao_grad_U.prod(Ao,grad_U,-1,1,-2,-1,-2);
	    //tmp1= G - dHdt - A_grad_U
	    tmp1_old.rs().set(tmp10).rest(Ao_grad_U);
	  }

	  adv_diff_ff->set_state(Un);//para obtener el Cp^{n+1}
	  adv_diff_ff->get_Cp(Cp);//para obtener el Cp^{n+1}
	  if (!lumped_mass) {
	    adv_diff_ff->enthalpy_fun
	      ->comp_W_Cp_N(N_Cp_N,SHAPE,SHAPE,
			    detJaco_new*WPG*rec_Dt_m/ALPHA);
	    matlocf.add(N_Cp_N);
	  }
	  adv_diff_ff->set_state(U); //vuelvo a t_{n+alpha}
	  adv_diff_ff->get_Cp(Cp); //vuelvo a t_{n+alpha}

	  // A_grad_N.prod(dshapex,A_jac,-1,1,-1,2,3);
	  adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);

	  // add ALE Galerkin terms
	  v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	  //	adv_diff_ff->get_Cp(Cp); //Cp in t_{n+alpha}
	  tmp_ALE_01.prod(v_mesh,dshapex,-1,-1,1);
	  tmp_ALE_02.prod(v_mesh,grad_U,-1,-1,1);
	
	  if (!weak_form) {
	    tmp_ALE_03.prod(SHAPE,Cp,1,2,3);
	    tmp_ALE_04.prod(tmp_ALE_03,tmp_ALE_02,1,2,-1,-1);
	    veccontr.axpy(tmp_ALE_04,wpgdet);
	    tmp_ALE_05.prod(tmp_ALE_03,tmp_ALE_01,1,2,4,3);
	    matlocf.axpy(tmp_ALE_05,-wpgdet);
	    // tmp2.prod(SHAPE,tmp1,1,2); // tmp2= SHAPE' * (G - dHdt - A_grad_U)
	    tmp2.prod(SHAPE,A_grad_U,1,2); // tmp2= SHAPE' * A_grad_U
	    veccontr.axpy(tmp2,-wpgdet);
	  
	    tmp11.set(0.).rest(fluxd); // tmp11 = - flux_d
	  
	    tmp23.set(SHAPE).scale(wpgdet);
	    tmp14.prod(A_grad_N,tmp23,3,2,4,1);
	    matlocf.add(tmp14);
	  } else {
	    // tmp11 = flux_c - flux_d
	    tmp11.set(flux).rest(fluxd); 
	    // Add ALE correction to flux
	    // tmp11 = flux_c - v_mesh*H^{n+alpha} - flux_d
	    tmp_ALE_flux.prod(Halpha,v_mesh,1,2);
	    tmp11.rest(tmp_ALE_flux);
	  
	    tmp23.set(SHAPE).scale(-wpgdet);
	    v_mesh_grad_N.prod(v_mesh,dshapex,-1,-1,1);
	    tmp_ALE_jac.prod(v_mesh_grad_N,Cp,1,2,3);// el Cp tiene que estar en t_{n+alpha}
	    tmp14b.set(A_grad_N).rest(tmp_ALE_jac);
	    tmp14.prod(tmp14b,tmp23,1,2,4,3);
	    matlocf.add(tmp14);
	  }	  
	  // tmp8= DSHAPEX * (w*flux_c - flux_d - v_mesh*H)
	  // w = weak_form
	  // tmp8.prod(dshapex,tmp11,-1,1,2,-1); // non-averaged shape function gradients
	  tmp8.prod(dshapex_gcl,tmp11,-1,1,2,-1);
	  tmp9.prod(SHAPE,tmp10,1,2); // tmp9 = SHAPE' * (G - dHdt)
#if 0
	  FMSHV(tmp8);
	  FMSHV(tmp9);
	  FMSHV(veccontr);
#endif
	  tmp8.add(tmp9);		// tmp8 = DSHAPEX * tmp11
	  veccontr.axpy(tmp8,wpgdet);
	  // printf("dHdt part: %f\n",veccontr.sum_all());
#if 0
	  FMSHV(veccontr);
#endif
	
	  // Diffusive term in matrix
#if 0
	  tmp17.set(dshapex).scale(wpgdet);
	  adv_diff_ff->comp_D_grad_N(D_grad_N,tmp17);
	  tmp18.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
#endif
	
	  adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
					    dshapex,wpgdet);
	  matlocf.add(grad_N_D_grad_N);

	  // add non linear contributions of
	  // diffusive matrix, Djac(U) ==> d Djac /dU
	  if (compute_dDdU_term) {
	    int flag = adv_diff_ff
	      ->comp_grad_N_dDdU_N(grad_N_dDdU_N,grad_U,
				   dshapex,SHAPE,wpgdet);
	    if (flag) matlocf.add(grad_N_dDdU_N);
	  }

	  // Reactive term in matrix (Galerkin part)
#if 0
	  tmp23.set(SHAPE).scale(wpgdet);
	  tmp24.prod(SHAPE,tmp23,1,2); // tmp24 = SHAPE' * SHAPE
	  tmp25.prod(tmp24,C_jac,1,3,2,4); // tmp25 = SHAPE' * SHAPE * C_jac
#endif

	  if(compute_reactive_terms){
	    if (!lumped_mass) {
	      adv_diff_ff->comp_N_N_C(N_N_C,SHAPE,wpgdet);
	      matlocf.add(N_N_C);
	    }
	  }

	  // adding shock-capturing term
	  delta_sc_v.set(0.0);
	  if (shocap > 0. || shocap_const > 0.) {
	    delta_sc_v.set(NAN);
	    adv_diff_ff->compute_delta_sc_v(delta_sc_v);
	    if (isnan(delta_sc_v.get(1))) {
	      if (shocap_const > 0.0) {
		delta_sc_v.set(0.0);
	      } else {
		delta_sc_v.set(0.0);
#if 0
		printf("Shock-capturing term NOT implemented for"
		       " this flux function\n");
		throw GenericError("not implemented shocap error");
#endif
	      }
	    }
	    tmp_shc_grad_U.prod(Cp_old,grad_U,2,-1,1,-1);
	    for (int jdf = 1; jdf <= ndof; jdf++) {
	      delta_sc_v.addel(delta_sc_old,jdf);
	    }
	    // tmp_sc.prod(dshapex,dshapex,-1,1,-1,2).scale(shocap*wpgdet);
	    tmp_sc.prod(dshapex,dshapex,-1,1,-1,2).scale(wpgdet);
	    // tmp_sc_v.prod(dshapex,grad_U,-1,1,-1,2);
	    tmp_sc_v.prod(dshapex,tmp_shc_grad_U,-1,1,-1,2);
	    for (int jdf = 1; jdf <= ndof; jdf++) {
	      double delta        = (double)delta_sc_v.get(jdf);
	      double delta_sc_eff = shocap * delta + shocap_const;
	      for (int kdf = 1; kdf <= ndof; kdf++) {
		//	      double tmp_shc_1=Cp.get(jdf,kdf);
		double tmp_shc_1 = Cp_old.get(jdf,kdf);
		matlocf.ir(2,jdf).ir(4,kdf).axpy(tmp_sc,delta_sc_eff*tmp_shc_1).rs();
		//       matlocf.ir(2,jdf).ir(4,jdf).axpy(tmp_sc,delta).rs();
	      }
	      tmp_sc_v.ir(2,jdf)
		.scale(-delta_sc_eff*wpgdet).rs();
	      delta_sc_v.rs();
	    }
	    veccontr.add(tmp_sc_v);
	  }

	  // adding ANISOTROPIC shock-capturing term
	  jvec.set(0.);
	  if (shocap_aniso > 0. || shocap_aniso_const > 0.) {
	    delta_aniso = 0.0;
	    jvec.set(NAN);
	    adv_diff_ff->compute_shock_cap_aniso(delta_aniso,jvec);
	    if (isnan(jvec.get(1))) {
	      if (shocap_aniso_const > 0.0) {
		jvec.set(0.0);
	      } else {
		jvec.set(0.0);
#if 0
		printf("Anisotropic Shock-capturing term NOT implemented for"
		       " this flux function\n");
		throw GenericError("not implemented shocap_aniso error");
#endif
	      }
	    }
	    
	    tmp_shc_grad_U.prod(Cp_old,grad_U,2,-1,1,-1);
	    tmp_j_grad_U.prod(jvec_old,tmp_shc_grad_U,-1,-1,1);
	    tmp_j_gradN.prod(jvec_old,dshapex,-1,-1,1);
	    
	    for (int jdf = 1; jdf <= ndof; jdf++) {
	      jvec.addel(delta_aniso_old,jdf);
	    }
	    
	    tmp_sc_aniso.prod(tmp_j_gradN,tmp_j_gradN,1,2)
	      .scale(wpgdet);
	    
	    tmp_sc_v_aniso.prod(tmp_j_gradN,tmp_j_grad_U,1,2);
	    
	    for (int jdf = 1; jdf <= ndof; jdf++) {
	      double delta_sc_aniso     = (double)jvec.get(jdf);
	      double delta_sc_aniso_eff = shocap_aniso * delta_sc_aniso + 
		shocap_aniso_const;
	      for (int kdf = 1; kdf <= ndof; kdf++) {
		double tmp_shc_aniso_1 = Cp_old.get(jdf,kdf);
		matlocf.ir(2,jdf).ir(4,kdf)
		  .axpy(tmp_sc_aniso,delta_sc_aniso_eff*tmp_shc_aniso_1).rs();
	      }
	      tmp_sc_v_aniso.ir(2,jdf)
		.scale(-delta_sc_aniso_eff*wpgdet).rs();
	      jvec.rs();
	    }
	    veccontr.add(tmp_sc_v_aniso);
	  }
	  
#ifndef USE_OLD_STATE_FOR_P_SUPG
	  // This computes either the standard `P_supg' perturbation
	  // function or other written by the user in the
	  // flux-function.
	  adv_diff_ff->comp_P_supg(P_supg);
#endif
	  
	  for (int jel = 1; jel <= nel; jel++) {
	    P_supg.ir(1,jel);
	    veccontr.ir(1,jel);
	    matlocf.ir(1,jel);
	    if (use_Ajac_old) 
	      tmp4.prod(tmp1_old,P_supg,-1,1,-1);
	    else tmp4.prod(tmp1,P_supg,-1,1,-1);
	    
	    veccontr.axpy(tmp4,wpgdet);
	    
	    tmp19.set(P_supg).scale(wpgdet);
	    if (use_Ajac_old)
	      tmp20.prod(tmp19,Ao_grad_N,1,-1,2,-1,3);
	    else tmp20.prod(tmp19,A_grad_N,1,-1,2,-1,3);
	    matlocf.add(tmp20);
	    
	    if(!lumped_mass) {
	      if(compute_reactive_terms){
		// Reactive term in matrix (SUPG term)
		adv_diff_ff->comp_N_P_C(N_P_C,P_supg,SHAPE,wpgdet);
		matlocf.add(N_P_C);
	      }
	      
	      tmp21.set(SHAPE).scale(wpgdet*rec_Dt_m/ALPHA);
	      adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
	      tmp22.prod(P_Cp,tmp21,1,3,2);
	      matlocf.add(tmp22);
	      
	      tmp_ALE_07.prod(P_Cp,tmp_ALE_01,1,3,2);
	      matlocf.axpy(tmp_ALE_07,-wpgdet);
	      tmp_ALE_06.prod(P_Cp,tmp_ALE_02,1,-1,-1);
	      veccontr.axpy(tmp_ALE_06,wpgdet);
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
      
      }// ipg loop
      
      volume_flag=0;
      if (comp_res) {
	veccontr.export_vals(element.ret_vector_values(*retval));
	// veccontr.print(nel,"veccontr:");
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
