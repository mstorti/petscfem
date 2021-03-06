//__INSERT_LICENSE__
//$Id: advdife-bubbly.cpp,v 1.3.56.1 2007/02/19 20:23:56 mstorti Exp $
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

#include "nwadvdif.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void detj_error(double &detJaco,int elem) {
  printf("Jacobian of element %d is negative or null\n"
	 " Jacobian: %f\n",elem,detJaco);
  detJaco = -detJaco;
  if (detJaco==0.) detJaco = 1.0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int NewAdvDif::ask(const char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

IdentityEF identity_ef;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "void AdvDifFF::get_log_vars(int,const int*)"
void NewAdvDifFF::get_log_vars(int &nlog_vars,const int *& log_vars) {
  const char *log_vars_entry;
  const int &ndof=ndof;
  elemset->get_entry("log_vars_list",log_vars_entry);
  VOID_IT(log_vars_v);
  string s;
  if (log_vars_entry) {
    s=string(log_vars_entry);	// Save local copy
    read_int_array(log_vars_v,log_vars_entry);
  }
  nlog_vars=log_vars_v.size();
  log_vars = &*log_vars_v.begin();
  int ierr=0;
  for (int j=0; j<nlog_vars; j++) {
    if (log_vars_v[j]<=0) {
      PetscPrintf(PETSCFEM_COMM_WORLD,"Non positive dof in "
		  "\"log_vars_list\" entry: dof %d\n",
		  log_vars_v[j]);
      ierr=1;
    } else if (log_vars_v[j]>ndof) {
      PetscPrintf(PETSCFEM_COMM_WORLD,"Dof grater that ndof in "
		  "\"log_vars_list\" entry: dof %d, ndof %d\n",
		  log_vars_v[j]);
      ierr=1;
    }
    if (ierr) {
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "Errors while reading \"log_vars_list\"\n");
      if (log_vars_entry)
	PetscPrintf(PETSCFEM_COMM_WORLD,
		    "In line \"%s\"\n",s.c_str());
      exit(1);
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "void log_transf()"
void log_transf(FastMat2 &true_lstate,const FastMat2 &lstate,
		const int nlog_vars,const int *log_vars) {
  // Copy to log_state
  true_lstate.set(lstate);
  // Transform only those fields in log_vars
  for (int k=0; k<nlog_vars; k++) {
    int dof=log_vars[k];
    true_lstate.ir(2,dof);
    true_lstate.fun(&exp);
  }
  true_lstate.ir(2);
}

NewAdvDifFF::NewAdvDifFF(const NewElemset *elemset_)
    : elemset(elemset_), enthalpy_fun(NULL) {
  // This is ugly!!
  // assert(new_adv_dif_elemset);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "NewAdvDif::before_assemble"
void NewAdvDif::
before_assemble(arg_data_list &arg_datav,Nodedata *nodedata,
		Dofmap *dofmap, const char *jobinfo,int myrank,
		int el_start,int el_last,int iter_mode,
		const TimeData *time_data) { }

void NewAdvDifFF::get_C(FastMat2 &C) {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "Using lumped needs definition for get_C() virtual function\n"
	      "in the flux function object.\n");
  assert(0);
}

void NewAdvDifFF::get_Cp(FastMat2 &Cp) {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "Using lumped needs definition for get_Cp() virtual function\n"
	      "in the flux function object.\n");
  assert(0);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "advective::assemble"
void NewAdvDif::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			     const Dofmap *dofmap,const char *jobinfo,
			     const ElementList &elemlist,
			     const TimeData *time_data) {

extern const char * jobinfo_fields;

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  int locdof,kldof,lldof;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  assert(npg>0);
  assert(ndim>0);

  int nelprops,nel,ndof;
  elem_params(nel,ndof,nelprops);
  int nen = nel*ndof;

  // Unpack Dofmap
  int neq,nnod;
  neq = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  // H is a generalized local property passed per node with the nodal
  // coordinates. In shallow water nH =1 and H is the depth. It is
  // needed in order to compute the source term. In 1D Euler it may be
  // the area section of the tube. Its gradient is needed for the
  // source term in the momentum eqs.
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),H(nH);

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  double *retvalt;

  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step.
  vector<double> *dtmin;
  double lambda_max;
  int jdtmin;
  GlobParam *glob_param;
  // The trapezoidal rule integration parameter
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)
  arg_data *staten,*stateo,*retval,*fdj_jac,*jac_prof,*Ajac;
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
    rec_Dt_m = 1./DT;
    if (glob_param->steady) rec_Dt_m = 0.;
    if (ADVDIF_CHECK_JAC)
      fdj_jac = &arg_data_v[++j];
  }

  FastMat2 matlocf(4,nel,ndof,nel,ndof),matlocf_mass(4,nel,ndof,nel,ndof);
  FastMat2 prof_nodes(2,nel,nel),prof_fields(2,ndof,ndof),matlocf_fix(4,nel,ndof,nel,ndof);
  FastMat2 Id_ndf(2,ndof,ndof),Id_nel(2,nel,nel),prof_fields_diag_fixed(2,ndof,ndof);

  //o Use the weak form for the Galerkin part of the advective term.
  NSGETOPTDEF(int,weak_form,1);
  //o Use lumped mass (used mainly to avoid oscillations for small time steps).
  NSGETOPTDEF(int,lumped_mass,0);
  //o Factor to weight the influence of the temporal term in the continuous phase
  NSGETOPTDEF(int,factor_dalpha_dt,0);
  //o Key for Using lumping only for gas phase
  NSGETOPTDEF(int,use_lumping_only_for_gas,0);

  int lumped_mass_phase = 0;
  int comp_liq  = !strcmp(jobinfo_fields,"liq");
  int comp_gas  = !strcmp(jobinfo_fields,"gas");
//  if (comp_gas && lumped_mass) lumped_mass_phase = 1;
  if (lumped_mass) lumped_mass_phase = 1;

  //o Add axisymmetric version for this particular elemset.
  NSGETOPTDEF(string,axisymmetric,"none");
  assert(axisymmetric.length()>0);
  if (axisymmetric=="none") axi=0;
  else if (axisymmetric=="x") axi=1;
  else if (axisymmetric=="y") axi=2;
  else if (axisymmetric=="z") axi=3;
  else {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Invalid value for \"axisymmetric\" option\n"
		"axisymmetric=\"%s\"\n",axisymmetric.c_str());
    PetscFinalize();
    exit(0);
  }

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

  // Not implemented yet:= not lumped_mass_phase + log-vars
  assert(!use_log_vars || lumped_mass_phase);
  // lumped_mass_phase:= If this options is activated then all the inertia
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
  double lambda_max_pg;

  dshapex.resize(2,ndimel,nel);
  FMatrix Jaco(ndimel,ndim),Jaco_av(ndimel,ndim),
    iJaco(ndimel,ndimel),
    flux(ndof,ndimel),fluxd(ndof,ndimel),mass(nel,nel),
    grad_U(ndimel,ndof),A_grad_U(ndof),
    G_source(ndof), dUdt(ndof), Un(ndof),
    Ho(ndof),Hn(ndof);
  // These are edclared but not used
  FMatrix nor,lambda,Vr,Vr_inv,U(ndof),Ualpha(ndof),
    lmass(nel),Id_ndof(ndof,ndof),
    tmp1,tmp2,tmp3,tmp4,tmp5,hvec(ndimel),tmp6,tmp7,
    tmp8,tmp9,tmp10,tmp11(ndof,ndimel),tmp12,tmp14,
    tmp15,tmp17,tmp19,tmp20,tmp21,tmp22,tmp23,
    tmp24;
  FMatrix tmp30;

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

  //  FastMat2 Cp(2,ndof,ndof),Cr(2,ndof,ndof);
  FastMat2 Cr(2,ndof,ndof);
  FastMat2 Cp_bis(2,ndof,ndof);

  if (axi) assert(ndim==3);

  //#define COMPUTE_FD_ADV_JACOBIAN
#ifdef COMPUTE_FD_ADV_JACOBIAN
  FastMat2 A_fd_jac(3,ndim,ndof,ndof),U_pert(1,ndof),flux_pert(2,ndof,ndimel);
#endif
  FastMat2 U0(1,ndof);
  double u0[] = {1.,1.,1.,2.,0.,0.,0.1,0.1};
  U0.set(u0);

  Uo.resize(1,ndof);
  FastMat2 Uo_mod;
  Uo_mod.resize(1,ndof);

  Id_ndof.set(0.);
  for (int j=1; j<=ndof; j++) Id_ndof.setel(1.,j,j);

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int indx_l1[] = {1,3,4,5};
  int indx_g1[] = {2,6,7,8};
  int indx_l2[] = {1,3,4};
  int indx_g2[] = {2,5,6};
  int *indx_l, *indx_g;

  if(ndim==3) {
	indx_l=indx_l1;
	indx_g=indx_g1;
  } else{
	indx_l=indx_l2;
	indx_g=indx_g2;
	  }

  // printf("[%d] %s start: %d last: %d\n",MY_RANK,jobinfo,el_start,el_last);
  for (ElementIterator element = elemlist.begin();
       element!=elemlist.end(); element++) {

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
      // log_transf(true_lstaten,lstaten,nlog_vars,log_vars);
      // log_transf(true_lstateo,lstateo,nlog_vars,log_vars);
    }

    // State at time t_{n+\alpha}
    lstate.set(0.).axpy(lstaten,ALPHA).axpy(lstateo,(1-ALPHA));
    log_transf(true_lstate ,lstate ,nlog_vars,log_vars);
    log_transf(true_lstateo,lstateo,nlog_vars,log_vars);

    veccontr.set(0.);
    mass.set(0.);
    lmass.set(0.);
    matlocf.set(0.);
    //    if (lumped_mass_phase) matlocf_mass.set(0.);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    if (0){
    // DEBUG
	int kk,ielhh;
	element.position(kk,ielhh);
	printf("Element %d \n",kk);
        lstate.print("Estado :");
    // END DEBUG
    }

    // loop over Gauss points

    Jaco_av.set(0.);
    for (ipg=0; ipg<npg; ipg++) {

      //      Matrix xpg = SHAPE * xloc;
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

      if (detJaco <= 0.) {
	int k,ielh;
	element.position(k,ielh);
	printf("advdife: Jacobian of element %d is negative or null\n"
	       " Jacobian: %f\n",k,detJaco);
	PetscFinalize();
	exit(0);
      }
      wpgdet = detJaco*WPG;

      // Set volume of element. This is somewhat incorrect
      // because we should loop over Gauss points.
      // Correct for simplices (tris and tetras) and not
      // deformed quads and hexas.
      if (axi >0){
        ind_axi_1 = (  axi   % 3)+1;
        ind_axi_2 = ((axi+1) % 3)+1;
        Jaco_axi.setel(Jaco.get(ind_axi_1,ind_axi_1),1,1);
        Jaco_axi.setel(Jaco.get(ind_axi_1,ind_axi_2),1,2);
        Jaco_axi.setel(Jaco.get(ind_axi_2,ind_axi_1),2,1);
        Jaco_axi.setel(Jaco.get(ind_axi_2,ind_axi_2),2,2);
        detJaco_axi = Jaco_axi.det();
        double wpgdet_axi = detJaco_axi*WPG;
        Volume = 0.5*double(npg)*fabs(wpgdet_axi);

      } else {
        Volume = double(npg)*wpgdet;
      }

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

	Uo_mod.set(Uo);
	if(comp_liq){
		Uo_mod.ir(1,2);
		Uo.ir(1,2);
		Un.ir(1,2);
		Uo_mod.set(0.).axpy(Uo,1.0-factor_dalpha_dt).axpy(Un,factor_dalpha_dt);
		Uo.rs();
		Un.rs();
		Uo_mod.rs();
		}


	adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo_mod);

/*
	Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
	dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);
	for (int k=0; k<nlog_vars; k++) {
	  int jdof=log_vars[k];
	  double UU=exp(Ualpha.get(jdof));
	  dUdt.ir(1,jdof).scale(UU);
	}
	dUdt.rs();
*/

	// Pass to the flux function the true positive values
	U.prod(SHAPE,true_lstate,-1,-1,1);
	grad_U.prod(dshapex,true_lstate,1,-1,-1,2);
	grad_Uo.prod(dshapex,true_lstateo,1,-1,-1,2);

	delta_sc=0;
	//	double lambda_max_pg;

	// Compute A_grad_U in the `old' state
	//	adv_diff_ff->set_state(Uo,grad_Uold);
	adv_diff_ff->set_state(Uo,grad_Uo);
	adv_diff_ff->compute_flux(Uo,iJaco,H,grad_H,flux,fluxd,
				  A_grad_U,grad_Uo,G_source,
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

    if (weak_form==1 || weak_form==-1) {
	Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
	dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);
	for (int k=0; k<nlog_vars; k++) {
	  int jdof=log_vars[k];
	  double UU=exp(Ualpha.get(jdof));
	  dUdt.ir(1,jdof).scale(UU);
	}
	dUdt.rs();
} else {
	adv_diff_ff->get_Cp(Cp_bis);
	// usamos Ualpha como Delta U aqui
	Ualpha.set(Un).minus(Uo);
	dUdt.prod(Cp_bis,Ualpha,1,-1,-1).scale(rec_Dt_m);
}

	// Set the state of the fluid so that it can be used to
 	// compute matrix products
	adv_diff_ff->set_state(U,grad_U);
	adv_diff_ff->enthalpy_fun->set_state(U);
	if (!lumped_mass_phase) {
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

#ifdef COMPUTE_FD_ADV_JACOBIAN
	double eps_fd=1e-8;
	for (int jdof=1; jdof<=ndof; jdof++) {
	  U_pert.set(U).is(1,jdof).add(eps_fd).rs();
	  adv_diff_ff->set_state(U_pert,grad_U);
	  adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux_pert,fluxd,
				    A_grad_U,grad_U,G_source,
				    tau_supg,delta_sc,
				    lambda_max_pg, nor,lambda,Vr,Vr_inv,0);
	  flux_pert.minus(flux).scale(1./eps_fd);
	  flux_pert.t();
	  A_fd_jac.ir(3,jdof).set(flux_pert).rs();
	  flux_pert.rs();
	}
#endif

	if (lambda_max_pg>lambda_max) lambda_max=lambda_max_pg;

	if (!lumped_mass_phase) {
	  tmp10.set(G_source);	// tmp10 = G - dUdt
	  tmp10.minus(dUdt);
	} else {
	  tmp10.set(G_source);	// tmp10 = G (only the non lumped part of source)

	  if (use_lumping_only_for_gas) {
	    // anulo la parte temporal de la fase gas que la introduzco con lumping
	    for (int i=0; i<=ndim; i++) {
	      dUdt.is(1,indx_g[i]);
	    }
	    dUdt.set(0.).rs();
	    tmp10.minus(dUdt);
	  }
	}

	tmp1.rs().set(tmp10).minus(A_grad_U); //tmp1= G - dUdt - A_grad_U

	if (!lumped_mass_phase) {
	  adv_diff_ff->enthalpy_fun->comp_W_Cp_N(N_Cp_N,SHAPE,SHAPE,wpgdet*rec_Dt_m);
	  matlocf.add(N_Cp_N);
	} else {

	  if (use_lumping_only_for_gas) {
	    // anulo la parte temporal de la fase gas que la introduzco sin lumping

	    adv_diff_ff->get_Cp(Cp_bis);

	    for (int ii=0; ii<=ndim; ii++) Cp_bis.is(1,indx_g[ii]);
	    for (int jj=0; jj<ndof; jj++) Cp_bis.is(2,jj+1);
	    Cp_bis.set(0.).rs();

	    for (int ii=0; ii<=ndim; ii++) Cp_bis.is(1,indx_l[ii]);
	    for (int ii=0; ii<=ndim; ii++) Cp_bis.is(2,indx_g[ii]);
	    Cp_bis.set(0.).rs();
	    //	adv_diff_ff->enthalpy_fun->comp_W_Cp_N(N_Cp_N,SHAPE,SHAPE,wpgdet*rec_Dt_m);
	    tmp30.prod(SHAPE,Cp_bis,1,2,3);
	    N_Cp_N.prod(tmp30,SHAPE,1,2,4,3).scale(wpgdet*rec_Dt_m);;
	    matlocf.add(N_Cp_N);

	  }
	}

	// A_grad_N.prod(dshapex,A_jac,-1,1,-1,2,3);
	adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);

	// Termino Galerkin
	if (weak_form) {
	  // weak version
	  tmp11.set(flux).minus(fluxd); // tmp11 = flux_c - flux_d

	  tmp23.set(SHAPE).scale(-wpgdet*ALPHA);
	  tmp14.prod(A_grad_N,tmp23,1,2,4,3);
	  matlocf.add(tmp14);
	} else {
	  // tmp2.prod(SHAPE,tmp1,1,2); // tmp2= SHAPE' * (G - dUdt - A_grad_U)
	  tmp2.prod(SHAPE,A_grad_U,1,2); // tmp2= SHAPE' * A_grad_U
	  veccontr.axpy(tmp2,-wpgdet);
	  tmp11.set(0.).minus(fluxd); // tmp11 = - flux_d

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

	// Debug and fixme later
#define COMPUTE_C_D_TERMS

#ifdef COMPUTE_C_D_TERMS

	// Diffusive term in matrix
#if 0
	tmp17.set(dshapex).scale(wpgdet*ALPHA);
	adv_diff_ff->comp_D_grad_N(D_grad_N,tmp17);
	tmp18.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
#endif

	adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
					  dshapex,wpgdet*ALPHA);
	matlocf.add(grad_N_D_grad_N);
#endif

        // Reactive term in matrix (Galerkin part)
#if 0
        tmp23.set(SHAPE).scale(wpgdet*ALPHA);
	tmp24.prod(SHAPE,tmp23,1,2); // tmp24 = SHAPE' * SHAPE
	tmp25.prod(tmp24,C_jac,1,3,2,4); // tmp25 = SHAPE' * SHAPE * C_jac
#endif


#ifdef COMPUTE_C_D_TERMS
	if (!lumped_mass_phase) {
	adv_diff_ff->comp_N_N_C(N_N_C,SHAPE,wpgdet*ALPHA);
	matlocf.add(N_N_C);
	}
#endif

#ifndef USE_OLD_STATE_FOR_P_SUPG
	// This computes either the standard `P_supg' perturbation
	// function or other written by the user in the
	// flux-function.
	adv_diff_ff->comp_P_supg(P_supg);
#endif

	    if (0){
	      // DEBUG
	      int kk,ielhh;
	      element.position(kk,ielhh);
	      printf("Element %d \n",kk);
	      dUdt.print("dUdt :");
	      G_source.print("G_source :");
	      // END DEBUG
	    }

// DEBUG
//	tmp1.rs().set(tmp10); //tmp1= G - dUdt - A_grad_U
// END DEBUG

	for (int jel=1; jel<=nel; jel++) {

   	  P_supg.ir(1,jel);

	  tmp4.prod(tmp1,P_supg,-1,1,-1);
	  //	  veccontr.ir(1,jel).axpy(tmp4,wpgdet).ir(1);
	  veccontr.ir(1,jel).axpy(tmp4,wpgdet);

	  matlocf.ir(1,jel);

	  tmp19.set(P_supg).scale(ALPHA*wpgdet);
	  tmp20.prod(tmp19,A_grad_N,1,-1,2,-1,3);
	  matlocf.add(tmp20);

#ifdef COMPUTE_C_D_TERMS
	  if(!lumped_mass_phase) {
	    // Reactive term in matrix (SUPG term)
	    adv_diff_ff->comp_N_P_C(N_P_C,P_supg,SHAPE,wpgdet*ALPHA);
	    matlocf.add(N_P_C);

	    tmp21.set(SHAPE).scale(wpgdet*rec_Dt_m);
	    adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
	    tmp22.prod(P_Cp,tmp21,1,3,2);
	    matlocf.add(tmp22);

	  }
#endif

	  if(lumped_mass_phase && use_lumping_only_for_gas) {
	    tmp21.set(SHAPE).scale(wpgdet*rec_Dt_m);
	    P_Cp.prod(P_supg,Cp_bis,1,-1,-1,2);
	    tmp22.prod(P_Cp,tmp21,1,3,2);
	    matlocf.add(tmp22);
	  }

	  /*
		    if(lumped_mass_phase) {
		    if (use_lumping_only_for_gas) {
		    // anulo la parte temporal de la fase gas que la introduzco sin lumping
		    tmp21.set(SHAPE).scale(wpgdet*rec_Dt_m);
		    P_Cp.prod(P_supg,Cp_bis,1,-1,-1,2);
		    tmp22.prod(P_Cp,tmp21,1,3,2);
		    matlocf.add(tmp22);
		    } else {

		    tmp21.set(SHAPE).scale(wpgdet*rec_Dt_m);
		    adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
		    tmp22.prod(P_Cp,tmp21,1,3,2);
		    matlocf.add(tmp22);
		    }
		    }
	  */

	  matlocf.rs();
	  veccontr.rs();
	}
	P_supg.rs();

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

      }
      volume_flag=0;

      // modif March 26
      // lumped mass
      lmass.axpy(SHAPE,wpgdet);

    }  // ipg loop

    if (lumped_mass_phase) {
      // lumped terms treatment
      // temporal and source terms are lumped out of the gauss points loop
      for (int j=1; j<=nel; j++) {


	lstaten.ir(1,j);
	Un.set(lstaten);
	lstaten.rs();

	lstateo.ir(1,j);
	Uo.set(lstateo);
	lstateo.rs();
	//      adv_diff_ff->set_state(Uo,grad_U);
	adv_diff_ff->set_state(Uo);

	Uo_mod.set(Uo);
	if(comp_liq){
	  Uo_mod.ir(1,2);
	  Uo.ir(1,2);
	  Un.ir(1,2);
	  Uo_mod.set(0.).axpy(Uo,1.0-factor_dalpha_dt).axpy(Un,factor_dalpha_dt);
	  Uo.rs();
	  Un.rs();
	  Uo_mod.rs();
	}

	// adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
	adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo_mod);

	adv_diff_ff->set_state(Un,grad_U);
	adv_diff_ff->compute_flux(Un,iJaco,H,grad_H,flux,fluxd,
				  A_grad_U,grad_U,G_source,
				  tau_supg,delta_sc,
				  lambda_max_pg, nor,lambda,Vr,Vr_inv,
				  COMP_SOURCE_LUMPED);
	adv_diff_ff->enthalpy_fun->enthalpy(Hn,Un);

	Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);

	//      adv_diff_ff->enthalpy_fun->comp_P_Cp(Cp,Id_ndf);
	adv_diff_ff->get_Cp(Cp_bis);
	Cp_bis.scale(rec_Dt_m);

	if (weak_form==1 || weak_form==-1) {
	  dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);
	}
	else {
	  dUdt.prod(Cp_bis,Un,1,-1,-1);
	  Ho.prod(Cp_bis,Uo,1,-1,-1);
          dUdt.minus(Ho);
	}

	for (int k=0; k<nlog_vars; k++) {
	  int jdof=log_vars[k];
	  double UU=exp(Ualpha.get(jdof));
	  dUdt.ir(1,jdof).scale(UU);
	}
	dUdt.rs();

	if (use_lumping_only_for_gas){
	  // anulo la parte temporal de la fase liquida que la introduzco sin lumping
	  for (int ii=0; ii<=ndim; ii++) {
	    dUdt.is(1,indx_l[ii]);
	  }
	  dUdt.set(0.).rs();
	}

	tmp10.set(G_source).minus(dUdt);	// tmp10 = G - dUdt

	adv_diff_ff->get_C(Cr);

	if (use_lumping_only_for_gas){

	for (int ii=0; ii<=ndim; ii++) Cp_bis.is(1,indx_l[ii]);
	for (int jj=0; jj<ndof; jj++) Cp_bis.is(2,jj+1);
	Cp_bis.set(0.).rs();

	//	for (int ii=0; ii<=ndim; ii++) Cp_bis.is(1,indx_g[ii]);
	Cp_bis.is(1,1,ndof);
	for (int ii=0; ii<=ndim; ii++) Cp_bis.is(2,indx_l[ii]);
	Cp_bis.set(0.).rs();

	}

      tmp12.set(Cr).add(Cp_bis);

      double mass_node;
      mass_node = double(lmass.get(j));
      veccontr.ir(1,j).axpy(tmp10,mass_node).rs();
      matlocf.ir(1,j).ir(3,j).axpy(tmp12,mass_node).rs();

      }
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

#if 0 // Lumping + log var to be implemented in the future ... ???
      if (lumped_mass_phase) {
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
	dUloc.set(lstaten).minus(lstateo);
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
	veccontr.minus(veccontr_mass);
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
#endif

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

// DEBUG matrix
    if (0){
	int kk,ielhh;
	element.position(kk,ielhh);
	printf("Element %d \n",kk);
        lstate.print("Estado :");
        matlocf.print("System Matrix :");
    }
// END DEBUG

  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "NewAdvDif::after_assemble"
void NewAdvDif::
after_assemble(const char *jobinfo) { }

double NewAdvDif::volume() const {
  if (volume_flag) {
    return Volume;
  } else {
    assert(0); // Not valid at this moment;
  }
}

const FastMat2 *NewAdvDif::grad_N() const {
  return &dshapex;
}

#if 0
void NewAdvDif::comp_P_supg(int is_tau_scalar) {
  if (is_tau_scalar) {
    double tau_supg_d = tau_supg.get(1,1);
    P_supg.set(Ao_grad_N).scale(tau_supg_d);
  } else {
    P_supg.prod(Ao_grad_N,tau_supg,1,2,-1,-1,3);
  }
}
#endif

void NewAdvDifFF::comp_P_supg(FastMat2 &P_supg) {
  assert(new_adv_dif_elemset);
  const NewAdvDif *e = new_adv_dif_elemset;
  if (e->ff_options & SCALAR_TAU) {
    double tau_supg_d = e->tau_supg.get(1,1);
    P_supg.set(e->Ao_grad_N).scale(tau_supg_d);
  } else {
    P_supg.prod(e->Ao_grad_N,e->tau_supg,1,2,-1,-1,3);
  }
}

void NewAdvDifFF::set_profile(FastMat2 &seed) {
  seed.set(1.);
}

void NewAdvDifFF::
Riemann_Inv(const FastMat2 &U, const FastMat2 &normal,
	    FastMat2 &Rie, FastMat2 &drdU, FastMat2 &C_) { }

#undef SHAPE
#undef DSHAPEXI
#undef WPG
