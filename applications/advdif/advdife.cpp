//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;

#include <vector>
#include <string>
#include <typeinfo>
#ifdef __GNUC__
#include <cxxabi.h>
#endif

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>
#include <src/generror.h>

#include "nwadvdif.h"

#define MAXPROP 100

extern GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#if 0
void detj_error(double &detJaco,int elem) {
  printf("Jacobian of element %d is negative or null\n"
	 " Jacobian: %f\n",elem,detJaco);
  detJaco = -detJaco;
  if (detJaco==0.) detJaco = 1.0;
}
#endif

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
  int nel,ndof,nelprops;
  elemset->elem_params(nel,ndof,nelprops);
  const char *log_vars_entry;
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
		  log_vars_v[j], ndof);
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
  : elemset(elemset_), new_adv_dif_elemset(NULL), enthalpy_fun(NULL) {
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
		const TimeData *time_data) {
  int ierr;
  //o Compute finite difference jacobian of fluxes for checking the
  //  analytical one. For each element the following norms are printed:
  //  analytical jacobian #|A_a|# , numerical jacobian #|A_n|# and the
  //  difference #|A_a-A_n|# . Incrementing #compute_fd_adv_jacobian==1#
  //  increases the verbosity. If #=1# the maximum values over all the
  //  elemset are printed. If #=2# the errors for all elements are
  //  reported. Finally, if #=3# also the jacobians themselves are
  //  printed. For 2 and 3, if #compute_fd_adv_jacobian_elem_list# is
  //  set, then only those elements are printed. If
  //  #compute_fd_adv_jacobian_rel_err_threshold# is set then only those
  //  elements for which the error is greater than the given threshold
  //  are reported.  Also, be warned that when run in parallel, printing
  //  for a lot of elements in different processors may be messy.
  NSGETOPTDEF_ND(int,compute_fd_adv_jacobian,0);
  //o The perturbation scale for computing the numerical jacobian
  //  (see #compute_fd_adv_jacobian# ).
  NSGETOPTDEF_ND(double,compute_fd_adv_jacobian_eps,1e-4);
  assert(compute_fd_adv_jacobian_eps > 0.);
  //o Report elements whose relative error in computing
  //  flux jacobians exceed these value.
  NSGETOPTDEF_ND(double,compute_fd_adv_jacobian_rel_err_threshold,0.);
  //o Use the GCL compliant versin of the algorithm 
  NSGETOPTDEF_ND(int,use_GCL_compliant,0);
  //o Use the ALE+GCL formulation versin of the algorithm 
  NSGETOPTDEF_ND(int,use_ALE_form,0);
  //o Use the BDF time integrator
  NSGETOPTDEF_ND(int,use_BDF_advdife,0);

#if 0
  //o Flags whether we are solving a precondioned
  // system with the dual time strategy
  NSGETOPTDEF_ND(int,precoflag,0);
#else
  precoflag = GLOB_PARAM->precoflag;
#endif

  A_fd_jac_norm_max = 0.;
  A_fd_jac_norm_min = DBL_MAX;
  A_jac_norm_max = 0.;
  A_jac_norm_min = DBL_MAX;
  A_jac_err_norm_max = 0.;
  A_jac_err_norm_min = DBL_MAX;
  A_rel_err_max = 0.;
  A_rel_err_min = DBL_MAX;
  comp_checked=0;
  comp_total=0;

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This function prints the name of the
// underlying flux function class. I read somewhere
// that GCC uses the C++ ABI interface for RTTI, and then
// I protect this with the __GNUC__ macro. 
static void print_flux_fun_name(NewAdvDifFF *fun) {
#ifdef __GNUC__
  size_t len;
  int s;
  char* p=abi::__cxa_demangle(typeid(*fun).name(), 
                              0, &len, &s);
  PetscPrintf(PETSCFEM_COMM_WORLD,
              "Flux function class \"%s\"\n",p);
  free(p);
#else
  PetscPrintf(PETSCFEM_COMM_WORLD,
              "Can not determine flux function class.\n");
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void NewAdvDifFF::get_C(FastMat2 &C) {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "Not defined get_C() virtual function\n"
	      "in the flux function object.\n");
  print_flux_fun_name(this);
  assert(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void NewAdvDifFF::get_Cp(FastMat2 &Cp) {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "Not defined get_Cp() virtual function\n"
	      "in the flux function object.\n");
  print_flux_fun_name(this);
  assert(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void NewAdvDifFF::get_Ajac(FastMat2 &Ajac) {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "Not defined get_Ajac() virtual function\n"
	      "in the flux function object.\n");
  print_flux_fun_name(this);
}

void NewAdvDifFF::compute_delta_sc_v(FastMat2 &delta_sc_v) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "NewAdvDif::assemble"
void NewAdvDif::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			     const Dofmap *dofmap,const char *jobinfo,
			     const ElementList &elemlist,
			     const TimeData *time_data) try {

  if (use_ALE_form) {
    new_assemble_ALE_formulation(arg_data_v,nodedata,
				 dofmap,jobinfo,elemlist,time_data);
    return;
  }
  if (use_BDF_advdife) {
    new_assemble_BDF(arg_data_v,nodedata,
                     dofmap,jobinfo,elemlist,time_data);
    return;
  }
  if (use_GCL_compliant) {
    new_assemble_GCL_compliant(arg_data_v,nodedata,
                               dofmap,jobinfo,elemlist,time_data);
    return;
  } else if (precoflag) {
    new_assemble_preco(arg_data_v,nodedata,
                       dofmap,jobinfo,elemlist,time_data);
    return;
  }

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF(int,ndim,0); //nd
  assert(npg>0);
  assert(ndim>0);

  int nelprops,nel,ndof;
  elem_params(nel,ndof,nelprops);
  // int nen = nel*ndof;

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
  //  FMatrix  Hloc(nel,nH),H(nH),vloc_mesh(nel,ndim),v_mesh(ndim);
  FMatrix  Hloc(nel,nH),H(nH),vloc_mesh(nel,ndim);
  //  FastMat2 v_mesh;

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  time_m = double(* (const Time *) time_data);

  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step.
  double lambda_max=NAN;
  int PFUNUSED jdtmin;
  GlobParam *glob_param=NULL;
  // The trapezoidal rule integration parameter
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)
  arg_data *staten=NULL,*stateo=NULL,*retval=NULL,
    *jac_prof=NULL, *Ajac=NULL, *fdj_jac=NULL;
  if (comp_res) {
    int j=-1;
    stateo = &arg_data_v[++j]; //[0]
    staten = &arg_data_v[++j]; //[1]
    retval  = &arg_data_v[++j];//[2]
    j++; // jdtmin = ++j;//[3]
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];//[4]
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
    rec_Dt_m = 1./DT;
    if (glob_param->steady) rec_Dt_m = 0.;

    // apply ALPHA to time step in order to take
    // into account higher order temporal integration
    rec_Dt_m = rec_Dt_m/ALPHA;

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
  //o Flag to turn on ALE (Arbitrary Lagrangian-Eulerian) computation. 
  NSGETOPTDEF_ND(int,ALE_flag,0);
  //o Pointer to old coordinates in
  //  #nodedata# array excluding the first "ndim" values
  NSGETOPTDEF(int,indx_ALE_xold,1);
  //o Compute the term non-symmetric term
  //  correspoding to nonlinearities in the
  //  diffusive Jacobian.
  NSGETOPTDEF(int,compute_dDdU_term,1);
#if 0
  // NO: THIS IS ALREADY IN use_Ajac_old
  //o If activated, the stability term is computed
  //  in the previous time step. The advantage is that
  //  you have true Newton convergence (quadratic), but
  //  the parameters are not updated while the field is moving. 
  NSGETOPTDEF(int,use_old_state_for_p_supg,1);
#endif

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
    dUloc(nel,ndof),matloc;
  FastMat2 true_lstate(2,nel,ndof),
    true_lstateo(2,nel,ndof),true_lstaten(2,nel,ndof);

  FastMat2 true_lstate_abs(2,nel,ndof);

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg,GP_FASTMAT2);
  GPdata gp_data_low(geometry.c_str(),ndimel,nel,1,GP_FASTMAT2);

  double detJaco, wpgdet, delta_sc, delta_sc_old;
  int ipg;
  double lambda_max_pg;

  dshapex.resize(2,ndimel,nel);
  FMatrix Jaco(ndimel,ndim),Jaco_av(ndimel,ndim),
    iJaco(ndimel,ndimel),
    flux(ndof,ndimel),fluxd(ndof,ndimel),mass(nel,nel),
    grad_U(ndimel,ndof), A_grad_U(ndof),Ao_grad_U(ndof),
    G_source(ndof), dUdt(ndof), Un(ndof),
    Ho(ndof),Hn(ndof);
  // FMatrix grad_U_norm(ndimel,ndof);
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

  // FastMat2 vaux(2,ndim,nel);

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
      }

      // State at time t_{n+\alpha}
      lstate.set(0.).axpy(lstaten,ALPHA).axpy(lstateo,(1-ALPHA));
      log_transf(true_lstate ,lstate ,nlog_vars,log_vars);
      log_transf(true_lstateo,lstateo,nlog_vars,log_vars);

      veccontr.set(0.);
      mass.set(0.);
      lmass.set(0.);
      matlocf.set(0.);

      v_mesh.set(0.);

      if (lumped_mass) matlocf_mass.set(0.);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE	 (*gp_data.FM2_shape[ipg])
#define WPG	 (gp_data.wpg[ipg])

#define DSHAPEXI_LOW (*gp_data_low.FM2_dshapexi[0])
#define SHAPE_LOW    (*gp_data_low.FM2_shape[0])
#define WPG_LOW	     (gp_data_low.wpg[0])

#undef SHV
#define SHV(mess,v)                                     \
      printf("[%d] elem %d, %s " #v " %f\n",            \
             MY_RANK,k_elem,mess,v.sum_square_all(),    \
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

      if (0){
        // DEBUG
        int kk,ielhh;
        element.position(kk,ielhh);
        printf("Element %d \n",kk);
        lstate.print("Estado :");
        // END DEBUG
      }
    
      // nodal computation of mesh velocity
      if (ALE_flag) {
        PETSCFEM_ASSERT(nH >= ndim,"This element requires the old mesh "
                        "position to be passed as an H field. nH %d, ndim %d",
                        nH,ndim);
        PETSCFEM_ASSERT(indx_ALE_xold>=1 && nH-indx_ALE_xold+1 >= ndim,
                        "bad indx_ALE_xold, not remaining enough columns. "
                        "indx_ALE_xold %d, nH %d, ndim %d",
                        indx_ALE_xold,nH,ndim);  
        Hloc.is(2,indx_ALE_xold,indx_ALE_xold+ndim-1);
        vloc_mesh.set(xloc).minus(Hloc).scale(rec_Dt_m*ALPHA).rs();
        Hloc.rs();
      }

      if (0){
        int kk,ielhh;
        element.position(kk,ielhh);
        if (rand() % 50 == 0){
          printf("Element %d \n",kk);
          xloc.print(" xloc^(n+1) :");
          Hloc.print(" Hloc^(n+1) :");
          vloc_mesh.print(" vloc_mesh^(n+1) :");
        }
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

          if (ALE_flag) v_mesh.prod(SHAPE_LOW,vloc_mesh,-1,-1,1);

          adv_diff_ff->compute_flux(Uo,iJaco_low,H,grad_H,flux,fluxd,
                                    A_grad_U,grad_Uo,G_source,
                                    tau_supg,delta_sc_old,
                                    lambda_max_pg, nor,lambda,Vr,Vr_inv,
                                    COMP_SOURCE | COMP_UPWIND);
          adv_diff_ff->comp_A_grad_N(Ao_grad_N,dshapex);
	
          if (use_Ajac_old) adv_diff_ff->get_Ajac(Ao);
	
          // SHV("antes de comp_P_supg",P_supg);
          adv_diff_ff->comp_P_supg(P_supg);
          // SHV("despues de comp_P_supg",P_supg);
	
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
          adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);
	
          if (ALE_flag) {
            //	  v_mesh.prod(SHAPE_LOW,vloc_mesh,-1,-1,1);
            adv_diff_ff->get_Cp(Cp_bis);
            tmp_ALE_01.prod(v_mesh,dshapex,-1,-1,1);
            tmp_ALE_02.prod(v_mesh,grad_U,-1,-1,1);
          }

          //	  if (use_Ajac_old) Ao_grad_U.prod(Ao,grad_U,-1,1,-2,-1,-2);

          tmp11.set(0.).minus(fluxd);
          tmp8.prod(dshapex,tmp11,-1,1,2,-1);
          veccontr.axpy(tmp8,wpgdet_low);
	
          adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
                                            dshapex,wpgdet_low);
          matlocf.add(grad_N_D_grad_N);
	
	  dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);
	  
	  tmp10.set(G_source);	// tmp10 = G - dUdt
	  if (!lumped_mass) tmp10.minus(dUdt);
	  
          tmp1.rs().set(tmp10).minus(A_grad_U); //tmp1= G - dUdt - A_grad_U
	  
          if (use_Ajac_old) {
            Ao_grad_U.prod(Ao,grad_U,-1,1,-2,-1,-2);
            tmp1_old.rs().set(tmp10).minus(Ao_grad_U); //tmp1= G - dUdt - A_grad_U
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
	  
            tmp19.set(P_supg).scale(wpgdet_low);
            if (use_Ajac_old)
              tmp20.prod(tmp19,Ao_grad_N,1,-1,2,-1,3);
            else tmp20.prod(tmp19,A_grad_N,1,-1,2,-1,3);
            matlocf.add(tmp20);
	  
            if(!lumped_mass) {
	    
              if(compute_reactive_terms){
                // Reactive term in matrix (SUPG term)
                adv_diff_ff->comp_N_P_C(N_P_C,P_supg,SHAPE_LOW,wpgdet_low);
                matlocf.add(N_P_C);
              }
	    
              tmp21.set(SHAPE_LOW).scale(wpgdet_low*rec_Dt_m);
              adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
              tmp22.prod(P_Cp,tmp21,1,3,2);
              matlocf.add(tmp22);
	    
              if (ALE_flag) {
#if 0
                tmp_ALE_07.prod(P_supg,tmp_ALE_01,1,3,2);
                matlocf.axpy(tmp_ALE_07,-wpgdet_low);
                tmp_ALE_06.prod(P_supg,tmp_ALE_02,1,-1,-1);
                veccontr.axpy(tmp_ALE_06,wpgdet_low);
#else
                tmp_ALE_07.prod(P_Cp,tmp_ALE_01,1,3,2);
                matlocf.axpy(tmp_ALE_07,-wpgdet_low);
                tmp_ALE_06.prod(P_Cp,tmp_ALE_02,1,-1,-1);
                veccontr.axpy(tmp_ALE_06,wpgdet_low);
#endif
              }
	    
            }
          }
	
          matlocf.rs();
          P_supg.rs();
          veccontr.rs();
	
          //	volume_flag = 0;
	
        }	      
      }
#if 0
      SHV("despues de use_low_gpdata block",matlocf);
      SHV("despues de use_low_gpdata block",veccontr);
#endif

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

        if (detJaco<=0.) {
          int k,ielh;
          element.position(k,ielh);
          detj_error(detJaco,k);
          set_error(1);
        }
        wpgdet = detJaco*WPG;

        /*
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
        */

        dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);

        if (nH>0) {
          H.prod(SHAPE,Hloc,-1,-1,1);
          grad_H.prod(dshapex,Hloc,1,-1,-1,2);
        }

        if (comp_res) {

          // state variables and gradient
          //	Un.prod(SHAPE,lstaten,-1,-1,1);
          Un.prod(SHAPE,lstate,-1,-1,1);
          adv_diff_ff->enthalpy_fun->enthalpy(Hn,Un);
          Uo.prod(SHAPE,lstateo,-1,-1,1);
          adv_diff_ff->enthalpy_fun->enthalpy(Ho,Uo);
          Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
          dUdt.set(Hn).minus(Ho).scale(rec_Dt_m);

          /*
            adv_diff_ff->get_Cp(Cp_bis);
            dUdt.set(Un).minus(Uo).scale(rec_Dt_m);
            dUdt.prod(Cp_bis,dUdt,1,-1,-1);
          */

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

          if (ALE_flag) v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);

          adv_diff_ff->compute_flux(Uo,iJaco,H,grad_H,flux,fluxd,
                                    A_grad_U,grad_Uo,G_source,
                                    tau_supg,delta_sc_old,
                                    lambda_max_pg, nor,lambda,Vr,Vr_inv,
                                    COMP_SOURCE | COMP_UPWIND);
          adv_diff_ff->comp_A_grad_N(Ao_grad_N,dshapex);
          Ao_grad_U.set(A_grad_U);
          if (shocap>0. || shocap_aniso>0.)
            adv_diff_ff->get_Cp(Cp_bis_old);
          if (use_Ajac_old) adv_diff_ff->get_Ajac(Ao);
          adv_diff_ff
            ->compute_shock_cap_aniso(delta_aniso_old,jvec_old);

          // This computes either the standard `P_supg' perturbation
          // function or other written by the user in the
          // flux-function.
          if (use_Ajac_old)
            adv_diff_ff->comp_P_supg(P_supg);

          // Set the state of the fluid so that it can be used to
          // compute matrix products
          adv_diff_ff->set_state(U,grad_U);
          adv_diff_ff->enthalpy_fun->set_state(U);

          // MODIF BETO 8/6
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
          /*
            adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux,fluxd,
            A_grad_U,grad_U,G_source,
            tau_supg,delta_sc,
            lambda_max_pg, nor,lambda,Vr,Vr_inv,
            COMP_SOURCE | COMP_UPWIND);
          */
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

          tmp10.set(G_source);	// tmp10 = G - dUdt
          if (!lumped_mass) tmp10.minus(dUdt);

          tmp1.rs().set(tmp10).minus(A_grad_U); //tmp1= G - dUdt - A_grad_U

          if (use_Ajac_old) {
            Ao_grad_U.prod(Ao,grad_U,-1,1,-2,-1,-2);
            tmp1_old.rs().set(tmp10).minus(Ao_grad_U); //tmp1= G - dUdt - A_grad_U
          }

	// MODIF BETO 8/6
	if (!lumped_mass) {
	  adv_diff_ff->enthalpy_fun
	    ->comp_W_Cp_N(N_Cp_N,SHAPE,SHAPE,
			  wpgdet*rec_Dt_m);
	  matlocf.add(N_Cp_N);
	}
	
	// A_grad_N.prod(dshapex,A_jac,-1,1,-1,2,3);
	adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);

	// add ALE Galerkin terms
	if (ALE_flag) {
	  // v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	  adv_diff_ff->get_Cp(Cp_bis);
	  tmp_ALE_01.prod(v_mesh,dshapex,-1,-1,1);
	  tmp_ALE_02.prod(v_mesh,grad_U,-1,-1,1);
	  tmp_ALE_03.prod(SHAPE,Cp_bis,1,2,3);

	  tmp_ALE_04.prod(tmp_ALE_03,tmp_ALE_02,1,2,-1,-1);
	  veccontr.axpy(tmp_ALE_04,wpgdet);

	  tmp_ALE_05.prod(tmp_ALE_03,tmp_ALE_01,1,2,4,3);
	  matlocf.axpy(tmp_ALE_05,-wpgdet);
	}

	// Termino Galerkin
	if (weak_form) {
	  // assert(!lumped_mass && beta_supg==1.); 
	  // Not implemented yet!!
	  // weak version

	  //	  tmp11.set(flux).minus(fluxd); // tmp11 = flux_c - flux_d
	  if (use_low_gpdata){
	  tmp11.set(flux); // tmp11 = flux_c (viscous part is integrated in 1 PG)
	  } else {
	  tmp11.set(flux).minus(fluxd); // tmp11 = flux_c - flux_d
	  }

            tmp23.set(SHAPE).scale(-wpgdet);
            tmp14.prod(A_grad_N,tmp23,1,2,4,3);
            matlocf.add(tmp14);
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

            tmp23.set(SHAPE).scale(wpgdet);
            tmp14.prod(A_grad_N,tmp23,3,2,4,1);
            matlocf.add(tmp14);
          }
          // tmp8= DSHAPEX * (w*flux_c - flux_d)
          //	      w = weak_form
          tmp8.prod(dshapex,tmp11,-1,1,2,-1);
          tmp9.prod(SHAPE,tmp10,1,2); // tmp9 = SHAPE' * (G - dUdt)
          tmp8.add(tmp9);		// tmp8 = DSHAPEX * tmp11
          veccontr.axpy(tmp8,wpgdet);

          // Diffusive term in matrix
#if 0
          tmp17.set(dshapex).scale(wpgdet);
          adv_diff_ff->comp_D_grad_N(D_grad_N,tmp17);
          tmp18.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
#endif
          if(use_low_gpdata==0){
            adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
                                              dshapex,wpgdet);
            matlocf.add(grad_N_D_grad_N);
          }

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

          // MODIF BETO 8/6
          if(compute_reactive_terms){
            if (!lumped_mass) {
              adv_diff_ff->comp_N_N_C(N_N_C,SHAPE,wpgdet);
              matlocf.add(N_N_C);
            }
          }

          // adding shock-capturing term
          delta_sc_v.set(0.0);
          if (shocap>0. ) {
            adv_diff_ff->compute_delta_sc_v(delta_sc_v);

            /*
              adv_diff_ff->get_Cp(Cp_bis);
              tmp_shc_grad_U.prod(Cp_bis,grad_U,2,-1,1,-1);
              for (int jdf=1; jdf<=ndof; jdf++) {
              delta_sc_v.addel(delta_sc,jdf);
            */

            //	  Cp_bis.set(Cp_bis_old);
            tmp_shc_grad_U.prod(Cp_bis_old,grad_U,2,-1,1,-1);
            for (int jdf=1; jdf<=ndof; jdf++) {
              //	    delta_sc_v.addel(delta_sc,jdf);
              delta_sc_v.addel(delta_sc_old,jdf);
            }

            //		  }

            tmp_sc.prod(dshapex,dshapex,-1,1,-1,2).scale(shocap*wpgdet);
            //	  tmp_sc_v.prod(dshapex,grad_U,-1,1,-1,2);
            tmp_sc_v.prod(dshapex,tmp_shc_grad_U,-1,1,-1,2);
            for (int jdf=1; jdf<=ndof; jdf++) {
              double delta = (double)delta_sc_v.get(jdf);
              for (int kdf=1; kdf<=ndof; kdf++) {
                //	      double tmp_shc_1=Cp_bis.get(jdf,kdf);
                double tmp_shc_1=Cp_bis_old.get(jdf,kdf);
                matlocf.ir(2,jdf).ir(4,kdf).axpy(tmp_sc,delta*tmp_shc_1).rs();
                //       matlocf.ir(2,jdf).ir(4,jdf).axpy(tmp_sc,delta).rs();
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
            //	  double delta_aniso = 0.0;
            delta_aniso = 0.0;
            adv_diff_ff
              ->compute_shock_cap_aniso(delta_aniso,jvec);
#if 0
            double jvec_norm = sqrt(jvec.sum_square_all());
            assert(jvec_norm>0.);// fixme:= lanzar execpcion
            jvec.scale(1./jvec_norm);
#endif

#if 0
            adv_diff_ff->get_Cp(Cp_bis);
            tmp_shc_grad_U.prod(Cp_bis,grad_U,2,-1,1,-1);
            tmp_j_grad_U.prod(jvec,tmp_shc_grad_U,-1,-1,1);
            tmp_j_gradN.prod(jvec,dshapex,-1,-1,1);

            tmp_sc_aniso.prod(tmp_j_gradN,tmp_j_gradN,1,2)
              .scale(shocap_aniso*wpgdet);

            tmp_sc_v_aniso.prod(tmp_j_gradN,tmp_j_grad_U,1,2);
#else
            tmp_shc_grad_U.prod(Cp_bis_old,grad_U,2,-1,1,-1);
            tmp_j_grad_U.prod(jvec_old,tmp_shc_grad_U,-1,-1,1);
            tmp_j_gradN.prod(jvec_old,dshapex,-1,-1,1);

            tmp_sc_aniso.prod(tmp_j_gradN,tmp_j_gradN,1,2)
              .scale(shocap_aniso*wpgdet);

            tmp_sc_v_aniso.prod(tmp_j_gradN,tmp_j_grad_U,1,2);

#endif

#if 0
            for (int jdf=1; jdf<=ndof; jdf++) {
              double delta = (double)delta_sc_v.get(jdf);
              for (int kdf=1; kdf<=ndof; kdf++) {
                double tmp_shc_1=Cp_bis.get(jdf,kdf);
                matlocf.ir(2,jdf).ir(4,kdf).axpy(tmp_sc,delta*tmp_shc_1).rs();
                //       matlocf.ir(2,jdf).ir(4,jdf).axpy(tmp_sc,delta).rs();
              }
              tmp_sc_v.ir(2,jdf).scale(-shocap*delta*wpgdet).rs();
              delta_sc_v.rs();
            }

            veccontr.add(tmp_sc_v);
#else
#if 0
            tmp_matloc_aniso.prod(tmp_sc_aniso,Cp_bis,1,3,2,4);
            matlocf.axpy(tmp_matloc_aniso,delta_aniso);
            veccontr.axpy(tmp_sc_v_aniso,
                          -shocap_aniso*delta_aniso*wpgdet);
#else
            tmp_matloc_aniso.prod(tmp_sc_aniso,Cp_bis_old,1,3,2,4);
            matlocf.axpy(tmp_matloc_aniso,delta_aniso_old);
            veccontr.axpy(tmp_sc_v_aniso,
                          -shocap_aniso*delta_aniso_old*wpgdet);
#endif
#endif
          }

          // This computes either the standard `P_supg' perturbation
          // function or other written by the user in the
          // flux-function.
          if (!use_Ajac_old)
            adv_diff_ff->comp_P_supg(P_supg);
	  
          for (int jel=1; jel<=nel; jel++) {
            P_supg.ir(1,jel);
	  
            veccontr.ir(1,jel);
            matlocf.ir(1,jel);
	  
            if (use_Ajac_old) {
              tmp4.prod(tmp1_old,P_supg,-1,1,-1);
	    } else {
	      tmp4.prod(tmp1,P_supg,-1,1,-1);
	    }
	    
	    // veccontr.ir(1,jel).axpy(tmp4,wpgdet);
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
	      
	      tmp21.set(SHAPE).scale(wpgdet*rec_Dt_m);
	      adv_diff_ff->enthalpy_fun->comp_P_Cp(P_Cp,P_supg);
	      tmp22.prod(P_Cp,tmp21,1,3,2);
	      matlocf.add(tmp22);
	      
	      if (ALE_flag) {
#if 0
		tmp_ALE_07.prod(P_supg,tmp_ALE_01,1,3,2);
		matlocf.axpy(tmp_ALE_07,-wpgdet);
		tmp_ALE_06.prod(P_supg,tmp_ALE_02,1,-1,-1);
		veccontr.axpy(tmp_ALE_06,wpgdet);
#else
		tmp_ALE_07.prod(P_Cp,tmp_ALE_01,1,3,2);
		matlocf.axpy(tmp_ALE_07,-wpgdet);
		tmp_ALE_06.prod(P_Cp,tmp_ALE_02,1,-1,-1);
		veccontr.axpy(tmp_ALE_06,wpgdet);
#endif
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
        //	volume_flag=0;
      
        lmass.axpy(SHAPE,wpgdet);
      
      }  // ipg loop
      
      volume_flag=0;

      if (comp_res) {

	// MODIF BETO 8/6
	if (lumped_mass) {
	  
	  // lumped terms treatment
	  // temporal and source terms are lumped out of the gauss points loop
	  for (int j=1; j<=nel; j++) {
	    
	    /*
	      lstaten.ir(1,j);
	      Un.set(lstaten);
	      lstaten.rs();
	    */
	    
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
	    
	    adv_diff_ff->get_C(Cr);

            tmp12.set(Cr).add(Cp_bis);
            
            double mass_node;
            mass_node = double(lmass.get(j));
            veccontr.ir(1,j).axpy(tmp10,mass_node).rs();
            matlocf.ir(1,j).ir(3,j).axpy(tmp12,mass_node).rs();
            
          }
        }

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

    } catch (GenericError e) {
      set_error(1);
      return;
    }

  FastMat2::void_cache();
  FastMat2::deactivate_cache();
 } catch (GenericError e) {
  set_error(1);
 }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef __FUNC__
#define __FUNC__ "NewAdvDif::after_assemble"
void NewAdvDif::
after_assemble(const char *jobinfo) {
  if (compute_fd_adv_jacobian && !MY_RANK
      && !strcmp(jobinfo,"comp_res")) {
    assert(SIZE==1); // should code after the MPI_Reduce's
    printf("Flux jacobian comps: total %d, checked %d (%5.2f%%)\n",
	   comp_total,comp_checked,double(comp_checked)/
	   double(comp_total)*100.0);
#undef PRINT
#define PRINT(label,val) \
printf(label ": min %g, max %g\n",val##_min,val##_max)
    PRINT("|A_a|",A_jac_norm);
    PRINT("|A_n|",A_fd_jac_norm);
    PRINT("|A_a - A_n|",A_jac_err_norm);
    PRINT("|A_a - A_n|/|A_n|",A_rel_err);
#undef PRINT
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
double NewAdvDif::volume() const {
  // Check call is valid at this moment;
  assert(volume_flag);
  return Volume;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
const FastMat2 *NewAdvDif::grad_N() const {
  return &dshapex;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void NewAdvDifFF::comp_P_supg(FastMat2 &P_supg) {
  assert(new_adv_dif_elemset);
  const NewAdvDif *e = new_adv_dif_elemset;
  if (e->ff_options & SCALAR_TAU) {
    double tau_supg_d = e->tau_supg.get(1,1);
    P_supg.set(e->Ao_grad_N).scale(tau_supg_d);
  } else {
    P_supg.prod(e->Ao_grad_N,e->tau_supg,1,2,-1,-1,3);
  }
  if (e->use_ALE() || e->ALE_form()) {
    // FMSHV(e->v_mesh);
    tmp_P_supg_ALE_1.prod(*(e->grad_N()),e->v_mesh,-1,1,-1);
    tmp_P_supg_ALE_2.prod(e->Cp,tmp_P_supg_ALE_1,2,3,1);
    // tmp_P_supg_ALE_3.prod(tmp_P_supg_ALE_2,tau_supg_c,1,2,-1,-1,3);
    if (e->ff_options & SCALAR_TAU) {
      double tau_supg_d = e->tau_supg.get(1,1);
      tmp_P_supg_ALE_3.set(tmp_P_supg_ALE_2).scale(tau_supg_d);
    } else {
      tmp_P_supg_ALE_3
        .prod(tmp_P_supg_ALE_2,e->tau_supg,1,2,-1,-1,3);
    }
    P_supg.minus(tmp_P_supg_ALE_3);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void NewAdvDifFF::get_bcconv_factor(FastMat2 &bcconv_factor) {

  int nelprops,nel,ndof;
  elemset->elem_params(nel,ndof,nelprops);

  bcconv_factor.resize(1,ndof);

  bcconv_factor.set(1.);

  int ierr = elemset->get_double("bcconv_factor",
		     *bcconv_factor.storage_begin(),1,ndof);

  /*
  const char *line;
  elemset->get_entry("bcconv_factor",line);
  if(line) {
    read_double_array(bcconv_factor,line);
    PETSCFEM_ASSERT0(bcconv_factor.size() == ndof,
		     "bcconv_factor needs ndof values \n");
  } else bcconv_factor.resize(ndof,1.0);
  */

}

/*
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void NewAdvDifFF::set_profile(FastMat2 &seed) {

  seed.set(1.);

  int nelprops,nel,ndof;
  elemset->elem_params(nel,ndof,nelprops);

  vector<double> bcconv_factor;
  const char *line;
  elemset->get_entry("bcconv_factor",line);
  if(line) {
    read_double_array(bcconv_factor,line);
    PETSCFEM_ASSERT0(bcconv_factor.size() == ndof,
		     "bcconv_factor needs ndof values \n");
  } else bcconv_factor.resize(ndof,1.0);

  for (int j=0; j<ndof; j++) {
	 seed.ir(1,j+1).scale(bcconv_factor[j]).rs();
	  }

}

*/

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void NewAdvDifFF::set_profile(FastMat2 &seed) {
  seed.set(1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void NewAdvDifFF::Riemann_Inv(const FastMat2 &U, const FastMat2 &normal,
			      FastMat2 &Rie, FastMat2 &drdU, FastMat2 &C_){
  PETSCFEM_ERROR0("Not implemented Riemman_Inv() method\n"
		  "for this flux function\n");
#undef SHAPE_LOW
#undef DSHAPEXI_LOW
#undef WPG_LOW

}

#undef SHAPE
#undef DSHAPEXI
#undef WPG
