//__INSERT_LICENSE__
//$Id: diff.cpp,v 1.5 2002/03/15 12:53:32 mstorti Exp $
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

#include "advective.h"
#include "diff.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int Diff::ask(char *,int &)"
int Diff::ask(const char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
DiffFF::~DiffFF() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "advective::assemble"
void Diff::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			const Dofmap *dofmap,const char *jobinfo,
			const ElementList &elemlist,
			const TimeData *time_data) {

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
  FMatrix  Hloc(nel,nH),H(nH),grad_H(ndim,nH);

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
  double alpha,rec_Dt;

  // The trapezoidal rule integration parameter 
  arg_data *staten_a,*stateo_a,*retval,*fdj_jac,*jac_prof,*Ajac;
  if (comp_res) {
    int j=-1;
    stateo_a = &arg_data_v[++j];
    staten_a = &arg_data_v[++j];
    retval  = &arg_data_v[++j];
    jdtmin = ++j;
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    Ajac = &arg_data_v[++j];
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
#ifdef CHECK_JAC
    fdj_jac = &arg_data_v[++j];
#endif
    // make local copies of global params
    alpha = glob_param->alpha;
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  }

  FastMat2 matlocf(4,nel,ndof,nel,ndof),
    matlocf_mass(4,nel,ndof,nel,ndof);
  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matlocf.set(1.);
  }

  //o Use lumped mass (used mainly to avoid oscillations for small time steps).
  NSGETOPTDEF(int,lumped_mass,0);

  diff_ff->start_chunk(); 

  // lumped_mass:= If this options is activated then all the inertia
  // term matrix comtributions are added to 'matlocf_mass' and the
  // vector contribution terms are discarded. Then at the last moment
  // matlocf_mass*(Un-Uo)/Dt is added.

  // Allocate local vecs
  FMatrix veccontr(nel,ndof),veccontr_mass(nel,ndof),
    xloc(nel,ndim),state(nel,ndof), 
    stateo(nel,ndof),staten(nel,ndof),dUloc_c(nel,ndof),
    dUloc(nel,ndof),matloc;

  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  double detJaco, wpgdet;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FastMat2 dshapex(2,ndim,nel),Jaco(2,ndim,ndim),iJaco(2,ndim,ndim),
    fluxd(2,ndim,ndof),mass(2,nel,nel),Ualpha(1,ndof),U(1,ndof),
    grad_U(2,ndim,ndof),G_source(1,ndof),dUdt(1,ndof), Uo(1,ndof),Un(1,ndof), 
    Ho(1,ndof),Hn(1,ndof),N_Cp_N(4,nel,ndof,nel,ndof);

  // These are edclared but not used
  FMatrix lmass(nel),Id_ndof(ndof,ndof),
    tmp1,tmp2,tmp3,tmp4,tmp5,hvec(ndim),tmp6,tmp7,
    tmp8,tmp9,tmp10,tmp11(ndof,ndim),tmp12,tmp14,
    tmp15,tmp17,tmp19,tmp20,tmp21,tmp22,tmp23,
    tmp24;
  FastMat2 grad_N_D_grad_N(4,nel,ndof,nel,ndof);

  Id_ndof.set(0.);
  for (int j=1; j<=ndof; j++) Id_ndof.setel(1.,j,j);

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  // printf("[%d] %s start: %d last: %d\n",MY_RANK,jobinfo,el_start,el_last);
  for (ElementIterator element = elemlist.begin(); 
       element!=elemlist.end(); element++) {

    FastMat2::reset_cache();

    // Initialize element
    diff_ff->element_hook(element); 
    // Get nodedata info (coords. etc...)
    element.node_data(nodedata,xloc.storage_begin(),
		       Hloc.storage_begin());

    if (comp_prof) {
      matlocf.export_vals(element.ret_mat_values(*jac_prof));
      continue;
    }

    if (comp_res) {
      stateo.set(element.vector_values(*stateo_a));
      staten.set(element.vector_values(*staten_a));
    }

    // State at time t_{n+\alpha}
    state.set(0.).axpy(staten,alpha).axpy(stateo,(1-alpha));

    veccontr.set(0.);
    mass.set(0.);
    lmass.set(0.);
    matlocf.set(0.);
    if (lumped_mass) matlocf_mass.set(0.);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      //      Matrix xpg = SHAPE * xloc;
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);

      detJaco = Jaco.det();
      if (detJaco <= 0.) {
	int k,ielh;
	element.position(k,ielh);
	printf("Jacobian of element %d is negative or null\n"
	       " Jacobian: %f\n",k,detJaco);
	PetscFinalize();
	exit(0);
      }
      wpgdet = detJaco*WPG;
      iJaco.inv(Jaco);
      dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);

      if (nH>0) {
	H.prod(SHAPE,Hloc,-1,-1,1);
	grad_H.prod(dshapex,Hloc,1,-1,-1,2);
      }

      if (comp_res) {
	// state variables and gradient
	Un.prod(SHAPE,staten,-1,-1,1);
	diff_ff->enthalpy(Hn,Un);
	Uo.prod(SHAPE,stateo,-1,-1,1);
	diff_ff->enthalpy(Ho,Uo);

	dUdt.set(Hn).rest(Ho).scale(rec_Dt);
	dUdt.rs();

	// U:= the state at time $t^n+\alpha*\Dt$
	U.prod(SHAPE,state,-1,-1,1);
	grad_U.prod(dshapex,state,1,-1,-1,2) ;

	diff_ff->gp_hook(ipg,U,grad_U);

	diff_ff->compute_flux(U,grad_U,fluxd,G_source,H,grad_H);

	tmp10.set(G_source);	// tmp10 = G - dUdt
	if (!lumped_mass) tmp10.rest(dUdt);

	diff_ff->comp_N_Cp_N(N_Cp_N,SHAPE,wpgdet*rec_Dt);
	if (lumped_mass) {
	  matlocf_mass.add(N_Cp_N);
	} else {
	  matlocf.add(N_Cp_N);
	}

	// Termino Galerkin
	tmp8.prod(dshapex,fluxd,-1,1,-1,2).scale(-1.); 
	tmp9.prod(SHAPE,tmp10,1,2); // tmp9 = SHAPE' * (G - dUdt)
	tmp8.add(tmp9);		// tmp8 = DSHAPEX * tmp11
	veccontr.axpy(tmp8,wpgdet);

	diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
				      dshapex,wpgdet*alpha);
	matlocf.add(grad_N_D_grad_N);

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

      }

    }

    if (comp_res) {

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
	dUloc.set(staten).rest(stateo);
	state.rs();

	// Compute inertia term with lumped mass
	veccontr_mass.prod(matlocf_mass,dUloc,1,2,-1,-2,-1,-2);
	// Add (rest) to vector contribution to be returned
	veccontr.rest(veccontr_mass);
	// Add to matrix contribution to be returned
	matlocf.add(matlocf_mass);
      }

      veccontr.export_vals(element.ret_vector_values(*retval));
#ifdef CHECK_JAC
      veccontr.export_vals(element.ret_fdj_values(*fdj_jac));
#endif
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
