/*
  This file belongs to he PETSc - FEM package a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
extern int MY_RANK,SIZE;
  
#include <vector>
#include <string>

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"
#include "../../src/fastmat2.h"

#include "nwadvdif.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int NewAdvDif::ask(const char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Returns the list of variables that are 
    logarithmic transformed.
    @param nlog_vars (input) the number of logarithmic variables
    @param log_vars (input) the list of logarithmic variables
    @author M. Storti
*/ 
#undef __FUNC__
#define __FUNC__ "void AdvDifFF::get_log_vars(int,const int*)"
void NewAdvDifFF::get_log_vars(int &nlog_vars,const int *& log_vars) {
  const char *log_vars_entry;
  const int &ndof=elemset->ndof;
  elemset->get_entry("log_vars_list",log_vars_entry); 
  VOID_IT(log_vars_v);
  string s;
  if (log_vars_entry) {
    s=string(log_vars_entry);	// Save local copy
    read_int_array(log_vars_v,log_vars_entry); 
  }
  nlog_vars=log_vars_v.size();
  log_vars = log_vars_v.begin();
  int ierr=0;
  for (int j=0; j<nlog_vars; j++) {
    if (log_vars_v[j]<=0) {
      PetscPrintf(PETSC_COMM_WORLD,"Non positive dof in "
		  "\"log_vars_list\" entry: dof %d\n",
		  log_vars_v[j]);
      ierr=1;
    } else if (log_vars_v[j]>ndof) {
      PetscPrintf(PETSC_COMM_WORLD,"Dof grater that ndof in "
		  "\"log_vars_list\" entry: dof %d, ndof %d\n",
		  log_vars_v[j]);
      ierr=1;
    }
    if (ierr) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Errors while reading \"log_vars_list\"\n");
      if (log_vars_entry)
	PetscPrintf(PETSC_COMM_WORLD,
		    "In line \"%s\"\n",s.c_str());
      exit(1);
    }
  }
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Transforms state vector from logarithmic. The indices of fields
    logarithmically tranformed are listed in \verb+log_vars+. 
    @author M. Storti
    CORREGIR:=
    @param true_lstate (output) Transformed from logarithm
    to positive variable. 
    @param lstate (ouput) input state logarithmically transformed
    (only those fields in \verb+log_vars+).
    @param nlog_vars (input) number of fields logarithmically
    transformed
    @param log_vars (input) list of fields logarithmically
    transformed.
*/ 
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "advective::assemble"
void NewAdvDif::new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
			     const Dofmap *dofmap,const char *jobinfo,
			     const ElementList &elemlist,
			     const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

  int ierr=0;

  int locdof,kldof,lldof;

  NSGETOPTDEF(int,npg,0); //nd
  NSGETOPTDEF_ND(int,ndim,0); //nd
  assert(npg>0);
  assert(ndim>0);
  
  int nelprops;
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
#ifdef CHECK_JAC
    fdj_jac = &arg_data_v[++j];
#endif
  }

  FastMat2 matlocf(4,nel,ndof,nel,ndof),matlocf_mass(4,nel,ndof,nel,ndof);
  if (comp_prof) {
    jac_prof = &arg_data_v[0];
    matlocf.set(1.);
  }

  //o Use the weak form for the Galerkin part of the advective term. 
  NSGETOPTDEF(int,weak_form,1);
  //o Weights the temporal term with $N+\beta P$, i.e.
  // $\beta=0$ is equivalent to waight the temporal term a la 
  // Galerkin and $\beta=1$ is equivalent to do the consistent SUPG weighting.
  NSGETOPTDEF(double,beta_supg,1.);
  //o Use lumped mass (used mainly to avoid oscillations for small time steps).
  NSGETOPTDEF(int,lumped_mass,0);

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
    dUloc(nel,ndof),matloc,eye_ndof(ndof,ndof);
  FastMat2 true_lstate(2,nel,ndof), 
    true_lstateo(2,nel,ndof),true_lstaten(2,nel,ndof);

  eye_ndof.set(0.).eye(1.);
  
  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  NGETOPTDEF_S(string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  double detJaco, wpgdet, delta_sc;
  int elem, ipg,node, jdim, kloc,lloc,ldof,ret_options;

  FMatrix dshapex(ndim,nel),Jaco(ndim,ndim),Jaco_av(ndim,ndim),iJaco(ndim,ndim),
    flux(ndof,ndim),fluxd(ndof,ndim),mass(nel,nel),
    grad_U(ndim,ndof), P_supg(ndof,ndof), A_grad_U(ndof),
    G_source(ndof), dUdt(ndof), Uo(ndof),Un(ndof), tau_supg(ndof,ndof);
  // These are edclared but not used
  FMatrix nor,lambda,Vr,Vr_inv,U(ndof),Ualpha(ndof),
    lmass(nel),Id_ndof(ndof,ndof),
    tmp1,tmp2,tmp3,tmp4,tmp5,hvec(ndim),tmp6,tmp7,
    tmp8,tmp9,tmp10,tmp11(ndof,ndim),tmp12,tmp13,tmp14,
    tmp15,tmp17,tmp19,
    tmp20,tmp21,tmp22,tmp23,
    tmp24;
  FastMat2 A_grad_N(3,nel,ndof,ndof),
    grad_N_D_grad_N(4,nel,ndof,nel,ndof),N_N_C(4,nel,ndof,nel,ndof),
    N_P_C(3,ndof,nel,ndof);

  Id_ndof.set(0.);
  for (int j=1; j<=ndof; j++) Id_ndof.setel(1.,j,j);

  // Initialize flux functions
  adv_diff_ff->start_chunk(ret_options); 

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int start_chunk=1;
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

      //      Matrix xpg = SHAPE * xloc;
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      Jaco_av.add(Jaco);

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
	Un.prod(SHAPE,lstaten,-1,-1,1);
	Uo.prod(SHAPE,lstateo,-1,-1,1);
	Ualpha.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
	dUdt.set(Un).rest(Uo).scale(1./DT);
	for (int k=0; k<nlog_vars; k++) {
	  int jdof=log_vars[k];
	  double UU=exp(Ualpha.get(jdof));
	  dUdt.ir(1,jdof).scale(UU);
	}
	dUdt.rs();

	// Pass to the flux function the true positive values
	U.prod(SHAPE,true_lstate,-1,-1,1);
	grad_U.prod(dshapex,true_lstate,1,-1,-1,2) ;

	delta_sc=0;
	double lambda_max_pg;

	adv_diff_ff->compute_flux(U,iJaco,H,grad_H,flux,fluxd,
				  A_grad_U,grad_U,G_source,
				  tau_supg,delta_sc,
				  lambda_max_pg, nor,lambda,Vr,Vr_inv,
				  COMP_SOURCE | COMP_UPWIND);

	if (lambda_max_pg>lambda_max) lambda_max=lambda_max_pg;

	tmp10.set(G_source);	// tmp10 = G - dUdt
	if (!lumped_mass) tmp10.rest(dUdt);
	if (beta_supg==1.) {
	  tmp1.rs().set(tmp10).rest(A_grad_U); //tmp1= G - dUdt - A_grad_U
	} else {
	  tmp1.set(dUdt).scale(-beta_supg).add(G_source);
	}

	tmp15.set(SHAPE).scale(wpgdet/DT);
	tmp12.prod(SHAPE,tmp15,1,2); // tmp12 = SHAPE' * SHAPE
	tmp13.prod(tmp12,eye_ndof,1,3,2,4); // tmp13 = SHAPE' * SHAPE * I
	if (lumped_mass) {
	  matlocf_mass.add(tmp13);
	} else {
	  matlocf.add(tmp13);
	}

	// A_grad_N.prod(dshapex,A_jac,-1,1,-1,2,3);
	adv_diff_ff->comp_A_grad_N(A_grad_N,dshapex);

	// Termino Galerkin
	if (weak_form) {
	  assert(!lumped_mass && beta_supg==1.); // Not implemented yet!!
	  // weak version
	  tmp11.set(flux).rest(fluxd); // tmp11 = flux_c - flux_d

	  tmp23.set(SHAPE).scale(-wpgdet*ALPHA);
#if 1     // La verdad es que no se cual de estos dos es!!!
	  // Parece que este es el correcto
	  tmp14.prod(A_grad_N,tmp23,1,2,4,3);
	  matlocf.add(tmp14);
#else
	  tmp14.prod(A_grad_N,tmp23,3,2,4,1);
	  matlocf.rest(tmp14);
#endif
	} else {
	  // tmp2.prod(SHAPE,tmp1,1,2); // tmp2= SHAPE' * (G - dUdt - A_grad_U)
	  tmp2.prod(SHAPE,A_grad_U,1,2); // tmp2= SHAPE' * A_grad_U
	  veccontr.axpy(tmp2,-wpgdet);
	  tmp11.set(0.).rest(fluxd); // tmp11 = - flux_d

	  tmp23.set(SHAPE).scale(wpgdet*ALPHA);
#if 1     // La verdad es que no se cual de estos dos es!!!
	  // Parece que este es el correcto
	  tmp14.prod(A_grad_N,tmp23,3,2,4,1); 
	  matlocf.add(tmp14);
#else
	  tmp14.prod(A_grad_N,tmp23,1,2,4,3);
	  matlocf.rest(tmp14);
#endif

	}
	// tmp8= DSHAPEX * (w*flux_c - flux_d)
	//            w = weak_form
	tmp8.prod(dshapex,tmp11,-1,1,2,-1); 
	tmp9.prod(SHAPE,tmp10,1,2); // tmp9 = SHAPE' * (G - dUdt)
	tmp8.add(tmp9);		// tmp8 = DSHAPEX * tmp11
	veccontr.axpy(tmp8,wpgdet);

	// Diffusive term in matrix
#if 0	
	tmp17.set(dshapex).scale(wpgdet*ALPHA);
	adv_diff_ff->comp_D_grad_N(D_grad_N,tmp17);
	tmp18.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
#endif
	adv_diff_ff->comp_grad_N_D_grad_N(grad_N_D_grad_N,
					  dshapex,wpgdet*ALPHA);
	matlocf.add(grad_N_D_grad_N);

        // Reactive term in matrix (Galerkin part)
#if 0
        tmp23.set(SHAPE).scale(wpgdet*ALPHA);
	tmp24.prod(SHAPE,tmp23,1,2); // tmp24 = SHAPE' * SHAPE
	tmp25.prod(tmp24,C_jac,1,3,2,4); // tmp25 = SHAPE' * SHAPE * C_jac
#endif
	adv_diff_ff->comp_N_N_C(N_N_C,SHAPE,wpgdet*ALPHA);
	matlocf.add(N_N_C);

	for (int jel=1; jel<=nel; jel++) {

	  // Should we branch? Not if we take always
	  // the same path...
	  if (ret_options & SCALAR_TAU) {
	    double tau_supg_d = tau_supg.get(1,1);
	    P_supg.set(A_grad_N.ir(1,jel)).scale(tau_supg_d);
	    A_grad_N.rs();
	  } else {
	    P_supg.prod(A_grad_N.ir(1,jel),tau_supg,1,-1,-1,2);
	    A_grad_N.rs();
	    // P_supg = A_grad_N * tau_supg;
	  }
	  
	  tmp4.prod(tmp1,P_supg,-1,1,-1);
	  veccontr.ir(1,jel).axpy(tmp4,wpgdet).ir(1);

	  matlocf.ir(1,jel);
	  tmp19.set(P_supg).scale(ALPHA*wpgdet);
	  tmp20.prod(tmp19,A_grad_N,1,-1,2,-1,3);
	  matlocf.add(tmp20);

          // Reactive term in matrix (SUPG term)
#if 0
          tmp26.set(P_supg).scale(wpgdet*ALPHA);
          tmp27.prod(SHAPE,C_jac,1,2,3); // tmp27 = SHAPE * C_jac 
	  tmp28.prod(tmp26,tmp27,1,-1,2,-1,3); // tmp28 = P_supg * C_jac * SHAPE
#endif
	  adv_diff_ff->comp_N_P_C(N_P_C,P_supg,SHAPE,wpgdet*ALPHA);
	  matlocf.add(N_P_C);

	  tmp21.set(SHAPE).scale(beta_supg*wpgdet/DT);
	  tmp22.prod(P_supg,tmp21,1,3,2);
	  if (lumped_mass) {
	    // I think that if 'lumped_mass' is used then the
	    // contribution to the mass matrix from the SUPG
	    // perturbation term is null, but I include it. 
	    matlocf_mass.ir(1,jel).add(tmp22).rs();
	  } else {
	    matlocf.add(tmp22);
	  }
	  matlocf.rs();
	}

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

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
      // PetscPrintf(PETSC_COMM_WORLD,
      // "On element %d, hloc %f, lambda_max %f, dtloc %f\n",
      // k,hloc,lambda_max,dtloc);
	
      if (dtloc<DTMIN || !WAS_SET) {
	DTMIN = dtloc;
	WAS_SET = 1;
//   	PetscPrintf(PETSC_COMM_WORLD,
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
