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
  
#include <vector>

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"
#include "../../src/fastmat2.h"

#include "advective.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int AdvDif::ask(char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "advective::assemble"
int AdvDif::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			Dofmap *dofmap,char *jobinfo,int myrank,
			int el_start,int el_last,int iter_mode,
			const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCSTO(iele,j,k) VEC3(locsto,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)
#define RETVALMATT(iele,j,k,p,q) VEC5(retvalt,iele,j,nel,k,ndof,p,nel,q,ndof)
#define RETVAL_FDJ(iele,j,k) VEC3(retval_fdj,iele,j,nel,k,ndof)

  int ierr=0;

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)
  
  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
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

  double *locst,*locsto,*retval,*retvalt;
#ifdef CHECK_JAC
  double *retval_fdj;
#endif
  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step. 
  vector<double> *dtmin;
  double lambda_max;
  int jdtmin;
  GlobParam *glob_param;
  // The trapezoidal rule integration parameter 
#define ALPHA (glob_param->alpha)
#define DT (glob_param->Dt)
  
  if (comp_res) {
    int j=-1;
    locsto = arg_data_v[++j].locst;
    locst = arg_data_v[++j].locst;
    retval = arg_data_v[++j].retval;
    jdtmin = ++j;
#define DTMIN ((*(arg_data_v[jdtmin].vector_assoc))[0])
#define WAS_SET arg_data_v[jdtmin].was_set
    if (comp_mat_each_time_step_g) retvalt = arg_data_v[++j].retval;
    glob_param = (GlobParam *)arg_data_v[++j].user_data;;
#ifdef CHECK_JAC
    retval_fdj = arg_data_v[++j].retval;
#endif
  }

  FastMat2 matlocf(4,nel,ndof,nel,ndof);
  if (comp_prof) {
    retvalt = arg_data_v[0].retval;
    matlocf.set(1.);
  }

  //o Use the weak form for the Galerkin part of the advective term. 
  SGETOPTDEF(int,weak_form,1);
  //o Use lumped mass.
  SGETOPTDEF(int,lumped_mass,1);
  //o Parameter to control the amount of SUPG perturbation 
  //     added to the mass matrix to be consistent SUPG
  //     \verb+beta_supg+=0 implies consistent Galerkin and
  //     \verb+beta_supg+=1 implies full consistent SUPG. 
  SGETOPTDEF(double,beta_supg,0.8);
 
  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof), 
    locstateo(nel,ndof),locstaten(nel,ndof),
    matloc,eye_ndof(ndof,ndof);

  eye_ndof.set(0.).eye(1.);
//    if (ndof != ndim+1) {
//      PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
//    }
  
  nen = nel*ndof;

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  // char *geom;
  // thash->get_entry("geometry",geom);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);
  // GPdata gp_data(geom,ndim,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco, wpgdet, delta_sc;
  int elem, ipg,node, jdim, kloc,lloc,ldof,ret_options;

  FMatrix dshapex(ndim,nel),Jaco(ndim,ndim),Jaco_av(ndim,ndim),iJaco(ndim,ndim),
    flux(ndof,ndim),fluxd(ndof,ndim),mass(nel,nel),
    grad_U(ndim,ndof), P_supg(ndof,ndof), A_grad_U(ndof),
    G_source(ndof), dUdt(ndof), Uo(ndof),Un(ndof), tau_supg(ndof,ndof);
  // These are edclared but not used
  FMatrix nor,lambda,Vr,Vr_inv,U(ndof),lmass(nel),Id_ndof(ndof,ndof),
    tmp1,tmp2,tmp3,tmp4,tmp5,hvec(ndim),tmp6,tmp7,
    tmp8,tmp9,tmp10,tmp11(ndof,ndim),tmp12,tmp13,tmp14,
    tmp15,tmp16,tmp17,tmp18,tmp19,tmp20,tmp21,tmp22,tmp23;
  FastMat2 A_jac(3,ndim,ndof,ndof),D_jac(4,ndim,ndim,ndof,ndof),
    A_grad_N(3,nel,ndof,ndof);

  Id_ndof.set(0.);
  for (int j=1; j<=ndof; j++) Id_ndof.setel(1.,j,j);

#ifdef USE_FASTMAT2_CACHE
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
#endif

  int ielh=-1;
  int start_chunk=1;
  ElementList elemlist(this,el_start,el_last+1,DO_NOT_INCLUDE_GHOST);
  for (ElementIterator element = elemlist.begin(); 
       element!=elemlist.end(); element++) {
    // if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    int k,ielh;
    element.position(k,ielh);

#ifdef USE_FASTMAT2_CACHE
    FastMat2::reset_cache();
#endif
    //    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    //    printf("element %d, prop %f\n",k,ELEMPROPS(k,0));
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) 
	xloc.setel(NODEDATA(node-1,jdim),kloc+1,jdim+1);
      for (int ih=1; ih<=nH; ih++) 
	Hloc.setel(NODEDATA(node-1,ndim+ih-1),kloc+1,ih);
    }

    if (comp_prof) {
      matlocf.export(&(RETVALMATT(ielh,0,0,0,0)));
      continue;
    }

    if (comp_res) {
      lambda_max=0;
      locstateo.set(&(LOCSTO(ielh,0,0))); // State at time t_n
      locstaten.set(&(LOCST(ielh,0,0))); // State at time t_{n+1}
    }
    
    // State at time t_{n+\alpha}
    locstate.set(0.).axpy(locstaten,ALPHA).axpy(locstateo,(1-ALPHA));

    veccontr.set(0.);
    mass.set(0.);
    lmass.set(0.);
    matlocf.set(0.);

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

	Un.prod(SHAPE,locstaten,-1,-1,1);
	Uo.prod(SHAPE,locstateo,-1,-1,1);
	U.set(0.).axpy(Uo,1-ALPHA).axpy(Un,ALPHA);
	dUdt.set(Un).rest(Uo).scale(1./DT);

	grad_U.prod(dshapex,locstate,1,-1,-1,2) ;

	delta_sc=0;
	double lambda_max_pg;

	ierr =  (*flux_fun)(U,ndim,iJaco,H,grad_H,flux,fluxd,A_jac,
			    A_grad_U,grad_U,G_source,D_jac,tau_supg,delta_sc,
			    lambda_max_pg,
			    thash,nor,lambda,Vr,Vr_inv,
			    &(ELEMPROPS(k,0)),NULL,COMP_SOURCE |
			    COMP_UPWIND,start_chunk,ret_options);

	if (lambda_max_pg>lambda_max) lambda_max=lambda_max_pg;

	tmp10.set(G_source).rest(dUdt);	// tmp10 = G - dUdt
	tmp1.rs().set(tmp10).rest(A_grad_U); //tmp1= G - dUdt - A_grad_U

	tmp15.set(SHAPE).scale(wpgdet/DT);
	tmp12.prod(SHAPE,tmp15,1,2); // tmp12 = SHAPE' * SHAPE
	tmp13.prod(tmp12,eye_ndof,1,3,2,4); // tmp13 = SHAPE' * SHAPE * I
	matlocf.add(tmp13);

	A_grad_N.prod(dshapex,A_jac,-1,1,-1,2,3);

	// Termino Galerkin
	if (weak_form) {
	  // weak version
	  tmp11.set(flux).rest(fluxd); // tmp11 = flux_c - flux_d

	  tmp23.set(SHAPE).scale(-wpgdet*ALPHA);
	  // tmp14.prod(A_grad_N,tmp23,1,2,4,3); // debug:=
	  // matlocf.add(tmp14);
	  tmp14.prod(A_grad_N,tmp23,3,2,4,1);
	  matlocf.rest(tmp14);
	} else {
	  // tmp2.prod(SHAPE,tmp1,1,2); // tmp2= SHAPE' * (G - dUdt - A_grad_U)
	  tmp2.prod(SHAPE,A_grad_U,1,2); // tmp2= SHAPE' * A_grad_U
	  veccontr.axpy(tmp2,-wpgdet);
	  tmp11.set(0.).rest(fluxd); // tmp11 = - flux_d

	  tmp23.set(SHAPE).scale(wpgdet*ALPHA);
	  // tmp14.prod(A_grad_N,tmp23,3,2,4,1);  // debug:= 
	  // matlocf.add(tmp14);
	  tmp14.prod(A_grad_N,tmp23,1,2,4,3);
	  matlocf.rest(tmp14);

	}
	// tmp8= DSHAPEX * (w*flux_c - flux_d)
	//            w = weak_form
	tmp8.prod(dshapex,tmp11,-1,1,2,-1); 
	tmp9.prod(SHAPE,tmp10,1,2); // tmp9 = SHAPE' * (G - dUdt)
	tmp8.add(tmp9);		// tmp8 = DSHAPEX * tmp11
	veccontr.axpy(tmp8,wpgdet);

	// Diffusive term in matrix
	tmp17.set(dshapex).scale(wpgdet*ALPHA);
	tmp16.prod(D_jac,tmp17,1,-1,2,3,-1,4);
	tmp18.prod(tmp16,dshapex,-1,2,4,1,-1,3);
	matlocf.add(tmp18);

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

	  tmp21.set(SHAPE).scale(wpgdet/DT);
	  tmp22.prod(P_supg,tmp21,1,3,2);
	  matlocf.add(tmp22);
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

      veccontr.export(&(RETVAL(ielh,0,0)));
#ifdef CHECK_JAC
      veccontr.export(&(RETVAL_FDJ(ielh,0,0)));
#endif
      if (comp_mat_each_time_step_g) 
	matlocf.export(&(RETVALMATT(ielh,0,0,0,0)));
    } else if (comp_prof) {
      matlocf.export(&(RETVALMATT(ielh,0,0,0,0)));
    }

  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();

}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      

