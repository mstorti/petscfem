//__INSERT_LICENSE__
//$Id: advecfm2.cpp,v 1.8 2001/12/20 21:58:55 mstorti Exp $

extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
  
#include <vector>

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>

#include "advective.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int AdvectiveFM2::ask(const char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_diag_mat_mass);
   DONT_SKIP_JOBINFO(comp_mat_mass);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "advective::assemble"
int AdvectiveFM2::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			Dofmap *dofmap,const char *jobinfo,int myrank,
			int el_start,int el_last,int iter_mode,
			const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_mat_mass);
  GET_JOBINFO_FLAG(comp_diag_mat_mass);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)
#define RETVALMATT(iele,j,k,p,q) VEC5(retvalt,iele,j,nel,k,ndof,p,nel,q,ndof)


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

  double *locst,*retval,*retvalt;
  if (comp_mat_mass || comp_diag_mat_mass) {
    retval = arg_data_v[0].retval;
  }
  // lambda_max:= the maximum eigenvalue of the jacobians.
  // used to compute the critical time step. 
  vector<double> *dtmin;
  double lambda_max;
  if (comp_res) {
    locst = arg_data_v[0].locst;
    retval = arg_data_v[1].retval;
#define DTMIN ((*(arg_data_v[2].vector_assoc))[0])
#define WAS_SET arg_data_v[2].was_set
    if (comp_mat_each_time_step_g) retvalt = arg_data_v[3].retval;
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
         matloc(nen,nen);

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
    flux(ndof,ndim),mass(nel,nel),
    grad_U(ndim,ndof), P_supg(ndof,ndof), A_grad_U(ndof),
    G_source(ndof), tau_supg(ndof,ndof);
  // These are edclared but not used
  FMatrix nor,lambda,Vr,Vr_inv,U(ndof),lmass(nel),Id_ndof(ndof,ndof),
    tmp1,tmp2,tmp3,tmp4,tmp5,hvec(ndim),tmp6,tmp7,
    tmp8,tmp9;
  FastMat2 A_jac(3,ndim,ndof,ndof),A_grad_N(3,nel,ndof,ndof);

  Id_ndof.set(0.);
  for (int j=1; j<=ndof; j++) Id_ndof.setel(1.,j,j);

#ifdef USE_FASTMAT2_CACHE
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
#endif

  int ielh=-1;
  int start_chunk=1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
#ifdef USE_FASTMAT2_CACHE
    FastMat2::reset_cache();
#endif
    ielh++;
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

    if (comp_res) {
      lambda_max=0;
      locstate.set(&(LOCST(ielh,0,0)));
    }

    veccontr.set(0.);
    mass.set(0.);
    lmass.set(0.);
    matloc.set(0.);

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

	U.prod(SHAPE,locstate,-1,-1,1);
	grad_U.prod(dshapex,locstate,1,-1,-1,2) ;

	delta_sc=0;
	double lambda_max_pg;
	//#define TEST_LOOP
#ifdef TEST_LOOP
//  	Chrono chrono;
//  	chrono.start();
//  	printf("start loop...\n");
//  	int ntimes=2000000;
//  	for (int kkk=0; kkk<ntimes; kkk++) {
#endif
	ierr =  (*flux_fun)(U,ndim,iJaco,H,grad_H,flux,A_jac,
			    A_grad_U,grad_U,G_source,tau_supg,delta_sc,
			    lambda_max_pg,
			    thash,nor,lambda,Vr,Vr_inv,
			    &(ELEMPROPS(k,0)),NULL,COMP_SOURCE |
			    COMP_UPWIND,start_chunk,ret_options);
#if 0
	SHV(U);
	SHV(flux);
	SHV(A_grad_U);
	SHV(G_source);
	SHV(tau_supg);  
	SHV(lambda);    
	SHV(Vr);        
	SHV(Vr_inv);    
#endif

#ifdef TEST_LOOP
//  	}
//  	double elaps=chrono.elapsed();
//  	printf("Elapsed: %f [secs], %d times, rate: %g [evals/sec]\n",
//  	       elaps,ntimes,ntimes/elaps);
//  	exit(0);
#endif
	if (lambda_max_pg>lambda_max) lambda_max=lambda_max_pg;

	tmp1.rs().set(G_source).rest(A_grad_U);
	// Termino Galerkin
	if (weak_form) {
	  // version debilitada
	  tmp8.prod(dshapex,flux,-1,1,2,-1);
	  tmp9.prod(SHAPE,G_source,1,2);
	  tmp8.add(tmp9);
	  veccontr.axpy(tmp8,wpgdet);
//  	  veccontr += wpgdet * (dshapex.t() * flux.t() +SHAPE.t() *
//  				G_source.t());
	} else {
	  // version sin debilitar
	  tmp2.prod(SHAPE,tmp1,1,2);
	  veccontr.axpy(tmp2,wpgdet);
	}

	A_grad_N.prod(dshapex,A_jac,-1,1,-1,2,3);
	for (int jel=1; jel<=nel; jel++) {

	  // Should we branch? not if we take ever the same path
	  if (ret_options & SCALAR_TAU) {
	    double tau_supg_d = tau_supg.get(1,1);
	    P_supg.set(A_grad_N.ir(1,jel)).scale(tau_supg_d);
	    A_grad_N.rs();
	  } else {
	    assert(0); // not coded yet...
	    // P_supg = A_grad_N * tau_supg;
	  }
	  
	  tmp4.prod(tmp1,P_supg,-1,1,-1);
	  veccontr.ir(1,jel).axpy(tmp4,wpgdet).ir(1);
	  if (consistent_supg_matrix_g) {
	    tmp7.set(Id_ndof).scale(SHAPE.get(jel)).axpy(P_supg,beta_supg);
	    int fr,lr,fc,lc;
	    for ( int iel=1 ; iel<= nel ; iel++) {
	      fr=(jel-1)*ndof+1;
	      lr=(jel-1)*ndof+ndof;
	      fc=(iel-1)*ndof+1;
	      lc=(iel-1)*ndof+ndof;
	      matloc.is(1,fr,lr).is(2,fc,lc).axpy(tmp7,wpgdet*SHAPE.get(iel)).rs();
	    }
	    // FMSHV(matloc);
	    matloc.rs();
	  } else if (local_time_step_g) {
	    // for (int iel=1; iel<=nel; iel++) 
	      mass.addel(wpgdet*SHAPE.get(jel),jel,jel);
	  }
	    
	}

	// Shock capturing term
	FastMat2::branch();
	if (delta_sc!=0.) {
#ifdef USE_FASTMAT2_CACHE
	  FastMat2::choose(0);
#endif
	  tmp5.prod(dshapex,grad_U,-1,1,-1,2);
	  veccontr.axpy(tmp5,-wpgdet*delta_sc);
	}
	FastMat2::leave();

	// This flag is defined in "advective.h"
#if 0
	// Consistent SUPG mass matrix
        mass = wpgdet * SHAPE.t() * SHAPE;
        matloc += kron(mass,Id_ndof);
        for ( int iel=1 ; iel<= nel ; iel++) {
	  mass.Column(iel) = wpgdet * SHAPE(iel);
        }
	matloc += kron(mass,P_supg);
#endif
	// 

      } else if (comp_mat_mass) {
	
	if (lumped_mass) {
	  for (int iel=1; iel<=nel; iel++) 
	    mass.addel(wpgdet*SHAPE.get(iel),iel,iel);
	} else {
	  tmp5.prod(SHAPE,SHAPE,1,2);
	  mass.axpy(tmp5,wpgdet);
	}

      } else if (comp_diag_mat_mass) {
	
	assert(0);
	// not coded yet...
	// lmass += wpgdet * SHAPE.t();

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

      }

    }

    if (comp_mat_mass) {

      matloc.kron(mass,Id_ndof);
      matloc.export_vals(&(RETVALMAT(ielh,0,0,0,0)));

    } else if (comp_diag_mat_mass) {

      assert(0); // not coded yet...
//        for (int kk=1; kk<=ndof; kk++)
//  	veccontr.Column(kk) = lmass;
//        veccontr >> &(RETVAL(ielh,0,0));

    } else if (comp_res) {

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

      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      if (comp_mat_each_time_step_g) 
	matloc.export_vals(&(RETVALMATT(ielh,0,0,0,0)));
#if 0
      FMSHV(veccontr);
      FMSHV(matloc);
      SHV(dtloc);
#endif
      
    }
  }
#ifdef USE_FASTMAT2_CACHE
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
#endif
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
