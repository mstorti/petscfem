//__INSERT_LICENSE__
//$Id: advective.nm.cpp,v 1.4 2001/05/02 00:08:54 mstorti Exp $

extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
  
#include <vector>

#include <newmat.h>

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"

#include "advective.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int Advective::ask(char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_diag_mat_mass);
   DONT_SKIP_JOBINFO(comp_mat_mass);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "advective::assemble"
int Advective::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			Dofmap *dofmap,char *jobinfo,int myrank,
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
  Matrix  Hloc(nel,nH),H(1,nH),grad_H(ndim,nH);

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

  // weak_form = 1 implies that we use the weak form for the galerkin
  // aprt of the advective term. 
  SGETOPTDEF(int,weak_form,0);
  SGETOPTDEF(int,lumped_mass,1);
  
  // allocate local vecs
  int kdof;
  Matrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof), 
         matloc(nen,nen);

//    if (ndof != ndim+1) {
//      PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
//    }
  
  nen = nel*ndof;

  // Gauss Point data
  char *geom;
  thash->get_entry("geometry",geom);
  GPdata gp_data(geom,ndim,nel,npg);

  // Definiciones para descargar el lazo interno
  double detJaco, wpgdet, delta_sc;
  int elem, ipg,node, jdim, kloc,lloc,ldof,ret_options;

  Matrix dshapex(ndim,nel),Jaco(ndim,ndim),iJaco(ndim,ndim),
    flux(ndof,ndim),mass(nel,nel),
    grad_U(ndim,ndof), P_supg(ndof,ndof), A_grad_U(ndof,1),
    A_grad_N(ndof,ndof), G_source(ndof,1), tau_supg(ndof,ndof);
  // These are edclared but not used
  Matrix nor,lambda,Vr,Vr_inv;
  RowVector U(ndof);
  ColumnVector lmass(nel);

  vector<Matrix *> A_jac;
  for (int jd=1; jd<=ndim; jd++) {
    A_jac.push_back(new Matrix(ndof,ndof));
  }

  DiagonalMatrix Id_ndof(ndof);
  Id_ndof = 1;

  int ielh=-1;
  int start_chunk=1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    ielh++;
    //    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    //    printf("element %d, prop %f\n",k,ELEMPROPS(k,0));
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc(kloc+1,jdim+1) = NODEDATA(node-1,jdim);
      }
      for (int ih=1; ih<=nH; ih++) {
	Hloc(kloc+1,ih) = NODEDATA(node-1,ndim+ih-1);
      }
    }
    if (comp_res) {
      lambda_max=0;
      locstate << &(LOCST(ielh,0,0));
      // cout << "elem: " << k << "\n locstate: \n" << locstate << endl;
    }

    veccontr = 0;
    mass     = 0;
    lmass     = 0;
    matloc=0;

#define DSHAPEXI (gp_data.dshapexi[ipg])
#define SHAPE    (gp_data.shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      //      Matrix xpg = SHAPE * xloc;
      Jaco = DSHAPEXI * xloc;

      detJaco = mydet(Jaco);
      if (detJaco <= 0.) {
	cout << "Jacobian of element " << k << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	exit(1);
      }
      wpgdet = detJaco*WPG;
      iJaco = Jaco.i();
      dshapex = iJaco * DSHAPEXI;

      H = SHAPE * Hloc;
      grad_H = dshapex * Hloc;

      if (comp_res) {
	// state variables and gradient

	U = SHAPE * locstate;
	grad_U = dshapex * locstate ;

	delta_sc=0;
	double lambda_max_pg;
#define TEST_LOOP
#ifdef TEST_LOOP
	Chrono chrono;
	chrono.start();
	printf("start loop...\n");
	int ntimes=2000000;
	for (int kkk=0; kkk<ntimes; kkk++) {
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
	}
	double elaps=chrono.elapsed();
	printf("Elapsed: %f [secs], %d times, rate: %g [evals/sec]\n",
	       elaps,ntimes,ntimes/elaps);
	exit(0);
#endif
	if (lambda_max_pg>lambda_max) lambda_max=lambda_max_pg;

	// Termino Galerkin
	if (weak_form) {
	  // version debilitada
	  veccontr += wpgdet * (dshapex.t() * flux.t() +SHAPE.t() *
				G_source.t());
	} else {
	  // version sin debilitar
	  veccontr += wpgdet * SHAPE.t() * (G_source - A_grad_U).t();

	}

	for (int jel=1; jel<=nel; jel++) {
	  A_grad_N = 0;
	  for (int jd=1; jd<=ndim; jd++) {
	    A_grad_N += dshapex(jd,jel)* AJAC(jd);
	  }
	  // tau_supg may be a 1 x 1 (scalar) or ndof x ndof Matrix. 
	  // if (tau_supg.Nrows() == 1) {
	  if (ret_options & SCALAR_TAU) {
	    double tau_supg_d = tau_supg(1,1);
	    P_supg = tau_supg_d * A_grad_N;
	  } else {
	    P_supg = A_grad_N * tau_supg;
	  }

	  veccontr.Row(jel)
	    += wpgdet * (G_source - A_grad_U).t() * P_supg.t();
	  if (consistent_supg_matrix_g) {
	    int fr,lr,fc,lc;
	    for ( int iel=1 ; iel<= nel ; iel++) {
	      fr=(jel-1)*ndof+1;
	      lr=(jel-1)*ndof+ndof;
	      fc=(iel-1)*ndof+1;
	      lc=(iel-1)*ndof+ndof;
	      matloc.SubMatrix(fr,lr,fc,lc)+=
		wpgdet*SHAPE(iel)*(SHAPE(jel)*Id_ndof+P_supg);
	    }
	  } else if (local_time_step_g) {
	    for (int iel=1; iel<=nel; iel++) 
	      mass(iel,iel) += wpgdet * SHAPE(iel);
	  }
	    
	}
	veccontr +=  - wpgdet * delta_sc * dshapex.t() * grad_U;

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
	    mass(iel,iel) += wpgdet*SHAPE(iel);
	} else {
	  mass += wpgdet * SHAPE.t() * SHAPE;
	}

      } else if (comp_diag_mat_mass) {
	
	lmass += wpgdet * SHAPE.t();

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

      }

    }

    if (comp_mat_mass) {

      matloc = kron(mass,Id_ndof);
      matloc >> &(RETVALMAT(ielh,0,0,0,0));

    } else if (comp_diag_mat_mass) {

      for (int kk=1; kk<=ndof; kk++)
	veccontr.Column(kk) = lmass;
      veccontr >> &(RETVAL(ielh,0,0));

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
      for (int jd=1; jd<=ndim; jd++) {
	hhh = 2.*sqrt(Jaco.Column(jd).SumSquare());
	if (jd==1 || hhh<hloc) hloc=hhh;
      }

      dtloc = hloc/lambda_max;
      // PetscPrintf(PETSC_COMM_WORLD,
      // "On element %d, hloc %f, lambda_max %f, dtloc %f\n",
      // k,hloc,lambda_max,dtloc);
	
      if (dtloc<DTMIN || !WAS_SET) {
	DTMIN = dtloc;
	WAS_SET = 1;
// 	PetscPrintf(PETSC_COMM_WORLD,
// 		    "setting dtmin: %f, hloc %f, lambda_max: %f\n",
// 		    DTMIN,hloc,lambda_max);
      }

      if (local_time_step_g) {
	mass /= dtloc;
	matloc = kron(mass,Id_ndof);
      }

      veccontr >> &(RETVAL(ielh,0,0));
      if (comp_mat_each_time_step_g) 
	matloc >> &(RETVALMATT(ielh,0,0,0,0));
      
    }
  }

  for (int jd=1; jd<=ndim; jd++) delete A_jac[jd-1];
  
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
