/* $Id: fracstep.cpp,v 1.1 2000/12/28 12:54:43 mstorti Exp $ */
 
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"

#include "fracstep.h"

#define MAXPROP 100


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "fracstep::assemble"
int fracstep::assemble(double *retval,Nodedata *nodedata,double *locst,
		       double *locst2,Dofmap *dofmap,int ijob,
		       char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode) {

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;

//    if (ijob==COMP_MAT) {
//      printf("This elemset routine doesn't computes matrices\n");
//      exit(1);
//    }

//    if (!(ijob==COMP_VEC || ijob==COMP_FDJ)) {
//      printf("Don't know how to compute ijob: %d\n",ijob);
//      exit(1);
//    }
    
  //#define NODEDATA(j,k) (nodedata[ndim*(j)+(k)])
  //#define NODEDATA(j,k) VEC2(nodedata,j,k,ndim)
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define IDENT(j,k) (ident[ndof*(j)+(k)]) 
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)
  
  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim,couple_velocity=0;
  ierr = get_int(thash,"couple_velocity",&couple_velocity,1); CHKERRA(ierr);
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  int nen = nel*ndof;

  // Unpack Dofmap
  int *ident,neq,nnod;
  //  double *load;
  ident = dofmap->ident;
  neq = dofmap->neq;
  nnod = dofmap->nnod;
  //  load = dofmap->load;
  // fixme:= esto por ahora lo dejo asi. Que eventualmente pueda
  // ser diferente el ndof de 
  //  ndof = dofmap->ndof;

  // Unpack nodedata
  int nu=nodedata->nu;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  int comp_res_mom  = !strcmp(jobinfo,"comp_res_mom");
  int comp_mat_mom  = !strcmp(jobinfo,"comp_mat_mom");
  int comp_res_poi  = !strcmp(jobinfo,"comp_res_poi");
  int comp_mat_poi  = !strcmp(jobinfo,"comp_mat_poi");
  int comp_res_prj  = !strcmp(jobinfo,"comp_res_prj");
  int comp_mat_prj  = !strcmp(jobinfo,"comp_mat_prj");
  int comp_mat_mom_prof  = !strcmp(jobinfo,"comp_mat_mom_prof");
  int comp_mat_poi_prof  = !strcmp(jobinfo,"comp_mat_poi_prof");

  // allocate local vecs
  int kdof;
  //jdofloc = new int[nel*ndof];
  //      Matrix matloc(nen,nen), xloc(nel,ndim), veccontr(nel,ndof),
  //        locstate(nel,ndof);
  Matrix veccontr(nel,ndof),xloc(nel,ndim),
    locstate(nel,ndof),locstate2(nel,ndof),tmp(nel,ndof),
    ustate2(nel,ndim);

  if (ndof != ndim+1) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  Matrix matloc(nen,nen), matlocmom(nel,nel), masspg(nel,nel),
    matlocmom2(nen,nen),grad_u_ext(ndof,ndof);
  grad_u_ext = 0;
  Matrix seed;
  if (comp_mat_mom || comp_mat_mom_prof || comp_mat_prj) {
    seed= Matrix(ndof,ndof);
    seed=0;
    for (int j=1; j<=ndim; j++) {
      seed(j,j)=1;
    }
  } else if (comp_mat_poi) {
    seed= Matrix(ndof,ndof);
    seed=0;
    seed(ndof,ndof)=1;
  }
    

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  double Dt,alpha=0.5,alphap=1;
  ierr = get_double(thash,"Dt",&Dt); CHKERRA(ierr);
  ierr = get_double(thash,"alpha",&alpha,1); CHKERRA(ierr);
  ierr = get_double(thash,"alpha_presion",&alphap,1); CHKERRA(ierr);
  int weak_poisson = 1;
  ierr = get_int(thash,"weak_poisson",&weak_poisson,1); CHKERRA(ierr);


  // Factor para la estabilizacion
  double taufac=1.0;


  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))

  int nprops=iprop;
  
  double rho=1.;

  // Gauss Point data
  char *geom;
  thash->get_entry("geometry",geom);
  GPdata gp_data(geom,ndim,nel,npg);

  // Definiciones para descargar el lazo interno
  Matrix grad_fi,Uintri;
  RowVector Pert, W;
  //  LogAndSign L;
  double detJaco, UU, u2, Peclet, psi, tau, div_u_star,
    wpgdet;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  Matrix dshapex(ndim,nel),Jaco(ndim,ndim),iJaco(ndim,ndim),
    grad_u(ndim,ndim),grad_u_star(ndim,ndim),dshapext(nel,ndim),resmom(nel,ndim);
  RowVector fi(ndof);
  ColumnVector grad_p(ndim);
  ColumnVector u(ndim),u_star(ndim),uintri(ndim),rescont(nel);
  
  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    //if (epart[k] != myrank+1) continue;
    ielh++;
    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    //    printf("element %d, prop %f\n",k,ELEMPROPS(k,0));
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc(kloc+1,jdim+1) = NODEDATA(node-1,jdim);
      }
    }
    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    locstate << &(LOCST(ielh,0,0));
    locstate2 << &(LOCST2(ielh,0,0));
    matlocmom = 0;
    matlocmom2 = 0;
    veccontr = 0;
    resmom = 0;
    rescont = 0;

#define DSHAPEXI (gp_data.dshapexi[ipg])
#define SHAPE    (gp_data.shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      //      Matrix xpg = SHAPE * xloc;
      Jaco = DSHAPEXI * xloc;

      elem = k;

      detJaco = mydet(Jaco);
      if (detJaco <= 0.) {
	cout << "Jacobian of element " << elem << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	exit(1);
      }
      wpgdet = detJaco*WPG;
      iJaco = Jaco.i();
      dshapex = iJaco * DSHAPEXI;
      // vector de estado componentes de velocidad
      
      if (comp_res_mom) {
	// state variables and gradient
	
	ustate2 = locstate2.Columns(1,ndim);

	u = (SHAPE * ustate2).t();
	u_star = (SHAPE * locstate.Columns(1,ndim)).t();

	grad_u = dshapex * ustate2;
	grad_u_star = dshapex * locstate.Columns(1,ndim);
	grad_p = dshapex * locstate2.Column(ndim+1);

	double u2 = u.SumSquare();
	uintri = iJaco * u;
	double Uh = sqrt(uintri.SumSquare())/2.;

	if(u2<=1e-6*(2. * Uh * VISC)) {
	  Peclet=0.;
	  psi=0.;
	  tau=0.;
	} else {
	  Peclet = u2 / (2. * Uh * VISC);
	  psi = 1./tanh(Peclet)-1/Peclet;
	  tau = psi/(2.*Uh);
	}
	Pert  = (taufac * tau) * u.t() * dshapex;
	W = SHAPE + Pert;

	// Termino - (u.grad) u - 1/rho grad p (Sin debilitar) 

	// explicit version
	//resmom += wpgdet * W.t() * (- u.t() * grad_u - (1/rho)* grad_p.t());

	// implicit version - Backward Euler
	//	resmom += wpgdet * W.t() * (- u.t() * grad_u_star
	//		    - (1/rho)* grad_p.t());

	// implicit version - Jack Nicholson
	// resmom += wpgdet * W.t() * (- 0.5 * u.t() * (grad_u + grad_u_star)
	//		    - (1/rho)* grad_p.t());

	// implicit version - Crank-Nicholson !!!! := debug
	resmom += wpgdet * W.t() * 
	  (- ((1-alpha) * u.t() * grad_u
	      + alpha * u_star.t() * grad_u_star)
	   - ((1-alphap)/rho)* grad_p.t());
	// sacamos gradiente de presion en la ec. de momento (conf. Codina)
	//- (1/rho)* grad_p.t());

	// Parte difusiva
	// version implicita
	// resmom -= wpgdet * VISC * dshapex.t() * grad_u_star;
	// version Crank-Nicholson 
	resmom -= wpgdet * VISC * dshapex.t() *
	  ((1-alpha) * grad_u+ alpha * grad_u_star);
	
	// Parte temporal
	resmom -= (wpgdet/Dt) * W.t() * (u_star-u).t();
       
      } else if (comp_mat_mom_prof) {

	masspg=1;
	grad_u_ext=0;
	if (couple_velocity) {
	  grad_u_ext.SubMatrix(1,ndim,1,ndim) = 1;
	} else {
	  for (jdim=1; jdim<=ndim; jdim++) {
	    grad_u_ext(jdim,jdim)=1;
	  }
	}
	matlocmom2 += kron(masspg,grad_u_ext.t());

      } else if (comp_mat_poi_prof) {

	masspg=1;
	grad_u_ext=0;
	grad_u_ext(ndim+1,ndim+1) = 1;
	matlocmom2 += kron(masspg,grad_u_ext);

      } else if (comp_mat_mom) {
	// state variables and gradient
	
	ustate2 = locstate2.Columns(1,ndim);

	u = (SHAPE * ustate2).t();
	u_star = (SHAPE * locstate.Columns(1,ndim)).t();

	grad_u = dshapex * ustate2;
	grad_u_star = dshapex * locstate.Columns(1,ndim);
	grad_p = dshapex * locstate2.Column(ndim+1);

	double u2 = u.SumSquare();
	uintri = iJaco * u;
	double Uh = sqrt(uintri.SumSquare())/2.;

	if(u2<=1e-6*(2. * Uh * VISC)) {
	  Peclet=0.;
	  psi=0.;
	  tau=0.;
	} else {
	  Peclet = u2 / (2. * Uh * VISC);
	  psi = 1./tanh(Peclet)-1/Peclet;
	  tau = psi/(2.*Uh);
	}
	Pert  = (taufac * tau) * u.t() * dshapex;
	W = SHAPE + Pert;

	// Termino - (u.grad) u - 1/rho grad p (Sin debilitar) 
	//	resmom += wpgdet * W.t() * (- u.t() * grad_u
	//                           - (1/rho)* grad_p.t());

	// implicit version - Backward Euler
	// matlocmom += wpgdet * W.t() * u.t() * dshapex;

	// implicit version - Jack - Nicholson
	// matlocmom += ( 0.5 * wpgdet) * W.t() * u.t() * dshapex;

	// implicit version - Crank - Nickolson
	matlocmom += ( alpha * wpgdet) * W.t() * u_star.t() * dshapex;
	masspg = W.t() * SHAPE;
	grad_u_ext.SubMatrix(1,ndim,1,ndim) = grad_u_star;
	if (couple_velocity)
	  matlocmom2 += ( alpha * wpgdet) * kron(masspg,grad_u_ext.t());
//  	cout << "masspg:" << endl << masspg << endl;
//  	cout << "grad_u_ext:" << endl << grad_u_ext << endl;

	// Parte difusiva
	matlocmom += (alpha *wpgdet) * VISC * dshapex.t() * dshapex ; // grad_u_star;
	
	// Parte temporal
	matlocmom += (wpgdet/Dt) * W.t() * SHAPE;
       
      } else if (comp_res_poi) {

	//double h_max = 2*Jaco.Norm1();
	//double Dt_crit = h_max*h_max/VISC;

	//double Dt_art= ( Dt<Dt_crit ? Dt_crit : Dt);
	double Dt_art=Dt;

	// printf("Dt %f, h_max %f, Dt_crit %f, Dt_art %f\n",
	// Dt, h_max, Dt_crit, Dt_art);

	ustate2 = locstate2.Columns(1,ndim);
	grad_p = dshapex * locstate.Column(ndim+1);
	// fixme:= Debilito el termino de la divergencia en el miembro
	// derecho de Poisson para que me de un r.h.s. con suma nula y
	// asi el Poisson sea independiente del nodo que se fija. 

	// version debilitada
	if (weak_poisson) {
	  u_star = (SHAPE * ustate2).t();
	  rescont -= wpgdet * (-(rho/Dt_art) * dshapex.t() * u_star
			       + dshapex.t() * grad_p);
	} else {
	  // version no debilitada
	  div_u_star = (dshapex * ustate2).Trace();
	  rescont -= wpgdet * ((rho/Dt_art) * SHAPE.t() * div_u_star
			       + dshapex.t() * grad_p);
	}

      } else if (comp_mat_poi) {

	double Dt_art=Dt;
	matlocmom += wpgdet * dshapex.t() * dshapex; 
 
      } else if (comp_res_prj) {

	grad_p = dshapex * locstate2.Column(ndim+1);

	u_star = (SHAPE * locstate2.Columns(1,ndim)).t();
	u = (SHAPE * locstate.Columns(1,ndim)).t();

	resmom += wpgdet * SHAPE.t() *
	  (-(alphap/rho) *grad_p - (u - u_star)/Dt).t();

      } else if (comp_mat_prj) {

	// fixme:= esto me parece que deberia ir con signo - !!
	matlocmom += wpgdet/Dt * SHAPE.t() * SHAPE ;

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

      }

    }
    if (comp_res_mom || comp_res_prj) {
      veccontr.Columns(1,ndim) = resmom;
    } else if (comp_res_poi) {
      veccontr.Column(ndim+1) = rescont;
    } else if (comp_mat_poi || comp_mat_prj ) {
      matloc = kron(matlocmom,seed);
    } else if (comp_mat_mom || comp_mat_mom_prof ) {
      matloc = kron(matlocmom,seed) + matlocmom2;
    } else if (comp_mat_poi_prof ) {
      matloc = matlocmom2;
    }

    
    if (ijob==COMP_VEC || ijob==COMP_FDJ || ijob==COMP_FDJ_PROF)
      veccontr >> &(RETVAL(ielh,0,0));
    if (ijob==COMP_MAT || ijob==COMP_MAT_PROF)
      matloc >> &(RETVALMAT(ielh,0,0,0,0));

  }
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      

  /*
    # Local Variables: $
    # mode: c++ $
    # End: $
  */

