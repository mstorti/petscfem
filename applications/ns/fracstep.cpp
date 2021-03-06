//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include "nsi_tet.h"
#include "fracstep.h"

#define MAXPROP 100

void nmprint(Matrix &A) { cout << A << endl; }

const double FIX = 0.1;
void fix_null_diagonal_entries(Matrix &A,Matrix &B,int n) {
  B = 0;
  for (int j=1; j<=n; j++) {
    if (A(j,j)==0.) {
      A(j,j) = FIX;
      B(j,j) = FIX;
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int fracstep::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_mat_prof);
  DONT_SKIP_JOBINFO(comp_res_mom);
  DONT_SKIP_JOBINFO(comp_mat_poi);
  DONT_SKIP_JOBINFO(comp_res_poi);
  DONT_SKIP_JOBINFO(comp_mat_prj);
  DONT_SKIP_JOBINFO(comp_res_prj);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "fracstep::assemble"
int fracstep::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		       Dofmap *dofmap,const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_) {

  assert(fractional_step);	// do not mix with the TET (monolithic) solver
  assert(reuse_mat);		// not coded yet Newmat version
				// wirhout reuse mat option
  int ierr=0;

  GET_JOBINFO_FLAG(comp_mat_prof);
  GET_JOBINFO_FLAG(comp_res_mom);
  GET_JOBINFO_FLAG(comp_mat_poi);
  GET_JOBINFO_FLAG(comp_res_poi);
  GET_JOBINFO_FLAG(comp_mat_prj);
  GET_JOBINFO_FLAG(comp_res_prj);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele) VEC2(retval,iele,0,nel*ndof)
#define RETVALMAT(iele)     VEC2(retvalmat,    iele,0,nen*nen)
#define RETVALMAT_MOM(iele) VEC2(retvalmat_mom,iele,0,nen*nen)
#define RETVALMAT_POI(iele) VEC2(retvalmat_poi,iele,0,nen*nen)
#define RETVALMAT_PRJ(iele) VEC2(retvalmat_prj,iele,0,nen*nen)

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

  // Unpack nodedata
  int nu=nodedata->nu;
  int nnod = dofmap->nnod;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  double *locst=NULL,*locst2=NULL,*retval=NULL,
    *retvalmat=NULL,*retvalmat_mom=NULL,*retvalmat_poi=NULL,
    *retvalmat_prj=NULL;

  // rec_Dt is the reciprocal of Dt (i.e. 1/Dt)
  // for steady solutions it is set to 0. (Dt=inf)
  GlobParam *glob_param=NULL;
  double Dt=NAN;
  arg_data *A_mom_arg=NULL,*A_poi_arg=NULL,*A_prj_arg=NULL;
  if (comp_mat_prof) {
    int ja=0;
    retvalmat_mom = arg_data_v[ja++].retval;
    retvalmat_poi = arg_data_v[ja++].retval;
    retvalmat_prj = arg_data_v[ja++].retval;
  } else if (comp_res_mom) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    A_mom_arg = &arg_data_v[ja];
    retvalmat = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else if (comp_mat_poi) {
    int ja=0;
    A_poi_arg = &arg_data_v[ja];
    retvalmat_poi = arg_data_v[ja].retval;
  } else if (comp_res_poi) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else if (comp_mat_prj) {
    int ja=0;
    A_prj_arg = &arg_data_v[ja];
    retvalmat_prj = arg_data_v[ja].retval;
  } else if (comp_res_prj) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else assert(0); // Not implemented yet!!

  Matrix veccontr(nel,ndof),xloc(nel,ndim),
    locstate(nel,ndof),locstate2(nel,ndof),tmp(nel,ndof),
    ustate2(nel,ndim);

  if (ndof != ndim+1) {
    PetscPrintf(PETSCFEM_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  Matrix matloc(nen,nen), matlocmom(nel,nel), masspg(nel,nel),
    matlocmom2(nen,nen), grad_u_ext(ndof,ndof),
    mom_profile(nen,nen), poi_profile(nen,nen),
    mom_mat_fix(nen,nen), poi_mat_fix(nen,nen);
  grad_u_ext = 0;

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  double alpha=0.5,alphap=1;
  ierr = get_double(thash,"alpha",&alpha,1); CHKERRA(ierr);
  ierr = get_double(thash,"alpha_presion",&alphap,1); CHKERRA(ierr);
  int weak_poisson = 1;
  ierr = get_int(thash,"weak_poisson",&weak_poisson,1); CHKERRA(ierr);

  // Factor para la estabilizacion
  double taufac=1.0;

  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))

  int nprops=iprop;
  
  //o Density
  TGETOPTDEF(thash,double,rho,1.);
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  //GPdata gp_data(geom,ndim,nel,npg);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg);

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
  
  masspg=1;
  grad_u_ext=0;
  if (couple_velocity) grad_u_ext.SubMatrix(1,ndim,1,ndim) = 1;
  else for (jdim=1; jdim<=ndim; jdim++) grad_u_ext(jdim,jdim) = 1;
  mom_profile = kron(masspg,grad_u_ext.t());
  fix_null_diagonal_entries(mom_profile,mom_mat_fix,nen);
  
  grad_u_ext=0;
  grad_u_ext(ndim+1,ndim+1) = 1;
  poi_profile = kron(masspg,grad_u_ext);
  fix_null_diagonal_entries(poi_profile,poi_mat_fix,nen);

  if (comp_mat_prof) {
    mom_profile >> arg_data_v[0].profile; // A_mom
    poi_profile >> arg_data_v[1].profile; // A_poi
    mom_profile >> arg_data_v[2].profile;	// A_prj
  } else if (comp_res_mom) {
    mom_profile >> A_mom_arg->profile;
  } else if (comp_mat_poi) {
    poi_profile >> A_poi_arg->profile;
  } else if (comp_res_poi) {
  } else if (comp_mat_prj) {
    mom_profile >> A_prj_arg->profile;
  } else if (comp_res_prj) {
  } else assert(0);

  Matrix seed;
  if (comp_res_mom || comp_mat_prj) {
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

  int ielh=-1;
  int SHV_debug=0;
#undef SHV
#define SHV(pp) { if (SHV_debug) cout << #pp ": " << endl << pp << endl; }
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    //if (epart[k] != myrank+1) continue;
    ielh++;

    if (comp_mat_prof) {
      // We export anything, because in fact `upload_vector'
      // doesn't even look at the `retvalmat' values. It looks at the
      // .profile member in the argument value
      matlocmom >> &(RETVALMAT_MOM(ielh));
      matlocmom >> &(RETVALMAT_PRJ(ielh));
      matlocmom >> &(RETVALMAT_POI(ielh));
      continue;
    } 

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
    if (comp_res_mom || comp_res_poi || comp_res_prj) {
      locstate << &(LOCST(ielh,0,0));
      locstate2 << &(LOCST2(ielh,0,0));
    }
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
      if (detJaco<=0.) {
	detj_error(detJaco,elem);
	set_error(1);
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
	SHV(grad_p);

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

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// RESIDUE CALCULATION
	resmom += wpgdet * W.t() * 
	  (- ((1-alpha) * u.t() * grad_u
	      + alpha * u_star.t() * grad_u_star)
	   - ((1-alphap)/rho)* grad_p.t());
	SHV(resmom);
	// sacamos gradiente de presion en la ec. de momento (conf. Codina)
	//- (1/rho)* grad_p.t());

	// Parte difusiva
	// version implicita
	// resmom -= wpgdet * VISC * dshapex.t() * grad_u_star;
	// version Crank-Nicholson 
	resmom -= wpgdet * VISC * dshapex.t() *
	  ((1-alpha) * grad_u+ alpha * grad_u_star);
	SHV(resmom);
	
	// Parte temporal
	resmom -= (wpgdet/Dt) * W.t() * (u_star-u).t();
	SHV(resmom);

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// JACOBIAN CALCULATION
	matlocmom += ( alpha * wpgdet) * W.t() * u_star.t() * dshapex;
	masspg = W.t() * SHAPE;
	grad_u_ext.SubMatrix(1,ndim,1,ndim) = grad_u_star;
	if (couple_velocity)
	  matlocmom2 += ( alpha * wpgdet) * kron(masspg,grad_u_ext.t());

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

	matlocmom += wpgdet * dshapex.t() * dshapex; 
 
      } else if (comp_res_prj) {

	grad_p = dshapex * locstate2.Column(ndim+1);

	u_star = (SHAPE * locstate2.Columns(1,ndim)).t();
	u = (SHAPE * locstate.Columns(1,ndim)).t();

	resmom += wpgdet * SHAPE.t() *
	  (-Dt*(alphap/rho) *grad_p - (u - u_star)).t();
	SHV(resmom);
	SHV(SHAPE);
	SHV(u);
	SHV(u_star);

      } else if (comp_mat_prj) {

	// fixme:= esto me parece que deberia ir con signo - !!
	matlocmom += wpgdet * SHAPE.t() * SHAPE ;

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

      }

    }
    if (comp_res_mom) {
      veccontr.Columns(1,ndim) = resmom;
      veccontr >> &(RETVAL(ielh));
      matloc = kron(matlocmom,seed) + mom_mat_fix;
      matloc >> &(RETVALMAT(ielh));
    } else if (comp_mat_poi) {
      matloc = kron(matlocmom,seed) + poi_mat_fix;
      matloc >> &(RETVALMAT_POI(ielh));
    } else if (comp_res_poi) {
      veccontr.Column(ndim+1) = rescont;
      veccontr >> &(RETVAL(ielh));
    } else if (comp_mat_prj) {
      matloc = kron(matlocmom,seed) + mom_mat_fix;
      matloc >> &(RETVALMAT_PRJ(ielh));
    } else if (comp_res_prj) {
      veccontr.Columns(1,ndim) = resmom;
      veccontr >> &(RETVAL(ielh));
    } else assert(0);

  }
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      

  /*
    # Local Variables: $
    # mode: c++ $
    # End: $
  */

