/* $Id: nsitetles.cpp,v 1.2 2001/01/04 20:06:18 mstorti Exp $ */

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/fastmat.h"

#include "nsi_tet.h"

#define ADD_GRAD_DIV_U_TERM
#define STANDARD_UPWIND
#define USE_FASTMAT

#ifndef USE_FASTMAT
#error "Solo anda con FastMat por ahora"
#endif

extern TextHashTable *GLOBAL_OPTIONS;

#define STOP {PetscFinalize(); exit(0);}
   
#ifndef USE_FASTMAT //---:---<*>---:---<*>---:---<*>---:---<*>---:
#define MATRIX_LIB Matrix
#undef SHV_OPTIONS
#define SHV_OPTIONS << setprecision(7) << setw(12)
//#define MSHV SHV
#define MSHV(s) {cout << #s << endl;  newmat_print(s);}
#else               //---:---<*>---:---<*>---:---<*>---:---<*>---:
#define MATRIX_LIB FastMat
#define MSHV FMSHV
#endif              //---:---<*>---:---<*>---:---<*>---:---<*>---:
#define MAXPROP 100

// Anular impresiones
#if 0
#undef SHV
#undef MSHV
#undef STOP
#define SHV(a) {}
#define MSHV(a) {}
#define STOP {}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "nsi_tet_les::assemble"
int nsi_tet_les::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		      Dofmap *dofmap,const char *jobinfo,int myrank,
		      int el_start,int el_last,int iter_mode,
		      const TimeData *time_) {

  double time = *(Time *)time_;
  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsi_tet\n");

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define IDENT(j,k) (ident[ndof*(j)+(k)]) 
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)


  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  int nen = nel*ndof;

  // Unpack Dofmap
  int *ident,neq,nnod;
  neq = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;
  if (comp_mat) {
    retvalmat = arg_data_v[0].retval;
  }

  double *hmin,Dt; 
  if (comp_mat_res) {
    locst = arg_data_v[0].locst;
    locst2 = arg_data_v[1].locst;
    retval = arg_data_v[2].retval;
    retvalmat = arg_data_v[3].retval;
    hmin = (arg_data_v[4].vector_assoc)->begin();
#define WAS_SET arg_data_v[4].was_set
    Dt = *(double *)(arg_data_v[5].user_data);
  }

  // Is =1 if we use a weak form for the gradient of pressure term
  SGETOPTDEF(int,weak_form,1);
  SGETOPTDEF(double,shock_capturing_factor,1);

  // allocate local vecs
  int kdof;
  MATRIX_LIB veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof), 
         locstate2(nel,ndof),xpg,G_body(ndim,1);

  if (ndof != ndim+1) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  MATRIX_LIB matloc(nen,nen), matlocmom(nel,nel), masspg(nel,nel),
    grad_u_ext(ndof,ndof);

#ifndef USE_FASTMAT
  grad_u_ext = 0;
#else
  grad_u_ext.set(0.);
#endif

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  GGETOPTDEF(int,LES,0);    // Trapezoidal method parameter. 
  GGETOPTDEF(double,C_smag,0.18); // Dijo Beto
  GGETOPTDEF(double,alpha,1.);    // Trapezoidal method parameter. 
  SGETOPTDEF(double,tau_fac,1.);  // Scale upwind

  // This is the volumetric force. It's a constant vector,
  // set to zero if not defined 
  G_body.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"G_body",G_body.store,1,ndim);

  assert(LES==0);
  double pi = 4*atan(1.0);

  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))

  int nprops=iprop;
  
  double rho=1.;

  // Gauss Point data
  char *geom;
  thash->get_entry("geometry",geom);
  
#ifndef USE_FASTMAT
  GPdata gp_data(geom,ndim,nel,npg);
#else
  //GPdata gp_data(geom,ndim,nel,npg);
  GPdata gp_data(geom,ndim,nel,npg,GP_FASTMAT);
#endif

  // Definiciones para descargar el lazo interno
  double detJaco, UU, u2, Peclet, psi, tau_supg, tau_pspg, div_u_star,
    p_star,wpgdet,velmod,tol,h_supg,fz,delta_supg,Uh;

#ifndef USE_FASTMAT

  RowVector P_supg, dmatw;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  Matrix dshapex(ndim,nel),Jaco(ndim,ndim),iJaco(ndim,ndim),
    grad_u(ndim,ndim),grad_u_star(ndim,ndim),resmom(nel,ndim),
    matij(ndof,ndof),Uintri,P_pspg;

  ColumnVector grad_p_star(ndim),u(ndim),u_star(ndim),
    uintri(ndim),rescont(nel),dmatu(ndim),svec(ndim);;

  DiagonalMatrix eye(ndim);
  eye = 1;

#else

  FastMat P_supg, W_supg, W_supg_t,dmatw;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FastMat dshapex,dshapext,Jaco(ndim,ndim),iJaco(ndim,ndim),
    grad_u(ndim,ndim),grad_u_star,strain_rate(ndim,ndim),
    grad_u_star_t,resmom(nel,ndim),dresmom,
    matij(ndof,ndof),Uintri,P_pspg,P_pspg_t,svec;

  FastMat grad_p_star(ndim,1),ut,u,u_star,u_star_t,du,
    uintri(ndim,1),rescont(nel,1),dmatu(ndim,1),ucols,ucols_new,
    ucols_star,pcol_star,pcol_new,pcol,fm_p_star,tmp1,tmp2,tmp3,tmp4,tmp5,
    massm,tmp7,tmp8,tmp9,tmp10;

  double tmp12;

  FastMat eye;
  eye.eye(ndim);

#endif
  MATRIX_LIB seed,one_nel,matloc_prof(nen,nen);

  if (comp_mat) {
#ifndef USE_FASTMAT

#ifdef ADD_GRAD_DIV_U_TERM
    matloc_prof=1;
#else
    seed.ReSize(ndof,ndof);
    one_nel.ReSize(nel,nel);
    matloc_prof.ReSize(nen,nen);
    for (int jj=1; jj<=ndim; jj++) {
      seed(jj,jj)=1.;
      seed(jj,ndof)=1.;
      seed(ndof,jj)=1.;
    }
    seed(ndof,ndof)=1.;
    one_nel=1;
    matloc_prof = kron(one_nel,seed);
#endif

#else

#ifdef ADD_GRAD_DIV_U_TERM
    matloc_prof.set(1.);
#else
    seed.set_size(ndof,ndof);
    seed.set(0.);
    one_nel.set_size(nel,nel);
    one_nel.set(0.);
    matloc_prof.set_size(nen,nen);
    for (int jj=1; jj<=ndim; jj++) {
      seed.set(jj,jj,1.);
      seed.set(jj,ndof,1.);
      seed.set(ndof,jj,1.);
    }
    seed.set(ndof,ndof,1.);
    one_nel.set(1.);
    kron(matloc_prof,one_nel,seed);
#endif

#endif

  }

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    ielh++;
    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    elem = k;

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
#ifndef USE_FASTMAT
	xloc(kloc+1,jdim+1) = NODEDATA(node-1,jdim);
#else
	xloc.set(kloc+1,jdim+1,NODEDATA(node-1,jdim));
#endif
      }
    }

    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_mat_res) {
#ifndef USE_FASTMAT
      locstate << &(LOCST(ielh,0,0));
      locstate2 << &(LOCST2(ielh,0,0));
#else
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
#endif
    }

#ifndef USE_FASTMAT
    matlocmom = 0;
    matloc = 0;
    veccontr = 0;
    resmom = 0;
    rescont = 0;
#else
    matlocmom.set(0.);
    matloc.set(0.);
    veccontr.set(0.);
    resmom.set(0.);
    rescont.set(0.);
#endif


#ifdef USE_FASTMAT
    ucols.columns(locstate2,1,ndim);
    ucols_new.columns(locstate,1,ndim);
    pcol.columns(locstate2,ndof,ndof);
    pcol_new.columns(locstate,ndof,ndof);

    ucols_star.set(ucols_new);
    ucols_star.scale(alpha);
    FMaxpy(ucols_star,1-alpha,ucols);

    pcol_star.set(pcol_new);
    pcol_star.scale(alpha);
    FMaxpy(pcol_star,1-alpha,pcol);

#endif

#ifndef USE_FASTMAT  //----------------------------
#define DSHAPEXI (gp_data.dshapexi[ipg])
#define SHAPE    (gp_data.shape[ipg])
#define WPG      (gp_data.wpg[ipg])
#else                //----------------------------
#define DSHAPEXI (*gp_data.FM_dshapexi[ipg])
#define SHAPE    (*gp_data.FM_shape[ipg])
#define WPG      (gp_data.wpg[ipg])
#endif               //----------------------------
    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

#ifndef USE_FASTMAT
      Jaco = DSHAPEXI * xloc;
      detJaco = mydet(Jaco);
#else
      FMp(xpg,SHAPE,xloc);
      FMp(Jaco,DSHAPEXI,xloc);
      FMdet(detJaco,Jaco);
#endif

      if (detJaco <= 0.) {
	cout << "Jacobian of element " << elem << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	exit(1);
      }
      wpgdet = detJaco*WPG;

#ifndef USE_FASTMAT
      iJaco = Jaco.i();
      dshapex = iJaco * DSHAPEXI;
#else
      FMinv(iJaco,Jaco);
      FMp(dshapex,iJaco,DSHAPEXI);
      dshapext.transpose(dshapex);
#endif

      double Area   = npg*wpgdet;
      double h_pspg,Delta;
      if (ndim==2) {
	h_pspg = sqrt(4.*Area/pi);
	Delta = sqrt(Area);
      } else if (ndim==3) {
	// h_pspg = pow(6*Area/pi,1./3.);
	// El pow() da segmentation violation cuando corro con -O !!
	h_pspg = cbrt(6*Area/pi);
	Delta = cbrt(Area);
      } else {
	PFEMERRQ("Only dimensions 2 and 3 allowed for this element.\n");
      }
      
      if (comp_mat_res) {
	// computes the minimum size of the mesh
	if (!WAS_SET || h_pspg<*hmin) {
	  WAS_SET = 1;
	  *hmin = h_pspg;
	}

	// state variables and gradient
#ifndef USE_FASTMAT
	u = (SHAPE * locstate2.Columns(1,ndim)).t();
#else
	FMp(ut,SHAPE,ucols);
	u.transpose(ut);
#endif

#ifndef USE_FASTMAT
	p_star = (SHAPE * locstate.Column(ndim+1)).AsScalar();
	u_star = (SHAPE * locstate.Columns(1,ndim)).t();
#else
	FMp(fm_p_star,SHAPE,pcol_star);
	fm_p_star.as_scalar(p_star);
	FMp(u_star_t,SHAPE,ucols_star);
	u_star.transpose(u_star_t);
#endif

#ifndef USE_FASTMAT
	grad_u      = dshapex * locstate2.Columns(1,ndim);
	grad_u_star = dshapex * locstate.Columns(1,ndim);
	grad_p_star = dshapex * locstate.Column(ndim+1);
#else
	FMp(grad_u,dshapex,ucols);
	FMp(grad_u_star,dshapex,ucols_star);
	grad_u_star_t.transpose(grad_u_star);
	FMp(grad_p_star,dshapex,pcol_star);

	strain_rate.set(grad_u_star);
	FMa(strain_rate,strain_rate,grad_u_star_t);
	strain_rate.scale(0.5);

	// Smagorinsky turbulence model
	double nu_eff;
	if (LES) {
	  double tr;
	  strain_rate.trace_of_product(strain_rate,tr);
	  double nu_t = SQ(C_smag*Delta)*sqrt(2*tr);
	  nu_eff = VISC + nu_t;
	} else {
	  nu_eff = VISC;
	}
#endif

#ifndef USE_FASTMAT
	double u2 = u.SumSquare();
	uintri = iJaco * u;
	Uh = sqrt(uintri.SumSquare())/2.;
#else
	u.sum_square(u2);
	FMp(uintri,iJaco,u);
	uintri.sum_square(Uh);
	Uh = sqrt(Uh)/2;
#endif

#ifdef STANDARD_UPWIND
	velmod = sqrt(u2);
        tol=1.0e-16;
        h_supg=0;
        if(velmod>tol) {
#ifndef USE_FASTMAT
          svec = u/velmod; 
#else
	  svec.set(u);
	  svec.scale(1./velmod);
#endif
#ifndef USE_FASTMAT
	  h_supg = (dshapex.t() * svec).Norm1();
#else
	  FMp(tmp9,dshapext,svec);
	  h_supg = tmp9.norm1();
#endif
          h_supg = (h_supg < tol ? tol : h_supg);
          h_supg = 2./h_supg;
        } else {
          h_supg = h_pspg;
        }

	Peclet = velmod * h_supg / (2. * nu_eff);
//	psi = 1./tanh(Peclet)-1/Peclet;
//	tau_supg = psi*h_supg/(2.*velmod);

#if 0 // pow() crashes when optimization is enabled!!
        tau_supg = pow(2./Dt,2.)+pow(2.*velmod/h_supg,2.)
	  +9.*pow(4.*nu_eff/pow(h_supg,2),2.);
        tau_supg = 1./sqrt(tau_supg);

        tau_pspg = pow(2./Dt,2.)+pow(2.*velmod/h_pspg,2.)
	  +9.*pow(4.*nu_eff/pow(h_pspg,2),2.);
#else
        tau_supg = SQ(2./Dt)+SQ(2.*velmod/h_supg)
	  +9.*SQ(4.*nu_eff/SQ(h_supg));
        tau_supg = 1./sqrt(tau_supg);

        tau_pspg = SQ(2./Dt)+SQ(2.*velmod/h_pspg)
	  +9.*SQ(4.*nu_eff/SQ(h_pspg));
#endif
        tau_pspg = 1./sqrt(tau_pspg);

        fz = (Peclet < 3 ? Peclet/3 : 1);
        delta_supg = 0.5*h_supg*velmod*fz;
	
	if (tau_fac != 1.) {
	  tau_pspg *= tau_fac;
	  tau_supg *= tau_fac;
	}

#else
	if(u2<=1e-6*(2. * Uh * nu_eff)) {
	  Peclet=0.;
	  psi=0.;
	  tau_supg=0.;
	} else {
	  Peclet = u2 / (2. * Uh * nu_eff);
	  psi = 1./tanh(Peclet)-1/Peclet;
	  tau_supg = psi/(2.*Uh);
	}

	// PSPG parameter for the stabilization of the
	// incompressibility conditions
	tau_pspg = h_pspg*h_pspg/ 2 / ( 6*nu_eff + sqrt(u2)*h_pspg ) ;
	PFEMERRQ("Not implemented yet shock capturing with standard upwind\n");
	double delta_supg=1e-8;
#endif
	
	delta_supg *= shock_capturing_factor;
	
	// P_supg es un vector fila
#ifndef USE_FASTMAT
	P_supg  = tau_supg * u.t() * dshapex;
#else
	FMp(P_supg,ut,dshapex);
	P_supg.scale(tau_supg);

	// Weight function and transposed
	FMa(W_supg,SHAPE,P_supg);
	W_supg_t.transpose(W_supg);
#endif

	// Pressure stabilizing term
#ifndef USE_FASTMAT
	P_pspg  = tau_pspg /rho *dshapex;
#else
	P_pspg.set(dshapex);
	P_pspg.scale(tau_pspg/rho);
	P_pspg_t.transpose(P_pspg);
#endif

#ifdef USE_FASTMAT
#endif

	// implicit version - Backward Euler

#ifndef USE_FASTMAT
	dmatu = (1/Dt)* (u_star-u) + (u.t() * grad_u_star).t() ;
#else
	FMp(dmatu,grad_u_star_t,u_star);

	du.set(u_star);
	FMaxpy(du,-1.,u);
	FMaxpy(dmatu,1/(alpha*Dt),du);
	FMaxpy(dmatu,-1,G_body);
#endif

	
#ifndef USE_FASTMAT
	div_u_star = (dshapex * locstate.Columns(1,ndim)).Trace();
#else
	dshapex.trace_of_product(ucols_star,div_u_star);
#endif


	// Galerkin - momentum
	// resmom tiene que tener nel*ndim
#ifndef USE_FASTMAT
	resmom -= wpgdet * rho * SHAPE.t() * dmatu.t();	
#else
	FMp(dresmom,dmatu,SHAPE);
	FMaxpy_t(resmom,-wpgdet * rho,dresmom);
#endif

	if (weak_form) {
#ifndef USE_FASTMAT
	  resmom -= wpgdet * dshapex.t() *
	    (nu_eff * grad_u_star - p_star * eye);	
#else
	  tmp1.set(strain_rate);
	  tmp1.scale(nu_eff);
	  FMaxpy(tmp1,-p_star,eye);
	  FMp(tmp2,dshapext,tmp1);
	  FMaxpy(resmom,-wpgdet,tmp2);
#endif
	} else {
#ifndef USE_FASTMAT
	  resmom -= wpgdet * dshapex.t() * (nu_eff * grad_u_star );	
	  resmom -= wpgdet * SHAPE.t() * grad_p_star.t() ;	
#else
	  FMp(tmp1,dshapext,strain_rate);
	  FMaxpy(resmom,-wpgdet*nu_eff,tmp1);
	  FMp(tmp10,grad_p_star,SHAPE);
	  FMaxpy_t(resmom,-wpgdet,tmp10);
#endif

	}

	// SUPG perturbation - momentum
#ifndef USE_FASTMAT
	resmom -= wpgdet * P_supg.t() * (rho *
					 dmatu+grad_p_star).t();
#else
	tmp3.set(grad_p_star);
	FMaxpy(tmp3,rho,dmatu);
	FMp(tmp4,tmp3,P_supg);
	FMaxpy_t(resmom,-wpgdet,tmp4);
#endif

        // shock capturing term - momentum
#ifndef USE_FASTMAT
        resmom -= wpgdet * delta_supg * rho * dshapex.t() * div_u_star;
#else
	FMaxpy(resmom, -wpgdet*delta_supg*rho*div_u_star, dshapext);
#endif

	// Galerkin - continuity
#ifndef USE_FASTMAT
	rescont += wpgdet *  SHAPE.t() * div_u_star;
#else
	FMaxpy_t(rescont,wpgdet*div_u_star,SHAPE);
#endif

	// PSPG perturbation - continuity
#ifndef USE_FASTMAT
	rescont += wpgdet * P_pspg.t() * (rho * dmatu +
					  grad_p_star);
#else
	FMp(tmp5,P_pspg_t,tmp3);
	FMaxpy(rescont,wpgdet,tmp5);
#endif

	// Parte temporal + convectiva (Galerkin)
#ifndef USE_FASTMAT
	matlocmom = (SHAPE+P_supg).t() * rho * ((1/Dt)*SHAPE + u.t() *
						dshapex);
#else
	FMp(massm,u_star_t,dshapex);
	FMaxpy(massm,1/(alpha*Dt),SHAPE);
	FMp(matlocmom,W_supg_t,massm);
	matlocmom.scale(rho);
#endif

	// Parte difusiva
#ifndef USE_FASTMAT
	matlocmom += nu_eff * dshapex.t() * dshapex ;
#else
	FMp(tmp7,dshapext,dshapex);
	FMaxpy(matlocmom,0.5 * nu_eff,tmp7);
#endif

#ifndef USE_FASTMAT
	dmatw =  rho * ((1/Dt)*SHAPE + u.t() * dshapex);
#else
	dmatw.set(massm);
	dmatw.scale(rho);
#endif

	for (int iloc=1; iloc<=nel; iloc++) {
	  for (int jloc=1; jloc<=nel; jloc++) {
#ifndef USE_FASTMAT
	    matij=0;
#else
	    matij.set(0.);
#endif
	    for (jdim=1; jdim<=ndim; jdim++) {
#ifndef USE_FASTMAT
	      matij(jdim,jdim) += matlocmom(iloc,jloc);
#else
	      double c;
	      matlocmom.get(iloc,jloc,c);
	      matij.add(jdim,jdim, c);
#endif
	      if (!weak_form) {
#ifndef USE_FASTMAT
		matij(jdim,ndof) +=  SHAPE(iloc) * dshapex(jdim,jloc) +
		  P_supg(iloc) * dshapex(jdim,jloc);
#else
		matij.add(jdim,ndof,SHAPE.get(1,iloc) * dshapex.get(jdim,jloc) +
		  P_supg.get(1,iloc) * dshapex.get(jdim,jloc));
#endif
              } else {
#ifndef USE_FASTMAT
		matij(jdim,ndof) += - dshapex(jdim,iloc) * SHAPE(jloc) +
		  P_supg(iloc) * dshapex(jdim,jloc);
#else
		matij.add(jdim,ndof,-dshapex.get(jdim,iloc) * SHAPE.get(1,jloc)
			  + P_supg.get(1,iloc) * dshapex.get(jdim,jloc));
#endif
              }
#ifndef USE_FASTMAT
	      matij(ndof,jdim) -= P_pspg(jdim,iloc)*dmatw(jloc)
		+ SHAPE(iloc)*dshapex(jdim,jloc);
	      matij(ndof,ndof) -=
		P_pspg(jdim,iloc)*dshapex(jdim,jloc);
#else
	      matij.add(ndof,jdim, -(P_pspg.get(jdim,iloc)*dmatw.get(1,jloc)
		+ SHAPE.get(1,iloc)*dshapex.get(jdim,jloc)));
	      matij.add(ndof,ndof,-(P_pspg.get(jdim,iloc)*dshapex.get(jdim,jloc)));
#endif
#ifdef USE_FASTMAT
	      tmp12 = (delta_supg * rho + 0.5*nu_eff) * 
		  dshapex.get(jdim,iloc);
#endif
	      // shock capturing term - momentum
	      for (int kdim=1; kdim<=ndim; kdim++) {
#ifndef USE_FASTMAT
                matij(jdim,kdim) += delta_supg * rho * 
		  dshapex(jdim,iloc) * dshapex(kdim,jloc);
#else
		matij.add(jdim,kdim, tmp12 * dshapex.get(kdim,jloc));
		//			  + 0.5*nu_eff*dshapex.get(jdim,jloc)*dshapex.get(kdim,kloc));
#endif
              }
	    }
	    int il1=(iloc-1)*ndof+1;
	    int il2=il1+ndof-1;
	    int jl1=(jloc-1)*ndof+1;
	    int jl2=jl1+ndof-1;
#ifndef USE_FASTMAT
	    matloc.SubMatrix(il1,il2,jl1,jl2) += wpgdet * matij;
#else
	    tmp8.set(matij);
	    tmp8.scale(wpgdet);
	    matloc.add(il1,il2,jl1,jl2,tmp8);
#endif
	    
	  }
	}

      } else if (comp_mat) {
	// don't make anything here !!
      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Don't know how to compute jobinfo: %s\n",jobinfo);
	CHKERRQ(ierr);
      }

    }

    if(comp_mat) {
#ifndef USE_FASTMAT
      matloc_prof >> &(RETVALMAT(ielh,0,0,0,0));
#else
      matloc_prof.copy(&(RETVALMAT(ielh,0,0,0,0)));
#endif
    }      

    if (comp_mat_res) {
#ifndef USE_FASTMAT
      veccontr.Columns(1,ndim) = resmom;
      veccontr.Column(ndim+1) = rescont;
      veccontr >> &(RETVAL(ielh,0,0));
      matloc >> &(RETVALMAT(ielh,0,0,0,0));
#else
      veccontr.set(1,nel,1,ndim,resmom);
      veccontr.set(1,nel,ndof,ndof,rescont);
      veccontr.copy(&(RETVAL(ielh,0,0)));
      matloc.copy(&(RETVALMAT(ielh,0,0,0,0)));
#endif
    }
  }
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
