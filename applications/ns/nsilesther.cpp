//__INSERT_LICENSE__
//$Id: nsilesther.cpp,v 1.11 2001/05/30 18:21:50 mstorti Exp $

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/fastmat2.h"

#include "nsi_tet.h"

//#define ADD_GRAD_DIV_U_TERM
#define STANDARD_UPWIND
#define USE_FASTMAT

extern TextHashTable *GLOBAL_OPTIONS;

#define STOP {PetscFinalize(); exit(0);}
   
#define MAXPROP 100

#define SQ(n) ((n)*(n))

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "nsi_tet_les_fm2::assemble"
int nsi_tet_les_ther::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			       Dofmap *dofmap,const char *jobinfo,int myrank,
			       int el_start,int el_last,int iter_mode,
			       const TimeData *time_) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(get_nearest_wall_element);

// added for thermal equation
  GET_JOBINFO_FLAG(comp_mat_th);
  GET_JOBINFO_FLAG(comp_mat_res_th);
  GET_JOBINFO_FLAG(comp_res_th);

  comp_mat_res_th=comp_mat_res;
  comp_mat_th=comp_mat;

// end added

  // Essentially treat comp_res as comp_mat_res but
  // with the side effect of update_jacobian=1
  int update_jacobian=1;	
  if (comp_res) {
    comp_mat_res=1;
    update_jacobian=0;
  }

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsilesther\n");

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)
#define NN_IDX(j) ELEMIPROPS_ADD(j,0)
#define IDENT(j,k) (ident[ndof*(j)+(k)]) 
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)


  int locdof,kldof,lldof;
  char *value;

  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  // ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  TGETOPTNDEF(thash,int,ndim,none); //nd
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
  WallData *wall_data;
  if (comp_mat || comp_mat_th) {
    retvalmat = arg_data_v[0].retval;
  } else if (get_nearest_wall_element) {
    wall_data = (WallData *)arg_data_v[0].user_data;
  }

  double *hmin,Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res || comp_mat_res_th) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    if (update_jacobian) retvalmat = arg_data_v[ja++].retval;
#ifdef RH60    // fixme:= STL vector compiler bug??? see notes.txt
    hmin = (arg_data_v[ja++].vector_assoc)->begin();
#else
    ja++;
#endif
    ja_hmin=ja;
    Dt = *(double *)(arg_data_v[ja++].user_data);
    wall_data = (WallData *)arg_data_v[ja++].user_data;
  } 

  //o Use a weak form for the gradient of pressure term.
  SGETOPTDEF(int,weak_form,1);
  //o Add shock-capturing term.
  SGETOPTDEF(double,shock_capturing_factor,0);

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof), 
         locstate2(nel,ndof),xpg,G_body(ndim),gravity(ndim);

  if (ndof != ndim+2) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+2\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FMatrix matloc(nen,nen), matlocmom(nel,nel), masspg(nel,nel),
    grad_u_ext(ndof,ndof);
  FastMat2 matlocf(4,nel,ndof,nel,ndof);

  grad_u_ext.set(0.);

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  //o Density
  SGETOPTDEF(double,rho,1.);
  //o Thermal conductivity
  SGETOPTDEF(double,kappa,1.);
  //o Specific heat - constant pressure 
  SGETOPTDEF(double,Cp,1.);
  //o buoyancy coefficient for Boussinesq term
  SGETOPTDEF(double,betath,0.);
  //o Reference temperature for thermal coupling Boussinesq term
  SGETOPTDEF(double,Tinfty,0.);
  //o Gravity acceleration for buoyancy terms
  gravity.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"gravity",gravity.storage_begin(),1,ndim);

  //o Add LES for this particular elemset.
  GGETOPTDEF(int,LES,0);
  //o Cache \verb+grad_div_u+ matrix
  SGETOPTDEF(int,cache_grad_div_u,0);
  //o Smagorinsky constant.
  SGETOPTDEF(double,C_smag,0.18); // Dijo Beto
  //o van Driest constant for the damping law.
  SGETOPTDEF(double,A_van_Driest,26); 
  //o Parameter for the trapezoidal rule time integration method. 
  GGETOPTDEF(double,alpha,1.);
  //o Scale the SUPG upwind term. 
  SGETOPTDEF(double,tau_fac,1.);  // Scale upwind
  //o Adjust the stability parameters, taking into account
  // the time step. 
  SGETOPTDEF(double,temporal_stability_factor,0.);  // Scale upwind

  //o _T: double[ndim] _N: G_body _D: null vector 
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"G_body",G_body.storage_begin(),1,ndim);

  SGETOPTDEF(double,Q_body,0.);  // body thermal flux

  double pi = 4*atan(1.0);

  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))

// for thermal equation
#if 0
  DEFPROP(density);
#define rho (*(propel+density_indx))
  DEFPROP(specific_heat);
#define CP (*(propel+specific_heat_indx))
  DEFPROP(conductivity);
#define KAPPA (*(propel+conductivity_indx))
#endif

  int nprops=iprop;
  
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  
  //GPdata gp_data(geom,ndim,nel,npg);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco, UU, u2, Peclet, psi, tau_supg, tau_pspg, div_u_star,
    p_star,wpgdet,velmod,tol,h_supg,fz,delta_supg,Uh;
  double tau_supg_th,T_star,dT,tmp1_th,tmp2_th;

  FastMat2 P_supg, W_supg, W_supg_t, dmatw, W_supg_th,P_supg_th,
           grad_div_u(4,nel,ndim,nel,ndim);
  double *grad_div_u_cache;
  int grad_div_u_was_cached;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FMatrix dshapex,dshapext,Jaco(ndim,ndim),iJaco(ndim,ndim),
    grad_u(ndim,ndim),grad_u_star,strain_rate(ndim,ndim),resmom(nel,ndim),
    dresmom(nel,ndim),matij(ndof,ndof),Uintri,P_pspg,svec;

  FMatrix grad_p_star(ndim),u,u_star,du,
    uintri(ndim),rescont(nel),dmatu(ndim),ucols,ucols_new,
    ucols_star,pcol_star,pcol_new,pcol,fm_p_star,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,
    massm,tmp7,tmp8,tmp9,tmp10,tmp11,tmp13,tmp14,tmp15,dshapex_c,xc,
    wall_coords(ndim),dist_to_wall,tmp16,tmp17,tmp18,tmp19,tmp20;

  // FMatrix T,Tcol,Tcol_new,Tcol_star,grad_T_star(ndim),massm_th,
  double T;
  FMatrix Tcol,Tcol_new,Tcol_star,grad_T_star(ndim),massm_th,
    tmp3_th,tmp4_th,tmp5_th,tmp6_th,tmp7_th,tmp8_th,tmp9_th,
    tmp10_th,tmp11_th,tmp12_th,tmp13_th,tmp14_th,tmp15_th,tmp16_th,
    tmp17_th,tmp18_th,tmp19_th;
  FMatrix resther(nel),matlocther(nel,nel);
  FastMat2 matlocmomther(3,nel,ndim,nel);

  double tmp12;
  double nu_eff;
  double tsf = temporal_stability_factor;

  FMatrix eye(ndim,ndim),seed,one_nel,matloc_prof(nen,nen);;
  eye.eye();

  if (comp_mat) {

    matloc_prof.set(1.);

#if 0
#ifdef ADD_GRAD_DIV_U_TERM
#else
    seed.resize(2,ndof,ndof);
    seed.set(0.);
    one_nel.resize(2,nel,nel);
    one_nel.set(0.);
    matloc_prof.resize(2,nen,nen);
    for (int jj=1; jj<=ndim; jj++) {
      seed.setel(jj,jj,1.);
      seed.setel(jj,ndof,1.);
      seed.setel(ndof,jj,1.);
    }
    seed.setel(ndof,ndof,1.);
    one_nel.set(1.);
    matloc_prof.kron(one_nel,seed);
#endif
#endif

  }

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    elem = k;

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    if (get_nearest_wall_element) {
      assert(LES);
      xc.sum(xloc,-1,1).scale(1./double(nel));
      int nn;
      wall_data->nearest(xc.storage_begin(),nn);
      NN_IDX(k) = nn;
      continue;
    }

    double grad_div_u_coef=0.;	// multiplies grad_div_u term
    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_mat_res || comp_res || comp_mat_res_th || comp_res_th) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));

      if (cache_grad_div_u) {
	grad_div_u_cache = (double *)local_store_address(ielh);
	grad_div_u_was_cached = (grad_div_u_cache!=NULL);
	if (!grad_div_u_was_cached) {
	  local_store_address(ielh) 
	    = grad_div_u_cache 
	    = new double[ndim*ndim*nel*nel];
	}
      }
    }

    matlocf.set(0.);
    veccontr.set(0.);
    matlocmom.set(0.);
    matlocther.set(0.);
    matlocmomther.set(0.);
    resmom.set(0.);
    rescont.set(0.);
    resther.set(0.);

    if (comp_mat_res && cache_grad_div_u) {
      if (grad_div_u_was_cached) {
	grad_div_u.set(grad_div_u_cache);
      } else {
	grad_div_u.set(0.);
      }
    }

    ucols.set(locstate2.is(2,1,ndim));
    pcol.set(locstate2.rs().ir(2,ndim+1));
    Tcol.set(locstate2.rs().ir(2,ndof));
    locstate2.rs();

    ucols_new.set(locstate.is(2,1,ndim));
    pcol_new.set(locstate.rs().ir(2,ndim+1));
    Tcol_new.set(locstate.rs().ir(2,ndof));
    locstate.rs();

    ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
    pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
    Tcol_star.set(Tcol_new).scale(alpha).axpy(Tcol,1-alpha);
    
    double shear_vel;
    int wall_elem;
    if (LES && comp_mat_res) {
      Elemset *wall_elemset;
      const double *wall_coords_;
      wall_data->nearest_elem_info(NN_IDX(k),wall_elemset,wall_elem,wall_coords_);
      wall_coords.set(wall_coords_);
      shear_vel = wall_elemset->elemprops_add[wall_elem];
    }

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);

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
      dshapex_c.set(dshapex);

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
      
      if (comp_mat_res || comp_mat_res_th) {
	// computes the minimum size of the mesh
	if (!WAS_SET || h_pspg<*hmin) {
	  WAS_SET = 1;
#ifdef RH60    // fixme:= STL vector compiler bug??? see notes.txt
	  *hmin = h_pspg;
#endif
	}

	// state variables and gradient
	u.prod(SHAPE,ucols,-1,-1,1);
	T = double(tmp8.prod(SHAPE,Tcol,-1,-1));
	//T.prod(SHAPE,Tcol,-1,-1);

	p_star = double(tmp8.prod(SHAPE,pcol_star,-1,-1));
	T_star = double(tmp8.prod(SHAPE,Tcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);

	// grad_u.prod(dshapex,ucols,1,-1,-1,2);

	grad_u_star.prod(dshapex,ucols_star,1,-1,-1,2);
	grad_p_star.prod(dshapex,pcol_star,1,-1,-1);
	grad_T_star.prod(dshapex,Tcol_star,1,-1,-1);

	u2 = u.sum_square_all();
	velmod = sqrt(u2);

        if (comp_mat_res || comp_res) {
 	  strain_rate.set(grad_u_star);
	  grad_u_star.t();
	  strain_rate.add(grad_u_star).scale(0.5);
	  grad_u_star.rs();

	  // Smagorinsky turbulence model
	  if (LES) {
	    double tr = (double) tmp15.prod(strain_rate,strain_rate,-1,-2,-1,-2);

	    dist_to_wall.prod(SHAPE,xloc,-1,-1,1).rest(wall_coords);
	    double ywall = sqrt(dist_to_wall.sum_square_all());
	    double y_plus = ywall*shear_vel/VISC;
	    double van_D = 1.-exp(-y_plus/A_van_Driest);
	  
	    double nu_t = SQ(C_smag*Delta*van_D)*sqrt(2*tr);
	    nu_eff = VISC + nu_t;
	  } else {
	    nu_eff = VISC;
	  }
        }

	// u2 = u.sum_square_all();
	uintri.prod(iJaco,u,1,-1,-1);
	Uh = uintri.sum_square_all();
	Uh = sqrt(Uh)/2;

	// velmod = sqrt(u2);
        tol=1.0e-16;
        h_supg=0;
	FastMat2::branch();
        if(velmod>tol) {
	  FastMat2::choose(0);
	  svec.set(u).scale(1./velmod);
	  h_supg = tmp9.prod(dshapex,svec,-1,1,-1).sum_abs_all();
          h_supg = (h_supg < tol ? tol : h_supg);
          h_supg = 2./h_supg;
        } else {
          h_supg = h_pspg;
        }
	FastMat2::leave();

        if (comp_mat_res || comp_res) {
	  Peclet = velmod * h_supg / (2. * nu_eff);
          tau_supg = tsf*SQ(2./Dt)+SQ(2.*velmod/h_supg)
	    +9.*SQ(4.*nu_eff/SQ(h_supg));
          tau_supg = 1./sqrt(tau_supg);

          tau_pspg = tsf*SQ(2./Dt)+SQ(2.*velmod/h_pspg)
	    +9.*SQ(4.*nu_eff/SQ(h_pspg));
          tau_pspg = 1./sqrt(tau_pspg);

          fz = (Peclet < 3 ? Peclet/3 : 1);
          delta_supg = 0.5*h_supg*velmod*fz;
	
	  if (tau_fac != 1.) {
	    tau_pspg *= tau_fac;
	    tau_supg *= tau_fac;
	  }
	  delta_supg *= shock_capturing_factor;

        } 
        if (comp_mat_res_th || comp_res_th) {
	  if (velmod>tol) {
	    double diff_coe=kappa/rho/Cp;
	    Peclet = velmod * h_supg / (2. * diff_coe);
	    psi = 1./tanh(Peclet)-1/Peclet;
	    tau_supg_th = psi*h_supg/(2.*velmod);
	  } else {
	    tau_supg_th = 0;
	  }
        }

        if (comp_mat_res || comp_res) {	
	  // P_supg es un vector fila
	  P_supg.prod(u,dshapex,-1,-1,1).scale(tau_supg);

	  // Weight function 
	  W_supg.set(P_supg).add(SHAPE);

	  // Pressure stabilizing term
	  P_pspg.set(dshapex).scale(tau_pspg/rho);
        }

        if (comp_mat_res_th || comp_res_th) {	

	  // P_supg stabilizing term for thermal equation
	  P_supg_th.prod(u,dshapex,-1,-1,1).scale(tau_supg_th);
	  W_supg_th.set(P_supg_th).add(SHAPE);

        }

	// implicit version - General Trapezoidal rule - parameter alpha

        if (comp_mat_res || comp_res) {
#if ADD_GRAD_DIV_U_TERM
	  dmatu.prod(u_star,grad_u_star,-1,-1,1);
#else
	  dmatu.prod(u,grad_u_star,-1,-1,1);
#endif

	  du.set(u_star).rest(u);
	  dmatu.axpy(du,1/(alpha*Dt)).rest(G_body);
	
	  div_u_star = double(tmp10.prod(dshapex,ucols_new,-1,-2,-2,-1));

	  // Galerkin - momentum
	  // resmom tiene que tener nel*ndim
	  dresmom.t().prod(dmatu,SHAPE,1,2).rs();
	  resmom.axpy(dresmom,-wpgdet * rho);

	  if (weak_form) {
	    tmp1.set(strain_rate).scale(2*nu_eff).axpy(eye,-p_star);
	    tmp2.prod(dshapex,tmp1,-1,1,-1,2);
	    resmom.axpy(tmp2,-wpgdet);
	  } else {
	    tmp6.prod(dshapex,strain_rate,-1,1,-1,2).scale(2*nu_eff);
	    tmp11.prod(SHAPE,grad_p_star,1,2).add(tmp6);
	    resmom.axpy(tmp11,-wpgdet);
	  }

	  // SUPG perturbation - momentum
	  tmp3.set(grad_p_star).axpy(dmatu,rho);
	  tmp4.prod(P_supg,tmp3,1,2);
	  resmom.axpy(tmp4,-wpgdet);

          // shock capturing term - momentum
	  resmom.axpy(dshapex.t(),-wpgdet*delta_supg*rho*div_u_star);
	  dshapex.rs();

	  // Galerkin - continuity
	  rescont.axpy(SHAPE,wpgdet*div_u_star);

	  // PSPG perturbation - continuity
	  tmp5.prod(P_pspg,tmp3,-1,1,-1);
	  rescont.axpy(tmp5,wpgdet);

	  // Parte temporal + convectiva (Galerkin)
	  massm.prod(u_star,dshapex,-1,-1,1);
	  massm.axpy(SHAPE,1/(alpha*Dt));
	  matlocmom.prod(W_supg,massm,1,2).scale(rho);
        }

        if (comp_mat_res_th || comp_res_th) {

	  // thermal residue
	  // dT.set(T_star).rest(T).scale(rho*Cp);
	  dT = (T_star-T)/(rho*Cp);
	  // Guarda que pasa con esto? Aparentemente deberia ser la
	  // linea de arriba
#if 0   // Este es el corregido por mi
	  // tmp1_th = dT/(alpha*Dt)-Q_body;
	  tmp1_th.set(dT).scale(1./(alpha*Dt)).add(-Q_body);
#else  // Este es el original (Beto). No compila
	  // (esta escrito tipo Newmat)
	  //tmp1_th.set(dT).scale(1/(alpha*Dt)).rest(Q_body);
	  tmp1_th = dT/(alpha*Dt)-Q_body;
#endif
	  //resther.axpy(SHAPE,-wpgdet*tmp1_th);
	  //resther.axpy(P_supg_th,-wpgdet*tmp1_th);
	  resther.axpy(W_supg_th,-wpgdet*tmp1_th);

          if (weak_form) {
            tmp3_th.set(u_star).scale(rho*Cp*T_star);
            tmp4_th.prod(tmp3_th,dshapex,-1,-1,1);
	    resther.axpy(tmp4_th,wpgdet);
	    tmp20.prod(u_star,grad_T_star,-1,-1);
	    tmp2_th = double(tmp20)*(rho*Cp);
	    resther.axpy(P_supg_th,-wpgdet*tmp2_th);
          } else {
	    tmp20.prod(u_star,grad_T_star,-1,-1);
	    tmp2_th = double(tmp20)*(rho*Cp);
	    resther.axpy(W_supg_th,-wpgdet*tmp2_th);
          }

	  tmp8_th.prod(dshapex,grad_T_star,-1,1,-1).scale(kappa);
	  resther.axpy(tmp8_th,-wpgdet);

	  // thermal matrix
	  // massm_th.set(0.).axpy(SHAPE,rho*Cp/(alpha*Dt));
	  massm_th.set(SHAPE).scale(rho*Cp/(alpha*Dt));
	  matlocther.prod(W_supg_th,massm_th,1,2);

          if (weak_form) {
            tmp3_th.set(u_star).scale(-rho*Cp);
            tmp5_th.prod(tmp3_th,SHAPE,1,2);
            tmp6_th.prod(dshapex,tmp5_th,-1,1,-1,2);
            matlocther.add(tmp6_th);
	    massm_th.prod(u_star,dshapex,-1,-1,1).scale(rho*Cp);
            tmp7_th.prod(P_supg_th,massm_th,1,2);
            matlocther.add(tmp7_th);
          } else {
	    massm_th.prod(u_star,dshapex,-1,-1,1).scale(rho*Cp);
            tmp7_th.prod(W_supg_th,massm_th,1,2);
            matlocther.add(tmp7_th);
          }

	  tmp9_th.prod(dshapex,dshapex,-1,1,-1,2);
	  matlocther.axpy(tmp9_th,kappa);

        }

        // adding flow and thermal coupling (Boussinesq)
        if (comp_mat_res && comp_mat_res_th) {

	  tmp10_th.set(gravity).scale(rho*betath);
          // Fixme .... W_supg eld SHAPE
          tmp11_th.prod(W_supg,tmp10_th,1,2);
          //tmp11_th.prod(SHAPE,tmp10_th,1,2);
          resmom.axpy(tmp11_th,-(T_star-Tinfty)*wpgdet);

	  tmp12_th.prod(tmp11_th,SHAPE,1,2,3).scale(wpgdet);
	  matlocf.is(2,1,ndim).ir(4,ndof).add(tmp12_th).rs();

	  tmp13_th.prod(P_pspg,tmp10_th,-1,1,-1);
          rescont.axpy(tmp13_th,(T_star-Tinfty)*wpgdet);

	  tmp14_th.prod(tmp13_th,SHAPE,1,2).scale(-wpgdet);
	  matlocf.ir(2,1+ndim).ir(4,ndof).add(tmp14_th).rs();

          if (weak_form) {
            tmp15_th.set(dshapex).scale(rho*Cp*T_star);
	    tmp16_th.prod(tmp15_th,SHAPE,3,1,2).scale(-wpgdet);
	    matlocf.is(4,1,ndim).ir(2,ndof).add(tmp16_th).rs();
         
            tmp17_th.set(grad_T_star).scale(rho*Cp);
	    tmp18_th.prod(tmp17_th,SHAPE,1,2);
	    tmp19_th.prod(P_supg,tmp18_th,1,3,2).scale(wpgdet);
	    matlocf.is(4,1,ndim).ir(2,ndof).add(tmp19_th).rs();
            
          } else {
            tmp17_th.set(grad_T_star).scale(rho*Cp);
	    tmp18_th.prod(tmp17_th,SHAPE,1,2);
	    tmp19_th.prod(W_supg,tmp18_th,1,3,2).scale(wpgdet);
	    matlocf.is(4,1,ndim).ir(2,ndof).add(tmp19_th).rs();
          }
          
        }

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

        if (comp_mat_res || comp_res) {

	// Parte difusiva
	tmp7.prod(dshapex,dshapex,-1,1,-1,2);
	matlocmom.axpy(tmp7,nu_eff);

	// dmatw =  rho * ((1/Dt)*SHAPE + u * dshapex);
	dmatw.set(massm).scale(rho);

	// assert(delta_supg==0); // para no poner shock-capturing despues
	tmp13.prod(P_pspg,dshapex,-1,1,-1,2);
	// W_pspg.set(SHAPE).add(P_pspg);

	if (update_jacobian) {
	  for (int iloc=1; iloc<=nel; iloc++) {
	    for (int jloc=1; jloc<=nel; jloc++) {
	      double c = wpgdet*matlocmom.get(iloc,jloc);
	      for (int ii=1; ii<=ndim; ii++) {
		matlocf.addel(c,iloc,ii,jloc,ii);
	      }
	    }
	  }

	  tmp16.prod(W_supg,dshapex,1,2,3).scale(wpgdet);
	  //matlocf.is(2,1,ndim).ir(4,ndof).add(tmp16).rs();
	  matlocf.is(2,1,ndim).ir(4,ndim+1).add(tmp16).rs();

	  //matlocf.ir(2,ndof).is(4,1,ndim);
	  matlocf.ir(2,ndim+1).is(4,1,ndim);
	  tmp17.prod(P_pspg,dmatw,3,1,2).scale(wpgdet);
	  matlocf.rest(tmp17);
	  tmp17.prod(dshapex,SHAPE,3,2,1).scale(wpgdet);
	  matlocf.rest(tmp17).rs();

	  //matlocf.ir(2,ndof).ir(4,ndof).axpy(tmp13,-wpgdet).rs();
	  matlocf.ir(2,ndim+1).ir(4,ndim+1).axpy(tmp13,-wpgdet).rs();

	  if (!cache_grad_div_u) {
	    tmp19.set(dshapex).scale((delta_supg*rho+nu_eff)*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,1,4,3);
	    matlocf.is(2,1,ndim).is(4,1,ndim).add(tmp18).rs();
	  } else {
	    grad_div_u_coef += (delta_supg*rho+nu_eff)*wpgdet;
	    if (!grad_div_u_was_cached) {
	      tmp18.prod(dshapex,dshapex,2,1,4,3);
	      grad_div_u.add(tmp18);
	    }
	  }
	} 
        } 
        if (comp_mat_res_th || comp_res_th) {
	  for (int iloc=1; iloc<=nel; iloc++) {
	    for (int jloc=1; jloc<=nel; jloc++) {
	      double c = wpgdet*matlocther.get(iloc,jloc);
	      matlocf.addel(c,iloc,ndof,jloc,ndof);
	    }
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
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }      

    if (comp_mat_res) {
      if (cache_grad_div_u) {
	grad_div_u_coef /= double(npg);
	matlocf.is(2,1,ndim).is(4,1,ndim)
	  .axpy(grad_div_u,grad_div_u_coef).rs();
	if (!grad_div_u_was_cached) {
	  grad_div_u.export_vals(grad_div_u_cache);
	}
      }
      veccontr.is(2,1,ndim).set(resmom)
	.rs().ir(2,ndim+1).set(rescont).rs();

      if (comp_mat_res_th) {
	veccontr.ir(2,ndof).set(resther).rs();
      }
      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      if (update_jacobian) matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
