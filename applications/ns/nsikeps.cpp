//__INSERT_LICENSE__
/* $Id: nsikeps.cpp,v 1.27 2003/07/26 00:57:56 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"

#define ADD_GRAD_DIV_U_TERM
#define STANDARD_UPWIND
#define USE_FASTMAT

extern TextHashTable *GLOBAL_OPTIONS;

#define STOP {PetscFinalize(); exit(0);}
  
#define MAXPROP 100

//#define SQ(n) ((n)*(n))
#define SQ(n) square(n)

/** Cutoff function. It is very near to ${\rm ctff(x)}\approx \rm tol$ for
    $x<0$ and ${\rm ctff}(x)=x$ for $x\gg \rm tol$. 
*/ 
double ctff(double x, double & diff_ctff, double tol) {
  double r=x/tol-1.; 
  double ee,vaux,ret;
  if (fabs(r)<1e-7) {
    ret = (1.+0.5*exp(r)/(1+(1./6.)*r*r))*tol;
    ee  = tol*tol;
//    vaux = 7.0*ee+x*x-2.0*x*tol;
//    diff_ctff = 3.0*ee*exp(r)*(9.0*ee+x*x-4.0*x*tol)/vaux/vaux;
    vaux = (1.+1./6.*r*r);
    diff_ctff  = 0.5*exp(r)*(vaux-1./6.*2*r)/vaux/vaux;
    // dfx1dx  = 0.5*exp(r).*(1+1/6*r.^2-1/6*2*r)./(1+1/6*r.^2).^2;
  } else if (r>0) {
    ret =  (x-tol)/(1.-exp(-2.*r))+tol;
    vaux  = exp(-2.*r); 
    diff_ctff = 1.0/(1.0-vaux)*(1.0-2.*r*vaux/(1.0-vaux));
    // dfx21dx = (1-2*r.*exp(-2*r)-exp(-2*r))./(1-exp(-2*r)).^2;
  } else {
    ee = exp(2.*r);
    ret = (x-tol)*ee/(ee-1.)+tol;
    vaux = ee-1.0;
    diff_ctff = ee/vaux*(1.0-2.0*r/vaux);
    // dfx22dx = (exp(2*r)./(exp(2*r)-1).^2).*(exp(2*r)-1-2*r);
    // printf("ctff(%g,%g) = %g\n",x,tol,ret);
  }
  return ret;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// modif nsi_tet
#undef __FUNC__
#define __FUNC__ "nsi_tet_les_fm2::assemble"
int nsi_tet_keps::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			       Dofmap *dofmap,const char *jobinfo,int myrank,
			       int el_start,int el_last,int iter_mode,
			       const TimeData *time_) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(get_nearest_wall_element);

// added for kappa-epsilon equation
  GET_JOBINFO_FLAG(comp_mat_ke);
  GET_JOBINFO_FLAG(comp_mat_res_ke);
  GET_JOBINFO_FLAG(comp_res_ke);

  comp_mat_res_ke=comp_mat_res;
  comp_mat_ke=comp_mat;

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

  int ierr=0, axi;
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsikeps\n");

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
  if (comp_mat || comp_mat_ke) {
    retvalmat = arg_data_v[0].retval;
  } else if (get_nearest_wall_element) {
    wall_data = (WallData *)arg_data_v[0].user_data;
  }

  GlobParam *glob_param;
  double *hmin,rec_Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res || comp_mat_res_ke) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    if (update_jacobian) retvalmat = arg_data_v[ja++].retval;
    hmin = &*(arg_data_v[ja++].vector_assoc)->begin();
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
    wall_data = (WallData *)arg_data_v[ja++].user_data;
  } 

  //o Use a weak form for the gradient of pressure term.
  SGETOPTDEF(int,weak_form,1);
  //o Add shock-capturing term.
  SGETOPTDEF(double,shock_capturing_factor,1);

  // allocate local vecs
  int kdof;
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof), 
         locstate2(nel,ndof),xpg,G_body(ndim),gravity(ndim);

  if (ndof != ndim+3) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+3\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FMatrix matloc(nen,nen), matlocmom(nel,nel), masspg(nel,nel),
    grad_u_ext(ndof,ndof);
  FastMat2 matlocf(4,nel,ndof,nel,ndof);

  grad_u_ext.set(0.);

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

  //o Do not add turbulent viscosity for the momentum eq. (for
  //   debugging). 
  SGETOPTDEF(double,turbulence_coef,1.);
  //o Mask to the production terms in the k and epsilon
  // equations. This terms are then scaled by
  // \verb+turbulence_coef*turb_prod_coef+. 
  SGETOPTDEF(double,turb_prod_coef,1.);
  //o Scales the term proportional to $h/u$ in the pspg stabilization term
  SGETOPTDEF(double,pspg_advection_factor,1.);
  //o Scales the PSPG stabilization term 
  SGETOPTDEF(double,pspg_factor,1.);
  //o Density
  SGETOPTDEF(double,rho,1.);
  //o C_mu
  SGETOPTDEF(double,C_mu,0.09);
  //o C_1
  SGETOPTDEF(double,C_1,1.44);
  //o C_2
  SGETOPTDEF(double,C_2,1.92);
  //o sigma_k 
  SGETOPTDEF(double,sigma_k,1.);
  //o sigma_e 
  SGETOPTDEF(double,sigma_e,1.3);
  //o Cutoff value for $k$
  SGETOPTDEF(double,kap_ctff_val,1e-20);
  //o Cutoff value for $\epsilon$
  SGETOPTDEF(double,eps_ctff_val,1e-20);

  //o Add axisymmetric version for this particular elemset.
  TGETOPTDEF_S(thash,string,axisymmetric,none);
  assert(axisymmetric.length()>0);
  if (axisymmetric=="none") axi=0;
  else if (axisymmetric=="x") axi=1;
  else if (axisymmetric=="y") axi=2;
  else if (axisymmetric=="z") axi=3;
  else {
    PetscPrintf(PETSC_COMM_WORLD,
		"Invalid value for \"axisymmetric\" option\n"
		"axisymmetric=\"%s\"\n",axisymmetric.c_str());
    PetscFinalize();
    exit(0);
  }

  //o Add LES for this particular elemset.
  GGETOPTDEF(int,LES,0);
  //o Cache \verb+grad_div_u+ matrix
  SGETOPTDEF(int,cache_grad_div_u,0);
  // alpha is taken from the global parameters
  double &alpha = glob_param->alpha;
  //o Scale the SUPG upwind term. 
  SGETOPTDEF(double,tau_fac,1.);  // Scale upwind
  //o Adjust the stability parameters, taking into account
  // the time step. 
  SGETOPTDEF(double,temporal_stability_factor,0.);  // Scale upwind

  //o _T: double[ndim] _N: G_body _D: null vector 
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"G_body",G_body.storage_begin(),1,ndim);

  double pi = 4*atan(1.0);

  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))

  int nprops=iprop;
  
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  
  //GPdata gp_data(geom,ndim,nel,npg);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco, UU, u2, Peclet, psi, tau_supg, tau_pspg, div_u_star,
    p_star,wpgdet,velmod,tol,h_supg,fz,delta_supg,Uh;

  FastMat2 P_supg, W_supg, W_supg_t, dmatw, 
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
    wall_coords(ndim),dist_to_wall,tmp16,tmp162,tmp17,tmp18,tmp19,tmp20;

  double tau_supg_k,kap,kap_star,dkap,tmp1_ke,tmp2_ke;
  double tau_supg_e,eps,eps_star,deps;
  FastMat2 W_supg_k,P_supg_k,W_supg_e,P_supg_e;

  FMatrix kapcol,kapcol_new,kapcol_star,grad_kap_star(ndim),massm_ke;
  FMatrix epscol,epscol_new,epscol_star,grad_eps_star(ndim);

  FMatrix tmp3_ke,tmp4_ke,tmp5_ke,tmp6_ke,tmp7_ke,
    tmp8_ke,tmp9_ke,tmp71_ke,tmp81_ke,tmp10_ke,tmp11_ke;

  FMatrix reskap(nel),reseps(nel);

  double diff_coe_kap,diff_coe_eps,Pkap,Peps,Peps_2,eps_over_kap,kap_over_eps;
  double strain_rate_scalar,Jaco_kk,Jaco_ke,Jaco_ek,Jaco_ee;
  FMatrix Jaco_k(2),Jaco_e(2);

  double tmp12;
  double nu_eff;
  double tsf = temporal_stability_factor;
  double mix_length_inv,delta_supg_k,delta_supg_e;
  double fk,fe,fv,eps_before_ctff,dfkdk,dfvdv,dfede,dfedk,vaux,GG,dflidli;
  double lambda_1,lambda_2,lambda_max;

  FMatrix eye(ndim,ndim),seed,one_nel,matloc_prof(nen,nen);;

  FMatrix Jaco_axi(2,2),u_axi;
  int ind_axi_1, ind_axi_2;
  double detJaco_axi;
         
  if (axi) assert(ndim==3);

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

  double slope = 1./2./C_2*(3.*C_2-C_1+sqrt(pow(C_1,2.) - 
                 10.*C_1*C_2+9.*pow(C_2,2.)));

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
#ifdef USE_ANN
      xc.sum(xloc,-1,1).scale(1./double(nel));
      int nn;
      wall_data->nearest(xc.storage_begin(),nn);
      NN_IDX(k) = nn;
      continue;
#else
    PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
    }

    double grad_div_u_coef=0.;	// multiplies grad_div_u term
    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_mat_res || comp_res || comp_mat_res_ke || comp_res_ke) {
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
    matlocmom.set(0.);
    veccontr.set(0.);
    resmom.set(0.);
    rescont.set(0.);
    reskap.set(0.);
    reseps.set(0.);

    if (comp_mat_res && cache_grad_div_u) {
      if (grad_div_u_was_cached) {
	grad_div_u.set(grad_div_u_cache);
      } else {
	grad_div_u.set(0.);
      }
    }

    if (comp_mat_res) {
      ucols.set(locstate2.is(2,1,ndim));
      pcol.set(locstate2.rs().ir(2,ndim+1));
      kapcol.set(locstate2.rs().ir(2,ndof-1));
      epscol.set(locstate2.rs().ir(2,ndof));
      locstate2.rs();

      ucols_new.set(locstate.is(2,1,ndim));
      pcol_new.set(locstate.rs().ir(2,ndim+1));
      kapcol_new.set(locstate.rs().ir(2,ndof-1));
      epscol_new.set(locstate.rs().ir(2,ndof));
      locstate.rs();
      
      ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
      pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
      kapcol_star.set(kapcol_new).scale(alpha).axpy(kapcol,1-alpha);
      epscol_star.set(epscol_new).scale(alpha).axpy(epscol,1-alpha);
    }
    
    double shear_vel;
    int wall_elem;
    if (LES && comp_mat_res) {
#ifdef USE_ANN
      Elemset *wall_elemset;
      const double *wall_coords_;
      wall_data->nearest_elem_info(NN_IDX(k),wall_elemset,wall_elem,wall_coords_);
      wall_coords.set(wall_coords_);
      shear_vel = wall_elemset->elemprops_add[wall_elem];
#else
      PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
    }

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])
#define WPG_SUM  (gp_data.wpg_sum)

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

      // Modificado x Beto
      // double Area   = npg*wpgdet;
      double Area = detJaco*WPG_SUM;
      // fin modificado x Beto

      double wpgdet_c = wpgdet * turbulence_coef * turb_prod_coef;
      // Combined factor that affects the turbulence
      // production (for debugging)
      double tpf = turbulence_coef*turb_prod_coef;
      double h_pspg,Delta;
      if (ndim==2) {
	h_pspg = sqrt(4.*Area/pi);
	Delta = sqrt(Area);
      } else if (ndim==3 && axi==0) {
	// h_pspg = pow(6*Area/pi,1./3.);
	// El pow() da segmentation violation cuando corro con -O !!
	h_pspg = cbrt(6*Area/pi);
	Delta = cbrt(Area);
      } else if (ndim==3 && axi>0) {
        ind_axi_1 = (  axi   % 3)+1;
        ind_axi_2 = ((axi+1) % 3)+1;

	//        Jaco.is(1,ind_axi_1,ind_axi_2).is(2,ind_axi_1,ind_axi_2);
	//        Jaco_axi.set(Jaco);
	//        Jaco.rs();

        Jaco_axi.setel(Jaco.get(ind_axi_1,ind_axi_1),1,1);
        Jaco_axi.setel(Jaco.get(ind_axi_1,ind_axi_2),1,2);
        Jaco_axi.setel(Jaco.get(ind_axi_2,ind_axi_1),2,1);
        Jaco_axi.setel(Jaco.get(ind_axi_2,ind_axi_2),2,2);

        detJaco_axi = Jaco_axi.det();
	// Modificado x Beto July 9 2003
	//        double wpgdet_axi = detJaco_axi*WPG;
	//        double Area_axi = 0.5*npg*fabs(wpgdet_axi);
        double Area_axi = 0.5*detJaco_axi*WPG_SUM;
	h_pspg = sqrt(4.*Area_axi/pi);
	Delta = sqrt(Area_axi);
	// fin modificado x Beto

        } else {
	PFEMERRQ("Only dimensions 2 and 3 allowed for this element.\n");
      }
      
      if (comp_mat_res || comp_mat_res_ke) {
	// computes the minimum size of the mesh
	if (!WAS_SET || h_pspg<*hmin) {
	  WAS_SET = 1;
	  *hmin = h_pspg;
	}

	// state variables and gradient
	u.prod(SHAPE,ucols,-1,-1,1);
	kap = double(tmp8.prod(SHAPE,kapcol,-1,-1));
	eps = double(tmp8.prod(SHAPE,epscol,-1,-1));
	grad_u.prod(dshapex,ucols,1,-1,-1,2);

	p_star = double(tmp8.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);
	kap_star = double(tmp8.prod(SHAPE,kapcol_star,-1,-1));
	eps_star = double(tmp8.prod(SHAPE,epscol_star,-1,-1));

	grad_u_star.prod(dshapex,ucols_star,1,-1,-1,2);
	grad_p_star.prod(dshapex,pcol_star,1,-1,-1);
	grad_kap_star.prod(dshapex,kapcol_star,1,-1,-1);
	grad_eps_star.prod(dshapex,epscol_star,1,-1,-1);

        // cut off of kappa & epsilon
	kap = ctff(kap,dfkdk,kap_ctff_val);
	eps = ctff(eps,dfede,eps_ctff_val);

 	strain_rate.set(grad_u_star);
	grad_u_star.t();
	strain_rate.add(grad_u_star).scale(0.5);
        grad_u_star.rs();

        strain_rate_scalar = strain_rate.sum_square_all();

	kap_star = ctff(kap_star,dfkdk,kap_ctff_val);
	eps_star = ctff(eps_star,dfede,eps_ctff_val);

	// if (dfvdv != 1. | dfkdk != 1.) {
        //  printf(" cutoff acting \n");
        // }

        double nu_t = C_mu*kap*kap/eps;
	nu_eff = VISC + turbulence_coef*nu_t;

	//	uintri.prod(iJaco,u,1,-1,-1);
	//	Uh = uintri.sum_square_all();
	//	Uh = sqrt(Uh)/2;

	u2 = u.sum_square_all();

        if(axi>0){
          u_axi.set(u);
          u_axi.setel(0.,axi);
          u2 = u_axi.sum_square_all();
        }

	velmod = sqrt(u2);
        tol=1.0e-16;
        h_supg=0;
	FastMat2::branch();
        if(velmod>tol) {
	  FastMat2::choose(0);
	  // svec:= a streamline oriented unit vector
	  //	  svec.set(u).scale(1./velmod);
	  if(axi>0){
            svec.set(u_axi).scale(1./velmod);
          } else {
            svec.set(u).scale(1./velmod);
          }
	  h_supg = tmp9.prod(dshapex,svec,-1,1,-1).sum_abs_all();
          h_supg = (h_supg < tol ? tol : h_supg);
          h_supg = 2./h_supg;
        } else {
          h_supg = h_pspg;
        }
	FastMat2::leave();

        if (comp_mat_res || comp_res) {
	  Peclet = velmod * h_supg / (2. * nu_eff);
          tau_supg = tsf*SQ(2.*rec_Dt)+SQ(2.*velmod/h_supg)
	    +9.*SQ(4.*nu_eff/SQ(h_supg));
          tau_supg = 1./sqrt(tau_supg);

          tau_pspg = tsf*SQ(2.*rec_Dt)
	    +SQ(pspg_advection_factor*2.*velmod/h_pspg)
	    +9.*SQ(4.*nu_eff/SQ(h_pspg));
          tau_pspg = pspg_factor/sqrt(tau_pspg);

          fz = (Peclet < 3. ? Peclet/3. : 1.);
          delta_supg = 0.5*h_supg*velmod*fz;
	
	  if (tau_fac != 1.) {
	    tau_pspg *= tau_fac;
	    tau_supg *= tau_fac;
	  }
	  delta_supg *= shock_capturing_factor;

        } 
        if (comp_mat_res_ke || comp_res_ke) {
	  diff_coe_kap=VISC+turbulence_coef*nu_t/sigma_k;
          diff_coe_kap = (diff_coe_kap < tol ? tol : diff_coe_kap);
	  diff_coe_eps=VISC+turbulence_coef*nu_t/sigma_e;
          diff_coe_eps = (diff_coe_eps < tol ? tol : diff_coe_eps);
	  FastMat2::branch();
	  if(velmod>tol) {
	    FastMat2::choose(0);
	    Peclet = velmod * h_supg / (2. * diff_coe_kap);
	    psi = 1./tanh(Peclet)-1/Peclet;
	    tau_supg_k = psi*h_supg/(2.*velmod);
            delta_supg_k = 0.5*h_supg*velmod*psi;

	    Peclet = velmod * h_supg / (2. * diff_coe_eps);
	    psi = 1./tanh(Peclet)-1/Peclet;
	    tau_supg_e = psi*h_supg/(2.*velmod);
            delta_supg_e = 0.5*h_supg*velmod*psi;

	  } else {
	    FastMat2::choose(1);
	    tau_supg_k = 0;
	    tau_supg_e = 0;
            delta_supg_k = 0;
            delta_supg_e = 0;
	  }
	  FastMat2::leave();
        }

        if (comp_mat_res || comp_res) {	
	  // P_supg es un vector fila
	  P_supg.prod(u,dshapex,-1,-1,1).scale(tau_supg);

	  // Weight function 
	  W_supg.set(P_supg).add(SHAPE);

	  // Pressure stabilizing term
	  P_pspg.set(dshapex).scale(tau_pspg/rho);
        }

        if (comp_mat_res_ke || comp_res_ke) {	

	  // P_supg stabilizing term for turbulent transport equations
	  P_supg_k.prod(u,dshapex,-1,-1,1).scale(tau_supg_k);
	  W_supg_k.set(P_supg_k).add(SHAPE);

	  P_supg_e.prod(u,dshapex,-1,-1,1).scale(tau_supg_e);
	  W_supg_e.set(P_supg_e).add(SHAPE);

        }

	// implicit version - General Trapezoidal rule - parameter alpha

        if (comp_mat_res || comp_res) {
	  // dmatu := material derivative of u, also
	  // includes G_body, i.e. the external force field. 
#ifdef ADD_GRAD_DIV_U_TERM
	  dmatu.prod(u_star,grad_u_star,-1,-1,1);
#else
	  dmatu.prod(u,grad_u_star,-1,-1,1);
#endif

	  du.set(u_star).rest(u);
	  dmatu.axpy(du,rec_Dt/alpha).rest(G_body);
	
	  div_u_star = double(tmp10.prod(dshapex,ucols_star,-1,-2,-2,-1));

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
#ifdef ADD_GRAD_DIV_U_TERM
	  massm.prod(u_star,dshapex,-1,-1,1);
#else
	  massm.prod(u,dshapex,-1,-1,1);
#endif
	  massm.axpy(SHAPE,rec_Dt/alpha);
	  matlocmom.prod(W_supg,massm,1,2).scale(rho);
        }

        if (comp_mat_res_ke || comp_res_ke) {

	  // kappa-epsilon residue

          // temporal terms
 #if 0
	  dkap = (kap_star-kap);
	  tmp1_ke = dkap*rec_Dt/alpha;
	  reskap.axpy(W_supg_k,-wpgdet*tmp1_ke);

	  deps = (eps_star-eps);
	  tmp1_ke = deps*rec_Dt/alpha;
	  reseps.axpy(W_supg_e,-wpgdet*tmp1_ke);
#endif

          // convective terms 
	  tmp20.prod(u_star,grad_kap_star,-1,-1);
	  tmp2_ke = double(tmp20);
	  reskap.axpy(W_supg_k,-wpgdet*tmp2_ke);

	  tmp20.prod(u_star,grad_eps_star,-1,-1);
	  tmp2_ke = double(tmp20);
	  reseps.axpy(W_supg_e,-wpgdet*tmp2_ke);

          // diffusive terms
	  tmp3_ke.prod(dshapex,grad_kap_star,-1,1,-1).scale(diff_coe_kap);
	  reskap.axpy(tmp3_ke,-wpgdet);
	  tmp3_ke.prod(dshapex,grad_eps_star,-1,1,-1).scale(diff_coe_eps);
	  reseps.axpy(tmp3_ke,-wpgdet);
 

        // adding production terms to the turbulent transport equations

          //kap_over_eps = (abs(eps_star) < tol ? kap_star/tol : kap_star/eps_star);
          //eps_over_kap = (abs(kap_star) < tol ? eps_star/tol : eps_star/kap_star);

          eps_over_kap = eps_star/kap_star;
          kap_over_eps = kap_star/eps_star;
       
          Pkap = 2.*strain_rate_scalar*C_mu*kap_star*kap_over_eps;
          Peps = 2.*C_1*C_mu*kap_star*strain_rate_scalar;
          Peps_2 = C_2*eps_over_kap*eps_star;

	  // kappa-epsilon matrix

          // using cutoff derivatives
          GG = 2.0*C_mu*strain_rate_scalar;
          fk = kap_star;
          fe = eps_star;
          Jaco_kk = rec_Dt/alpha -  tpf*GG*2.0*fk/fe*dfkdk; 
          Jaco_ke = -tpf*(-GG*pow((fk/fe),2)-1.0)*dfede;      
          Jaco_ek = -tpf*(C_1*GG+C_2*pow((fe/fk),2.))*dfkdk;
          Jaco_ee = rec_Dt/alpha-tpf*(-2.0*C_2*fe/fk*dfede);

          Jaco_k.ir(1,1).set(Jaco_kk).rs();
          Jaco_k.ir(1,2).set(Jaco_ke).rs();
          Jaco_e.ir(1,1).set(Jaco_ek).rs();
          Jaco_e.ir(1,2).set(Jaco_ee).rs();

          // consistent matrix
          tmp7_ke.prod(W_supg_k,SHAPE,1,2);
          tmp8_ke.prod(W_supg_e,SHAPE,1,2);
          //tmp7_ke.prod(SHAPE,SHAPE,1,2);
          //tmp8_ke.prod(SHAPE,SHAPE,1,2);

          tmp9_ke.prod(tmp7_ke,Jaco_k,1,2,3);
          matlocf.ir(2,ndof-1).is(4,ndof-1,ndof).axpy(tmp9_ke,wpgdet).rs();
          
#if 1
          tmp6_ke.set(W_supg_k).scale((Pkap-eps_star));
          reskap.axpy(tmp6_ke,wpgdet_c);

          // adding temporal terms here due to the lumped mass matrix option
	  tmp10_ke.set(kapcol_star).rest(kapcol);
	  tmp11_ke.prod(tmp7_ke,tmp10_ke,1,-1,-1);
	  reskap.axpy(tmp11_ke,-wpgdet*rec_Dt/alpha);
#else
	  tmp10_ke.set(kapcol_star).scale(Jaco_kk);
          tmp10_ke.axpy(epscol_star,Jaco_ke);
	  tmp6_ke.prod(tmp7_ke,tmp10_ke,1,-1,-1);
          reskap.axpy(tmp6_ke,-wpgdet);
          
	  tmp10_ke.set(kapcol).scale(-1.);
	  tmp11_ke.prod(tmp7_ke,tmp10_ke,1,-1,-1);
	  reskap.axpy(tmp11_ke,-wpgdet*rec_Dt/alpha);
#endif

          tmp9_ke.prod(tmp8_ke,Jaco_e,1,2,3);
          matlocf.ir(2,ndof).is(4,ndof-1,ndof).axpy(tmp9_ke,wpgdet).rs();
#if 1
          tmp6_ke.set(W_supg_e).scale(Peps-Peps_2);
          reseps.axpy(tmp6_ke,wpgdet_c);

	  tmp10_ke.set(epscol_star).rest(epscol);
	  tmp11_ke.prod(tmp8_ke,tmp10_ke,1,-1,-1);
	  reseps.axpy(tmp11_ke,-wpgdet*rec_Dt/alpha);
#else
	  tmp10_ke.set(kapcol_star).scale(Jaco_ek);
          tmp10_ke.axpy(epscol_star,Jaco_ee);
	  tmp6_ke.prod(tmp8_ke,tmp10_ke,1,-1,-1);
          reseps.axpy(tmp6_ke,-wpgdet);

	  tmp10_ke.set(epscol).scale(-1.);
	  tmp11_ke.prod(tmp8_ke,tmp10_ke,1,-1,-1);
	  reseps.axpy(tmp11_ke,-wpgdet*rec_Dt/alpha);
#endif

          // convective terms
	  massm_ke.prod(u_star,dshapex,-1,-1,1);
          tmp4_ke.prod(W_supg_k,massm_ke,1,2);
          matlocf.ir(2,ndof-1).ir(4,ndof-1).axpy(tmp4_ke,wpgdet).rs();

          tmp4_ke.prod(W_supg_e,massm_ke,1,2);
          matlocf.ir(2,ndof).ir(4,ndof).axpy(tmp4_ke,wpgdet).rs();

          // diffusive terms
	  tmp5_ke.prod(dshapex,dshapex,-1,1,-1,2);
          matlocf.ir(2,ndof-1).ir(4,ndof-1).axpy(tmp5_ke,wpgdet*diff_coe_kap).rs();
          matlocf.ir(2,ndof).ir(4,ndof).axpy(tmp5_ke,wpgdet*diff_coe_eps).rs();

          // shock capturing terms

          double aa = Jaco_kk/diff_coe_kap;
          double bb = Jaco_ke/diff_coe_kap;
          double cc = Jaco_ek/diff_coe_eps;
          double dd = Jaco_ee/diff_coe_eps;

          double discri = pow(0.5*(aa+dd),2.)-(aa*dd-cc*bb);

          //if(discri<0) {
          //   printf("discriminante negativo \n");
          //}

          if(discri<0) {
            lambda_1 = 0.5*(aa+dd);
            lambda_max = ( lambda_1 > tol ? lambda_1 : tol );
          } else {
            lambda_1 = 0.5*(aa+dd)+sqrt(discri);
            lambda_2 = 0.5*(aa+dd)-sqrt(discri);
            lambda_max = ( lambda_1 > lambda_2 ? lambda_1 : lambda_2 );
            lambda_max = ( lambda_max > tol ? lambda_max : tol );
          };

          lambda_max = lambda_max*pow(h_pspg,2.);
	  psi = 1./tanh(lambda_max)-1./lambda_max;
          double coef_delta = 1.0;
          delta_supg_k = delta_supg_k + coef_delta*lambda_max*psi*diff_coe_kap;
          delta_supg_e = delta_supg_e + coef_delta*lambda_max*psi*diff_coe_eps;

	  tmp3_ke.prod(dshapex,grad_kap_star,-1,1,-1).scale(delta_supg_k);
	  reskap.axpy(tmp3_ke,-wpgdet);
	  tmp3_ke.prod(dshapex,grad_eps_star,-1,1,-1).scale(delta_supg_e);
	  reseps.axpy(tmp3_ke,-wpgdet);

          matlocf.ir(2,ndof-1).ir(4,ndof-1).axpy(tmp5_ke,wpgdet*delta_supg_k).rs();
          matlocf.ir(2,ndof).ir(4,ndof).axpy(tmp5_ke,wpgdet*delta_supg_e).rs();

#if 0
       // acoplamiento de las ecs de transporte turbulento con las de momento
       // por ahora es solo una idea preliminar

          tmp7_ke.set(grad_kap_star);
	  tmp8_ke.prod(tmp7_ke,SHAPE,1,2);
	  tmp9_ke.prod(W_supg_k,tmp8_ke,1,3,2).scale(wpgdet);
	  matlocf.is(4,1,ndim).ir(2,ndof-1).add(tmp9_ke).rs();

          tmp7_ke.set(grad_eps_star);
	  tmp8_ke.prod(tmp7_ke,SHAPE,1,2);
	  tmp9_ke.prod(W_supg_e,tmp8_ke,1,3,2).scale(wpgdet);
	  matlocf.is(4,1,ndim).ir(2,ndof).add(tmp9_ke).rs();
#endif          
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

// old version but wrong
	 //  tmp16.prod(W_supg,dshapex,1,2,3).scale(wpgdet);
	 //  matlocf.is(2,1,ndim).ir(4,ndim+1).add(tmp16).rs();
 
// new version (Mario) I hope it is OK
        if (weak_form) {
           tmp16.prod(P_supg,dshapex,1,2,3).scale(wpgdet);
           tmp162.prod(dshapex,SHAPE,2,1,3).scale(-wpgdet);
           matlocf.is(2,1,ndim).ir(4,ndim+1).add(tmp16)
                                            .add(tmp162).rs();
       } else {
           tmp16.prod(W_supg,dshapex,1,2,3).scale(wpgdet);
           matlocf.is(2,1,ndim).ir(4,ndim+1).add(tmp16).rs();
       }
 
	  //matlocf.ir(2,ndof).is(4,1,ndim);
	  matlocf.ir(2,ndim+1).is(4,1,ndim);
	  tmp17.prod(P_pspg,dmatw,3,1,2).scale(wpgdet);
	  matlocf.rest(tmp17);
	  tmp17.prod(dshapex,SHAPE,3,2,1).scale(wpgdet);
	  matlocf.rest(tmp17).rs();

	  //matlocf.ir(2,ndof).ir(4,ndof).axpy(tmp13,-wpgdet).rs();
	  matlocf.ir(2,ndim+1).ir(4,ndim+1).axpy(tmp13,-wpgdet).rs();

	  if (!cache_grad_div_u) {
	    tmp19.set(dshapex).scale(nu_eff*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,3,4,1);
	    matlocf.is(2,1,ndim).is(4,1,ndim).add(tmp18).rs();
	    tmp19.set(dshapex).scale(delta_supg*rho*wpgdet);
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

      if (comp_mat_res_ke) {
	veccontr.ir(2,ndof-1).set(reskap).rs();
	veccontr.ir(2,ndof).set(reseps).rs();
      }
      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      // matlocf.set(0.);
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
