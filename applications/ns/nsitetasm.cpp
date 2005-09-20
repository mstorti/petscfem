//__INSERT_LICENSE__
//$Id: nsitetasm.cpp,v 1.1.2.1 2005/09/20 00:58:38 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "nsifunaux.h"

#define ADD_GRAD_DIV_U_TERM
#define STANDARD_UPWIND
#define USE_FASTMAT

extern TextHashTable *GLOBAL_OPTIONS;

#define STOP {PetscFinalize(); exit(0);}
   
#define MAXPROP 100

void vel_axi(FastMat2 &u,FastMat2 &u_axi, const int axi) {
  u_axi.set(u);
  if(axi>0){
    //  u_axi.set(u);
    u_axi.setel(0.,axi);
    //  u2 = u_axi.sum_square_all();
  }	
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
double compute_h_supg(FastMat2 &u_axi,FastMat2 &dshapex, double velmod, double h_pspg) {
  
  static FastMat2 tmp9,svec;
  double tol=1.0e-16, h_supg;
  
  h_supg=0;
  FastMat2::branch();
  if(velmod>tol) {
    FastMat2::choose(0);
    svec.set(u_axi).scale(1./velmod);
    h_supg = tmp9.prod(dshapex,svec,-1,1,-1).sum_abs_all();
    h_supg = (h_supg < tol ? tol : h_supg);
    h_supg = 2./h_supg;
  } else {
    h_supg = h_pspg;
  }
  FastMat2::leave();
  
  return h_supg;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double compute_rho_m(FastMat2 &rho_g_vp,FastMat2 &arho_g_vp, FastMat2 &alpha_g_vp, 
		     double &alpha_l, double rho_l, int nphases) {

  static FastMat2 Id_vp(2,nphases,nphases); 
  
  double alpha_d_sum = (double) alpha_g_vp.sum_all();
  alpha_l = 1.0 - alpha_d_sum;
  
  Id_vp.set(0.).d(1,2);
  Id_vp.set(rho_g_vp).rs();	
  arho_g_vp.prod(Id_vp,alpha_g_vp,1,-1,-1);
  double arho_l = rho_l*alpha_l;
  double rho_m = arho_l+arho_g_vp.sum_all();

  return rho_m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void compute_vel_g(const FastMat2 &v_m, FastMat2 &vslip_vp, const FastMat2 &vslip_user_vp, 
		   FastMat2 &vslip_m_vp, FastMat2 &v_g_vp, const FastMat2 &rho_g_vp, 
		   const double rho_m, const int g_dir, const FastMat2 &d_bubble_vp, 
		   const int nphases, const int use_modified_slag_vslip,
		   const FastMat2 &alpha_g_vp) {
  
  
  //  double vslip;
  
  if (vslip_user_vp.sum_abs_all()>0.0) {
    vslip_vp.set(vslip_user_vp);
  } else {
    for (int j=1; j<=nphases; j++) {
      double rb = d_bubble_vp.get(j);
      double vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
		      rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));
      
      vslip = (g_dir > 0 ? vslip : -vslip);
      vslip_vp.setel(vslip,j);
    }
  }
  vslip_m_vp.set(vslip_vp);
  // modifico velocidad slip de la escoria por la diferencia de densidades con la mezcla
  double factor,tmp;
  if (use_modified_slag_vslip) {
    if(1){      
      double rho_g = rho_g_vp.get(nphases);
      tmp = double(vslip_m_vp.get(nphases))*(1.-rho_g/rho_m);
    } else {
      double alpha_g = alpha_g_vp.get(nphases);
      if(0){
	factor = 1.0-tanh(alpha_g/0.3);
      } else {
	factor = pow((1.0-alpha_g),3.0);
      }
      tmp = double(vslip_m_vp.get(nphases))*factor;
    }
    vslip_m_vp.setel(tmp,nphases);
  }
  // velocidad del gas	
  for (int j=1; j<=nphases; j++) {
    v_g_vp.ir(2,j);
    v_g_vp.set(v_m).addel(vslip_m_vp.get(j),abs(g_dir));
    v_g_vp.rs();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "nsi_tet_asm::assemble"

int nsi_tet_asm::
assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
	 Dofmap *dofmap,const char *jobinfo,int myrank,
	 int el_start,int el_last,int iter_mode,
	 const TimeData *time_) {

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(get_nearest_wall_element);

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
  // PetscPrintf(PETSC_COMM_WORLD,"entrando a nsi_tet\n");

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)
#define NN_IDX(j) ELEMIPROPS_ADD(j,0)
#define IDENT(j,k) (ident[ndof*(j)+(k)]) 

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
  if (comp_mat) {
    retvalmat = arg_data_v[0].retval;
  } else if (get_nearest_wall_element) {
    wall_data = (WallData *)arg_data_v[0].user_data;
    if(!wall_data) {
      printf("Null 'wall_data' object found.\n");
      set_error(2);
      return 1;
    }
  }

  // rec_Dt is the reciprocal of Dt (i.e. 1/Dt)
  // for steady solutions it is set to 0. (Dt=inf)
  GlobParam *glob_param;
  double *hmin,Dt,rec_Dt;
  int ja_hmin;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res) {
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
  assert(weak_form==1);

  //o Solve coupled
  SGETOPTDEF(int,coupled,0);

  //o Add shock-capturing term.
  SGETOPTDEF(double,shock_capturing_factor,0.);

  //o Add shock-capturing term.
  SGETOPTDEF(double,shocap,0.);

  if(shocap>0.) shock_capturing_factor = shocap;

  //o Add pressure controlling term. 
  SGETOPTDEF(double,pressure_control_coef,0.);
  assert(pressure_control_coef>=0.);
  //o number of diferent phases
  SGETOPTDEF(int,nphases,1);
  assert(nphases>0);
  //o Sato coefficient
  SGETOPTDEF(double,Sato_model_coef,0.);
  // Schmidt number
  SGETOPTDEF(double,Sc_number,1.0);
  //o use_modified_slag_vslip : key to modify slag slip velocity 
  //				as a function of slag void fraction 
  SGETOPTDEF(int,use_modified_slag_vslip,0);

  // allocate local vecs
  int kdof;
  FastMat2 veccontr(2,nel,ndof),xloc(2,nel,ndim),locstate(2,nel,ndof), 
	 locstate2(2,nel,ndof),xpg,G_body(1,ndim);

  if (ndof != ndim+1+nphases) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1+nphases\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FMatrix matloc(nen,nen), matlocmom(nel,nel), masspg(nel,nel);
  FastMat2 matlocf(4,nel,ndof,nel,ndof);

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

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
  SGETOPTDEF(int,LES,0);
  //o Cache  #grad_div_u#  matrix
  SGETOPTDEF(int,cache_grad_div_u,0);
  //o Smagorinsky constant.
  SGETOPTDEF(double,C_smag,0.18); // Dijo Beto
  //o van Driest constant for the damping law.
  SGETOPTDEF(double,A_van_Driest,0); 
  assert(A_van_Driest>=0.);
  //o Scale the SUPG and PSPG stabilization term. 
  SGETOPTDEF(double,tau_fac,1.);  // Scale upwind
  //o Scales the PSPG stabilization term. 
  SGETOPTDEF(double,tau_pspg_fac,1.);  // Scale upwind
  //o Scale the residual term. 
  SGETOPTDEF(double,residual_factor,1.);
  //o Scale the jacobian term. 
  SGETOPTDEF(double,jacobian_factor,1.);
  //o Adjust the stability parameters, taking into account
  // the time step. If the  #steady#  option is in effect,
  // (which is equivalent to $\Dt=\infty$) then
  //  #temporal_stability_factor#  is set to 0.
  SGETOPTDEF(double,temporal_stability_factor,0.);  // Scale upwind
  if (comp_mat_res && glob_param->steady) temporal_stability_factor=0;

  //o Add to the  #tau_pspg#  term, so that you can stabilize with a term
  //  independently of $h$. (Mainly for debugging purposes). 
  SGETOPTDEF(double,additional_tau_pspg,0.);  // Scale upwind
  double &alpha = glob_param->alpha;

  //o _T: double[ndim] _N: G_body _D: null vector 
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"G_body",G_body.storage_begin(),1,ndim);

  double pi = 4*atan(1.0);

  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))

  int nprops=iprop;

  //o Density (should be changed for rho_m)
  TGETOPTDEF(thash,double,rho,1.);
  
  //o Density
  TGETOPTDEF(thash,double,rho_l,0.);
  assert(rho_l>0);

  //o viscosity
  TGETOPTDEF(thash,double,visco_l,VISC);

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  //GPdata gp_data(geom,ndim,nel,npg);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco, UU, u2, Peclet, psi, tau_supg, tau_pspg, div_u_star,
    p_star,wpgdet,velmod,tol,h_supg,fz,delta_supg,Uh;

  FastMat2 P_supg, W_supg, W_supg_t, dmatw,
    grad_div_u(4,nel,ndim,nel,ndim),P_pspg(2,ndim,nel),dshapex(2,ndim,nel);
  double *grad_div_u_cache;
  int grad_div_u_was_cached;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FMatrix Jaco(ndim,ndim),iJaco(ndim,ndim),
    grad_u(ndim,ndim),grad_u_star,strain_rate(ndim,ndim),resmom(nel,ndim),
    dresmom(nel,ndim),matij(ndof,ndof),Uintri,svec(ndim);

  FMatrix grad_p_star(ndim),u,u_star,du,
    uintri(ndim),rescont(nel),dmatu(ndim),ucols,ucols_new,
    ucols_star,pcol_star,pcol_new,pcol,fm_p_star,tmp1,tmp2,tmp3,tmp6,
    massm,tmp8,tmp9,tmp10,tmp11,tmp13,tmp14,tmp15,dshapex_c,xc,
    wall_coords(ndim),dist_to_wall,tmp16,tmp162,tmp17,tmp171,tmp172,
    tmp173,tmp174,tmp18,tmp19;

  FastMat2 tmp20(2,nel,nel),tmp21;

  FastMat2 vfcols,vfcols_new,vfcols_star,Id_vp(2,nphases,nphases);
  FastMat2 vf,vf_star,grad_vf,grad_vf_star,grad_rho_m(1,ndim);
  FastMat2 res_alpha_g(2,nel,nphases),du_g(1,nphases), //dmatu_g(1,nphases),
    dresmom_g(2,nel,nphases),matloc_alpha_g(4,nel,nphases,nel,nphases);
  FastMat2 tmp1_g,tmp11_g,tmp12_g,tmp2_g,tmp3_g,tmp4_g,tmp5_g,tmp6_g,tmp7_g,
    tmp9_g,tmp10_g,Ni_Nj,gNi_Nj;
  FastMat2 rescont_gal,rescont_pspg,resmom_prime,resmom_supg;
  FastMat2 gNi_gNj;

  double tmp12;
  double tsf = temporal_stability_factor;
  
  FMatrix eye(ndim,ndim),seed,one_nel,matloc_prof(nen,nen);
  
  FMatrix Jaco_axi(2,2),u_axi;
  int ind_axi_1, ind_axi_2;
  double detJaco_axi;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Definicion de parametros de las diferentes fases
  FastMat2 rho_g_vp,visco_g_vp,d_bubble_vp,vslip_user_vp,alpha_source_vp,visco_g_eff_vp;
  FastMat2 velmod_g_vp,h_supg_vp,P_supg_vp,u_axi_g_vp,delta_supg_vp;

  double tau_supg_g,h_supg_g,velmod_g,u2_g,delta_supg_g;

  //  phases density vector
  rho_g_vp.resize(1,nphases);
  // rho_g_vp.set(rho_g);
  ierr = get_double(GLOBAL_OPTIONS,"rho_phases",rho_g_vp.storage_begin(),1,nphases);

  //  phases viscosity vector
  visco_g_vp.resize(1,nphases);
  // visco_g_vp.set(visco_g);
  ierr = get_double(GLOBAL_OPTIONS,"visco_phases",visco_g_vp.storage_begin(),1,nphases);

  //  Bubble diameter vector
  d_bubble_vp.resize(1,nphases);
  d_bubble_vp.set(0.0);
  ierr = get_double(GLOBAL_OPTIONS,"d_bubble_phases",d_bubble_vp.storage_begin(),1,nphases);

  //  slip velocity vector
  vslip_user_vp.resize(1,nphases);
  vslip_user_vp.set(0.0);
  ierr = get_double(GLOBAL_OPTIONS,"vslip_user_phases",vslip_user_vp.storage_begin(),1,nphases);

  //  source term vector
  alpha_source_vp.resize(1,nphases);
  // alpha_source_vp.set(alpha_source);
  ierr = get_double(GLOBAL_OPTIONS,"alpha_source_phases",alpha_source_vp.storage_begin(),1,nphases);

  //o Direction of gravity
  TGETOPTDEF(thash,int,g_dir,ndim);

  double rho_g,vslip,rho_m,rho_m_old,arho_l,arho_g,vslip_m,alpha_l,alpha_g;
  double d_bubble,visco_m_eff,visco_t,visco_g,visco_g_eff;
  vector<int> alpha_indx_vp;
  int vl_indx = 1;
  int vl_indxe = vl_indx+ndim-1;
  int alpha_indx = ndim+2;
  alpha_indx_vp.resize(nphases);
  for (int j=0; j<nphases; j++) alpha_indx_vp[j] = vl_indxe+j+1;

  FastMat2 alpha_g_vp,arho_g_vp,v_m,v_mix,v_g,v_g_vp,v_g_vp_old;
  FastMat2 v_rel,v_rel_vp,vslip_vp,vslip_m_vp;

  alpha_g_vp.resize(1,nphases);
  arho_g_vp.resize(1,nphases);
  v_m.resize(1,ndim);
  v_mix.resize(1,ndim);
  v_g.resize(1,ndim);
  v_g_vp.resize(2,ndim,nphases);
  v_g_vp_old.resize(2,ndim,nphases);
  v_rel.resize(1,ndim);
  v_rel_vp.resize(2,ndim,nphases);
  vslip_vp.resize(1,nphases);
  visco_g_eff_vp.resize(1,nphases);
  velmod_g_vp.resize(1,nphases);
  h_supg_vp.resize(1,nphases);
  P_supg_vp.resize(2,nel,nphases);
  u_axi_g_vp.resize(2,ndim,nphases);
  delta_supg_vp.resize(1,nphases);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	 
  if (axi) assert(ndim==3);
  
  eye.eye();
  
  if (comp_mat) {
    
    if (coupled) {
      matloc_prof.set(1.);
    } else {      
      
      seed.resize(2,ndof,ndof);
      seed.set(0.);
      one_nel.resize(2,nel,nel);
      one_nel.set(0.);
      matloc_prof.resize(2,nen,nen);
      
#ifndef ADD_GRAD_DIV_U_TERM
      for (int jj=1; jj<=ndim; jj++) {
	seed.setel(1.,jj,jj);
	seed.setel(1.,jj,ndim+1);
	seed.setel(1.,ndim+1,jj);
      }
      seed.setel(1.,ndim+1,ndim+1);
#else
      seed.is(1,1,ndim+1).is(2,1,ndim+1).set(1.0).rs();
#endif

      // disperse phases
      for (int j=ndim+2; j<=ndim+1+nphases; j++) {
	seed.setel(1.,j,j);
      }
      one_nel.set(1.);
      matloc_prof.kron(one_nel,seed);
    }
    
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
      seed.setel(jj,ndim+1,1.);
      seed.setel(ndim+1,jj,1.);
    }
    seed.setel(ndim+1,ndim+1,1.);
    // disperse phases
    for (int j=ndim+2; j<=ndim+1+nphases; j++) {
      seed.setel(1.,j,j);
    }
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

    if (get_nearest_wall_element && A_van_Driest>0.) {
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
    //			 locstate  <- u^*
    if (comp_mat_res || comp_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));

      if (cache_grad_div_u) {
	grad_div_u_cache = (double *)local_store_address(k);
	grad_div_u_was_cached = (grad_div_u_cache!=NULL);
	if (!grad_div_u_was_cached) {
	  local_store_address(k) 
	    = grad_div_u_cache 
	    = new double[ndim*ndim*nel*nel];
	}
	//#define DEBUG_CACHE_GRAD_DIV_U
#ifdef	DEBUG_CACHE_GRAD_DIV_U	// debug:=
	if (k<2 && grad_div_u_was_cached) {
	  printf("element %d, cached grad_div_u: ",k);
	  for (int kkkk=0; kkkk<ndim*ndim*nel*nel; kkkk++)
	    printf("%f	",grad_div_u_cache[kkkk]);
	  printf("\n");
	}
	grad_div_u_was_cached = 0; // In order to recompute always
				   // the grad_div_u operator
#endif
      }
    }

    matlocmom.set(0.);
    matloc.set(0.);
    matlocf.set(0.);
    veccontr.set(0.);
    resmom.set(0.);
    rescont.set(0.);
    res_alpha_g.set(0.);
    matloc_alpha_g.set(0.);
    if (comp_mat_res && cache_grad_div_u) {
      if (grad_div_u_was_cached) {
	grad_div_u.set(grad_div_u_cache);
      } else {
	grad_div_u.set(0.);
      }
    }

    if (comp_res || comp_mat_res) {
      ucols.set(locstate2.is(2,1,ndim));
      pcol.set(locstate2.rs().ir(2,ndim+1));
      vfcols.set(locstate2.rs().is(2,ndim+2,ndim+1+nphases));
      locstate2.rs();
      
      ucols_new.set(locstate.is(2,1,ndim));
      pcol_new.set(locstate.rs().ir(2,ndim+1));
      vfcols_new.set(locstate.rs().is(2,ndim+2,ndim+1+nphases));
      locstate.rs();
      
      ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
      pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
      vfcols_star.set(vfcols_new).scale(alpha).axpy(vfcols,1-alpha);

      //#define PRINT_ELEM_DEBUG
#ifdef PRINT_ELEM_DEBUG
      if (k==0) {
	locstate2.print("locstate2 (t_n):");
	locstate.print("locstate (t_n+1):");
      }
#endif
    }
    
    double shear_vel;
    int wall_elem;
    if (LES && comp_mat_res && A_van_Driest>0.) {
#ifdef USE_ANN
      if (!wall_data) { set_error(2); return 1; }
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
#define SHAPE	 (*gp_data.FM2_shape[ipg])
#define WPG	 (gp_data.wpg[ipg])
#define WPG_SUM	 (gp_data.wpg_sum)

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);

      detJaco = Jaco.det();
      if (detJaco<=0.) {
	detj_error(detJaco,elem);
	set_error(1);
      }
      wpgdet = detJaco*WPG;
      iJaco.inv(Jaco);
      dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);

      // Modificado x Beto
      // double Area   = npg*wpgdet;
      double Area = detJaco*WPG_SUM;
      // fin modificado x Beto

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

	Jaco_axi.setel(Jaco.get(ind_axi_1,ind_axi_1),1,1);
	Jaco_axi.setel(Jaco.get(ind_axi_1,ind_axi_2),1,2);
	Jaco_axi.setel(Jaco.get(ind_axi_2,ind_axi_1),2,1);
	Jaco_axi.setel(Jaco.get(ind_axi_2,ind_axi_2),2,2);

	detJaco_axi = Jaco_axi.det();
	double Area_axi = 0.5*detJaco_axi*WPG_SUM;
	h_pspg = sqrt(4.*Area_axi/pi);
	Delta = sqrt(Area_axi);
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
	u.prod(SHAPE,ucols,-1,-1,1);
	vf.prod(SHAPE,vfcols,-1,-1,1);

	p_star = double(tmp8.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);
	vf_star.prod(SHAPE,vfcols_star,-1,-1,1);

	grad_u.prod(dshapex,ucols,1,-1,-1,2);
	grad_vf.prod(dshapex,vfcols,1,-1,-1,2);

	grad_u_star.prod(dshapex,ucols_star,1,-1,-1,2);
	grad_vf_star.prod(dshapex,vfcols_star,1,-1,-1,2);
	grad_p_star.prod(dshapex,pcol_star,1,-1,-1);

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// compute mixture density rho_m

      	alpha_g_vp.set(vf);
	rho_m_old = compute_rho_m(rho_g_vp,arho_g_vp,alpha_g_vp,alpha_l,rho_l,nphases);

	alpha_g_vp.set(vf_star);
	rho_m = compute_rho_m(rho_g_vp,arho_g_vp,alpha_g_vp,alpha_l,rho_l,nphases);

	grad_rho_m.set(0.);	
	for (int j=1; j<=nphases; j++) {
 	  rho_g = rho_g_vp.get(j);
	  grad_vf_star.ir(2,j);
	  tmp5_g.set(grad_vf_star).scale(rho_g-rho_l);
	  grad_rho_m.add(tmp5_g);
	}
	grad_vf_star.rs();

	v_m.set(u);   // probamos con poner la velocidad del paso anterior
	compute_vel_g(v_m,vslip_vp,vslip_user_vp,vslip_m_vp,v_g_vp_old,rho_g_vp,
		      rho_m_old,g_dir,d_bubble_vp,nphases,use_modified_slag_vslip,
		      vf);

	v_m.set(u_star);
	compute_vel_g(v_m,vslip_vp,vslip_user_vp,vslip_m_vp,v_g_vp,rho_g_vp,
		      rho_m,g_dir,d_bubble_vp,nphases,use_modified_slag_vslip,
		      vf_star);

	strain_rate.set(grad_u_star);
	grad_u_star.t();
	strain_rate.add(grad_u_star).scale(0.5);
	grad_u_star.rs();

	// Smagorinsky turbulence model
	if (LES) {
	  double tr = (double) tmp15.prod(strain_rate,strain_rate,-1,-2,-1,-2);
	  double van_D;
	  if (0 && A_van_Driest>0.) {
	    dist_to_wall.prod(SHAPE,xloc,-1,-1,1).rest(wall_coords);
	    double ywall = sqrt(dist_to_wall.sum_square_all());
	    double y_plus = ywall*shear_vel/VISC;
	    van_D = 1.-exp(-y_plus/A_van_Driest);
	    if (k % 250==0) printf("van_D: %f\n",van_D);
	  } else van_D = 1.;
	  
	  visco_t = SQ(C_smag*Delta*van_D)*sqrt(2*tr);
	} else {
	  visco_t = 0.0;
	}
	
	double visco_sato = 0.;
	if (nphases==1 && Sato_model_coef>0) {
	  for (int j=1; j<=nphases; j++) {
	    d_bubble = d_bubble_vp.get(j);
	    vslip = vslip_m_vp.get(j);
	    alpha_g = alpha_g_vp.get(j);
	    visco_sato = Sato_model_coef*rho_l*alpha_g*d_bubble*vslip;
	  }
	}
	
	visco_m_eff = alpha_l*visco_l+visco_sato + rho_m*visco_t;
	for (int j=1; j<=nphases; j++) {
	  alpha_g = alpha_g_vp.get(j);
	  visco_g = visco_g_vp.get(j);
	  visco_m_eff += alpha_g*visco_g;
	}
	visco_m_eff = (visco_m_eff <= 0 ? 1.0e-15 : visco_m_eff);
	
	for (int j=1; j<=nphases; j++) {
 	  //	  rho_g = rho_g_vp.get(j);
	  //	  visco_g_eff = rho_g*visco_t/Sc_number;
	  visco_g_eff = visco_t/Sc_number;
	  visco_g_eff = (visco_g_eff <= 0 ? 1.0e-15 : visco_g_eff);
	  visco_g_eff_vp.setel(visco_g_eff,j);
	}
	
	vel_axi(u,u_axi,axi);
	u2 = u_axi.sum_square_all();
	velmod = sqrt(u2);
	
	h_supg = compute_h_supg(u_axi,dshapex,velmod,h_pspg);

	Peclet = velmod * h_supg / (2. * visco_m_eff);
	
	tau_supg = tsf*SQ(2.*rec_Dt)+SQ(2.*velmod/h_supg)
	  +9.*SQ(4.*visco_m_eff/SQ(h_supg));
	tau_supg = 1./sqrt(tau_supg);

	tau_pspg = tsf*SQ(2.*rec_Dt)+SQ(2.*velmod/h_pspg)
	  +9.*SQ(4.*visco_m_eff/SQ(h_pspg));
	tau_pspg = 1./sqrt(tau_pspg);

	fz = (Peclet < 3. ? Peclet/3. : 1.);
	delta_supg = 0.5*h_supg*velmod*fz;
	
	if (tau_fac != 1.) {
	  tau_pspg *= tau_fac;
	  tau_supg *= tau_fac;
	}
	tau_pspg *= tau_pspg_fac;

	tau_pspg += additional_tau_pspg;
	
	delta_supg *= shock_capturing_factor;
	
	// P_supg es un vector fila
	P_supg.prod(u,dshapex,-1,-1,1).scale(tau_supg);

	// Weight function 
	W_supg.set(P_supg).add(SHAPE);

	// Pressure stabilizing term
	// use rho_l as a characteristic density to avoid numerical problems
	// for pressure stabilization
	P_pspg.set(dshapex).scale(tau_pspg/rho_l);  //debug:=

	velmod_g_vp.set(0.);
	for (int j=1; j<=nphases; j++) {
	  v_g_vp_old.ir(2,j);
	  u_axi_g_vp.ir(2,j);
	  vel_axi(v_g_vp_old,u_axi_g_vp,axi);
	  u2_g = u_axi_g_vp.sum_square_all();
	  velmod_g_vp.setel(sqrt(u2_g),j);
	}
	u_axi_g_vp.rs();
	v_g_vp_old.rs();
	
	// upwind for disperse phases
	for (int j=1; j<=nphases; j++) {
	  u_axi_g_vp.ir(2,j);
	  velmod_g = velmod_g_vp.get(j);
	  h_supg_g = compute_h_supg(u_axi_g_vp,dshapex,velmod_g,h_pspg);
	  h_supg_vp.setel(h_supg_g,j);
	  u_axi_g_vp.rs();
	}
	
	// SUPG perturbation function for disperse phases
	// velocidad del gas	
	P_supg_vp.set(0.);
	for (int j=1; j<=nphases; j++) {
	  v_g_vp_old.ir(2,j);
	  P_supg_vp.ir(2,j);
	  h_supg_g = h_supg_vp.get(j);
	  visco_g = visco_g_eff_vp.get(j);
	  velmod_g = velmod_g_vp.get(j);
	  tau_supg_g = SQ(2.*velmod_g/h_supg_g)+9.*SQ(4.*visco_g/SQ(h_supg_g));
	  tau_supg_g = 1./sqrt(tau_supg_g);	  

	  tau_supg_g *= tau_fac;

	  P_supg_vp.prod(v_g_vp_old,dshapex,-1,-1,1).scale(tau_supg_g);

	  Peclet = velmod_g * h_supg_g / (2. * visco_g);

	  fz = (Peclet < 3. ? Peclet/3. : 1.);
	  delta_supg_g = 0.5*h_supg_g*velmod_g*fz;

	  delta_supg_g *= shock_capturing_factor;

	  delta_supg_vp.setel(delta_supg_g,j);

	}
	v_g_vp_old.rs();
	P_supg_vp.rs();

	// implicit version - General Trapezoidal rule - parameter alpha
#ifdef ADD_GRAD_DIV_U_TERM
	dmatu.prod(u_star,grad_u_star,-1,-1,1);
#else
	dmatu.prod(u,grad_u_star,-1,-1,1);
#endif

	// momentum equation residue	
	du.set(u_star).rest(u);
	dmatu.axpy(du,rec_Dt/alpha).rest(G_body);
	
	resmom_prime.set(grad_p_star).axpy(dmatu,rho_m);

	div_u_star = double(tmp10.prod(dshapex,ucols_star,-1,-2,-2,-1));

	// Galerkin - momentum
	// resmom tiene que tener nel*ndim
	dresmom.t().prod(dmatu,SHAPE,1,2).rs();
	resmom.axpy(dresmom,-wpgdet * rho_m);

	if (weak_form) {
	  tmp1.set(strain_rate).scale(2*visco_m_eff).axpy(eye,-p_star);
	  tmp2.prod(dshapex,tmp1,-1,1,-1,2);
	  resmom.axpy(tmp2,-wpgdet);
	} else {
	  tmp6.prod(dshapex,strain_rate,-1,1,-1,2).scale(2*visco_m_eff);
	  tmp11.prod(SHAPE,grad_p_star,1,2).add(tmp6);
	  resmom.axpy(tmp11,-wpgdet);
	}

	// SUPG perturbation - momentum
	resmom_supg.prod(P_supg,resmom_prime,1,2);
	resmom.axpy(resmom_supg,-wpgdet);

	// shock capturing term - momentum
	resmom.axpy(dshapex.t(),-wpgdet*delta_supg*rho_m*div_u_star);
	dshapex.rs();
	tmp6_g.prod(grad_rho_m,u_star,-1,-1);
	resmom.axpy(dshapex.t(),-wpgdet*delta_supg*double(tmp6_g));
	dshapex.rs();

	// Galerkin - continuity
	rescont_gal.prod(dshapex,u_star,-1,1,-1);
	rescont.axpy(rescont_gal,-wpgdet*rho_m);

	// PSPG perturbation - continuity
	rescont_pspg.prod(P_pspg,resmom_prime,-1,1,-1);
	rescont.axpy(rescont_pspg,wpgdet);
	
	// Penalization term?
	if (pressure_control_coef) {
	  double p_star = tmp21.prod(SHAPE,pcol_star,-1,-1).get();
	  rescont.axpy(SHAPE,pressure_control_coef*p_star*wpgdet);
	}

	// temporal part + convective 
#ifdef ADD_GRAD_DIV_U_TERM
	massm.prod(u_star,dshapex,-1,-1,1);
#else
	massm.prod(u,dshapex,-1,-1,1);
#endif
	massm.axpy(SHAPE,rec_Dt/alpha);
	matlocmom.prod(W_supg,massm,1,2).scale(rho_m);

	gNi_gNj.prod(dshapex,dshapex,-1,1,-1,2);
	Ni_Nj.prod(SHAPE,SHAPE,1,2);
	gNi_Nj.prod(dshapex,SHAPE,1,2,3);

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// disperse phases part
	du_g.set(vf_star).rest(vf).scale(rec_Dt/alpha);

	tmp11_g.prod(dshapex,grad_vf_star,-1,1,-1,2);

	for (int j=1; j<=nphases; j++) {
	  res_alpha_g.ir(2,j);
	  P_supg_vp.ir(2,j);
	  v_g_vp.ir(2,j);
	  grad_vf_star.ir(2,j);
	  tmp11_g.ir(2,j);
	  dresmom_g.ir(2,j);

	  // SUPG term
	  double dmatu_g = double(tmp12_g.prod(v_g_vp,grad_vf_star,-1,-1));
	  dmatu_g += double(du_g.get(j));
	  dresmom_g.set(P_supg_vp).scale(dmatu_g);
	  // Galerkin terms
	  dresmom_g.axpy(SHAPE,double(du_g.get(j)));
	  tmp1_g.prod(dshapex,v_g_vp,-1,1,-1).scale(double(alpha_g_vp.get(j)));
	  dresmom_g.rest(tmp1_g);
	  // Diffusion term
	  dresmom_g.axpy(tmp11_g,double(visco_g_eff_vp.get(j)));
	  // Shock capturing term
	  dresmom_g.axpy(tmp11_g,double(delta_supg_vp.get(j)));
	  // Addition
	  res_alpha_g.axpy(dresmom_g,-wpgdet);

	  matloc_alpha_g.ir(2,j).ir(4,j);
	  // Ni Nj / Dt
	  matloc_alpha_g.set(Ni_Nj).scale(rec_Dt/alpha);
	  // P_supg_i Nj / Dt
	  tmp3_g.prod(P_supg_vp,SHAPE,1,2).scale(rec_Dt/alpha);
	  matloc_alpha_g.add(tmp3_g);
	  // tmp2_g = gNi v_g
	  tmp2_g.prod(dshapex,v_g_vp,-1,1,-1);
	  // tmp3_g = P_supg_i v_g gNj
	  tmp3_g.prod(P_supg_vp,tmp2_g,1,2);
	  // Addition
	  matloc_alpha_g.add(tmp3_g);
	  // - gNi v_g Nj ( weak form term)
	  tmp4_g.prod(tmp2_g,SHAPE,1,2);
	  matloc_alpha_g.rest(tmp4_g);
	  // diffusion term
	  matloc_alpha_g.axpy(gNi_gNj,double(visco_g_eff_vp.get(j)));
	  // shock-capturing term
	  matloc_alpha_g.axpy(gNi_gNj,double(delta_supg_vp.get(j)));

	  if(use_modified_slag_vslip && j==nphases){
	    vslip = double(vslip_user_vp.get(j));
	    rho_g = rho_g_vp.get(j);
	    // Galerkin term
	    // tmp2_g = gNi alpha_g * dVslip_s/dalpha_s Nj
	    dshapex.ir(1,abs(g_dir));
	    tmp3_g.prod(dshapex,SHAPE,1,2).scale(-double(alpha_g_vp.get(j)));
	    double dvslip_dalpha_s = vslip*(rho_g-rho_l)*rho_g/rho_m/rho_m;
	    matloc_alpha_g.axpy(tmp3_g,dvslip_dalpha_s);
	    dshapex.rs();
	    // SUPG term
	    // tmp3_g = P_supg_i v_g gNj	    
	    tmp3_g.prod(P_supg_vp,SHAPE,1,2);
	    tmp3_g.scale(double(grad_vf_star.get(abs(g_dir))));
	    matloc_alpha_g.axpy(tmp3_g,dvslip_dalpha_s);
	  }
	  
	}
	tmp11_g.rs();
	res_alpha_g.rs();
	P_supg_vp.rs();
	v_g_vp.rs();
	grad_vf_star.rs();
	matloc_alpha_g.rs();
	dresmom_g.rs();

	if (update_jacobian) {
	  matlocf.is(2,ndim+2,ndim+1+nphases).is(4,ndim+2,ndim+1+nphases);
	  matlocf.axpy(matloc_alpha_g,wpgdet).rs();
	}

 	for (int j=1; j<=nphases; j++) {
	    // Galerkin contribution --> dR_cont/d_alpha_g
	    matlocf.rs().ir(2,ndim+1).ir(4,ndim+1+j);
	    tmp10_g.prod(rescont_gal,SHAPE,1,2);
	    matlocf.axpy(tmp10_g,(rho_g-rho_l)*wpgdet);

	    // termino agregado por la variacion de rho_m con alpha_g
	    // en las ecuaciones de momento
	    tmp171.prod(W_supg,dmatu,1,2);
	    tmp172.prod(tmp171,SHAPE,1,2,3);
	    matlocf.rs().is(2,1,ndim).ir(4,ndim+1+j);
	    matlocf.axpy(tmp172,(rho_g-rho_l)*wpgdet).rs();
	    
	    // termino agregado por la variacion de rho_m con alpha_g
	    // en la ecuacion de continuidad
	    tmp173.prod(P_pspg,dmatu,-1,1,-1);
	    tmp174.prod(tmp173,SHAPE,1,2);
	    matlocf.rs().ir(2,ndim+1).ir(4,ndim+1+j);
	    matlocf.axpy(tmp174,-(rho_g-rho_l)*wpgdet).rs();	    
	}

	matlocf.rs();
	double vaux;
	// extra term in continuity and momentum equations
	resmom.ir(2,abs(g_dir));
 	for (int j=1; j<=nphases; j++) {

	  P_pspg.ir(1,abs(g_dir));

	  rho_g = rho_g_vp.get(j);
	  alpha_g = alpha_g_vp.get(j);
	  if (1) {
	  // version incluyendo la modificacion de la vslip x escoria
	  vslip_m = vslip_m_vp.get(j);
	  vaux = rho_g*vslip_m*vslip_m;
	  } else {
	  // version sin incluir esta modificacion 
	  vslip = vslip_vp.get(j);
	  vaux = rho_g*vslip*vslip;
	  }
	  dshapex.ir(1,abs(g_dir));
	  // Galerkin contribution (weak form)
	  resmom.axpy(dshapex,wpgdet*vaux*alpha_g);
	  // SUPG contribution
	  grad_vf_star.ir(1,abs(g_dir)).ir(2,j);
	  resmom.axpy(P_supg,-wpgdet*vaux*double(grad_vf_star));
	  // PSPG contribution
	  rescont.axpy(P_pspg,wpgdet*vaux*double(grad_vf_star));

	  // agregado termino por variacion de la vslip (slag)
	  if(use_modified_slag_vslip && (j==nphases) && 0 ) {
	    /*
	    double g_rho_m = double(grad_rho_m.get(abs(g_dir)));
	    double grad_vslip = double(vslip_user_vp.get(j))*rho_g/rho_m/rho_m*g_rho_m;
	    */
	    grad_vf_star.ir(1,abs(g_dir)).ir(2,nphases);
	    double grad_alpha = double(grad_vf_star);
	    grad_vf_star.rs();
	    double grad_vslip = -double(vslip_user_vp.get(j))*grad_alpha;
	    double grad_TE_alpha_2 = rho_g*alpha_g*vslip_m*grad_vslip;	    
	    // SUPG contribution
	    resmom.axpy(P_supg,-wpgdet*grad_TE_alpha_2);
	    // PSPG contribution
	    rescont.axpy(P_pspg,wpgdet*grad_TE_alpha_2);
	  }

	  if (update_jacobian && coupled) {

	    matlocf.rs().ir(2,abs(g_dir)).ir(4,ndim+1+j);
	    // Galerkin contribution = -gNi Nj
	    gNi_Nj.ir(1,abs(g_dir));
	    matlocf.axpy(gNi_Nj,-wpgdet*vaux);
	    // P_supg_i_gNj 
	    tmp9_g.prod(P_supg,dshapex,1,2);
	    matlocf.axpy(tmp9_g,wpgdet*vaux);

	    // PSPG contribution --> dR_cont(PSPG)/d_alpha_g
	    dshapex.ir(1,abs(g_dir));
	    tmp10_g.prod(P_pspg,dshapex,1,2);
	    matlocf.rs().ir(2,ndim+1).ir(4,ndim+1+j);
	    matlocf.axpy(tmp10_g,-wpgdet*vaux);

	  }
	}

	resmom.rs();
	dshapex.rs();
	P_pspg.rs();
	grad_vf_star.rs();
	matlocf.rs();
	gNi_Nj.rs();

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	
	// diffusive part
	matlocmom.axpy(gNi_gNj,visco_m_eff);

	// dmatw =  rho * ((1/Dt)*SHAPE + u * dshapex);
	dmatw.set(massm).scale(rho_m);

	if (update_jacobian) {
	  for (int iloc=1; iloc<=nel; iloc++) {
	    for (int jloc=1; jloc<=nel; jloc++) {
	      double c = wpgdet*matlocmom.get(iloc,jloc);
	      for (int ii=1; ii<=ndim; ii++) {
		// matloc.ir(1,iloc).ir(2,ii).ir(3,jloc).ir(4,ii);
		matlocf.addel(c,iloc,ii,jloc,ii);
	      }
	    }
	  }

	  if (weak_form) {
	    tmp16.prod(P_supg,dshapex,1,2,3).scale(wpgdet);
	    tmp162.prod(dshapex,SHAPE,2,1,3).scale(-wpgdet);
	    matlocf.is(2,1,ndim).ir(4,ndim+1).add(tmp16)
	      .add(tmp162).rs();
	  } else {
	    tmp16.prod(W_supg,dshapex,1,2,3).scale(wpgdet);
	    matlocf.is(2,1,ndim).ir(4,ndim+1).add(tmp16).rs();
	  }
	  
	  
	  tmp13.prod(P_pspg,dshapex,-1,1,-1,2);
	  if (pressure_control_coef) {
	    tmp20.prod(SHAPE,SHAPE,1,2);
	    tmp13.axpy(tmp20,pressure_control_coef);
	  }	  

	  // d_R(cont)/d_vel
	  matlocf.ir(2,ndim+1).is(4,1,ndim);
	  tmp17.prod(P_pspg,dmatw,3,1,2).scale(wpgdet);
	  matlocf.rest(tmp17);
	  tmp17.prod(dshapex,SHAPE,3,1,2).scale(-rho_m*wpgdet);
	  matlocf.rest(tmp17).rs();
	  // d_R(cont)/d_p
	  matlocf.ir(2,ndim+1).ir(4,ndim+1).axpy(tmp13,-wpgdet).rs();

	  matlocf.rs();

	  if (!cache_grad_div_u) {
	    tmp19.set(dshapex).scale(visco_m_eff*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,3,4,1);
	    matlocf.is(2,1,ndim).is(4,1,ndim).add(tmp18).rs();

	    tmp19.set(dshapex).scale(delta_supg*rho_m*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,1,4,3);
	    matlocf.is(2,1,ndim).is(4,1,ndim).add(tmp18).rs();

	    tmp19.set(dshapex).scale(delta_supg*wpgdet);
	    tmp7_g.prod(grad_rho_m,SHAPE,2,1);
	    tmp18.prod(tmp7_g,tmp19,1,2,4,3);
	    matlocf.is(2,1,ndim).is(4,1,ndim).add(tmp18).rs();

	  } else {
	    grad_div_u_coef += (delta_supg*rho_m+visco_m_eff)*wpgdet;
	    if (!grad_div_u_was_cached) {
	      tmp18.prod(dshapex,dshapex,2,1,4,3);
	      grad_div_u.add(tmp18);
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
#ifdef	DEBUG_CACHE_GRAD_DIV_U	// debug:=
	  if (k<2) {
	    printf("element %d, computed grad_div_u: ",k);
	    for (int kkkk=0; kkkk<ndim*ndim*nel*nel; kkkk++)
	      printf("%f	",grad_div_u_cache[kkkk]);
	    printf("\n");
	  }
#endif
	}
      }

      //      veccontr.is(2,1,ndim).set(resmom)
      //.rs().ir(2,ndim+1).set(rescont).rs();

      veccontr.is(2,1,ndim).set(resmom)
	.rs().ir(2,ndim+1).set(rescont)
	.rs().is(2,ndim+2,ndim+1+nphases).set(res_alpha_g).rs();


      if (residual_factor!=1.)
	veccontr.scale(residual_factor);
      veccontr.export_vals(&(RETVAL(ielh,0,0)));

      if (update_jacobian) {
	if (jacobian_factor!=1.)
	  matlocf.scale(jacobian_factor);
	matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      }
#ifdef PRINT_ELEM_DEBUG
      if (k==0) {
	veccontr.print("veccontr:");
	matlocf.print("matlocf:");
      }
#endif
      // Fixme : later para debugear por que se va a la mierda
      //	printf("element %d, residuo : ",k);
      //	veccontr.print("veccontr:");
	// fin Fixme

    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
}

#undef SHAPE	
#undef DSHAPEXI 
#undef WPG	
#undef SQ
