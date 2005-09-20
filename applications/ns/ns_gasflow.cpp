//__INSERT_LICENSE__
//$Id: ns_gasflow.cpp,v 1.2 2005/09/20 01:56:43 mstorti Exp $

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

// auxiliary functions

// set state
#define set_state {						\
	vel_supg.set(vel).rest(v_mesh).rs();			\
        if(axi>0){						\
          vel_supg.setel(0.,axi);				\
        }							\
	u2 = vel_supg.sum_square_all();				\
	velmod = sqrt(u2);					\
	g1=ga-1;						\
	rho_ene = p/g1+0.5*rho*square(velmod);	\
	ene = rho_ene/rho;					\
	entalpy=rho_ene+p;					\
        Cv=Rgas/g1;						\
        int_ene=p/g1/rho;				\
        Cp = ga*Cv;						\
	grad_T.set(grad_p).axpy(grad_rho,-p/rho).scale(1./Rgas/rho).rs(); \
}  

// ======================================================================================
#define compute_Cp_matrix {							\
  Cpjac.set(0.);								\
  Cpjac.setel(1.0,1,1);							\
  Cpjac.is(1,vl_indx,vl_indxe).ir(2,1).set(vel).rs();			\
  Cpjac.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).eye(rho).rs();	\
  Cpjac.ir(1,ndof).is(2,vl_indx,vl_indxe).set(vel).scale(rho).rs();	\
  Cpjac.setel(1.0/g1,ndof,ndof);						\
  Cpjac.setel(0.5*square(velmod),ndof,1);					\
}

// ======================================================================================
// Advective fluxes							
#define compute_Fa {					\
  flux.set(0.);						\
  flux.ir(1,1).set(vel).scale(rho).rs();		\
  Amom.prod(vel,vel,1,2);				\
  Y.set(0.).add(Amom).scale(rho).axpy(Id_ndim,p);	\
  flux.is(1,vl_indx,vl_indxe).add(Y).rs();		\
  flux.ir(1,ndof).set(vel).scale(entalpy).rs();		\
}

// ======================================================================================
// Diffusive fluxes									
#define compute_Fd {									\
  fluxd.set(0.);									\
  fluxd.is(1,vl_indx,vl_indxe).set(sigma).scale(visco_eff).rs();			\
  viscous_work.prod(sigma,vel,1,-1,-1).scale(visco_l);					\
  heat_flux.set(grad_T).scale(-cond_eff);						\
  fluxd.ir(1,ndof).set(viscous_work).rest(heat_flux).rs();				\
}

// ======================================================================================
// Shock capturing fluxes								
#define compute_shocap {							\
  flux_sc.set(0.);		                                                \
  flux_sc.prod(Cpjac_old,grad_U,1,-1,2,-1);                                     \
  for (int j=1; j<=ndim; j++) {							\
    flux_sc.ir(2,j).mult(delta_sc_v);    					\
  }										\
  flux_sc.rs();                  						\
}

// ======================================================================================
// Smagorinsky turbulence model						
#define compute_visco_LES {								\
	if (LES) {									\
	  double tr = (double) tmp15.prod(strain_rate,strain_rate,-1,-2,-1,-2);		\
	  double van_D;									\
	  if (A_van_Driest>0.) {							\
	    dist_to_wall.prod(SHAPE,xloc,-1,-1,1).rest(wall_coords);			\
	    double ywall = sqrt(dist_to_wall.sum_square_all());				\
	    double y_plus = ywall*shear_vel/VISC;					\
	    van_D = 1.-exp(-y_plus/A_van_Driest);					\
	    if (k % 250==0) printf("van_D: %f\n",van_D);				\
	  } else van_D = 1.;								\
	  nu_tur = SQ(C_smag*Delta*van_D)*sqrt(2*tr);					\
 	  nu_eff = VISC + nu_tur;							\
	} else {									\
	  nu_eff = VISC;								\
	}										\
        cond_t = Cp*rho*nu_tur/Pr_t;							\
        sutherland_factor = 1.0;        						\
        Tem = p/Rgas/rho;								\
											\
        if (sutherland_law !=0 )							\
           sutherland_factor = pow(Tem/Tem_infty,1.5)*                                  \
                                    ((Tem_infty+Tem_ref)/(Tem_ref+Tem));	        \
											\
        visco_l = sutherland_factor * VISC;						\
        visco_eff = (visco_l + rho*nu_tur);						\
        cond_eff = sutherland_factor * cond + cond_t;					\
        visco_bar = visco_l - cond_eff/Cv;						\
											\
	grad_u.d(1,2);									\
	tmp10.sum(grad_u,-1);								\
	grad_u.rs();									\
	double divu = double(tmp10);							\
											\
	sigma.set(0.);									\
	sigma.eye(divu).scale(-2./3.).axpy(strain_rate,2.0).rs();			\
}

// ======================================================================================
// Adjective Jacobians in primitive basis						      
#define compute_Ajac {								\
  Ajac.set(0.);									\
  Ajac.ir(2,1).ir(3,1).set(vel).rs();						\
  Ajac.ir(2,1).is(3,vl_indx,vl_indxe).set(Id_ndim).scale(rho).rs();		\
  Ajac.ir(2,ndof).ir(3,1).set(vel).scale(0.5*square(velmod)).rs();		\
  Ajac.ir(2,ndof).ir(3,ndof).set(vel).scale(ga/g1).rs();			\
  Ajac.ir(2,ndof).is(3,vl_indx,vl_indxe).prod(vel,vel,1,2).scale(rho).rs();     \
  Ajac.ir(2,ndof).is(3,vl_indx,vl_indxe).axpy(Id_ndim,ga/g1*p+0.5*rho*square(velmod)).rs();   \
  Ajac.is(2,vl_indx,vl_indxe).ir(3,1).prod(vel,vel,1,2).rs();				      \
  Ajac.is(2,vl_indx,vl_indxe).ir(3,ndof).set(Id_ndim).rs();				      \
  Ajac.is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe).prod(vel,Id_ndim,1,2,3).scale(rho).rs(); \
  Y3.prod(vel,Id_ndim,2,1,3).scale(rho).rs();				\
  Ajac.is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe).add(Y3).rs();	\
}

// ======================================================================================

// Diffusive jacobians in primitive variables
#define compute_Djac {								\
  if(ndim==3) {									\
    ip[0] = 1;									\
    ip[1] = 2;									\
    ip[2] = 3;									\
    ip[3] = 1;									\
  } else {									\
    ip[0] = 1;									\
    ip[1] = 2;									\
    ip[2] = 1;									\
  }										\
  Djac.set(0.);									\
  Djac.is(2,vl_indx,vl_indxe).is(4,vl_indx,vl_indxe).axpy(IdId,visco_eff).rs();	\
  tmp00 = 1./3.*visco_eff;							\
  for (int j=1; j<=ndim; j++) {							\
    Djac.addel(tmp00,j,j+1,j,j+1);						\
  }										\
  tmp01 = -cond_eff*Tem/rho;							\
  tmp02 = cond_eff/rho/Cv/g1;							\
  Djac.ir(2,ndof).ir(4,1).axpy(Id_ndim,tmp01).rs();				\
  Djac.ir(2,ndof).ir(4,ndof).axpy(Id_ndim,tmp02).rs();				\
  tmp_vel_ndim.prod(Id_ndim,vel,1,2,3);						\
  Djac.ir(2,ndof).is(4,vl_indx,vl_indxe).axpy(tmp_vel_ndim,visco_eff).rs();	\
  tmp_vel.set(vel).scale(tmp00);						\
  for (int j=1; j<=ndim; j++) {							\
    Djac.addel(tmp_vel.get(j),j,ndof,j,j+1);					\
  }										\
  tmp00 = 2./3.*visco_eff;							\
  for (int j=0; j<ndim; j++) {							\
    Djac.addel(-2./3.*visco_eff,ip[j],ip[j]+1,ip[j+1],ip[j+1]+1);		\
    Djac.addel(visco_eff,ip[j],ip[j+1]+1,ip[j+1],ip[j]+1);			\
    Djac.addel(-2./3.*visco_eff,ip[j+1],ip[j+1]+1,ip[j],ip[j]+1);		\
    Djac.addel(visco_eff,ip[j+1],ip[j]+1,ip[j],ip[j+1]+1);			\
    Djac.addel(visco_eff*vel.get(ip[j+1]),ip[j],ndof,ip[j+1],ip[j]+1);		\
    Djac.addel(-2./3.*visco_eff*vel.get(ip[j]),ip[j],ndof,ip[j+1],ip[j+1]+1);	\
    Djac.addel(-2./3.*visco_eff*vel.get(ip[j+1]),ip[j+1],ndof,ip[j],ip[j]+1);	\
    Djac.addel(visco_eff*vel.get(ip[j]),ip[j+1],ndof,ip[j],ip[j+1]+1);		\
  }										\
}

// ======================================================================================
// Stabilization parameters
#define compute_tau_supg {						\
    double tol=1.0e-16;							\
    h_supg = compute_h_supg(vel_supg,dshapex,velmod,h_pspg);		\
    double Peclet;							\
    double sonic_speed = sqrt(ga*p/rho);				\
    double velmax = velmod+sonic_speed;					\
    tau_supg_a = h_supg/2./velmax;					\
    Peclet = velmod * h_supg / (2. * visco_supg);			\
    double tol_shoc = 1e-010;						\
    double h_shoc, grad_rho_mod = sqrt(grad_rho.sum_square_all());	\
    FastMat2::branch();							\
    if(grad_rho_mod>tol_shoc) {						\
      FastMat2::choose(0);						\
      jvec.set(grad_rho).scale(1.0/grad_rho_mod);			\
      h_shoc = tmp9.prod(dshapex,jvec,-1,1,-1).sum_abs_all();		\
      h_shoc = (h_shoc < tol ? tol : h_shoc);				\
      h_shoc = 2./h_shoc;						\
    } else {								\
      FastMat2::choose(1);						\
      jvec.set(0.);							\
      h_shoc = h_supg;							\
    }									\
    FastMat2::leave();							\
    double fz = grad_rho_mod*h_shoc/rho;				\
    fz = pow(fz,shocap_beta);						\
    delta_sc_aniso = 0.5*h_shoc*velmax*fz;				\
    double fz2 = (Peclet < 3. ? Peclet/3. : 1.);			\
    delta_sc = 0.5*h_supg*velmax*fz2;					\
    Cpi.set(0.);                                                                  \
    Cpi.setel(1.0,1,1);                                                           \
    Cpi.is(1,vl_indx,vl_indxe).ir(2,1).set(vel).scale(-1.0/rho).rs();             \
    Cpi.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).eye(1.0/rho).rs();             \
    Cpi.ir(1,ndof).is(2,vl_indx,vl_indxe).set(vel).scale(-g1).rs();             \
    Cpi.setel(g1,ndof,ndof);             \
    Cpi.setel(0.5*square(velmod)*g1,ndof,1);             \
}

// ======================================================================================

#undef __FUNC__
#define __FUNC__ "ns_gasflow::assemble"
int ns_gasflow::
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

  // Hloc stores the old mesh coordinates
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),vloc_mesh(nel,ndim),v_mesh(ndim);

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
  //o Add shock-capturing term.
  SGETOPTDEF(double,shock_capturing_factor,0);
  //o Add pressure controlling term. 
  SGETOPTDEF(double,pressure_control_coef,0.);
  assert(pressure_control_coef>=0.);

  //o ALE_flag : flag to ON ALE computation
  SGETOPTDEF(int,ALE_flag,0);
  //o indx_ALE_xold : pointer to old coordinates in
  //  NODEDATA array excluding the first "ndim" values
  SGETOPTDEF(int,indx_ALE_xold,1);

  // allocate local vecs
  int kdof;
  FastMat2 veccontr(2,nel,ndof),xloc(2,nel,ndim),locstate(2,nel,ndof), 
    locstate2(2,nel,ndof),xpg,G_body(1,ndim),vrel,Id_ndim(2,ndim,ndim);

  if (ndof != ndim+2) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FastMat2 res(2,nel,ndof); //,reslocf(2,nel,ndof);
  FastMat2 mat(4,nel,ndof,nel,ndof),matlocf(4,nel,ndof,nel,ndof);

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

  // gamma coeficient
  SGETOPTDEF(double,ga,1.4);
  // constant of a particular gas for ideal gas law (state equation for the gas)
  SGETOPTDEF(double,Rgas,287.);

  // tau scheme 
  // [0] standard [default] tau = max(0, tau_a-tau_d-tau_delta) 
  // [1] see Tezduyar & Senga WCCM VI (2004) 
  //     inv(tau[jdof,jdof])^r_switch = 
  //                    inv(tau_adv[jdof])^r_switch+inv(tau_dif[jdof])^r_switch 
  // 
  SGETOPTDEF(int,tau_scheme,0);

  // r-switch : exponent for stabilization switching 
  //            see Tezduyar & Senga WCCM VI (2004) 
  SGETOPTDEF(double,r_switch,2.0);

  // beta parameter for shock capturing, see Tezduyar & Senga WCCM VI (2004)
  SGETOPTDEF(double,shocap_beta,1.0);
  //o Add a shock capturing term 
  // (copy from advdife.cpp)
  SGETOPTDEF(double,shocap,0.0);

  // Sutherland law key
  SGETOPTDEF(int,sutherland_law,0);
  // Sutherland law infinity Temperature
  SGETOPTDEF(double,Tem_infty,0.0);
  // Sutherland law reference Temperature
  SGETOPTDEF(double,Tem_ref,0.0);
  // Sutherland law implicitly
  SGETOPTDEF(int,sutherland_law_implicit,0);
  //o Turbulent Prandtl number
  SGETOPTDEF(double,Pr_t,1.);
  //o molecular dynamic viscosity 
  SGETOPTDEF(double,visco,0);
  //o molecular conductivity
  SGETOPTDEF(double,cond,0);

  //o _T: double[ndim] _N: G_body _D: null vector 
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"G_body",G_body.storage_begin(),1,ndim);

  double pi = 4*atan(1.0);

  /*
  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))
  */

#define VISC (visco)

  int nprops=iprop;
  
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

  FMatrix Jaco(ndim,ndim),iJaco(ndim,ndim),grad_u_old,
    grad_u(ndim,ndim),grad_u_star,strain_rate(ndim,ndim);

  FMatrix grad_p_star(ndim),u_star,du,uintri(ndim),dmatu(ndim),
    ucols,ucols_new,ucols_star,rhocol,rhocol_new,rhocol_star,grad_rho_star(ndim),
    pcol_star,pcol_new,pcol,fm_p_star,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,
    massm,tmp7,tmp8,tmp9,tmp10,tmp11,tmp13,tmp14,tmp15,dshapex_c,xc,
    wall_coords(ndim),dist_to_wall,tmp16,tmp162,tmp17,tmp18,tmp19;
  FastMat2 tmp20(2,nel,nel),tmp21,vel_supg;
  FastMat2 Up_old(1,ndof),Uc_old(1,ndof),Up(1,ndof),Uc(1,ndof);
  FastMat2 Ajac,Djac,Ao_grad_N,A_grad_N,sigma(2,ndim,ndim),grad_rho(1,ndim),grad_rho_old(1,ndim),
    jvec(1,ndim);
  FastMat2 vel,tau_supg_m(2,ndof,ndof),A_grad_U,grad_U(2,ndim,ndof),G_source(1,ndof),
    flux(2,ndof,ndim),fluxd(2,ndof,ndim),flux_sc(2,ndof,ndim),tau_supg_c(2,ndof,ndof),
    delta_sc_v(1,ndof);
  FastMat2 tmp1_G,tmp1_S,tmp2_G,tmp2_S,tmp3_G,tmp4_G,tmp4_S,tmp5_S;
  FastMat2 tmp1_G_M,tmp1_S_M,tmp2_G_M,tmp2_S_M,tmp3_G_M,tmp5_S_M;
  FastMat2 Ni_Nj,Pi_Nj,gNi_Nj,gNi_gNj,gNi_dot_gNj;
  FastMat2 Cpjac,Cpjac_old,Cpjac_old_delta,Cpi,dUcdt,dUpdt,grad_T(1,ndim),
    heat_flux(1,ndim),viscous_work(1,ndim);
  FastMat2 IdId,tmp_vel_ndim,tmp_vel;

  double tmp00,tmp01,tmp02;

  FastMat2 grad_p_old(1,ndim),grad_p(1,ndim);	
  FastMat2 Amom,Y(2,ndim,ndim),Y3;

  Id_ndim.eye();
  IdId.resize(4,ndim,ndim,ndim,ndim);					
  IdId.prod(Id_ndim,Id_ndim,1,3,2,4);					
  Ajac.resize(3,ndim,ndof,ndof);
  Cpjac.resize(2,ndof,ndof);
  Cpjac_old.resize(2,ndof,ndof);
  Cpjac_old_delta.resize(2,ndof,ndof);
  Cpi.resize(2,ndof,ndof);
  Djac.resize(4,ndim,ndof,ndim,ndof);

 vector<int> ip;
  ip.resize(ndim+1);

  double tmp12;
  double tsf = temporal_stability_factor;

  int vl_indx=2, vl_indxe = vl_indx+ndim-1;

  double rho,p,rho_old,p_old,rho_star;
  double visco_supg,cond_supg,nu_eff,cond_eff,tau_supg_a,visco_eff,nu_tur,
    visco_l,visco_bar,cond_t;
  double sutherland_factor = 1.0,Tem,delta_sc,delta_sc_aniso;

  FastMat2 u_old;

  FMatrix eye(ndim,ndim),seed,one_nel,matloc_prof(nen,nen);

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
      if(nH>0) Hloc.ir(1,kloc+1).set(&NODEDATA(node-1,0)+ndim);
    }
    xloc.rs();
    Hloc.rs();

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
    //                   locstate  <- u^*
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
#ifdef 	DEBUG_CACHE_GRAD_DIV_U	// debug:=
	if (k<2 && grad_div_u_was_cached) {
	  printf("element %d, cached grad_div_u: ",k);
	  for (int kkkk=0; kkkk<ndim*ndim*nel*nel; kkkk++)
	    printf("%f  ",grad_div_u_cache[kkkk]);
	  printf("\n");
	}
	grad_div_u_was_cached = 0; // In order to recompute always
				   // the grad_div_u operator
#endif
      }
    }

    mat.set(0.);
    matlocf.set(0.);
    veccontr.set(0.);
    res.set(0.);
    //    reslocf.set(0.);

    if (comp_mat_res && cache_grad_div_u) {
      if (grad_div_u_was_cached) {
	grad_div_u.set(grad_div_u_cache);
      } else {
	grad_div_u.set(0.);
      }
    }

    if (comp_res || comp_mat_res) {
      ucols.set(locstate2.is(2,vl_indx,vl_indxe));
      pcol.set(locstate2.rs().ir(2,ndof));
      rhocol.set(locstate2.rs().ir(2,1));
      locstate2.rs();
      
      ucols_new.set(locstate.is(2,vl_indx,vl_indxe));
      pcol_new.set(locstate.rs().ir(2,ndof));
      rhocol_new.set(locstate.rs().ir(2,1));
      locstate.rs();
      
      ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
      pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
      rhocol_star.set(rhocol_new).scale(alpha).axpy(rhocol,1-alpha);

      //#define PRINT_ELEM_DEBUG
#ifdef PRINT_ELEM_DEBUG
      if (k==0) {
	locstate2.print("locstate2 (t_n):");
	locstate.print("locstate (t_n+1):");
      }
#endif
    }

    // nodal computation of mesh velocity
    if (ALE_flag) {
      assert(nH >= ndim);
      assert(indx_ALE_xold >= nH+1-ndim);
      Hloc.is(2,indx_ALE_xold,indx_ALE_xold+ndim-1);
      vloc_mesh.set(xloc).rest(Hloc).scale(rec_Dt).rs();
      Hloc.rs();
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
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])
#define WPG_SUM  (gp_data.wpg_sum)

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

      double Area = detJaco*WPG_SUM;
      double h_pspg,Delta;

      double g1,rho_ene,ene,entalpy,Cv,int_ene,Cp;

      if (ndim==2) {
	h_pspg = sqrt(4.*Area/pi);
	Delta = sqrt(Area);
      } else if (ndim==3 && axi==0) {
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
	// old state
	rho_old = double(tmp8.prod(SHAPE,rhocol,-1,-1));
	p_old   = double(tmp8.prod(SHAPE,pcol,-1,-1));
	u_old.prod(SHAPE,ucols,-1,-1,1);

	// current state
	rho_star = double(tmp8.prod(SHAPE,rhocol_star,-1,-1));
	p_star = double(tmp8.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);

	grad_u_old.prod(dshapex,ucols,1,-1,-1,2);
	grad_p_old.prod(dshapex,pcol,1,-1,-1);
	grad_rho_old.prod(dshapex,rhocol,1,-1,-1);

	grad_u_star.prod(dshapex,ucols_star,1,-1,-1,2);
	grad_p_star.prod(dshapex,pcol_star,1,-1,-1);
	grad_rho_star.prod(dshapex,rhocol_star,1,-1,-1);
  
	v_mesh.set(0.0);
	if (ALE_flag) {
	  v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	}

	// set state to old state
	rho=rho_old; p=p_old;vel.set(u_old);
	grad_rho.set(grad_rho_old);grad_p.set(grad_p_old);
	set_state;

	// U(old) (ndof)
	Up_old.setel(rho,1).setel(p,ndof);
	Up_old.is(1,vl_indx,vl_indxe).set(vel).rs();

	// Uc(old) (ndof)
	Uc_old.setel(rho,1).setel(rho_ene,ndof);
	Uc_old.is(1,vl_indx,vl_indxe).set(vel).scale(rho).rs();

	// A(old)  (ndim,ndof,ndof)
	compute_Ajac;

	// A(old) grad N (nel,ndof,ndof)
	Ao_grad_N.prod(Ajac,dshapex,-1,2,3,-1,1);

	// Cp(old) (ndof,ndof)
	compute_Cp_matrix;
	Cpjac_old.set(Cpjac);

	// strain rate tensor computation at old state 
	grad_u.set(grad_u_star);
	strain_rate.set(grad_u);
	grad_u.t();
	strain_rate.add(grad_u).scale(0.5);
	grad_u.rs();

	// LES correction of diffusion parameters ( !!! at old state !!!! for upwind and much more)
	double nu_eff,nu_tur=0;
	compute_visco_LES;

	// Tau_supg (ndof,ndof)
	visco_supg = nu_eff;
	cond_supg  = cond_eff;
	grad_rho.set(grad_rho_old);	
	compute_tau_supg;	
	tau_supg_m.eye(tau_supg_a).scale(tau_fac);	
	tau_supg_c.prod(Cpi,tau_supg_m,1,-1,-1,2);
	delta_sc_v.set(delta_sc).scale(shocap);

	// P_supg (nel,ndof,ndof)
	P_supg.prod(Ao_grad_N,tau_supg_c,1,2,-1,-1,3);

	// set state to current state
	rho=rho_star;p=p_star;vel.set(u_star);
	grad_rho.set(grad_rho_star);grad_p.set(grad_p_star);
	set_state;
	grad_U.ir(2,1).set(grad_rho).ir(2,ndof).set(grad_p).rs();
	grad_U.is(2,vl_indx,vl_indxe).set(grad_u_star).rs();

	// Up current (ndof)
	Up.setel(rho,1).setel(p,ndof);
	Up.is(1,vl_indx,vl_indxe).set(vel).rs();

	// Uc current (ndof)
	Uc.setel(rho,1).setel(rho_ene,ndof);
	Uc.is(1,vl_indx,vl_indxe).set(vel).scale(rho).rs();

	// A  (ndim,ndof,ndof)
	compute_Ajac;

	// A grad U (ndof)
	A_grad_N.prod(Ajac,dshapex,-1,2,3,-1,1);
	A_grad_U.prod(Ajac,grad_U,-1,1,-2,-1,-2);

	// Cp (ndof,ndof)
	compute_Cp_matrix;

	// Fa : advective fluxes (ndof,ndim)
	compute_Fa;
	
	// Djac (ndim,ndof,ndim,ndof)
	compute_Djac;
	
	// Fd : diffusive fluxes (ndof,ndim)
	compute_Fd;

	// shock capturing term
	compute_shocap;
	
	// dUc / dt (ndof)
	if (0) {
	// DEBUG
	Cpjac.eye();
	dUpdt.set(Up).rest(Up_old).scale(rec_Dt);
	dUcdt.prod(Cpjac,dUpdt,1,-1,-1);
	} else {
	dUcdt.set(Uc).rest(Uc_old).scale(rec_Dt);
	}
	// G (source) (ndof)
	G_source.set(0.0);

	// residual temporal term (Galerkin)
	tmp1_G.prod(SHAPE,dUcdt,1,2);

	// residual temporal term (SUPG)
	tmp1_S.prod(P_supg,dUcdt,1,2,-1,-1);

	// residual convective term (Galerkin)
	tmp2_G.prod(dshapex,flux,-1,1,2,-1);

	// residual convective term (SUPG)
	tmp2_S.prod(P_supg,A_grad_U,1,2,-1,-1);

	// residual diffusive term (Galerkin)

	if (0) {
	// DEBUG
	Djac.set(0.0);
	Djac.ir(1,1).ir(3,1).eye().rs();
	Djac.ir(1,2).ir(3,2).eye().rs();
	fluxd.prod(Djac,grad_U,2,1,-1,-2,-1,-2);
	// END DEBUG
	}

	tmp3_G.prod(dshapex,fluxd,-1,1,2,-1);

	// residual for external forces (Galerkin)
	tmp4_G.prod(SHAPE,G_source,1,2);

	// residual for external forces (Galerkin)
	tmp4_S.prod(P_supg,G_source,1,2,-1,-1);

	// residual for shock capturing term
	tmp5_S.prod(dshapex,flux_sc,-1,1,2,-1);

	res.set(tmp1_G).axpy(tmp1_S,1.0).axpy(tmp2_G,-1.0).axpy(tmp2_S,1.0);
	res.add(tmp3_G).rest(tmp4_G).rest(tmp4_S).add(tmp5_S);
	veccontr.axpy(res,-wpgdet);

	// matrices computations
	Ni_Nj.prod(SHAPE,SHAPE,1,2);      // (nel,nel)
	Pi_Nj.prod(P_supg,SHAPE,1,2,4,3); // (nel,ndof,nel,ndof)
	gNi_Nj.prod(dshapex,SHAPE,1,2,3); // (ndim,nel,nel)
	gNi_gNj.prod(dshapex,dshapex,1,2,3,4); // (ndim,nel,ndim,nel)
	gNi_dot_gNj.prod(dshapex,dshapex,-1,1,-1,2); // (nel,nel)

	// temporal Galerkin
	tmp1_G_M.prod(Ni_Nj,Cpjac,1,3,2,4).scale(rec_Dt);  // (nel,ndof,nel,ndof)

	// temporal SUPG
	tmp1_S_M.prod(Pi_Nj,Cpjac,1,2,3,-1,-1,4).scale(rec_Dt);  // (nel,ndof,nel,ndof)

	// convective Galerkin
	tmp2_G_M.prod(gNi_Nj,Ajac,-1,1,3,-1,2,4);  // (nel,ndof,nel,ndof)

	// convective SUPG  (nel,ndof,ndof) x (nel,ndof,ndof)
	tmp2_S_M.prod(P_supg,A_grad_N,1,2,-1,3,-1,4);  // (nel,ndof,nel,ndof)

	// diffusive Galerkin (ndim,nel,ndim,nel) x (ndim,ndof,ndim,ndof)
	tmp3_G_M.prod(gNi_gNj,Djac,-1,1,-2,3,-1,2,-2,4);  // (nel,ndof,nel,ndof)

	// shock capturing term
	Cpjac_old_delta.set(Cpjac_old);
	for (int jdf=1; jdf<=ndof; jdf++) {
	  Cpjac_old_delta.ir(1,jdf).scale(double(delta_sc_v.get(jdf)));
	}
	Cpjac_old_delta.rs();
	tmp5_S_M.prod(gNi_dot_gNj,Cpjac_old_delta,1,3,2,4);

	mat.set(tmp1_G_M).axpy(tmp1_S_M,1.0).axpy(tmp2_G_M,-1.0).axpy(tmp2_S_M,1.0);
	mat.add(tmp3_G_M).add(tmp5_S_M);
	matlocf.axpy(mat,wpgdet);


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
      
      // DEBUG
      if (0 && k % 17==0) {
	printf(" elemento: %d \n",k);   
	locstate.print("Estado :");	
	Cpjac.print("Cp: ");
	Ajac.print("Ajac: ");
	flux.print("Flux_A:");
	fluxd.print("Flux_D:");
	P_supg.print("P_supg: ");
      }
      
      // FIN DEBUG
      
      
      //      veccontr.set(res);
      
      //      veccontr.is(2,1,ndim).set(resmom)
      //	.rs().ir(2,ndof).set(rescont).rs();
      
      if (residual_factor!=1.)veccontr.scale(residual_factor);      
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
