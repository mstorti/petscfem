//__INSERT_LICENSE__
//$Id: nsitetlesf.cpp,v 1.1.4.1 2007/03/08 19:22:16 dalcinl Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "./nsi_tet.h"
#include "./nsitetlesf.h"

#define STANDARD_UPWIND
#define USE_FASTMAT

extern TextHashTable *GLOBAL_OPTIONS;

#define STOP {PetscFinalize(); exit(0);}
   
#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "nsi_tet_les_full::ask(char *,int &)"
int nsi_tet_les_full::
ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  // default
  DONT_SKIP_JOBINFO(comp_prof);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  DONT_SKIP_JOBINFO(get_nearest_wall_element);
  // special
  DONT_SKIP_JOBINFO(comp_mat_mass);
  DONT_SKIP_JOBINFO(comp_mat_advec);
  DONT_SKIP_JOBINFO(comp_mat_poisson);
  return 0;
}

#undef __FUNC__
#define __FUNC__ "nsi_tet_les_full::assemble"
int nsi_tet_les_full::
assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
	 Dofmap *dofmap,const char *jobinfo,int myrank,
	 int el_start,int el_last,int iter_mode,
	 const TimeData *time_) {
  // default
  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(get_nearest_wall_element);
  // special
  GET_JOBINFO_FLAG(comp_mat_mass);
  GET_JOBINFO_FLAG(comp_mat_advec);
  GET_JOBINFO_FLAG(comp_mat_poisson);

  int update_matrix=0;
  if(comp_mat_mass) {
    comp_res = 1;
    update_matrix = 1;
  }
  if(comp_mat_advec) {
    comp_res = 1;
    update_matrix = 1;
  }
  if(comp_mat_poisson) {
    comp_res = 1;
    update_matrix = 1;
  }

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
  // PetscPrintf(PETSCFEM_COMM_WORLD,"entrando a nsi_tet\n");

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)
#define NN_IDX(j) ELEMIPROPS_ADD(j,0)

  int locdof,kldof,lldof;
  char *value;

  // Unpack Dofmap
  int neq,nnod;
  neq  = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  
  //o Number of Dimensions.
  SGETOPTDEF(int,ndim,nodedata->ndim);

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,default);
  if(geometry == "default") {
    GPdata::get_default_geom(ndim, nel, geometry);
  }
  //o Number of Gauss points.
  TGETOPTDEF(thash,int,npg,-1);
  if (npg == -1) {
    GPdata::get_default_npg(geometry, npg);
  }

  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  int nen = nel*ndof;

  // Hloc stores the old mesh coordinates
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),vloc_mesh(nel,ndim),v_mesh(ndim);

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  // Get arguments from arg_list
  double *locst=NULL,*locst2=NULL,*retval=NULL,*retvalmat=NULL;
  WallData *wall_data=NULL;
  if (comp_prof) {
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
  GlobParam *glob_param=NULL;
  double *hmin=NULL,Dt=INFINITY,rec_Dt=0.0;
  int ja_hmin=INT_MAX;
#define WAS_SET arg_data_v[ja_hmin].was_set
  if (comp_mat_res) {
    int ja=0;
    locst = arg_data_v[ja++].locst;
    locst2 = arg_data_v[ja++].locst;
    retval = arg_data_v[ja++].retval;
    if (update_jacobian) 
      retvalmat = arg_data_v[ja++].retval;
    else if (update_matrix)   
      retvalmat = arg_data_v[ja++].retval;
    hmin = &*(arg_data_v[ja++].vector_assoc)->begin();
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
    wall_data = (WallData *)arg_data_v[ja++].user_data;
  } 

  //o Use Stokes form of NS equations. 
  SGETOPTDEF(int,stokes_form,0);
  //o Explicit treatement of advection velocity. 
  SGETOPTDEF(int,oseen_form,0);
  //o Use Laplace form of NS equations. 
  SGETOPTDEF(int,laplace_form,0);
  //o Use a weak form for the gradient of pressure term.
  SGETOPTDEF(int,weak_form,1);

  //o Add pressure controlling term. 
  SGETOPTDEF(double,pressure_control_coef,0.);
  assert(pressure_control_coef>=0.);

  //o Use full jacobian
  SGETOPTDEF(int,use_full_jacobian,1);

  //o Flag to *turn on* ALE computation
  SGETOPTDEF(int,ALE_flag,0);
  //o Pointer to old coordinates in
  //  #nodedata# array excluding the first #ndim# values
  SGETOPTDEF(int,indx_ALE_xold,1);
  //o Assert `fractional_step' is not used. 
  SGETOPTDEF(int,fractional_step,0);
  PETSCFEM_ASSERT0(!fractional_step,
                   "This elemset is to be used only \n"
                   "with the monolithic version. ");  

  // allocate local vecs
  int kdof;
  FastMat2 veccontr(2,nel,ndof),xloc(2,nel,ndim),locstate(2,nel,ndof), 
    locstate2(2,nel,ndof),xpg,G_body(1,ndim),F_body(1,ndim),
    vrel(1,ndim);

  if (ndof != ndim+1) {
    PetscPrintf(PETSCFEM_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
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
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Invalid value for \"axisymmetric\" option\n"
		"axisymmetric=\"%s\"\n",axisymmetric.c_str());
    PetscFinalize();
    exit(0);
  }

  //o Add LES for this particular elemset.
  SGETOPTDEF(int,LES,0);
  //o Smagorinsky constant.
  SGETOPTDEF(double,C_smag,0.18); // Dijo Beto
  //o van Driest constant for the damping law.
  SGETOPTDEF(double,A_van_Driest,0); 
  assert(A_van_Driest>=0.);
  //o print Van Driest factor
  SGETOPTDEF(int,print_van_Driest,0); 

  SGETOPTDEF(int,update_wall_data,0);
  assert(!((A_van_Driest>0)&&(update_wall_data==0)));

  //o Add LSIC stabilization
  SGETOPTDEF(int,use_lsic,1);
  //o Explicit treatement of SUPG and PSPG stabilization terms. 
  SGETOPTDEF(int,use_explicit_upwind,0);
  if (oseen_form) use_explicit_upwind = 1;

  //o Adjust the stability parameters, taking into account
  // the time step. If the  #steady#  option is in effect,
  // (which is equivalent to $\Dt=\infty$) then
  //  #tst_factor#  is set to 0.
  SGETOPTDEF(int,use_tst, 1);  // Scale upwind
  SGETOPTDEF(double,tst_factor,1.);  // Scale upwind
  if (comp_mat_res && glob_param->steady) use_tst=0;

  //o Scale the all the stabilization term. 
  SGETOPTDEF(double,tau_fac,1.);  // Scale upwind
  //o Scales the SUPG stabilization term. 
  SGETOPTDEF(double,tau_supg_fac,1.);  // Scale upwind
  //o Scales the PSPG stabilization term. 
  SGETOPTDEF(double,tau_pspg_fac,1.);  // Scale upwind
  //o Scales the LSIC stabilization term. 
  SGETOPTDEF(double,tau_lsic_fac,1.);  // Scale upwind
  //o Add to the  #tau_pspg#  term, so that you can stabilize with a term
  //  independently of $h$. (Mainly for debugging purposes). 
  SGETOPTDEF(double,additional_tau_supg,0.);  // Scale upwind
  SGETOPTDEF(double,additional_tau_pspg,0.);  // Scale upwind
  SGETOPTDEF(double,additional_tau_lsic,0.);  // Scale upwind
  double tau_supg_add = additional_tau_supg;
  double tau_pspg_add = additional_tau_pspg;
  double tau_lsic_add = additional_tau_lsic;

  //o Scale the residual term. 
  SGETOPTDEF(double,residual_factor,1.);
  //o Scale the jacobian term. 
  SGETOPTDEF(double,jacobian_factor,1.);

  //o XXX Write me
  SGETOPTDEF(int,body_force,0);
  

  double alpha=NAN;
  if (comp_mat_res) {
    alpha = glob_param->alpha;
    if (glob_param->steady) alpha=1.0;
  }

  //o _T: double[ndim] _N: G_body _D: null vector 
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.set(0.);
  ierr = get_double(thash,"G_body",
		    G_body.storage_begin(),1,ndim);

  double pi = 4*atan(1.0);

  //o Density
  TGETOPTDEF(thash,double,rho,1.);
  //o Viscosity
  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))
  int nprops=iprop;
  
  // Definiciones para descargar el lazo interno
  double detJaco, UU, u2, Peclet, psi, 
    tau_supg, tau_pspg, tau_lsic,
    div_u_star,p_star,wpgdet,velmod,tol,h_supg,fz,Uh;

  FastMat2 P_supg, W_supg, W_supg_t, dmatw,
    P_pspg(2,ndim,nel),dshapex(2,ndim,nel);
  FMatrix respert;

  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FMatrix Jaco(ndim,ndim),iJaco(ndim,ndim),
    grad_u(ndim,ndim),grad_u_star,strain_rate(ndim,ndim),resmom(nel,ndim),
    dresmom(nel,ndim),matij(ndof,ndof),Uintri,svec;

  FMatrix grad_p_star(ndim),u,u_star,du,
    uintri(ndim),rescont(nel),dmatu(ndim),ucols,ucols_new,
    ucols_star,pcol_star,pcol_new,pcol,fm_p_star,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,
    massm,tmp7,tmp8,tmp9,tmp10,tmp11,tmp13,tmp14,tmp15,dshapex_c,xc,
    wall_coords(ndim),dist_to_wall,tmp16,tmp162,tmp17,tmp18,tmp19,
    tmp23,tmp24,tmp25,
    tmp26,tmp27,tmp28,tmp29,tmp30;
  FastMat2 tmp20(2,nel,nel),tmp21,vel_supg(1,ndim);

  double pdyn; FMatrix grad_pdyn(ndim); 
  FMatrix tmp40,tmp41,tmp42,tmp43,tmp44,tmp45;

  double tmp12;
  double tsf = use_tst ? tst_factor : 0.0;

  FMatrix eye(ndim,ndim),seed,one_nel,matloc_prof(nen,nen);

  FMatrix Jaco_axi(2,2),u_axi;
  int ind_axi_1, ind_axi_2;
  double detJaco_axi;
         
  if (axi) assert(ndim==3);

  eye.eye();

  if (comp_prof) {

    matloc_prof.set(1.);

  }

  if (body_force) this->bf_init(nodedata);

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

    if ( A_van_Driest>0.&& get_nearest_wall_element ) {
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

    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_mat_res || comp_res) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }

    matlocmom.set(0.);
    matloc.set(0.);
    matlocf.set(0.);
    veccontr.set(0.);
    resmom.set(0.);
    rescont.set(0.);

    if (comp_res || comp_mat_res) {
      ucols.set(locstate2.is(2,1,ndim));
      pcol.set(locstate2.rs().ir(2,ndof));
      locstate2.rs();
      
      ucols_new.set(locstate.is(2,1,ndim));
      pcol_new.set(locstate.rs().ir(2,ndof));
      locstate.rs();
      
      ucols_star.set(ucols_new).scale(alpha).axpy(ucols,1-alpha);
      pcol_star.set(pcol_new).scale(alpha).axpy(pcol,1-alpha);
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
      vloc_mesh.set(xloc).minus(Hloc).scale(rec_Dt).rs();
      Hloc.rs();
    }
    
    double shear_vel=NAN;
    int wall_elem;
    if (LES && comp_mat_res && A_van_Driest>0.) {
#ifdef USE_ANN
      if (!wall_data) { set_error(2); return 1; }
      Elemset *wall_elemset;
      const double *wall_coords_;
      wall_data->nearest_elem_info(NN_IDX(k),wall_elemset,wall_elem,wall_coords_);
      wall_coords.set(wall_coords_);
      shear_vel = wall_elemset->elemprops_add[wall_elem];

      xc.sum(xloc,-1,1).scale(1./double(nel));
      dist_to_wall.set(xc).minus(wall_coords);

#else
      PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
    }

    if (body_force) this->bf_eval_el(k);

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

// #define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
// #define SHAPE    (*gp_data.FM2_shape[ipg])
// #define WPG      (gp_data.wpg[ipg])
// #define WPG_SUM  (gp_data.wpg_sum)

      FastMat2& DSHAPEXI = *gp_data.FM2_dshapexi[ipg];
      FastMat2& SHAPE    = *gp_data.FM2_shape[ipg];
      double    WPG      = gp_data.wpg[ipg];
      double    WPG_SUM  = gp_data.wpg_sum;

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
      if (ndim==1) {
	h_pspg = Delta = Area;
	//PFEMERRQ("Only dimensions 2 and 3 allowed for this element.\n");
      } else if (ndim==2) {
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
      
      if (comp_mat_res) {
	// computes the minimum size of the mesh
	if (!WAS_SET || h_pspg<*hmin) {
	  WAS_SET = 1;
	  *hmin = h_pspg;
	}

	// state variables and gradient
	u.prod(SHAPE,ucols,-1,-1,1);
	p_star = double(tmp8.prod(SHAPE,pcol_star,-1,-1));
	u_star.prod(SHAPE,ucols_star,-1,-1,1);

	grad_u.prod(dshapex,ucols,1,-1,-1,2);
	grad_u_star.prod(dshapex,ucols_star,1,-1,-1,2);
	grad_p_star.prod(dshapex,pcol_star,1,-1,-1);

	if (laplace_form) {
	  strain_rate.set(grad_u_star);
	  strain_rate.scale(0.5);
	} else {
	  strain_rate.set(grad_u_star);
	  grad_u_star.t();
	  strain_rate.add(grad_u_star);
	  grad_u_star.rs();
	  strain_rate.scale(0.5);
	}

	v_mesh.set(0.0);
	if (ALE_flag) {
	  v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	}

	// Smagorinsky turbulence model
	double mu_eff=NAN;
	if (LES) {
	  double van_D=NAN,ywall=NAN;
	  double tr = (double) tmp15.prod(strain_rate,strain_rate,-1,-2,-1,-2);
	  //	  double van_D;
	  if (A_van_Driest>0.) {
	    //	    dist_to_wall.prod(SHAPE,xloc,-1,-1,1).minus(wall_coords);
	    ywall = sqrt(dist_to_wall.sum_square_all());
	    double y_plus = ywall*shear_vel/(VISC/rho);
	    van_D = 1.-exp(-y_plus/A_van_Driest);
	  } else van_D = 1.;

	  double nu_t = SQ(C_smag*Delta*van_D)*sqrt(2*tr);
	  mu_eff = VISC + rho * nu_t;

	  if (print_van_Driest && (k % 1==0))
	    printf("element %d , y: %f,van_D: %f, nu_eff: %f\n",
		   ielh, ywall, van_D, mu_eff/rho);

	} else {
	  mu_eff = VISC;
	}

#if 1
	// XXX tau_supg and tau_pspg computed in u^n
	vel_supg.set(u).minus(v_mesh).rs();
#else
	if (use_explicit_upwind)
	  vel_supg.set(u).minus(v_mesh).rs();
	else
	  vel_supg.set(u_star).minus(v_mesh).rs();
#endif
        if (axi>0) { vel_supg.setel(0.,axi); }
	
	u2 = vel_supg.sum_square_all();
	
	velmod = sqrt(u2);
        tol=1.0e-16;
        h_supg=0;
	FastMat2::branch();
        if(velmod>tol) {
	  FastMat2::choose(0);
	  svec.set(vel_supg).scale(1./velmod);
	  h_supg = tmp9.prod(dshapex,svec,-1,1,-1).sum_abs_all();
          h_supg = ((h_supg < tol) ? tol : h_supg);
          h_supg = 2./h_supg;
        } else {
          h_supg = h_pspg;
        }
	FastMat2::leave();

	double nu_eff = mu_eff/rho;

        tau_supg = 
	  tsf*SQ(2.*rec_Dt) +
	  SQ(2.*velmod/h_supg) +
	  9.*SQ(4.*nu_eff/SQ(h_supg));
        tau_supg = 1./sqrt(tau_supg);

        tau_pspg = 
	  tsf*SQ(2.*rec_Dt) +
	  SQ(2.*velmod/h_pspg) +
	  9.*SQ(4.*nu_eff/SQ(h_pspg));
        tau_pspg = 1./sqrt(tau_pspg);

	Peclet = velmod * h_supg / (2. * nu_eff);
        fz = (Peclet < 3.) ? Peclet/3. : 1.;
        tau_lsic = 0.5*velmod*h_supg*fz;
	
	if (stokes_form) { 
	  tau_supg = 0.0;
	  tau_pspg = 0.0;
	  tau_pspg += tsf*SQ(2.*rec_Dt);
	  tau_pspg += 9.*SQ(4.*nu_eff/SQ(h_pspg));
	  tau_pspg = 1./sqrt(tau_pspg);
	}

	tau_supg *= (tau_fac * tau_supg_fac); tau_supg += tau_supg_add;
	tau_pspg *= (tau_fac * tau_pspg_fac); tau_pspg += tau_pspg_add;
	tau_lsic *= (tau_fac * tau_lsic_fac); tau_lsic += tau_lsic_add;

#if 1
	// XXX tau_supg and tau_pspg computed in u^n
	if (use_explicit_upwind)
	  vel_supg.set(u).minus(v_mesh).rs();
	else
	  vel_supg.set(u_star).minus(v_mesh).rs();
	if (axi>0) { vel_supg.setel(0.,axi); }
#endif

	if (stokes_form) { 
	  vrel.set(0.); 
	  vel_supg.set(0.);
	} else if (oseen_form) { 
	  vrel.set(u).minus(v_mesh).rs(); 
	} else { 
	  vrel.set(u_star).minus(v_mesh).rs(); 
	}
	
	// P_supg es un vector fila
	P_supg.prod(vel_supg,dshapex,-1,-1,1).scale(tau_supg);

	// Weight function, W = N + P_supg
	W_supg.set(P_supg).add(SHAPE);

	// Pressure stabilizing term
	P_pspg.set(dshapex).scale(tau_pspg/rho);  //debug:=

	// implicit version - General Trapezoidal rule - parameter alpha
	du.set(u_star).minus(u);
	dmatu.prod(vrel,grad_u_star,-1,-1,1);
	dmatu.axpy(du,rec_Dt/alpha).minus(G_body);

	if (body_force) {
	  F_body.set(0.);
	  this->bf_eval_pg(SHAPE,dshapex,F_body);
	  dmatu.minus(F_body);
	}
	
	
	// Galerkin - momentum
	dresmom.prod(SHAPE,dmatu,1,2);
	resmom.axpy(dresmom,-rho*wpgdet);

	if (weak_form) {
	  tmp1.set(strain_rate).scale(2*mu_eff).axpy(eye,-p_star);
	  tmp2.prod(dshapex,tmp1,-1,1,-1,2);
	  resmom.axpy(tmp2,-wpgdet);
	} else {
	  tmp6.prod(dshapex,strain_rate,-1,1,-1,2).scale(2*mu_eff);
	  tmp11.prod(SHAPE,grad_p_star,1,2).add(tmp6);
	  resmom.axpy(tmp11,-wpgdet);
	}

	// Galerkin - continuity
	div_u_star = double(tmp10.prod(dshapex,ucols_star,-1,-2,-2,-1));
	rescont.axpy(SHAPE,wpgdet*div_u_star);

	// SUPG/PSPG perturbation - residual
	tmp3.set(grad_p_star).axpy(dmatu,rho);
	respert.set(tmp3).scale(wpgdet);

	// SUPG perturbation - momentum
	tmp4.prod(P_supg,respert,1,2);
	resmom.minus(tmp4);

	// PSPG perturbation - continuity
	tmp5.prod(P_pspg,respert,-1,1,-1);
	rescont.add(tmp5);
	
        // LSIC perturbation - momentum
	if (use_lsic) {
	  resmom.axpy(dshapex.t(),-tau_lsic*rho*div_u_star*wpgdet);
	  dshapex.rs();
	}

	// Penalization term?
	if (pressure_control_coef) {
	  //double p_star = tmp21.prod(SHAPE,pcol_star,-1,-1).get();
	  rescont.axpy(SHAPE,pressure_control_coef*p_star*wpgdet);
	}

	if (update_jacobian) {

	  // temporal part + convective (Galerkin)
	  massm.prod(vrel,dshapex,-1,-1,1);
	  massm.axpy(SHAPE,rec_Dt/alpha);
	  matlocmom.prod(W_supg,massm,1,2).scale(rho*wpgdet);
	  
	  // diffusive part
	  tmp7.prod(dshapex,dshapex,-1,1,-1,2);
	  matlocmom.axpy(tmp7,mu_eff*wpgdet);

	  // dmatw =  rho * ((1/Dt)*SHAPE + u * dshapex);
	  dmatw.set(massm).scale(rho);

	  for (int ii=1; ii<=ndim; ii++) {
	    matlocf.ir(2,ii).ir(4,ii)
	      .add(matlocmom).rs();
	  }

	  if (!laplace_form) {
	    tmp19.set(dshapex).scale(mu_eff*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,3,4,1);
	    matlocf.is(2,1,ndim).is(4,1,ndim)
	      .add(tmp18).rs();
	  }

	  if (use_lsic) {
	    tmp19.set(dshapex).scale(tau_lsic*rho*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,1,4,3);
	    matlocf.is(2,1,ndim).is(4,1,ndim)
	      .add(tmp18).rs();
	  }

	  if (weak_form) {
	    tmp16.prod(P_supg,dshapex,1,2,3).scale(wpgdet);
	    tmp162.prod(dshapex,SHAPE,2,1,3).scale(-wpgdet);
	    matlocf.is(2,1,ndim).ir(4,ndof)
	      .add(tmp16).add(tmp162).rs();
	  } else {
	    tmp16.prod(W_supg,dshapex,1,2,3).scale(wpgdet);
	    matlocf.is(2,1,ndim).ir(4,ndof)
	      .add(tmp16).rs();
	  }

	  matlocf.ir(2,ndof).is(4,1,ndim);
	  tmp17.prod(P_pspg,dmatw,3,1,2).scale(wpgdet);
	  matlocf.minus(tmp17);
	  tmp17.prod(dshapex,SHAPE,3,2,1).scale(wpgdet);
	  matlocf.minus(tmp17);
	  matlocf.rs();

	  tmp13.prod(P_pspg,dshapex,-1,1,-1,2);
	  matlocf.ir(2,ndof).ir(4,ndof)
	    .axpy(tmp13,-wpgdet).rs();
	  // Penalization term?
	  if (pressure_control_coef) {
	    tmp20.prod(SHAPE,SHAPE,1,2);
	    tmp13.axpy(tmp20,pressure_control_coef);
	    matlocf.ir(2,ndof).ir(4,ndof)
	      .axpy(tmp13,-wpgdet).rs();
	  }


	  if (use_full_jacobian) {
	    // Jacobian term for advective term
	    if (!oseen_form && !stokes_form) {
	      tmp23.prod(SHAPE,grad_u_star,2,3,1);
	      // W * N  \dep u/ \dep x
	      tmp24.prod(W_supg,tmp23,1,2,3,4);
	      matlocf.is(2,1,ndim).is(4,1,ndim)
		.axpy(tmp24,rho*wpgdet).rs();
	      // P * N  \dep u/ \dep x
	      tmp25.prod(P_pspg,tmp23,-1,1,-1,2,3);
	      matlocf.ir(2,ndof).is(4,1,ndim)
		.axpy(tmp25,-rho*wpgdet).rs();
	    }
            // Jacobian term for SUPG perturbation
	    if (!use_explicit_upwind && !stokes_form) {
#if 1         // XXX tau_supg computed in u^n
	      tmp26.set(eye).scale(tau_supg);
#else         /* XXX we should try to compute der_tau_supg*/ 
	      double der_tau_supg = 0.0;
	      double mod_vel_supg = velmod?velmod:1.0;
	      tmp26.prod(vel_supg,vel_supg,1,2)
		.scale(der_tau_supg/mod_vel_supg);
	      tmp26.axpy(eye,tau_supg);
#endif
	      tmp27.prod(tmp26,dshapex,1,-1,-1,2);
	      tmp28.set(respert);
	      tmp29.prod(tmp27,tmp28,1,2,3);
	      tmp30.prod(SHAPE,tmp29,3,4,1,2);
	      matlocf.is(2,1,ndim).is(4,1,ndim)
	      	.add(tmp30).rs();
	    }
	  }
	}
	
	if (comp_mat_mass) {

	  tmp23.set(SHAPE).scale(wpgdet);
	  tmp24.prod(SHAPE,tmp23,1,2);
	  matlocmom.set(tmp24);

	  //double tau = 1./3.*SQ(h_pspg)/4.;
	  //tmp25.set(dshapex).scale(tau*wpgdet);
	  //tmp26.prod(dshapex,tmp25,-1,1,-1,2);
	  //matlocmom.add(tmp26);

	  for (int ii=1; ii<=ndof; ii++) {
	    matlocf.ir(2,ii).ir(4,ii)
	      .add(matlocmom).rs();
	  }
	}

	if (comp_mat_advec) {

	  // temporal part + convective (Galerkin)
	  massm.prod(vrel,dshapex,-1,-1,1);
	  massm.axpy(SHAPE,rec_Dt/alpha);
	  matlocmom.prod(SHAPE,massm,1,2).scale(rho*wpgdet);
	  // diffusive part
	  tmp7.prod(dshapex,dshapex,-1,1,-1,2);
	  matlocmom.axpy(tmp7,mu_eff*wpgdet);

	  for (int ii=1; ii<=ndof; ii++) {
	    matlocf.ir(2,ii).ir(4,ii)
	      .add(matlocmom).rs();
	  }
	}

	if (comp_mat_poisson) {

	  tmp25.set(dshapex).scale(wpgdet);
	  tmp26.prod(dshapex,tmp25,-1,1,-1,2);
	  matlocmom.set(tmp26);
	  
	  for (int ii=1; ii<=ndof; ii++) {
	    matlocf.ir(2,ii).ir(4,ii)
	      .add(matlocmom).rs();
	  }
	}

      } else if (comp_prof) {
	// don't make anything here !!
      } else {
	PetscPrintf(PETSCFEM_COMM_WORLD,
		    "Don't know how to compute jobinfo: %s\n",jobinfo);
	CHKERRQ(ierr);
      }

    }

    if(comp_prof) {
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
    }      

    if (comp_mat_res) {

      veccontr.is(2,1,ndim).set(resmom)
	.rs().ir(2,ndof).set(rescont).rs();

      if (residual_factor!=1.)
	veccontr.scale(residual_factor);
      veccontr.export_vals(&(RETVAL(ielh,0,0)));

      if (update_jacobian) {
	if (jacobian_factor!=1.)
	  matlocf.scale(jacobian_factor);
	matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      }

      if (update_matrix) {
	matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      }
      
#ifdef PRINT_ELEM_DEBUG
      if (k==0) {
	veccontr.print("veccontr:");
	matlocf.print("matlocf:");
      }
#endif
    }
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();

  return 0;
}

// #undef SHAPE    
// #undef DSHAPEXI 
// #undef WPG      
// #undef WPG_SUM
