//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "fracstep.h"

#define MAXPROP 100

// extern WallData wall_data;

const double FIX = 1.0;
/** Fixes all diagonal terms in matrix #A# so that they are #>0#, i.e.,
    makes #A(i,i)=FIX# if #A(i,i)=0# and leaves #A(i,i)# unaltered otherwise. 
    At the same time makes #B# a diagonal matrix which is the mask added to #A#. 
    In brief: B(i,i)=0; if (A(i,i)==0) A(i,i)=B(i,i)=FIX;
    @param A (input) matrix to be fixed.
    @param B (input) mask for A
*/ 
static void 
fix_null_diagonal_entries(FastMat2 &A,FastMat2 &B,int &nc) {
  int n = A.dim(1);
  B.set(0.);
  for (int j=1; j<=n; j++) {
    if (A.get(j,j)==0.) {
      A.setel((nc==-1? 1.0 : ++nc),j,j);
      B.setel(FIX,j,j);
    }
  }
}

static void 
fix_null_diagonal_entries(FastMat2 &A,FastMat2 &B) {
  int nc=-1;
  fix_null_diagonal_entries(A,B,nc);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define fracstep fracstep_fm2
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This is a function that approximates the absolute function
// and has the nice properties of being infinite differentiable
// and has monotone increasing derivative. It has a removable
// singularity at x=0
double smabs(double x,double &dydx) {
  double y, tol=1e-7;
  if (fabs(x)<tol) {
    // Remove singularity using a polynomial expansion
    y = 1.0+x*x/3.0;
    dydx = (2.0/3.0)*x;
  } else {
    // Use the complete expression
    double 
      sh = sinh(x),
      th = tanh(x);
    y = x/th;
    dydx = 1.0/th - x/(sh*sh);
  }
  return y;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This is a reduced version that only returns the function
// not the derivativae
double smabs(double x) {
  double yd;
  return smabs(x,yd);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "fracstep::assemble"
int fracstep::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		       Dofmap *dofmap,const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_) {

  assert(fractional_step);
  int ierr=0, axi;

  GET_JOBINFO_FLAG(comp_mat_prof);
  GET_JOBINFO_FLAG(comp_res_mom);
  GET_JOBINFO_FLAG(comp_mat_poi);
  GET_JOBINFO_FLAG(comp_res_poi);
  GET_JOBINFO_FLAG(comp_mat_prj);
  GET_JOBINFO_FLAG(comp_res_prj);
  GET_JOBINFO_FLAG(get_nearest_wall_element);

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
#define ELEMIPROPS_ADD(j,k) VEC2(elemiprops_add,j,k,neliprops_add)
#define NN_IDX(j) ELEMIPROPS_ADD(j,0)
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

  // Unpack Dofmap
  int nnod = dofmap->nnod;

  // Hloc stores:
  // for ALE: the old mesh coordinates
  // or scalar fields for other couplings
  int nH = nu-ndim;
  FMatrix  Hloc(nel,nH),vloc_mesh(nel,ndim),v_mesh(ndim),Tloc(nel),T_star;

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  double *locst=NULL,*locst2=NULL,*retval=NULL,
    *retvalmat=NULL,*retvalmat_mom=NULL,*retvalmat_poi=NULL,
    *retvalmat_prj=NULL;

  WallData *wall_data_p=NULL;
  if (get_nearest_wall_element) {
    wall_data_p = (WallData *)arg_data_v[0].user_data;
    if(!wall_data_p) {
      printf("Null 'wall_data' object found.\n");
      set_error(2);
      return 1;
    }
  }
    
  // rec_Dt is the reciprocal of Dt (i.e. 1/Dt)
  // for steady solutions it is set to 0. (Dt=inf)
  GlobParam *glob_param=NULL;
  double Dt=NAN,rec_Dt=NAN;
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
    rec_Dt = 1./Dt;
    if (glob_param->steady) rec_Dt=0.;
    wall_data_p = (WallData *)(arg_data_v[ja++].user_data);
  } else if (comp_mat_poi) {
    int ja=0;
    A_poi_arg = &arg_data_v[ja];
    retvalmat_poi = arg_data_v[ja++].retval;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
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
    if (!reuse_mat) {
      A_prj_arg = &arg_data_v[ja];
      retvalmat_prj = arg_data_v[ja].retval;
      ja++;
    }
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    Dt = glob_param->Dt;
  } else assert(0); // Not implemented yet!!

  FastMat2 veccontr(2,nel,ndof),xloc(2,nel,ndim),
    locstate(2,nel,ndof),locstate2(2,nel,ndof),tmp(2,nel,ndof),
    ustate2(2,nel,ndim),G_body(1,ndim),vrel,gravity(1,ndim);

  if (ndof != ndim+1) {
    PetscPrintf(PETSCFEM_COMM_WORLD,"ndof != ndim+1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  FastMat2 matloc(4,nel,ndof,nel,ndof), matlocmom(2,nel,nel), masspg(2,nel,nel),
    matlocmom2(4,nel,ndof,nel,ndof), grad_u_ext(2,ndof,ndof),
    mom_profile(4,nel,ndof,nel,ndof), poi_profile(4,nel,ndof,nel,ndof),
    mom_mat_fix(2,nen,nen), poi_mat_fix(2,nen,nen);

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


  //o ALE_flag : flag to ON ALE computation
  SGETOPTDEF(int,ALE_flag,0);
  //o indx_ALE_xold : pointer to old coordinates in
  //  NODEDATA array excluding the first "ndim" values
  SGETOPTDEF(int,indx_ALE_xold,1);

  //o Bouyancy_flag : coupling fractional step with thermal equation
  SGETOPTDEF(int,BUOYANCY_flag,0);
  //o indx_buoyancy_T : pointer to temperatures in
  //  NODEDATA array excluding the first "ndim" values
  SGETOPTDEF(int,indx_buoyancy_T,0);

  if(BUOYANCY_flag>0)assert(indx_buoyancy_T>0);

  //o buoyancy coefficient for Boussinesq term
  SGETOPTDEF(double,betath,0.);
  //o Reference temperature for thermal coupling Boussinesq term
  SGETOPTDEF(double,Tinfty,0.);
  //o Gravity acceleration for buoyancy terms
  gravity.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"gravity",gravity.storage_begin(),1,ndim);

  //o Axis for selective Darcy term (damps incoming flow
  //at outlet bdry's)
  SGETOPTDEF(int,darcy_axi,0); 
  double axi_sign = 1.0;
  if (darcy_axi<0) {
    axi_sign = -1.0;
    darcy_axi = -darcy_axi;
  }
  assert(darcy_axi<=ndim);

  //o Reference velocity for selectiv Darcy term. 
  SGETOPTDEF(double,darcy_uref,-1.0); 
  if (darcy_axi) { assert(darcy_uref>=0.); }

  //o Reference velocity for selectiv Darcy term. 
  SGETOPTDEF(double,darcy_factor_global,-1.0); 
  if (darcy_axi) { assert(darcy_factor_global>=0.); }

#if 0
  //o Coefficient affecting selective Darcy term. 
  SGETOPTDEF(double,darcy_coef,0.0); 
  assert(darcy_coef>=0. && darcy_axi>0);
#endif

  double alpha=0.5, gammap=0.0;
  ierr = get_double(thash,"alpha",&alpha,1); CHKERRA(ierr);
  ierr = get_double(thash,"gamma_pressure",&gammap,1); CHKERRA(ierr);
  int weak_poisson = 1;
  ierr = get_int(thash,"weak_poisson",&weak_poisson,1); CHKERRA(ierr);

  // Factor para la estabilizacion
  double taufac=1.0;
  ierr = get_double(thash,"taufac",&taufac,1); CHKERRA(ierr);

  DEFPROP(viscosity);
#define VISC (*(propel+viscosity_indx))
  DEFPROP(darcy_factor);
#define DARCY (*(propel+darcy_factor_indx))

  int nprops=iprop;

  //o _T: double[ndim] _N: G_body _D: null vector 
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.set(0.);
  ierr = get_double(GLOBAL_OPTIONS,"G_body",G_body.storage_begin(),1,ndim);
  
  //o Density
  TGETOPTDEF(thash,double,rho,1.);
  //o Factor masking the fractional step stabilization term. 
  TGETOPTDEF(thash,double,dt_art_fac,1.);
  //o characteristic velocity for pressure stabilization 
  TGETOPTDEF(thash,double,U_char,0.);
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);

  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  FastMat2 grad_fi,P_supg,W;
  //  LogAndSign L;
  double detJaco, UU, u2, Peclet, psi, tau, div_u_star,
    wpgdet;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  FastMat2 dshapex(2,ndim,nel),Jaco(2,ndim,ndim),iJaco(2,ndim,ndim),
    grad_u(2,ndim,ndim),grad_u_star(2,ndim,ndim),dshapext(2,nel,ndim),
    resmom(2,nel,ndim), fi(1,ndof), grad_p(1,ndim), grad_p_star(1,ndim),
    u(1,ndim),u_star(1,ndim),uintri(1,ndim),rescont(1,nel);
  FastMat2 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp71,tmp72,tmp73,tmp8,tmp9,tmp10,
    tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,xc,
    wall_coords(1,ndim),dist_to_wall;
  FastMat2 vel_supg;

  FMatrix Jaco_axi(2,2),u_axi,strain_rate(ndim,ndim);
  int ind_axi_1, ind_axi_2;
  double detJaco_axi;         
  if (axi) assert(ndim==3);

  masspg.set(1.);
  grad_u_ext.set(0.);
  grad_u_ext.is(1,1,ndim).is(2,1,ndim);
  if (couple_velocity) grad_u_ext.set(1.);
  else grad_u_ext.eye();
  grad_u_ext.rs();

  // grad_u_ext was transposed in the Newmat version. Why? if it's symmetric
#if 0
  mom_profile.prod(masspg,grad_u_ext,1,3,4,2);
#else
  mom_profile.set(0.);
  for (int j=1; j<=ndim; j++) {
    double w = j;
    mom_profile.ir(2,j).ir(4,j).set(w);
  }
  mom_profile.rs();
#endif

  int nc = ndim;
  mom_profile.reshape(2,nen,nen);
  fix_null_diagonal_entries(mom_profile,mom_mat_fix,nc);
  mom_profile.reshape(4,nel,ndof,nel,ndof);
  mom_mat_fix.reshape(4,nel,ndof,nel,ndof);
  
  grad_u_ext.set(0.);
  grad_u_ext.setel(1.,ndim+1,ndim+1);
  poi_profile.prod(masspg,grad_u_ext,1,3,4,2);
  poi_profile.reshape(2,nen,nen);
  fix_null_diagonal_entries(poi_profile,poi_mat_fix);
  poi_profile.reshape(4,nel,ndof,nel,ndof);
  poi_mat_fix.reshape(4,nel,ndof,nel,ndof);

  if (comp_mat_prof) {
    mom_profile.export_vals(arg_data_v[0].profile); // A_mom
    poi_profile.export_vals(arg_data_v[1].profile); // A_poi
    mom_profile.export_vals(arg_data_v[2].profile); // A_prj
  } else if (comp_res_mom) {
    mom_profile.export_vals(A_mom_arg->profile);
    matloc.set(0.).add(mom_mat_fix);
  } else if (comp_mat_poi) {
    poi_profile.export_vals(A_poi_arg->profile);
  } else if (comp_res_poi) {
  } else if (comp_mat_prj) {
    mom_profile.export_vals(A_prj_arg->profile);
  } else if (comp_res_prj) {
    if (!reuse_mat) mom_profile.export_vals(A_prj_arg->profile);
  } else assert(0);

  FastMat2 seed;
  if (comp_res_mom || comp_mat_prj || (comp_res_prj && !reuse_mat)) {
    seed.resize(2,ndof,ndof).set(0.)
      .is(1,1,ndim).is(2,1,ndim).eye().rs();
  } else if (comp_mat_poi) {
    seed.resize(2,ndof,ndof).set(0.).
      setel(1.,ndof,ndof).rs();
  }

  double PFUNUSED max_u_neg = 0., max_force=0.;
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  int SHV_debug=0;
#undef SHV
#define SHV(pp) { if (SHV_debug) pp.print(#pp); }
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    //if (epart[k] != myrank+1) continue;
    ielh++;

    if (comp_mat_prof) {
      // We export anything, because in fact `upload_vector'
      // doesn't even look at the `retvalmat' values. It looks at the
      // .profile member in the argument value
      matlocmom.export_vals(&(RETVALMAT_MOM(ielh)));
      matlocmom.export_vals(&(RETVALMAT_PRJ(ielh)));
      matlocmom.export_vals(&(RETVALMAT_POI(ielh)));
      continue;
    } 

    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    //    printf("element %d, prop %f\n",k,ELEMPROPS(k,0));
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
      wall_data_p->nearest(xc.storage_begin(),nn);
      NN_IDX(k) = nn;
      continue;
#else
      PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
    }

    // tenemos el estado locstate2 <- u^n
    //                   locstate  <- u^*
    if (comp_res_mom || comp_res_poi || comp_res_prj) {
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));
    }
    matlocmom.set(0.);
    matlocmom2.set(0.);
    veccontr.set(0.);
    resmom.set(0.);
    rescont.set(0.);

    if (comp_res_mom) {
      // nodal computation of mesh velocity
      if (ALE_flag) {
	assert(nH >= ndim);
	assert(indx_ALE_xold >= nH+1-ndim);
	Hloc.is(2,indx_ALE_xold,indx_ALE_xold+ndim-1);
	vloc_mesh.set(xloc).minus(Hloc).scale(rec_Dt).rs();
	Hloc.rs();
      } else if (BUOYANCY_flag) {
	assert(nH >= 1); // at least one field for temperature. 
	assert(indx_buoyancy_T >= 1);
	Hloc.ir(2,indx_buoyancy_T);
	Tloc.set(Hloc).rs();
	Hloc.rs();
      }
    }

    double shear_vel=NAN;
    if (comp_res_mom) {
      int wall_elem;
      if (LES && A_van_Driest>0.) {
#ifdef USE_ANN
	if (!wall_data_p) { set_error(2); return 1; }
	Elemset *wall_elemset;
	const double *wall_coords_;
	wall_data_p->nearest_elem_info(NN_IDX(k),wall_elemset,wall_elem,wall_coords_);
	wall_coords.set(wall_coords_);
	shear_vel = wall_elemset->elemprops_add[wall_elem];

	xc.sum(xloc,-1,1).scale(1./double(nel));
	dist_to_wall.set(xc).minus(wall_coords);
	
#else
	PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
      }
      
    }
    

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])
#define WPG_SUM  (gp_data.wpg_sum)

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      //      Matrix xpg = SHAPE * xloc;
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);

      elem = k;

      detJaco = Jaco.det();
      if (detJaco<=0.) {
	detj_error(detJaco,elem);
	set_error(1);
      }
      wpgdet = detJaco*WPG;
      iJaco.inv(Jaco);
      dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);

      double Area = detJaco*WPG_SUM, Delta;
      if (ndim==2) {
	Delta = sqrt(Area);
      } else if (ndim==3 && axi==0) {
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
	Delta = sqrt(Area_axi);
      } else {
	PFEMERRQ("Only dimensions 2 and 3 allowed for this element.\n");
      }

      double u2 ,Uh;      

      if (comp_res_mom) {
	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// PREDICTOR STEP
	// state variables and gradient
	locstate2.is(2,1,ndim);
	u.prod(SHAPE,locstate2,-1,-1,1);
	grad_u.prod(dshapex,locstate2,1,-1,-1,2);
	locstate2.rs();

	locstate.is(2,1,ndim);
	u_star.prod(SHAPE,locstate,-1,-1,1);
	grad_u_star.prod(dshapex,locstate,1,-1,-1,2);
	locstate.rs();

	locstate2.ir(2,ndof);
	grad_p.prod(dshapex,locstate2,1,-1,-1);
	locstate2.rs();

	v_mesh.set(0.0);
	if (ALE_flag) {
	  v_mesh.prod(SHAPE,vloc_mesh,-1,-1,1);
	}

	if (BUOYANCY_flag) {
	  T_star.prod(SHAPE,Tloc,-1,-1);
	}

	// Smagorinsky turbulence model
	double nu_eff;
	double van_D=NAN,ywall=NAN;
	if (LES) {

	strain_rate.set(grad_u_star);
	grad_u_star.t();
	strain_rate.add(grad_u_star).scale(0.5);
	grad_u_star.rs();

#if 0
	double tr = (double) tmp15.prod(strain_rate,strain_rate,-1,-2,-1,-2);
	//	  double van_D;
	if (A_van_Driest>0.) {
	  //	  dist_to_wall.prod(SHAPE,xloc,-1,-1,1).minus(wall_coords);
	  //	    double ywall = sqrt(dist_to_wall.sum_square_all());
	  ywall = sqrt(dist_to_wall.sum_square_all());
	  double y_plus = ywall*shear_vel/VISC;
	  van_D = 1.-exp(-y_plus/A_van_Driest);
	  //	    if (k % 250==0) printf("van_D: %f\n",van_D);
	  //	    if (k % 499==0) printf("van_D: %f\n",van_D);
	} else van_D = 1.;
	  
	double nu_t = SQ(C_smag*Delta*van_D)*sqrt(2*tr);
	nu_eff = VISC + nu_t;
#else
	double tr = (double) tmp15.prod(strain_rate,strain_rate,-1,-2,-1,-2);
	//	double van_D;
	if (A_van_Driest>0.) {
	  // dist_to_wall.prod(SHAPE,xloc,-1,-1,1).minus(wall_coords);
	  ywall = sqrt(dist_to_wall.sum_square_all());
	  double y_plus = ywall*shear_vel/VISC;
	  van_D = 1.-exp(-y_plus/A_van_Driest);
	  //	  if (k % 250==0) printf("van_D: %f\n",van_D);
	} else van_D = 1.;
	
	double nu_t = SQ(C_smag*Delta*van_D)*sqrt(2*tr);
	nu_eff = VISC + nu_t;
#endif
	} else {
	  nu_eff = VISC;
	}

	//	if (print_van_Driest && (ywall<0.02)) 
	if (print_van_Driest) 
	  printf("element %d , y: %f,van_D: %f, nu_eff: %f\n",
		 ielh, ywall, van_D,nu_eff);

	if (0) {
	  u2 = u.sum_square_all(); 
	} else {
	  vel_supg.set(u).minus(v_mesh).rs();
	  if(axi>0){
	    vel_supg.setel(0.,axi);
	  }
	  u2 = vel_supg.sum_square_all();
	}

	if(0) {	
	  if(axi>0){
	    u_axi.set(u);
	    u_axi.setel(0.,axi);
	    u2 = u_axi.sum_square_all();
	    uintri.prod(iJaco,u_axi,1,-1,-1);
	    Uh = sqrt(uintri.sum_square_all())/2.;
	  } else {
	    uintri.prod(iJaco,u,1,-1,-1);
	    Uh = sqrt(uintri.sum_square_all())/2.;
	  }
	  
	  if(u2<=1e-6*(2.*Uh*nu_eff)) {
	    Peclet=0.;
	    psi=0.;
	    tau=0.;
	  } else {
	    //	  Peclet = u2 / (2. * Uh * VISC);
	    Peclet = u2 / (2. * Uh * nu_eff);
	    psi = 1./tanh(Peclet)-1/Peclet;
	    tau = psi/(2.*Uh);
	  }
	  P_supg.prod(u,dshapex,-1,-1,1).scale(taufac * tau);
	} else {

	  u_axi.set(vel_supg);
	  uintri.prod(iJaco,u_axi,1,-1,-1);
	  Uh = sqrt(uintri.sum_square_all())/2.;
	  
	  if(u2<=1e-6*(2.*Uh*nu_eff)) {
	    Peclet=0.;
	    psi=0.;
	    tau=0.;
	  } else {
	    Peclet = u2 / (2. * Uh * nu_eff);
	    psi = 1./tanh(Peclet)-1/Peclet;
	    tau = psi/(2.*Uh);
	  }
	  P_supg.prod(vel_supg,dshapex,-1,-1,1).scale(taufac * tau);
	}

	W.set(SHAPE).add(P_supg);
	
	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// RESIDUE CALCULATION
	if(0) {
	tmp1.prod(u,grad_u,-1,-1,1).scale(-(1-alpha));
	tmp2.prod(u_star,grad_u_star,-1,-1,1).scale(-alpha);
	} else {
	vrel.set(u).minus(v_mesh).rs();
	tmp1.prod(vrel,grad_u,-1,-1,1).scale(-(1-alpha));
	vrel.set(u_star).minus(v_mesh).rs();
	tmp2.prod(vrel,grad_u_star,-1,-1,1).scale(-alpha);
	}
	tmp1.add(tmp2).axpy(grad_p,-(gammap/rho));
	tmp3.prod(W,tmp1,1,2);
	resmom.axpy(tmp3,wpgdet);
	SHV(resmom);
	// sacamos gradiente de presion en la ec. de momento (conf. Codina)
	//- (1/rho)* grad_p.t());

	// Parte difusiva
	// version implicita
	// resmom -= wpgdet * VISC * dshapex.t() * grad_u_star;
	// version Crank-Nicholson 
	tmp4.set(grad_u).scale(1-alpha).axpy(grad_u_star,alpha);
	tmp5.prod(dshapex,tmp4,-1,1,-1,2);
	//	resmom.axpy(tmp5,-wpgdet*VISC);
	resmom.axpy(tmp5,-wpgdet*nu_eff);
	SHV(resmom);
	
	// Parte temporal
	tmp6.set(u_star).minus(u);
	tmp7.prod(W,tmp6,1,2);
	resmom.axpy(tmp7,-wpgdet/Dt);
	SHV(resmom);

	// Fuerza de cuerpo
	tmp71.prod(W,G_body,1,2);
	resmom.axpy(tmp71,wpgdet);

	// Bouyancy force
	if (BUOYANCY_flag) {	  
	  tmp72.set(gravity).scale(rho*betath);
          tmp73.prod(W,tmp72,1,2);
          resmom.axpy(tmp73,-(double(T_star)-Tinfty)*wpgdet);
	}
	
	// Selective Darcy term
	if (darcy_axi) {
	  // Velocity along `axi' direction
	  double uu = u.get(darcy_axi)*axi_sign;
	  // Smoothed u^+
	  double au = smabs(uu/darcy_uref)*darcy_uref;
	  // Force acting in direction positive when
	  // velocity comes in negative direction. 
	  double darcy = DARCY*darcy_factor_global;
	  double force = axi_sign*darcy*(au-uu)/2.0;
	  resmom.ir(2,darcy_axi).axpy(W,wpgdet*force).rs();
	  if (force>max_force) {
	    max_force = force;
	    max_u_neg = uu;
	  }
#if 0
	  if (ipg==0 && k%100==0) 
	    printf("k %d, uu %f, au %f, darcy %f, force %f, axi_sign %f\n",
		   k,uu,au,darcy,force,axi_sign);
#endif
	}

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// JACOBIAN CALCULATION
	if(0){
	  tmp17.prod(u_star,dshapex,-1,-1,1);
	} else {
	  tmp17.prod(vrel,dshapex,-1,-1,1);
	}
	tmp8.prod(W,tmp17,1,2);
	matlocmom.axpy(tmp8,alpha*wpgdet);

	masspg.prod(W,SHAPE,1,2);
	tmp9.prod(masspg,grad_u_ext,1,3,2,4);
	if (couple_velocity)
	  matlocmom2.is(2,1,ndim).is(4,1,ndim).axpy(tmp9,alpha * wpgdet).rs();

	// Parte difusiva
	tmp10.prod(dshapex,dshapex,-1,1,-1,2);
	//	matlocmom.axpy(tmp10,alpha*wpgdet*VISC);
	matlocmom.axpy(tmp10,alpha*wpgdet*nu_eff);
	
	// Parte temporal
	matlocmom.axpy(masspg,wpgdet/Dt);
       
      } else if (comp_res_poi || comp_mat_poi) {

	double dt_fac = 1.0;
	if (dt_art_fac > 0.0) {
	  double area = detJaco * gp_data.wpg_sum;
	  double h = pow(area,1.0/double(ndim));
	  dt_fac += dt_art_fac/Dt*(1. /(4.*VISC/h/h+2.*U_char/h)); 
	  //	  dt_fac += dt_art_fac*h*h/VISC;
	} 

	if (comp_res_poi) {
	  locstate.ir(2,ndof);
	  grad_p.prod(dshapex,locstate,1,-1,-1);
	  locstate.rs();
	  
	  // The Poisson equation is treated in weak form. This
	  // guarantees that the solution of the Poisson equation is
	  // independent of which node is fixed its pressure. 
	  locstate2.is(2,1,ndim);
	  u_star.prod(SHAPE,locstate2,-1,-1,1);
	  locstate2.rs();
	  tmp12.set(u_star).scale(-rho/Dt).axpy(grad_p,dt_fac-gammap);
	  tmp11.prod(dshapex,tmp12,-1,1,-1);
	  rescont.axpy(tmp11,-wpgdet);
#if 0
	  // No-weak form
	  div_u_star = (dshapex * ustate2).Trace();
	  rescont -= wpgdet * ((rho/Dt_art) * SHAPE.t() * div_u_star
			       + dshapex.t() * grad_p);
#endif
	}

	if (comp_mat_poi) {
	  tmp13.prod(dshapex,dshapex,-1,1,-1,2);
	  matlocmom.axpy(tmp13,dt_fac*wpgdet);
	}
 
      } else if (comp_res_prj) {

	locstate.ir(2,ndof);
	grad_p_star.prod(dshapex,locstate,1,-1,-1);
	locstate.rs();

	locstate2.ir(2,ndof);
	grad_p.prod(dshapex,locstate2,1,-1,-1);

	locstate2.rs().is(2,1,ndim);
	u_star.prod(SHAPE,locstate2,-1,-1,1);
	locstate2.rs();

	locstate.is(2).is(2,1,ndim);
	u.prod(SHAPE,locstate,-1,-1,1);
	locstate.rs();

	tmp14.set(u_star).minus(u).axpy(grad_p_star,-Dt/rho);
	if (gammap) tmp14.axpy(grad_p,+gammap*Dt/rho);
	tmp15.prod(SHAPE,tmp14,1,2);
	resmom.axpy(tmp15,wpgdet);

	SHV(resmom);
	SHV(SHAPE);
	SHV(u);
	SHV(u_star);

      } else if (comp_mat_prj || (comp_res_prj && !reuse_mat) ) {

	tmp16.prod(SHAPE,SHAPE,1,2);
	matlocmom.axpy(tmp16,wpgdet);

      } else {

	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);

      }

    }
    if (comp_res_mom) {
      veccontr.is(2,1,ndim).set(resmom);
      veccontr.rs().export_vals(&(RETVAL(ielh)));
      // matloc.prod(matlocmom,seed,1,3,2,4).add(mom_mat_fix);
      for (int j=1; j<=ndim; j++)
	matloc.ir(2,j).ir(4,j).set(matlocmom);
      matloc.rs();
      matloc.export_vals(&(RETVALMAT(ielh)));
    } else if (comp_mat_poi) {
      matloc.prod(matlocmom,seed,1,3,2,4).add(poi_mat_fix);
      matloc.export_vals(&(RETVALMAT_POI(ielh)));
    } else if (comp_res_poi) {
      veccontr.ir(2,ndof).set(rescont);
      veccontr.rs().export_vals(&(RETVAL(ielh)));
    } else if (comp_mat_prj) {
      matloc.prod(matlocmom,seed,1,3,2,4).add(mom_mat_fix);
      matloc.export_vals(&(RETVALMAT_PRJ(ielh)));
    } else if (comp_res_prj) {
      veccontr.is(2,1,ndim).set(resmom);
      veccontr.rs().export_vals(&(RETVAL(ielh)));
      if (comp_res_prj && !reuse_mat) {
        matloc.prod(matlocmom,seed,1,3,2,4).add(mom_mat_fix);
        matloc.export_vals(&(RETVALMAT_PRJ(ielh)));
      }
    } else assert(0);
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
#if 0
  if (comp_res_mom) {
    vector<double> maxunv(SIZE), maxfv(SIZE);
    MPI_Gather(&max_force,1, MPI_DOUBLE,
	       &maxfv[0],1,MPI_DOUBLE,0,PETSCFEM_COMM_WORLD);
    MPI_Gather(&max_u_neg,1,MPI_DOUBLE,
	       &maxunv[0],1,MPI_DOUBLE,0,PETSCFEM_COMM_WORLD);
    if (!myrank) {
      int jmax = 0;
      for (int j=0; j<SIZE; j++) 
	if (maxfv[j]>maxfv[jmax]) jmax=j;
      printf("DARCY stat: max(-u) %f, max(force) %f\n",
	     maxunv[jmax],maxfv[jmax]);
    }
  }
#endif
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
