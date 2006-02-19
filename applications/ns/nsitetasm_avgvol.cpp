//__INSERT_LICENSE__
//$Id: nsitetasm_avgvol.cpp,v 1.3 2006/02/19 23:59:47 mstorti Exp $

//
// This elemset uses the volume average velocity model 
// instead of mass average one. We keep v_m as the velocity
// but it should be taken as volume average.
//

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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double compute_rho_m_2(FastMat2 &rho_g_vp,FastMat2 &arho_g_vp, FastMat2 &alpha_g_vp, 
		       double &alpha_l, double rho_l, int nphases, FastMat2 &grad_alpha_g_vp, 
		       FastMat2 &grad_rho_m,const int ndim) {

  static FastMat2 Id_vp(2,nphases,nphases),grad_alpha_l(1,ndim),tmp(2,ndim,nphases); 
  
  double alpha_d_sum = (double) alpha_g_vp.sum_all();
  alpha_l = 1.0 - alpha_d_sum;
  
  Id_vp.set(0.).d(1,2);
  Id_vp.set(rho_g_vp).rs();	
  arho_g_vp.prod(Id_vp,alpha_g_vp,1,-1,-1);
  double arho_l = rho_l*alpha_l;
  double rho_m = arho_l+arho_g_vp.sum_all();

  tmp.prod(grad_alpha_g_vp,Id_vp,1,-1,-1,2);
  grad_rho_m.sum(tmp,1,-1);

  grad_alpha_l.sum(grad_alpha_g_vp,1,-1).scale(-1.0);
  grad_rho_m.axpy(grad_alpha_l,rho_l);

  return rho_m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double compute_rho_m_2_bcc(FastMat2 &rho_g_vp,FastMat2 &arho_g_vp, FastMat2 &alpha_g_vp, 
		       double &alpha_l, double rho_l, int nphases,const int ndim) {

  //  static FastMat2 Id_vp(2,nphases,nphases);
  
  double alpha_d_sum = (double) alpha_g_vp.sum_all();
  alpha_l = 1.0 - alpha_d_sum;

  arho_g_vp.set(alpha_g_vp).mult(rho_g_vp);
  
  double arho_l = rho_l*alpha_l;
  double rho_m = arho_l+arho_g_vp.sum_all();

  return rho_m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void compute_vel_g_2(const FastMat2 &V_v, FastMat2 &vslip_vp, const FastMat2 &vslip_user_vp, 
		     FastMat2 &vslip_m_vp, const FastMat2 &rho_g_vp,const double rho_m,
		     const double rho_l,const int g_dir,const FastMat2 &d_bubble_vp,
		     const int nphases,const int use_modified_slag_vslip,const double alpha_l,
		     const FastMat2 &alpha_g_vp,FastMat2 &v_g_vp,FastMat2 &grad_v_g_vp,
		     const FastMat2 &grad_rho_m,FastMat2 &grad_alpha_g_vp,const int ndim,
		     const FastMat2 &grad_V_v) {

  double rho_g,alpha_g,vslip;
  static FastMat2 grad_vslip_vp(2,ndim,nphases),V_m(1,ndim),grad_V_m(2,ndim,ndim),
    grad_Vm_Vv(2,ndim,ndim);
  static FastMat2 grad_vslip_l(1,ndim),grad_alpha_l(1,ndim);
  static FastMat2 tmp1_g_vp(1,nphases),tmp2_g_vp(1,nphases),tmp3_g_vp(1,nphases),
    tmp_g_vp(1,nphases),tmp0_g_vp(1,nphases);

  if (vslip_user_vp.sum_abs_all()>0.0) {
    vslip_vp.set(vslip_user_vp);
  } else {
    for (int j=1; j<=nphases; j++) {
      double rb = d_bubble_vp.get(j);
      vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
	       rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));
      
      vslip = (g_dir > 0 ? vslip : -vslip);
      vslip_vp.setel(vslip,j);
    }
  }

  // modifico velocidad slip de la escoria por la diferencia de densidades con la mezcla

  grad_vslip_vp.set(0.0);
  vslip_m_vp.set(vslip_vp);
  if (use_modified_slag_vslip) {
    rho_g = rho_g_vp.get(nphases);
    vslip = double(vslip_m_vp.get(nphases));
    double tmp = vslip*(1.-rho_g/rho_m);
    vslip_m_vp.setel(tmp,nphases);
    grad_vslip_vp.ir(2,nphases).set(grad_rho_m).scale(rho_g*vslip/rho_m/rho_m);
  }

  // Vm , Vl_slip 

  V_m.set(V_v);
  double vslip_l, tmp_vslip_l=0;

  for (int j=1; j<=nphases; j++) {
    rho_g = rho_g_vp.get(j);
    alpha_g = alpha_g_vp.get(j);
    grad_alpha_g_vp.ir(2,j);
    grad_vslip_vp.ir(2,j);
    vslip = vslip_m_vp.get(j);

    V_m.addel(-alpha_g*vslip*(1.0-rho_g/rho_l),abs(g_dir));
    tmp_vslip_l = tmp_vslip_l + rho_g*alpha_g*vslip;


  }
  grad_alpha_g_vp.rs();
  grad_vslip_vp.rs();
  grad_V_m.rs();

  // disperse phase velocity
  v_g_vp.set(0.);
  for (int j=1; j<=nphases; j++) {
    v_g_vp.ir(2,j);
    v_g_vp.set(V_m).addel(vslip_m_vp.get(j),abs(g_dir));
    v_g_vp.rs();
  }
  // continuum phase velocity
  v_g_vp.ir(2,nphases+1);
  vslip_l = -tmp_vslip_l/rho_l/alpha_l; 
  v_g_vp.set(V_m).addel(vslip_l,abs(g_dir)).rs();

  // continuous phase

  tmp0_g_vp.set(1.0).axpy(rho_g_vp,-1.0/rho_l);
  tmp1_g_vp.set(alpha_g_vp).mult(rho_g_vp);
  tmp2_g_vp.set(vslip_m_vp).mult(rho_g_vp);

  grad_vslip_l.set(0.);
  grad_Vm_Vv.set(0.);

  grad_Vm_Vv.ir(1,abs(g_dir));

  for (int j=1; j<=ndim; j++) {
    grad_vslip_vp.ir(1,j);
    grad_alpha_g_vp.ir(1,j);

    tmp_g_vp.set(tmp1_g_vp).mult(grad_vslip_vp);		   
    grad_vslip_l.addel(-tmp_g_vp.sum_all()/rho_l/alpha_l,j);

    tmp_g_vp.set(tmp2_g_vp).mult(grad_alpha_g_vp);		   
    grad_vslip_l.addel(-tmp_g_vp.sum_all()/rho_l/alpha_l,j);

    tmp_g_vp.set(tmp0_g_vp).mult(alpha_g_vp).mult(grad_vslip_vp);		   
    grad_Vm_Vv.addel(-tmp_g_vp.sum_all(),j);

    tmp_g_vp.set(tmp0_g_vp).mult(vslip_m_vp).mult(grad_alpha_g_vp);		   
    grad_Vm_Vv.addel(-tmp_g_vp.sum_all(),j);

  }
  grad_alpha_g_vp.rs();
  grad_vslip_vp.rs();
  grad_Vm_Vv.rs();

  grad_alpha_l.sum(grad_alpha_g_vp,1,-1).scale(-1.0);
  //  tmp3_g_vp.set(alpha_g_vp).mult(rho_g_vp).mult(vslip_m_vp);
  grad_vslip_l.axpy(grad_alpha_l,tmp_vslip_l/rho_l/alpha_l/alpha_l);

  // gradient of phase velocities
  grad_V_m.set(grad_V_v).add(grad_Vm_Vv);
  // disperse phases
  for (int j=1; j<=nphases; j++) {
    grad_vslip_vp.ir(2,j);
    grad_v_g_vp.ir(3,j).set(grad_V_m);
    grad_v_g_vp.ir(2,abs(g_dir)).add(grad_vslip_vp).rs();
  }
  grad_vslip_vp.rs();

  // continuum phase
  grad_v_g_vp.ir(3,nphases+1).set(grad_V_m);
  grad_v_g_vp.ir(2,abs(g_dir)).add(grad_vslip_l).rs();
  
}

/*
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void compute_vel_g_2(const FastMat2 &V_v, FastMat2 &vslip_vp, const FastMat2 &vslip_user_vp, 
		     FastMat2 &vslip_m_vp, const FastMat2 &rho_g_vp,const double rho_m,
		     const double rho_l,const int g_dir,const FastMat2 &d_bubble_vp,
		     const int nphases,const int use_modified_slag_vslip,const double alpha_l,
		     const FastMat2 &alpha_g_vp,FastMat2 &v_g_vp,FastMat2 &grad_v_g_vp,
		     const FastMat2 &grad_rho_m,FastMat2 &grad_alpha_g_vp,const int ndim,
		     const FastMat2 &grad_V_v) {

  double rho_g,alpha_g,vslip;
  static FastMat2 grad_vslip_vp(2,ndim,nphases),V_m(1,ndim),grad_V_m(2,ndim,ndim); 
  static FastMat2 grad_vslip_l(1,ndim),grad_alpha_l(1,ndim);
  static FastMat2 tmp1_g_vp(1,nphases),tmp2_g_vp(1,nphases),tmp3_g_vp(1,nphases),
    tmp_g_vp(1,nphases);

  if (vslip_user_vp.sum_abs_all()>0.0) {
    vslip_vp.set(vslip_user_vp);
  } else {
    for (int j=1; j<=nphases; j++) {
      double rb = d_bubble_vp.get(j);
      vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
	       rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));
      
      vslip = (g_dir > 0 ? vslip : -vslip);
      vslip_vp.setel(vslip,j);
    }
  }
  grad_vslip_vp.set(0.0);
  vslip_m_vp.set(vslip_vp);
  // modifico velocidad slip de la escoria por la diferencia de densidades con la mezcla
  if (use_modified_slag_vslip) {
    rho_g = rho_g_vp.get(nphases);
    vslip = double(vslip_m_vp.get(nphases));
    double tmp = vslip*(1.-rho_g/rho_m);
    vslip_m_vp.setel(tmp,nphases);
    grad_vslip_vp.ir(2,nphases).set(grad_rho_m).scale(rho_g*vslip/rho_m/rho_m);
  }
  // mass averaged mixture velocity
  V_m.set(V_v);
  grad_V_m.set(grad_V_v);
  double vslip_l=0;
  for (int j=1; j<=nphases; j++) {
    rho_g = rho_g_vp.get(j);
    alpha_g = alpha_g_vp.get(j);
    grad_alpha_g_vp.ir(2,j);
    grad_vslip_vp.ir(2,j);
    vslip = vslip_m_vp.get(j);
    V_m.addel(-alpha_g*vslip*(1.0-rho_g/rho_l),abs(g_dir));
    vslip_l = vslip_l - rho_g*alpha_g*vslip;
    //    grad_V_m.ir(2,abs(g_dir)).axpy(grad_alpha_g_vp,-vslip*(1.0-rho_g/rho_l));
    //    grad_V_m.ir(2,abs(g_dir)).axpy(grad_vslip_vp,-alpha_g*(1.0-rho_g/rho_l));
  }
  grad_alpha_g_vp.rs();
  grad_vslip_vp.rs();
  grad_V_m.rs();

  // disperse phase velocity
  v_g_vp.set(0.);
  for (int j=1; j<=nphases; j++) {
    v_g_vp.ir(2,j);
    v_g_vp.set(V_m).addel(vslip_m_vp.get(j),abs(g_dir));
    v_g_vp.rs();
  }
  // continuum phase velocity
  v_g_vp.ir(2,nphases+1);
  vslip_l = vslip_l/rho_l/alpha_l; 
  v_g_vp.set(V_m).addel(vslip_l,abs(g_dir)).rs();

  // gradient of phase velocities
  // disperse phases
  for (int j=1; j<=nphases; j++) {
    grad_vslip_vp.ir(2,j);
    grad_v_g_vp.ir(3,j).set(grad_V_m);
    grad_v_g_vp.ir(2,abs(g_dir)).add(grad_vslip_vp).rs();
  }
  grad_vslip_vp.rs();

  // continuous phase
  tmp1_g_vp.set(alpha_g_vp).mult(rho_g_vp);
  tmp2_g_vp.set(vslip_m_vp).mult(rho_g_vp);
  grad_vslip_l.set(0.);
  for (int j=1; j<=ndim; j++) {
    grad_vslip_vp.ir(1,j);
    grad_alpha_g_vp.ir(1,j);
    tmp_g_vp.set(tmp1_g_vp).mult(grad_vslip_vp);		   
    grad_vslip_l.addel(-tmp_g_vp.sum_all()/rho_l/alpha_l,j);
    tmp_g_vp.set(tmp2_g_vp).mult(grad_alpha_g_vp);		   
    grad_vslip_l.addel(-tmp_g_vp.sum_all()/rho_l/alpha_l,j);
  }
  grad_alpha_g_vp.rs();
  grad_vslip_vp.rs();

  grad_alpha_l.sum(grad_alpha_g_vp,1,-1).scale(-1.0);
  tmp3_g_vp.set(alpha_g_vp).mult(rho_g_vp).mult(vslip_m_vp);
  grad_vslip_l.axpy(grad_alpha_l,tmp3_g_vp.sum_all()/rho_l/alpha_l/alpha_l);

  grad_v_g_vp.ir(3,nphases+1).set(grad_V_m);
  grad_v_g_vp.ir(2,abs(g_dir)).add(grad_vslip_l).rs();
  
}
*/

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void compute_vel_g_2_bcc(const FastMat2 &V_v, FastMat2 &vslip_vp, const FastMat2 &vslip_user_vp, 
		     FastMat2 &vslip_m_vp, const FastMat2 &rho_g_vp,const double rho_m,
		     const double rho_l,const int g_dir,const FastMat2 &d_bubble_vp,
		     const int nphases,const int use_modified_slag_vslip,const double alpha_l,
		     const FastMat2 &alpha_g_vp,FastMat2 &v_g_vp,const int ndim) {

  double rho_g,alpha_g,vslip;
  static FastMat2 V_m(1,ndim);

  if (vslip_user_vp.sum_abs_all()>0.0) {
    vslip_vp.set(vslip_user_vp);
  } else {
    for (int j=1; j<=nphases; j++) {
      double rb = d_bubble_vp.get(j);
      vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
	       rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));
      
      vslip = (g_dir > 0 ? vslip : -vslip);
      vslip_vp.setel(vslip,j);
    }
  }
  vslip_m_vp.set(vslip_vp);
  // modifico velocidad slip de la escoria por la diferencia de densidades con la mezcla
  if (use_modified_slag_vslip) {
    rho_g = rho_g_vp.get(nphases);
    vslip = double(vslip_m_vp.get(nphases));
    double tmp = vslip*(1.-rho_g/rho_m);
    vslip_m_vp.setel(tmp,nphases);
  }
  // mass averaged mixture velocity
  V_m.set(V_v);
  double vslip_l=0;
  for (int j=1; j<=nphases; j++) {
    rho_g = rho_g_vp.get(j);
    alpha_g = alpha_g_vp.get(j);
    vslip = vslip_m_vp.get(j);
    V_m.addel(-alpha_g*vslip*(1.0-rho_g/rho_l),abs(g_dir));
    vslip_l = vslip_l - rho_g*alpha_g*vslip;
  }

  // disperse phase velocity
  v_g_vp.set(0.);
  for (int j=1; j<=nphases; j++) {
    v_g_vp.ir(2,j);
    v_g_vp.set(V_m).addel(vslip_m_vp.get(j),abs(g_dir));
    v_g_vp.rs();
  }
  // continuum phase velocity
  v_g_vp.ir(2,nphases+1);
  vslip_l = vslip_l/rho_l/alpha_l; 
  v_g_vp.set(V_m).addel(vslip_l,abs(g_dir)).rs();

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "nsi_tet_asm_avgvol::assemble"

int nsi_tet_asm_avgvol::
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
    dresmom(nel,ndim),matij(ndof,ndof),Uintri,svec(ndim),stress_tensor(ndim,ndim),
    grad_u_tmp(ndim,ndim);

  FMatrix grad_p_star(ndim),u,u_star,du,
    uintri(ndim),rescont(nel),dmatu(ndim),ucols,ucols_new,
    ucols_star,pcol_star,pcol_new,pcol,fm_p_star,tmp1,tmp2,tmp3,tmp6,
    massm,tmp8,tmp9,tmp10,tmp11,tmp13,tmp14,tmp15,dshapex_c,xc,
    wall_coords(ndim),dist_to_wall,tmp16,tmp162,tmp17,tmp171,tmp172,
    tmp173,tmp174,tmp18,tmp19;

  FastMat2 tmp20(2,nel,nel),tmp21;

  FastMat2 vfcols,vfcols_new,vfcols_star,Id_vp(2,nphases,nphases);

  FastMat2 vf,vf_star,grad_vf,grad_vf_star,grad_rho_m(1,ndim),grad_u_vp,
    grad_u_star_vp;
  FastMat2 res_alpha_g(2,nel,nphases),du_g(1,nphases),dmatu2,du_vp,dmatu_vp,
    dresmom_g(2,nel,nphases),matloc_alpha_g(4,nel,nphases,nel,nphases),T_extra_1;
  FastMat2 tmp1_g,tmp11_g,tmp12_g,tmp13_g,tmp14_g,tmp2_g,tmp3_g,tmp4_g,tmp6_g,tmp7_g,
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
  double d_bubble,visco_m_eff,visco_t,visco_g,visco_g_eff,visco_l_eff;
  vector<int> alpha_indx_vp;
  int vl_indx = 1;
  int vl_indxe = vl_indx+ndim-1;
  int alpha_indx = ndim+2;
  alpha_indx_vp.resize(nphases);
  for (int j=0; j<nphases; j++) alpha_indx_vp[j] = vl_indxe+j+1;

  FastMat2 alpha_g_vp,arho_g_vp,u_star_vp,u_vp;
  FastMat2 v_rel,v_rel_vp,vslip_vp,vslip_m_vp;

  alpha_g_vp.resize(1,nphases);
  arho_g_vp.resize(1,nphases);

  u_vp.resize(2,ndim,nphases+1);
  u_star_vp.resize(2,ndim,nphases+1);
  du_vp.resize(2,ndim,nphases+1);
  dmatu_vp.resize(2,ndim,nphases+1);

  grad_u_star_vp.resize(3,ndim,ndim,nphases+1);
  grad_u_vp.resize(3,ndim,ndim,nphases+1);

  T_extra_1.resize(1,ndim);

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

      /*
      // compute continuum and disperse phase velocities nodally 
      // at current time step ===> ucols_star_vp(nel,ndim,nphases+1)
      compute_vel_g(ucols_star,vfcols_star,rho_g_vp,rho_l,vslip_vp,
		    vslip_user_vp,vslip_m_vp,g_dir,d_bubble_vp,nphases,
		    use_modified_slag_vslip,ucols_star_vp);
      
      // compute continuum and disperse phase velocities nodally 
      // at old time step ===> ucols_vp(nel,ndim,nphases+1)
      compute_vel_g(ucols,vfcols,rho_g_vp,rho_l,vslip_vp,
		    vslip_user_vp,vslip_m_vp,g_dir,d_bubble_vp,nphases,
		    use_modified_slag_vslip,ucols_vp);

      */

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

	if (detJaco_axi<=0.) {
	  detj_error(detJaco_axi,elem);
	  set_error(1);
	}
	
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
	grad_p_star.prod(dshapex,pcol_star,1,-1,-1);
	grad_vf_star.prod(dshapex,vfcols_star,1,-1,-1,2);

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 


	// compute mixture density rho_m_old (old state)
	rho_m_old = compute_rho_m_2(rho_g_vp,arho_g_vp,vf,alpha_l,rho_l,nphases,grad_vf,
				grad_rho_m,ndim);

	// compute continuum and disperse phase velocities nodally 
	// at old time step ===> ucols_vp(nel,ndim,nphases+1)
	compute_vel_g_2(u,vslip_vp,vslip_user_vp,vslip_m_vp,rho_g_vp,rho_m,
			rho_l,g_dir,d_bubble_vp,nphases,use_modified_slag_vslip,alpha_l,
			vf,u_vp,grad_u_vp,grad_rho_m,grad_vf,ndim,grad_u);       

	alpha_g_vp.set(vf_star);

	// compute mixture density rho_m
	rho_m = compute_rho_m_2(rho_g_vp,arho_g_vp,alpha_g_vp,alpha_l,rho_l,nphases,grad_vf_star,
				grad_rho_m,ndim);

	// compute continuum and disperse phase velocities nodally 
	// at current time step ===> ucols_star_vp(nel,ndim,nphases+1)
	compute_vel_g_2(u_star,vslip_vp,vslip_user_vp,vslip_m_vp,rho_g_vp,rho_m,
			rho_l,g_dir,d_bubble_vp,nphases,use_modified_slag_vslip,alpha_l,
			vf_star,u_star_vp,grad_u_star_vp,grad_rho_m,grad_vf_star,ndim,grad_u_star);

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
	
	visco_l_eff = visco_l+rho_l*visco_t;
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
	
	velmod_g_vp.set(0.);
	for (int j=1; j<=nphases; j++) {
	  u_vp.ir(2,j);
	  u_axi_g_vp.ir(2,j);
	  vel_axi(u_vp,u_axi_g_vp,axi);
	  u2_g = u_axi_g_vp.sum_square_all();
	  velmod_g_vp.setel(sqrt(u2_g),j);
	}
	u_axi_g_vp.rs();
	u_vp.rs();

	h_supg = compute_h_supg(u_axi,dshapex,velmod,h_pspg);
	
	// upwind for disperse phases
	for (int j=1; j<=nphases; j++) {
	  u_axi_g_vp.ir(2,j);
	  velmod_g = velmod_g_vp.get(j);
	  h_supg_g = compute_h_supg(u_axi_g_vp,dshapex,velmod_g,h_pspg);
	  h_supg_vp.setel(h_supg_g,j);
	  u_axi_g_vp.rs();
	}
	
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

	// SUPG perturbation function for disperse phases
	// velocidad del gas	
	P_supg_vp.set(0.);
	for (int j=1; j<=nphases; j++) {
	  u_vp.ir(2,j);
	  P_supg_vp.ir(2,j);
	  h_supg_g = h_supg_vp.get(j);
	  visco_g = visco_g_eff_vp.get(j);
	  velmod_g = velmod_g_vp.get(j);
	  tau_supg_g = SQ(2.*velmod_g/h_supg_g)+9.*SQ(4.*visco_g/SQ(h_supg_g));
	  tau_supg_g = 1./sqrt(tau_supg_g);	  

	  tau_supg_g *= tau_fac;

	  P_supg_vp.prod(u_vp,dshapex,-1,-1,1).scale(tau_supg_g);

	  Peclet = velmod_g * h_supg_g / (2. * visco_g);

	  fz = (Peclet < 3. ? Peclet/3. : 1.);
	  delta_supg_g = 0.5*h_supg_g*velmod_g*fz;

	  delta_supg_g *= shock_capturing_factor;

	  delta_supg_vp.setel(delta_supg_g,j);

	}
	u_vp.rs();
	P_supg_vp.rs();

	// implicit version - General Trapezoidal rule - parameter alpha
#ifdef ADD_GRAD_DIV_U_TERM
	dmatu.prod(u_star,grad_u_star,-1,-1,1);
#else
	dmatu.prod(u,grad_u_star,-1,-1,1);
#endif

	// momentum equation residue	
	du.set(u_star).rest(u);
	dmatu.axpy(du,rec_Dt/alpha);
	dmatu2.set(dmatu).scale(rho_l).axpy(G_body,-rho_m);
       	resmom_prime.set(grad_p_star).add(dmatu2);
	// divergence of  velocity
	div_u_star = double(tmp10.prod(dshapex,ucols_star,-1,-2,-2,-1));

	// Galerkin - momentum
	// resmom tiene que tener nel*ndim
	dresmom.prod(SHAPE,dmatu2,1,2).rs();
	resmom.axpy(dresmom,-wpgdet);


	if(1) {
	stress_tensor.set(0.0);
	// disperse phase stress tensor
	for (int j=1; j<=nphases; j++) {
	  grad_u_tmp.set(grad_u_star_vp.ir(3,j));
	  strain_rate.set(grad_u_tmp);
	  grad_u_tmp.t();
	  strain_rate.add(grad_u_tmp).scale(0.5);
	  grad_u_star_vp.rs();
	  visco_g = visco_g_vp.get(j);
	  alpha_g = alpha_g_vp.get(j);
	  stress_tensor.axpy(strain_rate,2*visco_g*alpha_g);
	}
	// continuum phase stress tensor
	grad_u_tmp.set(grad_u_star_vp.ir(3,nphases+1));
	strain_rate.set(grad_u_tmp);
	grad_u_tmp.t();
	strain_rate.add(grad_u_tmp).scale(0.5);
	grad_u_star_vp.rs();
	stress_tensor.axpy(strain_rate,2*visco_l_eff*alpha_l);
	} else {
	  stress_tensor.set(strain_rate).scale(2*visco_m_eff);
	}

	if (weak_form) {
	  //	  tmp1.set(strain_rate).scale(2*visco_m_eff).axpy(eye,-p_star);
	  tmp1.set(stress_tensor).axpy(eye,-p_star);
	  tmp2.prod(dshapex,tmp1,-1,1,-1,2);
	  resmom.axpy(tmp2,-wpgdet);
	} else {
	  //	  tmp6.prod(dshapex,strain_rate,-1,1,-1,2).scale(2*visco_m_eff);
	  tmp6.prod(dshapex,stress_tensor,-1,1,-1,2);
	  tmp11.prod(SHAPE,grad_p_star,1,2).add(tmp6);
	  resmom.axpy(tmp11,-wpgdet);
	}

	// SUPG perturbation - momentum
	resmom_supg.prod(P_supg,resmom_prime,1,2);
	resmom.axpy(resmom_supg,-wpgdet);

	// shock capturing term - momentum
	resmom.axpy(dshapex.t(),-wpgdet*delta_supg*rho_l*div_u_star);
	dshapex.rs();

	// Galerkin - continuity
	rescont_gal.prod(dshapex,u_star,-1,1,-1);
	rescont.axpy(rescont_gal,-wpgdet);

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
	matlocmom.prod(W_supg,massm,1,2).scale(rho_l);

	gNi_gNj.prod(dshapex,dshapex,-1,1,-1,2);
	Ni_Nj.prod(SHAPE,SHAPE,1,2);
	gNi_Nj.prod(dshapex,SHAPE,1,2,3);

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
	// disperse phases part
	matloc_alpha_g.set(0.);
	du_g.set(vf_star).rest(vf).scale(rec_Dt/alpha);
	tmp11_g.prod(dshapex,grad_vf_star,-1,1,-1,2);

	for (int j=1; j<=nphases; j++) {

	  res_alpha_g.ir(2,j);
	  P_supg_vp.ir(2,j);
	  u_star_vp.ir(2,j);
	  grad_vf_star.ir(2,j);
	  tmp11_g.ir(2,j);
	  dresmom_g.ir(2,j);

	  // SUPG term
	  double dmatu_g = double(tmp13_g.prod(u_star_vp,grad_vf_star,-1,-1));
	  dmatu_g += double(du_g.get(j));
	  dresmom_g.set(P_supg_vp).scale(dmatu_g);

	  // Galerkin terms
	  dresmom_g.axpy(SHAPE,double(du_g.get(j)));
	  tmp1_g.prod(dshapex,u_star_vp,-1,1,-1).scale(double(alpha_g_vp.get(j)));
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
	  tmp2_g.prod(dshapex,u_star_vp,-1,1,-1);
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
	  
	}
	tmp11_g.rs();
	res_alpha_g.rs();
	P_supg_vp.rs();
	u_star_vp.rs();
	grad_vf_star.rs();
	matloc_alpha_g.rs();
	dresmom_g.rs();

	if (update_jacobian) {
	  matlocf.is(2,ndim+2,ndim+1+nphases).is(4,ndim+2,ndim+1+nphases);
	  matlocf.axpy(matloc_alpha_g,wpgdet).rs();
	}

	// TE_1 : first extra term (propto Dvk/Dt)

	// compute the material derivative of each disperse phase velocity
	dmatu_vp.set(0.);
	for (int j=1; j<=nphases; j++) {
	  //	  u_vp.ir(2,j);
	  u_star_vp.ir(2,j);
	  grad_u_star_vp.ir(3,j);
	  dmatu_vp.ir(2,j);

 	  dmatu_vp.prod(u_star_vp,grad_u_star_vp,-1,-1,1);
	}

	//	u_vp.rs();
	u_star_vp.rs();
	grad_u_star_vp.rs();
	dmatu_vp.rs();

	du_vp.set(u_star_vp).rest(u_vp);
	dmatu_vp.axpy(du_vp,rec_Dt/alpha);

	for (int j=1; j<=nphases; j++) {

	  rho_g = rho_g_vp.get(j);
	  alpha_g = alpha_g_vp.get(j);
	  dmatu_vp.ir(2,j);

	  T_extra_1.set(dmatu_vp).scale((rho_g-rho_l)*alpha_g);

	  // Galerkin and SUPG contribution 
          tmp14_g.prod(W_supg,T_extra_1,1,2);
	  resmom.axpy(tmp14_g,-wpgdet);
	  // PSPG contribution
          tmp12_g.prod(P_pspg,T_extra_1,-1,1,-1);
	  rescont.axpy(tmp12_g,wpgdet);
	  
	  if (update_jacobian && coupled) {

	    matlocf.is(2,1,ndim).ir(4,ndim+1+j);
	    // Galerkin and SUPG contribution
	    tmp9_g.prod(W_supg,SHAPE,1,2);
	    tmp10_g.prod(tmp9_g,dmatu_vp,1,3,2);
	    matlocf.axpy(tmp10_g,wpgdet*(rho_g-rho_l)).rs();

	    // PSPG contribution --> dR_cont(PSPG)/d_alpha_g
	    tmp12_g.prod(P_pspg,dmatu_vp,-1,1,-1);
	    tmp9_g.prod(tmp12_g,SHAPE,1,2);
	    matlocf.ir(2,ndim+1).ir(4,ndim+1+j);
	    matlocf.axpy(tmp9_g,-wpgdet*(rho_g-rho_l)).rs();
	    
	  }
	}

	/*
	// TE_2 : second extra term 
	T_extra_2.set(0.);
	div_u_star_vp = double(tmpTE2_0.prod(dshapex,ucols_star_vp,-1,-2,-2,-1,1));

	for (int j=1; j<=nphases; j++) {
	  u_star_vp.ir(2,j);
	  grad_u_star_vp.ir(3,j);
	  div_u_star_vp.ir(1,j);
	  alpha_g = alpha_g_vp.get(j);
	  grad_vf_star.ir(2,j);

 	  tmpTE2_1.prod(u_star_vp,grad_u_star_vp,-1,-1,1);
	  tmpTE2_1.axpy(u_star_vp,double(div_u_star_vp)).scale(alpha_g);

 	  tmpTE2_2.prod(u_star_vp,grad_vf_star,-1,-1);
 	  tmpTE2_1.axpy(u_star_vp,double(tmpTE2_2));

	  T_extra_2.add(tmpTE2_1);

	}
	u_star_vp.rs();
	grad_u_star_vp.rs();
	div_u_star_vp.rs();
	grad_vf_star.rs();
	
	// for continuum phase
	
	u_star_vp.ir(2,nphases+1);
	grad_u_star_vp.ir(3,nphases+1);
	div_u_star_vp.ir(1,nphases+1);
	grad_alpha_l.set(sum(grad_vf_star,1,-1));
	
	tmpTE2_1.prod(u_star_vp,grad_u_star_vp,-1,-1,1);
	tmpTE2_1.axpy(u_star_vp,double(div_u_star_vp)).scale(alpha_l);
	
	tmpTE2_2.prod(u_star_vp,grad_vf_star,-1,-1);
	tmpTE2_1.axpy(u_star_vp,double(tmpTE2_2));
	
	T_extra_2.add(tmpTE2_1);
	
	// Galerkin and SUPG contribution 
	tmp11_g.prod(W_supg,T_extra_1,1,2);
	resmom.axpy(tmp11_g,-wpgdet);
	// PSPG contribution
	tmp12_g.prod(P_pspg,T_extra_1,-1,1,-1);
	rescont.axpy(tmp12_g,wpgdet);
	
	if (update_jacobian && coupled) {
	  
	  matlocf.is(2,1,ndim).ir(4,ndim+1+j);
	  // Galerkin and SUPG contribution
	  tmp9_g.prod(W_supg,SHAPE,1,2);
	  tmp10_g.prod(tmp9_g,dmatu_vp,1,3,2);
	  matlocf.axpy(tmp10_g,wpgdet*(rho_g-rho_l));
	  
	  // PSPG contribution --> dR_cont(PSPG)/d_alpha_g
	  tmp12_g.prod(P_pspg,dmatu_vp,-1,1,-1);
	  matlocf.rs().ir(2,ndim+1).ir(4,ndim+1+j);
	  matlocf.axpy(tmp12_g,-wpgdet*(rho_g-rho_l));
	  
	}
      }
      
      resmom.rs();
      dshapex.rs();
      P_pspg.rs();
      grad_vf_star.rs();
      matlocf.rs();
      gNi_Nj.rs();
 
	*/
     
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      
      // diffusive part
      matlocmom.axpy(gNi_gNj,visco_m_eff);
      
      dmatw.set(massm).scale(rho_l);
      
      if (update_jacobian) {
	for (int iloc=1; iloc<=nel; iloc++) {
	  for (int jloc=1; jloc<=nel; jloc++) {
	      double c = wpgdet*matlocmom.get(iloc,jloc);
	      for (int ii=1; ii<=ndim; ii++) {
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
	  tmp17.prod(dshapex,SHAPE,3,1,2).scale(-wpgdet);
	  matlocf.rest(tmp17).rs();
	  // d_R(cont)/d_p
	  matlocf.ir(2,ndim+1).ir(4,ndim+1).axpy(tmp13,-wpgdet).rs();

	  matlocf.rs();

	  if (!cache_grad_div_u) {
	    tmp19.set(dshapex).scale(visco_m_eff*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,3,4,1);
	    matlocf.is(2,1,ndim).is(4,1,ndim).add(tmp18).rs();

	    tmp19.set(dshapex).scale(delta_supg*rho_l*wpgdet);
	    tmp18.prod(dshapex,tmp19,2,1,4,3);
	    matlocf.is(2,1,ndim).is(4,1,ndim).add(tmp18).rs();

	  } else {
	    grad_div_u_coef += (delta_supg*rho_l+visco_m_eff)*wpgdet;
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
