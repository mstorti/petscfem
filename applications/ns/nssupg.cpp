//__INSERT_LICENSE__
//$Id: nssupg.cpp,v 1.1 2002/07/30 16:13:27 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>

extern TextHashTable *GLOBAL_OPTIONS;

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Linearized free surface boundary condition, but with a
// disipative term and Galerkin (not collocation as in `nssup.cpp'). 
// The free surface equation is $-\dot\eta+\cdump\Delta\eta+w = 0$
// The $\cdump$ gives some dumping avoiding spurious waves. 

#undef __FUNC__
#define __FUNC__ "ns_sup_g::assemble"
int ns_sup_g::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		       Dofmap *dofmap,const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_) {

  int ierr;
  //o Add LES for this particular elemset.
  GGETOPTDEF(int,LES,0);
  assert(!LES);
  //o This factor scales the temporal derivative term
  // int the free surface equation $w = f \dot p$ equation,
  // so that for $f=0$ we recover the slip boundary
  // condition $w=0$. 
  TGETOPTDEF(thash,double,fs_eq_factor,1.);
  assert(nel % 2==0); // one layer of nodes is for \eta the other for w
  nel2 = nel/2;

  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_mat_res);
  GET_JOBINFO_FLAG(comp_res);

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

  TGETOPTNDEF(thash,int,ndim,none); //nd
  // Spatial coordinate along the gravity field. Normally this is
  // the #z# direction in 3D (normal_dir=3) and #y# in 2D (normal_dir=2). 
  int normal_dir = ndim; // hardwired to z in 3D, y in 2D
  int nen = nel*ndof;

  // Get arguments from arg_list
  double *locst,*locst2,*retval,*retvalmat;
  if (comp_mat) {
    retvalmat = arg_data_v[0].retval;
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
    hmin = (arg_data_v[ja++].vector_assoc)->begin();
    ja_hmin=ja;
    glob_param = (GlobParam *)(arg_data_v[ja++].user_data);
    rec_Dt = 1./glob_param->Dt;
    if (glob_param->steady) rec_Dt=0.;
  } 
  double &alpha = glob_param->alpha;

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    locstate(2,nel,ndof), locstate2(2,nel,ndof),
    veccontr(2,nel,ndof),
    eta(1,nel2), eta_new(1,nel2), 
    w(1,nel2), w_new(1,nel2), w_star(1,nel2), 
    matlocf(4,nel,ndof,nel,ndof),xloc(2,nel2,ndim);

  double w,w_new,w_star,eta,eta_new,eta_star,res;

  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

  veccontr.set(0.);
  matlocf.set(0.);

  if (comp_mat) {
    matloc_prof.set(1.);
  }

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    if(comp_mat) {

      // this is irrelevant, the profile is passed now
      // via the .profile entry
      matloc_prof.export_vals(&(RETVALMAT(ielh,0,0,0,0)));
      
    } else if (comp_mat_res || comp_res) {

      // Load local node coordinates in local vector
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
      for (kloc=0; kloc<nel2; kloc++) {
	node = ICONE(k,kloc);
	xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
      }
      xloc.rs();
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));

      // Surface elevations and normal velocities at t^{n+1}
      locstate.is(1,nel2+1,nel).ir(2,1);
      eta.set(locstate);
      locstate.rs();
      locstate.is(1,1,nel2).ir(2,normal_dir);
      w.set(locstate);
      locstate.rs();

      // Surface elevations and normal velocities at t^{n}
      locstate2.is(1,nel2+1,nel).ir(2,1);
      eta.set(locstate2);
      locstate2.rs();
      locstate2.is(1,1,nel2).ir(2,normal_dir);
      w.set(locstate2);
      locstate2.rs();


      locstate.is(1,nel2+1,nel).ir(2,1);
      eta_new = locstate.get(2,1);
      // Normal velocities at t^n, t^{n+1}, t^*
      w = locstate2.get(1,normal_dir);
      w_new = locstate.get(1,normal_dir);
      w_star = alpha * w_new + (1.-alpha) * w;
      // Residual of equation (eta^{n+1}-eta^n)/Dt - w = 0
      res = -(eta_new-eta)*rec_Dt*fs_eq_factor + w_star;
      veccontr.setel(res,2,1);
      // alpha's here are not clear. 
      matlocf.setel(rec_Dt/alpha*fs_eq_factor,2,1,2,1)
	.setel(-alpha,2,1,1,ndim);
      veccontr.export_vals(&(RETVAL(ielh,0,0)));
      if (update_jacobian) matlocf.export_vals(&(RETVALMAT(ielh,0,0,0,0)));

    }
  }

  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
