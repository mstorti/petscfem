//__INSERT_LICENSE__
//$Id: nssupg.cpp,v 1.9 2002/09/30 21:45:32 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>

extern TextHashTable *GLOBAL_OPTIONS;

#ifdef ROSI_COUPLING_MODULE
double AVERAGE_ELEVATION=0.;
#else
const double AVERAGE_ELEVATION=0.;
#endif

#define MAXPROP 100
extern int fractional_step, reuse_mat;

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

  assert(!fractional_step);	// Not implemented with `fractional step'
  int ierr;
  //o Add LES for this particular elemset.
  GGETOPTDEF(int,LES,0);
  assert(!LES);
  //o $\Cnst{eq}$=\alltt{fs_eq_factor} (see doc for {\tt
  // free\_surface\_damp} option) is a factor that scales the free
  // surface ``rigidity''. $\Cnst{eq}=1$ (which is the default) means
  // no scaling, a zero value means infinitely rigid (as for an
  // inifinite gravity).
  TGETOPTDEF(thash,double,fs_eq_factor,1.);
  //o $\Cnst{lf}=$\altt{free\_surface\_set\_level\_factor} tries to
  // keep the free surface level constant by adding a term $\propto
  // \bar\eta$ to the free surface level.  (see doc //for {\tt
  // free\_surface\_damp}})
  TGETOPTDEF(thash,double,free_surface_damp,0.);
  //o This adds a $\Cnst{lf}\eta$ term in the free surface equation
  // in order to have the total meniscus volume constant. 
  TGETOPTDEF(thash,double,free_surface_set_level_factor,0.);
#define C_FSL free_surface_set_level_factor
#define C_EQ fs_eq_factor
#define C_DAMP free_surface_damp

  assert(nel % 2==0); // one layer of nodes is for \eta the other for w
  int nel2 = nel/2;

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

  int ndimel = ndim-1;
  // Unpack nodedata
  int nu=nodedata->nu;
  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel2,npg,GP_FASTMAT2);

  FastMat2 matloc_prof(4,nel,ndof,nel,ndof),
    locstate(2,nel,ndof), locstate2(2,nel,ndof),
    veccontr(2,nel,ndof),
    eta(1,nel2), eta_new(1,nel2), 
    w(1,nel2), w_new(1,nel2), w_star(1,nel2), 
    res(1,nel2), matlocf(4,nel,ndof,nel,ndof),xloc(2,nel2,ndim),
    res_pg, res_pgg, n(1,ndim), tmp, tmp1, grad_eta, tmp2, 
    mass_mat(2,nel2,nel2), lap_mat(2,nel2,nel2),
    tmp3, tmp4, eta_star(1,nel2), Jaco, iJaco, dshapex(2,ndimel,nel2);

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
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
      for (int kloc=0; kloc<nel2; kloc++) {
	int node = ICONE(k,kloc);
	xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
      }
      xloc.rs();
      locstate.set(&(LOCST(ielh,0,0)));
      locstate2.set(&(LOCST2(ielh,0,0)));

      // Surface elevations and normal velocities at t^{n+1}
      // `_new' values
      locstate.is(1,nel2+1,nel).ir(2,1);
      eta_new.set(locstate);
      locstate.rs();
      locstate.is(1,1,nel2).ir(2,normal_dir);
      w_new.set(locstate);
      locstate.rs();

      // Surface elevations and normal velocities at t^{n}
      locstate2.is(1,nel2+1,nel).ir(2,1);
      eta.set(locstate2);
      locstate2.rs();
      locstate2.is(1,1,nel2).ir(2,normal_dir);
      w.set(locstate2);
      locstate2.rs();

      res.set(0.);
      mass_mat.set(0.);
      lap_mat.set(0.);
      matlocf.set(0.);
      veccontr.set(0.);

      // Normal velocities at t^{*}
      w_star.set(w_new).scale(alpha).axpy(w,1-alpha);
      // elevation at t^*
      eta_star.set(eta_new).scale(alpha).axpy(eta,1-alpha);
      // Residual of free surface equation at t^* (nodal values).
      tmp.set(eta_new).rest(eta).scale(-rec_Dt*C_EQ)
	.add(-AVERAGE_ELEVATION*C_FSL*C_EQ)
	.add(w_star);

#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

      veccontr.is(1,nel2+1,nel).ir(2,1);
      // loop over Gauss points
      for (int ipg=0; ipg<npg; ipg++) {
	
	assert(normal_dir==3);
	xloc.is(2,1,2);
	Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
	xloc.rs();
	double detJaco = Jaco.det();
	iJaco.inv(Jaco);
	double wpgdet = detJaco*WPG;

	dshapex.prod(iJaco,DSHAPEXI,1,-1,-1,2);
	// Residual of equation (eta^{n+1}-eta^n)/Dt - w = 0
	res_pg.prod(tmp,SHAPE,-1,-1);
	res_pgg.set(SHAPE).scale(res_pg.get()*wpgdet);

#if 1
	grad_eta.prod(dshapex,eta_star,1,-1,-1);
	tmp2.prod(dshapex,grad_eta,-1,1,-1);
	res_pgg.axpy(tmp2,-C_DAMP*wpgdet);
#endif
	veccontr.add(res_pgg);

	// Jacobian term
	tmp3.prod(SHAPE,SHAPE,1,2);
	mass_mat.axpy(tmp3,wpgdet);
#if 1
	tmp4.prod(dshapex,dshapex,-1,1,-1,2);
	lap_mat.axpy(tmp4,wpgdet);
#endif
      }

      veccontr.rs().export_vals(&(RETVAL(ielh,0,0)));

      matlocf.is(1,nel2+1,nel).ir(2,1).is(3,nel2+1,nel).ir(4,1)
	.set(mass_mat)
	.scale(rec_Dt/alpha*C_EQ) // temporal term
	.axpy(lap_mat,C_DAMP) // surface diffusion term
	.rs();
      // fixme:= perhaps alpha doesn't go here... 
      matlocf.is(1,nel2+1,nel).ir(2,1).is(3,1,nel2).ir(4,normal_dir)
	.set(mass_mat).scale(-alpha).rs();
      
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
