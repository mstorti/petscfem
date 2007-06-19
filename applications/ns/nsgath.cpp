//__INSERT_LICENSE__
//$Id mstorti-v6-1-11-g23c052c Mon Jun 18 13:08:24 2007 -0300$

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/gatherer.h>

#include "./nsgath.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void force_integrator::init() {
  int ierr;
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  ndim_m = ndim;

  //o Dimension of the element
  TGETOPTNDEF(thash,int,ndimel,ndim-1); 

  //o Viscosity
  TGETOPTDEF_ND(thash,double,viscosity,NAN); 

  //o Density
  TGETOPTDEF_ND(thash,double,rho,1.0); 

  //o Distance to the wall
  TGETOPTDEF_ND(thash,double,y_wall,NAN); 

  //o Add wall-law contribution
  TGETOPTDEF_ND(thash,int,add_wall_law_contrib,0); 

  nu = NAN;
  if (add_wall_law_contrib) {
    PETSCFEM_ASSERT0(!add_wall_law_contrib || 
                     (!isnan(y_wall) && !isnan(rho) && !isnan(viscosity)),
                     "add_wall_law_contrib was set, but not physical "
                     "values were given");
    PETSCFEM_ASSERT0(viscosity>=0,"viscosity must be non-negative"); 
    PETSCFEM_ASSERT0(rho>=0,"Density must be non-negative"); 
    PETSCFEM_ASSERT0(y_wall>=0,"Wall-law must be non-negative"); 
    nu = viscosity/rho;
  }

  assert(ndimel==ndim-1);
  assert(gather_length==ndim || gather_length==2*ndim);
  compute_moment = (gather_length==2*ndim);
  force.resize(1,ndim);
  moment.resize(1,ndim);
  x_center.resize(1,ndim).set(0.);
  dx.resize(1,ndim);
  //o _T: double[ndim] _N: moment_center _D: null vector 
  // _DOC: Center of moments. _END
  get_double(thash,"moment_center",x_center.storage_begin(),1,ndim);  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
force_integrator::wall_law_solver_t
::wall_law_solver_t() {
  yp1 = 5;
  yp2 = 30;
  clog = 2.5;
  c1 = -yp1*log(yp1)+yp1;
  c2 = yp1*log(yp2)+c1 -clog*log(yp2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
double force_integrator::wall_law_solver_t
::residual(double yplus,void *user_data) {
  double f;
  if (yplus<yp1) f = yplus;
  else if (yplus<yp2) f = yp1*log(yplus)+c1;
  else f = clog*log(yplus)+c2;
  return Re_wall - yplus*f;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
double force_integrator::wall_law_solver_t
::fdot(double yplus) {
  double fd;
  if (yplus<yp1) fd = 1;
  else if (yplus<yp2) fd = yp1/yplus;
  else fd = clog/yplus;
  return fd;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void force_integrator::set_pg_values(vector<double> &pg_values,FastMat2 &u,
				     FastMat2 &uold,FastMat2 &xpg,FastMat2 &n,
				     double wpgdet,double time) {
  // Force contribution = normal * pressure * weight of Gauss point
  force.set(n).scale(-wpgdet*u.get(ndim_m+1));
  if (add_wall_law_contrib) {
    // Add contribution following wall-law
    // FIXME:= we should check here that the wall law
    // used here is the same used in the wall-law element!!
    u.is(1,1,ndim_m);
    double uu = u.norm_p_all();
    // Now we have to invert the relation Re_wall = y+ f(y+)
    // for y+, where f() is the universal law of the wall
    wall_law_solver.Re_wall = uu*y_wall/nu;
    double yplus = wall_law_solver.sol();
    double coeff;
    if (yplus>wall_law_solver.yp1) {
      double ustar = nu*yplus/y_wall;
      double tau_wall = rho*ustar*ustar;
      coeff = tau_wall/uu;
    } else {
      coeff = viscosity/y_wall;
    }
    force.axpy(u,-coeff);
    u.rs();
  }
  // export forces to return vector
  force.export_vals(&*pg_values.begin());
  if (compute_moment) {
    // Position offset of local point to center of moments
    dx.set(xpg).rest(x_center);
    // Moment contribution = force X dx
    moment.cross(force,dx);
    // export forces to return vector
    moment.export_vals(&*pg_values.begin()+ndim_m);
#if 0
#define SHM(name) name.print(#name ": ")
    SHM(xpg);
    SHM(dx);
    SHM(force);
    SHM(moment);
#endif
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void force_integrator::clean() {
  force.clear();
  x_center.clear();
  dx.clear();
  moment.clear();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void flow_rate_integrator::init() {
  int ierr;
  //o Dimension of the embedding space
  TGETOPTDEF(thash,int,ndim,0);
  assert(ndim>0);
  // o If Navier-Stokes compressible
  // TGETOPTNDEF(thash,int,nsc,1);
  ndim_m=ndim;
  //o Dimenson of the element
  TGETOPTDEF(thash,int,ndimel,ndim-1); 
  assert(ndimel==ndim-1); // not implemented ndimel!=ndim-1
  assert(gather_length==1);
  Q.resize(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void flow_rate_integrator::set_pg_values(vector<double> &pg_values,FastMat2 &u,
				     FastMat2 &uold,FastMat2 &xpg,FastMat2 &n,
				     double wpgdet,double time) {
  u.is(1,1,ndim_m);
  Q.prod(n,u,-1,-1).scale(wpgdet);;
  u.rs();
  Q.export_vals(&*pg_values.begin());
}
