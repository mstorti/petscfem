//__INSERT_LICENSE__
//$Id: gatherer.cpp,v 1.17 2003/01/25 15:28:58 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/gatherer.h>

#include "./gatherer.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void force_integrator::init() {
  int ierr;
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  ndim_m = ndim;
  //o Dimenson of the element
  TGETOPTNDEF(thash,int,ndimel,ndim-1); 
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void force_integrator::set_pg_values(vector<double> &pg_values,FastMat2 &u,
				     FastMat2 &uold,FastMat2 &xpg,FastMat2 &n,
				     double wpgdet,double time) {
  // Force contribution = normal * pressure * weight of Gauss point
  force.set(n).scale(-wpgdet*u.get(4));
  // export forces to return vector
  force.export_vals(pg_values.begin());
  if (compute_moment) {
    // Position offset of local point to center of moments
    dx.set(xpg).rest(x_center);
    // Moment contribution = force X dx
    moment.cross(force,dx);
    // export forces to return vector
    moment.export_vals(pg_values.begin()+ndim_m);
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
  TGETOPTNDEF(thash,int,ndim,none);
  //o If Navier-Stokes compressible
  TGETOPTNDEF(thash,int,nsc,1);
  ndim_m=ndim;
  //o Dimenson of the element
  TGETOPTNDEF(thash,int,ndimel,ndim-1); 
  assert(ndimel==ndim-1);
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
  Q.export_vals(pg_values.begin());
}
