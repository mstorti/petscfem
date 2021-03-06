//__INSERT_LICENSE__
//$Id: advdfgth.cpp,v 1.4 2005/08/22 19:42:41 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/gatherer.h>

#include "./advdfgth.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void flow_rate_integrator::init() {
  int ierr;
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  ndim_m=ndim;
  //o Dimenson of the element
  TGETOPTDEF(thash,int,ndimel,ndim-1); 
  assert(ndimel==ndim-1);
  assert(gather_length==1);
  Q.resize(0);

  //o Give mass rate or volumetric rate
  TGETOPTDEF(thash,int,use_mass_rate,1); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void flow_rate_integrator::set_pg_values(vector<double> &pg_values,FastMat2 &u,
				     FastMat2 &uold,FastMat2 &xpg,FastMat2 &n,
				     double wpgdet,double time) {
  double rho = u.get(1); 
  u.is(1,2,ndim_m+1); 
  double w = (use_mass_rate ? rho*wpgdet : wpgdet);
  Q.prod(n,u,-1,-1).scale(w);
  u.rs();
  Q.export_vals(&*pg_values.begin());
}
