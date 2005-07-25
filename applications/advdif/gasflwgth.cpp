//__INSERT_LICENSE__
//$Id: gasflwgth.cpp,v 1.2 2005/07/25 03:11:39 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/gatherer.h>

#include "./gasflwgth.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
gasflow_force_integrator
::init() {
  int ierr;
  //o Dimension of the embedding space
  TGETOPTNDEF(thash,int,ndim,none);
  ndim_m=ndim;
  //o Dimenson of the element
  TGETOPTDEF(thash,int,ndimel,ndim-1); 
  assert(ndimel==ndim-1);
  assert(gather_length==ndim);
  F.resize(1,ndim);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
gasflow_force_integrator
::set_pg_values(vector<double> &pg_values,FastMat2 &u,
		FastMat2 &uold,FastMat2 &xpg,FastMat2 &n,
		double wpgdet,double time) {
  double p = u.get(ndim_m+2); 
  F.set(n).scale(-p*wpgdet)
    .export_vals(&*pg_values.begin());
}
