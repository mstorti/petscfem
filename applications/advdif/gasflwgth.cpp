//__INSERT_LICENSE__
//$Id: gasflwgth.cpp,v 1.4 2005/07/26 18:41:30 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/gatherer.h>

#include "./gasflwgth.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void 
read_double_array2(const Elemset *e,
		  const char *key,
		  vector<double> &v) {
  v.clear();
  const char *line;
  e->thash->get_entry(key,line);
  if (line) read_double_array(v,line);
}

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
  if (ndim==2) {
    assert(gather_length==ndim || gather_length==3);
    comp_moments = (gather_length==3);
    dim_mom = 1;
  } else if (ndim==3) {
    assert(gather_length==ndim || gather_length==2*ndim);
    comp_moments = (gather_length==2*ndim);
    dim_mom = ndim;
  } else { assert(0); }
  F.resize(1,ndim);
  M.resize(1,dim_mom);
  dx.resize(1,ndim);

  // Reference position to take moments
  x0.resize(1,ndim).set(0.);

  vector<double> x0v;
  //o _T: vector<int>  _N: x0 _D: 0 x ndim vector
  // _DOC: Position from where to take moments. 
  read_double_array2(this,"x0",x0v);
  if (x0v.size()) {
    assert(x0v.size()==ndim);
    x0.set(&x0v[0]);
  }
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

  if (comp_moments) {
    dx.set(xpg).rest(x0);
    M.cross(F,dx)
      .export_vals(&pg_values[ndim_m]);
  }

#if 0
  F.print("F: ");
  M.print("M: ");
#endif

}
