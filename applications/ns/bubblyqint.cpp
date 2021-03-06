//__INSERT_LICENSE__
//$Id: bubblyqint.cpp,v 1.3 2006/09/04 20:18:04 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/srfgath.h>
#include "./bubblyqint.h"

void bubbly_flow_rate_integrator::init() {
  int ierr;
  //o Number of dimensions
  TGETOPTDEF_ND(thash,int,ndim,-1);
  assert(ndim>0);
  //o Number of disperse phases
  TGETOPTDEF_ND(thash,int,nphases,1);

  //o Direction of gravity
  TGETOPTDEF_ND(thash,int,g_dir,ndim);
  // Ojo puede ser negativo
  assert(g_dir>0);

  assert(ndof>=ndim+1+nphases);

  //  slip velocity vector
  vslip_user_vp.resize(1,nphases);
  vslip_user_vp.set(0.0);
  ierr = get_double(GLOBAL_OPTIONS,"vslip_user_phases",
		    vslip_user_vp.storage_begin(),1,nphases);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int bubbly_flow_rate_integrator::vals_per_plane() { 
  // Integral quantities are:
  // 0: cut area
  // 1: bulk flow rate
  // k+2: flow rate for `k' phase (0-base)
  return nphases+2; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_flow_rate_integrator
::set_ip_values(vector<double> &ip_values,FastMat2 &u,
		FastMat2 &xpg,FastMat2 &n,double time) {
  ip_values[0] = 1.0;
  u.is(1,1,ndim);
  tmp.prod(u,n,-1,-1);
  u.rs();
  double un = double(tmp);
  ip_values[1] = un;
  double gn = n.get(g_dir);
  for (int j=0; j<nphases; j++) {
    double u_phase_n = un + vslip_user_vp.get(1+j)*gn;
    ip_values[j+2] = u_phase_n * u.get(ndim+2+j);
  }
}
