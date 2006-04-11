//__INSERT_LICENSE__
//$Id: rhhook.cpp,v 1.1 2006/04/11 12:36:20 mstorti Exp $

#include <src/debug.h>
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/util3.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/rhhook.h>
#include <src/dvector.h>
#include <src/dvector2.h>

extern int MY_RANK, SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::
init(Mesh &mesh_a,Dofmap &dofmap_a,
     const char *name_a) {
#if 0
  nu = mesh_a.nodedata->nu;
  assert(nu==2*ndim);

  // Store original coords for reference
  asdp->coords0.a_resize(2,asdp->nnod,ndim);
  for (int j=0; j<asdp->nnod; j++)
    for (int k=0; k<ndim; k++) 
      asdp->coords0.e(j,k)=COORDS(j,k);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::
time_step_pre(double time,int step) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::close() {
}
