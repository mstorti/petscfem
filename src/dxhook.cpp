//__INSERT_LICENSE__
//$Id: dxhook.cpp,v 1.2 2003/02/04 13:32:01 mstorti Exp $
#ifdef USE_SSL

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>

#include <HDR/sockets.h>

#define PF_DX_PORT "5555"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::init(Mesh &mesh_a,Dofmap &dofmap,
			   const char *name_a) {
  PetscPrintf(PETSC_COMM_WORLD,
	      "dx_hook: starting socket at port: %s\n",PF_DX_PORT);
  srvr_root = Sopen("","s" PF_DX_PORT);
  assert(srvr_root);
  PetscPrintf(PETSC_COMM_WORLD,"Done.\n");
  mesh = &mesh_a;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::time_step_pre(double time,int step) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
  PetscPrintf(PETSC_COMM_WORLD,
	      "dx_hook: accepting connections\n");
  srvr = Saccept(srvr_root);
  assert(srvr);
  Nodedata *nodedata = mesh->nodedata;
  double *xnod = nodedata->nodedata;
  int ndim = nodedata->ndim;
  int nnod = nodedata->nnod;
  int nu = nodedata->nu;
  Sprintf(srvr,"nodes %d %d\n",ndim,nnod);
  for (int node=0; node<nnod; node++)
    Swrite(srvr,xnod+node*nu,ndim*sizeof(double));
  Sclose(srvr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::close() {
  Sclose(srvr_root);
}

#endif
