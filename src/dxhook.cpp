//__INSERT_LICENSE__
//$Id: dxhook.cpp,v 1.4 2003/02/05 02:38:11 mstorti Exp $
#ifdef USE_SSL

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>

#include <HDR/sockets.h>

#define PF_DX_PORT "5555"

extern int MY_RANK, SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::init(Mesh &mesh_a,Dofmap &dofmap_a,
			   const char *name_a) {
  if (!MY_RANK) {
    printf("dx_hook: starting socket at port: %s\n",PF_DX_PORT);
    srvr_root = Sopen("","s" PF_DX_PORT);
    assert(srvr_root);
    PetscPrintf(PETSC_COMM_WORLD,"Done.\n");
  }
  mesh = &mesh_a;
  dofmap = &dofmap_a;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::time_step_pre(double time,int step) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {

  Nodedata *nodedata = mesh->nodedata;
  double *xnod = nodedata->nodedata;
  int ndim = nodedata->ndim;
  int nnod = nodedata->nnod;
  int nu = nodedata->nu;
  PetscPrintf(PETSC_COMM_WORLD,
	      "dx_hook: accepting connections\n");
  if (!MY_RANK) {
    srvr = Saccept(srvr_root);
    assert(srvr);
    // Send node coordinates
    Sprintf(srvr,"nodes %d %d\n",ndim,nnod);
    for (int node=0; node<nnod; node++)
      Swrite(srvr,xnod+node*nu,ndim*sizeof(double));
  }
  // Send results
  int ndof = dofmap->ndof;
  
  double *fields = NULL;
  if (!MY_RANK) fields = new double[ndof*nnod];
  int ierr = state2fields(fields,state(),dofmap,time_data()); assert(!ierr);
  if (!MY_RANK) {
    Sprintf(srvr,"fields %d %d\n",ndof,nnod);
    Swrite(srvr,fields,ndof*nnod*sizeof(double));
  }
  delete[] fields;

  // Send connectivities for each elemset
  Darray *elist = mesh->elemsetlist;
  for (int j=0; j<da_length(elist); j++) {
    Elemset *e = *(Elemset **)da_ref(elist,j);
    if (!MY_RANK) {
      Sprintf(srvr,"icone %d %d %s %s\n",e->nelem,e->nel,e->name(),
	      e->type);
      printf("type %s\n",e->type);
      Swrite(srvr,e->icone,e->nelem*e->nel*sizeof(int));
    }
  }
  Sprintf(srvr,"end\n");

  Sclose(srvr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::close() {
  Sclose(srvr_root);
}

#endif
