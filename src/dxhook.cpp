//__INSERT_LICENSE__
//$Id: dxhook.cpp,v 1.8 2003/02/07 19:27:26 mstorti Exp $
#ifdef USE_SSL

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>

#include <HDR/sockets.h>

extern int MY_RANK, SIZE;

// I didn't found a definition for this
#ifndef IPPORT_MAX
#define IPPORT_MAX 65536
#endif

// It seems that DX hangs when using doubles for coordinates and
// data in the moment of making `DXEndField()' (it seems that internally
// it happens in th moment of doing `DXBoundingBox()'. So that
// I will use only floats. 
#define DX_USE_FLOATS

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::init(Mesh &mesh_a,Dofmap &dofmap_a,
			   const char *name_a) {
  int ierr;
  char skthost[10];
  if (!MY_RANK) {
    //o TCP/IP port for communicating with DX (5000 < dx_port < 65536). 
    TGETOPTDEF(mesh_a.global_options,int,dx_port,5314);
    PETSCFEM_ASSERT(dx_port>IPPORT_USERRESERVED && 
		    dx_port<IPPORT_MAX,"\"dx_port\" number must be in the range\n"
		    "IPPORT_USERRESERVED < dx_port < IPPORT_MAX, dx_port = %d\n"
		    "[current values are: IPPORT_USERRESERVED = %d\n, IPPORT_MAX  = %d]\n",
		    dx_port, IPPORT_USERRESERVED, IPPORT_MAX);
    printf("dx_hook: starting socket at port: %d\n",dx_port);
    sprintf(skthost,"S%d",dx_port);
    srvr_root = Sopen("",skthost);
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
    Sprintf(srvr,"nodes nodes %d %d\n",ndim,nnod);
    for (int node=0; node<nnod; node++) {
#ifndef DX_USE_FLOATS
      Swrite(srvr,xnod+node*nu,ndim*sizeof(double));
#else
      for (int j=0; j<ndim; j++) {
	float val = (float)*(xnod+node*nu+j);
	Swrite(srvr,&val,sizeof(float));
      }
#endif
    }
  }
  // Send results
  int ndof = dofmap->ndof;
  
  double *state_p = NULL;
  if (!MY_RANK) state_p = new double[ndof*nnod];
  int ierr = state2fields(state_p,state(),dofmap,time_data()); assert(!ierr);
  if (!MY_RANK) {
    Sprintf(srvr,"state state %d %d\n",ndof,nnod);
#ifndef DX_USE_FLOATS
    Swrite(srvr,state_p,ndof*nnod*sizeof(double));
#else
    for (int j=0; j<ndof*nnod; j++) {
      float val = (float)*(state_p+j);
      Swrite(srvr,&val,sizeof(float));
    }
#endif
  }
  delete[] state_p;

  // Send connectivities for each elemset
  Darray *elist = mesh->elemsetlist;
  for (int j=0; j<da_length(elist); j++) {
    Elemset *e = *(Elemset **)da_ref(elist,j);
    e->dx(srvr,nodedata,state_p);
  }

  // Send termination signal
  Sprintf(srvr,"end\n");

  Sclose(srvr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::close() {
  Sclose(srvr_root);
}

#endif
