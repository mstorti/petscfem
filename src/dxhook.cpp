//__INSERT_LICENSE__
//$Id: dxhook.cpp,v 1.15 2003/02/09 14:50:57 mstorti Exp $
#ifdef USE_SSL

#include <src/debug.h>
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/util3.h>
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
  //o TCP/IP port for communicating with DX (5000 < dx_port < 65536). 
  TGETOPTDEF(mesh_a.global_options,int,dx_port,5314);
  TGETOPTDEF_ND(mesh_a.global_options,int,steps,1);
  if (!MY_RANK) {
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

#define PF_DBG(name,format)					\
PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] " 		\
                        #name "  " #format "\n",MY_RANK,name);	\
PetscSynchronizedFlush(PETSC_COMM_WORLD); 
#define PF_DBG_INT(name)  PF_DBG(name,%d) 
#define PF_DBG_DBL(name)  PF_DBG(name,%f) 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static void *wait_connection(void *arg) {
  dx_hook *hook = (dx_hook *)arg;
  return hook->wait_connection();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void *dx_hook::wait_connection() {
  srvr = Saccept(srvr_root);
  connection_state_master = connected;
  return NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
dx_hook::connection_state_t 
dx_hook::connection_state() {
  if (!MY_RANK) connection_state_m = connection_state_master;
  ierr = MPI_Bcast (&connection_state_m, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  return connection_state_m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::set_connection_state(connection_state_t s) {
  connection_state_m = connection_state_master = s;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::re_launch_connection() {
  set_connection_state(not_connected);
  if (!MY_RANK) {
    ierr = pthread_create(&thread,NULL,&::wait_connection,this);
    assert(!ierr);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
  // this is for CHECK_COOKIE
#define sock srvr
  int cookie, cookie2, dx_step;
  string state_file;
  if (steps && step_cntr--) return;
  if (!steps) {
    if(connection_state() == not_launched) {
      re_launch_connection();
      sleep(1);
    }
    if (connection_state()==not_connected) return;
    else if (connection_state()==connected) {
      void *retval;
      if (!MY_RANK) {
	ierr = pthread_join(thread,&retval);
	assert(!ierr);
      }
      set_connection_state(not_launched);
    } else assert(0);
  }
#define BUFSIZE 512
  static char *buf = (char *)malloc(BUFSIZE);
  static size_t Nbuf = BUFSIZE;

  vector<string> tokens;
  Nodedata *nodedata = mesh->nodedata;
  double *xnod = nodedata->nodedata;
  int ndim = nodedata->ndim;
  int nnod = nodedata->nnod;
  int nu = nodedata->nu;
  PetscPrintf(PETSC_COMM_WORLD,
	      "dx_hook: accepting connections, step %d\n",step);

  // Process DX options. 
  if (!MY_RANK) {
    srvr = Saccept(srvr_root);
    assert(srvr);

    Sgetline(&buf,&Nbuf,srvr);
    tokenize(buf,tokens);

    // Parse DX options
    int j=0;
    while (1) {
      if (j>=tokens.size()) break;
      if (tokens[j]=="steps") {
	int stepso;
	assert(!string2int(tokens[++j],stepso));
	if (stepso>=0) {
	  if (stepso!=steps) 
	    printf("dx_hook: changed \"steps\" %d -> %d from DX\n",
		   steps,stepso);
	  steps=stepso;
	}
      } else if (tokens[j]=="step") {
	assert(!string2int(tokens[++j],dx_step));
      } else if (tokens[j]=="state_file") {
	state_file = tokens[++j];
      } else {
	printf("Unknown option \"%s\"\n",tokens[j].c_str());
      }
      j++;
    }
    printf("dx_hook: Got steps %d, dx_step %d, state_file %s\n",
	   steps,dx_step,state_file.c_str());
  }
  // Options are read in master and
  // each option is sent to the slaves with MPI_Bcast
  ierr = MPI_Bcast (&steps, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  ierr = MPI_Bcast (&dx_step, 1, MPI_INT, 0,PETSC_COMM_WORLD);
  ierr = string_bcast(state_file,0,PETSC_COMM_WORLD);
  
  step_cntr = steps-1;
  
  if (!MY_RANK) {
    // Send node coordinates
    cookie = rand();
    Sprintf(srvr,"nodes nodes %d %d %d\n",ndim,nnod,cookie);
    // printf("sending nodes nodes %d %d %d\n",ndim,nnod,cookie);
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
    CHECK_COOKIE(nodes);
  }

  // Send results
  int ndof = dofmap->ndof;
  
  double *state_p = NULL;
  if (!MY_RANK) state_p = new double[ndof*nnod];
  ierr = state2fields(state_p,state(),dofmap,time_data()); assert(!ierr);
  if (!MY_RANK) {
    cookie = rand();
    Sprintf(srvr,"state state %d %d %d\n",ndof,nnod,cookie);
#ifndef DX_USE_FLOATS
    Swrite(srvr,state_p,ndof*nnod*sizeof(double));
#else
    for (int j=0; j<ndof*nnod; j++) {
      float val = (float)*(state_p+j);
      Swrite(srvr,&val,sizeof(float));
    }
#endif
    CHECK_COOKIE(state);
  }
  delete[] state_p;

  // Send connectivities for each elemset
  Darray *elist = mesh->elemsetlist;
  for (int j=0; j<da_length(elist); j++) {
    Elemset *e = *(Elemset **)da_ref(elist,j);
    e->dx(srvr,nodedata,state_p);
  }

  // Send termination signal
  if (!MY_RANK) {
    Sprintf(srvr,"end\n");
    Sclose(srvr);
  }
  if(connection_state() == not_launched) re_launch_connection();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dx_hook::close() {
  // I'm not too much sure how to do this correctly
  if (connection_state()==connected) {
    // Here we should shut down the connection
    // sending some message to the client
    set_connection_state(not_connected);
  }
  if (connection_state()==not_connected) {
    pthread_cancel(thread);
    set_connection_state(not_launched);
  }
  assert(connection_state()==not_launched);
  if (!MY_RANK) Sclose(srvr_root); 
}

#endif
