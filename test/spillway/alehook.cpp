//__INSERT_LICENSE__
//$Id: alehook.cpp,v 1.6 2003/03/25 21:23:22 mstorti Exp $
#define _GNU_SOURCE

#include <cstdio>
#include <cassert>

#include <map>

#include <src/vecmacros.h>
#include <src/texthash.h>
#include <src/texthf.h>
#include <src/fem.h>
#include <src/util3.h>
#include <src/hook.h>
#include <src/dlhook.h>
#include <src/autostr.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/ampli.h>
#include "../plate/fifo.h"

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class ale_hook {
private:
  Mesh *mesh;
  int ndim;
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

void ale_hook::init(Mesh &mesh_a,Dofmap &dofmap,
	  TextHashTableFilter *options,const char *name) {
  int ierr;
  mesh = &mesh_a;
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,ndim,0);
  assert(ndim>0);
}

void ale_hook::time_step_pre(double time,int step) {}

void ale_hook::time_step_post(double time,int step,
			      const vector<double> &gather_values) {

  int ierr;
  PetscPrintf(PETSC_COMM_WORLD,"ALEHOOK: begins time_step_post() ....\n");
  int nnod = mesh->nodedata->nnod;
  int nu = mesh->nodedata->nu;
  AutoString command;
  command.sprintf("/usr/bin/make petscfem_step=%d spillway_mesh",step);
  int stat = system(command.str());
  command.clear();
  if (!MY_RANK) {
    // char * const argv[] = {"/usr/bin/make","spillway_mesh",NULL};
    // int stat = execv(argv[0],argv);
    if (stat==-1) {
      printf("ALEHOOK: Couldn't launch octave process...\n");
      abort();
    }
    FILE *fid = fopen("spillway.nod.tmp","r");
    // double *nodedata = mesh->nodedata->nodedata;
#define NODEDATA(j,k) VEC2(mesh->nodedata->nodedata,j,k,nu)
    for (int j=0; j<nnod; j++) {
      for (int k=0; k<ndim; k++) {
	int nread = fscanf(fid,"%lf",&NODEDATA(j,k));
	assert(nread==1);
      }
      double dummy;
      for (int k=ndim; k<nu; k++) {
	int nread = fscanf(fid,"%f",&dummy);
	assert(nread==1);
      }
    }
    fclose(fid);
  }
  ierr = MPI_Bcast(mesh->nodedata->nodedata, nnod*nu, MPI_DOUBLE, 0,PETSC_COMM_WORLD);
  assert(!ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"ALEHOOK: ends time_step_post() ....\n");
}

void ale_hook::close() {}

DL_GENERIC_HOOK(ale_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class ale_hook2 {
private:
  FILE *ns2mmv,*mmv2ns;
  Mesh *mesh;
  // The reference coordinates
  dvector<double> xnod0;
  int nnod,ndim,nu;
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

void ale_hook2::init(Mesh &mesh_a,Dofmap &dofmap,
	  TextHashTableFilter *options,const char *name) { 
  if (!MY_RANK) {
    if (1) {
      printf("ALE_HOOK2_INIT: Starting ALE_HOOK...\n");
      pid_t pid = fork();
      if (pid==-1) {
	printf("ALE_HOOK2_INIT: Couldn't fork MESH_MOVE process...\n");
	abort();
      }
      if (pid==0) {
	char * const argv[] = {"/usr/bin/make","mesh_move",NULL};
	int stat = execv(argv[0],argv);
	if (stat==-1) {
	  printf("ALE_HOOK2_INIT: Couldn't \"execv\" MESH_MOVE  process...\n");
	  abort();
	}
      } else printf("ALE_HOOK2_INIT: MESH_MOVE pid is %d\n",pid);
      printf("ALE_HOOK2_INIT: Done.\n");
    }

    printf("ALE_HOOK2_INIT: Opening fifos for communicating with MESH_MOVE.\n");

    ns2mmv = fopen("ns2mmv.fifo","w");
    assert(ns2mmv);
    setvbuf(ns2mmv,NULL,_IOLBF,0);

    mmv2ns = fopen("mmv2ns.fifo","r");
    assert(mmv2ns);

    printf("ALE_HOOK2: Done.\n");
    
  }
  mesh = &mesh_a;
  nnod = mesh->nodedata->nnod;
  ndim = mesh->nodedata->ndim;
  nu = mesh->nodedata->nu;
  if (!MY_RANK) {
    xnod0.a_resize(2,nnod,ndim);
    
    for (int k=0; k<nnod; k++) 
      for (int j=0; j<ndim; j++) 
	xnod0.e(k,j) = mesh->nodedata->nodedata[k*nu+j];
  }
}

void ale_hook2::time_step_pre(double time,int step) {}

void ale_hook2::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  // Displacements are read in server and sent to slaves
  if (!MY_RANK) {
    int ierr;
    fprintf(ns2mmv,"step %d\n",step);
    // Here goes reading the data from  mesh_move and assigning
    // to the new mesh
    int mmv_step = int(read_doubles(mmv2ns,"mmv_step_ok"));
    assert(step==mmv_step);
    
    // Reads displacements computed by `mesh_move' and add to nodedata
    FILE *fid = fopen("spillway_mmv.state.tmp","r");
    int nnod = mesh->nodedata->nnod;
    int ndim = mesh->nodedata->ndim;
    int nu = mesh->nodedata->nu;
    double d;
    double *nodedata = mesh->nodedata->nodedata;
    for (int k=0; k<nnod; k++) {
      for (int j=0; j<ndim; j++) {
	fscanf(fid,"%lf",&d);
	*nodedata++ = xnod0.e(k,j) + d;
      }
      for (int j=ndim; j<nu; j++) fscanf(fid,"%lf",&d);
    }
  }
  int ierr = MPI_Bcast(mesh->nodedata->nodedata, nnod*nu, MPI_DOUBLE, 0,PETSC_COMM_WORLD);
  assert(!ierr);
}

void ale_hook2::close() {
  if (!MY_RANK) {
    xnod0.clear();
    fprintf(ns2mmv,"step %d\n",-1);
  }
}

DL_GENERIC_HOOK(ale_hook2);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// spines:= directions (normalized) along which the displacements
// are imposed
dvector<double> spines;
// displ:= current displacements (they cumulate during time steps)
dvector<double> displ;
// fs2indx_t:= type of fs2indx
typedef map<int,int> fs2indx_t;
// fs2indx:= map (FS -> index FS list)
fs2indx_t fs2indx;

/** Hook executed in the mesh move run */ 
class ale_mmv_hook {
private:
  /// Fifos for cummunicating with the NS run
  FILE *ns2mmv,*mmv2ns;
  /// Number of nodes, dimension, number of nodes on the FS
  int nnod, ndim, nfs;
  // time step, relaxation coefficient
  double Dt, fs_relax;
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ale_mmv_hook::init(Mesh &mesh,Dofmap &dofmap,
	  TextHashTableFilter *options,const char *name) { 

  if (!MY_RANK) {
    printf("MESH_MOVE: Opening fifos for communicating with ALE_HOOK2.\n");
    
    ns2mmv = fopen("ns2mmv.fifo","r");
    assert(ns2mmv);
    
    mmv2ns = fopen("mmv2ns.fifo","w");
    assert(mmv2ns);
    setvbuf(mmv2ns,NULL,_IOLBF,0);
    
    printf("ALE_HOOK2: Done.\n");
  }
  int ierr = MPI_Barrier(PETSC_COMM_WORLD);
  assert(!ierr);

  nnod = mesh.nodedata->nnod;
  ndim = 2;
  // Read list of nodes on the FS ad spines (all of size nfs).
  FILE *fid = fopen("spillway.nod_fs.tmp","r");
  FILE *fid2 = fopen("spillway.spines.tmp","r");
  int indx=0;
  while(1) {
    // printf("indx %d\n",indx);
    int node;
    int nread = fscanf(fid,"%d",&node);
    if (nread==EOF) break;
    double n;
    for (int j=0; j<ndim; j++) {
      nread = fscanf(fid2,"%lf",&n);
      // printf("read %f\n",n);
      assert(nread==1);
      spines.push(n);
    }
    assert(fs2indx.find(node)==fs2indx.end());
    // load map
    fs2indx[node] = indx++;
  }
  fclose(fid);
  fclose(fid2);
  // Number of nodes on the FS
  nfs = fs2indx.size();
  assert(spines.size()==ndim*nfs);
  spines.reshape(2,nfs,ndim);
  
  displ.set_chunk_size(ndim*nfs);
  displ.a_resize(2,nfs,ndim);
  displ.set(0.);
  
  //o Time step.
  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,Dt,0.);
  assert(Dt>0.);
  //o Relaxation factor for the update of the free surface position. 
  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,fs_relax,1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ale_mmv_hook::time_step_pre(double time,int step) {
  int ierr;
  // verify step sent by NS run
  int step_sent = int(read_doubles(ns2mmv,"step"));
  if (step_sent==-1) {
    PetscPrintf(PETSC_COMM_WORLD, 
		"MESH_MOVE: step==-1 received, stopping myself.\n");
    PetscFinalize();
    exit(0);
  }
  assert(step=step_sent);
  // Open state file. Will read velocities on the FS
  // and update displ += v * Dt
  FILE *fid = fopen("spillway.state.tmp","r");
  // string buffer 
  AutoString line;
  vector<string> tokens;
  fs2indx_t::iterator qe = fs2indx.end();
  // Reads the whole file
  for (int node=1; node<=nnod; node++) {
    // Reads line (even if it is not in the FS)
    line.getline(fid);
    // if node is not on the FS, then skip
    fs2indx_t::iterator q = fs2indx.find(node);
    if (q==qe) continue;
    // index in the containers
    int indx = q->second;
    double v, vn_tmp=0., n2=0.;
    // Compute normal component of velocity
    // printf("node %d, vel, nor ",node);
    tokens.clear();
    tokenize(line.str(),tokens);
    assert(tokens.size()==ndim+1);
    for (int j=0; j<ndim; j++) {
      string2dbl(tokens[j],v);
      double n = spines.e(indx,j);
      // printf(" %f %f",v,n);
      vn_tmp += v*n;
      n2 += n*n;
    }
    // printf("\n");
    double tol = 1e-5;
    assert(fabs(n2-1.0)<tol);
    for (int j=0; j<ndim; j++) displ.e(indx,j) += fs_relax*Dt*vn_tmp*spines.e(indx,j);
  }
  fclose(fid);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ale_mmv_hook::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  fprintf(mmv2ns,"mmv_step_ok %d\n",step);
}

void ale_mmv_hook::close() {}

DL_GENERIC_HOOK(ale_mmv_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class fs_coupling : public DLGenericTmpl {
public:
  fs_coupling() { }
  void init(TextHashTable *thash) { }
  double eval(double);
  ~fs_coupling() { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double fs_coupling::eval(double) { 
  int ndim=2;
  int f = field();
  assert(f<=ndim);
  int node_c = node();
  fs2indx_t::iterator q = fs2indx.find(node_c);
  PETSCFEM_ASSERT(q != fs2indx.end(),
		  "Can't find node in FS node list\n"
		  "node %d\n",node_c);
  int indx = q->second;
  double val = displ.e(indx,f-1);
  // printf("fs_coupling: node %d, field %d, indx %d -> %f\n",node_c,f,indx,val);
  return val;
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(fs_coupling);
