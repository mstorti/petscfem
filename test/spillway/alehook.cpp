//__INSERT_LICENSE__
//$Id: alehook.cpp,v 1.3 2003/03/23 17:29:31 mstorti Exp $
#define _GNU_SOURCE

#include <cstdio>
#include <cassert>

#include <map>

#include <src/vecmacros.h>
#include <src/texthash.h>
#include <src/texthf.h>
#include <src/fem.h>
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
    if (0) {
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
}

void ale_hook2::time_step_pre(double time,int step) {}

void ale_hook2::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  if (!MY_RANK) {
    int ierr;
    fprintf(ns2mmv,"step %d\n",step);
  }
}

void ale_hook2::close() {}

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
  // time step
  double Dt;
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

    printf("MESH_MOVE: Opening fifos for communicating with ALE_HOOK2.\n");

    ns2mmv = fopen("ns2mmv.fifo","r");
    assert(ns2mmv);

    mmv2ns = fopen("mmv2ns.fifo","w");
    assert(mmv2ns);
    setvbuf(mmv2ns,NULL,_IOLBF,0);

    printf("ALE_HOOK2: Done.\n");

    nnod = mesh.nodedata->nnod;
    ndim = 2;
    // Read list of nodes on the FS ad spines (all of size nfs).
    FILE *fid = fopen("spillway.nod_fs.tmp","r");
    FILE *fid2 = fopen("spillway.spines.tmp","r");
    int indx=0;
    while(1) {
      printf("indx %d\n",indx);
      int node;
      int nread = fscanf(fid,"%d",&node);
      if (nread==EOF) break;
      double n;
      for (int j=0; j<ndim; j++) {
	nread = fscanf(fid2,"%lf",&n);
	printf("read %f\n",n);
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
    
    int ierr;
    TGETOPTDEF_ND(GLOBAL_OPTIONS,double,Dt,0.);
    assert(Dt>0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ale_mmv_hook::time_step_pre(double time,int step) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ale_mmv_hook::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  int ierr;
  // verify step sent by NS run
  int step_sent = int(read_doubles(ns2mmv,"step"));
  assert(step=step_sent);
  // Open state file. Will read velocities on the FS
  // and update displ += v * Dt
  FILE *fid = fopen("spillway.state.tmp","r");
  // string buffer 
  AutoString line;
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
    for (int j=0; j<ndim; j++) {
      int nread = sscanf(line.str(),"%lf",&v);
      double n = spines.e(indx,j);
      vn_tmp += v*n;
      n2 += n*n;
    }
    double tol = 1e-5;
    assert(fabs(n2-1.0)<tol);
    for (int j=0; j<ndim; j++) displ.e(indx,j) += Dt*spines.e(indx,j);
  }
  fclose(fid);
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
  return displ.e(indx,f-1);
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(fs_coupling);
