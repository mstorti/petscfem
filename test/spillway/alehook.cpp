//__INSERT_LICENSE__
//$Id: alehook.cpp,v 1.11 2003/03/31 17:07:10 mstorti Exp $
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
#include <src/fastmat2.h>
#include "../plate/fifo.h"

extern int MY_RANK,SIZE;
#define CASE_NAME "wave"
//#define CASE_NAME "spillway"

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
  command.sprintf("/usr/bin/make petscfem_step=%d " CASE_NAME "_mesh",step);
  int stat = system(command.str());
  command.clear();
  if (!MY_RANK) {
    if (stat==-1) {
      printf("ALEHOOK: Couldn't launch octave process...\n");
      abort();
    }
    FILE *fid = fopen(CASE_NAME ".nod.tmp","r");
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
    printf("ALE_HOOK2: starting init()\n");
    int ierr;
    //o Flag fo launching the `mesh_move' process.
    TGETOPTDEF(GLOBAL_OPTIONS,int,launch_mesh_move,1);
    if (launch_mesh_move) {
      printf("ALE_HOOK2_INIT: Starting ALE_HOOK...\n");
      pid_t pid = fork();
      if (pid==-1) {
	printf("ALE_HOOK2_INIT: Couldn't fork MESH_MOVE process...\n");
	abort();
      }
      if (pid==0) {
	char * const argv[] = {"/usr/bin/make",CASE_NAME "_mmv",NULL};
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
    printf("ALE_HOOK2: ending init()\n");
  }
}

void ale_hook2::time_step_pre(double time,int step) {}

void ale_hook2::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  // Displacements are read in server and sent to slaves
  if (!MY_RANK) {
    printf("ALE_HOOK2: starting time_step_post()\n");
    int ierr;
    fprintf(ns2mmv,"step %d\n",step);
    // Here goes reading the data from  mesh_move and assigning
    // to the new mesh
    int mmv_step = int(read_doubles(mmv2ns,"mmv_step_ok"));
    assert(step==mmv_step);
    
    // Reads displacements computed by `mesh_move' and add to nodedata
    FILE *fid = fopen(CASE_NAME "_mmv.state.tmp","r");
    double d;
    double *nodedata = mesh->nodedata->nodedata;
    int debug=0;
    for (int k=0; k<nnod; k++) {
      for (int j=0; j<ndim; j++) {
	fscanf(fid,"%lf",&d);
	*nodedata++ = xnod0.e(k,j) + d;
	if (debug) printf("x0 %f, d %f",xnod0.e(k,j),d);
      }
      if (debug) printf("\n");
    }
#if 0
    nodedata = mesh->nodedata->nodedata;
    for (int k=0; k<nnod; k++)
	printf("node %d, x,y: %f %f\n",k+1,nodedata[k*nu+0],nodedata[k*nu+1]);
#endif
    printf("ALE_HOOK2: ending time_step_post()\n");
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
// fs:= vector of nodes on the free surface (FS)
dvector<int> fs;
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
  // Pointer to the mesh, in order to obtain the nodedata
  Mesh *mesh;
  // Flags whether the FS is assumed to be cyclic 
  int cyclic_fs;
  // If the problem is cyclic then this is the period 
  double cyclic_length;
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ale_mmv_hook::init(Mesh &mesh_a,Dofmap &dofmap,
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

  mesh = &mesh_a;
  // fixme:= some things here will not run in parallel
  nnod = mesh->nodedata->nnod;
  ndim = 2;
  // Read list of nodes on the FS ad spines (all of size nfs).
  FILE *fid = fopen(CASE_NAME ".nod_fs.tmp","r");
  FILE *fid2 = fopen(CASE_NAME ".spines.tmp","r");
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
    fs.push(node);
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
  //o Assume problem is periodic 
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,cyclic_fs,0);
  //o Assume problem is periodic 
  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,cyclic_length,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ale_mmv_hook::time_step_pre(double time,int step) {
  PetscPrintf(PETSC_COMM_WORLD,"ALE_MMV_HOOK: step %d, time_step_pre() starts.\n",step);
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
  FILE *fid = fopen(CASE_NAME ".state.tmp","r");
  // string buffer 
  AutoString line;
  vector<string> tokens;
  fs2indx_t::iterator qe = fs2indx.end();
  // Get coordinates pointer
  double *nodedata = mesh->nodedata->nodedata;
  int nu = mesh->nodedata->nu;
  int nnod = mesh->nodedata->nnod;
  FastMat2 x0(1,ndim),x1(1,ndim),x2(1,ndim),
    x01(1,ndim),x12(1,ndim),normal(1,ndim);
  // Reads the whole file
  for (int node=1; node<=nnod; node++) {
    // Reads line (even if it is not in the FS)
    line.getline(fid);
    // if node is not on the FS, then skip
    fs2indx_t::iterator q = fs2indx.find(node);
    if (q==qe) continue;
    // index in the containers
    int indx = q->second;
    // Compute normal component of velocity
    int debug = 0;
    if (debug) printf("node %d, vel, nor, spin ",node);
    tokens.clear();
    tokenize(line.str(),tokens);
    assert(tokens.size()==ndim+1);
    // -- Compute normal --
    // FS nodes at the end of the FS must be treated specially.
    // Currently assumes cyclic FS
    assert(cyclic_fs);
    // Compute coordinates of three nodes on the free surface
    // Node at the center
    x1.set(&nodedata[nu*(node-1)]);
    // Previous node
    int div0,div2;
    int indx0 = modulo(indx-1,nfs,&div0);
    int node0 = fs.e(indx0);
    x0.set(&nodedata[nu*(node0-1)]);
    x0.addel(div0*cyclic_length,1);

    // Next node
    int indx2 = modulo(indx+1,nfs,&div2);
    int node2 = fs.e(indx2);
    x2.set(&nodedata[nu*(node2-1)]);
    x2.addel(div2*cyclic_length,1);

#if 1
    printf("indx %d, div0 %d, div2 %d\n",indx,div0,div2);
    x0.print("x0:");
    x1.print("x1:");
    x2.print("x2:");
#endif

    // Segements 0-1 1-2
    x01.set(x1).rest(x0);
    x12.set(x2).rest(x1);
    double l01 = x01.norm_p_all(2);
    double l12 = x12.norm_p_all(2);

    // weighted tangent at 1
    normal.set(x01).scale(l01).axpy(x12,l12);
    double l = normal.norm_p_all(2);
    normal.scale(1./l);

    // rotate to get normal
    assert(ndim==2);		// not implemented yet 3D
    double nx = normal.get(2);
    double ny = -normal.get(1);
    normal.setel(nx,1);
    normal.setel(ny,2);

    // normalize
    double n2 = normal.norm_p_all(2);
    normal.scale(1./n2);

    // Compute displacement
    double v, vn=0., sn=0.;
    for (int j=0; j<ndim; j++) {
      string2dbl(tokens[j],v);
      double s = spines.e(indx,j);
      double n = normal.get(j+1);
      if (debug) printf(" %f %f %f",v,n,s);
      vn += v*n;
      sn += s*n;
    }
    if (debug) printf(", vn, sn: %f\n",vn,sn);
    double d = fs_relax*Dt*vn/sn;
    for (int j=0; j<ndim; j++) {
      displ.e(indx,j) += d*spines.e(indx,j);
    }
    if (debug) printf("%f %f\n",displ.e(indx,0),displ.e(indx,1));
  }
  fclose(fid);
  PetscPrintf(PETSC_COMM_WORLD,"ALE_MMV_HOOK: time_step_pre() ends.\n");
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
