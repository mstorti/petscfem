//__INSERT_LICENSE__
//$Id: godunov_hook.cpp,v 1.1.2.2 2004/07/07 15:17:53 mstorti Exp $
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
#include "../../test/plate/fifo.h"

extern int MY_RANK,SIZE;
#define CASE_NAME "reflection"
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

const vector<double> *gather_values = NULL;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Hook executed in the advdif run */ 
class godunov_hook {
private:
  FILE *adv2god,*god2adv;
  int nnod,ndim,nu;
  Mesh *mesh;
  dvector<double> xnod0;
  void read_mesh();
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

void godunov_hook::init(Mesh &mesh_a,Dofmap &dofmap,
	  TextHashTableFilter *options,const char *name) { 
  if (!MY_RANK) {
    printf("ADV_HOOK: starting init()\n");
    int ierr;
    //o Flag fo launching the godunov process.
    TGETOPTDEF(GLOBAL_OPTIONS,int,launch_godunov,1);
    if (launch_godunov) {
      printf("ADV_HOOK_INIT: Starting ADV_HOOK...\n");
      pid_t pid = fork();
      if (pid==-1) {
	printf("ADV_HOOK_INIT: Couldn't fork GODUNOV process...\n");
	abort();
      }
      if (pid==0) {
	char * const argv[] = {"/usr/bin/make",CASE_NAME "_god",NULL};
	int stat = execv(argv[0],argv);
	if (stat==-1) {
	  printf("ADV_HOOK_INIT: Couldn't \"execv\" GODUNOV  process...\n");
	  abort();
	}
      } else printf("ADV_HOOK_INIT: GODUNOV pid is %d\n",pid);
      printf("ADV_HOOK_INIT: Done.\n");
    }
    
    printf("ADV_HOOK_INIT: Opening fifos for communicating with GODUNOV.\n");

    adv2god = fopen("adv2god.fifo","w");
    assert(adv2god);
    setvbuf(adv2god,NULL,_IOLBF,0);
    
    god2adv = fopen("god2adv.fifo","r");
    assert(god2adv);
    
    printf("ADV_HOOK: Done.\n");
    
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

    printf("ADV_HOOK: ending init()\n");
  }
  int ierr;
  //o Restart previous run
  TGETOPTDEF(GLOBAL_OPTIONS,int,restart,0);
  if (restart) read_mesh();
}

void godunov_hook::time_step_pre(double time,int step) {}

void godunov_hook::read_mesh() {
  // Reads states computed by godunov and add to nodedata
  FILE *fid = fopen(CASE_NAME "_god.state.tmp","r");
  double d;
  double *nodedata = mesh->nodedata->nodedata;
  for (int k=0; k<nnod; k++) {
    for (int j=0; j<ndim; j++) {
      fscanf(fid,"%lf",&d);
      *nodedata++ = xnod0.e(k,j) + d;
      // if (fs_debug) printf("x0 %f, d %f",xnod0.e(k,j),d);
    }
    // if (fs_debug) printf("\n");
  }
#if 0
  nodedata = mesh->nodedata->nodedata;
  for (int k=0; k<nnod; k++)
    printf("node %d, x,y: %f %f\n",k+1,nodedata[k*nu+0],nodedata[k*nu+1]);
#endif
}

void godunov_hook::time_step_post(double time,int step,
			      const vector<double> &gather_values_a) {
  // Pass to global for computing the bottom flux
  gather_values = &gather_values_a;
  // Displacements are read in server and sent to slaves
  if (!MY_RANK) {
    printf("ADV_HOOK: starting time_step_post()\n");
    int ierr;
    fprintf(adv2god,"step %d\n",step);
    int god_step = int(read_doubles(god2adv,"god_step_ok"));
    assert(step==god_step);
    read_mesh();
    printf("ADV_HOOK: ending time_step_post()\n");
  }
  int ierr = MPI_Bcast(mesh->nodedata->nodedata, nnod*nu, MPI_DOUBLE, 0,PETSC_COMM_WORLD);
  assert(!ierr);

  // Volume correction
  int ntime = time_bottom_vel.size();
  double tol=1e-10;
}

void godunov_hook::close() {
  if (!MY_RANK) {
    xnod0.clear();
    fprintf(adv2god,"step %d\n",-1);
  }
}

DL_GENERIC_HOOK(godunov_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// normals:= directions (normalized) along which the reflections
// are imposed
dvector<double> normals;
// states:= current states (they cumulate during time steps)
dvector<double> states;
// fs:= vector of nodes on the ficticious boundary (FB)
dvector<int> fb;
// fb2indx_t:= type of fb2indx
typedef map<int,int> fb2indx_t;
// fb2indx:= map (FB -> index FB list)
fb2indx_t fb2indx;

/** Hook executed in the godunov run */ 
class adv_god_hook {
private:
  /// Fifos for cummunicating with the advdif run
  FILE *adv2god,*god2adv;
  /// Number of nodes, dimension, number of nodes on the FB
  int nnod, ndim, nfb;
  /// time step
  double Dt;
  // Pointer to the mesh, in order to obtain the nodedata
  Mesh *mesh;
  /// Print some values related to the update of the FB
  int fb_debug;
  /// Restart a previous run
  int restart;
  /// Number of times the filter is applied
  int nfilt;
  /// temporary buffer
  dvector<double> states_old, dn, sn, ds, ds2;
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adv_god_hook::init(Mesh &mesh_a,Dofmap &dofmap,
	  TextHashTableFilter *options,const char *name) { 

  if (!MY_RANK) {
    printf("ADV_GOD_HOOK: Opening fifos for communicating with ADV_HOOK.\n");
    
    adv2god = fopen("adv2god.fifo","r");
    assert(adv2god);
    
    god2adv = fopen("god2adv.fifo","w");
    assert(god2adv);
    setvbuf(god2adv,NULL,_IOLBF,0);
    
    printf("ADV_GOD_HOOK: Done.\n");
  }
  int ierr = MPI_Barrier(PETSC_COMM_WORLD);
  assert(!ierr);
  printf("ADV_GOD_HOOK: trace -2\n");

  mesh = &mesh_a;
  // fixme:= some things here will not run in parallel
  nnod = mesh->nodedata->nnod;
  ndim = 2;
  // Read list of nodes on the FB ad normals (all of size nfb).
  FILE *fid = fopen(CASE_NAME ".nod_fb.tmp","r");
  FILE *fid2 = fopen(CASE_NAME ".normals.tmp","r");
  int indx=0;
  printf("ADV_GOD_HOOK: trace -1\n");
  while(1) {
    // printf("indx %d\n",indx);
    int node;
    int nread = fscanf(fid,"%d",&node);
    if (nread==EOF) break;
    double s;
    // Normals are normalized
    double s2 = 0.;
    int pos = normals.size();
    for (int j=0; j<ndim; j++) {
      nread = fscanf(fid2,"%lf",&s);
      // printf("read %f\s",s);
      assert(nread==1);
      s2 += s*s;
      normals.push(s);
    }
    // Normalize normals
    s2 = sqrt(s2);
    for (int j=0; j<ndim; j++) normals.e(pos+j) /= s2;

    // Verify that normals nodes are unique
    assert(fb2indx.find(node)==fb2indx.end());
    fb.push(node);
    // load map
    fb2indx[node] = indx++;
  }
  printf("ADV_GOD_HOOK: trace 0\n");

  fclose(fid);
  fclose(fid2);
  // Number of nodes on the FB
  nfb = fb2indx.size();
  assert(normals.size()==ndim*nfb);
  normals.reshape(2,nfb,ndim);
  
  states.set_chunk_size(ndim*nfb);
  states.a_resize(2,nfb,ndim);
  states.set(0.);

  states_old.set_chunk_size(ndim*nfb);
  states_old.a_resize(2,nfb,ndim);
  states_old.set(0.);
  
  dn.set_chunk_size(nfb);
  dn.a_resize(1,nfb);
  dn.set(0.);
  
  ds.set_chunk_size(nfb);
  ds.a_resize(1,nfb);
  ds.set(0.);
  
  ds2.set_chunk_size(nfb);
  ds2.a_resize(1,nfb);
  ds2.set(0.);
  
  sn.set_chunk_size(nfb);
  sn.a_resize(1,nfb);
  sn.set(0.);
  
  //o Time step.
  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,Dt,0.);
  assert(Dt>0.);
  //o Print some values related to the update of the FB
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,fs_debug,0);
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,restart,0);

  printf("ADV_GOD_HOOK: trace 1\n");
  if (!MY_RANK && !restart) {
    // This is to rewind the file
    FILE *fid = fopen(CASE_NAME ".fsh.tmp","w");
    fclose(fid);
  } 

s  printf("ADV_GOD_HOOK: trace 2\n");
  // Read the last computed mesh
  if (!MY_RANK && restart) {
    fid = fopen(CASE_NAME "_god.state.tmp","r");
    double *xnod = mesh->nodedata->nodedata;
    int nu = mesh->nodedata->nu;
    
    fb2indx_t::iterator q,qe = fb2indx.end();
    for (int node=1; node<=nnod; node++) {
      indx = -1;
      q = fb2indx.find(node);
      if(q!=qe) indx = q->second;
      double d;
      for (int j=0; j<ndim; j++) {
	int nread = fscanf(fid,"%lf",&d);
	assert(nread==1);
	// if (indx>=0) states.e(indx,j) = d - xnod[(node-1)*nu+j];
	if (indx>=0) states.e(indx,j) = d;
      }
    }
    fclose(fid);
  }

  if (!MY_RANK) printf("ADV_GOD_HOOK: ending init().\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adv_god_hook::time_step_pre(double time,int step) {
  PetscPrintf(PETSC_COMM_WORLD,"ADV_GOD_HOOK: step %d, time_step_pre() starts.\n",step);
  int ierr;
  // verify step sent by NS run
  int step_sent = int(read_doubles(adv2god,"step"));
  if (step_sent==-1) {
    PetscPrintf(PETSC_COMM_WORLD, 
		"ADV_GOD_HOOK: step==-1 received, stopping myself.\n");
    PetscFinalize();
    exit(0);
  }
  assert(step=step_sent);
  // Open state file. Will read velocities on the FB
  // and update states += v * Dt
  FILE *fid = fopen(CASE_NAME ".state.tmp","r");
  // string buffer 
  AutoString line;
  vector<string> tokens;
  fs2indx_t::iterator qe = fb2indx.end();
  // Get coordinates pointer
  double *nodedata = mesh->nodedata->nodedata;
  int nu = mesh->nodedata->nu;
  int nnod = mesh->nodedata->nnod;
  FastMat2 x0(1,ndim),x1(1,ndim),x2(1,ndim),
    dx0(1,ndim),dx1(1,ndim),dx2(1,ndim),
    x01(1,ndim),x12(1,ndim),normal(1,ndim);
  // Reads the whole file
  states_old.set(states.buff());
  for (int node=1; node<=nnod; node++) {
    // Reads line (even if it is not in the FB)
    line.getline(fid);
    // if node is not on the FB, then skip
    fb2indx_t::iterator q = fb2indx.find(node);
    if (q==qe) continue;
    // index in the containers
    int indx = q->second;
    if (indx==0 || indx==nfb) continue;
    // Compute normal component of velocity
    int fb_debug = 0;
    tokens.clear();
    tokenize(line.str(),tokens);
    assert(tokens.size()==ndim+1);
    // -- Compute normal --
    // FB nodes at the end of the FB must be treated specially.
    // Currently assumes cyclic FB
    // assert(cyclic_fs);
    // Compute coordinates of three nodes on the free surface
    // Node at the center
    dx1.set(&states_old.e(indx,0));
    x1.set(&nodedata[nu*(node-1)]).add(dx1);

    // Previous node
    int div0=0,div2=0;
    int indx0, node0;
    if (cyclic_fs) {
      indx0 = modulo(indx-1,nfs,&div0);
      if (indx0<0) indx0 = (cyclic_fs ? modulo(indx0,nfs,&div2) : 0);
      node0 = fs.e(indx0);
      x0.set(&nodedata[nu*(node0-1)]);
      dx0.set(&states_old.e(indx0,0));
      x0.addel(div0*cyclic_length,1).add(dx0);
    } else {
      indx0 = indx-1;
      assert(indx0>=0);
    }

    // Next node
    int indx2,node2;
    if (cyclic_fs) {
      indx2 = indx+1;
      if (indx2 >= nfs) indx2 = (cyclic_fs ? modulo(indx2,nfs,&div2) : nfs-1);
      node2 = fs.e(indx2);
      x2.set(&nodedata[nu*(node2-1)]);
      dx2.set(&states_old.e(indx2,0));
      x2.addel(div2*cyclic_length,1).add(dx2);
    } else {
      indx2 = indx+1;
      assert(indx2 >= 0);
    }

    // Segments 0-1 1-2
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
    if (fs_debug) normal.print("normal: ");

    // Compute statesacement
    double v, vn=0., ssn=0., s2=0.;
    n2 = 0.;
    if (fs_debug) printf("node %d, vel, nor, spin ",node);
    for (int j=0; j<ndim; j++) {
      string2dbl(tokens[j],v);
      double s = spines.e(indx,j);
      double n = normal.get(j+1);
      if (fs_debug) printf(" %f %f %f",v,n,s);
      vn += v*n;
      ssn += s*n;
      n2 += n*n;
      s2 += s*s;
    }
    double tol = 1e-10;
    assert(fabs(n2-1.0)<tol);
    assert(fabs(s2-1.0)<tol);
    if (fs_debug) printf(", vn, sn: %f\n",vn,ssn);
    dn.e(indx) = fs_relax * Dt * vn;
    sn.e(indx) = ssn;
  }

  // Compute statesacements along the spines
  for (int indx=0; indx<nfs; indx++) {
    if (!cyclic_fs && indx==0) continue;
    double d=0.;
    for (int j=0; j<ndim; j++) 
      d += spines.e(indx,j)*states.e(indx,j);
    d += dn.e(indx)/sn.e(indx); // This is the correction due to non alignement
				// between the normal and the spine
    ds.e(indx) = d;
  }

  // Smoothing and updating
  if (cyclic_fs) {
    for (int k=0; k<nfilt; k++) {
      for (int indx1=0; indx1<nfs; indx1++) {
	int indx0 = modulo(indx1-1,nfs);
	int indx2 = modulo(indx1+1,nfs);
	ds2.e(indx1) = ds.e(indx1) 
	  + fs_smoothing_coef*(ds.e(indx2)-2*ds.e(indx1)+ds.e(indx0));
      }
      for (int indx=0; indx<nfs; indx++) ds.e(indx) = ds2.e(indx);
    }
  } else {
    // Not implemented yet, do nothing
  }
  
  // Compute vector statesacements from amplitude
  // statesacements projected on the spines
  for (int indx1=0; indx1<nfs; indx1++)
    for (int j=0; j<ndim; j++)
      states.e(indx1,j) = ds.e(indx1)*spines.e(indx1,j);
  
  fclose(fid);
  if (!MY_RANK) {
    fid = fopen(CASE_NAME ".fsh.tmp","a");
    for (int indx=0; indx<nfs; indx++) {
      int node = fs.e(indx);
      for (int j=0; j<ndim; j++) 
	fprintf(fid,"%f ",nodedata[nu*(node-1)+j]+states.e(indx,j));
      fprintf(fid,"\n");
    }
    fclose(fid);
  }
  PetscPrintf(PETSC_COMM_WORLD,"ADV_GOD_HOOK: time_step_pre() ends.\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adv_god_hook::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  fprintf(god2adv,"god_step_ok %d\n",step);
}

void adv_god_hook::close() {}

DL_GENERIC_HOOK(adv_god_hook);
