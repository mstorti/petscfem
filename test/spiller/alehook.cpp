//__INSERT_LICENSE__
//$Id: alehook.cpp,v 1.2 2003/03/19 19:06:41 mstorti Exp $
#define _GNU_SOURCE

#include <cstdio>
#include <cassert>

#include <src/vecmacros.h>
#include <src/texthash.h>
#include <src/texthf.h>
#include <src/fem.h>
#include <src/hook.h>
#include <src/dlhook.h>

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

  PetscPrintf(PETSC_COMM_WORLD,"ALEHOOK: begins time_step_post() ....\n");
  // char * const argv[] = {"/usr/bin/make","spillway_mesh",NULL};
  // int stat = execv(argv[0],argv);
  int stat = system("/usr/bin/make spillway_mesh");
  if (stat==-1) {
    printf("ALEHOOK: Couldn't launch octave process...\n");
    abort();
  }
  FILE *fid = fopen("spiller.nod.tmp","r");
  int nnod = mesh->nodedata->nnod;
  int nu = mesh->nodedata->nu;
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
  PetscPrintf(PETSC_COMM_WORLD,"ALEHOOK: ends time_step_post() ....\n");
}

void ale_hook::close() {}

DL_GENERIC_HOOK(ale_hook);
