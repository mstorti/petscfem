/*__INSERT_LICENSE__*/
// $Id: pfmat4.cpp,v 1.3 2004/07/19 16:38:55 rodrigop Exp $

// Tests for the `PFMat' class

#include <src/debug.h>

#include <stdlib.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/graph.h>

// Runs a simple example for testing the PFMat matrix classes.

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int processor(int k) const { return 0; }
  ~Part() {}
} part;

int main(int argc,char **args) {
  PFMat *A_p;
  Vec x,b;
  double tol=1e-10,fill=1.;
  // debug.activate();
  int myrank,size;
  PetscInitialize(&argc,&args,NULL,NULL);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  MY_RANK = myrank;
  SIZE = size;
  FILE *fid = fopen("sw_matrix.dat","r");
  int N;
  int nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  IISDMat AA(N,N,part,PETSC_COMM_WORLD);
  PFMat &A = AA;
  A.set_option("iisd_subpart",1);
  A.set_option("rtol",1e-13);
  A.set_option("atol",0);
  A.set_option("print_fsm_transition_info",0);
  A.set_option("use_compact_profile",0);
  while (1) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if (j < 0) break;
    A.set_profile(j,k);
  }
  fclose(fid);
  A.create();
  fid = fopen("sw_matrix.dat","r");
  nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  while (1) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if (j < 0) break;
    A.set_value(j,k,v);
  }
  fclose(fid);
  A.assembly(MAT_FINAL_ASSEMBLY);
  int ierr;
  ierr = VecCreateMPI(PETSC_COMM_WORLD,
			  N,N,&b); CHKERRQ(ierr); 
  ierr = VecDuplicate(b,&x); CHKERRQ(ierr); 
  fid = fopen("b.dat","r");
  for (int j=0; j<N; j++) {
    double v;
    nread = fscanf(fid,"%lf",&v);
    ierr = VecSetValues(b,1,&j,&v,INSERT_VALUES); CHKERRQ(ierr); 
  }
  fclose(fid);
  ierr = VecAssemblyBegin(b); CHKERRA(ierr);
  ierr = VecAssemblyEnd(b); CHKERRA(ierr);
  ierr = A.solve(b,x); CHKERRA(ierr); 
  ierr = VecView(x,PETSC_VIEWER_STDOUT_SELF);
  A.clear();
  PetscFinalize();
}
