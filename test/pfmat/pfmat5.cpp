/*__INSERT_LICENSE__*/
// $Id: pfmat5.cpp,v 1.3 2004/07/19 16:38:55 rodrigop Exp $

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
#include <src/petscmat.h>

// Runs a simple example for testing the PFMat matrix classes.

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
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
  FILE *fid = fopen("/u/rodrigop/PETSC/petsc-2.1.6/src/contrib/oberman/laplacian_q1/lap_A.dat","r");
  int ierr;
  int N;
  int nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  part.comm_size = size;
  part.N = N;
  MY_RANK = myrank;
  SIZE = size;
  IISDMat AA(N,N,part,PETSC_COMM_WORLD);
  PFMat &A = AA;
  A.set_option("iisd_subpart",1);
  A.set_option("use_interface_full_preco",1);
  A.set_option("use_interface_full_preco_nlay",5);
  A.set_option("print_interface_full_preco_conv",0);
  A.set_option("rtol",1e-9);
  A.set_option("atol",1.e-14);
  A.set_option("dtol",1.e5);
  A.set_option("maxits",200);
  A.set_option("print_fsm_transition_info",0);
  A.set_option("use_compact_profile",0);
  A.set_option("print_internal_loop_conv",1);
  //  A.set_option("preco_type","jacobi");
  int nonzeros=0;
  while (1) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if (j < 0) break;
    nonzeros++;
  }
  fclose(fid);
  fid = fopen("/u/rodrigop/PETSC/petsc-2.1.6/src/contrib/oberman/laplacian_q1/lap_A.dat","r");
  nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  for (int i=0;i<nonzeros;i++){
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if ((j*size)/N != myrank) continue; 
    A.set_profile(j,k);
  }
  fclose(fid);
  A.create();
  fid = fopen("/u/rodrigop/PETSC/petsc-2.1.6/src/contrib/oberman/laplacian_q1/lap_A.dat","r");
  nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Setting matrix\n");CHKERRA(ierr);}
  for (int i=0;i<nonzeros;i++) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if (((j)*size)/N == myrank) A.set_value((j),(k),(v));
    //    A.set_value(j,k,v);
  }
  fclose(fid);
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Assembling matrix\n");CHKERRA(ierr);}
  A.assembly(MAT_FINAL_ASSEMBLY);
  int nhere=0;
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->setting vector\n");CHKERRA(ierr);}
  for (int k=0; k<N; k++) 
    if (part.processor(k)==myrank) nhere++;
  ierr = VecCreateMPI(PETSC_COMM_WORLD,
		      nhere,PETSC_DETERMINE,&b); CHKERRA(ierr); 
  ierr = VecDuplicate(b,&x); CHKERRA(ierr); 
  fid = fopen("/u/rodrigop/PETSC/petsc-2.1.6/src/contrib/oberman/laplacian_q1/lap_b.dat","r");
  for (int j=0; j<N; j++) {
    double v;
    nread = fscanf(fid,"%lf",&v);
    ierr = VecSetValues(b,1,&j,&v,INSERT_VALUES); CHKERRQ(ierr); 
  }
  fclose(fid);
  int its;
  ierr = VecAssemblyBegin(b); CHKERRA(ierr);
  ierr = VecAssemblyEnd(b); CHKERRA(ierr);
  if (!myrank)  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"data-->read\n");CHKERRA(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Solving system\n");CHKERRA(ierr);
  }
  ierr = A.solve(b,x); CHKERRA(ierr); 
  //  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);
  A.clear();
  PetscFinalize();
}
