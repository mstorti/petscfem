/*__INSERT_LICENSE__*/
// $Id: pfmat7.cpp,v 1.2 2004/10/24 16:48:43 mstorti Exp $

// Tests for the `PFMat' class

#include <src/debug.h>

#include <stdlib.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/graph.h>
#include <src/petscmat.h>
#include <src/spetscmat.h>

// Runs a simple example for testing the PFMat matrix classes.

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
  ~Part() {}
} part;

int main(int argc,char **args) {
  Vec x,b;
  double tol=1e-10,fill=1.;
  // debug.activate();
  int myrank,size;
  PetscInitialize(&argc,&args,NULL,NULL);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  MY_RANK = myrank;
  SIZE = size;
  FILE *fid = fopen("symm.dat","r");
  int ierr;
  int N;
  int nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  part.comm_size = size;
  part.N = N;
  MY_RANK = myrank;
  SIZE = size;
  PFSymmPETScMat AA(N,N,part,PETSC_COMM_WORLD);
  PFMat &A = AA;

  int nhere=0;
  if (!myrank)  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
		       "now-->setting vector\n");
    CHKERRQ(ierr);
  }
  for (int k=0; k<N; k++) 
    if (part.processor(k)==myrank) nhere++;
  ierr = VecCreateMPI(PETSC_COMM_WORLD,
		      nhere,PETSC_DETERMINE,&b); 
  CHKERRQ(ierr); 
  ierr = VecDuplicate(b,&x); CHKERRQ(ierr); 

  for (int k=0; k<N; k++) {
    double v;
    nread = fscanf(fid,"%lf",&v);
    assert(nread==1);
    if (part.processor(k)==myrank) 
    ierr = VecSetValues(b,1,&k,&v,INSERT_VALUES); 
    CHKERRQ(ierr); 
  }
  ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

  while(1) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    if (nread == EOF) break;
    assert(nread==3);
    A.set_profile(j,k);
  }
  fclose(fid);
  A.create();

#if 0
  nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Setting matrix\n");CHKERRQ(ierr);}
  for (int i=0;i<nonzeros;i++) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if (((j)*size)/N == myrank) A.set_value((j),(k),(v));
    //    A.set_value(j,k,v);
  }
  fclose(fid);
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Assembling matrix\n");CHKERRQ(ierr);}
  A.assembly(MAT_FINAL_ASSEMBLY);
  int nhere=0;
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->setting vector\n");CHKERRQ(ierr);}
  for (int k=0; k<N; k++) 
    if (part.processor(k)==myrank) nhere++;
  ierr = VecCreateMPI(PETSC_COMM_WORLD,
		      nhere,PETSC_DETERMINE,&b); CHKERRQ(ierr); 
  ierr = VecDuplicate(b,&x); CHKERRQ(ierr); 
  fid = fopen("/u/rodrigop/PETSC/petsc-2.1.6/src/contrib/oberman/laplacian_q1/lap_b.dat","r");
  for (int j=0; j<N; j++) {
    double v;
    nread = fscanf(fid,"%lf",&v);
    ierr = VecSetValues(b,1,&j,&v,INSERT_VALUES); CHKERRQ(ierr); 
  }
  fclose(fid);
  int its;
  ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
  if (!myrank)  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"data-->read\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Solving system\n");CHKERRQ(ierr);
  }
  ierr = A.solve(b,x); CHKERRQ(ierr); 
  //  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);
  A.clear();
  PetscFinalize();
#endif
}
