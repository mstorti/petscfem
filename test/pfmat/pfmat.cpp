/*__INSERT_LICENSE__*/
// $Id: pfmat.cpp,v 1.1.2.3 2001/12/24 12:53:52 mstorti Exp $

// Tests for the `PFMat' class
#include <src/debug.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/graph.h>

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
  ~Part() {}
} part;

int main(int argc,char **args) {

  debug.activate();
  int myrank,size;
  PetscInitialize(&argc,&args,NULL,NULL);
  const int N=10;

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  part.comm_size = size;
  part.N = N;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  MY_RANK = myrank;
  SIZE = size;
  IISDMat A(N,N,part,PETSC_COMM_WORLD);
  
  for (int j=0; j<N; j++) {
    if (j % size != myrank) continue; // Load periodically
    if (j+1<N) A.set_profile(j,j+1);
    if (j-1>=0) A.set_profile(j,j-1);
  }
  debug.trace("main 0");
  A.create();
  A.print();
  PetscFinalize();
  exit(0);
 
}
