#include <src/debug.h>
#include <cstdio>
#include <vector>

//#include <mpi.h>
#include <petsc.h>

#include <src/iisdgraph.h>

int  N;
extern int SIZE,MY_RANK;

class MyPart : public DofPartitioner {
  //  int processor(int j) { return j >= MY_RANK*N/2 || j < (MY_RANK+1)*N/2; }
  int processor(int j) { return j*SIZE/N; }
} d_part;

int main(int argc,char **argv) {

  PetscInitialize(&argc,&argv,NULL,NULL);

  N=10;
  if (MY_RANK==0 && argc>1) sscanf(argv[1],"%d",&N);
  MPI_Bcast (&N, 1, MPI_INT, 0,MPI_COMM_WORLD);

  // MPI_Init(&argc,&argv);
  MPI_Comm_rank (PETSC_COMM_WORLD, &MY_RANK);
  MPI_Comm_size (PETSC_COMM_WORLD, &SIZE);
  Debug debug(0,PETSC_COMM_WORLD);

  debug.trace("after init");
  GPartitioner g_part(&d_part);
  debug.trace("0");

  // debug.activate();
  Debug::init();
  debug.trace("1");
  
  StoreGraph g(N,&g_part,PETSC_COMM_WORLD);
  debug.trace("2");
  for (int k=0; k<N; k++) {
    if (k % SIZE != MY_RANK) continue;
    if (k+1<N) g.add(k,k+1);
    if (k-1>=0) g.add(k,k-1);
  }
  debug.trace("before scatter");
  g.scatter();
  g.print();
  PetscFinalize();
  exit(0);
 
}
