#include <src/debug.h>
#include <cstdio>
#include <vector>

//#include <mpi.h>
#include <petsc.h>

#include <src/iisdgraph.h>
#include <src/linkgraph.h>

int  N;
extern int SIZE,MY_RANK;

class MyPart : public DofPartitioner {
  //  int processor(int j) { return j >= MY_RANK*N/2 || j < (MY_RANK+1)*N/2; }
  int processor(int j) const { return j*SIZE/N; }
} d_part;

int main(int argc,char **argv) {

  PetscInitialize(&argc,&argv,NULL,NULL);

  N=10;
  if (MY_RANK==0 && argc>1) sscanf(argv[1],"%d",&N);
  MPI_Bcast (&N, 1, MPI_INT, 0,MPI_COMM_WORLD);

  int ntime=1;
  if (MY_RANK==0 && argc>2) sscanf(argv[2],"%d",&ntime);
  MPI_Bcast (&ntime, 1, MPI_INT, 0,MPI_COMM_WORLD);

  // MPI_Init(&argc,&argv);
  MPI_Comm_rank (PETSC_COMM_WORLD, &MY_RANK);
  MPI_Comm_size (PETSC_COMM_WORLD, &SIZE);
  Debug debug(0,PETSC_COMM_WORLD);

  // debug.activate();
  Debug::init();
  debug.trace("1");

  for (int j=0; j<ntime; j++) {
    if (j % 100 ==0) 
      PetscPrintf(PETSC_COMM_WORLD,"j %d\n",j);
    LinkGraphWrapper g(0,&d_part,PETSC_COMM_WORLD);
    g.init(N);
    // StoreGraph g(N,&d_part,PETSC_COMM_WORLD);
    debug.trace("2");
    for (int k=0; k<N; k++) {
      if (k % SIZE != MY_RANK) continue;
      if (k+1<N) g.add(k,k+1);
      if (k-1>=0) g.add(k,k-1);
    }
    debug.trace("before scatter");
    g.scatter();
    if (j==ntime-1) g.print();
    g.clear();
  }
  PetscFinalize();
  exit(0);
 
}
