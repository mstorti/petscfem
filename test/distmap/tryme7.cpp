#include <cstdio>
#include <vector>

#include <mpi.h>

#include <src/iisdgraph.h>

const int  N=100;
int SIZE,RANK;

class MyPart : public DofPartitioner {
  //  int processor(int j) { return j >= RANK*N/2 || j < (RANK+1)*N/2; }
  int processor(int j) { return j*SIZE/N; }
} d_part;

int main(int argc,char **argv) {

  GPartitioner g_part(&d_part);

  MPI_Init(&argc,&argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &RANK);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);

  StoreGraph g(N,&g_part);
  for (int k=0; k<N; k++) {
    if (k % SIZE != RANK) continue;
    if (k+1<N) g.add(k,k+1);
    if (k-1>=0) g.add(k,k-1);
  }
  g.scatter();
  g.print();
}
