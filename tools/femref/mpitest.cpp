#include <cassert>
#include <cstdio>
#include <mpi.h>


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main (int argc, char **argv) {
  MPI_Init(&argc,&argv);
  int myrank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
#define N 5
  double v[N];
  if (!myrank) {
    for (int j=0; j<N; j++) v[j] = 23.0;
    MPI_Send(v,N,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
  } else {
    MPI_Status status;
    MPI_Recv(v,N,MPI_DOUBLE,0,0,
	     MPI_COMM_WORLD,&status);
    printf("[%d] received ");
    for (int j=0; j<N; j++) 
      printf("%lg ",v[j]);
    printf("\n");
  }
  MPI_Finalize();
}

