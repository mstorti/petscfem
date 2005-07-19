//__INSERT_LICENSE__
//$Id: mpicookie.cpp,v 1.1 2005/07/19 10:59:55 mstorti Exp $

// This is a test program for testing if setting an
// environment variable from within a program is seen
// from _outside_ the pogram, i.e. lookin at
// /proc/<pid>/environ in Linux. 

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

int main(int argc, char **argv) {
  MPI_Init(&argc,&argv);

  srand(time(NULL));
  unsigned int cookie = rand();

  int  myrank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  MPI_Bcast (&cookie,1, MPI_INT,0,MPI_COMM_WORLD);
  char line[100];
  sprintf(line,"%x",cookie);
  setenv("cookie",line,1);
  
  if (!myrank) {
    printf("press any char to continue > ");
    getchar();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize(); 
}
