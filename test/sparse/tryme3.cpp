/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.3 2002/07/18 20:01:44 mstorti Exp $

#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <src/utils.h>
#include <src/iisdgraph.h>
#include <sles.h>

TextHashTable GLOBAL_OPTIONS;

int main(int argc, char **argv) {
  PetscInitialize(&argc,&argv,NULL,NULL);

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);
  MPI_Comm_rank(PETSC_COMM_WORLD,&MY_RANK);

  int M=1000000;
  StoreGraph g(M,&seq_partitioner);
  for (int j=0; j<M; j++) { g.add(irand(1,M*M),irand(1,M*M)); }
  g.clear();
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  delete[] a;

  PetscFinalize();
  exit(0);
}
