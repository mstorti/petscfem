/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.4 2002/07/18 21:43:37 mstorti Exp $

#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <src/utils.h>
#include <src/iisdgraph.h>
#include <sles.h>

TextHashTable GLOBAL_OPTIONS;

#define ALLOC malloc_alloc

typedef set< int, less<int>, ALLOC > SET;

class MAP : public map< int, SET, less<int>, ALLOC > {
public:
  void add(int i,int j) { (*this)[i].insert(j); }
};

int main(int argc, char **argv) {
  PetscInitialize(&argc,&argv,NULL,NULL);

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);
  MPI_Comm_rank(PETSC_COMM_WORLD,&MY_RANK);

  int M=1000000;
  // StoreGraph g(M,&seq_partitioner);
  MAP g;
  for (int j=0; j<M; j++) { g.add(irand(1,M),irand(1,M)); }
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  g.clear();
  delete[] a;

  PetscFinalize();
  exit(0);
}
