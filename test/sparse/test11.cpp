/*__INSERT_LICENSE__*/
// $Id: test11.cpp,v 1.2 2006/04/08 21:28:58 mstorti Exp $

#include <cstdio>

#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/dvecpar.h>
#include <src/dvecpar2.h>

struct A {
  double x;
  int y;
};

extern Debug *GLOBAL_DEBUG;
extern int MY_RANK,SIZE;

int main(int argc, char **argv) {
  PetscInitialize(&argc,&argv,NULL,NULL);

  // MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank(MPI_COMM_WORLD, &MY_RANK);

  Debug debug(1,PETSC_COMM_WORLD);
  Debug::init();
  GLOBAL_DEBUG = &debug;

  dvector<double> v;
  dvector_read_parallel("file.dat",v);
  double sum = 0.0;
  for (int k=0; k<v.size(); k++) 
    sum += v.ref(k);
  printf("sum %f\n",sum);

  int n=20;
  dvector<A> av;
  double xsum=0.0;
  int ysum=0;
  if (!MY_RANK) {
    av.resize(n);
    for (int j=0; j<n; j++) {
      av.ref(j).x = 0.1*j;
      av.ref(j).y = 100+j;
    }
  }
  debug.trace("antes de clone");
  dvector_clone_parallel(av);
  debug.trace("despues de clone");

  for (int j=0; j<n; j++) {
    xsum += av.ref(j).x;
    ysum += av.ref(j).y;
  }
  printf("[%d] xsum %f, ysum %d\n",MY_RANK,xsum,ysum);
  PetscFinalize();
}
