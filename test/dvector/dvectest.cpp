//__INSERT_LICENSE__
// $Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$

#include <cstdio>
#include <mpi.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/dvecpar.h>
#include <src/dvecpar2.h>

void print_shapes(int myrank,dvector<double> &a) {
  vector<int> shape;
  a.get_shape(shape);
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
                          "[%d] rank %d, shape: ",myrank,a.rank());
  for (unsigned int j=0; j<shape.size(); j++)
    PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
                            "%d ",shape[j]);
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,"\n");
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD);
}

int main(int argc,char **args) {
  PetscFemInitialize(&argc,&args,NULL,NULL);
  int myrank,size;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  dvector<double> a;
  int ok = 1;

  if (!myrank) a.a_resize(2,10,3);
  a.set(1.0);
  dvector_clone_parallel(a);
  // print_shapes(myrank,a);
  int ok1,ok1l=1;
  ok1l &= (a.rank()==2);
  ok1l &= a.size(0)==10;
  ok1l &= a.size(1)==3;
  MPI_Reduce(&ok1l,&ok1,1,MPI_INT,MPI_LAND,0,MPI_COMM_WORLD);
  PetscPrintf(PETSCFEM_COMM_WORLD,"shape (3,10) OK? %d\n",ok1);
  ok &= ok1;

  dvector<double> b;
  if (!myrank) b.resize(1000);
  b.set(1.0);
  dvector_clone_parallel(b);
  // print_shapes(myrank,b);
  ok1l = 1;
  ok1l &= (b.rank()==0);
  ok1l &= (b.size()==1000);
  MPI_Reduce(&ok1l,&ok1,1,MPI_INT,MPI_LAND,0,MPI_COMM_WORLD);
  PetscPrintf(PETSCFEM_COMM_WORLD,"not shaped dvector OK? %d\n",ok1);
  ok &= ok1;
  
  if (!myrank) b.clear();
  b.set(1.0);
  dvector_clone_parallel(b);
  // print_shapes(myrank,b);
  ok1l = 1;
  ok1l &= (b.rank()==0);
  ok1l &= (b.size()==0);
  MPI_Reduce(&ok1l,&ok1,1,MPI_INT,MPI_LAND,0,MPI_COMM_WORLD);
  PetscPrintf(PETSCFEM_COMM_WORLD,"empty vector  OK? %d\n",ok1);
  ok &= ok1;

  PetscPrintf(PETSCFEM_COMM_WORLD,"All tests OK? %d\n",ok);

  PetscFinalize();
  return 0;
}
