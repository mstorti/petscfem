/*__INSERT_LICENSE__*/
// $Id: distmap.cpp,v 1.1 2001/07/30 00:12:19 mstorti Exp $
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../../src/distmap.h"
#include <petsc.h>

int SIZE, MYRANK;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistMap<int,double>
::size_of_pack(const map<int,double>::iterator iter) const {
  return sizeof(int)+sizeof(double);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMap<int,double>::
pack(const int &k,const double &v,char **buff) const {
  memcpy(*buff,&k,sizeof(int));
  *buff += sizeof(int);
  memcpy(*buff,&v,sizeof(double));
  *buff += sizeof(double);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMap<int,double>::
unpack(int &k,double &v,const char *buff) {
  memcpy(&k,buff,sizeof(int));
  buff += sizeof(int);
  memcpy(&v,buff,sizeof(double));
  buff += sizeof(double);
}

#define M 100
int proc(int l) { return int((l*SIZE)/M);}

int
DistMap<int,double>::processor(const map<int,double>::iterator k) const {
  // return int((k->first*size)/M);
  return proc(k->first);
};

int main(int argc,char **argv) {
  int row;

  srand (time (0));
  /// Initializes MPI
  PetscInitialize(&argc,&argv,0,0);
  // MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MYRANK);

  DistMap<int,double> S;
  int N=100;
  for (int j=0; j<N; j++) {
    row = int(double(rand())/double(RAND_MAX)*double(M));
    S[row] = double(rand())/double(RAND_MAX);
    // PetscSynchronizedPrintf(PETSC_COMM_WORLD,
    // "[%d] row %d, proc %d\n",MYRANK,row,proc(row));
  }
  //  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  S.scatter();
  MPI_Finalize();
}
