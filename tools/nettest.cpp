//__INSERT_LICENSE__
//$Id: nettest.cpp,v 1.5 2004/07/12 17:42:21 mstorti Exp $
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <cassert>
#include <unistd.h>
#include <ctype.h>

#include <petsc.h>

/** Computes random numbers between 0 and 1
*/ 
double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

/** Computes an intger random number in the
    range `imin'-`imax'
*/ 
int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

//#define CHECK

int main(int argc,char **argv) {

  int size,rank,p1=0,p2=1;
  int chunk_size=100000,ntimes=3;
  double tol=1e-10,sum_check;

  MPI_Status stat;
  
  /// time related quantities
  double start;

  /// Initializes MPI
  PetscInitialize(&argc,&argv,NULL,NULL);
  /// Initializes random
  srand(time (0));

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#define MAXNAME 100
  char name[MAXNAME];
  int resultlen, ierr;
  ierr =  MPI_Get_processor_name(name,&resultlen);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] name: \"%s\"\n",rank,name);
  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  while (chunk_size>1) {
  
    PetscPrintf(MPI_COMM_WORLD,
		"Results for chunk_size = %d double, (%d bytes, %d bits) ----------\n",
		chunk_size,chunk_size*8,chunk_size*64);
    double *buff = new double[chunk_size];
    for (int j=0; j<chunk_size; j++) buff[j]=(double)j;
    
    sum_check=double(chunk_size)*double(chunk_size-1)/2.;
    
    for (p1=0; p1<size-1; p1++) {
      for (p2=p1+1; p2<size; p2++) {

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	for (int jtime=0; jtime < ntimes; jtime++) {
	  if (rank==p1) {
	    
	    MPI_Send(buff,chunk_size,MPI_DOUBLE,p2,p1,MPI_COMM_WORLD);
	    
	    MPI_Recv(buff,chunk_size,MPI_DOUBLE,p2,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
#ifdef CHECK
	    double sum=0;
	    for (int j=0; j<chunk_size; j++) sum += buff[j];
	    assert(fabs(sum-sum_check)<tol);
#endif
	    
	  } else if (rank==p2) {
	    
	    MPI_Recv(buff,chunk_size,MPI_DOUBLE,p1,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
#ifdef CHECK
	    double sum=0;
	    for (int j=0; j<chunk_size; j++) sum += buff[j];
	    assert(fabs(sum-sum_check)<tol);
#endif
	    MPI_Send(buff,chunk_size,MPI_DOUBLE,p1,p2,MPI_COMM_WORLD);
	    
	  }
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank==0) {
	  double elapsed = MPI_Wtime()-start;
	  double speed = 2.*double(ntimes)*double(chunk_size)*64./elapsed/1e6;
	  printf("proc %d <-> proc %d connection bandwidth: %g Mbit/sec\n",p1,p2,speed);
	}
      }
    }
    delete[] buff;
    chunk_size /= 10;
    ntimes *= 10;
  }

  PetscFinalize();
  exit(0);

}
