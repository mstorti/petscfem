/*__INSERT_LICENSE__*/
// $Id: distmat.cpp,v 1.7 2001/08/13 01:33:28 mstorti Exp $
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <petsc.h>

#include <utils.h>
#include <maximizr.h>
#include <buffpack.h>
#include <distmap.h>
#include <distmat.h>

int SIZE, MYRANK, M;

class TrivialPartitioner  {
public:
  int processor(int dof) { return int((dof*SIZE)/M); };
  int processor(map<int,Row>::iterator k) {
    return processor(k->first);
  };
};

#define MAT(i,j) VEC2(mat,i,j,M)
int main(int argc,char **argv) {
  TrivialPartitioner part;
    
  int j,k,N,row,col,root=0,debug_v;
  double d,e,w,err,errb,tol;
  double *mat,*matc;
  Maximizer< pair<int,int> > maxerr;
  pair<int,int> p;

  srand (time (0));
  /// Initializes MPI
  PetscInitialize(&argc,&argv,0,0);

  DistMatrix<TrivialPartitioner> S(&part);
  // MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MYRANK);

  if (argc!=5) {
    PetscPrintf(PETSC_COMM_WORLD,"argc: %d\n",argc);
    for (j=0; j < argc; j++) {
      PetscPrintf(PETSC_COMM_WORLD,"argv[%d]: \"%s\"\n",j,argv[j]);
    }
    PetscPrintf(PETSC_COMM_WORLD,"usage: distmap.bin N M tol debug_v\n");
    PetscFinalize();
    exit(0);
  }

  sscanf(argv[1],"%d",&M);
  sscanf(argv[2],"%d",&N);
  sscanf(argv[3],"%lf",&tol);
  sscanf(argv[4],"%d",&debug_v);

  PetscPrintf(PETSC_COMM_WORLD,"Args: M %d, N %d, tol %g, debug_v %d\n",
	      M,N,tol,debug_v);

  MPI_Bcast (&M, 1, MPI_INT, root,MPI_COMM_WORLD);
  MPI_Bcast (&N, 1, MPI_INT, root,MPI_COMM_WORLD);
  MPI_Bcast (&tol, 1, MPI_DOUBLE, root,MPI_COMM_WORLD);
  MPI_Bcast (&debug_v, 1, MPI_INT, root,MPI_COMM_WORLD);
  
  mat = new double[M*M];
  matc = new double[M*M];
  for (j=0; j<M*M; j++) mat[j]=0.;

  for (int j=0; j<N; j++) {
    row = int(double(rand())/double(RAND_MAX)*double(M));
    col = int(double(rand())/double(RAND_MAX)*double(M));
    e = double(rand())/double(RAND_MAX);
    if (debug_v)
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[%d]  (%d,%d) -> %f\n",
			      MYRANK,row,col,e);
    S.insert_val(row,col,e);
    MAT(row,col) += e;
  }
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  S.scatter();
  MPI_Allreduce(mat,matc,
		M*M,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  maxerr.reset();
  for (j=0; j<M; j++) {
    if (part.processor(j)==MYRANK) {
      err = 0;
      for (k=0; k<M; k++) {
	e = S.val(j,k);
	w = MAT(j,k);
	err = maxd(2,err,fabs(e-w));
	p.first=j;
	p.second=k;
	maxerr.scan(p,fabs(e-w));
//  	if (w!=0) 
//  	  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
//  				  "(%d,%d) desired %f, found %f\n",j,k,w,e);
      }
    }
  }
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "[%d] max error at (%d,%d) -> %g\n",MYRANK,
			  maxerr.t.first,maxerr.t.second,err);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  MPI_Reduce(&err,&errb,1,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD,
	      "max error over all processors -> %g\n"
	      "error < tol OK ? > %d\n",
	      errb,errb<tol);

  PetscFinalize();
  delete[] mat;
  delete[] matc;
}
