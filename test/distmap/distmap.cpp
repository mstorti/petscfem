/*__INSERT_LICENSE__*/
// $Id: distmap.cpp,v 1.3 2001/07/31 18:10:42 mstorti Exp $
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../../src/distmap.h"
#include <petsc.h>

int SIZE, MYRANK, M;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Test for the distributed map class
// A distributed map<int,double> class

// Simply returns the size of the int+ double
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistMap<int,double>
::size_of_pack(const map<int,double>::iterator iter) const {
  return sizeof(int)+sizeof(double);
}

// Copy the int and double to the buffer. Update pointer *buff
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMap<int,double>::
pack(const int &k,const double &v,char *&buff) const {
  memcpy(buff,&k,sizeof(int));
  buff += sizeof(int);
  memcpy(buff,&v,sizeof(double));
  buff += sizeof(double);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMap<int,double>::
unpack(int &k,double &v,const char *& buff) {
  memcpy(&k,buff,sizeof(int));
  buff += sizeof(int);
  memcpy(&v,buff,sizeof(double));
  buff += sizeof(double);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMap<int,double>::
combine(const pair<int,double> &p) {
  map<int,double>::iterator iter = find(p.first);
  if (find(p.first) != end()) {
    iter->second += p.second;
  } else {
    insert(p);
  }
};

int proc(int l) { return int((l*SIZE)/M);}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double maxd(int n,...) {
  va_list list;
  va_start(list,n);
  double max,item;
  for (int kk=0; kk<n; kk++) {
    item = va_arg(list,double);
    max = ( kk==0 ? item : ( max > item ? max : item));
  }
  va_end(list);
  return max;
}

int
DistMap<int,double>::processor(const map<int,double>::iterator k) const {
  // return int((k->first*size)/M);
  return proc(k->first);
};

int main(int argc,char **argv) {
  int j,N,row,root=0;
  double d,e,err,errb,tol;
  
  map<int,double>::iterator k;
  vector<double> vec,vecc;
  srand (time (0));
  /// Initializes MPI
  PetscInitialize(&argc,&argv,0,0);

  // MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MYRANK);

  if (argc!=4) {
    PetscPrintf(PETSC_COMM_WORLD,"argc: %d\n",argc);
    for (j=0; j < argc; j++) {
      PetscPrintf(PETSC_COMM_WORLD,"argv[%d]: \"%s\"\n",j,argv[j]);
    }
    PetscPrintf(PETSC_COMM_WORLD,"usage: distmap.bin N M tol\n");
    PetscFinalize();
    exit(0);
  }

  sscanf(argv[1],"%d",&M);
  sscanf(argv[2],"%d",&N);
  sscanf(argv[3],"%lf",&tol);

  PetscPrintf(PETSC_COMM_WORLD,"Args: M %d, N %d, tol %g\n",
	      M,N,tol);

  MPI_Bcast (&M, 1, MPI_INT, root,MPI_COMM_WORLD);
  MPI_Bcast (&N, 1, MPI_INT, root,MPI_COMM_WORLD);
  MPI_Bcast (&tol, 1, MPI_DOUBLE, root,MPI_COMM_WORLD);
  
  vec.resize(M,0);
  vecc.resize(M,0);

  DistMap<int,double> S;
  for (int j=0; j<N; j++) {
    row = int(double(rand())/double(RAND_MAX)*double(M));
    e = double(rand())/double(RAND_MAX);
    S[row] += e;
    vec[row] += e;
  }

  S.scatter();
  MPI_Allreduce(vec.begin(),vecc.begin(),M,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

//    if (MYRANK==0) 
//      for (j=0; j<M; j++) printf("%d -> %f\n",j,vecc[j]);

  err = 0;
  for (j=0; j<N; j++) {
    if (proc(j)==MYRANK) {
      k = S.find(j);
      d = (k!=S.end() ? k->second : 0);
      err = maxd(2,err,fabs(d-vecc[j]));
    }
  }
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "[%d] max error -> %g\n",MYRANK,err);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  MPI_Reduce(&err,&errb,1,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD,
	      "max error over all processors -> %g\n"
	      "error < tol OK ? > %d\n",
	      errb,errb<tol);

  MPI_Finalize();
}
