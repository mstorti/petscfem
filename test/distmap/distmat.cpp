/*__INSERT_LICENSE__*/
// $Id: distmat.cpp,v 1.1 2001/07/31 20:06:58 mstorti Exp $
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../../src/distmap.h"
#include <petsc.h>

int SIZE, MYRANK, M;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Test for the distributed map class
// A distributed map<int,int,double> class
class Row : public  map<int,double> {
  int size_of_pack() const;
  void pack(char *&buff) const;
}

typedef DistMap<int,Row> DistMat;
typedef map<int,Row> BasMap;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistMat::
size_of_pack(const DistMat::iterator iter) const {
  int n = size();
  // size + row number + size*(int+double)
  return (n+2)*sizeof(int)+n*sizeof(double);
}

template<class T>
void pack<T>(const T &t,char *& buff) {
  memcpy(buff,&t,sizeof(T));
  buff += sizeof(T);
}
  
template<class T>
void unpack<T>(T &t,const char *& buff) {
  memcpy(&t,buff,sizeof(T));
  buff += sizeof(T);
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
pack(const int &k,const Row &row,char *&buff) const {
  int n = size();
  Row::iterator iter;
  pack<int>(n,buff);
  pack<int>(k,buff);
  for (iter=row.begin(); iter!=row.end(); iter++) {
    pack<int>(iter->first,buff);
    pack<double>(iter->second,buff);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
pack(int &k,Row &row,const char *&buff) const {
  int n,j,key,val;
  unpack<int>(n,buff);
  unpack<int>(k,buff);
  for (j=0; j<n; j++) {
    unpack<int>(key,buff);
    unpack<double>(val,buff);
    row[key] += val;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
combine(const pair<int,Row> &p) {
  map<int,Row>::iterator iter = find(p.first);
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
  MPI_Allreduce(vec.begin(),vecc.begin(),
		M,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

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

  PetscFinalize();
}
