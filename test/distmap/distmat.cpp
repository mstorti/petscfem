/*__INSERT_LICENSE__*/
// $Id: distmat.cpp,v 1.2 2001/08/01 02:25:43 mstorti Exp $
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
};

typedef DistMap<int,Row> DistMat;
typedef map<int,Row> BasMap;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistMat::
size_of_pack(const DistMat::iterator iter) const {
  int n = iter->second.size();
  // size + row number + size*(int+double)
  return (n+2)*sizeof(int)+n*sizeof(double);
}

namespace BufferPack {

  template<class T>
  void pack(const T &t,char *& buff) {
    memcpy(buff,&t,sizeof(T));
    buff += sizeof(T);
  }
  
  template<class T>
  void unpack(T &t,const char *& buff) {
    memcpy(&t,buff,sizeof(T));
    buff += sizeof(T);
  }
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
pack(const int &k,const Row &row,char *&buff) const {
  int n = row.size();
  Row::const_iterator iter;
  BufferPack::pack(n,buff);
  BufferPack::pack(k,buff);
  for (iter=row.begin(); iter!=row.end(); iter++) {
    BufferPack::pack(iter->first,buff);
    BufferPack::pack(iter->second,buff);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
unpack(int &k,Row &row,const char *&buff) {
  int n,j,key;
  double val;
  BufferPack::unpack(n,buff);
  BufferPack::unpack(k,buff);
  for (j=0; j<n; j++) {
    BufferPack::unpack(key,buff);
    BufferPack::unpack(val,buff);
    row[key] += val;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistMat::
combine(const pair<int,Row> &p) {
  int n,j;
  map<int,Row>::iterator iter = find(p.first);
  Row::iterator r;
  Row::const_iterator q;
  if (iter == end()) {
    printf("[%d] combining row %d\n",myrank,p.first);
    insert(p);
  } else {
    Row &oldr = iter->second;
    const Row &newr = p.second;
    for (q = newr.begin(); q!=newr.end(); q++) {
      r = oldr.find(q->first);
      if (r == oldr.end()) {
	r->second += q->second;
      } else {
	oldr.insert(*q);
      }
    }
  }
}

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
DistMat::processor(const DistMat::iterator k) const {
  return proc(k->first);
};

class DistMatrix : public DistMat {
public:
  void insert_val(int i,int j,double v);
  double val(int i,int j);
};

void DistMatrix::insert_val(int i,int j,double v) {
  DistMat::iterator I = find(i);
  Row::iterator J;
  if (I == end()) {
    insert(pair<int,Row>(i,Row()));
    I = find(i);
  }
  Row &row = I->second;
  J = row.find(j);
  if (J == row.end()) {
    row.insert(pair<int,double>(j,v));
  } else {
    J->second += v;
  }
}

double DistMatrix::val(int i,int j) {
  Row::iterator J;
  DistMat::iterator I = find(i);
  if (I == end()) return 0.;
  Row &row = I->second;
  J = row.find(j);
  if (J == row.end()) {
    return 0.;
  } else {
    return J->second;
  }
}  

template<class T> 
class Maximizer {
private:
  int ini;
  T t_;
  double max_;
public:
  void reset() {ini=0;};
  const T &t;
  const double &max;
  Maximizer() : ini(0), max_(0), max(max_), t(t_) {};
  void scan(T &tt,double val);
};

template<class T> 
void Maximizer<T>::scan(T &tt,double val) {
  if (!ini || val>max_) {
    t_=tt; 
    max_=val;
    ini=1;
  }
};

#define MAT(i,j) VEC2(mat,i,j,M)
int main(int argc,char **argv) {
#if 0 // Test for the maximizer class
  double e;
  Maximizer<int> maxerr;
  maxerr.reset();
  for (int j=0; j<100; j++) {
    e = double(rand())/double(RAND_MAX);
    printf("%d -> %f\n",j,e);
    maxerr.scan(j,e);
  }
  printf("max val: j %d, val %f\n",maxerr.t,maxerr.max);
#endif
    
  int j,k,N,row,col,root=0,debug_v;
  double d,e,w,err,errb,tol;
  double *mat,*matc;
  Maximizer< pair<int,int> > maxerr;

  srand (time (0));
  /// Initializes MPI
  PetscInitialize(&argc,&argv,0,0);

  DistMatrix S;
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
    if (proc(j)==MYRANK) {
      err = 0;
      for (k=0; k<M; k++) {
	e = S.val(j,k);
	w = MAT(j,k);
	err = maxd(2,err,fabs(e-w));
	maxerr.scan(pair<int,int>(j,k),fabs(e-w));
      }
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
  delete[] mat;
  delete[] matc;
}
