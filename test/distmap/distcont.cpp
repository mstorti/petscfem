/*__INSERT_LICENSE__*/
// $Id: distcont.cpp,v 1.6 2003/07/03 04:32:11 mstorti Exp $
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include <map>

#include <petsc.h>

extern int SIZE, MY_RANK;
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/utils.h>

int M;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
int wait_from_console(char *s=NULL) {
  static int deac=0;
  if (deac) return 0;
  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  int ierr;
  char ans;
  ierr = MPI_Barrier(PETSC_COMM_WORLD);
  if (myrank==0) {
    if (s!=NULL) printf("%s --- ",s);
    printf("Continue? (n/RET=y) > ");
    fflush(stdout);
    scanf("%c",&ans);
  }
  ierr = MPI_Bcast (&ans, 1, MPI_CHAR, 0,PETSC_COMM_WORLD);
  CHKERRQ(ierr); 
  if (ans=='n') {
    PetscFinalize();
    exit(0);
  } else if (ans=='d') {
    deac = 1;
  } 
  ierr = MPI_Barrier(PETSC_COMM_WORLD);
  CHKERRQ(ierr);  
}
#endif

typedef map<int,double> Map_id;
typedef pair<int,double> VT; // ValueType

class TrivialPartitioner {
public:
  int processor(int j) { return int((j*SIZE)/M);};
  void processor(const VT &k,int &nproc,int *plist);
};

void 
TrivialPartitioner::processor(const VT &k,int &nproc,int *plist) {
  nproc=1;
  plist[0] = processor(k.first);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Test for the distributed map class
// A distributed map<int,double> class

// Simply returns the size of the int+ double
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int DistCont<Map_id,VT,TrivialPartitioner>
::size_of_pack(const VT &p) const {
  return sizeof(int)+sizeof(double);
}

// Copy the int and double to the buffer. Update pointer *buff
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<Map_id,VT,TrivialPartitioner>::
pack(const VT &p,char *&buff) const {
  memcpy(buff,&p.first,sizeof(int));
  buff += sizeof(int);
  memcpy(buff,&p.second,sizeof(double));
  buff += sizeof(double);
}
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<Map_id,VT,TrivialPartitioner>::
unpack(VT &p,const char *& buff) {
  memcpy(&p.first,buff,sizeof(int)); // debug:=
  buff += sizeof(int);
  memcpy(&p.second,buff,sizeof(double)); // debug:=
  buff += sizeof(double);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DistCont<Map_id,VT,TrivialPartitioner>::
combine(const VT &p) {
  Map_id::iterator iter = find(p.first);
  if (iter != end()) {
    iter->second += p.second;
  } else {
    insert(p);
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
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
#endif

typedef DistCont<Map_id,VT,TrivialPartitioner> Mapp;

int main(int argc,char **argv) {
  int j,N,row,root=0;
  double d,e,err,errb,tol;
  TrivialPartitioner part;
  
  Map_id::iterator k;
  vector<double> vec,vecc;
  srand (time (0));
  /// Initializes MPI
  PetscInitialize(&argc,&argv,0,0);

  // wait_from_console("starting..."); 

  // MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

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

  Mapp S(&part);
  for (int j=0; j<N; j++) {
    row = int(double(rand())/double(RAND_MAX)*double(M));
    e = double(rand())/double(RAND_MAX);
    S[row] += e;
    vec[row] += e;
    // printf("[%d] loading S[%d] += %f\n",MY_RANK,row,e);
  }

  S.scatter();
  MPI_Allreduce(&*vec.begin(),&*vecc.begin(),M,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  err = 0;
  for (j=0; j<N; j++) {
    if (part.processor(j)==MY_RANK) {
      k = S.find(j);
      d = (k!=S.end() ? k->second : 0);
      e = vecc[j];
      err = maxd(2,err,fabs(d-e));
//        if (d!=0 || e!=0) printf("[%d]  S[%d] = %f, (expected %f)\n",
//  			       MY_RANK,j,d,vecc[j]);
      if (fabs(d-e)>tol ) printf("[%d]  S[%d] = %f, (expected %f, err %f)\n",
			       MY_RANK,j,d,e,fabs(d-e));
    }
  }
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "[%d] max error -> %g\n",MY_RANK,err);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  MPI_Reduce(&err,&errb,1,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD,
	      "max error over all processors -> %g\n"
	      "error < tol OK ? > %d\n",
	      errb,errb<tol);

  MPI_Finalize();
}
