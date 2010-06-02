//__INSERT_LICENSE__
//$Id: tryme7.cpp,v 1.5 2003/01/08 15:49:04 mstorti Exp $
#include <src/debug.h>
#include <cstdio>
#include <vector>

//#include <mpi.h>
#include <petsc.h>

#include <src/distmap.h>
#include <src/distmap2.h>
#include <src/distcont.h>
#include <src/distcont2.h>
#include <src/iisdgraph.h>
#include <src/linkgraph.h>

int  N;
extern int SIZE,MY_RANK;

class MyPart : public DofPartitioner {
public:
  //  int processor(int j) { return j >= MY_RANK*N/2 || j < (MY_RANK+1)*N/2; }
  int processor(int j) const { return j*SIZE/N; }
} d_part;

void check_fun(int &ok,int k,int e,int N,
               set<int> &ngbrs, int rank) {
  if (e>=0 && e<N) {
    int erase = ngbrs.erase(e);
    if (erase!=1) {
      printf("[%d] key %d, elem %d, erase %d\n",
             rank,k,e,erase);
      ok = 0;
    }
  }
}

int main(int argc,char **argv) {

  PetscFemInitialize(&argc,&argv,NULL,NULL);

  N=10;
  if (MY_RANK==0 && argc>1) sscanf(argv[1],"%d",&N);
  MPI_Bcast (&N, 1, MPI_INT, 0,MPI_COMM_WORLD);

  int ntime=1;
  if (MY_RANK==0 && argc>2) sscanf(argv[2],"%d",&ntime);
  MPI_Bcast (&ntime, 1, MPI_INT, 0,MPI_COMM_WORLD);

  // MPI_Init(&argc,&argv);
  MPI_Comm_rank (PETSC_COMM_WORLD, &MY_RANK);
  MPI_Comm_size (PETSC_COMM_WORLD, &SIZE);
  Debug debug(0,PETSC_COMM_WORLD);

  // debug.activate();
  Debug::init();
  debug.trace("1");

  for (int j=0; j<ntime; j++) {
    // if (j % 100 ==0) 
    //   PetscPrintf(PETSC_COMM_WORLD,"j %d\n",j);
    LinkGraphWrapper g(0,&d_part,PETSC_COMM_WORLD);
    g.init(N);
    // StoreGraph g(N,&d_part,PETSC_COMM_WORLD);
    debug.trace("2");
    for (int k=0; k<N; k++) {
      if (k % SIZE != MY_RANK) continue;
      if (k+1<N) g.add(k,k+1);
      if (k-1>=0) g.add(k,k-1);
    }
    debug.trace("before scatter");
    g.scatter();
    int ok = 1;
    for (int k=0; k<N; k++) {
      set<int> ngbrs;
      g.set_ngbrs(k,ngbrs);
      int proc = d_part.processor(k);
      // printf("[%d] k %d, proc %d\n",MY_RANK,k,proc);
      if (proc!=MY_RANK) {
        if (!ngbrs.empty()) {
          printf("[%d] key %d belongs to proc %d and non-empty\n",
                 MY_RANK,k,proc);
          ok=0;
        }
      } else {
        check_fun(ok,k,k-1,N,ngbrs,MY_RANK);
        check_fun(ok,k,k+1,N,ngbrs,MY_RANK);
        ok &= ngbrs.empty();
      }
#if 0
      set<int>::iterator q = ngbrs.begin();
      while (q != ngbrs.end()) {
        int proc = d_part.processor(*q);
        if (proc!=MY_RANK) {
          printf("[%d] k %d, element %d not OK, belongs to %d\n",
                 MY_RANK,k,*q,proc);
          ok=0;
        }
        q++;
      }
#endif
    }
    if (j==ntime-1) {
      // g.print();
      // printf("[%d] OK? %d\n",MY_RANK,ok);
      int okglob;
      MPI_Reduce(&ok,&okglob,1,MPI_INT,MPI_LAND,0,MPI_COMM_WORLD);
      if (!MY_RANK) printf("ok on all procs ? %d\n",okglob);
    }
    g.clear();
  }
  PetscFinalize();
  exit(0);
}
