/*__INSERT_LICENSE__*/
// $Id: pfmat.cpp,v 1.1.2.4 2001/12/24 18:18:25 mstorti Exp $

// Tests for the `PFMat' class
#include <src/debug.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/graph.h>

extern int MY_RANK,SIZE;

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
  ~Part() {}
} part;

int main(int argc,char **args) {

  Vec x,b;
  int ierr;
  debug.activate();
  int myrank,size;
  PetscInitialize(&argc,&args,NULL,NULL);
  const int N=10;

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  part.comm_size = size;
  part.N = N;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  MY_RANK = myrank;
  SIZE = size;
  IISDMat A(N,N,part,PETSC_COMM_WORLD);
  
  for (int j=0; j<N; j++) {
    if (j % size != myrank) continue; // Load periodically
    if (j+1<N) A.set_profile(j,j+1);
    if (j-1>=0) A.set_profile(j,j-1);
  }
  debug.trace("main 0");
  A.set_option("iisd_subpart","2");
  A.set_option("print_internal_loop_conv","1");
  A.create();
  A.print();

  int nhere=0;
  for (int k=0; k<N; k++) 
    if (part.processor(k)==myrank) nhere++;
  
  ierr = VecCreateMPI(PETSC_COMM_WORLD,
		      nhere,PETSC_DETERMINE,&b); CHKERRA(ierr); 
  ierr = VecDuplicate(b,&x); CHKERRA(ierr); 

  A.zero_entries();
  for (int j=0; j<N; j++) {
    if (j % size != myrank) continue; // Load periodically
    A.set_value(j,j,2.);
    if (j+1<N) A.set_value(j,j+1,-1.);
    if (j-1>=0) A.set_value(j,j-1,-1.);
    double val=1.;
    ierr = VecSetValues(b,1,&j,&val,INSERT_VALUES); CHKERRA(ierr);
  }
  ierr = VecAssemblyBegin(b); CHKERRA(ierr);
  ierr = VecAssemblyEnd(b); CHKERRA(ierr);
  A.assembly_begin(MAT_FINAL_ASSEMBLY); 
  A.assembly_end(MAT_FINAL_ASSEMBLY); 

  ierr = A.solve(b,x); CHKERRA(ierr); 
  ierr = VecView(x,VIEWER_STDOUT_WORLD);

  PetscFinalize();
  exit(0);
 
}
