/*__INSERT_LICENSE__*/
// $Id: sw_nn_preco.cpp,v 1.1.2.1 2004/07/07 15:47:57 mstorti Exp $

// Tests for the `PFMat' class

#include <src/debug.h>

#include <stdlib.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/petscmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/graph.h>
#include <petscsles.h>

// Runs a simple example for testing the PFMat matrix classes.

extern int MY_RANK,SIZE;
//extern int PCCreate_NN(PC);

class Part : public DofPartitioner {
public:
  int N,comm_size;
  int processor(int k) const { return k*comm_size/N; }
  ~Part() {}
} part;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args) {
  PFMat *A_p;
  Vec x,b;
  double tol=1e-10,fill=1.;
  // debug.activate();
  int myrank,size;
  PetscInitialize(&argc,&args,NULL,NULL);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  Debug debug(0,PETSC_COMM_WORLD);
  MY_RANK = myrank;
  SIZE = size;
  FILE *fid = fopen("sw_matrix_nn.dat","r");
  int N;
  int pdebug_wait=1;
  //  while (pdebug_wait);
  int nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  part.comm_size = size;
  part.N = N;
  MY_RANK = myrank;
  SIZE = size;
  PETScMat AA(N,N,part,PETSC_COMM_WORLD);
  PFMat &A = AA;
  A.set_option("rtol",1e-9);
  A.set_option("atol",1.e-13);
  A.set_option("dtol",1.e5);
  A.set_option("maxits",200);
  A.set_option("Krylov_dim",200);
  A.set_option("print_fsm_transition_info",0);
  A.set_option("use_compact_profile",0);
  A.set_option("print_internal_loop_conv","1");
  A.set_option("solver","petsc");
  A.set_option("KSP_method","gmres");
  A.set_option("preco_type","nn");
  A.set_option("mat_type","is");
  int nonzeros=0;
  int ierr;
  Debug::init();
  debug.activate();
  while (1) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if (j < 0) break;
    nonzeros++;
  }
  fclose(fid);
  fid = fopen("sw_matrix_nn.dat","r");
  nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Setting profile\n");CHKERRA(ierr);}
  for (int i=0;i<nonzeros;i++){
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if ((j*size)/N != myrank) continue; 
    A.set_profile(j,k);
  }
  fclose(fid);
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Creating matrix\n");CHKERRA(ierr);}
  debug.trace("creating matrix");
  A.create();
  fid = fopen("sw_matrix_nn.dat","r");
  nread = fscanf(fid,"%d",&N);
  assert(nread==1);
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Setting matrix\n");CHKERRA(ierr);}
  for (int i=0;i<nonzeros;i++) {
    int j,k;
    double v;
    nread = fscanf(fid,"%d %d %lf",&j,&k,&v);
    assert(nread==3);
    if (((j)*size)/N == myrank) A.set_value((j),(k),(v));
    //    if (!myrank)  printf("i:= %d\n",i);
  }
  fclose(fid);
  int nhere=0;
  debug.trace("vector create");
  for (int k=0; k<N; k++) 
    if (part.processor(k)==myrank) nhere++;
  ierr = VecCreateMPI(PETSC_COMM_WORLD,nhere,PETSC_DETERMINE,&x); CHKERRA(ierr); 
  ierr = VecCreateGhost(PETSC_COMM_WORLD,nhere,PETSC_DETERMINE,A.ghosts,A.g_ghosts,&b); CHKERRA(ierr); 
  //  ierr = VecDuplicate(x,&b); CHKERRA(ierr); 
  ierr = VecSetLocalToGlobalMapping(b,A.ltog);CHKERRQ(ierr);
  fid = fopen("b.dat","r");
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Setting rhs\n");CHKERRA(ierr);}
  for (int j=0; j<N; j++) {
    double v;
    nread = fscanf(fid,"%lf",&v);
    ierr = A.set_value_v(b,j,v,INSERT_VALUES); CHKERRQ(ierr);
  }
  fclose(fid);
  ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

  int its;
  if (!myrank)  {ierr = PetscPrintf(PETSC_COMM_WORLD,"now-->Assembling matrix\n");CHKERRA(ierr);}
  A.assembly(MAT_FINAL_ASSEMBLY);
  if (!myrank)  {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------\n");CHKERRA(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"data-->Read\n");CHKERRA(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"now -->Solving system\n");CHKERRA(ierr);
  }
  debug.trace("solving system");
  ierr = A.solve(b,x); CHKERRA(ierr); 
  A.clear();
  PetscFinalize();
}
