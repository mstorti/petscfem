//__INSERT_LICENSE__
//$Id: spetscmat.cpp,v 1.6.40.1 2005/09/25 22:58:44 mstorti Exp $

#include <src/petscmat.h>
#include <src/pfmat.h>
#include <src/spetscmat.h>
#include <src/graph.h>
#include <src/pfgtpmacr.h>

static int spetscmat_mult(Mat A,Vec x,Vec y);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScSymmMat::build_ksp"
PETScSymmMat
::PETScSymmMat(int MM,int NN,
		 const DofPartitioner &part_a,
		 MPI_Comm comm_a) 
  : PETScMat(MM,NN,part_a,comm_a) {
  MPI_Comm_rank(comm,&myrank);
  MPI_Comm_size(comm,&size);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScSymmMat::build_ksp"
int PETScSymmMat::build_ksp() {

  int ierr;
  //o Absolute tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,double,atol,1e-6);
  //o Relative tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,double,rtol,1e-3);
  //o Divergence tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,double,dtol,1e+3);
  //o Maximum iteration number in solving the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,int,maxits,100);
  //o Prints convergence in the solution of the CG iteration. 
  TGETOPTDEF_ND_PF(thash,int,print_internal_loop_conv,0);
  //o Chooses the preconditioning operator. 
  TGETOPTDEF_S_PF(thash,string,preco_type,jacobi);

  int rstart,rend,nlocal;
  ierr = MatGetOwnershipRange(A,&rstart,&rend);
  nlocal = rend-rstart;

  assert(M==N); // this should go elsewhere
  ierr = MatCreateShell(comm,nlocal,nlocal,M,M,this,&Asymm);
  CHKERRQ(ierr); 
  P=A;

  MatShellSetOperation(Asymm,MATOP_MULT,
		       (void (*)(void))(&spetscmat_mult));
  MatShellSetOperation(Asymm,MATOP_MULT_TRANSPOSE,
		       (void (*)(void))(&spetscmat_mult));

  ierr = KSPDestroy_maybe(ksp); CHKERRQ(ierr);
  ierr = KSPCreate(comm,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,Asymm,P,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);

  set_preco(preco_type);

  ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr);

  ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits);CHKERRQ(ierr); 

  ierr = KSPSetMonitor(ksp,PFPETScMat_default_monitor,
		       this,NULL); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScSymmMat::mult"
int PETScSymmMat::mult(Vec x,Vec y) {
  ierr = MatMult(A,x,y); CHKERRQ(ierr); 
  ierr = MatMultTransposeAdd(A,x,y,y); CHKERRQ(ierr); 
  double scal=0.5;
  ierr = VecScale(y,scal); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int PETScSymmMat::insert_p(int row,int col) {
#if 0
  return col>=row;
#else
  return (col-row) % 2;
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int PETScSymmMat::insert_p2(int row,int col) {
  int x;
  if (row==col) x = 1;
  else if (col>row) x = insert_p(row,col);
  else x = !insert_p(col,row);
  return x;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScSymmMat::set_profile_a"
int PETScSymmMat::set_profile_a(int row,int col) {
  if (insert_p2(row,col)) lgraph->add(row,col);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScSymmMat::set_value_a"
int PETScSymmMat::set_value_a(int row,int col,PetscScalar value,
	      InsertMode mode) {
  int roww = row, coll = col;
  if (!insert_p2(row,col)) { roww=col; coll=row; }
  // printf("MatSetValue at %d, %d\n",roww,coll);
  ierr = MatSetValues(A,1,&roww,1,&coll,&value,mode); 
  CHKERRQ(ierr); 
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "spetscmat_mult"
static int spetscmat_mult(Mat A,Vec x,Vec y) {
  void *ctx;
  PETScSymmMat *pfsA;
  int ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr); 
  pfsA = (PETScSymmMat *) ctx;
  ierr = pfsA->mult(x,y); CHKERRQ(ierr); 
  return 0;
}
