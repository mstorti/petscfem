//__INSERT_LICENSE__
//$Id: spetscmat.cpp,v 1.3 2004/10/24 17:56:02 mstorti Exp $

#include <src/petscmat.h>
#include <src/pfmat.h>
#include <src/spetscmat.h>
#include <src/graph.h>
#include <src/pfgtpmacr.h>

extern int MY_RANK,SIZE;

static int spetscmat_mult(Mat A,Vec x,Vec y);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFSymmPETScMat::build_sles"
int PFSymmPETScMat::build_sles() {

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

  ierr = SLESDestroy_maybe(sles); CHKERRQ(ierr);
  ierr = SLESCreate(comm,&sles); CHKERRQ(ierr);
  ierr = SLESSetOperators(sles,Asymm,P,SAME_NONZERO_PATTERN); 
  CHKERRQ(ierr);
  ierr = SLESGetKSP(sles,&ksp); CHKERRQ(ierr);
  ierr = SLESGetPC(sles,&pc); CHKERRQ(ierr);

  set_preco(preco_type);

  ierr = KSPSetType(ksp,"cg"); CHKERRQ(ierr);

  ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); 
  CHKERRQ(ierr); 

  ierr = KSPSetMonitor(ksp,PFPETScMat_default_monitor,
		       this,NULL); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFSymmPETScMat::mult"
int PFSymmPETScMat::mult(Vec x,Vec y) {
  ierr = MatMult(A,x,y); CHKERRQ(ierr); 
  ierr = MatMultTransposeAdd(A,x,y,y); CHKERRQ(ierr); 
  double scal=0.5;
  ierr = VecScale(&scal,y); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "spetscmat_mult"
static int spetscmat_mult(Mat A,Vec x,Vec y) {
  void *ctx;
  PFSymmPETScMat *pfsA;
  int ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr); 
  pfsA = (PFSymmPETScMat *) ctx;
  ierr = pfsA->mult(x,y); CHKERRQ(ierr); 
  return 0;
}

