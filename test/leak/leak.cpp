/*__INSERT_LICENSE__*/
// $Id: leak.cpp,v 1.1 2002/10/06 20:24:01 mstorti Exp $

// Tests for the `PFMat' class
#include <petscsles.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {

  PetscInitialize(&argc,&args,NULL,NULL);
  Mat A;
  Vec b,x;
  SLES sles;
  PC pc;
  KSP ksp;

  const int N=10;
  double h = 1./double(N);
  
  int ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,N,N,N,N,
			     3,NULL,3,NULL,&A); CHKERRQ(ierr); 
  ierr = VecCreateMPI(PETSC_COMM_WORLD,N,N,&b); CHKERRQ(ierr); 
  ierr = VecDuplicate(b,&x); CHKERRQ(ierr); 

  ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRQ(ierr); 
  ierr = SLESGetKSP(sles,&ksp); CHKERRQ(ierr); 
  ierr = SLESGetPC(sles,&pc); CHKERRQ(ierr); 

  ierr = SLESSetOperators(sles,A,A,SAME_NONZERO_PATTERN); CHKERRQ(ierr); 
  ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr); 
  ierr = PCSetType(pc,PCJACOBI); CHKERRQ(ierr); 

  int niter=2;

  for (int iter=0; iter<=niter; iter++) {
    for (int k=0; k<N; k++) {
      ierr = MatSetValue(A,k,k,1.,INSERT_VALUES); CHKERRQ(ierr); 
      
      int l=k+1;
      if (l>=N) l-=N;
      ierr = MatSetValue(A,k,l,-0.1,INSERT_VALUES); CHKERRQ(ierr); 
      
      l=k-1;
      if (l<0) l+=N;
      ierr = MatSetValue(A,k,l,-0.1,INSERT_VALUES); CHKERRQ(ierr); 
      ierr = VecSetValue(b,k,1.,INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);
    
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
    
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); 

    int its;
    ierr = SLESSolve(sles,b,x,&its); CHKERRQ(ierr); 

    ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); 

  }
  ierr = MatDestroy(A);
  ierr = VecDestroy(x);
  ierr = VecDestroy(b);
  ierr = SLESDestroy(sles);

  PetscFinalize();
}
