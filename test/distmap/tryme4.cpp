#include <petsc.h>
#include <sles.h>

int main(int argc,char **args) {
  const int N=100;
  const int M=10;
  int ierr,m=3,its;
  double v, *a, sum;
  SLES sles;
  KSP ksp;
  PC pc;
  Mat A;
  Vec x,b;

  PetscInitialize(&argc,&args,NULL,NULL);
  ierr = VecCreateSeq(PETSC_COMM_SELF,N,&x);
  ierr = VecCreateSeq(PETSC_COMM_SELF,N,&b);
  
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,N,N,m,
 			 PETSC_NULL,&A); 

  ierr = SLESCreate(PETSC_COMM_SELF,&sles); CHKERRA(ierr); 
  ierr = SLESSetOperators(sles,A,
			  A,SAME_NONZERO_PATTERN); CHKERRA(ierr); 
  ierr = SLESGetKSP(sles,&ksp); CHKERRA(ierr); 
  ierr = SLESGetPC(sles,&pc); CHKERRA(ierr); 

  ierr = KSPSetType(ksp,KSPPREONLY); CHKERRA(ierr); 
  ierr = PCSetType(pc,PCLU); CHKERRA(ierr); 
//    double pc_lu_fill=2.;
//    ierr = PCLUSetFill(pc,pc_lu_fill); CHKERRA(ierr); 
  ierr = PCLUSetUseInPlace(pc); CHKERRA(ierr);

  for (int j=0; j<N; j++) {
    v=1.;
    VecSetValue(b,j,v,ADD_VALUES);
    v=1.;
    MatSetValue(A,j,j,v,ADD_VALUES);
    v=.1;
    MatSetValue(A,j,(j+M)%N,v,ADD_VALUES);
    MatSetValue(A,j,(j+N-M)%N,v,ADD_VALUES);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);

  ierr = SLESSolve(sles,b,x,&its);
  ierr = VecGetArray(x,&a); CHKERRQ(ierr); 
  sum = 0.;
  for (int j=0; j<N; j++) sum = sum += a[j];
  printf("sum(x) = %f (expected %f)\n",sum,double(N)/1.2);

  PetscFinalize();
}
