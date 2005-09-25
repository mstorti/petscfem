//__INSERT_LICENSE__
//$Id: tryme4.cpp,v 1.6.80.1 2005/09/25 22:58:58 mstorti Exp $
#include <petscksp.h>

int null_monitor(KSP ksp,int n, double rnorm,void *A_) {return 0;}

int main(int argc,char **args) {
  const int N=100000, M=100;
  int ierr,m=3,its, Nit=10;
  double v, *a, *a2, sum,sum2,dit;
  KSP ksp;
  PC pc;
  Mat A;
  Vec x,b,x2,b2;

  PetscInitialize(&argc,&args,NULL,NULL);
  ierr = VecCreateSeq(PETSC_COMM_SELF,N,&x);
  ierr = VecCreateSeq(PETSC_COMM_SELF,N,&b);
  ierr = VecCreateSeq(PETSC_COMM_SELF,N,&b2);
  ierr = VecCreateSeq(PETSC_COMM_SELF,N,&x2);

  dit=0.;
  while (1) {
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,N,N,m,
			   PETSC_NULL,&A); 
    ierr = MatZeroEntries(A); CHKERRA(ierr);

    ierr = KSPCreate(PETSC_COMM_SELF,&ksp); CHKERRA(ierr); 
    ierr = KSPSetOperators(ksp,A,
			    A,SAME_NONZERO_PATTERN); CHKERRA(ierr); 
    ierr = KSPGetPC(ksp,&pc); CHKERRA(ierr); 

#if 0
    ierr = KSPSetType(ksp,KSPPREONLY); CHKERRA(ierr); 
    ierr = PCSetType(pc,PCLU); CHKERRA(ierr); 
    double pc_lu_fill=1.7;
    ierr = PCLUSetFill(pc,pc_lu_fill); CHKERRA(ierr); 
//      ierr = PCLUSetUseInPlace(pc); CHKERRA(ierr);
#else
    ierr = KSPSetType(ksp,KSPGMRES); CHKERRA(ierr); 
    if (KSP_method=="gmres") {
      ierr = KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);
      CHKERRQ(ierr);
    }
    ierr = KSPSetTolerances(ksp,0,0,1e10,1); CHKERRA(ierr); 
    ierr = PCSetType(pc,PCLU); CHKERRA(ierr); 
#endif
    
    for (int j=0; j<N; j++) {
      v=dit;
      VecSetValue(b,j,v,ADD_VALUES);
      v=2.*dit;
      VecSetValue(b2,j,v,ADD_VALUES);
      v=1.;
      MatSetValue(A,j,j,v,ADD_VALUES);
      v=.1;
      MatSetValue(A,j,(j+M)%N,v,ADD_VALUES);
      MatSetValue(A,j,(j+N-M)%N,v,ADD_VALUES);
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);

    ierr = KSPSolve(ksp,b,x); CHKERRA(ierr);
    //    ierr = KSPGetIterationNumber(ksp,&its)

    ierr = KSPSolve(ksp,b2,x2); CHKERRA(ierr); 
    //    ierr = KSPGetIterationNumber(ksp,&its)

    ierr = VecGetArray(x,&a); CHKERRA(ierr); 
    ierr = VecGetArray(x2,&a2); CHKERRA(ierr); 
    sum = 0.;
    sum2 = 0.;
    for (int j=0; j<N; j++) {
      sum = sum += a[j];
      sum2 = sum2 += a2[j];
    }
    if (int(dit) % 100 ==0) {
      printf("sum(x) = %f (expected %f)\n",sum,dit*double(N)/1.2);
      printf("sum(x2) = %f (expected %f)\n",sum2,2*dit*double(N)/1.2);
    }
    ierr = VecRestoreArray(x,&a);
    ierr = VecRestoreArray(x2,&a2);
    ierr = KSPDestroy(ksp); CHKERRA(ierr); 
    ierr = MatDestroy(A); CHKERRA(ierr); 
    dit +=1.;
    if(Nit>0 && dit > double(Nit)) break;
  }
  PetscFinalize();
}
