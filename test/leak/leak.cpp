/*__INSERT_LICENSE__*/
// $Id: leak.cpp,v 1.3 2002/10/06 21:51:35 mstorti Exp $

#define _GNU_SOURCE

// Tests for the `PFMat' class
#ifdef USE_OLD_PETSC_VERSION
#include <sles.h>
#else
#include <petscsles.h>
#endif

extern int SIZE;

inline double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

void print_memory_usage(void) {
  int ierr;
  char *file,*line;
  int mem,mem_min,mem_max,mem_sum,mem_avrg;
  size_t n=0;
  assert(asprintf(&file,"/proc/%d/status",getpid()));
  FILE *fid = fopen(file,"r");
  while(1) {
    assert(getline(&line,&n,fid)!=-1);
    if (sscanf(line,"VmRSS: %d kB",&mem)) break;
  }
  fclose(fid);
  MPI_Allreduce(&mem,&mem_min,1,MPI_INT,MPI_MIN,PETSC_COMM_WORLD);
  MPI_Allreduce(&mem,&mem_max,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
  MPI_Allreduce(&mem,&mem_sum,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
  mem_avrg = int(ceil(double(mem_sum)/double(SIZE)));
  PetscPrintf(PETSC_COMM_WORLD,
	      "[Memory usage(kB): min %d, max %d, avrg %d]\n",
	      mem_min,mem_max,mem_avrg);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {

  PetscInitialize(&argc,&args,NULL,NULL);
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);

  Mat A;
  Vec b,x;
  SLES sles;
  PC pc;
  KSP ksp;
  Debug debug2;
  debug2.activate("memory_usage");

  const int N=100000;
  double h = 1./double(N);
  
  int niter=200;

  for (int iter=0; iter<=niter; iter++) {
    // debug2.trace("enter iter");
    print_memory_usage();
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

    for (int k=0; k<N; k++) {
      MatSetValue(A,k,k,1.,INSERT_VALUES); 
      
      int l=k+1;
      if (l>=N) l-=N;
      MatSetValue(A,k,l,-0.1,INSERT_VALUES); 
      
      l=k-1;
      if (l<0) l+=N;
      MatSetValue(A,k,l,-0.1,INSERT_VALUES); 
      VecSetValue(b,k,2.*drand()-1,INSERT_VALUES); 
    }
    ierr = VecAssemblyBegin(b);
    ierr = VecAssemblyEnd(b);
    
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
    
    int its;
    ierr = SLESSolve(sles,b,x,&its); CHKERRQ(ierr); 

    PetscPrintf(PETSC_COMM_WORLD,"iter %d, converged on %d CG iters\n",
		iter,its);
#if 0
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); 
    ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); 
#endif

    ierr = MatDestroy(A);
    ierr = VecDestroy(x);
    ierr = VecDestroy(b);
    ierr = SLESDestroy(sles);

  }
  PetscFinalize();
}
