#include <sles.h>

static char help[] = "Basic finite element program.\n\n";

void petscfem_initialize(void) {
#define ARGC 1
  char **args;
  int argc;
  //  args = new (char *)[1];
  args = (char **)malloc(sizeof(char *));
  args[0] = "/home/mstorti/PETSC/petscfem/run/../applications/advective/adv.bin";
  printf("hi in petscfem_initialize\n");
  argc=ARGC;
  PetscInitialize(&argc,&args,(char *)0,help);
}  

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
int make_rhs_vector(Vec *b) {
  int n,ierr,j;
  printf("hi in make_rhs_vector\n");
  b = (Vec *)malloc(sizeof(Vec));
  n=10;

  ierr =  VecCreateSeq(PETSC_COMM_SELF,n,b); CHKERRA(ierr); 

  for (j=0; j<n; j++) {
    double val;
    val=(double)(j+1);
    VecSetValue(*b,j,val,INSERT_VALUES);
  }
  ierr = VecAssemblyBegin(*b); CHKERRA(ierr);
  ierr = VecAssemblyEnd(*b); CHKERRA(ierr);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void print_vector(Vec *x) {
  printf("hi in print_vector\n");
  VecView(*x,VIEWER_STDOUT_WORLD);
}
