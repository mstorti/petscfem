//__INSERT_LICENSE__
//$Id: laplace.cpp,v 1.8 2001/12/20 21:58:55 mstorti Exp $
 
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/arglist.h>
#include <src/utils.h>
#include <src/getprop.h>

#include "genload.h"
#include "lapla.h"
#include <time.h>

int MY_RANK,SIZE;

static char help[] = "Basic finite element program.\n\n";

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

int MyKSPMonitor(KSP ksp,int n,double rnorm,void *dummy)
{
  PetscPrintf(PETSC_COMM_WORLD,
	      "iteration %d KSP Residual norm %14.12e \n",n,rnorm);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
    SET_ELEMSET_TYPE(lapla)
    SET_ELEMSET_TYPE(genload)
	{
	printf("not known elemset \"type\": %s\n",type);
	exit(1);
	}
}

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {
  Vec     x, res;
  Mat     A;                          /* linear system matrix */
  SLES    sles;     /* linear solver context */
  PC      pc;           /* preconditioner context */
  KSP     ksp;        /* Krylov subspace method context */
  double  norm, *sol, scal; /* norm of solution error */
  int     flg, ierr,  size,  ndim, nel, nen, neq, myrank, its;
  double tol=2e-6;
  char fcase[FLEN+1];
  Dofmap *dofmap;
  Mesh *mesh;

  PetscInitialize(&argc,&args,(char *)0,help);
  print_copyright();

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  ierr = OptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  // Read the mesh
  read_mesh(mesh,fcase,dofmap,neq,size,myrank);
  dofmap->create_MPI_vector(x);

  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  string save_file = string("save_file.out");
  get_string(mesh->global_options,"save_file",save_file,1,1);
  printf("retornado por get_string: \"%s\"\n",save_file.c_str());

  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);
  ierr = VecSet(&scal,res); CHKERRA(ierr);

  ierr = opt_read_vector(mesh,x,dofmap,myrank);

  arg_list argl;

  VOID_IT(argl);
  argl.arg_add(&A,PROFILE);
  ierr = assemble(mesh,argl,dofmap,"comp_prof"); CHKERRA(ierr);

  ierr = MatZeroEntries(A); CHKERRA(ierr);
      
  // SLES solver for the laplacian
  ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRA(ierr);
  ierr = SLESSetOperators(sles,A,
			  A,DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
  ierr = SLESGetKSP(sles,&ksp); CHKERRA(ierr);
  ierr = SLESGetPC(sles,&pc); CHKERRA(ierr);

  ierr = KSPSetType(ksp,KSPCG); CHKERRA(ierr);
  ierr = PCSetType(pc,PCJACOBI); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,
         PETSC_DEFAULT); CHKERRA(ierr);
  ierr = KSPSetMonitor(ksp,MyKSPMonitor,PETSC_NULL);

#if 0
  // Computes matrix by finite differences
  VOID_IT(argl);
  argl.arg_add(&x,IN_VECTOR);
  argl.arg_add(&res,OUT_VECTOR);
  ierr = assemble(mesh,argl,dofmap,"comp_res"); CHKERRA(ierr);

  VOID_IT(argl);
  argl.arg_add(&x,PERT_VECTOR);
  argl.arg_add(&A,OUT_MATRIX_FDJ);
  ierr = assemble(mesh,argl,dofmap,"comp_res"); CHKERRA(ierr);
#else
  // Computes analytic matrix 
  VOID_IT(argl);
  argl.arg_add(&x,IN_VECTOR);
  argl.arg_add(&res,OUT_VECTOR);
  argl.arg_add(&A,OUT_MATRIX);
  ierr = assemble(mesh,argl,dofmap,"comp_res_mat"); CHKERRA(ierr);
#endif

  ierr = SLESSolve(sles,res,x,&its); CHKERRA(ierr); 
  print_vector(save_file.c_str(),x,dofmap);

  PetscFinalize();
  exit(0);
}

