//__INSERT_LICENSE__
//$Id: spdirect.cpp,v 1.6.74.1 2006/05/20 21:11:20 dalcinl Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

#include <petscmat.h>

#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>
#include <src/spdirect.h>
#include <src/iisdmat.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "SparseDirect::create_a"
int SparseDirect::create_a() {
  int size;
  MPI_Comm_size (PETSCFEM_COMM_WORLD, &size);
  assert(size==1); // Only sequential use
  return 0;
}

int SparseDirect::factor_and_solve_a(Vec &res,Vec &dx) {
  double *res_p, *dx_p;
  int m,j,ierr;

  ierr = VecGetArray(res,&res_p); CHKERRQ(ierr); 
  ierr = VecGetArray(dx,&dx_p); CHKERRQ(ierr); 
  m = A_p->rows();
  for (j=0; j<m; j++) 
    dx_p[j] = res_p[j];
  A_p->solve(dx_p);
  ierr = VecRestoreArray(res,&res_p); CHKERRQ(ierr); 
  ierr = VecRestoreArray(dx,&dx_p); CHKERRQ(ierr); 
  return 0;
}

int SparseDirect::solve_only_a(Vec &res,Vec &dx) {
  return factor_and_solve_a(res,dx);
}

int SparseDirect::set_value_a(int row,int col,PetscScalar value,
			     InsertMode mode) {
  double val = value;
  if (mode==ADD_VALUES)
    val += A_p->get(row,col);
  A_p->set(row,col,val);
  return 0;
}

int SparseDirect::duplicate_a(MatDuplicateOption op,const PFMat &B) {
  const SparseDirect * BB = dynamic_cast<const SparseDirect *>(&B);
  assert(BB);
#if 0
  *A_p = *(BB->A_p);
#else
  A_p->duplicate(*BB->A_p);
  if (op==MAT_DO_NOT_COPY_VALUES)
    A_p->clear();
  // Duplicate non-factored data only
  return 0;
#endif
}
