//__INSERT_LICENSE__
//$Id: pfmatsp.cpp,v 1.2 2001/10/06 23:37:08 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

#include <typeinfo>
#ifdef RH60
#include "libretto.h"
#else
#include <libretto/libretto.h>
#endif
#include <mat.h>

#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/iisdmat.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "SparseDirect::create"
void SparseDirect::create(Darray *da,const Dofmap *dofmap,
		int debug_compute_prof) {
  int size;
  MPI_Comm_size (PETSC_COMM_WORLD, &size);
  assert(size==1); // Only sequential use
  const int &neq = dofmap->neq;
  A.resize(neq,neq);

}

int SparseDirect::solve(Vec res,Vec dx) {
  double *res_p, *dx_p;
  int m,j,ierr;

  ierr = VecGetArray(res,&res_p); CHKERRQ(ierr); 
  ierr = VecGetArray(dx,&dx_p); CHKERRQ(ierr); 
  m = A.rows();
  for (j=0; j<m; j++) 
    dx_p[j] = res_p[j];
  A.solve(dx_p);
  ierr = VecRestoreArray(res,&res_p); CHKERRQ(ierr); 
  ierr = VecRestoreArray(dx,&dx_p); CHKERRQ(ierr); 
  return 0;
}

void SparseDirect::set_value(int row,int col,Scalar value,
			     InsertMode mode=ADD_VALUES) {
  double val = value;
  if (mode==ADD_VALUES)
    val += A.get(row,col);
  A.set(row,col,val);
}

int SparseDirect::duplicate(MatDuplicateOption op,const PFMat &B) {
  const SparseDirect * BB = dynamic_cast<const SparseDirect *>(&B);
  assert(BB);
  A = BB->A;
  return 0;
}
