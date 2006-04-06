// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvecpar.h,v 1.1 2006/04/06 16:33:12 mstorti Exp $
#ifndef PETSCFEM_DVECPAR_H
#define PETSCFEM_DVECPAR_H

#include <src/dvector.h>
#include <mpi.h>

template<class T>
MPI_Datatype dvector_mpi_type();

template<class T>
void 
dvector_clone_parallel(dvector<T> &w,int root=0) {
  int size = w.size();
  int ierr = MPI_Bcast(&size,1,dvector_mpi_type<T>(),root,PETSC_COMM_WORLD);
  assert(!ierr);

  if (MY_RANK!=root) w.resize(size).defrag();
  ierr = MPI_Bcast(w.buff(), size, MPI_DOUBLE, 
		   root,PETSC_COMM_WORLD);
  assert(!ierr);
}

template<class T>
void 
dvector_read_parallel(const char *file,
		      dvector<T> &w,int root=0) {
  if (MY_RANK==root) {
    if (w.size()==0) {
      w.cat(file).defrag();
    } else {
      w.read(file).defrag();
    }
  }
  dvector_clone_parallel(w,root);
}

#endif
