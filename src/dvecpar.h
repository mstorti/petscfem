// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvecpar.h,v 1.3 2006/04/08 20:37:15 mstorti Exp $
#ifndef PETSCFEM_DVECPAR_H
#define PETSCFEM_DVECPAR_H

#include <cassert>
#include <mpi.h>
#include <petsc.h>
#include <src/dvector.h>

extern int MY_RANK,SIZE;

template<class T>
class dvector_mpi_type {
public:
  MPI_Datatype type();
};

template<class T>
void 
dvector_clone_parallel(dvector<T> &w,int root=0) {
  int size = w.size();
  MPI_Datatype t = dvector_mpi_type<T>().type();
  int ierr = MPI_Bcast(&size,1,t,
		       root,PETSC_COMM_WORLD);
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
