// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvecpar2.h,v 1.2 2006/04/08 21:37:35 mstorti Exp $
#ifndef PETSCFEM_DVECPAR2_H
#define PETSCFEM_DVECPAR2_H

#include <cassert>
#include <mpi.h>
#include <petsc.h>

#include <src/dvector.h>
#include <src/debug.h>

extern int MY_RANK,SIZE;

#define MPI_TYPE_UNDEFINED INT_MAX

template<class T>
class dvector_mpi_type {
public:
  MPI_Datatype type() { return MPI_TYPE_UNDEFINED; }
};

template<class T>
void 
dvector_clone_parallel(dvector<T> &w,int root=0) {
  int size = w.size();
  int ierr = MPI_Bcast(&size,1,MPI_INT,
		       root,PETSC_COMM_WORLD);
  assert(!ierr);
  
  if (MY_RANK!=root) w.resize(size);
  w.defrag();
  MPI_Datatype t = dvector_mpi_type<T>().type();
  if (t!=MPI_TYPE_UNDEFINED) {
    ierr = MPI_Bcast(w.buff(), size, MPI_DOUBLE, 
		     root,PETSC_COMM_WORLD);
    assert(!ierr);
  } else {
    ierr = MPI_Bcast(w.buff(), size*sizeof(T), MPI_CHAR, 
		     root,PETSC_COMM_WORLD);
    assert(!ierr);
  }
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
