// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvecpar2.h,v 1.2.20.1 2007/02/19 20:23:56 mstorti Exp $
#ifndef PETSCFEM_DVECPAR2_H
#define PETSCFEM_DVECPAR2_H

#include <cassert>
#include <mpi.h>
#include <petsc.h>

#include <src/dvector.h>
#include <src/debug.h>

extern int MY_RANK,SIZE;

#define MPI_TYPE_UNDEFINED INT_MAX

// In some installations, depending on the sepcifics of MPI
// it seems that MPI_DATATYPE is a pointer and then it
// complains when compared to an INT_MAX. I tried to use the
// following define and obvioulsy it complains but I didn't confirm
// if it really works. 
// #define MPI_TYPE_UNDEFINED NULL

template<class T>
class dvector_mpi_type {
public:
  MPI_Datatype type() { return MPI_TYPE_UNDEFINED; }
};

template<class T>
void 
dvector_clone_parallel(dvector<T> &w,int root) {
  int size = w.size();
  int PFUNUSED ierr = MPI_Bcast(&size,1,MPI_INT,
		       root,PETSCFEM_COMM_WORLD);
  assert(!ierr);
  
  if (MY_RANK!=root) w.resize(size);
  vector<int> shape;
  int rank;
  if (MY_RANK==root) {
    w.get_shape(shape);
    rank = shape.size();
  }
  ierr = MPI_Bcast(&rank,1,MPI_INT,root,PETSCFEM_COMM_WORLD);
  if (MY_RANK!=root) shape.resize(rank);
  ierr = MPI_Bcast(&shape[0],rank,MPI_INT,root,PETSCFEM_COMM_WORLD);
  
  if (MY_RANK!=root) w.reshape(shape);
  w.defrag();
  MPI_Datatype t = dvector_mpi_type<T>().type();
  if (t!=MPI_TYPE_UNDEFINED) {
    ierr = MPI_Bcast(w.buff(), size, MPI_DOUBLE, 
		     root,PETSCFEM_COMM_WORLD);
    assert(!ierr);
  } else {
    ierr = MPI_Bcast(w.buff(), size*sizeof(T), MPI_CHAR, 
		     root,PETSCFEM_COMM_WORLD);
    assert(!ierr);
  }
}

template<class T>
void 
dvector_read_parallel(const char *file,
		      dvector<T> &w,int root) {
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
