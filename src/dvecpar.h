// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvecpar.h,v 1.4 2006/04/08 21:35:17 mstorti Exp $
#ifndef PETSCFEM_DVECPAR_H
#define PETSCFEM_DVECPAR_H

#include <cassert>
#include <mpi.h>
#include <petsc.h>
#include <src/dvector.h>


/** Clones a dvector object from the #root# processor to all
    the others (like #MPI_Bcast()#). 
    @param w (input/output) the dvector to be broadcasted. As 
    a side effect it is defragmented in the root procesor and
    is created as defragmented in all others.
    @param root (input) the processor from which the dvector 
    is boradcasted */ 
template<class T>
void 
dvector_clone_parallel(dvector<T> &w,int root=0);

/** Reads a vector in the #root# processor and broadcasts it to 
    all the others. The final vectors are defragmented. 
    @param file (input) The name of the file
    @param w (input/output) the dvector to be broadcasted. As 
    a side effect it is defragmented in the root procesor and
    is created as defragmented in all others.
    @param root (input) the processor from which the dvector 
    is boradcasted */ 
template<class T>
void 
dvector_read_parallel(const char *file,
		      dvector<T> &w,int root=0);

#endif
