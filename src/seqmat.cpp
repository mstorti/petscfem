//__INSERT_LICENSE__
//$Id: seqmat.cpp,v 1.1 2001/09/19 21:39:04 mstorti Exp $

#include "seqmat.h"

void SeqMat::set_value(int row,int col,Scalar value,InsertMode mode=ADD_VALUES) {
  //  map 
  //store
}

void SeqMat::create(Darray *da,const Dofmap *dofmap_,int debug_compute_prof=0) {
  MPI_Comm_size (PETSC_COMM_WORLD, &size);
  assert(size==1);
}


/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/u/mstorti/PETSC/petscfem-beta-1.93/tools/pfcpp")
  End: 
*/
