// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: spetscmat.h,v 1.4 2004/10/25 02:09:32 mstorti Exp $
#ifndef PETSCFEM_SPETSCMAT_H
#define PETSCFEM_SPETSCMAT_H

// Symmetric PETSc Matrix only stores the upper diagonal part
#include <vector>

#include <petsc.h>
#include <petscsles.h>

#include <src/iisdgraph.h>
#include <src/graphdv.h>
#include <src/linkgraph.h>
#include <src/pfmat.h>

class PETScSymmMat : public PETScMat {

protected:
  Mat Asymm;
  int myrank, size;

  int insert_p(int row,int col);
  int insert_p2(int row,int col);

public:

  PETScSymmMat(int MM,int NN,const DofPartitioner &part_a,
		 MPI_Comm comm_a = PETSC_COMM_WORLD);

  int mult(Vec x,Vec y);

  /// Adds an element to the matrix profile
  int set_profile_a(int j,int k);

  int set_value_a(int row,int col,PetscScalar value,
		  InsertMode mode=ADD_VALUES);

  int duplicate_a(MatDuplicateOption op,const PFMat &A) { assert(0); }

  int build_sles();

};

#endif
