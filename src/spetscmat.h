// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: spetscmat.h,v 1.1 2004/10/24 16:25:21 mstorti Exp $
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

class PFSymmPETScMat : public PETScMat {

protected:
  Mat Asymm;

public:

  PFSymmPETScMat(int MM,int NN,const DofPartitioner &part_a,
		 MPI_Comm comm_a = PETSC_COMM_WORLD) 
    : PETScMat(MM,NN,part_a,comm_a) {}

  int mult(Vec x,Vec y);

  /// Adds an element to the matrix profile
  int set_profile_a(int j,int k) {
    if (k>=j) lgraph->add(j,k);
    return 0;
  }

  int set_value_a(int row,int col,PetscScalar value,
		  InsertMode mode=ADD_VALUES) {
    int roww=row, coll=col;
    if (col<row) { roww=col; coll=row; }
    ierr = MatSetValues(A,1,&row,1,&col,&value,mode); 
    CHKERRQ(ierr); 
    return 0;
  };

  int duplicate_a(MatDuplicateOption op,const PFMat &A) { assert(0); }

  int build_sles();

};

#endif
