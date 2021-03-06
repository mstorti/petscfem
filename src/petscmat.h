// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: petscmat.h,v 1.6.72.1 2007/02/19 20:23:56 mstorti Exp $
#ifndef PETSCMAT_H
#define PETSCMAT_H

#include <vector>

#include <src/fem.h>
#include <src/pfmat.h>
#include <src/pfptscmat.h>

// This is the OO wrapper to PETSc matrix
class PETScMat : public PFPETScMat {
protected:
  int M,N;
  int ierr;
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int factor_and_solve_a(Vec &res,Vec &dx);
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int solve_only_a(Vec &res,Vec &dx);

  /// Maps dofs in this processors to global dofs
  vector<int> dofs_proc;
  /// Poitner to the storage area in `dofs_proc'
  int *dofs_proc_v;

  /** Maps dof's in this processor to global
      ones. The inverse of `dofs_proc'. 
  */
  map<int,int> proc2glob;

public:
  /// Destructor (calls almost destructor)
  ~PETScMat();

  PETScMat(int MM,int NN,const DofPartitioner &pp,MPI_Comm comm_ =
	  PETSCFEM_COMM_WORLD) : 
    PFPETScMat(MM,pp,comm_), 
    M(MM), N(NN) {};

  // returns the j-th dimension
  int size(int j) { 
    assert(j==1 || j==2);
    return (j==1 ? M : N);
  }
  /** Call the assembly begin function for the underlying PETSc matrices
      @param type (input) PETSc assembly type
      @return A PETSc error code 
  */ 
  int assembly_begin_a(MatAssemblyType type) {
    return MatAssemblyBegin(A,type);};
  /** Call the assembly end function for the underlying PETSc matrices
      @param type (input) PETSc assembly type
      @return A PETSc error code 
  */ 
  int assembly_end_a(MatAssemblyType type) {
    return MatAssemblyEnd(A,type);};
  /** Sets individual values on the operator #A(row,col) = value#
      @param row (input) first index
      @param col (input) second index
      @param value (input) the value to be set
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
  */ 
  int set_value_a(int row,int col,PetscScalar value,InsertMode mode=ADD_VALUES) {
    ierr = MatSetValues(A,1,&row,1,&col,&value,mode); CHKERRQ(ierr); 
    return 0;
  };
  /// Sets all values of the operator to zero.
  int clean_mat_a();
  /** Creates the matrix from the profile computed in #da#
      @param da (input) dynamic array containing the adjacency matrix
      of the operator
      @param dofmap (input) the dofmap of the operator (contains
      information about range of dofs per processor. 
      @param debug_compute_prof (input) flag for debugging the process
      of building the operator.
  */ 
  int create_a();
  /// Duplicate matrices 
  int duplicate_a(MatDuplicateOption op,const PFMat &A);
  int view(PetscViewer viewer=PETSC_VIEWER_STDOUT_WORLD);
};

#endif
