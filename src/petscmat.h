// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: petscmat.h,v 1.1.2.1 2001/12/27 19:55:47 mstorti Exp $
#ifndef PETSCMAT_H
#define PETSCMAT_H

#include <vector>

#include <src/pfmat.h>
#include <src/pfptscmat.h>

// This is the OO wrapper to PETSc matrix
class PETScMat : public PFMat {
  int M,N;
  /// PETSc error code 
  int ierr;
  /// The graph storing the profile object
  StoreGraph lgraph;
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int factor_and_solve(Vec &res,Vec &dx);
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int solve_only(Vec &res,Vec &dx);
public:
  /// Destructor (calls almost destructor)
  ~PETScMat() {clear();};
  /// clear memory (almost destructor)
  void clear();

  PETScMat(int MM,int NN,const DofPartitioner &pp,MPI_Comm comm_ =
	  PETSC_COMM_WORLD) : 
    PFPETScMat(comm_), 
    M(MM), N(NN), 
    part(pp), lgraph(M,&part,comm), 
    A_LL_other(NULL), A_LL(NULL), 
    local_solver(PETSc), pf_part(pp) {};


  // returns the j-th dimension
  int size(int j) { 
    assert(j==1 || j==2);
    return (j==1 ? M : N);
  }

  /// Constructor 
  PETScMat() : PFMat() {};
  /// Adds an element to the matrix profile
  int set_profile(int j,int k) {
    lgraph.add(j,k);
  }

  /** Call the assembly begin function for the underlying PETSc matrices
      @param type (input) PETSc assembly type
      @return A PETSc error code 
  */ 
  int assembly_begin(MatAssemblyType type) {
    return MatAssemblyBegin(A,type);};
  /** Call the assembly end function for the underlying PETSc matrices
      @param type (input) PETSc assembly type
      @return A PETSc error code 
  */ 
  int assembly_end(MatAssemblyType type) {
    return MatAssemblyEnd(A,type);};
  /** Sets individual values on the operator #A(row,col) = value#
      @param row (input) first index
      @param col (input) second index
      @param value (input) the value to be set
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
  */ 
  void set_value(int row,int col,Scalar value,InsertMode mode=ADD_VALUES) {
    MatSetValues(A,1,&row,1,&col,&value,mode);};
  /// Sets all values of the operator to zero.
  int zero_entries();
  /** Creates the matrix from the profile computed in #da#
      @param da (input) dynamic array containing the adjacency matrix
      of the operator
      @param dofmap (input) the dofmap of the operator (contains
      information about range of dofs per processor. 
      @param debug_compute_prof (input) flag for debugging the process
      of building the operator.
  */ 
  void create(Darray *da,const Dofmap *dofmap_,int debug_compute_prof=0);
  /// Duplicate matrices 
  int duplicate(MatDuplicateOption op,const PFMat &A);
  int view(Viewer viewer);
};

#endif
