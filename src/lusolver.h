// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: lusolver.h,v 1.2 2002/09/05 18:34:18 mstorti Exp $
#ifndef LUSOLVER_H
#define LUSOLVER_H

// This is the OO wrapper to PETSc matrix
class LUPETScMat : public PETScMat {
  /// PETSc error code 
  int ierr;
public:
  /// Destructor (calls almost destructor)
  ~LUPETScMat() {clear();};
  /// Constructor 
  LUPETScMat() : PETScMat() {};
  /** Call the assembly begin function for the underlying PETSc matrices
      @param type (input) PETSc assembly type
      @return A PETSc error code 
  */ 
  void set_value(int row,int col,PetscScalar value,InsertMode mode=ADD_VALUES) {
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
  int view(PetscViewer viewer);
};


#endif
