// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfmat.h,v 1.2 2001/07/16 14:14:16 mstorti Exp $
#ifndef PFMAT_H
#define PFMAT_H

#include <vector>

// This is a wrapper to the PETSc Matrix class and allows us to define
// new types
class PFMat {
protected:
  /// We will have always a PETSc matrix for the system and for the preconditioner
  Mat A,P;
  /// The PETSc linear solver
  SLES sles;
  /// The PETSc preconditioner
  PC pc;
  /// Krylov subspace method context
  KSP ksp;
  /// PETSc absolute tolerance
  double atol;
  /// PETSc relative tolerance
  double rtol;
  /// PETSc divergence tolerance
  double dtol;
  /// Krylov space dimension
  int Krylov_dim;
  /// Maximum number of iterations
  int maxits;
  /// flags whether the internal convergence is print or not
  int print_internal_loop_conv;
  /// number of iterations done
  int its_;
  /// flags whether the system was built or not
  int sles_was_built;
public:
  /// Constructor, initialize variables
  PFMat() : sles_was_built(0), A(NULL), P(NULL) {};
  /// Virtual destructor
  virtual ~PFMat()=0;
  /// clear memory (almost destructor)
  virtual void clear();
  /// calls MatAssemblyBegin on internal matrices, see PETSc doc
  virtual int assembly_begin(MatAssemblyType type)=0;
  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
  virtual int assembly_end(MatAssemblyType type)=0;
  /** Creates the matrix from the profile computed in #da#
      @param da (input) dynamic array containing the adjacency matrix
      of the operator
      @param dofmap (input) the dofmap of the operator (contains
      information about range of dofs per processor. 
      @param debug_compute_prof (input) flag for debugging the process
      of building the operator.
  */ 
  virtual void create(Darray *da,const Dofmap *dofmap_,
		      int debug_compute_prof=0)=0;
  /** Sets individual values on the operator #A(row,col) = value#
      @param row (input) first index
      @param col (input) second index
      @param value (input) the value to be set
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
  */ 
  virtual void set_value(int row,int col,Scalar value,
			 InsertMode mode=ADD_VALUES)=0;
  /// Sets all values of the operator to zero.
  virtual int zero_entries()=0;
  /** Performs all operations needed before permorming the solution of
      the linear system (creating the PETSc SLES, etc...). Sets
      oprtions from the #thash# table. 
      @param name (input) a string to be prepended to all options that
      apply to the operator. (not implemented yet)
  */ 
  virtual void build_sles(TextHashTable *thash,char *name=NULL);
  /** Defines how to report convergence in the internal loop. 
      Derive this function to obtain a different effect from the
      default one. 
      @param ksp (input) the KSP of the linear solver
      @param n (input) current iteration number
      @param rnorm (input) norm of the residual for the current
      iteration 
      @return A PETSc error code 
  */ 
  virtual int monitor(KSP ksp,int n,double rnorm);
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  virtual void solve(Vec res,Vec dx);
  /// returns the number of iterations spent in the last solve
  virtual int its() {return its_;};
};

/** Wrapper monitor. You customize the monitor by deriving the
    #monitor# function member. 
*/
int PFMat_default_monitor(KSP ksp,int n,double rnorm,void *);

// This is the OO wrapper to PETSc matrix
class PETScMat : public PFMat {
  /// PETSc error code 
  int ierr;
public:
  /// Destructor (calls almost destructor)
  ~PETScMat() {clear();};
  /// clear memory (almost destructor)
  void clear();
  /// Constructor 
  PETScMat() : PFMat() {};
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
  int zero_entries() {ierr=MatZeroEntries(A); CHKERRQ(ierr);};
  /** Creates the matrix from the profile computed in #da#
      @param da (input) dynamic array containing the adjacency matrix
      of the operator
      @param dofmap (input) the dofmap of the operator (contains
      information about range of dofs per processor. 
      @param debug_compute_prof (input) flag for debugging the process
      of building the operator.
  */ 
  void create(Darray *da,const Dofmap *dofmap_,int debug_compute_prof=0);
};

#if 0
class LUsubdMat : public PFMat {
  const Dofmap *dofmap;
  int n_int,n_loc;
  vector<int> map;
  Mat A_LL,A_LI,A_IL,A_II;
public:
  void create(Darray *da,const Dofmap *dofmap_,
	      int debug_compute_prof=0);
  void clear();
};
#endif

#endif
