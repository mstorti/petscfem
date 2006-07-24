// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfmat.h,v 1.40 2006/07/24 04:28:15 mstorti Exp $
#ifndef PFMAT_H
#define PFMAT_H

#include <vector>

#include <src/part.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/distmap.h>
#include <src/distmat.h>
#include <src/graph.h>

class pfmatFSM;

/// This is the generic matrix class
class PFMat {
  friend class pfmatFSMContext;
  /// Pointers to pass args to solver routines through the FSM layer
  Vec *res_p,*dx_p;

protected:
  /// Print Finite State Machine transitions
  int print_fsm_transition_info;

  /// Allows to pass PETSc error codes through the FSM layer
  int ierr;
  pfmatFSM *fsm;

  ///@name Actions of the finite state machine.
  //@{
  /// The action corresponding to `set_profile'
  virtual int set_profile_a(int j,int k)=0;

  /// The action corresponding to `create'
  virtual int create_a()=0;

  /// The action corresponding to `set_value'
  virtual int set_value_a(int row,int col,PetscScalar value,
			   InsertMode mode=ADD_VALUES)=0;

  /// The action corresponding to `set_values'
  virtual int set_values_a(int nrows,int *idxr,int ncols,int *idxc,
			   PetscScalar *values, InsertMode mode=ADD_VALUES);

  /// calls MatAssemblyBegin on internal matrices, see PETSc doc
  virtual int assembly_begin_a(MatAssemblyType type)=0;

  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
  virtual int assembly_end_a(MatAssemblyType type)=0;

  /// Factorizes matrix and solves linear system. Args are passed via pointers.
  virtual int factor_and_solve_a(Vec &res,Vec &dx)=0;

  /// Solves (matrix should be factorized). Args are passed via pointers.
  virtual int solve_only_a(Vec &res,Vec &dx)=0;

  /// Cleans the factored part
  virtual int clean_factor_a()=0;

  /// Sets all values of the operator to zero.
  virtual int clean_mat_a()=0;

  /// clear profile memory 
  virtual int clean_prof_a()=0;

  /// duplicate matrix
  virtual int duplicate_a(MatDuplicateOption op,const PFMat &A) {
    printf("Not implemented yet\n"); assert(0);
  }
  //@}

  /// Factorizes matrix and solves linear system. Args are passed via pointers.
  int factor_and_solve_A() { factor_and_solve_a(*res_p,*dx_p); }

  /// Solves linear system. Args are passed via pointers.
  int solve_only_A() { solve_only_a(*res_p,*dx_p); }

public:
  /// This is the `factory' of `PFMat' matrices.
  static PFMat *dispatch(int N,DofPartitioner &part,const char *s);

  /// Returns the size of the j-th dimension
  virtual int size(int j)=0;

  /// Return both sizes
  void size(int &m, int&n) {
    m = size(1); n=size(2);
  }

  PFMat();

  /// Virtual destructor
  virtual ~PFMat();

  /// Adds an element to the matrix profile
  int set_profile(int row,int col);

  /// Creates the matrix from the profile graph entered with `profile'
  int create();

  /** Sets individual values on the operator #A(row,col) = value#
      @param row (input) first index
      @param col (input) second index
      @param value (input) the value to be set
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
  */ 
  int set_value(int row,int col,PetscScalar value,
		InsertMode mode=ADD_VALUES);

  /** Sets an array of values on the operator . For
      #mode=INSERT_VALUES# it is equivalent to #for (j=0; j<nrows;
      j++) for (k=0; k<ncols; k++) A(idxr[j],idxc[k]) =
      values[j*ncols+k]#. For #mode=ADD_VALUES# the assignment
      operator is replaced by #+=#.
      @param nrows (input) number of row indices
      @param idxr (input) array of row indices of size #nrows#
      @param ncols (input) number of column indices
      @param idxc (input) array of column indices of size #ncols#
      @param values (input) array of values of size #nrows*ncols# 
      stored by row (C style). 
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES# */
  int set_values(int nrows,int *idxr,int ncols,int *idxc,
		 PetscScalar *values, InsertMode mode=ADD_VALUES);

  /// calls MatAssemblyBegin on internal matrices, see PETSc doc
  int assembly_begin(MatAssemblyType type);

  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
  int assembly_end(MatAssemblyType type);

  /// This calls both.
  int assembly(MatAssemblyType type);

  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int solve(Vec &res,Vec &dx);

  /// Factorizes matrix and solves linear system. Args are passed via pointers.
  int factor_and_solve(Vec &res,Vec &dx);

  /// Solves (matrix should be factorized). Args are passed via pointers.
  int solve_only(Vec &res,Vec &dx);

  /// Cleans the factored part
  int clean_factor();

  /// Sets all values of the operator to zero.
  int clean_mat();

  /// Cleans the profile related part
  int clean_prof();

  /// clear memory (almost destructor)
  int clear();

  /** Defines how to report convergence in the internal loop. 
      Derive this function to obtain a different effect from the
      default one. 
      @param n (input) current iteration number
      @param rnorm (input) norm of the residual for the current
      iteration 
      @return A PETSc error code 
  */ 
  virtual int monitor(int n,double rnorm) {}

  /// returns the number of iterations spent in the last solve
  virtual int its() {return 0;};

  /// Prints the matrix to a PETSc viewer
  virtual int view(PetscViewer viewer=PETSC_VIEWER_STDOUT_WORLD)=0;

  /// Duplicate matrices (currently not implemented for IISDMat)
  int duplicate(MatDuplicateOption op,const PFMat &A);

  virtual void set_option(const char *key,const char *value)=0;
  virtual void set_option(const char *name,int *val,int n=1)=0;
  void set_option(const char *name,int val) {
    int val_ = val;
    return set_option(name,&val_);
  }
  virtual void set_option(const char *name,double *val,int n=1)=0;
  void set_option(const char *name,double val) {
    double val_ = val;
    return set_option(name,&val_);
  }

  virtual void get_option(const char *key,const char *&value) const=0;
  virtual void get_option(const char *name,int *retval,int defval=0,int
			 n=1)=0;
  virtual void get_option(const char *name,double *retval,int defval=0,int
			 n=1)=0;
  virtual void get_option(const char *name,string &retval,int defval=0,int
			 n=1 )=0;

  int print_fsm_transition_info_f() { 
    return print_fsm_transition_info; };
};

#endif
