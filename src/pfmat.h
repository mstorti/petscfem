// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfmat.h,v 1.28.2.10 2002/01/07 16:26:00 mstorti Exp $
#ifndef PFMAT_H
#define PFMAT_H

#include <vector>

#include <src/part.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/distmap.h>
#include <src/distmat.h>
#include <src/graph.h>

#define PF_ACTION_DECL(action) void action() 
#define PF_ACTION_DEF(action) \
        void pfmatFSMContext::action() { matrix_p->action(); }

#define PF_ACTION_LIST 				\
  PF_ACTION(create_a); 

#if 0
  PF_ACTION(factor_and_solve);			\
  PF_ACTION(solve_only);			\
  PF_ACTION(clean_factor);
#endif

class PFMat;

#define PF_ACTION(name) PF_ACTION_DECL(name) 
class pfmatFSMContext {
public:
  PFMat * matrix_p;
  // pfmatFSMContext(PFMat *p) : matrix_p(p) {};
  PF_ACTION_LIST;
  void FSMError(const char *e,const char *s) { 
    printf("PFMat: Not valid event \"%s\" in state \"%s\"\n",e,s);
  }
};
#undef PF_ACTION

#include "pfmatFSM.h"

class pfmatFSM;

/** This is a wrapper to the PETSc Matrix class and allows us to define
    new types
*/
class PFMat {
  friend class pfmatFSMContext;
protected:
  pfmatFSM fsm;

  ///@name Actions of the finite state machine.
  //@{
  /// The action corresponding to `set_profile'
  virtual int set_profile_a(int j,int k)=0;

  /// The action corresponding to `create'
  virtual void create_a()=0;

  /// The action corresponding to `set_value'
  virtual void set_value_a(int row,int col,Scalar value,
			   InsertMode mode=ADD_VALUES)=0;
  //@}

public:
  /// Returns the size of the j-th dimension
  virtual int size(int j)=0;

  /// Return both sizes
  void size(int &m, int&n) {
    m = size(1); n=size(2);
  }

  PFMat() { fsm.matrix_p = this; }

  /// Virtual destructor
  virtual ~PFMat()=0;

  /// Adds an element to the matrix profile
  int set_profile(int row,int col) {
    fsm.set_profile();
    set_profile_a(row,col);
  }

  /// Creates the matrix from the profile graph entered with `profile'
  void create() {
    fsm.create();
    create_a();
  }

  /** Sets individual values on the operator #A(row,col) = value#
      @param row (input) first index
      @param col (input) second index
      @param value (input) the value to be set
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
  */ 
  void set_value(int row,int col,Scalar value,
		 InsertMode mode=ADD_VALUES) {
    fsm.set_value();
    set_value_a(row,col, value, mode);
  }

  /// calls MatAssemblyBegin on internal matrices, see PETSc doc
  virtual int assembly_begin(MatAssemblyType type)=0;

  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
  virtual int assembly_end(MatAssemblyType type)=0;

  /// This calls both.
  int assembly(MatAssemblyType type) {
    int ierr;
    ierr = assembly_begin(type); CHKERRQ(ierr); 
    ierr = assembly_end(type); CHKERRQ(ierr); 
  }

  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  virtual int solve(Vec &res,Vec &dx)=0;

  /// Cleans the factored part
  virtual int clean_factor()=0;

  /// Sets all values of the operator to zero.
  virtual int zero_entries()=0;

  /// clear memory (almost destructor)
  virtual void clear()=0;

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
  virtual int view(Viewer viewer=VIEWER_STDOUT_WORLD)=0;

  /// Duplicate matrices (currently not implemented for IISDMat)
  virtual int duplicate(MatDuplicateOption op,const PFMat &A)=0;

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

  // virtual void print()=0;
};

#if 0
PFMat * PFMat_dispatch(const char *s);
#endif


#endif
