// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfmat.h,v 1.28.2.12 2002/01/09 16:31:07 mstorti Exp $
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

#define PF_ACTION_DEF(action)			\
void pfmatFSMContext::action() {		\
  matrix_p->ierr = matrix_p->action();		\
  if (matrix_p->ierr)				\
    printf("pfmatFSMContext::action ierr=%d\n",	\
	   matrix_p->ierr);			\
}

#define PF_ACTION_LIST				\
  PF_ACTION(clean_prof_a);			\
  PF_ACTION(clean_mat_a);			\
  PF_ACTION(clean_factor_a);			\
  PF_ACTION(factor_and_solve_A);		\
  PF_ACTION(solve_only_A);

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
  /// Pointers to pass args to solver routines through the FSM layer
  Vec *res_p,*dx_p;
protected:
  /// Print Finite State Machine transitions
  int print_fsm_transition_info;
  /// Allows to pass PETSc error codes through the FSM layer
  int ierr;
  pfmatFSM fsm;

  ///@name Actions of the finite state machine.
  //@{
  /// The action corresponding to `set_profile'
  virtual int set_profile_a(int j,int k)=0;

  /// The action corresponding to `create'
  virtual int create_a()=0;

  /// The action corresponding to `set_value'
  virtual int set_value_a(int row,int col,Scalar value,
			   InsertMode mode=ADD_VALUES)=0;

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
  //@}

  /// Factorizes matrix and solves linear system. Args are passed via pointers.
  int factor_and_solve_A() { factor_and_solve_a(*res_p,*dx_p); }

  /// Solves linear system. Args are passed via pointers.
  int solve_only_A() { solve_only_a(*res_p,*dx_p); }

public:
  /// Returns the size of the j-th dimension
  virtual int size(int j)=0;

  /// Return both sizes
  void size(int &m, int&n) {
    m = size(1); n=size(2);
  }

  PFMat() : ierr(0), 
    print_fsm_transition_info(0) { fsm.matrix_p = this; }

  /// Virtual destructor
  virtual ~PFMat()=0;

  /// Adds an element to the matrix profile
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_profile"
  int set_profile(int row,int col) {
    fsm.set_profile();
    CHKERRQ(ierr); 

    ierr = set_profile_a(row,col);
    CHKERRQ(ierr); 
  }

  /// Creates the matrix from the profile graph entered with `profile'
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::create"
  int create() {
    fsm.create();
    CHKERRQ(ierr); 

    ierr = create_a();
    CHKERRQ(ierr); 
  }

  /** Sets individual values on the operator #A(row,col) = value#
      @param row (input) first index
      @param col (input) second index
      @param value (input) the value to be set
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
  */ 
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_value"
  int set_value(int row,int col,Scalar value,
		InsertMode mode=ADD_VALUES) {
    fsm.set_value();
    CHKERRQ(ierr); 

    ierr = set_value_a(row,col,value,mode); CHKERRQ(ierr); 
    return 0;
  }

  /// calls MatAssemblyBegin on internal matrices, see PETSc doc
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::assembly_begin"
  int assembly_begin(MatAssemblyType type) {
    fsm.assembly_begin();
    CHKERRQ(ierr); 

    ierr = assembly_begin_a(type); CHKERRQ(ierr);
    return 0;
  }

  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::assembly_end"
  int assembly_end(MatAssemblyType type) {
    fsm.assembly_end();
    CHKERRQ(ierr); 

    ierr = assembly_end_a(type); CHKERRQ(ierr); 
    return 0;
  }

  /// This calls both.
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::assembly"
  int assembly(MatAssemblyType type) {
    int ierr;
    ierr = assembly_begin(type); CHKERRQ(ierr); 
    ierr = assembly_end(type); CHKERRQ(ierr); 
    return 0;
  }

  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::solve"
  int solve(Vec &res,Vec &dx) {
    res_p = &res;
    dx_p = &dx;
    fsm.solve(); CHKERRQ(ierr); 
    return 0;
  }

  /// Factorizes matrix and solves linear system. Args are passed via pointers.
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::factor_and_solve"
  int factor_and_solve(Vec &res,Vec &dx) {
    res_p = &res;
    dx_p = &dx;
    fsm.factor_and_solve(); CHKERRQ(ierr); 
    return 0;
  }

  /// Solves (matrix should be factorized). Args are passed via pointers.
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::solve_only"
  int solve_only(Vec &res,Vec &dx) {
    res_p = &res;
    dx_p = &dx;
    fsm.solve_only(); CHKERRQ(ierr); 
    return 0;
  }

  /// Cleans the factored part
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clean_factor"
  int clean_factor() { 
    fsm.clean_factor(); CHKERRQ(ierr); 
    return 0;
  }

  /// Sets all values of the operator to zero.
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clean_mat"
  int clean_mat() { 
    fsm.clean_mat(); CHKERRQ(ierr); 
    return 0;
  }

  /// Cleans the profile related part
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clean_prof"
  int clean_prof() { 
    fsm.clean_prof(); CHKERRQ(ierr); 
    return 0;
  }

  /// clear memory (almost destructor)
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clear"
  int clear() { 
    fsm.clear(); CHKERRQ(ierr); 
    return 0;
  }

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

  int print_fsm_transition_info_f() { 
    return print_fsm_transition_info; };
};

#if 0
PFMat * PFMat_dispatch(const char *s);
#endif


#endif
