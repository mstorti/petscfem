// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfmat.h,v 1.28.2.7 2001/12/27 10:12:37 mstorti Exp $
#ifndef PFMAT_H
#define PFMAT_H

#include <vector>

#include <src/part.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/distmap.h>
#include <src/distmat.h>
#include <src/graph.h>

//#define PFFSM
#ifdef PFFSM
#undef FSM_ACTION
#define FSM_ACTION(action) void pfmatFSMContext::action()	\
  {matrix_p->action();}

class pfmatFSMContext {
public:
  PFMat * matrix_p;
  pfmatFSMContext() {};
  FSM_ACTION(factor_and_solve);
  FSM_ACTION(solve_only);
  FSM_ACTION(clean_factor);
  void FSMError(const char *e,const char *s) { 
    printf("Not valid \"%s\" event in state \"%s\"",e,s);
  }
};

#include "pfmatFSM.h"
#endif

/** This is a wrapper to the PETSc Matrix class and allows us to define
    new types
*/
class PFMat {
protected:
  // virtual after
public:
  /// Returns the size of the j-th dimension
  virtual int size(int j)=0;

  /// Return both sizes
  void size(int &m, int&n) {
    m = size(1); n=size(2);
  }

  PFMat() ;
#if 0
  /// Constructor, initialize variables
  PFMat(int m, int n,DofPartitioner &pp);
#endif

  /// Virtual destructor
  virtual ~PFMat()=0;

  /// clear memory (almost destructor)
  virtual void clear()=0;

  /// calls MatAssemblyBegin on internal matrices, see PETSc doc
  virtual int assembly_begin(MatAssemblyType type)=0;

  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
  virtual int assembly_end(MatAssemblyType type)=0;

  /// Adds an element to the matrix profile
  virtual int set_profile(int j,int k)=0;

  /** Creates the matrix from the profile graph computed in #g#
  */ 
  virtual void create()=0;

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

  /// Destroy the SLES associated with the operator. 
  virtual int clean_factor()=0;

  /** Defines how to report convergence in the internal loop. 
      Derive this function to obtain a different effect from the
      default one. 
      @param n (input) current iteration number
      @param rnorm (input) norm of the residual for the current
      iteration 
      @return A PETSc error code 
  */ 
  virtual int monitor(int n,double rnorm) {}
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  virtual int solve(Vec &res,Vec &dx)=0;
  /// returns the number of iterations spent in the last solve
  virtual int its() {return 0;};
  /// Prints the matrix to a PETSc viewer
  virtual int view(Viewer viewer)=0;
  /// Derive this if you want to manage directly the preconditioning. 
  virtual int set_preco(const string & preco_type)=0;
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

  virtual void print()=0;
};

#if 0
PFMat * PFMat_dispatch(const char *s);
#endif


#endif
