// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfmat.h,v 1.28.2.4 2001/12/21 01:29:34 mstorti Exp $
#ifndef PFMAT_H
#define PFMAT_H

#include <vector>

#include <src/part.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/distmap.h>
#include <src/distmat.h>
#include <src/graph.h>

/// This partitioner is based on the dofmap of the mesh. 
class DofmapPartitioner : public IntRowPartitioner {
  /// Pointer to the dofmap 
  const Dofmap *dofmap;
public:
  /// Dof partitioning (currently based on ranges of the dof's).  
  int dofpart(int row);
  /// Constructor from the dofmap
  DofmapPartitioner(const Dofmap *dfm);
  /// Destructor
  ~DofmapPartitioner();
  /// interfaces with the DistMap class
  int processor(map<int,Row>::iterator k) {
    return dofpart(k->first);
  }
};

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
#ifdef PFFSM
  /// The Finite State Machine Object
  friend class pfmatFSM;
  pfmatFSM fsm;
#endif
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
  /// Defines the KSP method
  string KSP_method;
  /// The options database
  TextHashTable thash;
  /// These are the actions for the state machine
  int factored;
  virtual int factor_and_solve(Vec &res,Vec &dx)=0;
  virtual int solve_only(Vec &res,Vec &dx)=0;
  virtual int clean_factor() {assert(0);}; // fixme:= make it pure
  // virtual after
public:
  /// Constructor, initialize variables
  PFMat(DofPartitioner &pp);

  /// Virtual destructor
  virtual ~PFMat();

  /// clear memory (almost destructor)
  virtual void clear();

  /// calls MatAssemblyBegin on internal matrices, see PETSc doc
  virtual int assembly_begin(MatAssemblyType type)=0;

  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
  virtual int assembly_end(MatAssemblyType type)=0;

  /// Adds an element to the matrix profile
  set_profile(int j,int k);

  /** Creates the matrix from the profile graph computed in #g#
  */ 
  virtual void create();

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
  virtual int build_sles(TextHashTable *thash,char *name=NULL);

  /// Destroy the SLES associated with the operator. 
  virtual int destroy_sles();

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
  virtual int solve(Vec &res,Vec &dx);
  /// returns the number of iterations spent in the last solve
  virtual int its() {return its_;};
  /// Prints the matrix to a PETSc viewer
  virtual int view(Viewer viewer)=0;
  /// Derive this if you want to manage directly the preconditioning. 
  virtual int set_preco(const string & preco_type);
  /// Duplicate matrices (currently not implemented for IISDMat)
  virtual int duplicate(MatDuplicateOption op,const PFMat &A);
  void set_option(const char *key,const char *value) {
    thash.add_entry(key,value);
  }
  void get_option(const char *key,const char *&value) const {
    // Remove constness asuming this is OK.
    // (Should fix the `TextHashTable' class
    ((TextHashTable &)thash).get_entry(key,value);
  }
  int get_int(const char *name,int *retval,int defval=0,int n=1) {
    ::get_int(&thash,name,retval,defval,n);
  }
};

PFMat * PFMat_dispatch(const char *s);

/** Wrapper monitor. You customize the monitor by deriving the
    #monitor# function member. 
*/
int PFMat_default_monitor(KSP ksp,int n,double rnorm,void *);


#endif
