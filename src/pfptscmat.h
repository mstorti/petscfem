// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfptscmat.h,v 1.1.2.5 2001/12/27 21:35:27 mstorti Exp $
#ifndef PFPTSCMAT_H
#define PFPTSCMAT_H

#include <vector>

#include <petsc.h>
#include <sles.h>

#include <src/iisdgraph.h>
#include <src/pfmat.h>

class PFPETScMat : public PFMat {
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
  /// Defines the KSP method
  string KSP_method;
  /// The options database
  TextHashTable thash;
  /// These are the actions for the state machine
  int factored;

  /// The partitioner object
  const DofPartitioner &part;

  /// The graph storing the profile object
  StoreGraph lgraph;

  /// IntRowPartitioner based on a DofPartitioner
  class PFPETScPart : public IntRowPartitioner {
  public:
    const DofPartitioner &part;
    int processor(map<int,Row>::iterator k) {
      return part.processor(k->first);
    }
    PFPETScPart(const DofPartitioner & p) : part(p) {};
    ~PFPETScPart() {}
  } pf_part;

  /// The communicator
  MPI_Comm comm;

  virtual int factor_and_solve(Vec &res,Vec &dx)=0;
  virtual int solve_only(Vec &res,Vec &dx)=0;
  int clean_factor();

public:

  int solve(Vec &res,Vec &dx);

  PFPETScMat(int MM,const DofPartitioner &pp,MPI_Comm comm_) 
    : comm(comm_), part(pp), pf_part(part), lgraph(MM,&part,comm_), 
    sles_was_built(0), 
    A(NULL), P(NULL), factored(0) {}

  ~PFPETScMat();
  void clear();
  int duplicate(MatDuplicateOption op,const PFMat &A);
  virtual int build_sles();
  int set_preco(const string & preco_type);
  int monitor(int n,double rnorm);

  /// Adds an element to the matrix profile
  int set_profile(int j,int k) {
    lgraph.add(j,k);
  }

  /// returns the number of iterations spent in the last solve
  int its() {return its_;};

  void set_option(const char *key,const char *value) {
    thash.add_entry(key,value);
  }
  void get_option(const char *key,const char *&value) const {
    // Remove constness asuming this is OK.
    // (Should fix the `TextHashTable' class
    ((TextHashTable &)thash).get_entry(key,value);
  }
  void get_option(const char *name,int *retval,int defval=0,int n=1) {
    ::get_int(&thash,name,retval,defval,n);
  }
  void get_option(const char *name,double *retval,int defval=0,int n=1) {
    ::get_double(&thash,name,retval,defval,n);
  }
  void get_option(const char *name,string &retval,int defval=0,int
		  n=1 ) {
    get_string(&thash,name,retval,defval,n); }


  /// Not defined yet
  void set_option(const char *name,int *val,int n=1) {
    thash.add_entry(name,val,n); }
  void set_option(const char *name,double *val,int n=1) {
    thash.add_entry(name,val,n); }

};

/** Wrapper monitor. You customize the monitor by deriving the
    #monitor# function member. 
*/
int PFPETScMat_default_monitor(KSP ksp,int n,double rnorm,void *);

#endif
