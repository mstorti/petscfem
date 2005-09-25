// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfptscmat.h,v 1.16.70.1 2005/09/25 18:47:03 mstorti Exp $
#ifndef PFPTSCMAT_H
#define PFPTSCMAT_H

#include <vector>

#include <petsc.h>
#include <petscksp.h>

#include <src/iisdgraph.h>
#include <src/graphdv.h>
#include <src/linkgraph.h>
#include <src/pfmat.h>

class PFPETScMat : public PFMat {
protected:
  /// Number of rows/columns
  int mat_size;
  /// We will have always a PETSc matrix for the system and for the preconditioner
  Mat A,P;
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
  // int ksp_was_built; // now included in `factored'
  /// Defines the KSP method
  string KSP_method;
  /// The options database
  TextHashTable thash;
  /// These are the actions for the state machine
  int factored;
  // define subproblems in Additive Schwarz prec, default is 0.
  int asm_define_sub_problems;
  // preconditioning type on each sub-block, default is ilu(0)
  string asm_sub_preco_type;
  // number of blocks in ASM prec., default is one block per processor
  int asm_lblocks;
  // overlap between blocks in ASM prec., default is 1
  int asm_overlap;
  
  /// The partitioner object
  const DofPartitioner &part;

  /// The graph storing the profile object
  StoreGraph *lgraph;
  StoreGraph1 lgraph1;
  graphdv_dis lgraph_dv;
  LinkGraphWrapper lgraph_lkg;

  /// IntRowPartitioner based on a DofPartitioner
  class PFPETScPart : public IntRowPartitioner {
  public:
    const DofPartitioner &part;
    int processor(map<int,Row>::const_iterator k) {
      return part.processor(k->first);
    }
    PFPETScPart(const DofPartitioner & p) : part(p) {};
    ~PFPETScPart() {}
  } pf_part;

  /// The communicator
  MPI_Comm comm;

  /// Cleans linear system 'ksp'
  int clean_factor_a();

  /// Cleans linear matrix entries
  // int clean_mat_a();

  /// Cleans profile related stuff
  int clean_prof_a();

public:

  PFPETScMat(int MM,const DofPartitioner &pp,MPI_Comm comm_);

  ~PFPETScMat();

  int duplicate_a(MatDuplicateOption op,const PFMat &A);
  virtual int build_ksp();
  virtual int set_preco(const string & preco_type);
  int monitor(int n,double rnorm);

  /// Adds an element to the matrix profile
  int set_profile_a(int j,int k) {
    lgraph->add(j,k);
    return 0;
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
