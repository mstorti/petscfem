// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: pfmat.h,v 1.23 2001/10/06 23:37:08 mstorti Exp $
#ifndef PFMAT_H
#define PFMAT_H

#include <vector>

#include <src/distmap.h>
#include <src/distmat.h>

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


/// This is the basic distributed matrix class. 
//typedef DistMap<int,Row,IntPartitioner> DistMat;

/** This is a wrapper to the PETSc Matrix class and allows us to define
    new types
*/
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
  /// Defines the KSP method
  string KSP_method;
public:
  /// Constructor, initialize variables
  PFMat() : sles_was_built(0), A(NULL), P(NULL) {};
  /// Virtual destructor
  virtual ~PFMat();
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
  virtual int build_sles(TextHashTable *thash,char *name=NULL);
  /** Destroy the SLES associated with the operator. 
  */ 
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
  virtual int solve(Vec res,Vec dx);
  /// returns the number of iterations spent in the last solve
  virtual int its() {return its_;};
  /// Prints the matrix to a PETSc viewer
  virtual int view(Viewer viewer)=0;
  /// Derive this if you want to manage directly the preconditioning. 
  virtual int set_preco(const string & preco_type);
  /// Duplicate matrices (currently not implemented for IISDMat)
  virtual int duplicate(MatDuplicateOption op,const PFMat &A);
};

PFMat * PFMat_dispatch(const char *s);

/** Wrapper monitor. You customize the monitor by deriving the
    #monitor# function member. 
*/
int PFMat_default_monitor(KSP ksp,int n,double rnorm,void *);


#endif
