// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: iisdmat.h,v 1.7 2001/10/06 23:37:08 mstorti Exp $
#ifndef IISDMAT_H
#define IISDMAT_H

#include <vector>

#include <src/pfmat.h>
#include <src/sparse.h>

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
  /// Duplicate matrices 
  int duplicate(MatDuplicateOption op,const PFMat &A);
  int view(Viewer viewer);
};

int iisd_jacobi_pc_apply(void *ctx,Vec,Vec);

/** Solves iteratively on the `interface' (between subdomain) nodes
    and solving by a direct method in the internal nodes.
 */
class IISDMat : public PFMat {
  /** Type of nodes, L: local, I: internal.
      Type of block (PETSc sense) D: diagonal, O: off-diagonal. 
   */
  static const int D,O,L,I;
  /// The dofmap of the related mesh (this should be changed)!!
  const Dofmap *dofmap;
  /// Number of interface nodes 
  int n_int;
  /// Number of local nodes
  int n_loc;
  /// Dof's in this processor begin here (#k1<= k <=k2#)
  int k1;
  /// Dof's in this processor end here (#k1<= k <=k2#)
  int k2;
  /// Number of interface dof's in all processors
  int n_int_tot;
  /// Number of local dof's in all processors
  int n_loc_tot;
  /** Local dof's in this processor start at this position in the
      MPI local vector.
   */
  int n_locp;
  /** Interface nodes in this processor start at this position in the
      MPI local vector.
   */
  int n_intp;
  /// Number of dof's in this processor. 
  int neqp;
  /// Maps old numbering in new numbering 
  vector<int> map;
  /** Local dof's in this processor are in the range
      #n_loc_v[myrank] <= dof < #n_loc_v[myrank+1]
   */
  vector<int> n_loc_v;
  /** Interface dof's in this processor are in the range
      #n_loc_v[myrank] <= dof < #n_loc_v[myrank+1]
   */
  vector<int> n_int_v;
  /// The PETSc nnz vector for the local part
  vector<int> d_nnz_LL;
  /// Version SuperLU of the local matrix
  Sparse::Mat A_LL_SLU;
  /// Local-Local matrix (sequential matrix on each processor). 
  Mat A_LL;
  /// Local-Interface matrix (MPI matrix).
  Mat A_LI;
  /// Interface-Local matrix (MPI matrix).
  Mat A_IL;
  /// Interface-Interface matrix (MPI matrix).
  Mat A_II;
  /// Shortcuts to the #A_LL#, #A_IL#, #A_LI# and #A_II# matrices. 
  Mat *AA[2][2];
  /// Partitioner
  DofmapPartitioner *part;
  /** Here we put all non-local things that are in the loca-local
      block on other processors
  */
  DistMatrix *A_LL_other;
  /// The mode we are inserting values
  InsertMode insert_mode;
  /// Auxiliar MPI vector that contains all local dof's
  Vec x_loc;
  /** Auxiliar sequential vector that contains local dof's in this
      processor
  */
  Vec x_loc_seq;
  /** Auxiliar sequential vector that contains local dof's in this
      processor
  */
  Vec y_loc_seq;
  /// SLES for local solution (en each processor)
  SLES sles_ll;
  /// PC for local solution (en each processor)
  PC pc_ll;
  /// KSP for local solution (en each processor)
  KSP ksp_ll;
  /// Warn if not appropriate setting for preconditioning type
  static int warn_iisdmat;
  /** For a global dof #gdof# gives the #block# (`local' or
      `interface') and the number in that block. 
      @param gdof (input) dof in old numbering
      @param block (output) the block to which the dof belongs 
      @param ldof (input) number of dof in new ordering relative to
      its block 
  */ 
  void map_dof(int gdof,int &block,int &ldof);
  /** Solves the local problem #A_LL x_loc = s * y_loc# for #x_loc#. 
      #x_loc# and #y_loc# may be aliased. 
      @param x_loc (output) the solution vector
      @param y_loc (input) the right hand side
      @param trans (input) if non null solves with the traspose of #A_LL_# 
      @param s (input) scale factor (usually for changing sign)
  */ 
  int local_solve(Vec x_loc,Vec y_loc,int trans=0,double s=1.);
  /** Idem to #local_solve# but for SuperLU (Sparse::Mat)
      representation of the local problem.
  */
  int local_solve_SLU(Vec x_loc,Vec y_loc,int trans=0,double c=1.);
  /// Diagonal of Interface matrix to use as preconditioning
  Vec A_II_diag;
  /// PETSc LU fill parameter 
  double pc_lu_fill;
  /// Layers of nodes of the preconditioning
  vector< set<int> > int_layers;

public:

  /// Local solver type
  enum LocalSolver {PETSc, SuperLU} local_solver;

  /** Creates the matrix from the profile computed in #da#
      @param da (input) dynamic array containing the adjacency matrix
      of the operator
      @param dofmap (input) the dofmap of the operator (contains
      information about range of dofs per processor. 
      @param debug_compute_prof (input) flag for debugging the process
      of building the operator.
  */ 
  void create(Darray *da,const Dofmap *dofmap_,
	      int debug_compute_prof=0);
  /** Applies the Schur operator #y = S * x#
      @param x (input) a given interface vector
      @param y (output) the result of applying the Schur operator on
      #x#. Usually involves a solution (by a direct method) on each
      subdomain. 
  */ 
  virtual int mult(Vec x,Vec y);
  /** Applies the traspose of the Schur operator (see #mult#). 
      @param x (input) a given interface vector
      @param y (output) the result of applying the Schur operator on
      #x#. 
  */ 
  virtual int mult_trans(Vec x,Vec y);
  /** Sets individual values on the operator #A(row,col) = value#
      @param row (input) first index
      @param col (input) second index
      @param value (input) the value to be set
      @param mode (input) either #ADD_VALUES# (default) or #INSERT_VALUES#
  */ 
  void set_value(int row,int col,Scalar value,
		 InsertMode mode=ADD_VALUES);
  /// Clear the object (almost destructor)
  void clear();
  /// Sets the underlying matrices to zero
  int zero_entries();
  /// Calls MatAssemblyBegin on internal matrices, see PETSc doc
  int assembly_begin(MatAssemblyType type);
  /// calls MatAssemblyEnd on internal matrices, see PETSc doc
  int assembly_end(MatAssemblyType type);
  /// Prints the matrix to a PETSc viewer
  int view(Viewer viewer);
  /** Performs all operations needed before permorming the solution of
      the linear system (creating the PETSc SLES, etc...). Sets
      oprtions from the #thash# table. 
      @param name (input) a string to be prepended to all options that
      apply to the operator. (not implemented yet)
  */ 
  int build_sles(TextHashTable *thash,char *name=NULL);
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int solve(Vec res,Vec dx);
  /// Derive this if you want to manage directly the preconditioning. 
  int set_preco(const string & preco_type);
  /// Destroy the SLES associated with the operator. 
  virtual int destroy_sles();
  IISDMat() : A_LL_other(NULL), A_LL(NULL), part(NULL),
  local_solver(PETSc) {};
  /// The PETSc wrapper function calls this
  int jacobi_pc_apply(Vec x,Vec y); 
  /// Destructor
  ~IISDMat();
};

/// Direct SuperLU solver. 
class SparseDirect : public PFMat {
public:
  Sparse::Mat A;
  /// destructor
  ~SparseDirect() {A.clear();};
  /// clear memory (almost destructor)
  void clear() {A.clear();};
  /// does nothing here (only sequential use...)
  int assembly_begin(MatAssemblyType type) { return 0;};
  /// does nothing here (only sequential use...)
  int assembly_end(MatAssemblyType type) { return 0;};
  /// Resizes the underlying #A# matrix. 
  void create(Darray *da,const Dofmap *dofmap_,
		      int debug_compute_prof=0);
  void set_value(int row,int col,Scalar value,
		 InsertMode mode=ADD_VALUES);
  /// Sets all values of the operator to zero.
  int zero_entries() {A.clear(); return 0;};
  // Does nothing
  int build_sles(TextHashTable *thash,char *name=NULL) {return 0;};
  // Does nothing
  int destroy_sles() {return 0;};
  // solve system calling Sparse::Mat::solve()
  int solve(Vec res,Vec dx);
  /// returns the number of iterations spent in the last solve
  int its() {return 1;};
  /// Prints the matrix to a PETSc viewer
  int view(Viewer viewer=NULL) {A.print(); return 0;};
  /// Derive this if you want to manage directly the preconditioning. 
  int set_preco(const string & preco_type) {return 0;};
  /// Duplicate matrices (currently not implemented for IISDMat)
  int duplicate(MatDuplicateOption op,const PFMat &B);
};

#endif
