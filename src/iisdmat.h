// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: iisdmat.h,v 1.14.2.11 2001/12/30 20:00:25 mstorti Exp $
#ifndef IISDMAT_H
#define IISDMAT_H

#include <vector>

#include <src/pfmat.h>
#include <src/sparse.h>
#include <src/sparse2.h>
#include <src/iisdgraph.h>
#include <src/pfptscmat.h>

#define TGETOPTDEF_ND_PFMAT(thash,type,name,default)		\
        name = default;						\
        ierr = ::get_##type(thash,#name,&name,1);		\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

#define TGETOPTDEF_S_ND_PFMAT(thash,type,name,default)		\
        name=type(#default);					\
        ierr = ::get_##type(thash,#name,name,1);		\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

#define TGETOPTDEF_S_PFMAT(thash,type,name,default)		\
        type name=type(#default);				\
        ierr = get_##type(thash,#name,name,1);			\
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

int iisd_jacobi_pc_apply(void *ctx,Vec,Vec);

/** Solves iteratively on the `interface' (between subdomain) nodes
    and solving by a direct method in the internal nodes.
 */
class IISDMat : public PFPETScMat {
  int M,N;
  /** Type of dofs: L: local, I: interface
      Type of block (PETSc sense) D: diagonal, O: off-diagonal. 
  */
  static const int D,O,L,I;
  /// Number of interface nodes 
  int n_int;
  /// Number of local nodes
  int n_loc;
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
  /** The number of subdomains in which local nodes to each processor
      are subdivided 
  */
  int iisd_subpart;
  /// Maps old numbering in new numbering 
  vector<int> map_dof;
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
  /// Version  of the local matrix
  Sparse::SuperLUMat A_LL_SLU;
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
  void map_dof_fun(int gdof,int &block,int &ldof);
  /** Solves the local problem #A_LL x_loc = s * y_loc# for #x_loc#. 
      #x_loc# and #y_loc# may be aliased. 
      @param x_loc (output) the solution vector
      @param y_loc (input) the right hand side
      @param trans (input) if non null solves with the traspose of #A_LL_# 
      @param s (input) scale factor (usually for changing sign)
  */ 
  int local_solve(Vec x_loc,Vec y_loc,int trans=0,double s=1.);
  /** Idem to #local_solve# but for SuperLU (Sparse::Mat)
g      representation of the local problem.
  */
  int local_solve_SLU(Vec x_loc,Vec y_loc,int trans=0,double c=1.);
  /// Diagonal of Interface matrix to use as preconditioning
  Vec A_II_diag;
  /// PETSc LU fill parameter 
  double pc_lu_fill;
  /// Layers of nodes of the preconditioning
  vector< set<int> > int_layers;
  /** Performs all operations needed before permorming the solution of
      the linear system (creating the PETSc SLES, etc...). 
  */ 
  //  int build_sles();
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int factor_and_solve(Vec &res,Vec &dx);
  /** Solve the linear system 
      @param res (input) the rhs vector
      @param dx (input) the solution vector
  */ 
  int solve_only(Vec &res,Vec &dx);

  /// Clean all data related to factorization
  int clean_factor();

  /// Maps dofs in this processors to global dofs
  vector<int> dofs_proc;
  /// Poitner to the storage area in `dofs_proc'
  int *dofs_proc_v;

  /** Maps dof's in this processor to global
      ones. The inverse of `dofs_proc'. 
  */
  map<int,int> proc2glob;

public:

  // returns the j-th dimension
  int size(int j) { 
    assert(j==1 || j==2);
    return (j==1 ? M : N);
  }

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
  void create();

  /** Applies the Schur operator #y = S * x#
      @param x (input) a given interface vector
      @param y (output) the result of applying the Schur operator on
      #x#. Usually involves a solution (by a direct method) on each
      subdomain. 
  */ 
  int mult(Vec x,Vec y);
  /** Applies the traspose of the Schur operator (see #mult#). 
      @param x (input) a given interface vector
      @param y (output) the result of applying the Schur operator on
      #x#. 
  */ 
  int mult_trans(Vec x,Vec y);
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
  int view(Viewer viewer=VIEWER_STDOUT_WORLD);
  /// Derive this if you want to manage directly the preconditioning. 
  int set_preco(const string & preco_type);

  // void print(void);
  /// Constructor
  IISDMat(int MM,int NN,const DofPartitioner &pp,MPI_Comm comm_ =
	  PETSC_COMM_WORLD) : 
    PFPETScMat(MM,pp,comm_), 
    M(MM), N(NN), 
    A_LL_other(NULL), A_LL(NULL), 
    local_solver(PETSc), sles_ll(NULL) {};
  /// The PETSc wrapper function calls this
  int jacobi_pc_apply(Vec x,Vec y); 
  /// Destructor
  ~IISDMat();
};

#endif
