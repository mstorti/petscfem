#ifndef PETSCFEM_CHIMERA_H
#define PETSCFEM_CHIMERA_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This is (let's say) wrapper hook to define nodes that are
// bdry (internal and external). We put this as a hook but in
// fact the only think that we use is that it is linked
// dinamically. But in fact we don't use at as hook (properly speaking). 
class chimera_hook_t : public Hook {
public:
  virtual void
  mark_bdry_nodes(set<int> &ebdry,set<int> &ibdry,
                  double time,int step) { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Auxiliary struct to store coefficients in the sparse
// interpolator. It represents the row and column indices
// J,K and the coefficient AJK.  Normally the interpolated
// callue for node j is phi{j} = sum{k} a{jk}*phi{k}
struct ajk_t {
  // Row and column indices
  int j,k;
  // Interpolator coefficient
  double ajk;
  // Utility constructor
  ajk_t(int j_,int k_,double ajk_)
    : j(j_), k(k_), ajk(ajk_) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*> This
// This class is the context for the PETSc MatShell object.
// When the MatShell object will call the MatMult operation
// we will pass it to this object. 
class chimera_mat_shell_t {
public:
  chimera_mat_shell_t() : A(NULL), xnod(NULL) {}

  // Initializes the problem
  int init(Mat A);
  // Do some tasks before starting the iterative solver
  int before_solve(Vec x,Vec res,double time,int step);
  // This is the matrix vector product and is called in each
  // iteration of the solver loop. It computes the product
  // of the FEM part (matrix A) and then the matrix-free
  // part due to the interpolation and boundary conditions.
  int mat_mult(Vec x,Vec y);
  int mat_mult_transpose(Vec x,Vec y);
  // The underlying PETSc matrix (assmbled by FEM)
  Mat A;
  // Pointer to internal PF array of node coordinates
  double *xnod;
  // Number of columns in XNOD
  int nu;
  // Returns component K of coordinates of node J
  // double xcoords(int j,int k) { return XNOD(j,k); }
  // For the internal boundaries (those that must be interpolated
  // from the other domain, maps the equation number to the
  // node number.
  // WARNING:= right now I think that in mny places we are
  // using that the equation number and the node are the
  // same. This restricts that the boundary conditions are
  // imposed in the Chimera module itself and no
  // parallelism.
  map<int,int> eq2node,node2eq;
  // Options for this Chimera module
  Json::Value opts;
  // Number of nodes of the W1 and W2 domains
  int nnod,nnod1,nnod2,nelem1,nelem2;
  // The interpolators in format (ROW,COL,COEF)
  vector<ajk_t> z;
  // The interpolators are ordered by ROW so that we can use
  // CSR pointers and then we store the beginning and end
  // for a certain ROW in the correspondings Z12PTR array.
  vector<int> zptr;
  // Problem specific: marks external bdry nodes (nodes at
  // bdries of the subdomains not at external bdries)
  // int isexternal(double *x);

  // Mark which nodes are imposed. They may be external
  // boundary nodes (ebdry) and internal bdry nodes (ibdry),
  // i.e. nodes that are inside the physical domain but in
  // the boundary to other overlapping region and so their
  // values must be interpolated from the other domain. 
  set<int> ibdry,ebdry;
};

// The PETSc MatMult and MatMultTranspose operations for
// Chimera MatShell object
int chimera_mat_mult(Mat Ashell,Vec x,Vec y);
int chimera_mat_mult_transpose(Mat Ashell,Vec x,Vec y);

// FIXME:= is it needed??
void init_hooks();

extern chimera_hook_t *CHIMERA_HOOK_P;

#endif
