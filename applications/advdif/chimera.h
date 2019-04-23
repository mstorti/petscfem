#ifndef PETSCFEM_CHIMERA_H
#define PETSCFEM_CHIMERA_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Marks the nodes that are bdry (internal and external)
class chimera_hook_t : public Hook {
public:
  virtual void mark_bdry_nodes(set<int> &ebdry,set<int> &ibdry) { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
struct ajk_t {
  int j,k;
  double ajk;
  ajk_t(int j_,int k_,double ajk_)
    : j(j_), k(k_), ajk(ajk_) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class chimera_mat_shell_t {
public:
  chimera_mat_shell_t() : A(NULL), xnod(NULL) {}

  // Initializes the problem
  int init0(Mat A);
  // Initializes the problem before starting the iterative
  // solver
  int init1(Vec x,Vec res);
  // This is called in each iteration of the solver iteration loop
  int mat_mult(Vec x,Vec y);
  // The underlying PETSc matrix
  Mat A;
  // Pointer to internal PF array of node coordinates
  double *xnod;
  int nu;
  // Returns component K of coordinates of node J
  // double xcoords(int j,int k) { return XNOD(j,k); }
  // For the internal boundaries (those that must be interpolated
  // from the other domain, maps the equation number to the
  // node number
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
  // Read integer array from H5 double dataset
  void h5d2i(const char *file,const char *dset,dvector<int> &w);
  // Problem specific: marks external bdry nodes (nodes at
  // bdries of the subdomains not at external bdries)
  int isexternal(double *x);
#if 0  
  // Read the bdry from H5 and separate external and internal nodes
  void readbdry(const char *file,const char *dset,set<int> &ebdry,
                set<int> &ibdry,int domain);
#endif
  // Mark which nodes are imposed. They may be external
  // boundary nodes (ebdry) and internal bdry nodes (ibdry),
  // i.e. nodes that are inside the physical domain but in
  // the boundary to other overlapping region and so their
  // values must be interpolated from the other domain. 
  // This stores the rows that are fixed
  vector<int> rows;
  set<int> ibdry,ebdry;
};

extern chimera_hook_t *CHIMERA_HOOK_P;

#endif
