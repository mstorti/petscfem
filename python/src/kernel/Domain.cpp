// $Id$

#include <sstream>

#include "Domain.h"

#include <fem.h>
#include <dofmap.h>

PF4PY_NAMESPACE_BEGIN

Domain::~Domain()
{ }

Domain::Domain(const Domain& domain)
  : comm(domain.comm), 
    ndim(domain.ndim), 
    nnod(domain.nnod), 
    ndof(domain.ndof),
    mesh(domain.mesh),
    dofset(domain.dofset),
    appctx(domain.appctx)
{ }

Domain::Domain(int ndim, int nnod, int ndof)
  : comm(PETSC_COMM_WORLD),
    ndim(ndim), 
    nnod(nnod), 
    ndof(ndof)
{
  if (comm==MPI_COMM_NULL) throw Error("Domain: null communicator");
  if (ndim<1)              throw Error("Domain: ndim < 1");
  if (ndim>3)              throw Error("Domain: ndim > 3");
  if (nnod<0)              throw Error("Domain: nnod < 0");
  if (ndof<1)              throw Error("Domain: ndof < 1");
  
  this->mesh   = new Mesh   (this->comm, this->ndim, this->nnod);
  this->dofset = new Dofset (this->comm, this->nnod, this->ndof);
  this->appctx = new AppNS  ();
}

Domain::Domain(int ndim, int nnod, int ndof, MPI_Comm comm)
  : comm(comm), 
    ndim(ndim), 
    nnod(nnod), 
    ndof(ndof)
{
  if (comm==MPI_COMM_NULL) throw Error("Domain: null communicator");
  if (ndim<1)              throw Error("Domain: ndim < 1");
  if (ndim>3)              throw Error("Domain: ndim > 3");
  if (nnod<0)              throw Error("Domain: nnod < 0");
  if (ndof<1)              throw Error("Domain: ndof < 1");
  
  this->mesh   = new Mesh   (this->comm, this->ndim, this->nnod);
  this->dofset = new Dofset (this->comm, this->nnod, this->ndof);
  this->appctx = new AppNS  ();
}

Options&
Domain::getOptions() const
{
  return this->mesh->getOptions();
}

void
Domain::setOptions(const Options& options)
{
  this->mesh->setOptions(options);
}

DTable<double>& 
Domain::getNodedata() const
{
  return this->mesh->getNodedata();
}

void
Domain::setNodedata(const DTable<double>& nodetable)
{
  const std::pair<int,int>& shape = nodetable.getShape();
  if (this->nnod != shape.first)
    throw Error("Domain: rows != nnod in nodedata");
  if (this->ndim > shape.second)
    throw Error("Domain: cols < ndim in nodedata");
  this->mesh->setNodedata(nodetable);
}

DTable<double>& 
Domain::getField(const std::string& name) const
{
  return this->mesh->getField(name);
}

void  
Domain::setField(const std::string& name,
		 DTable<double>& data)
{
  const std::pair<int,int>& shape = data.getShape();
  if (this->nnod != shape.first)
    throw Error("Domain: rows != nnod in nodedata");
  this->mesh->setField(name,data);
}

Elemset& 
Domain::getElemset(int index) const
{
  return this->mesh->getElemset(index);
}

void
Domain::addElemset(const Elemset& elemset)
{
  return this->mesh->addElemset(elemset);
}

inline void 
check_nodes(const int nn, const int nodes[], const int nnod) {
  for(int n=0; n<nn; ++n) {
    if (nodes[n] < 0)
      throw Error("invalid node, out of range (node<0)");
    if (nodes[n] >= nnod)
      throw Error("invalid node, out of range (node>=nnod)");
  }
}
inline void 
check_fields(const int nf, const int fields[], const int ndof) {
  for(int f=0; f<nf; ++f) {
    if (fields[f] < 0)
      throw Error("invalid field, out of range (field<0)");
    if (fields[f] >= ndof) 
      throw Error("invalid field, out of range (field>=ndof)");
  }
}
inline void 
check_values(const int nv, const double values[]) {
  for(int v=0; v<nv; ++v) {
    if (values[v] != values[v])
      throw Error("invalid value entry, not a number (NaN)");
  }
}
inline void 
check_coeffs(const int nc, const double coeffs[]) {
  for(int c=0; c<nc; ++c) {
    //if (coeffs[c] == 0)
    //  throw Error("invalid coeff entry, coeff is 0.0");
    if (coeffs[c] != coeffs[c])
      throw Error("invalid coeff entry, not a number (NaN)");
  }
}

void 
Domain::setFixation(int nn, const int    nodes[],
		    int nf, const int    fields[],
		    int nv, const double values[]) 
{
  std::vector<int>    _f;
  std::vector<double> _v;

  if (nn == 0 || nf == 0) return;

  check_nodes (nn, nodes,  this->nnod);
  check_fields(nf, fields, this->ndof);
  check_values(nv, values);
  if (nn > 1 && nf == 1) { _f.resize(nn, fields[0]); nf = nn; fields = &_f[0]; }
  if (nn != nf) throw Error("incompatible array sizes, size(nodes) != size(fields)");
  if (nn > 1 && nv == 1) { _v.resize(nn, values[0]); nv = nn; values = &_v[0]; }
  if (nn != nv) throw Error("incompatible array sizes, size(nodes) != size(values)");

  this->dofset->add_fixations(nn, nodes, fields, values);
}

void 
Domain::setFixation(int nn, const int nodes[],
		    int nf, const int fields[],
		    Amplitude& amplitude) {
  std::vector<int>    _f;
  std::vector<double> _v;

  if (nn == 0 || nf == 0) return;
  check_nodes (nn, nodes,  this->nnod);
  check_fields(nf, fields, this->ndof);
  if (nn > 1 && nf == 1) { _f.resize(nn, fields[0]); nf = nn; fields = &_f[0]; }
  if (nn != nf) throw Error("incompatible array sizes, size(nodes) != size(fields)");

  this->dofset->add_fixations(nn, nodes, fields, NULL, &amplitude);
}


void 
Domain::setFixation(int nn, const int    nodes[],
		    int nf, const int    fields[],
		    int nv, const double values[],
		    Amplitude& amplitude) 
{
  std::vector<int>    _f;
  std::vector<double> _v;

  if (nn == 0 || nf == 0) return;

  check_nodes (nn, nodes,  this->nnod);
  check_fields(nf, fields, this->ndof);
  check_values(nv, values);
  if (nn > 1 && nf == 1) { _f.resize(nn, fields[0]); nf = nn; fields = &_f[0]; }
  if (nn != nf) throw Error("incompatible array sizes, size(nodes) != size(fields)");
  if (nn > 1 && nv == 1) { _v.resize(nn, values[0]); nv = nn; values = &_v[0]; }
  if (nn != nv) throw Error("incompatible array sizes, size(nodes) != size(values)");

  this->dofset->add_fixations(nn, nodes, fields, values, &amplitude);
}

void 
Domain::setPeriodic(int n1, const int nodes1[],
		    int n2, const int nodes2[]) 
{
  double C[2] = {1.0, -1.0};
  int    N[2];
  int    F[2];

  if ((n1 == 0 && n2 == 0)) return;
  
  if (n1 != n2) throw Error("incompatible array sizes, size(nodes1) != size(nodes2)");
  check_nodes(n1, nodes1, this->nnod);
  check_nodes(n2, nodes2, this->nnod);
  const int ndof = this->ndof;
  for(int n=0; n<n1; ++n) {
    N[0] = nodes1[n];
    N[1] = nodes2[n];
    for(int f=0; f<ndof; ++f) {
      F[0] = F[1] = f;
      this->dofset->add_constraints(2, C, N, F);
    }
  }
}

void 
Domain::setPeriodic(int n1, const int nodes1[],
		    int n2, const int nodes2[],
		    int nf, const int fields[]) {
  double C[2] = {1.0, -1.0};
  int    N[2];
  int    F[2];

  if ((n1 == 0 && n2 == 0) || nf == 0) return;
  
  if (n1 != n2) throw Error("incompatible array sizes, size(nodes1) != size(nodes2)");

  check_nodes (n1, nodes1, this->nnod);
  check_nodes (n2, nodes2, this->nnod);
  check_fields(nf, fields, this->ndof);

  for(int n=0; n<n1; n++) {
    N[0] = nodes1[n];
    N[1] = nodes2[n];
    for(int f=0; f<nf; f++) {
      F[0] = F[1] = fields[f]; 
      this->dofset->add_constraints(2, C, N, F);
    }
  }
}

void 
Domain::setConstraint(int nc, const double coeffs[],
		      int nn, const int    nodes[],
		      int nf, const int    fields[])
{
  if (nn == 0 && nf == 0 && nc == 0) return;
//   if (nc < 2)   throw Error("array size too small, size(coeffs) < 2");
//   if (nn < 2)   throw Error("array size too small, size(nodes) < 2 ");
//   if (nf < 2)   throw Error("array size too small, size(fields) < 2");
//   if (nn != nf) throw Error("incompatible array sizes, size(nodes) != size(fields)");
//   if (nn != nc) throw Error("incompatible array sizes, size(nodes) != size(coeffs)");
//   check_nodes (nn, nodes,  this->nnod);
//   check_fields(nf, fields, this->ndof);
//   check_coeffs(nc, coeffs);

//   this->dofset->add_constraints(nc, coeffs, nodes, fields);
  if (nc == nn * nf) {
    std::vector<int> tmp(nf); int *N = &tmp[0];
    for(int n=0; n<nn; n++) {
      for(int i=0; i<nf; i++) N[i] = nodes[n];
      this->dofset->add_constraints(nf, coeffs+nf*n, N, fields);
    }
  } else { 
    if (nc < 2)   throw Error("array size too small, size(coeffs) < 2");
    if (nn < 2)   throw Error("array size too small, size(nodes) < 2 ");
    if (nf < 2)   throw Error("array size too small, size(fields) < 2");
    if (nn != nf) throw Error("incompatible array sizes, size(nodes) != size(fields)");
    if (nn != nc) throw Error("incompatible array sizes, size(nodes) != size(coeffs)");
    this->dofset->add_constraints(nc, coeffs, nodes, fields);
  }
}

void 
Domain::setUp() 
{
  srand(1); // XXX Errors in sucessive setups (mesh partitioning) !!!
  this->mesh->setup(this);
  this->dofset->setup(this);
  this->appctx->setup(this);
}

Mesh&
Domain::getMesh() const
{
  if (!this->mesh) throw Error("Domain: mesh not created");
  return this->mesh;
}

Dofset&
Domain::getDofset() const
{
  if (!this->dofset) throw Error("Domain: dofset not created");
  return this->dofset;
}

AppCtx&
Domain::getAppCtx() const
{
  if (!this->appctx) throw Error("Domain: appctx not created");
  return this->appctx;
}

static void 
mk_vec(Vec& o, Comm comm, const std::pair<int,int>& sizes,
       const std::string& type, const std::string& name)
{
  PetscTruth valid;
  PetscInt n, N;
  const char *otype = (type.size() > 0) ? type.c_str() : NULL;
  if (o == PETSC_NULL) {
    try {
      // create object
      n = sizes.first; N = sizes.second;
      PF4PY_PETSC_CALL(VecCreate, (comm, &o));
      PF4PY_PETSC_CALL(VecSetSizes, (o, n, N));
      // set object type
      if (!otype) otype = (comm.getSize() == 1) ? VECSEQ : VECMPI;
      PF4PY_PETSC_CALL(VecSetType, (o, otype));
    } catch(...) {
      if (o != PETSC_NULL) { VecDestroy(o); o = PETSC_NULL; }
      throw;
    }
  } else {
    // check object handle
    PF4PY_PETSC_CALL(VecValid, (o, &valid));
    if (!valid)
      throw Error("invalid handle in " + name + " vector");
    // check object comm
    MPI_Comm ocomm;
    PF4PY_PETSC_CALL(PetscObjectGetComm,((PetscObject)o, &ocomm));
    if (comm != ocomm)
      throw Error("incompatible communicator in " + name + " vector");
    // check sizes
    PF4PY_PETSC_CALL(VecGetLocalSize, (o, &n));
    PF4PY_PETSC_CALL(VecGetSize,      (o, &N));
    if (n >= 0 && n != sizes.first)
      throw Error("incompatible local size in " + name + " vector");
    if (N >= 0 && N != sizes.second)
      throw Error("incompatible global size in " + name + " vector");
    if (n < 0 || N < 0) {
      n = sizes.first; N = sizes.second;
      PF4PY_PETSC_CALL(VecSetSizes, (o, n, N));
      if (!otype) otype = (comm.getSize() == 1) ? VECSEQ : VECMPI;
    }
    if (otype) PF4PY_PETSC_CALL(VecSetType, (o, otype));
  }
}

static void 
mk_mat(Mat& o, Comm comm, const std::pair<int,int>& sizes,
       const std::string& type, const std::string& name)
{
  PetscTruth valid;
  PetscInt m, n, M, N;
  const char *otype = (type.size() > 0) ? type.c_str() : NULL;
  if (o == PETSC_NULL) {
    try {
      // create object
      m = n = sizes.first; M = N = sizes.second;
      PF4PY_PETSC_CALL(MatCreate, (comm, &o));
      PF4PY_PETSC_CALL(MatSetSizes, (o, m, n, M, N));
      // set object type
      if (!otype) otype = (comm.getSize() == 1) ? MATSEQAIJ : MATMPIAIJ;
      PF4PY_PETSC_CALL(MatSetType, (o, otype));
    } catch(...) {
      if (o != PETSC_NULL) { MatDestroy(o); o = PETSC_NULL; }
      throw;
    }
  } else {
    // check object handle
    PF4PY_PETSC_CALL(MatValid, (o, &valid));
    if (!valid)
      throw Error("invalid handle in " + name + " matrix");
    // check object comm
    MPI_Comm ocomm;
    PF4PY_PETSC_CALL(PetscObjectGetComm,((PetscObject)o, &ocomm));
    if (comm != ocomm)
      throw Error("incompatible communicator in " + name + " matrix");
    // check sizes
    PF4PY_PETSC_CALL(MatGetLocalSize, (o, &m, &n));
    PF4PY_PETSC_CALL(MatGetSize,      (o, &M, &N));
    if (m >= 0 && m != sizes.first)
      throw Error("incompatible local row size in " + name + " matrix");
    if (n >= 0 && n != sizes.first)
      throw Error("incompatible local column size in " + name + " matrix");
    if (M >= 0 && M != sizes.second)
      throw Error("incompatible global row size in " + name + " matrix");
    if (N >= 0 && N != sizes.second)
      throw Error("incompatible global column size in " + name + " matrix");
    if (n < 0 || m < 0 || M < 0 || N < 0) {
      m = n = sizes.first; M = N = sizes.second;
      PF4PY_PETSC_CALL(MatSetSizes, (o, m, n, M, N));
      if (!otype) otype = (comm.getSize() == 1) ? MATSEQAIJ : MATMPIAIJ;
    }
    if (otype) PF4PY_PETSC_CALL(MatSetType, (o, otype));
  }
}

static PetscErrorCode 
mat_dummy_setvals(Mat mat,
		  PetscInt m,const PetscInt idxm[],
		  PetscInt n,const PetscInt idxn[],
		  const PetscScalar v[],InsertMode addv) { return 0; }

static PetscErrorCode 
mat_dummy_assem(Mat mat, MatAssemblyType type) { return 0; }

static void 
mk_mat_dummy(Mat& o, Comm comm, const std::pair<int,int>& sizes)
{
  PetscInt m, n, M, N;
  try {
    // create object
    m = n = sizes.first; M = N = sizes.second;
    PF4PY_PETSC_CALL(MatCreate, (comm, &o));
    PF4PY_PETSC_CALL(MatSetSizes, (o, m, n, M, N));
    // set object type
    PF4PY_PETSC_CALL(MatSetType, (o, MATSHELL));
    MatShellSetOperation(o,MATOP_SET_VALUES,
			 (void (*)(void)) mat_dummy_setvals);
    MatShellSetOperation(o,MATOP_ASSEMBLY_BEGIN,
			 (void (*)(void)) mat_dummy_assem);
    MatShellSetOperation(o,MATOP_ASSEMBLY_END,
			 (void (*)(void)) mat_dummy_assem);
  } catch(...) {
    if (o != PETSC_NULL) { MatDestroy(o); o = PETSC_NULL; }
    throw;
  }
}

void 
Domain::allocateSolution(Vec& u) const
{
  Comm     comm = PETSC_COMM_SELF;
  int      nnod = this->getNNod();
  int      ndof = this->getNDof();
  std::pair<int,int> sizes(nnod*ndof,nnod*ndof);
  mk_vec(u, comm, sizes, VECSEQ, "solution");
}

void 
Domain::allocateState(Vec& x, const std::string& type) const
{
  Dofset::Impl* dofmap = this->getDofset();
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  
  Comm               comm  = this->getComm();
  std::pair<int,int> sizes = this->dofset->getSizes();
  mk_vec(x, comm, sizes, type, "state");
}

void 
Domain::allocateResidual(Vec& r, const std::string& type) const
{
  Dofset::Impl* dofmap = this->getDofset();
  if (dofmap == NULL) throw Error("Domain: dofset not ready");

  Comm               comm  = this->getComm();
  std::pair<int,int> sizes = this->dofset->getSizes();
  mk_vec(r, comm, sizes, type, "residual");
}

void 
Domain::allocateJacobian(Mat& J, const std::string& type) const
{
  Dofset::Impl* dofmap = this->getDofset();
  if (dofmap == NULL) throw Error("Domain: dofset not ready");

  Comm               comm  = this->getComm();
  std::pair<int,int> sizes = this->dofset->getSizes();
  
  if (type == std::string("dummy")) {
    mk_mat_dummy(J, comm, sizes); return;
  }
  mk_mat(J, comm, sizes, type, "jacobian");
  if (!this->appctx) return;
  
  if (type == std::string(MATIS)) { 

    std::vector<int> dofs, xadj, adjncy;
    if (!this->appctx->sdgraph(this, dofs, xadj, adjncy)) return;
    MPI_Comm mcomm = MPI_COMM_NULL;
    ISLocalToGlobalMapping mapping = PETSC_NULL;
    Mat A = PETSC_NULL;
    PF4PY_PETSC_CALL(PetscObjectGetComm,((PetscObject)J, &mcomm));
    PF4PY_PETSC_CALL(ISLocalToGlobalMappingCreate,
		     (mcomm,dofs.size(),&dofs[0],&mapping));
    int nout, nin = adjncy.size();
    int *idx = &adjncy[0];
    PF4PY_PETSC_CALL(ISGlobalToLocalMappingApply,
		     (mapping,IS_GTOLM_DROP,nin,idx,&nout,idx));
    if (nin != nout) throw Error("Domain: internal error, "
				 "bad subdomain graph");
    PF4PY_PETSC_CALL(MatSetLocalToGlobalMapping,(J,mapping));
    PF4PY_PETSC_CALL(ISLocalToGlobalMappingDestroy,(mapping));
    PF4PY_PETSC_CALL(MatISGetLocalMat,(J,&A));
    PF4PY_PETSC_CALL(MatSeqAIJSetPreallocationCSR,
		     (A, &xadj[0], &adjncy[0], PETSC_NULL));
    return;

  } else {

    std::vector<int> xadj, adjncy;
    if (!this->appctx->profile(this, xadj, adjncy)) return;
    if (comm.getSize() == 1)
      PF4PY_PETSC_CALL(MatSeqAIJSetPreallocationCSR,
		       (J, &xadj[0], &adjncy[0], PETSC_NULL));
    else
      PF4PY_PETSC_CALL(MatMPIAIJSetPreallocationCSR,
		       (J, &xadj[0], &adjncy[0], PETSC_NULL));
  }
}

void 
Domain::assemble(double t, Vec x, Vec r, Mat J) const 
{
  assert(t == t);

  Mesh::Impl*   mesh   = this->getMesh();
  Dofset::Impl* dofmap = this->getDofset();

  if (mesh   == NULL) throw Error("Domain: mesh not ready");
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  if (!this->appctx)  throw Error("Domain: appctx not ready");
  
  this->appctx->assemble(this, t, x, r, J);
}

void
Domain::assemble(double t1, Vec x1, double t0, Vec x0, 
		 Vec r, Mat J, double alpha) const 
{
  assert(t1 == t1);
  assert(t0 == t0);

  Mesh::Impl*   mesh   = this->getMesh();
  Dofset::Impl* dofmap = this->getDofset();

  if (mesh   == NULL) throw Error("Domain: mesh not ready");
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  if (!this->appctx)  throw Error("Domain: appctx not ready");

  if (alpha != alpha) 
    throw Error("Domain: invalid value, 'alpha' is NaN");
  if (alpha < 0.0) 
    throw Error("Domain: invalid value, 'alpha' < 0.0");
  if (alpha > 1.0) 
    throw Error("Domain: invalid value, 'alpha' > 1.0");

  this->appctx->assemble(this, t1, x1, t0, x0, r, J, alpha);
}


void 
Domain::assemble(const std::string& jobname,
		 double t, Vec x, 
		 Vec r, Mat J) const 
{
  assert(t == t);

  Mesh::Impl*   mesh   = this->getMesh();
  Dofset::Impl* dofmap = this->getDofset();

  if (mesh   == NULL) throw Error("Domain: mesh not ready");
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  if (!this->appctx)  throw Error("Domain: appctx not ready");
  
  this->appctx->assemble(this, jobname, t, x, r, J);
}

void 
Domain::assemble(const std::string& jobname,
		 double t1, Vec x1, double t0, Vec x0,
		 Vec r, Mat J, double alpha) const
{
  assert(t1 == t1);
  assert(t0 == t0);

  Mesh::Impl*   mesh   = this->getMesh();
  Dofset::Impl* dofmap = this->getDofset();

  if (mesh   == NULL) throw Error("Domain: mesh not ready");
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  if (!this->appctx)  throw Error("Domain: appctx not ready");

  if (alpha != alpha) 
    throw Error("Domain: invalid value, 'alpha' is NaN");
  if (alpha < 0.0) 
    throw Error("Domain: invalid value, 'alpha' < 0.0");
  if (alpha > 1.0) 
    throw Error("Domain: invalid value, 'alpha' > 1.0");
  
  this->appctx->assemble(this, jobname, t1, x1, t0, x0, r, J, alpha);
}

#if 0
// This is a simple Richardon iteration with Jacobi preconditioning to
// solve the overdetermined problem Q*x=y.
// The problem is solved in a least-square sense: Q'*Q*x = Q*y.
// The Richardson/Jacobi iteration is then
// x^{n+1} = x^{n} + inv(D)*(y - Q'*Q*x^{n}), where D = diag(Q'*Q)
static std::string
dofmap_solve(Dofset::Impl* dofmap,
	     double *x, double *y,
	     int& iter, double& r0, double& r,
	     double rtol=1e-10, double atol=1e-10, int niter=100, 
	     double omega=1.0)
{
  assert(x != NULL);
  assert(y != NULL);
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  int nrow = nnod*ndof;
  int ncol = dofmap->neqtot;
  int m;
  vector<double> dv(ncol, 0.0); double *d = &dv[0];
  vector<double> zv(nrow, 0.0); double *z = &zv[0];
  vector<double> vv(nrow, 0.0); double *v = &vv[0];
  vector<double> wv(ncol, 0.0); double *w = &wv[0];
  // z <- Q'*y
  dofmap->qtxpy(z, y, 1.0);
  // d <- diag(Q'*Q)
  for (int j=0; j<nrow; j++) {
    int kdof = j % ndof + 1;
    int node = (j / ndof) + 1;
    const int *dofs;
    const double *coef;
    dofmap->get_row(node,kdof,m,&dofs,&coef);
    for (int k=0; k<m; k++) {
      int dof = dofs[k]-1;
      d[dof] += square(coef[k]);
    }
  }
  for (int j=0; j<ncol; j++) 
    if (d[j] == 0.0) d[j] = 1.0;

  for (iter=0; iter<niter; iter++) {
    // v <- Q*x
    for (int j=0; j<nrow; j++) v[j]=0.0; // memset(v, 0, sizeof(double)*nrow);
    dofmap->qxpy(v, x, 1.0);
    // w <- Q'*Q*x
    for (int j=0; j<ncol; j++) w[j]=0.0; // memset(w, 0, sizeof(double)*ncol);
    dofmap->qtxpy(w, v, 1.0);
    // x <- x + omega *D^(-1)(z - w) 
    // D = diag(Q), z = Q'*y; w = Q'*Q*x
    r = 0.0;
    for (int j=0; j<ncol; j++) {
      x[j] += omega*(z[j]-w[j])/d[j];
      r += square(z[j]-w[j]);
    }
    r = sqrt(r);
    //printf("i: %3d |r|: %f\n", iter, r);
    if (iter==0) r0 = r;
    if (r < atol || r < rtol*r0) break;
  }
  if (iter < niter) return "";
  std::stringstream message;
  message << "Domain: projection did not converge "
	  << "in "   << niter << " iterations. "
	  << "rtol: " << rtol << ", "
	  << "|r0|: " << r0   << ", "
	  << "atol: " << atol << ", "
	  << "|r|: "  << r    << "";
  return message.str();
}
#endif


// We have to solve the system Q*x = y. 
// In this version we solve it with the CG
// on the normal equations. Normally the matrix Q
// is a permutation and in other cases it is very well
// conditioned so that this is a good choice. So we
// solve Q'*Q*x = Q'*y, -> H*x = z
#undef __FUNC__
#define __FUNC__ "dofmap::mult"
static 
PetscErrorCode dofmap_multQtQ(Mat A,Vec x,Vec y) 
{
  PetscErrorCode ierr;
  Dofmap         *dofmap = NULL;
  PetscFunctionBegin;
  ierr = MatShellGetContext(A,(void**)&dofmap);CHKERRQ(ierr); 
  dofmap->mult(A,x,y);
  PetscFunctionReturn(0);
}
#undef __FUNC__
#define __FUNC__ "dofmap::solve"
static std::string
dofmap_solve(Dofset::Impl* dofmap,
	     double *state, double *solution,
	     int& iter, double& r0, double& r,
	     double rtol=1e-5, double atol=1e-10, int maxit=100)
{
  PetscInt its;
  Vec x, z;
  Mat H;
  KSP ksp;
  PC pc;
  PetscScalar *xp,*zp;
  PetscErrorCode ierr;

  assert(state    != NULL);
  assert(solution != NULL);
  //int nnod = dofmap->nnod;
  //int ndof = dofmap->ndof;
  //int nrow = nnod*ndof;
  int ncol = dofmap->neqtot;

  PetscFunctionBegin;

  // Define auxiliar Vec's
  ierr = VecCreateSeq(PETSC_COMM_SELF, ncol, &x);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ncol, &z);

  // Define auxiliar Mat (shell)
  ierr = MatCreateShell(PETSC_COMM_SELF,ncol,ncol,ncol,ncol,dofmap,&H);
  MatShellSetOperation(H,MATOP_MULT,(void(*)(void))dofmap_multQtQ);

  // Define auxiliar KSP and PC
  ierr = KSPCreate(PETSC_COMM_SELF,&ksp);
  ierr = KSPSetOptionsPrefix(ksp,"pf_dofmap_");
  ierr = KSPSetType(ksp,KSPCG);
  ierr = KSPGetPC(ksp,&pc);
  ierr = PCSetType(pc,PCNONE);
  ierr = KSPSetOperators(ksp,H,H,SAME_NONZERO_PATTERN);
  ierr = KSPSetTolerances(ksp,rtol,atol,PETSC_DEFAULT,maxit);
  ierr = KSPSetFromOptions(ksp);

  // Compute rhs, z = Q'*y
  ierr = VecSet(z,0);
  ierr = VecGetArray(z,&zp);
  dofmap->qtxpy(zp,solution,1.0);
  ierr = VecRestoreArray(z,&zp);
  ierr = VecNorm(z,NORM_2,&r0);

  // Solve problem
  ierr = KSPSolve(ksp,z,x);
  ierr = KSPGetResidualNorm(ksp,&r);
  ierr = KSPGetIterationNumber(ksp,&its);
  ierr = KSPGetTolerances(ksp,&rtol,&atol,0,0);

  // Copy solution
  ierr = VecGetArray(x,&xp);
  ierr = PetscMemcpy(state,xp,ncol*sizeof(PetscScalar));
  ierr = VecRestoreArray(x,&xp);
  
  // Destroy auxiliary objects
  ierr = VecDestroy(x);
  ierr = VecDestroy(z);
  ierr = MatDestroy(H);
  ierr = KSPDestroy(ksp);
  
  if (its < maxit) PetscFunctionReturn("");
  
  std::stringstream message;
  message << "Domain: projection did not converge "
	  << "in "    << maxit << " iterations. "
	  << "rtol: " << rtol << ", "
	  << "|r0|: " << r0   << ", "
	  << "atol: " << atol << ", "
	  << "|r|: "  << r    << "";

  PetscFunctionReturn(message.str());
}


int
Domain::buildState(double time, Vec solution, Vec state) const
{
  assert(time == time);

  Dofset::Impl* dofmap = *this->dofset;
  if (dofmap == NULL) throw Error("Domain: dofset not ready");

  int nnod = this->dofset->getNNod();
  int ndof = this->dofset->getNDof();
  std::pair<int,int> sizes = this->dofset->getSizes();
  int ldofs = sizes.first;
  int gdofs = sizes.second;

  Comm comm = this->getComm();
  int  rank = comm.getRank();

  int root = 0;

  /* check provided state vector */
  PetscTruth stt_valid; VecValid(state,&stt_valid);
  PF4PY_ASSERT(stt_valid == PETSC_TRUE, "provided state vector is not valid");
  PetscInt n, N; VecGetLocalSize(state,&n); VecGetSize(state,&N);
  PF4PY_ASSERT(ldofs == n,"provided state vector has wrong local size");
  PF4PY_ASSERT(gdofs == N,"provided state vector has wrong global size");

  /* check provided solution vector */
  if (rank == root) {
    PetscTruth sol_valid; VecValid(solution,&sol_valid);
    PF4PY_ASSERT(sol_valid == PETSC_TRUE, "provided solution vector is not valid");
    PetscInt sol_size; VecGetLocalSize(solution,&sol_size);
    PF4PY_ASSERT(sol_size == nnod*ndof,   "provided solution vector has wrong local size ");
  }
  
  /* create scatter */
  std::pair<VecScatter,Vec> scatter(PETSC_NULL,PETSC_NULL);
  if (!scatter.first) PF4PY_PETSC_CALL(VecScatterCreateToAll,
				       (state, &scatter.first, &scatter.second));
  
  std::string result;

  /* build state in root processor */
  int i=0; double r0=0.0; double r=0.0;
  PetscScalar* stt_array;
  VecGetArray(scatter.second,&stt_array);
  if (rank == root) {
    // get solution array
    PetscScalar* sol_array;
    VecGetArray(solution,&sol_array);
    // solve
    PetscScalar* sttbuff = PETSC_NULL;
    PetscMalloc(dofmap->neqtot*sizeof(PetscScalar), &sttbuff);
    result = dofmap_solve(dofmap, &sttbuff[0], sol_array, i, r0, r);
    PetscMemcpy(stt_array, &sttbuff[0], dofmap->neq*sizeof(PetscScalar));
    PetscFree(sttbuff);
    // restore solution array
    VecRestoreArray(solution,&sol_array);
    //
  } else {
    PetscMemzero(stt_array,gdofs*sizeof(PetscScalar));
  }
  VecRestoreArray(scatter.second,&stt_array);
  MPI_Bcast(&i, 1, MPI_INT, root, comm);

  /* scatter state values to all processors */
  VecZeroEntries(state);
  PF4PY_PETSC_CALL(VecScatterBegin, (scatter.first, scatter.second, state,
				     ADD_VALUES, SCATTER_REVERSE));
  PF4PY_PETSC_CALL(VecScatterEnd,   (scatter.first, scatter.second, state, 
				     ADD_VALUES, SCATTER_REVERSE));
  
  /* destroy scatter */
  PF4PY_PETSC_DESTROY(VecScatterDestroy, scatter.first);
  PF4PY_PETSC_DESTROY(VecDestroy,        scatter.second);

  if (result.size() > 0) throw Error(result);
  
  return i;
}

static void
dofmap_apply(Dofset::Impl* dofmap, 
	     double t, const double x[], double u[])
{
  assert(t == t);
  assert(x != NULL);
  assert(u != NULL);
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  //int neq  = dofmap->neq;
  Time time; time.set(t);
  for (int i=0; i<nnod; i++)
    for (int j=0; j<ndof; j++)
      dofmap->get_nodal_value(i+1, j+1, x, &time, u[i*ndof+j]);
}

void 
Domain::buildSolution(double time, Vec state, Vec solution) const
{
  
  Dofset::Impl* dofmap = this->getDofset();
  if (dofmap == NULL) throw Error("Domain: dofset not ready");

  int nnod = this->dofset->getNNod();
  int ndof = this->dofset->getNDof();
  std::pair<int,int> sizes = this->dofset->getSizes();
  int ldofs = sizes.first;
  int gdofs = sizes.second;

  /* check provided state vector */
  PetscTruth stt_valid; VecValid(state,&stt_valid);
  PF4PY_ASSERT(stt_valid == PETSC_TRUE, "provided state vector is not valid");
  PetscInt n, N; VecGetLocalSize(state,&n); VecGetSize(state,&N);
  PF4PY_ASSERT(ldofs == n, "provided state vector has wrong local size");
  PF4PY_ASSERT(gdofs == N, "provided state vector has wrong global size");

  /* check provided solution vector */
  PetscTruth sol_valid; PetscInt sol_size = 0; 
  if (solution != PETSC_NULL) {
    VecValid(solution,&sol_valid);
    PF4PY_ASSERT(sol_valid == PETSC_TRUE,
		 "provided solution vector is not valid");
    VecGetLocalSize(solution,&sol_size);
    PF4PY_ASSERT(sol_size == 0 || sol_size == nnod*ndof, 
		 "provided solution vector has wrong local size ");
  }
  
  /* create scatter */
  std::pair<VecScatter,Vec> scatter(PETSC_NULL,PETSC_NULL);
  if (!scatter.first) PF4PY_PETSC_CALL(VecScatterCreateToAll,
				       (state, &scatter.first, &scatter.second));
  /* scatter state values to all processors */
  PF4PY_PETSC_CALL(VecScatterBegin, (scatter.first, state, scatter.second,
				     INSERT_VALUES, SCATTER_FORWARD));
  PF4PY_PETSC_CALL(VecScatterEnd,   (scatter.first, state, scatter.second,
				     INSERT_VALUES,SCATTER_FORWARD));
  
  /* user does not want solution on this processor */
  if (solution == PETSC_NULL) goto exit;
  
  PetscScalar* sol_array;
  VecGetArray(solution, &sol_array);
  if (sol_size > 0 && sol_array != NULL) {
    PetscScalar* stt_array;
    VecGetArray(scatter.second, &stt_array);
    //
    if (sol_array != NULL)
      dofmap_apply(dofmap,time,stt_array,sol_array);
    /* array of state vector was not touched */
    VecRestoreArray(scatter.second, &stt_array);
  }
  VecRestoreArray(solution, &sol_array);

 exit:
  /* destroy scatter */
  PF4PY_PETSC_DESTROY(VecScatterDestroy, scatter.first);
  PF4PY_PETSC_DESTROY(VecDestroy,        scatter.second);
}


PF4PY_NAMESPACE_END
