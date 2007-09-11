// $Id$

#include <algorithm>
#include <numeric>

#include "Dofset.h"
#include "Domain.h"

#include <fem.h>
#include <dofmap.h>
#include <readmesh.h>


PF4PY_NAMESPACE_BEGIN
class Dofset::Proxy
{
private:
  Proxy();
  Proxy(const Proxy&);
  Proxy& operator=(const Proxy&);

protected:
  std::auto_ptr< Dofset::Impl > _ptr;
public:
  operator Dofset::Impl*() const { return this->_ptr.get(); }

protected:
  std::vector<int>   _neqproc;
  VecScatter         _ghost_scatter;
  VecScatter         _scatter_print;

public:
  ~Proxy()
  {
    Dofset::Impl* dofmap = *this;
    if (dofmap == NULL) return;
#if 0
    // private:
    PF4PY_DELETE_VCTR(dofmap->idmap2);
    PF4PY_DELETE_VCTR(dofmap->special_ptr);
    PF4PY_DELETE_VCTR(dofmap->sp_eq);
    PF4PY_DELETE_VCTR(dofmap->coefs);
#endif
    // public:
    PF4PY_DELETE_SCLR(dofmap->ident);
    PF4PY_DELETE_SCLR(dofmap->id);
    PF4PY_DELETE_FUNC(da_destroy, dofmap->fixa);
    PF4PY_DELETE_FUNC(da_destroy, dofmap->q);

    dofmap->startproc  = NULL;
    dofmap->neqproc    = NULL;
    dofmap->tpwgts     = NULL;
    dofmap->npart      = NULL;
    dofmap->ghost_dofs = NULL;

    PF4PY_PETSC_DESTROY(VecScatterDestroy, *(dofmap->ghost_scatter));
    PF4PY_PETSC_DESTROY(VecScatterDestroy, *(dofmap->scatter_print));
  }

  Proxy(Dofset* dofset, Dofset::Impl* impl)
    : _ptr(impl),
      _neqproc(),
      _ghost_scatter(0),
      _scatter_print(0)
  {
    Dofset::Impl* dofmap = *this;
    if (dofmap == NULL) return;
    
    Comm comm = dofset->getComm();
    int  nnod = dofset->getNNod();
    int  ndof = dofset->getNDof();
    int  size = comm.getSize();
    
    this->_neqproc.resize(size+1, 0);

    // preallocation
    std::vector<float>& tpwgts     = dofset->tpwgts;
    std::vector<int>&   startproc  = dofset->ownership;
    std::vector<int>&   neqproc    = this->_neqproc;
    std::vector<int>&   ghost_dofs = dofset->ghosts;
    //ghost_dofs.reserve(??);

    dofmap->comm          = MPI_COMM_NULL;
    dofmap->ident         = NULL;
    dofmap->ghost_dofs    = NULL;
    dofmap->id            = NULL;
    dofmap->startproc     = NULL;
    dofmap->neqproc       = NULL;
    dofmap->tpwgts        = NULL;
    dofmap->npart         = NULL;
    dofmap->ghost_scatter = NULL;
    dofmap->scatter_print = NULL;
    dofmap->ghost_dofs    = NULL;
    //
    dofmap->comm          = comm;
    dofmap->size          = size;
    dofmap->nnod          = nnod;
    dofmap->ndof          = ndof;
    dofmap->tpwgts        = &tpwgts[0];
    dofmap->startproc     = &startproc[0];
    dofmap->neqproc       = &neqproc[0];
    dofmap->neq           = 0;
    //dofmap->npart         = &npart[0];
    dofmap->ghost_dofs    = &ghost_dofs;
    dofmap->ghost_scatter = &this->_ghost_scatter;
    dofmap->scatter_print = &this->_scatter_print;
    dofmap->id            = new idmap(nnod*ndof, NULL_MAP);
  }

  void fill_idmap(Mesh::Impl* mesh)
  { 
    Dofset::Impl* dofmap = *this;
    assert(mesh != NULL);
    assert(mesh->elemsetlist != NULL);
    assert(dofmap->id != NULL);
    std::size_t nelemset = int(da_length(mesh->elemsetlist));
    for (std::size_t i=0; i<nelemset; i++) {
      Elemset::Impl** eptr = (Elemset::Impl**) da_ref(mesh->elemsetlist, i);
      assert (eptr != NULL);
      Elemset::Impl* elemset = *eptr;
      assert (elemset != NULL);
      int  ndof   = dofmap->ndof;
      int  nelem  = elemset->nelem;
      int  nel    = elemset->nel;
      int* icone  = elemset->icone;
      assert(icone != NULL);
      for (int j=0; j<nelem*nel; j++) {
	for (int k=1; k<=ndof; k++) {
	  int edof = dofmap->edof(icone[j], k);
	  dofmap->id->set_elem(edof, edof, 1.0);
	}
      }
    }
  }
  
  void add_f_entry(const Dofset::Fixation& f)
  {
    Dofset::Impl* dofmap = *this;
    // extract entries
    int node  = f.node+1;
    int field = f.field+1;
    double value = f.value;
    ::Amplitude* amp = NULL;
    if (f.amp != NULL) amp = *(f.amp);
    // add fixation
    row_t row; dofmap->get_row(node, field, row);
    if (row.size() != 1 || row.begin()->second != 1.0)
      throw Error("Dofset: bad fixation for node/field combination");
    std::vector<fixation_entry>& fixed      = dofmap->fixed;
    std::map<int,int>&           fixed_dofs = dofmap->fixed_dofs;
    int keq  = row.begin()->first;
    int edof = amp ? dofmap->edof(node, field) : -1;
    fixation_entry f_entry(value, amp, edof);
    if (fixed_dofs.find(keq) == fixed_dofs.end()) {
      fixed.push_back(f_entry);
      fixed_dofs[keq] = fixed.size()-1;
    } else {
      fixed[fixed_dofs[keq]] = f_entry;
    }
  }

  void add_c_entry(const Dofset::Constraint& constraint)
  {
    Dofset::Impl* dofmap = *this;
    // fill constraint
    ::Constraint C; C.reserve(constraint.size());
    for (std::size_t i=0; i<constraint.size(); i++) {
      const Dofset::constraint& c = constraint[i];
      C.add_entry(c.node+1, c.field+1, c.coeff);
    }
    // add constraint
    dofmap->set_constraint(C);
  }
  
}; // class Dofset::Proxy

Dofset::Impl* Dofset::getimpl() const
{ 
  Dofset::Proxy* proxy = this->proxy.get();
  if (proxy != NULL) return *proxy;
  return NULL;
}

PF4PY_NAMESPACE_END


PF4PY_NAMESPACE_BEGIN

Dofset::~Dofset()
{ }

Dofset::Dofset()
  : comm(MPI_COMM_SELF),
    nnod(0),
    ndof(0),
    tpwgts(1,1),
    ownership(2, 0),
    ghosts(),
    fixations(),
    amplitudes(),
    constraints(),
    proxy()
{ 
  throw Error("Dofset::Dofset()"); 
}


Dofset::Dofset(const Dofset& dofset)
  : comm(dofset.comm),
    nnod(dofset.nnod), 
    ndof(dofset.ndof),
    tpwgts(dofset.tpwgts),
    ownership(dofset.ownership),
    ghosts(dofset.ghosts),
    fixations(dofset.fixations),
    amplitudes(dofset.amplitudes),
    constraints(dofset.constraints),
    proxy()

{ 
  throw Error("Dofset::Dofset(const Dofset&)"); 
}

Dofset::Dofset(MPI_Comm comm, int nnod, int ndof)
  :  comm(comm),
     nnod(nnod),
     ndof(ndof),
     tpwgts(),
     ownership(),
     ghosts(),
     fixations(),
     amplitudes(),
     constraints(),
     proxy()
  
{ 
  int size = Comm(this->comm).getSize();
  this->tpwgts.resize(size, 1.0/float(size));
  this->ownership.resize(size+1, 0);
}

inline void
Dofset::add_fixa(int n, int f, double v, Amplitude* a)
{
  this->fixations.push_back(Dofset::Fixation(n, f, v, a));
}

void
Dofset::add_fixations(int n,
		      const int    node[],
		      const int    field[],
		      const double value[],
		      Amplitude*   amp)
{
  if (n == 0) return;
  assert(node  != NULL);
  assert(field != NULL);
  this->fixations.reserve(this->fixations.size() + n);
  if (value == NULL)
    for (int i=0; i<n; i++) 
      this->add_fixa(node[i], field[i], 1.0, amp);
  else
    for (int i=0; i<n; i++) 
      this->add_fixa(node[i], field[i], value[i], amp);
  this->amplitudes.add(amp);
}

void
Dofset::add_constraints(int n, 
			const double coeff[],
			const int    node[],
			const int    field[])
{
  if (n == 0) return;
  assert(coeff != NULL);
  assert(node  != NULL);
  assert(field != NULL);
  this->constraints.push_back(Dofset::Constraint());
  Dofset::Constraint& c = this->constraints.back();
  c.reserve(n);
  for (int i=0; i<n; i++) 
    c.push_back(Dofset::constraint(coeff[i], node[i], field[i]));
}

void
Dofset::del_fixations()
{
  this->fixations.clear();
  this->amplitudes.clear();
}

void
Dofset::del_constraints()
{
  this->constraints.clear();
}

void
Dofset::clear()
{
  int ownsz = this->ownership.size();
  this->ownership.resize(0);
  this->ownership.resize(ownsz, 0);
  this->ghosts.clear();
  this->fixations.clear();
  this->amplitudes.clear();
  this->constraints.clear();
  this->proxy.reset(NULL);
}

void 
setup_dofmap(MPI_Comm,
	     Mesh::Impl*,
	     Dofset::Impl*,
	     std::vector<int>&);

void
Dofset::setup(Domain* domain)
{
  Mesh::Impl*   mesh   = domain->getMesh();
  if (mesh == NULL) throw Error("Dofset: mesh is not ready");
  // create dofmap
  Dofset::Impl* dofmap = new Dofset::Impl;
  this->proxy.reset(new Dofset::Proxy(this, dofmap));
  // fill idmap
  this->proxy->fill_idmap(mesh);
  // add fixations
  Dofset::FixationList::const_iterator f  = this->fixations.begin();
  Dofset::FixationList::const_iterator fe = this->fixations.end();
  while (f != fe) this->proxy->add_f_entry(*f++);
  // add constraints
  Dofset::ConstraintList::const_iterator c  = this->constraints.begin();
  Dofset::ConstraintList::const_iterator ce = this->constraints.end();
  while (c != ce) this->proxy->add_c_entry(*c++);
  // setup dofmap
  MPI_Comm           comm  = domain->getComm();
  std::vector<int>&  npart = domain->mesh->nodepart;
  setup_dofmap(comm, mesh, dofmap, npart);
}

PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN

const std::vector<float>&
Dofset::getWeights() const
{
  return this->tpwgts;
}

void 
Dofset::setWeights(const std::vector<float>& weights)
{
  if (weights.size() != this->tpwgts.size())
    throw Error("Dofset: invalid size for weights");
  else {
    std::vector<float>::const_iterator i = weights.begin();
    std::vector<float>::const_iterator n = weights.end();
    while (i != n) 
      if (*i++ < 0.0)
	throw Error("Dofset: negative weight entry");
  }
  std::copy(weights.begin(), weights.end(), this->tpwgts.begin());
  std::vector<float>::iterator i = this->tpwgts.begin();
  std::vector<float>::iterator n = this->tpwgts.end();
  float accum = std::accumulate(i, n, 0.0);
  while (i != n) *i++ /= accum;
}

std::pair<int,int>       
Dofset::getSizes() const
{
  int size; MPI_Comm_size(this->comm, &size);
  int rank; MPI_Comm_rank(this->comm, &rank);
  int local  = this->ownership[rank+1] - this->ownership[rank];
  int global = this->ownership[size];
  return std::pair<int,int>(local, global);
}

std::pair<int,int>       
Dofset::getSizes(int rank) const
{
  int size; MPI_Comm_size(this->comm, &size);
  if (rank < 0) rank += size;
  if (rank < 0 || rank >= size)
    throw Error("Dofset: invalid rank, out of range");
  int local  = this->ownership[rank+1] - this->ownership[rank];
  int global = this->ownership[size];
  return std::pair<int,int>(local, global);
}

std::pair<int,int>
Dofset::getRange() const
{
  int size; MPI_Comm_size(this->comm, &size);
  int rank; MPI_Comm_rank(this->comm, &rank);
  int first = this->ownership[rank];
  int last  = this->ownership[rank+1];
  return std::pair<int,int>(first, last);
}

std::pair<int,int>
Dofset::getRange(int rank) const
{
  int size; MPI_Comm_size(this->comm, &size);
  if (rank < 0) rank += size;
  if (rank < 0 || rank >= size)
    throw Error("Dofset: invalid rank, out of range");
  int first = this->ownership[rank];
  int last  = this->ownership[rank+1];
  return std::pair<int,int>(first, last);
}

std::vector<int>
Dofset::getDist() const
{
  Dofset::Impl* dofmap = *this;
  if (dofmap == NULL) throw Error("Dofset: dofmap not ready");
  return this->ownership;
}

void
Dofset::getGhostDofs(std::vector<int>& gdofs) const
{
  Dofset::Impl* dofmap = *this;
  if (dofmap == NULL) throw Error("Dofset: dofmap not ready");
  const std::vector<int>& ghosts = this->ghosts;
  gdofs.reserve(gdofs.size() + ghosts.size());
  gdofs.insert(gdofs.end(), ghosts.begin(), ghosts.end());
}

void
Dofset::getLocalDofs(std::vector<int>& ldofs) const
{
  Dofset::Impl* dofmap = *this;
  if (dofmap == NULL) throw Error("Dofset: dofmap not ready");
  std::pair<int,int> range = this->getRange();
  const int start = range.first;
  const int end   = range.second; 
  const std::vector<int>& ghosts = this->ghosts;
  std::vector<int>::const_iterator middle =
    std::lower_bound(ghosts.begin(), ghosts.end(), start);
  ldofs.reserve(ldofs.size() + ghosts.size() + (end-start));
  ldofs.insert(ldofs.end(), ghosts.begin(), middle);
  for (int dof=start; dof<end; ) ldofs.push_back(dof++);
  ldofs.insert(ldofs.end(), middle, ghosts.end());
}

void 
Dofset::getFieldDofs(int field, std::vector<int>& fdofs) const
{
  Dofset::Impl* dofmap = *this;
  if (dofmap == NULL) throw Error("Dofset: dofmap not ready");
  if (field < 0 || field >= this->getNDof())
    throw Error("Dofset: invalid field, out of range");
  const int nnod = this->getNNod();
  const int ndof = this->getNDof();
  std::pair<int,int> range = this->getRange();
  const int first = range.first; 
  const int last  = range.second;
  fdofs.clear(); fdofs.reserve(1+(last-first)/ndof);
  for (int node=0; node<nnod; node++) {
    int n; const int *dofp; const double *coefp;
    dofmap->get_row(node+1,field+1,n,&dofp,&coefp);
    for (int i=0; i<n; i++) {
      const int dof = dofp[i]-1;
      if (dof < first || dof >= last) continue;
      fdofs.push_back(dof);
    }
  }
  std::vector<int>::iterator begin, end, uniq;
  begin = fdofs.begin();
  end   = fdofs.end();
  std::sort(begin, end);
  uniq = std::unique(begin, end);
  fdofs.erase(uniq, end);
}


PF4PY_NAMESPACE_END
