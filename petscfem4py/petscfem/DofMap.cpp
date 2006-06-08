// $Id: DofMap.cpp,v 1.1.2.10 2006/06/08 15:44:52 dalcinl Exp $

#include "DofMap.h"

#include <fem.h>
#include <dofmap.h>


PYPF_NAMESPACE_BEGIN


DofMap::~DofMap() 
{ 
  DofMap::Base* dofmap = *this;
  if (dofmap == NULL) return;
#if 0
  // private:
  PYPF_DELETE_VCTR(dofmap->idmap2);
  PYPF_DELETE_VCTR(dofmap->special_ptr);
  PYPF_DELETE_VCTR(dofmap->sp_eq);
  PYPF_DELETE_VCTR(dofmap->coefs);
#endif
  // public:
  PYPF_DELETE_SCLR(dofmap->ident);
  PYPF_DELETE_SCLR(dofmap->ghost_dofs);
  PYPF_DELETE_SCLR(dofmap->id);
  PYPF_DELETE(da_destroy,dofmap->fixa);
  PYPF_DELETE(da_destroy,dofmap->q);
  PYPF_DELETE_VCTR(dofmap->startproc);
  PYPF_DELETE_VCTR(dofmap->neqproc);
  PYPF_DELETE_VCTR(dofmap->tpwgts);
  PYPF_DELETE_VCTR(dofmap->npart);
  if (dofmap->ghost_scatter && *(dofmap->ghost_scatter))
    PYPF_PETSC_DESTROY(VecScatterDestroy, *(dofmap->ghost_scatter));
  PYPF_DELETE_SCLR(dofmap->ghost_scatter);
  if (dofmap->scatter_print && *(dofmap->scatter_print))
    PYPF_PETSC_DESTROY(VecScatterDestroy, *(dofmap->scatter_print));
  PYPF_DELETE_SCLR(dofmap->scatter_print);
  delete dofmap;
}


DofMap::DofMap(Mesh& mesh, Dofset& dofset)
  : Handle(new DofMap::Base),
    Object(dofset.getComm()),
    ampset(dofset.amplitudes)
{ 
  this->build(mesh, dofset); 
}

DofMap::DofMap(Mesh& mesh, Dofset& dofset, MPI_Comm comm)
  : Handle(new DofMap::Base),
    Object(comm),
    ampset(dofset.amplitudes)
{ 
  this->build(mesh, dofset); 
}

void 
DofMap::build(const Mesh& mesh, const Dofset& dofset)
{
  DofMap::Base* dofmap = *this;

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
  
  dofmap->ghost_dofs    = new vector<int>;
  dofmap->ghost_scatter = new VecScatter;
  dofmap->scatter_print = new VecScatter;
  *dofmap->ghost_scatter = PETSC_NULL;
  *dofmap->scatter_print = PETSC_NULL;

  // allocation
  int nnod, ndof, size;
  dofset.getSizes(&nnod,&ndof);
  MPI_Comm comm = this->getComm();
  MPI_Comm_size(comm, &size);
  dofmap->comm      = comm;
  dofmap->size      = size;
  dofmap->nnod      = nnod;
  dofmap->ndof      = ndof;
  dofmap->startproc = new int[size+1];
  dofmap->neqproc   = new int[size+1];
  dofmap->tpwgts    = new float[size];
  dofmap->npart     = new int[nnod];
  memset(dofmap->startproc, 0, sizeof(int)*(size+1));
  memset(dofmap->neqproc,   0, sizeof(int)*(size+1));
  memset(dofmap->npart,     0, sizeof(int)*nnod);
  for (int i=0; i<size; dofmap->tpwgts[i++] = 1.0/float(size));

  // create and fill idmap
  dofmap->id = new idmap(nnod*ndof, NULL_MAP);
  for (int i=0; i<mesh.getSize(); i++) {
    Elemset& elemset = mesh.getElemset(i);
    int  nelem, nel; const int* icone;
    elemset.getData(&nelem, &nel, &icone);
    for (int j=0; j<nelem*nel; j++) {
      for (int k=0; k<ndof; k++) {
	int edof = dofmap->edof(icone[j], k+1);
	dofmap->id->set_elem(edof, edof, 1.0);
      }
    }
  }
  // add fixations  
  const Dofset::FixationList& fixations = dofset.fixations;
  Dofset::FixationList::const_iterator f = fixations.begin();
  while (f != fixations.end()) this->add_fixation(*f++);
  // add constraints
  const Dofset::ConstraintList& constraints = dofset.constraints;
  Dofset::ConstraintList::const_iterator c = constraints.begin();
  while (c != constraints.end()) this->add_constraint(*c++);
}

void 
DofMap::add_fixation(const Dofset::Fixation& f)
{
  DofMap& dofmap = *this;

  int node  = f.node+1;
  int field = f.field+1;
  double value = f.value;
  Amplitude::Base* amp = f.amp;

  row_t row; dofmap->get_row(node, field, row);
  if (row.size()!=1)
    throw Error("bad fixation for node/field combination");
  if (row.begin()->second != 1.0)
    throw Error("bad fixation for node/field combination");
  std::vector<fixation_entry>& fixed      = dofmap->fixed;
  std::map<int,int>&           fixed_dofs = dofmap->fixed_dofs;
  int keq  = row.begin()->first;
  int edof = amp ? dofmap->edof(node, field) : -1;
  fixation_entry fix_entry(value, amp, edof);
  if (fixed_dofs.find(keq) == fixed_dofs.end()) {
    fixed.push_back(fix_entry);
    fixed_dofs[keq] = fixed.size()-1;
  } else {
    fixed[fixed_dofs[keq]] = fix_entry;
  }
}

void 
DofMap::add_constraint(const Dofset::Constraint& constraint)
{
  DofMap& dofmap = *this;
  Constraint C;
  for (int i=0; i<constraint.size(); i++) {
    const Dofset::constraint& c = constraint[i];
    C.add_entry(c.node+1, c.field+1, c.coeff);
  }
  dofmap->set_constraint(C);
}


int 
DofMap::getSize() const
{
  const DofMap& dofmap = *this;
  return dofmap->startproc[dofmap->size];
}

int 
DofMap::getLocalSize() const
{
  const DofMap& dofmap = *this;
  int rank; MPI_Comm_rank(dofmap->comm, &rank);
  const int* range = dofmap->startproc;
  return range[rank+1] - range[rank];
}

void
DofMap::getSizes(int* local, int* global) const
{
  if (local)  *local  = this->getLocalSize();
  if (global) *global = this->getSize();
}

void
DofMap::getRange(int* first, int* last) const
{
  const DofMap& dofmap =*this;
  int rank; MPI_Comm_rank(dofmap->comm, &rank);
  const int* range = dofmap->startproc;
  if (first) *first = range[rank];
  if (last)  *last  = range[rank+1];
}

void
DofMap::getDist(int* size, int* ranges[]) const
{
  const DofMap& dofmap = *this;
  if (size)   *size   = dofmap->size + 1;
  if (ranges) *ranges = dofmap->startproc;
}

int
DofMap::getNNod() const
{ 
  const DofMap& dofmap = *this;
  return dofmap->nnod;
};

int
DofMap::getNDof() const
{ 
  const DofMap& dofmap = *this;
  return dofmap->ndof;
};


void
DofMap::apply(int nstt, const double state[],
	      int nsol, double solution[],
	      double _time) const
{
  DofMap& dofmap = const_cast<DofMap&>(*this);

  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  int neq  = dofmap->neq;

  PYPF_ASSERT(nstt==neq,       "invalid size for state array");
  PYPF_ASSERT(nsol==nnod*ndof, "invalid size for solution array");

  Time time; time.set(_time);
  for (int i=0; i<nnod; i++)
    for (int j=0; j<ndof; j++)
      dofmap->get_nodal_value(i+1, j+1, state, &time, solution[i*ndof+j]);
}

void
DofMap::solve(int nsol, const double solution[],
	      int nstt, double state[]) const
{
  DofMap& dofmap = const_cast<DofMap&>(*this);

  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  int neq  = dofmap->neq;
  int neqt = dofmap->neqtot;

  PYPF_ASSERT(nstt==neq,       "invalid size for state array");
  PYPF_ASSERT(nsol==nnod*ndof, "invalid size for solution array");

  std::vector<double> solvec; solvec.resize(nnod*ndof);
  std::vector<double> sttvec; sttvec.resize(neqt);

  memcpy(&solvec[0], solution, nnod*ndof*sizeof(double));
  dofmap->solve(&sttvec[0], &solvec[0]);
  memcpy(state, &sttvec[0], neq*sizeof(double));
}


void
DofMap::view() const
{
  const DofMap& dofmap = *this;
  if (dofmap->id != NULL) dofmap->id->print();
}

PYPF_NAMESPACE_END
