// $Id: DofMap.cpp,v 1.1.2.5 2006/03/28 22:13:25 rodrigop Exp $

#include "Mesh.h"
#include "DofMap.h"


#include <fem.h>
#include <dofmap.h>


PYPF_NAMESPACE_BEGIN


DofMap::~DofMap() 
{ 
  DofMap::Base* dofmap = *this;
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

// DofMap::DofMap()
//   : Handle(new DofMap::Base), Object(),
//     nnod(0), ndof(0)
// { 
//   DofMap::Base* dofmap = *this;

//   dofmap->ident      = NULL;
//   dofmap->ghost_dofs = NULL;
//   dofmap->id         = NULL;
//   dofmap->startproc  = NULL;
//   dofmap->neqproc    = NULL;
//   dofmap->tpwgts     = NULL;
//   dofmap->npart      = NULL;
//   dofmap->ghost_scatter = NULL;
//   dofmap->scatter_print = NULL;
// }

DofMap::DofMap(DofMap::Base* dm)
  : Handle(dm), Object(),
    nnod(dm->nnod), ndof(dm->ndof),
    frozen(true)
{ 
//   if (dofmap->ghost_dofs == NULL)
//     dofmap->ghost_dofs    = new vector<int>;
//   if (dofmap->ghost_scatter == NULL)
//     dofmap->ghost_scatter = new VecScatter;
//   if (dofmap->scatter_print == NULL)
//     dofmap->scatter_print = new VecScatter;
//   if (dofmap->id == NULL)
//     dofmap->id = new idmap(this->nnod*this->ndof, NULL_MAP);
}


DofMap::DofMap(Mesh* mesh, int _ndof)
  : Handle(new DofMap::Base), Object(),
    nnod(mesh->nodedata->nnod), ndof(_ndof),
    frozen(false)
{
  DofMap::Base* dofmap   = *this;

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

  int nnod = dofmap->nnod = this->nnod;
  int ndof = dofmap->ndof = this->ndof;
  dofmap->id = new idmap(nnod*ndof, NULL_MAP);

  for (int i=0; i<mesh->elemsetlist.size(); i++) {
    Elemset& elemset = *mesh->elemsetlist[i];
    int  nelem, nel, *icone;
    elemset.getConnectivity(&nelem, &nel, &icone);
    for (int j=0; j<nelem*nel; j++) {
      for (int k=1; k<=ndof; k++) {
	int edof = dofmap->edof(icone[j], k);
	dofmap->id->set_elem(edof, edof, 1.);
      }
    }
  }

  int size, rank;
  MPI_Comm_size(this->comm, &size);
  MPI_Comm_rank(this->comm, &rank);
  int*   startproc = new int[size+1];
  int*   neqproc   = new int[size+1];
  float* tpwgts    = new float[size];
  int *  npart     = new int[nnod];
  memset(startproc, 0, sizeof(int)*(size+1));
  memset(neqproc,   0, sizeof(int)*(size+1));
  memset(npart,     0, sizeof(int)*nnod);
  for (int i=0; i<size;  tpwgts[i++] = 1.0/float(size));
  dofmap->size      = size;
  dofmap->startproc = startproc;
  dofmap->neqproc   = neqproc;
  dofmap->tpwgts    = tpwgts;
  dofmap->npart     = npart;
  
}

static void 
dofmap_set_fixation(Dofmap* dofmap, int node, int field, double value) {
  row_t row;
  dofmap->get_row(node, field, row);
  if (row.size()!=1) throw Error("bad fixation for node/field combination");
  std::vector<fixation_entry>& fixed      = dofmap->fixed;
  std::map<int,int>&           fixed_dofs = dofmap->fixed_dofs;
  int keq = row.begin()->first;
  if (fixed_dofs.find(keq) == fixed_dofs.end()) {
    fixed.push_back(fixation_entry(value));
    fixed_dofs[keq] = fixed.size()-1;
  } else {
    fixed[fixed_dofs[keq]] = fixation_entry(value);
  }
}

void
DofMap::addFixations(int n, int node[], int field[], double value[])
{
  if (this->frozen) throw Error("DofMap object is frozen");
  DofMap& dofmap = *this;
  for (int i=0; i<n; i++)
    dofmap_set_fixation(dofmap, node[i]+1, field[i]+1, value[i]);
}

void
DofMap::addConstraints(int n, int node[], int field[], double coef[])
{
  if (this->frozen) throw Error("DofMap object is frozen");
  DofMap& dofmap = *this;
  ::Constraint constraint;
  for (int i=0; i<n; i++) 
    constraint.add_entry(node[i]+1, field[i]+1, coef[i]);
  dofmap->set_constraint(constraint);
}

void
DofMap::getSizes(int* local, int* global) const
{
  int size, rank;
  MPI_Comm_size(this->comm, &size);
  MPI_Comm_rank(this->comm, &rank);
  const DofMap& dofmap = *this;
  *local  = dofmap->neqproc[rank];
  *global = dofmap->startproc[size];
}

void
DofMap::getRange(int* start, int* end) const
{
  DofMap& dofmap = const_cast<DofMap&>(*this);
  int rank, dof1, dof2;
  MPI_Comm_rank(this->comm, &rank);
  dofmap->dof_range(rank, dof1, dof2);
  *start = dof1;
  *end   = dof2+1;
}

void
DofMap::getRanges(int* size, int* ranges[]) const
{
  const DofMap& dofmap = *this;
  *size = dofmap->size + 1;
  *ranges = dofmap->startproc;
}

int
DofMap::getNnod() const
{ 
  const DofMap& dofmap = *this;
  return dofmap->nnod;
};

int
DofMap::getNdof() const
{ 
  const DofMap& dofmap = *this;
  return dofmap->ndof;
};

int
DofMap::getNfix() const
{
  const DofMap& dofmap = *this;
  return dofmap->neqf;
};


void
DofMap::setUp()
{ }

void
DofMap::view() const
{
  const DofMap& dofmap = *this;
  if (dofmap->id != NULL) dofmap->id->print();
}


PYPF_NAMESPACE_END
