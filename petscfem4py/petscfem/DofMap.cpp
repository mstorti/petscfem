// $Id: DofMap.cpp,v 1.1.2.4 2006/03/20 16:06:00 rodrigop Exp $

#include "DofMap.h"


#include <fem.h>
#include <dofmap.h>


PYPF_NAMESPACE_BEGIN


OptionTable*
DofMap::get_opt_table() const
{
  throw Error("no options for DofMap");
  return NULL;
}

DofMap::~DofMap() 
{ 
  DofMap::Base* dofmap = *this;

  // private:
#if 0
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
  if (dofmap->ghost_scatter)
    PYPF_PETSC_DESTROY(VecScatterDestroy, *(dofmap->ghost_scatter));
  if (dofmap->scatter_print)
    PYPF_PETSC_DESTROY(VecScatterDestroy, *(dofmap->scatter_print));
  PYPF_DELETE_SCLR(dofmap->ghost_scatter);
  PYPF_DELETE_SCLR(dofmap->scatter_print);
  PYPF_DELETE_SCLR(dofmap);
}

DofMap::DofMap()
  : Ptr(new DofMap::Base)
{ 
  DofMap::Base* dofmap = *this;

  dofmap->ident      = NULL;
  dofmap->ghost_dofs = NULL;
  dofmap->id         = NULL;
  dofmap->startproc  = NULL;
  dofmap->neqproc    = NULL;
  dofmap->tpwgts     = NULL;
  dofmap->npart      = NULL;
  dofmap->ghost_scatter = NULL;
  dofmap->scatter_print = NULL;
}

DofMap::DofMap(DofMap::Base* _dofmap)
  : Ptr(_dofmap)
{ 
  
}

#if 1
DofMap::DofMap(int nnod, int ndof)
  : Ptr(new DofMap::Base)
{ 
  MPI_Comm comm = PETSC_COMM_WORLD;

  int comm_size, comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  idmap* id = new idmap(nnod*ndof, NULL_MAP);

  int*   startproc = new int[comm_size+1];
  int*   neqproc   = new int[comm_size+1];
  int *  npart     = new int[nnod];
  float* tpwgts    = new float[comm_size];
  memset(startproc, 0, (comm_size+1) * sizeof(int));
  memset(neqproc,   0, (comm_size+1) * sizeof(int));
  memset(npart,   0,   (nnod)        * sizeof(int));
  for (int i=0; i<comm_size; i++) 
    tpwgts[i] = 1.0/float(comm_size);

  DofMap::Base* dofmap = *this;

  dofmap->comm = comm;
  dofmap->nnod = nnod;
  dofmap->ndof = ndof;
  dofmap->id   = id;

  dofmap->size      = comm_size;
  dofmap->startproc = startproc;
  dofmap->neqproc   = neqproc;
  dofmap->tpwgts    = tpwgts;
  dofmap->npart     = npart;
  
  for (int i=1; i<=nnod; i++)
    for (int j=1; j<=ndof; j++) {
      int idx = (*this)->edof(i,j);
      dofmap->id->set_elem(idx,idx,1.);
    }
}
#endif

void
DofMap::addFixations(int n, int node[], int field[], double value[])
{
  DofMap::Base* dofmap = *this;

  row_t row;
  for (int i=0; i<n; i++) {
    dofmap->get_row(node[i]+1, field[i]+1, row);
    if (row.size()!=1) throw Error("bad fixation for node/field combination");
    int keq = row.begin()->first;
#if 0
    dofmap->fixed.push_back(fixation_entry(value[i]));
    dofmap->fixed_dofs[keq] = dofmap->fixed.size()-1;
#else
    typedef fixation_entry fixentry;
    std::vector<fixentry>& fixed      = dofmap->fixed;
    std::map<int,int>&     fixed_dofs = dofmap->fixed_dofs;
    if (fixed_dofs.find(keq) == fixed_dofs.end()) {
      fixed.push_back(fixentry(value[i]));
      fixed_dofs[keq] = fixed.size()-1;
    } else {
      fixed[fixed_dofs[keq]] = fixentry(value[i]);
    }
#endif
  }
}

void
DofMap::addConstraints(int n, int node[], int field[], double coef[])
{
  /* test */
  if ((*this)->id == NULL) throw Error("null id map");

  DofMap::Base* dofmap = *this;

  ::Constraint constraint;
  for (int i=0; i<n; i++) {
    constraint.add_entry(node[i]+1, field[i]+1, coef[i]);
  }
  dofmap->set_constraint(constraint);
}

void
DofMap::getSizes(int* local, int* global) const
{
  DofMap::Base* dofmap = *this;

  int size, rank;
  MPI_Comm_size(dofmap->comm, &size);
  MPI_Comm_rank(dofmap->comm, &rank);
  *local  = dofmap->neqproc[rank];
  *global = dofmap->startproc[size];
}

void
DofMap::getRange(int* start, int* end) const
{
  DofMap::Base* dofmap = *this;

  int rank, dof1, dof2;
  MPI_Comm_rank(dofmap->comm, &rank);
  dofmap->dof_range(rank, dof1, dof2);
  *start = dof1;
  *end   = dof2+1;
}

void
DofMap::getRanges(int* size, int* ranges[]) const
{
  DofMap::Base* dofmap = *this;

  MPI_Comm_size(dofmap->comm, size);
  (*size)+=1;
  *ranges = dofmap->startproc;
}

void 
DofMap::getNnod(int* nnod) 
{ 
  *nnod = (*this)->nnod; 
};

void 
DofMap::getFixSize(int* neq_fix) {
  *neq_fix = (*this)->neqf;
};


void
DofMap::setUp()
{
  DofMap::Base* dofmap = *this;
  dofmap->freeze();
}

void
DofMap::view() const
{
  DofMap::Base* dofmap = *this;

  if (dofmap->id != NULL) dofmap->id->print();
}


PYPF_NAMESPACE_END
