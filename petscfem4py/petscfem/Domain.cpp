
// $Id: Domain.cpp,v 1.1.2.8 2006/07/26 23:31:34 dalcinl Exp $

#include <algorithm>

#include "Domain.h"

#include <fem.h>
#include <dofmap.h>

int domain_setup(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap);
extern int            SIZE, MY_RANK;
extern Mesh*          GLOBAL_MESH;
extern TextHashTable* GLOBAL_OPTIONS;


PYPF_NAMESPACE_BEGIN

Domain::~Domain() 
{ 
  PYPF_DECREF(this->dofset);
  PYPF_DECREF(this->mesh);
  PYPF_DECREF(this->dofmap);
}

static void
Domain_setUp(const Domain* domain)
{
  PETSCFEM_COMM_WORLD = domain->getComm();
  GLOBAL_MESH         = domain->getMesh();
  MPI_Comm_size(PETSCFEM_COMM_WORLD, &SIZE);
  MPI_Comm_rank(PETSCFEM_COMM_WORLD, &MY_RANK);
  domain_setup(PETSCFEM_COMM_WORLD, GLOBAL_MESH, domain->getDofMap());
  PETSCFEM_COMM_WORLD = MPI_COMM_NULL;
  GLOBAL_MESH         = NULL;
}

Domain::Domain(const Domain& domain)
  : Object(domain),
    dofset(domain.dofset),
    mesh(domain.mesh),
    dofmap(domain.dofmap)
{
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofset);
  PYPF_INCREF(this->dofmap);
}
    
Domain::Domain(Nodeset& nodeset,
	       const std::vector<Elemset*>& elemsets,
	       Dofset& dofset)
  : Object(dofset.getComm()),
    dofset(&dofset),
    mesh(new Mesh(nodeset, elemsets)),
    dofmap(new DofMap(*mesh, dofset))
{
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofset);
  PYPF_INCREF(this->dofmap);

  Domain_setUp(this);
}

Domain::Domain(Nodeset& nodeset,
	       const std::vector<Elemset*>& elemsets,
	       Dofset& dofset, MPI_Comm comm)
  : Object(comm),
    dofset(&dofset),
    mesh(new Mesh(nodeset, elemsets, comm)),
    dofmap(new DofMap(*mesh, dofset, comm))
{
  this->dofset->setComm(comm);
  
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofset);
  PYPF_INCREF(this->dofmap);

  Domain_setUp(this);
}


Nodeset&
Domain::getNodeset() const
{
  return this->mesh->getNodeset();
}

Dofset&
Domain::getDofset() const
{
  return *this->dofset;
}

Mesh&
Domain::getMesh() const
{
  return *this->mesh;
}

DofMap&
Domain::getDofMap() const
{
  return *this->dofmap;
}

int
Domain::getDim() const
{
  return this->mesh->getNodeset().getDim();
}

void
Domain::getSizes(int* nnod, int* ndof) const
{
  if (nnod) *nnod = this->dofmap->getNNod();
  if (ndof) *ndof = this->dofmap->getNDof();
}

void
Domain::getDofSizes(int* local, int* global) const
{
  this->dofmap->getSizes(local, global);
}

void 
Domain::getDofRange(int* first, int* last) const 
{
  this->dofmap->getRange(first, last);
}
void 
Domain::getDofDist(int* rsize, int* ranges[]) const
{
  this->dofmap->getDist(rsize, ranges);
}

void 
Domain::getOwnedDofs(int* start, int* end) const 
{
  this->dofmap->getRange(start, end);
}

void
Domain::getGhostDofs(std::vector<int>& gdofs) const
{
  const DofMap& dofmap = this->getDofMap();
  const std::vector<int>& ghosts = *dofmap->ghost_dofs;
  gdofs.reserve(gdofs.size() + ghosts.size());
  gdofs.insert(gdofs.end(), ghosts.begin(), ghosts.end());
}

void
Domain::getLocalDofs(std::vector<int>& ldofs) const
{
  const DofMap& dofmap = this->getDofMap();
  int start,end; this->getOwnedDofs(&start, &end);
  const std::vector<int>& ghosts = *dofmap->ghost_dofs;
  std::vector<int>::const_iterator middle =
    std::lower_bound(ghosts.begin(), ghosts.end(), start);
  ldofs.reserve(ldofs.size() + ghosts.size() + (end-start));
  ldofs.insert(ldofs.end(), ghosts.begin(), middle);
  for (int dof=start; dof<end; ) ldofs.push_back(dof++);
  ldofs.insert(ldofs.end(), middle, ghosts.end());
}

void 
Domain::getSplitDofs(std::vector<std::vector<int> >& split) const
{
  int nnod,  ndof; this->getSizes(&nnod, &ndof);
  int first, last; this->getDofRange(&first, &last);
  const DofMap& dofmap = this->getDofMap();
  split.clear(); split.resize(ndof);
  for (int field=0; field<ndof; field++) {
    std::vector<int>& dofs = split[field];
    for (int node=0; node<nnod; node++) {
      int n; const int *dofp; const double *coefp;
      dofmap->get_row(node+1,field+1,n,&dofp,&coefp);
      for (int i=0; i<n; i++) {
	int dof = dofp[i]-1;
	if (dof < first || dof >= last) continue;
	dofs.push_back(dof);
      }
    }
    std::vector<int>::iterator first, last, uniq;
    first = dofs.begin();
    last  = dofs.end();
    std::sort(first, last);
    uniq = std::unique(first, last);
    dofs.erase(uniq, last);
  }
}

PYPF_NAMESPACE_END
