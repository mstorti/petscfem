// $Id: Domain.cpp,v 1.1.2.2 2006/05/25 00:26:57 dalcinl Exp $

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
  PYPF_DECREF(this->mesh);
  PYPF_DECREF(this->dofset);
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
    mesh(domain.mesh),
    dofset(domain.dofset),
    dofmap(domain.dofmap)
{
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofset);
  PYPF_INCREF(this->dofmap);
}
    

Domain::Domain(Mesh& mesh,
	       Dofset& dofset)
  : Object(dofset.getComm()),
    mesh(&mesh),
    dofset(&dofset),
    dofmap(new DofMap(mesh, dofset))
{
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofset);
  PYPF_INCREF(this->dofmap);
  Domain_setUp(this);
}

Domain::Domain(Nodeset& nodeset,
	       const std::vector<Elemset*>& elemsets,
	       Dofset& dofset)
  : Object(dofset.getComm()),
    mesh(new Mesh(nodeset, elemsets)),
    dofset(&dofset),
    dofmap(new DofMap(*mesh, dofset))
{
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofset);
  PYPF_INCREF(this->dofmap);
  Domain_setUp(this);
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

int
Domain::getSize() const
{
  return this->dofmap->getNNod() * this->dofmap->getNDof();
}

void
Domain::getSizes(int* nnod, int* ndof) const
{
  if (nnod) *nnod = this->dofmap->getNNod();
  if (ndof) *ndof = this->dofmap->getNDof();
}

int
Domain::getDofSize() const
{
  return this->dofmap->getSize();
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


#include <set>

void 
Domain::getOwnedDofs(int* start, int* end) const 
{
  this->dofmap->getRange(start, end);
}

void
Domain::getGhostDofs(std::vector<int>& gdofs) const
{
  const DofMap& dofmap = this->getDofMap();
  const std::vector<int>& ghost_dofs = *dofmap->ghost_dofs;
  gdofs.reserve(gdofs.size() + ghost_dofs.size());
  gdofs.insert(gdofs.end(), ghost_dofs.begin(), ghost_dofs.end());
}

void
Domain::getLocalDofs(std::vector<int>& ldofs) const
{
  const DofMap& dofmap = this->getDofMap();
  // add ghost dofs
  std::vector<int>& ghost_dofs = *dofmap->ghost_dofs;
  std::set<int> dofset(ghost_dofs.begin(), ghost_dofs.end());
  // add ghost in local range
  int start,end; this->getDofRange(&start, &end);
  for (int dof=start; dof<end; dof++) dofset.insert(dof);
  // fill output vector
  ldofs.reserve(ldofs.size() + dofset.size());
  ldofs.insert(ldofs.end(), dofset.begin(), dofset.end());

#if 0
  // old version !!

  Mesh&   mesh   = this->getMesh();
  DofMap& dofmap = this->getDofMap();

  // iterate over all elemsets and build
  // the set of nodes of local elements
  std::set<int> nodset;
  int proc; MPI_Comm_rank(this->getComm(), &proc); proc++;
  for (int elset=0; elset<mesh.getSize(); elset++) {
    const Elemset& elemset = mesh.getElemset(elset);
    int nelem, nel;
    const int *icone;
    int       *part;
    elemset.getData(&nelem, &nel, &icone);
    elemset.getPart(NULL, &part);
    for (int i=0; i<nelem; i++)
      if (part[i] == proc)
	for (int j=0; j<nel; j++)
	  nodset.insert(icone[i*nel+j]);
  }
  
  // iterate over set of nodes of local elements
  // and build the set of associated dofs
  std::set<int> dofset;
  int           nnod = nodset.size();
  int           ndof = dofmap.getNDof();
  int           neqs = dofmap.getSize();
  idmap*        id   = dofmap->id;
  IdMapRow      row;
  set<int>::iterator node_it = nodset.begin();
  for (int n=0; n<nnod; n++, node_it++) {
    int node = *node_it;
    for (int field=1; field<=ndof; field++) {
      int edof = (node-1) * ndof + field; // dofmap->edof()
      id->get_row(edof, row);
      int nrows = row.size();
      for (int i=0; i<nrows; i++) {
	IdMapEntry* entry = &row[i];
	int dof = entry->j - 1;
	if (dof<neqs) dofset.insert(dof);
      }
    }
  }
  nodset.clear();
  // fill output vector
  ldofs.reserve(ldofs.size() + dofset.size());
  ldofs.insert(ldofs.end(), dofset.begin(), dofset.end());
#endif
}

void 
Domain::getDofGraph(std::vector<int>& xadj,
		    std::vector<int>& adjncy) const
{
  Mesh&   mesh   = this->getMesh();
  DofMap& dofmap = this->getDofMap();
  
}

PYPF_NAMESPACE_END
