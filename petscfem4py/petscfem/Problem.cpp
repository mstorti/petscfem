// $Id: Problem.cpp,v 1.1.2.7 2006/05/21 00:26:05 dalcinl Exp $

#include "Problem.h"

#include <fem.h>
#include <dofmap.h>


int my_mesh_part(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap);
extern TextHashTable* GLOBAL_OPTIONS;
extern Mesh*          GLOBAL_MESH;
extern int            SIZE, MY_RANK;


PYPF_NAMESPACE_BEGIN


Problem::~Problem() 
{ 
  PYPF_DECREF(this->mesh);
  PYPF_DECREF(this->dofmap);
}


static void
Problem_setUp(MPI_Comm comm, Mesh& mesh, DofMap& dofmap)
{
  PETSCFEM_COMM_WORLD = comm;
  GLOBAL_MESH         = mesh;
  MPI_Comm_size(PETSCFEM_COMM_WORLD, &SIZE);
  MPI_Comm_rank(PETSCFEM_COMM_WORLD, &MY_RANK);
  my_mesh_part(comm, mesh, dofmap);
  PETSCFEM_COMM_WORLD = MPI_COMM_NULL;
  GLOBAL_MESH         = NULL;
}

Problem::Problem(const Problem& P)
  : Object(P),
    mesh(P.mesh), 
    dofmap(P.dofmap)
{
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofmap);
}

Problem::Problem(Mesh& mesh, DofMap& dofmap)
  : Object(mesh.getComm()),
    mesh(&mesh),
    dofmap(&dofmap)
{
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofmap);
  Problem_setUp(this->getComm(), this->getMesh(), this->getDofMap());
}

Problem::Problem(Nodeset& nodeset,
		 const std::vector<Elemset*>& elemsets,
		 Dofset& dofset)
  : Object(),
    mesh(new Mesh(nodeset, elemsets)),
    dofmap(new DofMap(*mesh, dofset))
{
  this->setComm(this->mesh->getComm());
  PYPF_INCREF(this->mesh);
  PYPF_INCREF(this->dofmap);
  Problem_setUp(this->getComm(), this->getMesh(), this->getDofMap());
}

Mesh&
Problem::getMesh() const
{
  return *this->mesh;
}

DofMap&
Problem::getDofMap() const
{
  return *this->dofmap;
}

int
Problem::getDim() const
{
  return this->mesh->getNodeset().getDim();
}

int
Problem::getSize() const
{
  return this->dofmap->getNNod() * this->dofmap->getNDof();
}

void
Problem::getSizes(int* nnod, int* ndof) const
{
  if (nnod) *nnod = this->dofmap->getNNod();
  if (ndof) *ndof = this->dofmap->getNDof();
}

int
Problem::getDofSize() const
{
  return this->dofmap->getSize();
}

void 
Problem::getDofSizes(int* local, int* global) const
{
  this->dofmap->getSizes(local, global);
}

void 
Problem::getDofRange(int* first, int* last) const 
{
  this->dofmap->getRange(first, last);
}

static void
solution2state(DofMap& dofmap, const double solution[], double state[])
{
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  int neq  = dofmap->neq;
  int neqt = dofmap->neqtot;
  std::vector<double> solvec; solvec.resize(nnod*ndof);
  std::vector<double> sttvec; sttvec.resize(neqt);

  PetscMemcpy(&solvec[0], solution, nnod*ndof*sizeof(double));
  dofmap->solve(&sttvec[0], &solvec[0]);
  PetscMemcpy(state, &sttvec[0], neq*sizeof(double));
}

static void
state2solution(DofMap& dofmap,
	       const double state[], double solution[], double t)
{
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  Time time; time.set(t);
  for (int i=0; i<nnod; i++)
    for (int j=0; j<ndof; j++)
      dofmap->get_nodal_value(i+1, j+1, state, &time, solution[i*ndof+j]);
}

static void 
check_vec_sizes(DofMap& dofmap, Vec sol_vec, Vec stt_vec) 
{
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  int neq  = dofmap->neq;
  int sol_size; VecGetLocalSize(sol_vec, &sol_size);
  int stt_size; VecGetLocalSize(stt_vec, &stt_size);
  if (sol_size != nnod*ndof)
    throw Error("invalid local size for solution vector");
  if (stt_size != neq)
    throw Error("invalid local size for state vector");
}

void
Problem::buildSolution(Vec state, Vec solution) const
{
  DofMap& dofmap = *this->dofmap;
  check_vec_sizes(dofmap, solution, state);

  double* stt_buff;
  double* sol_buff;
  VecGetArray(state,    &stt_buff);
  VecGetArray(solution, &sol_buff);

  state2solution(dofmap, stt_buff, sol_buff, 0.0);

  VecRestoreArray(state,    &stt_buff);
  VecRestoreArray(solution, &sol_buff);
  PetscObjectStateDecrease(state);
}

void
Problem::buildState(Vec solution, Vec state) const
{

  DofMap& dofmap = *this->dofmap;
  check_vec_sizes(dofmap, solution, state);
  
  double* sol_buff;
  double* stt_buff;
  VecGetArray(solution, &sol_buff);
  VecGetArray(state, &stt_buff);
  
  solution2state(dofmap, sol_buff, stt_buff);

  VecRestoreArray(solution, &sol_buff);
  VecRestoreArray(state,    &stt_buff);
  PetscObjectStateDecrease(solution);
}

#include <set>

void
Problem::getLocalDofs(int* _ndofs, int* _dofs[]) const
{
  Mesh&   mesh = *this->mesh;
  DofMap& dofmap = *this->dofmap;

  // iterate over all elemsets and build
  // the set of nodes of local elements
  std::set<int> nodset;
  int           proc;
  MPI_Comm_rank(this->getComm(), &proc); proc++;
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
  

  // build array from set of dofs
  int ndofs = dofset.size();
  int *dofs = new int[ndofs];
  std::set<int>::iterator dof_it = dofset.begin();
  for (int i=0; i<ndofs; ) dofs[i++] = *dof_it++;
  dofset.clear();

  // output data
  *_ndofs = ndofs; *_dofs = dofs;
}

void
Problem::preAssemble()
{
  // MPI
  PETSCFEM_COMM_WORLD = this->getComm();
  MPI_Comm_size(PETSCFEM_COMM_WORLD, &SIZE);
  MPI_Comm_rank(PETSCFEM_COMM_WORLD, &MY_RANK);
  // Mesh
  GLOBAL_MESH = this->getMesh();
  // Options
  if (this->options.empty())
    GLOBAL_OPTIONS = OPTIONS::GLOBAL;
  else
    GLOBAL_OPTIONS = this->options;
}

void
Problem::postAssemble()
{
  PETSCFEM_COMM_WORLD = MPI_COMM_NULL;
  GLOBAL_MESH         = NULL;
  GLOBAL_OPTIONS      = OPTIONS::GLOBAL;
}

PYPF_NAMESPACE_END
