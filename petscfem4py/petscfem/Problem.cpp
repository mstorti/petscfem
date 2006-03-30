// $Id: Problem.cpp,v 1.1.2.5 2006/03/30 15:18:14 rodrigop Exp $

#include "Problem.h"

#include <fem.h>
#include <dofmap.h>
#include <readmesh.h>
//#include <getprop.h>
//#include <utils.h>
//#include <util2.h>
//#include <sttfilter.h>
//#include <pfmat.h>
//#include <hook.h>
extern Mesh* GLOBAL_MESH;
int my_mesh_part(MPI_Comm comm, Mesh* mesh, Dofmap* dofmap);


PYPF_NAMESPACE_BEGIN


Problem::~Problem() 
{ 
  this->mesh->decref();
  this->dofmap->decref();
}

// Problem::Problem()
//   : Object(),
//     nnod(0), ndim(0), ndof(0),
//     mesh(new Mesh), dofmap(new DofMap),
//     setupcalled(false)
// { 
//   this->mesh->incref();
//   this->dofmap->incref();
// }

// Problem::Problem(const Problem& p) 
//   : Object(p),
//     nnod(p.nnod), ndim(p.ndim), ndof(p.ndof),
//     mesh(p.mesh), dofmap(p.dofmap),
//     setupcalled(p.setupcalled)
// { 
//   this->mesh->incref();
//   this->dofmap->incref();
// }


Problem::Problem(Mesh* mesh, DofMap* dofmap)
  : Object(),
    nnod(mesh->nodedata->nnod), ndim(mesh->nodedata->ndim), ndof(dofmap->ndof),
    mesh(mesh), dofmap(dofmap),
    setupcalled(false)
{
  this->mesh->incref();
  this->dofmap->incref();
}


static void
problem_setup(MPI_Comm comm, Mesh& mesh, DofMap& dofmap)
{
  GLOBAL_MESH = mesh;
  my_mesh_part(comm, mesh, dofmap);
  GLOBAL_MESH = NULL;
}

void
Problem::setUp()
{
  if (this->setupcalled) return;
  
  problem_setup(this->comm, *this->mesh, *this->dofmap);
  this->dofmap->frozen = true;
  
  this->setupcalled = true;
}


Mesh*
Problem::getMesh() const
{
  return this->mesh;
}

DofMap*
Problem::getDofMap() const
{
  return this->dofmap;
}

void 
Problem::getDofSizes (int* local, int* global) const 
{
  this->dofmap->getSizes(local, global);
}

void 
Problem::getDofRange (int* first, int* last) const 
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
state2solution(DofMap& dofmap, const double state[], double solution[]) 
{
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  for (int i=0; i<nnod; i++)
    for (int j=0; j<ndof; j++)
      dofmap->get_nodal_value(i+1, j+1, state, NULL, solution[i*ndof+j]);
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

  state2solution(dofmap, stt_buff, sol_buff);

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


// void Problem::read(const std::string& filename)
// {

//   char* fcase = const_cast<char*>(filename.c_str());
//   int neq, size, rank;
//   MPI_Comm_size(PETSC_COMM_WORLD, &size);
//   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

//   Mesh::Base*   mesh   = NULL;
//   DofMap::Base* dofmap = NULL;;

//   read_mesh(mesh, fcase, dofmap, neq, size, rank);

//   this->comm = PETSC_COMM_WORLD;

//   this->nnod = mesh->nodedata->nnod;
//   this->ndim = mesh->nodedata->ndim;
//   this->ndof = dofmap->ndof;

//   this->mesh->decref();
//   this->mesh = new Mesh(mesh);
//   this->mesh->incref();

//   this->dofmap->decref();
//   this->dofmap = new DofMap(dofmap);
//   this->dofmap->incref();

//   this->setupcalled = true;

// }



PYPF_NAMESPACE_END
