// $Id: Problem.cpp,v 1.1.2.4 2006/03/28 22:13:25 rodrigop Exp $

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

Problem::Problem()
  : Object(),
    nnod(0), ndim(0), ndof(0),
    mesh(new Mesh), dofmap(new DofMap),
    setupcalled(false)
{ 
  this->mesh->incref();
  this->dofmap->incref();
}

Problem::Problem(const Problem& p) 
  : Object(p),
    nnod(p.nnod), ndim(p.ndim), ndof(p.ndof),
    mesh(p.mesh), dofmap(p.dofmap),
    setupcalled(p.setupcalled)
{ 
  this->mesh->incref();
  this->dofmap->incref();
}


Problem::Problem(Mesh* mesh, DofMap* dofmap)
  : Object(),
    nnod(0), ndim(0), ndof(0),
    mesh(mesh), dofmap(dofmap),
    setupcalled(false)
{
  this->mesh->incref();
  this->dofmap->incref();
  
  this->nnod = (*this->mesh)->nodedata->nnod;
  this->ndim = (*this->mesh)->nodedata->ndim;
  this->ndof = (*this->dofmap)->ndof;
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

void Problem::read(const std::string& filename)
{

  char* fcase = const_cast<char*>(filename.c_str());
  int neq, size, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  Mesh::Base*   mesh   = NULL;
  DofMap::Base* dofmap = NULL;;

  read_mesh(mesh, fcase, dofmap, neq, size, rank);

  this->comm = PETSC_COMM_WORLD;

  this->nnod = mesh->nodedata->nnod;
  this->ndim = mesh->nodedata->ndim;
  this->ndof = dofmap->ndof;

  this->mesh->decref();
  this->mesh = new Mesh(mesh);
  this->mesh->incref();

  this->dofmap->decref();
  this->dofmap = new DofMap(dofmap);
  this->dofmap->incref();

  this->setupcalled = true;

}


static void
solution2state(Dofmap *dofmap, double solution[], double state[]) 
{
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  int neq  = dofmap->neq;
  int neqt = dofmap->neqtot;

  std::vector<double> solvec; solvec.resize(nnod*ndof);
  std::vector<double> sttvec; sttvec.resize(neqt);

  PetscMemcpy(&solvec[0], solution, neqt*sizeof(double));
  dofmap->solve(&sttvec[0], &solvec[0]);
  PetscMemcpy(state, &sttvec[0], neq*sizeof(double));
}

void
Problem::buildState(Vec solution, Vec state)
{

  DofMap::Base* dofmap = *this->dofmap;
  
  double* sol_buff;
  double* stt_buff;
  VecGetArray(solution, &sol_buff);
  VecGetArray(state, &stt_buff);
  
  solution2state(dofmap, sol_buff, stt_buff);

  VecRestoreArray(solution, &sol_buff);
  VecRestoreArray(state,    &stt_buff);
  PetscObjectStateDecrease(solution);
}


static void
state2solution(Dofmap *dofmap, double state[], double solution[]) 
{
  int nnod = dofmap->nnod;
  int ndof = dofmap->ndof;
  for (int k=1; k<=nnod; k++) {
    for (int kldof=1; kldof<=ndof; kldof++) {
      double dval;
      dofmap->get_nodal_value(k, kldof, state, NULL, dval);
      solution[ndof*(k-1)+kldof-1] = dval;
    }
  }
}

void
Problem::buildSolution(Vec state, Vec solution)
{
  DofMap::Base* dofmap = *this->dofmap;

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
Problem::setUp()
{
#if 0
  /* setup */
  if (this->setupcalled) return;
  /* ------ */
  this->mesh->setUp();
  this->dofmap->setUp();
  /* ------ */
  this->setupcalled = true;
#else

  if (this->setupcalled) return;

  Mesh::Base*   mesh    = *this->mesh;
  DofMap::Base* dofmap  = *this->dofmap;

  GLOBAL_MESH = mesh;

  my_mesh_part(this->comm, mesh, dofmap);

  this->dofmap->frozen = true;

  this->setupcalled = true;
  //for (int i=0; i<this->mesh->elemsetlist.size(); 
  //     (*this->mesh->elemsetlist[i++])->initialize());
#endif

}

void
Problem::clear()
{
  /* mesh */
  this->mesh->decref();
  this->mesh = new Mesh();
  this->mesh->incref();
  /* dofmap */
  this->dofmap->decref();
  this->dofmap = new DofMap();
  this->dofmap->incref();
}



PYPF_NAMESPACE_END
