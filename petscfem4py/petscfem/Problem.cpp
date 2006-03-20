// $Id: Problem.cpp,v 1.1.2.3 2006/03/20 16:06:00 rodrigop Exp $

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


PYPF_NAMESPACE_BEGIN


OptionTable*
Problem::get_opt_table() const
{
  throw Error("no options for Problem");
  return NULL;
}

Problem::~Problem() 
{ 
  this->mesh->decref();
  this->dofmap->decref();
}

Problem::Problem() 
  : comm(MPI_COMM_NULL),
    nnod(0), ndim(0), ndof(0),
    mesh(new Mesh), dofmap(new DofMap),
    setupcalled(false)
{ 
  this->mesh->incref();
  this->dofmap->incref();
}

Problem::Problem(int nnod, int ndim, int ndof) 
  : comm(PETSC_COMM_WORLD), 
    nnod(0), ndim(0), ndof(0),
    mesh(new Mesh), dofmap(new DofMap),
    setupcalled(false)

{ 
//   MPI_Comm comm = this->comm;

//   // check arguments
//   if (comm == MPI_COMM_NULL) throw Error("null communicator");
//   if (nnod < 1)              throw Error("invalid number of nodes");
//   if (ndim < 1)              throw Error("invalid number of dimensions");
//   if (ndim > 3)              throw Error("invalid number of dimensions");
//   if (ndof < 1)              throw Error("invalid number of degree of fredom");
  
//   int comm_size, comm_rank;
//   MPI_Comm_size(comm, &comm_size);
//   MPI_Comm_rank(comm, &comm_rank);

//   // build mesh
//   Mesh::Base* mesh = new Mesh::Base;
//   mesh->global_options = new OptionTable;
//   mesh->nodedata       = new Nodedata::Base;
//   mesh->elemsetlist    = da_create(sizeof(Elemset::Base*));
  
//   // build nodedata
//   double* xyz     = new double[nnod*ndim];
//   memset(xyz, 0, nnod*ndim*sizeof(double));
//   Nodedata::Base* ndata = mesh->nodedata;
//   ndata->options  = new OptionTable;
//   ndata->nodedata = xyz;
//   ndata->nnod     = nnod;
//   ndata->ndim     = ndim;
//   ndata->nu       = ndim;

//   // build dofmap
//   float* tpwgts  = new float[comm_size];
//   for (int i=0; i<comm_size; i++) 
//     tpwgts[i] = 1.0/float(comm_size);
//   DofMap::Base* dofmap = new DofMap::Base;
//   dofmap->comm   = comm;
//   dofmap->nnod   = nnod;
//   dofmap->ndof   = ndof;
//   dofmap->size   = comm_size;
//   dofmap->tpwgts = tpwgts;
//   dofmap->id     = new ::idmap(nnod*ndof, NULL_MAP);
//   for (int i=1; i<=nnod; i++)
//     for (int j=1; j<=ndof; j++) {
//       int idx = dofmap->edof(i, j);
//       dofmap->id->set_elem(idx, idx, 1.0);
//     }
  
//   // set objects
//   this->nnod   = nnod;
//   this->ndim   = ndim;
//   this->ndof   = ndof;
//   this->comm   = comm;
//   this->mesh   = mesh;
//   this->dofmap = dofmap;

//   // setup flag
//   this->setupcalled = false;
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

void Problem::fromFile(const std::string& filename) 
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
  /* setup */
  if (this->setupcalled) return;
  /* ------ */

  /* ------ */
  this->setupcalled = true;
}

void
Problem::clear()
{
  /* mesh */
  this->mesh->decref();
  this->mesh = new Mesh(mesh);
  this->mesh->incref();
  /* dofmap */
  this->dofmap->decref();
  this->dofmap = new DofMap(dofmap);
  this->dofmap->incref();
}



PYPF_NAMESPACE_END
