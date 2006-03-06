// $Id: Problem.cpp,v 1.1.2.2 2006/03/06 16:56:04 rodrigop Exp $

#include "Problem.h"

#include <fem.h>
#include <dofmap.h>
#include <readmesh.h>
#include <getprop.h>
#include <utils.h>
#include <util2.h>
//#include <sttfilter.h>
//#include <pfmat.h>
//#include <hook.h>

PyPF::Problem::~Problem() 
{ 
  Mesh::Base* mesh = this->mesh;
  this->mesh = NULL;
  delete mesh;
  
  DofMap::Base* dofmap = this->dofmap;
  this->dofmap = NULL;
  delete dofmap;
}

PyPF::Problem::Problem() 
  : comm(MPI_COMM_NULL),
    nnod(0), ndim(0), ndof(0),
    mesh(), dofmap(),
    setupcalled(false)
{ }

PyPF::Problem::Problem(int nnod, int ndim, int ndof)
{ 
  MPI_Comm comm = PETSC_COMM_WORLD;

  // check arguments
  if (comm == MPI_COMM_NULL) throw Error("null communicator");
  if (nnod < 1)              throw Error("invalid number of nodes");
  if (ndim < 1)              throw Error("invalid number of dimensions");
  if (ndim > 3)              throw Error("invalid number of dimensions");
  if (ndof < 1)              throw Error("invalid number of degree of fredom");
  
  int comm_size, comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  // build mesh
  Mesh::Base* mesh = new Mesh::Base;
  mesh->global_options = new OptionTable;
  mesh->nodedata       = new Nodedata::Base;
  mesh->elemsetlist    = da_create(sizeof(Elemset::Base *));
  
  // build nodedata
  double* xyz     = new double[nnod*ndim];
  memset(xyz, 0, nnod*ndim*sizeof(double));
  Nodedata::Base* ndata = mesh->nodedata;
  ndata->options  = new OptionTable;
  ndata->nodedata = xyz;
  ndata->nnod     = nnod;
  ndata->ndim     = ndim;
  ndata->nu       = ndim;

  // build dofmap
  float* tpwgts  = new float[comm_size];
  for (int i=0; i<comm_size; i++) 
    tpwgts[i] = 1.0/float(comm_size);
  DofMap::Base* dofmap = new DofMap::Base;
  dofmap->nnod   = nnod;
  dofmap->ndof   = ndof;
  dofmap->size   = comm_size;
  dofmap->tpwgts = tpwgts;
  dofmap->id     = new ::idmap(nnod*ndof, NULL_MAP);
  for (int i=1; i<=nnod; i++)
    for (int j=1; j<=ndof; j++) {
      int idx = dofmap->edof(i, j);
      dofmap->id->set_elem(idx, idx, 1.0);
    }
  
  // set objects
  this->nnod   = nnod;
  this->ndim   = ndim;
  this->ndof   = ndof;
  this->comm   = comm;
  this->mesh   = mesh;
  this->dofmap = dofmap;

  // setup flag
  this->setupcalled = false;
}

void PyPF::Problem::fromFile(const std::string& filename) 
{
  if (this->mesh) { 
    Mesh::Base* mesh = this->mesh;
    this->mesh = NULL;
    delete mesh;
  }
  if (this->dofmap) { 
    DofMap::Base* dofmap = this->dofmap;
    this->dofmap = NULL;
    delete dofmap;
  }

  char* fcase = const_cast<char*>(filename.c_str());
  int neq, size, rank;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  read_mesh(this->mesh, fcase, this->dofmap, neq, size, rank);

  this->comm = PETSC_COMM_WORLD;
  this->nnod = this->mesh->nodedata->nnod;
  this->ndim = this->mesh->nodedata->ndim;
  this->ndof = this->dofmap->ndof;
  this->setupcalled = true;
}

MPI_Comm&
PyPF::Problem::getComm() 
{
  return this->comm;
}

// PyPF::Nodedata
// PyPF::Problem::getNodeData() 
// {
//   if (!this->mesh) throw Error("null mesh");
//   if (!this->mesh->nodedata) throw Error("null nodedata in mesh");
//   return this->mesh->nodedata;
// }

PyPF::Mesh
PyPF::Problem::getMesh()
{
  /* test */
  if (!this->mesh) throw Error("null Mesh");
  
  return this->mesh;
}

PyPF::DofMap
PyPF::Problem::getDofMap()
{
  /* test */
  if (!this->dofmap) throw Error("null DofMap");
  
  return this->dofmap;
}

void
PyPF::Problem::setUp()
{
  /* test */
  if (this->comm == MPI_COMM_NULL) throw Error("null communicator");
  if (!this->mesh)                 throw Error("null mesh");
  if (!this->mesh->nodedata)       throw Error("null node data");
  if (!this->mesh->elemsetlist)    throw Error("null elemset list");
  if (!this->dofmap)               throw Error("null dofmap");
  
  /* setup */
  if (this->setupcalled) return;
  /* ------ */

  /* ------ */
  this->setupcalled = true;
}

