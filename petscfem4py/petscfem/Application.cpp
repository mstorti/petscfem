// $Id: Application.cpp,v 1.1.2.3 2006/05/30 20:14:04 dalcinl Exp $

#include "Application.h"

#include <fem.h>
#include <dofmap.h>

extern TextHashTable* GLOBAL_OPTIONS;
extern Mesh*          GLOBAL_MESH;
extern int            SIZE, MY_RANK;

PYPF_NAMESPACE_BEGIN

Application::~Application() 
{ 
  PYPF_DECREF(this->domain);
}

Application::Application()
  : Object(),
    domain(NULL)
{ }

Application::Application(const Application& app)
  : Object(app),
    domain(app.domain)
{
  PYPF_INCREF(this->domain);
}

Application::Application(Domain& domain)
  : Object(domain.getComm()),
    domain(&domain)
{
  PYPF_INCREF(this->domain);
}

Domain&
Application::getDomain() const
{
  return *this->domain;
}

// #undef __FUNCT__  
// #define __FUNCT__ "Application_getGhostValues"
// static PetscErrorCode
// Application_getGhostValues(const Application& app,
// 			   const double state[],
// 			   std::vector<double>& buff)
// {
//   const DofMap& dofmap = app.getDomain().getDofMap();
//   int ldofsz, gdofsz; dofmap.getSizes(&ldofsz, &gdofsz);
//   int  nghost = dofmap->ghost_dofs->size();
//   buff.reserve(nghost+1); buff.resize(nghost);
//   double* ghosts = &buff[0];
//   Vec sttvec, ghtvec;

//   PetscFunctionBegin;
//   VecCreateMPIWithArray(dofmap.getComm(),ldofsz,gdofsz,state,&sttvec);
//   VecCreateSeqWithArray(PETSC_COMM_SELF,nghost,ghosts,&ghtvec);
//   VecScatter scatter = *(dofmap->ghost_scatter);
//   VecScatterBegin(sttvec,ghtvec,INSERT_VALUES,SCATTER_FORWARD,scatter);
//   VecScatterEnd  (sttvec,ghtvec,INSERT_VALUES,SCATTER_FORWARD,scatter);
//   VecDestroy(sttvec);
//   VecDestroy(ghtvec);
//   PetscFunctionReturn(0);
// }

void
Application::getNodalValues(const double state[], double time,
			    double values[]) const
{
  int nn, nf;
  this->domain->getSizes(&nn,&nf);
  const DofMap::Base* dm = this->domain->getDofMap();
  Time time_data = time;
  for (int i=0; i<nn; i++)
    for (int j=0; j<nf; j++)
      dm->get_nodal_value(i+1, j+1, 
			  state, &time_data, 
			  values[i*nf+j]);

}

void
Application::getNodalValues(const double state[], double time,
			    int nn, const int nodes[],
			    int nf, const int fields[],
			    double values[]) const
{
  const DofMap::Base* dm = this->domain->getDofMap();
  Time time_data = time;
  for (int i=0; i<nn; i++)
    for (int j=0; j<nf; j++)
      dm->get_nodal_value(nodes[i]+1, fields[j]+1,
			  state, &time_data, 
			  values[i*nf+j]);
}

void 
Application::assemble(const Application& app, const ArgList& args)
{
  const Domain& domain = app.getDomain();
  const Mesh&   mesh   = domain.getMesh();
  const DofMap& dofmap = domain.getDofMap();
  
  // pre-assemble
  PETSCFEM_COMM_WORLD = domain.getComm();
  MPI_Comm_size(PETSCFEM_COMM_WORLD, &SIZE);
  MPI_Comm_rank(PETSCFEM_COMM_WORLD, &MY_RANK);
  GLOBAL_MESH = mesh;
  if (app.options.empty())
    GLOBAL_OPTIONS = OPTIONS::GLOBAL;
  else
    GLOBAL_OPTIONS = app.options;

  // assemble
  int ierr = ::assemble(mesh, args, dofmap, args.job(), args.time());

  // post-assemble
  PETSCFEM_COMM_WORLD = MPI_COMM_NULL;
  GLOBAL_MESH         = NULL;
  GLOBAL_OPTIONS      = OPTIONS::GLOBAL;

  if (ierr) throw Error("PETScFEM assemble error");

}

PYPF_NAMESPACE_END
