// $Id: Application.cpp,v 1.1.2.4 2006/06/22 22:35:42 dalcinl Exp $

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
  PYPF_PETSC_DESTROY(VecScatterDestroy, this->scatter);
  PYPF_PETSC_DESTROY(VecDestroy,        this->state);
}

Application::Application()
  : Object(),
    domain(NULL),
    scatter(PETSC_NULL), state(PETSC_NULL)
{ }

Application::Application(const Application& app)
  : Object(app),
    domain(app.domain),
    scatter(PETSC_NULL), state(PETSC_NULL)
{
  PYPF_INCREF(this->domain);
}

Application::Application(Domain& domain)
  : Object(domain.getComm()),
    domain(&domain),
    scatter(PETSC_NULL), state(PETSC_NULL)
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
Application::buildSolution(double time, Vec state, Vec solution)
{
  int nnod,  ndof;  this->domain->getSizes(&nnod,&ndof);
  int ldofs, gdofs; this->domain->getDofSizes(&ldofs,&gdofs);

  /* check provided state vector */
  PetscTruth stt_valid; VecValid(state,&stt_valid);
  PYPF_ASSERT(stt_valid == PETSC_TRUE, "provided state vector is not valid");
  PetscInt n, N; VecGetLocalSize(state,&n); VecGetSize(state,&N);
  PYPF_ASSERT(ldofs == n,"provided state vector has wrong local size");
  PYPF_ASSERT(gdofs == N,"provided state vector has wrong global size");
  
  /* scatter state values to all processors */
  if (!this->scatter) PYPF_PETSC_CALL(VecScatterCreateToAll, (state,&this->scatter,&this->state));
  PYPF_PETSC_CALL(VecScatterBegin, (state,this->state,INSERT_VALUES,SCATTER_FORWARD,this->scatter));
  PYPF_PETSC_CALL(VecScatterEnd,   (state,this->state,INSERT_VALUES,SCATTER_FORWARD,this->scatter));
  
  /* user does not want solution on this processor */
  if (solution == PETSC_NULL) return;
  
  /* check provided solution vector */
  PetscTruth sol_valid; VecValid(solution,&sol_valid);
  PYPF_ASSERT(sol_valid == PETSC_TRUE, "provided solution vector is not valid");
  PetscInt sol_size; VecGetLocalSize(solution,&sol_size);
  PYPF_ASSERT(sol_size == 0 || sol_size == nnod*ndof, 
	      "provided solution vector has wrong local size ");
  
  PetscScalar* sol_array;
  VecGetArray(solution, &sol_array);
  if (sol_size > 0 && sol_array != NULL) {
    PetscScalar* stt_array;
    VecGetArray(this->state, &stt_array);
    if (sol_array != NULL) {
      const DofMap::Base* dm = this->domain->getDofMap();
      const Time time_data = time;
      for (int i=0; i<nnod; i++)
	for (int j=0; j<ndof; j++)
	  dm->get_nodal_value(i+1, j+1,stt_array,&time_data,sol_array[i*ndof+j]);
    }
    /* array of state vector was not touched */
    VecRestoreArray(this->state, &stt_array);
    PetscObjectStateDecrease((PetscObject)this->state);
  }
  VecRestoreArray(solution, &sol_array);
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
