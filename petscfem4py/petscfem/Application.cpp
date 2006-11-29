// $Id: Application.cpp,v 1.1.2.7 2006/11/29 22:35:09 dalcinl Exp $

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


void 
Application::allocateState(Vec& r, const std::string& vec_type) const
{
  MPI_Comm comm;
  comm = this->getComm();
  PetscInt n, N;
  this->domain->getDofSizes(&n, &N);

  if (r != PETSC_NULL) throw Error("state vector already allocated");
  // create
  PYPF_PETSC_CALL(VecCreate, (comm, &r));
  PYPF_PETSC_CALL(VecSetSizes, (r, n, N));
  // set type
  const char *vtype = NULL;
  if (vec_type.size() > 0)
    vtype = vec_type.c_str();
  else {
    int size; MPI_Comm_size(comm, &size);
    vtype = (size == 1) ? VECSEQ : VECMPI;
  }
  PYPF_PETSC_CALL(VecSetType, (r, vtype));
}

void 
Application::allocateResidual(Vec& r, const std::string& vec_type) const
{
  MPI_Comm comm;
  PetscInt n, N;
  comm = this->getComm();
  this->domain->getDofSizes(&n, &N);

  if (r != PETSC_NULL) throw Error("residual vector already allocated");
  // create
  PYPF_PETSC_CALL(VecCreate, (comm, &r));
  PYPF_PETSC_CALL(VecSetSizes, (r, n, N));
  // set type
  const char *vtype = NULL;
  if (vec_type.size() > 0)
    vtype = vec_type.c_str();
  else {
    int size; MPI_Comm_size(comm, &size);
    vtype = (size == 1) ? VECSEQ : VECMPI;
  }
  PYPF_PETSC_CALL(VecSetType, (r, vtype));
}

void 
Application::allocateJacobian(Mat& J, const std::string& mat_type) const
{
  MPI_Comm comm;
  comm = this->getComm();
  PetscInt n, N;
  this->domain->getDofSizes(&n, &N);
  
  if (J != PETSC_NULL) throw Error("jacobian matrix already allocated");
  // create
  PYPF_PETSC_CALL(MatCreate, (comm, &J));
  PYPF_PETSC_CALL(MatSetSizes, (J, n, n, N, N));
  // set type
  const char *mtype = NULL;
  if (mat_type.size() > 0)
    mtype = mat_type.c_str();
  else 
    mtype = MATAIJ;
  PYPF_PETSC_CALL(MatSetType,(J, mtype));
}

void 
Application::allocateSolution(Vec& r) const
{
  MPI_Comm comm = PETSC_COMM_SELF;
  PetscInt nnod, ndof, n, N;
  this->domain->getSizes(&nnod, &ndof);
  n = N = nnod*ndof;

  if (r != PETSC_NULL) throw Error("solution vector already allocated");
  // create
  PYPF_PETSC_CALL(VecCreate, (comm, &r));
  PYPF_PETSC_CALL(VecSetSizes, (r, n, N));
  // set type
  const char *vtype = VECSEQ;
  PYPF_PETSC_CALL(VecSetType, (r, vtype));
}

void 
Application::buildState(Vec solution, Vec state)
{
  int nnod,  ndof;  this->domain->getSizes(&nnod,&ndof);
  int ldofs, gdofs; this->domain->getDofSizes(&ldofs,&gdofs);

  /* check provided state vector */
  PetscTruth stt_valid; VecValid(state,&stt_valid);
  PYPF_ASSERT(stt_valid == PETSC_TRUE, "provided state vector is not valid");
  PetscInt n, N; VecGetLocalSize(state,&n); VecGetSize(state,&N);
  PYPF_ASSERT(ldofs == n,"provided state vector has wrong local size");
  PYPF_ASSERT(gdofs == N,"provided state vector has wrong global size");

  MPI_Comm comm;  PetscObjectGetComm((PetscObject)state,&comm);
  int      rank;  MPI_Comm_rank(comm,&rank);

  int root = 0;
  
  /* check provided solution vector */
  if (rank==root) {
    PetscTruth sol_valid; VecValid(solution,&sol_valid);
    PYPF_ASSERT(sol_valid == PETSC_TRUE, "provided solution vector is not valid");
    PetscInt sol_size; VecGetLocalSize(solution,&sol_size);
    PYPF_ASSERT(sol_size == nnod*ndof,   "provided solution vector has wrong local size ");
  }
  
  /* create scatter */
  if (!this->scatter) PYPF_PETSC_CALL(VecScatterCreateToAll, (state,&this->scatter,&this->state));
  
  /* build state in root processor */
  //VecZeroEntries(this->state);
  PetscScalar* stt_array;
  VecGetArray(this->state,&stt_array);
  if (rank==root) {
    PetscScalar* sol_array; 
    VecGetArray(solution,&sol_array);
#if 1
    DofMap::Base* dm = this->domain->getDofMap();
    //std::vector<PetscScalar> sttbuff(dm->neqtot);
    PetscScalar* sttbuff; 
    PetscMalloc(dm->neqtot*sizeof(PetscScalar), &sttbuff);
    dm->solve(&sttbuff[0],sol_array);
    PetscMemcpy(stt_array, &sttbuff[0], dm->neq*sizeof(PetscScalar));
    PetscFree(sttbuff);
#else
    const DofMap& dofmap = this->domain->getDofMap();
    dofmap.solve(nnod*ndof,sol_array,gdofs,stt_array);
#endif
    VecRestoreArray(solution,&sol_array);
    PetscObjectStateDecrease((PetscObject)solution);
  } else {
    PetscMemzero(stt_array,gdofs*sizeof(PetscScalar));
  }
  VecRestoreArray(this->state,&stt_array);

  /* scatter state values to all processors */
  VecZeroEntries(state);
  PYPF_PETSC_CALL(VecScatterBegin, (this->state,state,ADD_VALUES,SCATTER_REVERSE,this->scatter));
  PYPF_PETSC_CALL(VecScatterEnd,   (this->state,state,ADD_VALUES,SCATTER_REVERSE,this->scatter));


}

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
#if 0
      const DofMap::Base* dm = this->domain->getDofMap();
      const Time time_data = time;
      for (int i=0; i<nnod; i++)
	for (int j=0; j<ndof; j++)
	  dm->get_nodal_value(i+1, j+1,stt_array,&time_data,sol_array[i*ndof+j]);
#else
      const DofMap& dofmap = this->domain->getDofMap();
      dofmap.apply(time,gdofs,stt_array,nnod*ndof,sol_array);
#endif
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
    GLOBAL_OPTIONS = Options::GLOBAL;
  else
    GLOBAL_OPTIONS = app.options;

  // assemble
  int ierr = ::assemble(mesh, args, dofmap, args.job(), args.time());

  // post-assemble
  PETSCFEM_COMM_WORLD = MPI_COMM_NULL;
  GLOBAL_MESH         = NULL;
  GLOBAL_OPTIONS      = Options::GLOBAL;

  if (ierr) throw Error("PETScFEM assemble error");

}

PYPF_NAMESPACE_END
