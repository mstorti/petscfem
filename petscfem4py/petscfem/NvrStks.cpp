// $Id: NvrStks.cpp,v 1.1.2.1 2006/05/25 00:33:11 dalcinl Exp $

#include "NvrStks.h"

#include <fem.h>
#include <readmesh.h>
#include <arglist.h>
#include <sttfilter.h>

#include <applications/ns/nsi_tet.h>

extern TextHashTable* GLOBAL_OPTIONS;
extern Mesh*          GLOBAL_MESH;
extern int            SIZE, MY_RANK;

static const char COMP_RES[]     = "comp_res";
static const char COMP_RES_JAC[] = "comp_mat_res";

PYPF_NAMESPACE_BEGIN

struct State: public ::State {
  Vec _vs;

  State()                 : ::State(), _vs(PETSC_NULL)
  { this->vec = &this->_vs; }
  State(const State& stt) : ::State(), _vs(PETSC_NULL)
  { this->set(stt._vs, stt.time); }
  State(Vec v, double t)  : ::State(), _vs(PETSC_NULL)
  { this->set(v,t); }
  ~State()
  { if(this->vec == &this->_vs) this->vec = NULL; } // ::State::~State() calls VecDestroy()
  void set(Vec v, double t)
  { this->_vs = v; this->vec = &this->_vs; this->time = t; }
};

struct ArgsNS {
  arg_list       argl;
  State          state0, state1;
  Time           tstar;
  vector<double> hmin;
  WallData*      wall_data;
  GlobParam      glob_param;
  Vec            xst;
  Vec            res;
  Mat            Jac;
  const char*    jobinfo;
  ArgsNS() : argl(),
	     state0(), state1(), tstar(),
	     hmin(1), wall_data(NULL), glob_param(),
	     xst(PETSC_NULL), res(PETSC_NULL), Jac(PETSC_NULL),
	     jobinfo(NULL)
             { }
  ~ArgsNS()  
  {
    this->argl.clear();
    PYPF_PETSC_DESTROY(VecDestroy, this->xst);
    PYPF_DELETE_SCLR(this->wall_data);
  }
  
  void clear() { this->argl.clear(); }

  arg_list&       getArgList()  { return this->argl;    }
  const char*     getJobInfo()  { return this->jobinfo; }
  const TimeData* getTimeData() { return &this->tstar;  }

  void pack(Vec x0, double t0,
	    Vec x1, double t1,
	    Vec r, Mat J,
	    double alpha, bool steady)
  {
    if (x0 == PETSC_NULL) {
      if (!this->xst) PYPF_PETSC_CALL(VecDuplicate(x1, &this->xst));
      x0 = this->xst;
      PYPF_PETSC_CALL(VecCopy(x1, x0));
    } else if (x1 == PETSC_NULL) {
      if (!this->xst) PYPF_PETSC_CALL(VecDuplicate(x0, &this->xst));
      x1 = this->xst;
      PYPF_PETSC_CALL(VecCopy(x0, x1));
    }

    // tstart:  t^* = t^0 + alpha * (t^1 - t^0)
    this->tstar.set(t0 + alpha * (t1-t0));

    // clear arg list
    this->argl.clear();

    // add new state
    this->state1.set(x1, t1);
    this->argl.arg_add(&state1, IN_VECTOR | USE_TIME_DATA);

    // add old state
    this->state0.set(x0, t0);
    this->argl.arg_add(&state0, IN_VECTOR | USE_TIME_DATA);

    // add residual vector and jacobian matrix
    this->res = r;
    this->Jac = J;
    this->argl.arg_add(&this->res, OUT_VECTOR);
    this->jobinfo = COMP_RES;
    if (this->Jac != PETSC_NULL) {
      this->argl.arg_add(&this->Jac, OUT_MATRIX);
      this->jobinfo = COMP_RES_JAC;
    }

    // add hmin (XXX: What is this arg for???)
    this->argl.arg_add(&this->hmin, VECTOR_MIN);

    // setup and add global parameters
    this->glob_param.steady    = steady?1:0;
    this->glob_param.alpha     = alpha;
    this->glob_param.Dt        = t1-t0;;
    this->glob_param.inwt      = 0; // XXX: this is ok?
    this->glob_param.state_old = &this->state0;
    this->glob_param.state     = &this->state1;
    this->argl.arg_add(&this->glob_param, USER_DATA);
  
    // add wall_data
    this->argl.arg_add(this->wall_data, USER_DATA);
  }

};

PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN

NvrStks::~NvrStks()
{
  PYPF_DELETE_SCLR(this->args);
}

NvrStks::NvrStks()
  : Application(),
    alpha(1.0), steady(false),
    args(new ArgsNS)
{ }

NvrStks::NvrStks(const NvrStks& ns)
  : Application(ns),
    alpha(ns.alpha), steady(ns.steady),
    args(new ArgsNS)
{ }

NvrStks::NvrStks(Domain& domain, double alpha, bool steady) 
  : Application(domain),
    alpha(alpha), steady(steady),
    args(new ArgsNS)
{ }


static void
NvrStks_preAssemble(const NvrStks* ns)
{
  const Domain& domain = ns->getDomain();
  // MPI
  PETSCFEM_COMM_WORLD = domain.getComm();
  MPI_Comm_size(PETSCFEM_COMM_WORLD, &SIZE);
  MPI_Comm_rank(PETSCFEM_COMM_WORLD, &MY_RANK);
  // Mesh
  GLOBAL_MESH = domain.getMesh();
  // Options
//   if (this->options.empty())
//     GLOBAL_OPTIONS = OPTIONS::GLOBAL;
//   else
//     GLOBAL_OPTIONS = this->options;
  GLOBAL_OPTIONS      = OPTIONS::GLOBAL;
}

static void
NvrStks_Assemble(Mesh& mesh, DofMap& dofmap, ArgsNS& args) {
  int ierr = ::assemble(mesh, args.argl, dofmap, args.jobinfo, &args.tstar);
  if (ierr) throw Error("PETScFEM assemble error");
}

static void
NvrStks_postAssemble(const NvrStks* ns)
{
  PETSCFEM_COMM_WORLD = MPI_COMM_NULL;
  GLOBAL_MESH         = NULL;
  //GLOBAL_OPTIONS      = OPTIONS::GLOBAL;
}


void 
NvrStks::assemble(Vec x, double t,
		  Vec r, Mat J) const
{
  double alpha  = 1.0;
  bool   steady = true;

  const Mesh&   mesh   = this->domain->getMesh();
  const DofMap& dofmap = this->domain->getDofMap();
  ArgsNS*       args   = this->args;
  
  // XXX test vector and matrix sizes against dofmap !!!
  args->pack(PETSC_NULL, t, x, t, r, J, alpha, steady);
  NvrStks_preAssemble(this);
  int ierr = ::assemble(mesh, args->argl, dofmap, args->jobinfo, &args->tstar);
  NvrStks_postAssemble(this);
  if (ierr) throw Error("NvrStks::operator(): assemble error");

  PYPF_PETSC_CALL(VecScale(r, -1.0/alpha));
}


void
NvrStks::assemble(Vec x0, double t0, Vec x1, double t1, Vec r, Mat J) const
{
  double alpha  = this->alpha;
  bool   steady = this->steady;

  if (steady) alpha = 1.0;
  else if (alpha <= 0.0) throw Error("invalid value, 'alpha' <= 0.0");
  else if (alpha >  1.0) throw Error("invalid value, 'alpha' > 1.0");
  if (t1 < t0 )          throw Error("invalid value, 't1' < 't0'");
    
  const Mesh&   mesh   = this->domain->getMesh();
  const DofMap& dofmap = this->domain->getDofMap();
  ArgsNS*       args   = this->args;
  
  // XXX test vector and matrix sizes against dofmap !!!
  args->pack(x0, t0, x1, t1, r, J, alpha, steady);
  NvrStks_preAssemble(this);
  int ierr = ::assemble(mesh, args->argl, dofmap, args->jobinfo, &args->tstar);
  NvrStks_postAssemble(this);
  if (ierr) throw Error("NvrStks::operator(): assemble error");

  PYPF_PETSC_CALL(VecScale(r, -1.0/alpha));
 
}

PYPF_NAMESPACE_END
