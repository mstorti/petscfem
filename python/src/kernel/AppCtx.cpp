// $Id$

#include "AppCtx.h"
#include "Domain.h"

#include <fem.h>
#include <dofmap.h>

// ---------------------------------------------------------------- //

PF4PY_NAMESPACE_BEGIN

AppCtx::Args::~Args()
{ }

AppCtx::Args::Args()
  : impl(new AppCtx::Args::Impl)
{ }

PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN

AppCtx::~AppCtx()
{ }

AppCtx::AppCtx(const AppCtx& app)
  : args(0)
{ }

AppCtx::AppCtx()
  : args(0)
{ }

PF4PY_NAMESPACE_END


#include "gvars.h"

PF4PY_NAMESPACE_BEGIN

void 
AppCtx::assemble(const Domain* domain, const AppCtx::Args* args)
{ 
  assert(domain != NULL);
  assert(args   != NULL);

  Comm           C = domain->getComm();
  const Mesh&    M = domain->getMesh();
  const Dofset&  D = domain->getDofset();
  const Options& O = domain->getOptions();
  
  MPI_Comm       comm   = C;
  Mesh::Impl*    mesh   = M;
  Dofset::Impl*  dofmap = D;
  Options::Impl* thash  = NULL;
  if (comm   == MPI_COMM_NULL) throw Error("Domain: null comunicator");
  if (mesh   == NULL)          throw Error("Domain: mesh not ready");
  if (dofmap == NULL)          throw Error("Domain: dofmap not ready");
  if (!O.empty()) thash = O;

  GlobalVars gvars(comm, thash, mesh, dofmap);

  Args::Impl* argl = args->argl(); assert(argl != NULL);
  const char* job  = args->job();  assert(job  != NULL);
  const Time* time = args->time(); assert(time != NULL);
  
  int ierr = ::assemble(mesh, *argl, dofmap, job, time);
  if (ierr) throw Error(ierr, "AppCtx: "
			"PETSc-FEM error in assemble");
}

PF4PY_NAMESPACE_END

// ---------------------------------------------------------------- //

/* 
   the following is a vile hack to avoid residual uploading 
*/

typedef void (*VecOp)(void);
#define VECOP_SET_VALUES ((VecOperation)19)
static PetscErrorCode vec_setvals_empty
(Vec vec,PetscInt n,const PetscInt i[],const PetscScalar v[],InsertMode addv)
{ return 0; }
#undef __FUNCT__
#define __FUNCT__ "VecCreateEmpty"
static PetscErrorCode VecCreateEmpty(MPI_Comm comm, Vec *empty) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecCreateSeqWithArray(comm,0,PETSC_NULL,empty);CHKERRQ(ierr);
  ierr = VecSetOperation(*empty,VECOP_SET_VALUES,(VecOp)vec_setvals_empty);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* 
   the following is a vile hack to avoid jacobian uploading 
*/

typedef void (*MatOp)(void);
static PetscErrorCode mat_setvals_empty
(Mat mat,PetscInt m,const PetscInt i[],PetscInt n,const PetscInt j[],const PetscScalar v[],InsertMode addv)
{ return 0; }
#undef __FUNCT__
#define __FUNCT__ "MatCreateEmpty"
static PetscErrorCode MatCreateEmpty(MPI_Comm comm, Mat *empty) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatCreateShell(comm,0,0,0,0,PETSC_NULL,empty);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*empty,MATOP_SET_VALUES,(MatOp)mat_setvals_empty);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ---------------------------------------------------------------- //

#include "sttfilter.h"

PF4PY_NAMESPACE_BEGIN

struct State : public ::State {
  Vec _vs;
  ~State()  // ::State::~State() calls VecDestroy() !
  { if(this->vec == &this->_vs) this->vec = NULL; }
  State()                 
    : ::State(), _vs(PETSC_NULL)
  { this->vec = &this->_vs; }
  State(const State& stt) 
    : ::State(), _vs(PETSC_NULL)
  { this->set(stt.time, stt._vs); }
  State(double t, Vec x) 
    : ::State(), _vs(PETSC_NULL)
  { this->set(t,x); }
  State(const std::pair<double,Vec>& stt) 
    : ::State(), _vs(PETSC_NULL)
  { this->set(stt.first, stt.second); }
  void set(double t, Vec x)
  { this->time = t; this->_vs = x; this->vec = &this->_vs; }
};

PF4PY_NAMESPACE_END

// ---------------------------------------------------------------- //

#include <pfmat.h>
#include <pfptscmat.h>

PF4PY_NAMESPACE_BEGIN

struct ProfileBase : public PFPETScMat {
  ProfileBase(const Dofset &ds)
    : PFPETScMat(ds.getSizes().second,
		 *((::Dofmap*)ds),
		 ds.getComm()) { }
  int size(int j)                                    { return mat_size; }
  int set_value_a(int, int, PetscScalar, InsertMode) { return 0; }
  int assembly_begin_a(MatAssemblyType type)         { return 0; }
  int assembly_end_a(MatAssemblyType type)           { return 0; }
  int factor_and_solve_a(Vec&, Vec&)                 { return 0; }
  int solve_only_a(Vec&, Vec&)                       { return 0; }
  int clean_mat_a()                                  { return 0; }
  int clean_factor_a()                               { return 0; }
  int view(PetscViewer viewer)                       { return 0; }

};

struct Profile : public ProfileBase {
  const Dofset* dofset;
  std::vector<int> buf1;
  std::vector<int> buf2;
  std::vector<int> buf3;
  Profile(const Dofset &ds)
    : ProfileBase(ds),
      dofset(&ds),
      buf1(), buf2(), buf3() { }
  virtual void pre_build() {}
  virtual void build()=0;
  virtual void post_build() {}
  int create_a() 
  {
    this->buf1.clear();
    this->buf2.clear();
    this->buf3.clear();
    
    this->pre_build();
    lgraph->scatter();
    this->build();
    lgraph->clear();
    this->post_build();

    return 0;
  }
};

struct AdjGraph : public StoreGraph {
  std::pair<int,int>           range;
  std::map<int,std::set<int> > dmap;
  AdjGraph(const std::pair<int,int>& range)
    : range(range) { }
  void add(int i, int j) {
    if (i <  this->range.first)  return;
    if (i >= this->range.second) return;
    this->dmap[i].insert(j);
  }
  void set_ngbrs(int i, GSet& gset){
    gset.swap(this->dmap[i]);
  }
  void scatter() {}
  void clear() { this->dmap.clear(); }
};

struct ProfileAdj : public Profile {

  AdjGraph pgraph;

  ProfileAdj(const Dofset &ds)
    : Profile(ds), pgraph(ds.getRange())
  { 
    //this->lgraph = &(this->pgraph); 
  }

  void build() 
  {
    std::pair<int,int> range = this->dofset->getRange();
    std::vector<int>&  xadj   = this->buf1;
    std::vector<int>&  adjncy = this->buf2;
    int first = range.first;
    int last  = range.second;
    int dof = first;
    int neq = last - first;
    GSet adj;
    xadj.reserve(neq+1);
    xadj.push_back(0);
    for (int k=0; k<neq; k++) {
      adj.clear();
      lgraph->set_ngbrs(dof++, adj);
      adjncy.insert(adjncy.end(), adj.begin(), adj.end());
      xadj.push_back(adjncy.size());
    }
  }
};

struct ProfileNNZ : public Profile {
  ProfileNNZ(const Dofset &ds)
    : Profile(ds)
  {
    
  }
  void build() 
  {
    std::pair<int,int> range = this->dofset->getRange();
    std::vector<int>& d_nnz  = this->buf1;
    std::vector<int>& o_nnz  = this->buf2;
    int first = range.first;
    int last  = range.second;
    int dof = first;
    int neq = last - first;
    GSet adj;
    d_nnz.resize(neq, 0); 
    o_nnz.resize(neq, 0);
    for (int k=0; k<neq; k++) {
      adj.clear();
      lgraph->set_ngbrs(dof++,adj);
      GSet::iterator n  = adj.begin();
      GSet::iterator ne = adj.end();
      while (n!=ne) {
	int node = *n++;
	if (node < first || node >= last)
	  o_nnz[k]++;
	else 
	  d_nnz[k]++;
      }
    }
  }
};

struct AdjSdGraph : public StoreGraph {
  std::map<int,std::set<int> > sdmap;

  void add(int i, int j) {
    this->sdmap[i].insert(j);
  }
  void set_ngbrs(int i, GSet& gset){
    gset.swap(this->sdmap[i]);
  }
  void scatter() {}
  void clear() { this->sdmap.clear(); }
};

struct ProfileAdjSd : public Profile {

  AdjSdGraph pgraph;

  ProfileAdjSd(const Dofset &ds)
    : Profile(ds), pgraph()
  { 
    this->lgraph = &(this->pgraph); 
  }
  void build()
  {
    std::vector<int>& dofs   = this->buf1;
    std::vector<int>& xadj   = this->buf2;
    std::vector<int>& adjncy = this->buf3;

    this->dofset->getLocalDofs(dofs);
    xadj.reserve(dofs.size()+1);
    xadj.push_back(0);

    int neq = dofs.size();
    if (neq == 0) return;

    std::pair<int,int> range = this->dofset->getRange();
    std::vector<int> ghosts;
    this->dofset->getGhostDofs(ghosts);
    std::vector<int>::iterator gb = ghosts.begin();
    std::vector<int>::iterator ge = ghosts.end();

    GSet adj;
    for (int k=0; k<neq; k++) {
      adj.clear();
      lgraph->set_ngbrs(dofs[k], adj);
#if 1
      GSet::iterator d = adj.begin();
      GSet::iterator e = adj.end();
      while (d != e) {
	GSet::iterator p = d++;
	int node = *p;
	if (node < range.first || node >= range.second)
	  if (!std::binary_search(gb, ge, node)) adj.erase(p);
      }
#endif
      adjncy.insert(adjncy.end(), adj.begin(), adj.end());
      xadj.push_back(adjncy.size());
    }

  }
};

PF4PY_NAMESPACE_END

// ---------------------------------------------------------------- //


class WallData;

PF4PY_NAMESPACE_BEGIN

struct ArgsNS : public AppCtx::Args {
  Time           tstar;
  State          state1, state0;
  vector<double> hmin;
  WallData*      wall_data;
  GlobParam      glob_param;
  Vec            res, dummy_res;
  Mat            Jac, dummy_Jac;
  std::string    jobinfo;

  static const char COMP_PROFILE[];
  static const char COMP_RES[];
  static const char COMP_RES_JAC[];

  ArgsNS() : Args(),
	     tstar(), state1(), state0(),
	     hmin(1), wall_data(NULL), glob_param(),
	     res(PETSC_NULL), dummy_res(PETSC_NULL),
	     Jac(PETSC_NULL), dummy_Jac(PETSC_NULL),
	     jobinfo()
  { }
  ~ArgsNS() 
  {
    //PF4PY_DELETE_SCLR(this->wall_data); XXX
    PF4PY_PETSC_DESTROY(VecDestroy, this->dummy_res);
    PF4PY_PETSC_DESTROY(MatDestroy, this->dummy_Jac);
  }
  
  const char* job()  const { return this->jobinfo.c_str(); }
  const Time* time() const { return &this->tstar;  }
  
  void pack(const std::string& jobname,
	    double t1, Vec x1,
	    double t0, Vec x0,
	    Vec r, Mat J,
	    double alpha, bool steady)
  {
    Args::Impl* argl = this->argl();

    double Dt = t1 - t0;
    if (Dt == 0.0) {
      steady = true;
      alpha = 1.0;
      Dt = 1.0;
    }

    // clear arg list
    argl->clear();

    // tstart:  t^* = t^0 + alpha * (t^1 - t^0)
    this->tstar.set(t0 + alpha * (t1-t0));

    // add new state
    this->state1.set(t1, x1);
    argl->arg_add(&this->state1, IN_VECTOR | USE_TIME_DATA);

    // add old state
    this->state0.set(t0, (x0!=PETSC_NULL) ? x0 : x1);
    argl->arg_add(&this->state0, IN_VECTOR | USE_TIME_DATA);

    // add residual vector and jacobian matrix
    this->res = r;
    this->Jac = J;
    if (this->res == PETSC_NULL) {
      if (this->dummy_res == PETSC_NULL) 
	PF4PY_PETSC_CALL(VecCreateEmpty, (PETSC_COMM_SELF, &this->dummy_res));
      argl->arg_add(&this->dummy_res, OUT_VECTOR);
    } else {
      argl->arg_add(&this->res, OUT_VECTOR);
    }
    if (this->Jac == PETSC_NULL) {
      this->jobinfo = ArgsNS::COMP_RES;
    } else {
      argl->arg_add(&this->Jac, OUT_MATRIX);
      this->jobinfo = ArgsNS::COMP_RES_JAC;
    }

    if (!jobname.empty()) {
      this->jobinfo = std::string("comp_") + jobname;
    }

    // add hmin (XXX: What is this arg for???)
    argl->arg_add(&this->hmin, VECTOR_MIN);
    
    // setup and add global parameters
    this->glob_param.steady    = steady ? 1 : 0;
    this->glob_param.alpha     = alpha;
    this->glob_param.Dt        = Dt;
    this->glob_param.inwt      = 0; // XXX: this is ok?
    this->glob_param.time      = &this->tstar;
    this->glob_param.state     = &this->state1;
    this->glob_param.state_old = &this->state0;
    this->glob_param.x         = this->state1._vs;
    this->glob_param.xold      = this->state0._vs;
    argl->arg_add(&this->glob_param, USER_DATA);
  
    // add wall_data
    argl->arg_add(this->wall_data, USER_DATA);
  }
  
};

const char ArgsNS::COMP_PROFILE[] = "comp_mat";
const char ArgsNS::COMP_RES[]     = "comp_res";
const char ArgsNS::COMP_RES_JAC[] = "comp_mat_res";

PF4PY_NAMESPACE_END



PF4PY_NAMESPACE_BEGIN

AppNS::~AppNS()
{ }

AppNS::AppNS(const AppNS& ns)
  : AppCtx(ns)
{
  this->args.reset(new ArgsNS);
}

AppNS::AppNS()
  :  AppCtx()
{ 
  this->args.reset(new ArgsNS);
}

void
AppNS::assemble(const Domain* domain,
		const std::string& jobname,
		double t, Vec x,
		Vec r, Mat J) const
{
  ArgsNS* args = reinterpret_cast<ArgsNS*>(this->args.get());
  args->pack(jobname, t, x, t, PETSC_NULL, r, J, 1.0, true);

  AppCtx::assemble(domain, args);

  if (r != PETSC_NULL) PF4PY_PETSC_CALL(VecScale, (r, -1));
}

void
AppNS::assemble(const Domain* domain,
		const std::string& jobname,
		double t1, Vec x1,
		double t0, Vec x0,
		Vec r, Mat J, double alpha) const
{
  if (alpha == 0.0) 
    throw Error("AppNX: invalid value, 'alpha' == 0.0");

  ArgsNS* args = reinterpret_cast<ArgsNS*>(this->args.get());
  args->pack(jobname, t1, x1, t0, x0, r, J, alpha, false);

  AppCtx::assemble(domain, args);
  
  if (r != PETSC_NULL) PF4PY_PETSC_CALL(VecScale, (r, -1.0/alpha));
}

PF4PY_NAMESPACE_END


PF4PY_NAMESPACE_BEGIN

bool
AppNS::sdgraph(const Domain* domain,
	       std::vector<int>& dofs,
	       std::vector<int>& xadj,
	       std::vector<int>& adjncy) const
{
  const Dofset& dofset = domain->getDofset();
  Dofset::Impl* dofmap = dofset;
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  
  ProfileAdjSd A(dofset);

  ArgsNS* args = reinterpret_cast<ArgsNS*>(this->args.get());
  args->tstar = 0.0;
  args->jobinfo = ArgsNS::COMP_PROFILE;
  args->argl()->clear();
  args->argl()->arg_add(&A,PROFILE|PFMAT);

  AppCtx::assemble(domain, args);
  
  dofs.swap(A.buf1); xadj.swap(A.buf2); adjncy.swap(A.buf3);
  
  return true;
}

bool
AppNS::profile(const Domain* domain,
	       std::vector<int>& xadj,
	       std::vector<int>& adjncy) const
{
  const Dofset& dofset = domain->getDofset();
  Dofset::Impl* dofmap = dofset;
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  
  ProfileAdj A(dofset);

  ArgsNS* args = reinterpret_cast<ArgsNS*>(this->args.get());
  args->tstar = 0.0;
  args->jobinfo = ArgsNS::COMP_PROFILE;
  args->argl()->clear();
  args->argl()->arg_add(&A,PROFILE|PFMAT);

  AppCtx::assemble(domain, args);
  
  xadj.swap(A.buf1); adjncy.swap(A.buf2);
  
  return true;
}

PF4PY_NAMESPACE_END


// ---------------------------------------------------------------- //

extern int local_time_step_g;
extern int consistent_supg_matrix_g;
extern int comp_mat_each_time_step_g;

PF4PY_NAMESPACE_BEGIN

struct ArgsAD : public AppCtx::Args {
  Time           tstar;
  State          state1, state0;
  vector<double> hmin;
  WallData*      wall_data;
  GlobParam      glob_param;
  Vec            res, dummy_res;
  Mat            Jac, dummy_Jac;
  std::string    jobinfo;

  static const char COMP_PROFILE[];
  static const char COMP_RES[];
  static const char COMP_RES_JAC[];

  ArgsAD() : Args(),
	     tstar(), state1(), state0(),
	     hmin(1), wall_data(NULL), glob_param(),
	     res(PETSC_NULL), dummy_res(PETSC_NULL),
	     Jac(PETSC_NULL), dummy_Jac(PETSC_NULL),
	     jobinfo()
  { }
  ~ArgsAD() 
  {
    //PF4PY_DELETE_SCLR(this->wall_data); XXX
    PF4PY_PETSC_DESTROY(VecDestroy, this->dummy_res);
    PF4PY_PETSC_DESTROY(MatDestroy, this->dummy_Jac);
  }
  
  const char* job()  const { return this->jobinfo.c_str(); }
  const Time* time() const { return &this->tstar;  }
  
  void pack(const std::string& jobname,
	    double t1, Vec x1,
	    double t0, Vec x0,
	    Vec r, Mat J,
	    double alpha, bool steady)
  {
    Args::Impl* argl = this->argl();

    double Dt = t1 - t0;
    if (Dt == 0.0) {
      steady = true;
      alpha = 1.0;
      Dt = 1.0;
    }
    
    local_time_step_g = 0;

    // clear arg list
    argl->clear();

    // tstart:  t^* = t^0 + alpha * (t^1 - t^0)
    this->tstar.set(t0 + alpha * (t1-t0));

    // add old state
    this->state0.set(t0, (x0!=PETSC_NULL) ? x0 : x1);
    argl->arg_add(&this->state0, IN_VECTOR | USE_TIME_DATA);

    // add new state
    this->state1.set(t1, x1);
    argl->arg_add(&this->state1, IN_VECTOR | USE_TIME_DATA);

    // add residual vector
    this->res = r;
    if (this->res == PETSC_NULL) {
      if (this->dummy_res == PETSC_NULL)
	PF4PY_PETSC_CALL(VecCreateEmpty, (PETSC_COMM_SELF, &this->dummy_res));
      argl->arg_add(&this->dummy_res, OUT_VECTOR);
    } else {
      argl->arg_add(&this->res, OUT_VECTOR);
    }

    // add hmin
    argl->arg_add(&this->hmin, VECTOR_MIN);

    // add  jacobian matrix
    this->Jac = J;
    if (this->Jac == PETSC_NULL) {
      if (this->dummy_Jac == PETSC_NULL)
	PF4PY_PETSC_CALL(MatCreateEmpty, (PETSC_COMM_SELF, &this->dummy_Jac));
      argl->arg_add(&this->dummy_Jac, OUT_MATRIX);
      this->jobinfo = ArgsAD::COMP_RES;
      
      consistent_supg_matrix_g  = 0;
      comp_mat_each_time_step_g = 0;

    } else {
      argl->arg_add(&this->Jac, OUT_MATRIX);
      this->jobinfo = ArgsAD::COMP_RES_JAC;

      consistent_supg_matrix_g  = 1;
      comp_mat_each_time_step_g = 1;

    }

    // setup and add global parameters
    this->glob_param.steady    = steady?1:0;
    this->glob_param.alpha     = alpha;
    this->glob_param.Dt        = Dt;
    this->glob_param.inwt      = 0;
    this->glob_param.time      = &this->tstar;
    this->glob_param.state     = &this->state1;
    this->glob_param.state_old = &this->state0;
    this->glob_param.x         = this->state1._vs;
    this->glob_param.xold      = this->state0._vs;
    argl->arg_add(&this->glob_param, USER_DATA);
  
    // add wall_data
    argl->arg_add(this->wall_data, USER_DATA);
  }
  
};

const char ArgsAD::COMP_PROFILE[] = "comp_prof";
const char ArgsAD::COMP_RES[]     = "comp_res";
const char ArgsAD::COMP_RES_JAC[] = "comp_res";

PF4PY_NAMESPACE_END


PF4PY_NAMESPACE_BEGIN

AppAD::~AppAD()
{ }

AppAD::AppAD(const AppAD& ad)
  : AppCtx(ad)
{
  this->args.reset(new ArgsAD);
}

AppAD::AppAD()
  : AppCtx()
{ 
  this->args.reset(new ArgsAD);
}

void
AppAD::assemble(const Domain* domain,
		const std::string& jobname,
		double t, Vec x,
		Vec r, Mat J) const
{
  ArgsAD* args = reinterpret_cast<ArgsAD*>(this->args.get());
  args->pack(jobname, t, x, t, PETSC_NULL, r, J, 1.0, true);

  AppCtx::assemble(domain, args);

  if (r != PETSC_NULL) PF4PY_PETSC_CALL(VecScale, (r, -1));
}

void
AppAD::assemble(const Domain* domain,
		const std::string& jobname,
		double t1, Vec x1,
		double t0, Vec x0,
		Vec r, Mat J, double alpha) const
{
  ArgsAD* args = reinterpret_cast<ArgsAD*>(this->args.get());
  args->pack(jobname, t1, x1, t0, x0, r, J, alpha, false);

  AppCtx::assemble(domain, args);
  
  PetscScalar rscale = (alpha == 0.0) ? -1.0 : -1.0/alpha;
  if (r != PETSC_NULL) PF4PY_PETSC_CALL(VecScale, (r, rscale));
}

PF4PY_NAMESPACE_END


PF4PY_NAMESPACE_BEGIN

bool
AppAD::sdgraph(const Domain* domain,
	       std::vector<int>& ldofs,
	       std::vector<int>& xadj,
	       std::vector<int>& adjncy) const
{
  const Dofset& dofset = domain->getDofset();
  Dofset::Impl* dofmap = dofset;
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  
  ProfileAdjSd A(dofset);

  ArgsAD* args = reinterpret_cast<ArgsAD*>(this->args.get());
  args->tstar = 0.0;
  args->jobinfo = ArgsAD::COMP_PROFILE;
  args->argl()->clear();
  args->argl()->arg_add(&A,PROFILE|PFMAT);

  AppCtx::assemble(domain, args);
  
  ldofs.swap(A.buf1); xadj.swap(A.buf2); adjncy.swap(A.buf3);
  
  return true;
}

bool
AppAD::profile(const Domain* domain,
	       std::vector<int>& xadj,
	       std::vector<int>& adjncy) const
{
  const Dofset& dofset = domain->getDofset();
  Dofset::Impl* dofmap = dofset;
  if (dofmap == NULL) throw Error("Domain: dofset not ready");
  
  ProfileAdj A(dofset);

  ArgsAD* args = reinterpret_cast<ArgsAD*>(this->args.get());
  args->tstar = 0.0;
  args->jobinfo = ArgsAD::COMP_PROFILE;
  args->argl()->clear();
  args->argl()->arg_add(&A,PROFILE|PFMAT);

  AppCtx::assemble(domain, args);
  
  xadj.swap(A.buf1); adjncy.swap(A.buf2);
  
  return true;
}

PF4PY_NAMESPACE_END

// ---------------------------------------------------------------- //
