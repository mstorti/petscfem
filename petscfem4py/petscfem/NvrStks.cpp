// $Id: NvrStks.cpp,v 1.1.2.3 2006/05/30 20:22:23 dalcinl Exp $

#include "NvrStks.h"

#include <fem.h>
#include <sttfilter.h>

#include <applications/ns/nsi_tet.h>

PYPF_NAMESPACE_BEGIN

struct State: public ::State {
  Vec _vs;

  State()                 : ::State(), _vs(PETSC_NULL)
  { this->vec = &this->_vs; }
  State(const State& stt) : ::State(), _vs(PETSC_NULL)
  { this->set(stt._vs, stt.time); }
  State(Vec v, double t)  : ::State(), _vs(PETSC_NULL)
  { this->set(v,t); }
  ~State()  // ::State::~State() calls VecDestroy()
  { if(this->vec == &this->_vs) this->vec = NULL; }
  void set(Vec v, double t)
  { this->_vs = v; this->vec = &this->_vs; this->time = t; }
};

static const char COMP_RES[]     = "comp_res";
static const char COMP_RES_JAC[] = "comp_mat_res";

struct ArgsNS : public ArgList {
  State          state0, state1;
  Time           tstar;
  vector<double> hmin;
  WallData*      wall_data;
  GlobParam      glob_param;
  Vec            xst;
  Vec            res;
  Mat            Jac;
  const char*    jobinfo;
  ArgsNS() : ArgList(),
	     state0(), state1(), tstar(),
	     hmin(1), wall_data(NULL), glob_param(),
	     xst(PETSC_NULL), res(PETSC_NULL), Jac(PETSC_NULL),
	     jobinfo(NULL)
             { }
  ~ArgsNS()  
  {
    PYPF_PETSC_DESTROY(VecDestroy, this->xst);
    PYPF_DELETE_SCLR(this->wall_data);
  }

  const char* job()  const { return this->jobinfo; }

  const Time* time() const { return &this->tstar;  }
  
  void pack(Vec x0, double t0,
	    Vec x1, double t1,
	    Vec r, Mat J,
	    double alpha, bool steady)
  {

    ArgList::Base& argl = *this;

    if (x0 == PETSC_NULL) {
      if (!this->xst) PYPF_PETSC_CALL(VecDuplicate(x1, &this->xst));
      x0 = this->xst;
      PYPF_PETSC_CALL(VecCopy(x1, x0));
    } else if (x1 == PETSC_NULL) {
      if (!this->xst) PYPF_PETSC_CALL(VecDuplicate(x0, &this->xst));
      x1 = this->xst;
      PYPF_PETSC_CALL(VecCopy(x0, x1));
    }

    // clear arg list
    argl.clear();

    // tstart:  t^* = t^0 + alpha * (t^1 - t^0)
    this->tstar.set(t0 + alpha * (t1-t0));


    // add new state
    this->state1.set(x1, t1);
    argl.arg_add(&this->state1, IN_VECTOR | USE_TIME_DATA);

    // add old state
    this->state0.set(x0, t0);
    argl.arg_add(&this->state0, IN_VECTOR | USE_TIME_DATA);

    // add residual vector and jacobian matrix
    this->res = r;
    this->Jac = J;
    argl.arg_add(&this->res, OUT_VECTOR);
    if (this->Jac == PETSC_NULL) {
      this->jobinfo = COMP_RES;
    } else {
      argl.arg_add(&this->Jac, OUT_MATRIX);
      this->jobinfo = COMP_RES_JAC;
    }

    // add hmin (XXX: What is this arg for???)
    argl.arg_add(&this->hmin, VECTOR_MIN);

    // setup and add global parameters
    this->glob_param.steady    = steady?1:0;
    this->glob_param.alpha     = alpha;
    this->glob_param.Dt        = t1-t0;;
    this->glob_param.inwt      = 0; // XXX: this is ok?
    this->glob_param.state_old = &this->state0;
    this->glob_param.state     = &this->state1;
    argl.arg_add(&this->glob_param, USER_DATA);
  
    // add wall_data
    argl.arg_add(this->wall_data, USER_DATA);
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


void 
NvrStks::assemble(Vec x, double t, Vec r, Mat J) const
{
  // XXX test vector and matrix sizes against dofmap !!!

  double alpha  = 1.0;
  bool   steady = true;

  args->pack(PETSC_NULL, t, x, t, r, J, alpha, steady);
  Application::assemble(*this, *this->args);
  PYPF_PETSC_CALL(VecScale(r, -1.0/alpha));
}


void
NvrStks::assemble(Vec x0, double t0, Vec x1, double t1, Vec r, Mat J) const
{
  // XXX test vector and matrix sizes against dofmap !!!

  double alpha  = this->alpha;
  bool   steady = this->steady;
  if (steady) alpha = 1.0;
  else if (alpha <= 0.0) throw Error("invalid value, 'alpha' <= 0.0");
  else if (alpha >  1.0) throw Error("invalid value, 'alpha' > 1.0");
  if (t1 < t0 )          throw Error("invalid value, 't1' < 't0'");
    
  args->pack(x0, t0, x1, t1, r, J, alpha, steady);
  Application::assemble(*this, *this->args);
  PYPF_PETSC_CALL(VecScale(r, -1.0/alpha));
}

PYPF_NAMESPACE_END
