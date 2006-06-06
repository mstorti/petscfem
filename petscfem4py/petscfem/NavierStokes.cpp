// $Id: NavierStokes.cpp,v 1.1.2.7 2006/06/06 15:46:52 dalcinl Exp $

#include "NavierStokes.h"

#include <fem.h>
#include <readmesh.h>
#include <arglist.h>
#include <sttfilter.h>

#include <applications/ns/nsi_tet.h>


PYPF_NAMESPACE_BEGIN

static const char COMP_RES[]     = "comp_res";
static const char COMP_RES_JAC[] = "comp_mat_res";

class ArgList {

public:
  arg_list argl;
  
public:  
  ArgList() : argl() { }
  virtual ~ArgList()  { this->argl.clear(); }

  virtual void pack(Vec x0, double t0, 
		    Vec x1, double t1,
		    Vec r, Mat J,
		    double alpha, bool steady) = 0;

  virtual void pack(Vec x0, double t0,
		    Vec x1, double t1,
		    Vec r, Mat J,
		    double alpha=1.0)
  { this->pack(x0, t0, x1, t1, r, J, alpha, false); }

  virtual void pack(Vec x, double t, Vec r, Mat J) 
  { this->pack(PETSC_NULL, t, x, t, r, J, 1.0, true);}

  virtual void pack(Vec x, Vec r, Mat J) 
  { this->pack(x, 0.0, r, J); }

  virtual void pack(Vec r, Mat J)
  { this->pack(PETSC_NULL, r, J); }

  virtual void clear()
  { this->argl.clear(); }

};


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


class ArgsNS {

public:
  arg_list       argl;
  State          state0, state1;
  Time           tstar;
  vector<double> hmin;
  WallData*      wall_data;
  GlobParam      glob_param;
  Vec            stt;
  Vec            res;
  Mat            Jac;
  const char*    jobinfo;

public:  
  ArgsNS() : argl(),
	     state0(), state1(),
	     tstar(),
	     hmin(1), wall_data(NULL), glob_param(),
	     stt(PETSC_NULL), res(PETSC_NULL), Jac(PETSC_NULL),
	     jobinfo(NULL)
             { }
  ~ArgsNS() 
  {
    PYPF_PETSC_DESTROY(VecDestroy,this->stt);
    PYPF_DELETE_SCLR(this->wall_data); 
  }
  
  void pack(Vec x0, double t0, 
	    Vec x1, double t1,
	    Vec r, Mat J,
	    double alpha, bool steady)
  {
    if (x0 == PETSC_NULL) {
      if (!this->stt) VecDuplicate(x1, &this->stt);
      x0 = this->stt;
      VecCopy(x1, x0);
    } else if (x1 == PETSC_NULL) {
      if (!this->stt) VecDuplicate(x0, &this->stt);
      x1 = this->stt;
      VecCopy(x0, x1);
    }

    // tstart:  t^* = t^0 + alpha * (t^1 - t^0)
    this->tstar.set(t0 + alpha * (t1-t0));

    // clear arglist
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

  void pack(Vec x, double t, Vec r, Mat J)
  { this->pack(PETSC_NULL, t, x, t, r, J, 1.0, true);}

  void clear() 
  { this->argl.clear(); }

};

PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN

NavierStokes::~NavierStokes()
{ 
  PYPF_DELETE_SCLR(this->args);
}

NavierStokes::NavierStokes(Mesh& mesh, Dofset& dofset)
  : Problem(mesh, dofset),
    args(new ArgsNS)
{ }

NavierStokes::NavierStokes(Nodeset& nodeset,
			   const std::vector<Elemset*>& elemsets,
			   Dofset& dofset)
  : Problem(nodeset, elemsets, dofset),
    args(new ArgsNS)
{ }

void
NavierStokes::assemble(Vec x, double t, Vec r, Mat J)
{
  this->assemble(PETSC_NULL, t, x, t, r, J);
}

static inline void
check_time_parms(double t0, double t1, double alpha)
{
}

void
NavierStokes::assemble(Vec x0, double t0,
		       Vec x1, double t1,
		       Vec r, Mat J,
		       double alpha, bool steady)
{
  // XXX test vector and matrix sizes against dofmap !!!

  if (steady) alpha = 1.0;
  else if (alpha <= 0.0) throw Error("invalid value, 'alpha' <= 0.0");
  else if (alpha >  1.0) throw Error("invalid value, 'alpha' > 1.0");
  if (t1 < t0 )          throw Error("invalid value, 't1' < 't0'");
    
  const Mesh&   mesh   = this->getMesh();
  const DofMap& dofmap = this->getDofMap();
  ArgsNS*       args   = this->args;

  args->pack(x0, t0, x1, t1, r, J, alpha, steady);
  this->preAssemble();
  int ierr = ::assemble(mesh, args->argl, dofmap, args->jobinfo, &args->tstar);
  this->postAssemble();
  args->clear();
  if (ierr) throw Error("NavierStokes::assemble(...)");
  VecScale(r, -1.0/alpha);
}


// static PetscErrorCode
// matzero_unary(Mat mat) { return 0;}

// static PetscErrorCode
// matzero_setvalues
// (Mat mat,
//  PetscInt m,const PetscInt im[],
//  PetscInt n,const PetscInt in[],
//  const PetscScalar v[], InsertMode addv) { return 0; }

// static PetscErrorCode
// matzero_assembly(Mat mat, MatAssemblyType type) { return 0;}


// class MatZero {
//   Mat mat;
// public:
//   ~MatZero() 
//   { 
//     if (this->mat) MatDestroy(this->mat);
//     this->mat= PETSC_NULL; 
//   }
//   MatZero(MPI_Comm comm, PetscInt m, PetscInt n, PetscInt M, PetscInt N)
//     : mat(PETSC_NULL)
//   {
//     //Mat& mat = this->mat;
//     MatCreate(comm, &mat);
//     MatSetSizes(mat, m, n, M, N);
//     MatSetType(mat, MATSHELL);
//     typedef void (*MatOp)(void);
//     MatShellSetOperation(mat, MATOP_SET_VALUES,     (MatOp)matzero_setvalues);
//     MatShellSetOperation(mat, MATOP_ASSEMBLY_BEGIN, (MatOp)matzero_assembly);
//     MatShellSetOperation(mat, MATOP_ASSEMBLY_END,   (MatOp)matzero_assembly);
//     MatShellSetOperation(mat, MATOP_ZERO_ENTRIES,   (MatOp)matzero_unary);
//   }
//   inline operator Mat() { return this->mat; }
// };


// void
// NavierStokes::assembleResidual(Vec x0, double t0,
// 			       Vec x1, double t1,
// 			       Vec r,   
// 			       double alpha, bool steady)
// {
//   int local, global;
//   this->getDofSizes(&local, &global);
//   MatZero J(this->getComm(), local, local, global, global);
//   this->assemble(x0, t0, x1, t1, r, J, alpha, steady);
// }

PYPF_NAMESPACE_END
