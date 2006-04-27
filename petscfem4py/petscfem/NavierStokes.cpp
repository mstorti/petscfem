// $Id: NavierStokes.cpp,v 1.1.2.6 2006/04/27 19:09:17 rodrigop Exp $

#include "NavierStokes.h"

#include <fem.h>
#include <readmesh.h>
#include <arglist.h>
#include <sttfilter.h>

#include <applications/ns/nsi_tet.h>


//extern Mesh*          GLOBAL_MESH;

PYPF_NAMESPACE_BEGIN


static Vec  _v_dummy;
static Time _t_dummy;

struct NSArgs {

  /* XXX: ~::State() calls VecDestroy()!!! */
  struct MyState: public State {
    MyState() : State(_v_dummy, _t_dummy) { this->vec = NULL; }
    ~MyState()                            { this->vec = NULL; }
    void set_vec(Vec &v)                  { this->vec = &v;   }
  };

  arg_list argl;
  Time time0, time1, time_star;
  Vec vecdummy;
  MyState state0, state1;
  vector<double> hmin;
  WallData *wall_data;
  GlobParam glob_param;

  NSArgs() : argl(),
	     time0(), time1(), time_star(),
	     vecdummy(PETSC_NULL), 
	     state0(), state1(),
	     hmin(1),
	     wall_data(NULL),
	     glob_param() { }
  ~NSArgs() { 
    this->argl.clear();
    PYPF_DELETE_SCLR(this->wall_data);
  };
  
  void pack(Vec& x0, double t0,  
	    Vec& x1, double t1,
	    Vec& r, Mat& J,
	    double alpha, int steady)
  {
    // clear arglist
    this->argl.clear();
    
    // new state
    this->time1.set(t1);
    this->state1.set_vec(x1);
    this->state1.set_time(this->time1);
    argl.arg_add(&state1, IN_VECTOR | USE_TIME_DATA);
    // old state
    this->time0.set(t0);
    this->state0.set_vec(x0);
    this->state0.set_time(this->time0);
    argl.arg_add(&state0, IN_VECTOR | USE_TIME_DATA);
    // time_start:  t^* = t^0 + alpha * Dt, Dt = t^1 - t^0
    this->time_star.set(t0 + alpha * (t1-t0));
    
    // add residual vector and jacobian matrix
    argl.arg_add(&r, OUT_VECTOR);
    argl.arg_add(&J, OUT_MATRIX);

    // add hmin (XXX: What is that???)
    argl.arg_add(&(this->hmin), VECTOR_MIN);

    // setup and add global_param
    this->glob_param.steady    = steady;
    this->glob_param.alpha     = alpha;
    this->glob_param.Dt        = t1-t0;;
    this->glob_param.inwt      = 0; // XXX: this is ok?
    this->glob_param.state_old = &(this->state0);
    this->glob_param.state     = &(this->state1);
    argl.arg_add(&glob_param, USER_DATA);
  
    // add wall_data
    argl.arg_add(this->wall_data, USER_DATA);
  }

  void pack(Vec& x0, double t0,  
	    Vec& x1, double t1,
	    Vec& r,  
	    double alpha, int steady)
  {
    // clear arglist
    this->argl.clear();
    
    // new state
    this->time1.set(t1);
    this->state1.set_vec(x1);
    this->state1.set_time(this->time1);
    argl.arg_add(&state1, IN_VECTOR | USE_TIME_DATA);
    // old state
    this->time0.set(t0);
    this->state0.set_vec(x0);
    this->state0.set_time(this->time0);
    argl.arg_add(&state0, IN_VECTOR | USE_TIME_DATA);
    // time_start:  t^* = t^0 + alpha * Dt, Dt = t^1 - t^0
    this->time_star.set(t0 + alpha * (t1-t0));
    
    // add residual vector
    argl.arg_add(&r, OUT_VECTOR);

    // add hmin (XXX: What is that???)
    argl.arg_add(&(this->hmin), VECTOR_MIN);

    // setup and add global_param
    this->glob_param.steady    = steady;
    this->glob_param.alpha     = alpha;
    this->glob_param.Dt        = t1-t0;;
    this->glob_param.inwt      = 0; // XXX: this is ok?
    this->glob_param.state_old = &(this->state0);
    this->glob_param.state     = &(this->state1);
    argl.arg_add(&glob_param, USER_DATA);
  
    // add wall_data
    argl.arg_add(this->wall_data, USER_DATA);
  }


  void clear() { this->argl.clear(); }

};

NavierStokes::~NavierStokes() 
{ 
  NSArgs*& nsargs = reinterpret_cast<NSArgs*&>(this->nsargs);
  PYPF_DELETE_SCLR(nsargs);
}

NavierStokes::NavierStokes(Mesh& mesh, DofMap& dofmap)
  : Problem(mesh, dofmap),
    nsargs(reinterpret_cast<void*>(new NSArgs))
{ }

NavierStokes::NavierStokes(Nodeset& nodeset,
			   const std::vector<Elemset*>& elemsets,
			   Dofset& dofset)
  : Problem(nodeset, elemsets, dofset),
    nsargs(reinterpret_cast<void*>(new NSArgs))
{ }

void
NavierStokes::assemble(Vec x, double t, Vec r, Mat J)
{
  // XXX test vector and matrix sizes against dofmap !!!

  Mesh::Base*   mesh   = *this->mesh;
  DofMap::Base* dofmap = *this->dofmap;
  NSArgs* args = reinterpret_cast<NSArgs*>(this->nsargs);

  //GLOBAL_MESH = mesh;
  
  args->pack(x, t, x, t, r, J, 1.0, 1);
  this->preAssemble();
  int ierr = ::assemble(mesh, args->argl, dofmap, "comp_mat_res", &args->time_star);
  this->postAssemble();
  args->clear();
  if (ierr) throw Error("NavierStokes::assemble()");
  VecScale(r, -1.0);
}

static inline void
check_time_parms(double t0, double t1, double alpha)
{
  if (t1 < t0 )    throw Error("invalid value, 't1' < 't0'");
  if (alpha <= 0.0) throw Error("invalid value, 'alpha' <= 0.0");
  if (alpha >  1.0) throw Error("invalid value, 'alpha' > 1.0");
}

void
NavierStokes::assemble(Vec x0, double t0, 
		       Vec x1, double t1,
		       Vec r, Mat J, 
		       double alpha)
{
  // XXX test vector and matrix sizes against dofmap !!!

  /* test */ check_time_parms(t0, t1, alpha);
    
  Mesh::Base*   mesh   = *this->mesh;
  DofMap::Base* dofmap = *this->dofmap;
  NSArgs* args = reinterpret_cast<NSArgs*>(this->nsargs);

  args->pack(x0, t0, x1, t1, r, J, alpha, 0);
  this->preAssemble();
  int ierr = ::assemble(mesh, args->argl, dofmap, "comp_mat_res", &args->time_star);
  this->postAssemble();
  args->clear();
  if (ierr) throw Error("NavierStokes::assemble()");
  VecScale(r, -1.0/alpha);
}


static PetscErrorCode
matzero_unary(Mat mat) { return 0;}

static PetscErrorCode
matzero_setvalues
(Mat mat,
 PetscInt m,const PetscInt im[],
 PetscInt n,const PetscInt in[],
 const PetscScalar v[], InsertMode addv) { return 0; }

static PetscErrorCode
matzero_assembly(Mat mat, MatAssemblyType type) { return 0;}


class MatZero {
  Mat mat;
public:
  ~MatZero() 
  { 
    if (this->mat) MatDestroy(this->mat); 
    this->mat= PETSC_NULL; 
  }
  MatZero(MPI_Comm comm, PetscInt m, PetscInt n, PetscInt M, PetscInt N)
    : mat(PETSC_NULL)
  {
    //Mat& mat = this->mat;
    MatCreate(comm, &mat);
    MatSetSizes(mat, m, n, M, N);
    MatSetType(mat, MATSHELL);
    typedef void (*MatOp)(void);
    MatShellSetOperation(mat, MATOP_SET_VALUES,     (MatOp)matzero_setvalues);
    MatShellSetOperation(mat, MATOP_ASSEMBLY_BEGIN, (MatOp)matzero_assembly);
    MatShellSetOperation(mat, MATOP_ASSEMBLY_END,   (MatOp)matzero_assembly);
    MatShellSetOperation(mat, MATOP_ZERO_ENTRIES,   (MatOp)matzero_unary);
  }
  inline operator Mat() { return this->mat; }
};


void
NavierStokes::assembleResidual(Vec x0, double t0,
			       Vec x1, double t1,
			       Vec r,   
			       double alpha)
{
  int local, global;
  this->getDofSizes(&local, &global);
  MatZero J(this->getComm(), local, local, global, global);
  this->preAssemble();
  this->assemble(x0, t0, x1, t1, r, J, alpha);
  this->postAssemble();
  VecScale(r, -1.0/alpha);
}

PYPF_NAMESPACE_END
