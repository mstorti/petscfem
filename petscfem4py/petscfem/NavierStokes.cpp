// $Id: NavierStokes.cpp,v 1.1.2.5 2006/03/30 15:18:14 rodrigop Exp $

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
    if (this->wall_data) delete this->wall_data;
  };
  
  void pack(Vec& x0, double t0,  Vec& x1, double t1,
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

  void clear() { this->argl.clear(); }

};

NavierStokes::~NavierStokes() 
{ 
  NSArgs*& nsargs = reinterpret_cast<NSArgs*&>(this->nsargs);
  PYPF_DELETE_SCLR(nsargs);
}

// NavierStokes::NavierStokes()
//   : Problem(),
//     nsargs(reinterpret_cast<void*>(new NSArgs))
// { }

// NavierStokes::NavierStokes(const NavierStokes& ns)
//   : Problem(ns),
//     nsargs(reinterpret_cast<void*>(new NSArgs))
// { }

NavierStokes::NavierStokes(Mesh* mesh, DofMap* dofmap) 
  : Problem(mesh, dofmap),
    nsargs(reinterpret_cast<void*>(new NSArgs))
{ }

// NavierStokes* 
// NavierStokes::fromFile(const std::string& filename)
// {
//   NavierStokes* problem = new NavierStokes();
//   problem->read(filename);
//   return problem;
// }

void
NavierStokes::assemble(Vec x, double t, Vec r, Mat J)
{
  // XXX test vector and matrix sizes against dofmap !!!

  Mesh::Base*   mesh   = *this->mesh;
  DofMap::Base* dofmap = *this->dofmap;
  NSArgs* args = reinterpret_cast<NSArgs*>(this->nsargs);

  //GLOBAL_MESH = mesh;
  
  args->pack(x, t, x, t, r, J, 1.0, 1);
  int ierr = ::assemble(mesh, args->argl, dofmap, "comp_mat_res", &args->time_star);
  args->clear();
  if (ierr) throw Error("NavierStokes::assemble()");
  VecScale(r, -1.0);
}

void
NavierStokes::assemble(Vec x0, double t0, Vec x1, double t1,
		       Vec r, Mat J, double alpha)
{
  // XXX test vector and matrix sizes against dofmap !!!

  /* test */
  if (t1 < t0 )    throw Error("invalid value, 't1' < 't0'");
  if (alpha < 0.0) throw Error("invalid value, 'alpha' < 0.0");
  if (alpha > 1.0) throw Error("invalid value, 'alpha' > 1.0");

  Mesh::Base*   mesh   = *this->mesh;
  DofMap::Base* dofmap = *this->dofmap;
  NSArgs* args = reinterpret_cast<NSArgs*>(this->nsargs);

  //GLOBAL_MESH = mesh;
  
  args->pack(x0, t0, x1, t1, r, J, alpha, 0);
  int ierr = ::assemble(mesh, args->argl, dofmap, "comp_mat_res", &args->time_star);
  args->clear();
  if (ierr) throw Error("NavierStokes::assemble()");
  VecScale(r, -alpha);
}


PYPF_NAMESPACE_END
