//__INSERT_LICENSE__
// $Id: advabso.cpp,v 1.23 2007/01/30 19:03:44 mstorti Exp $
#include "./advabso.h"
#include "./gasflow.h"

static
double msign(double x) {
  return (x<0. ?  1. : 0.0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbsoWall::
lm_initialize() { 
  int ff_options=0;
  adv_diff_ff->start_chunk(ff_options);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
TurnWallFun::~TurnWallFun() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class MyTurnWallFun : public TurnWallFun {
private:
  AdvectiveAbsoWall *elemset;
public:
  void 
  init(AdvectiveAbsoWall *elemset_a) {
    elemset = elemset_a;
  }
  int 
  is_wall(int elem,int node,FastMat2 &x,double t) {
    double y = x.get(2);
    return y>0.5;
  }
  ~MyTurnWallFun() { }
} my_turn_wall_fun;


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbsoWall::
init() {
  int nelprops;
  elem_params(nel,ndof,nelprops);
  int ierr;

  //o Dimension of the space.
  NSGETOPTDEF_ND(int,ndim,0);

  //o Flags whether to use the old state ad the
  // boundary as reference for the linear absorbing
  // boundary condition.
  NSGETOPTDEF_ND(int,use_old_state_as_ref,0);

  //o Use special combination for choosing reference state.
  //  If velocity is outgoing ( #uref.n>0#, #n# is exterior
  //  normal, #u# is the flow velocity), then the
  //  #use_old_state_as_ref# strategy is used, otherwise the
  //  state reference is used. 
  NSGETOPTDEF_ND(int,switch_to_ref_on_incoming,0);

  //o Do correction for wave characteristic computation
  //  do to mesh velocity (ALE).
  NSGETOPTDEF_ND(int,ALE_flag,0);

  //o Activate the transition from absorbing boundary condition
  //  to `wall' boundary condition
  NSGETOPTDEF_ND(int,activate_turn_wall,1);

  //o Turn absorbing boundary condition to `wall' boundary condition
  NSGETOPTDEF_ND(int,vel_indx,-1);
  PETSCFEM_ASSERT(!activate_turn_wall 
                  || vel_indx>=1 && vel_indx+ndim<=ndof,
                  "vel_indx should point to a valid dof range"
                  " if turn_wall is activated."
                  "vel_indx %d, ndim %d, ndof %d\n",
                  vel_indx,ndim,ndof);  
  
  vector<double> urefv;
  const char *line;
  get_entry("Uref",line);
  use_uref_glob = 0;
  if(line) {
    use_uref_glob = 1;
    read_double_array(urefv,line);
    PETSCFEM_ASSERT0(urefv.size() == (unsigned int)ndof,
		     "Uref needs ndof values \n");
    Uref_glob.resize(1,ndof);
    Uref_glob.set(&*urefv.begin());
  } 
  if (use_uref_glob || use_old_state_as_ref) assert(nel>=2);
  else assert(nel==3);

  ndimel = adv_diff_ff->dim();
  if (ndimel<0) ndimel = ndim;

  flux.resize(2,ndof,ndimel);
  fluxd.resize(2,ndof,ndimel);
  A_grad_U.resize(1,ndof);
  grad_U.resize(2,ndimel,ndof).set(0.);
  normal.resize(1,ndimel);
  vmesh.resize(1,ndimel);
  A_jac.resize(2,ndof,ndof);
  S.resize(2,ndof,ndof);
  invS.resize(2,ndof,ndof);
  tmp1.resize(2,ndof,ndof);
  c.resize(2,2,ndof);
  Pi_m.resize(2,ndof,ndof);
  Pi_p.resize(2,ndof,ndof);
  Cp.resize(2,ndof,ndof);
  Uold.resize(2,nel,ndof);
  invCp.resize(2,ndof,ndof);
  Ufluid.resize(1,ndimel);
  x.resize(1,ndim);
  if (ALE_flag) {
    get_prop(vmesh_prop,"vmesh");
    assert(vmesh_prop.length == ndimel);
  }

  get_prop(normal_prop,"normal");
  assert(normal_prop.length == ndimel);
  // The state on the reference node
  Uref.resize(1,ndof);
  dU.resize(1,ndof);
  // The lagrange multipliers (state for the 2nd node)
  Ulambda.resize(1,ndof);
  // The state of the outlet node
  Uo.resize(1,ndof);
  mask.resize(1,ndof)
    .set(0.0).is(1,vel_indx,vel_indx+ndim-1).set(1.0).rs();
  rlam.resize(1,ndof);
  node_list.resize(nel);
  turn_wall_fun = &my_turn_wall_fun;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbsoWall::
res(int k,FastMat2 &U,FastMat2 &r,
    FastMat2 &w,FastMat2 &jac) {
  double delta_sc=0.0,
      lambda_max_pg=0.0;
  U.ir(1,1); Uo.set(U);
  U.ir(1,2); Ulambda.set(U); U.rs();

  get_xloc(xloc,Hloc);

  xloc.ir(1,1);
  x.set(xloc);
  xloc.rs();

  Hloc.ir(1,1);
  H.set(Hloc);
  Hloc.rs();

  get_connect(node_list);

  int is_wall = 
    turn_wall_fun->is_wall(k+1,node_list[0],
                           x,get_time());
  
  FastMat2::branch();
  if (!is_wall) {
    FastMat2::choose(0);
    // As `use_old_state_as_ref' but
    // for this element. 
    int use_old_state_as_ref_elem;
    if (switch_to_ref_on_incoming) {

      // Flux Functions writters must add setUfluid funct
      // in order to extract the fluid velocities from Uref
      // U(3,:) contains the reference value
      // or alternatively the global `Uref' option
    
      if (use_uref_glob) Uref.set(Uref_glob);
      else { U.ir(1,3); Uref.set(U); U.rs(); }
      adv_diff_ff->set_Ufluid(Uref,Ufluid);
      double urefn = unor.prod(Ufluid,normal,-1,-1);
      use_old_state_as_ref_elem = urefn>0;
    } else {
      use_old_state_as_ref_elem 
        = use_old_state_as_ref; 
    }

    FastMat2::branch();
    if (use_old_state_as_ref_elem) {
      FastMat2::choose(0);
      get_old_state(Uold);
      Uold.ir(1,1);
      Uref.set(Uold);
      Uold.rs();
    } else {
      FastMat2::choose(1);
      if (use_uref_glob) Uref.set(Uref_glob);
      else { U.ir(1,3); Uref.set(U); U.rs(); }
    }
    FastMat2::leave();
  
    adv_diff_ff->set_state(Uref,grad_U);
    adv_diff_ff
      ->compute_flux(Uref, dummy, dummy, dummy, flux, fluxd,
                     A_grad_U, grad_U, dummy,
                     dummy, delta_sc, lambda_max_pg, dummy,
                     dummy, dummy, dummy, 0);
    adv_diff_ff->comp_A_jac_n(A_jac,normal);
    adv_diff_ff->get_Cp(Cp);
    invCp.inv(Cp);
    if (ALE_flag) {
      vnor.prod(vmesh,normal,-1,-1);
      // A_jac.d(1,2).add(-double(vnor)).rs();
      double vn = double(vnor);
      A_jac.axpy(Cp,-vn);
    }
    // tmp1 = Cp \ A
    tmp1.prod(invCp,A_jac,1,-1,-1,2);
    c.eig(tmp1,S);
    invS.inv(S);
    double aimag = c.ir(1,2).sum_square_all();
    assert(aimag<1e-10);
    c.ir(1,1);
    // Pi_m = projector on negative eigenvalues space
    // Pi_p = projector on positive eigenvalues space
    // Pi_m + Pi_p = A_jac
    Pi_m.set(0.).d(1,2)
      .set(c).fun(msign).rs();
    c.rs();
    tmp1.prod(Pi_m,invS,1,-1,-1,2);
    Pi_m.prod(S,tmp1,1,-1,-1,2);
#if 0
    tmp1.prod(Pi_p,invS,1,-1,-1,2);
    Pi_p.prod(S,tmp1,1,-1,-1,2);
    tmp1.set(Pi_m).add(Pi_p);
#endif
    // residual is the projection of U-Uref
    // on to the space of incoming waves
    dU.set(Uo).rest(Uref);
    r.prod(Pi_m,dU,1,-1,-1);
    // The vector of reactions is the pojector on
    // to the incoming wave space: w = Cp * Pi_m
    tmp1.prod(Cp,Pi_m,1,-1,-1,2);
    w.set(0.).ir(1,1).set(tmp1).rs();
    jac.ir(2,1).set(Pi_m);
#if 0
    FastMat2::branch();
    if (nel>=3 && !use_old_state_as_ref_elem) {
      FastMat2::choose(0);
      jac.ir(2,3).set(Pi_m).scale(-1.0);
    }
    FastMat2::leave();
#endif
    jac.rs();
  } else {
    FastMat2::choose(1);
    r.set(Uo).mult(mask);
    rlam.set(1.0).rest(mask).mult(Ulambda);
    r.add(rlam);
    jac.ir(2,1).d(1,3).set(mask);
    jac.rs().ir(2,2).d(1,3).set(1.0).rest(mask).rs();
    w.ctr(jac,3,1,2);
  }
  FastMat2::leave();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbsoWall::
element_hook(ElementIterator &element) {
  normal.set(prop_array(element,normal_prop));
  if (ALE_flag)
    vmesh.set(prop_array(element,vmesh_prop));
  adv_diff_ff->element_hook(element);
}
