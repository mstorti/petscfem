// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
#ifndef PETSCFEM_ADAPTOR_H
#define PETSCFEM_ADAPTOR_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Generic class that allows the user not make explicitly the element
    loop, but instead he has to code the `element_connector' function
    that computes the residual vector and jacobian of the element. 
*/
class adaptor : public ns_volume_element { 
private:
  /// Flags whether the elements have been initialized or not
  int elem_init_flag, use_fastmat2_cache;
  arg_data_list *arg_data_vp;
  /// Index of element local to chunk
  int ielh;
  /// Perturbed index, if=-1: not in FDJ loop,
  /// if 0: reference state, if >0 
  int jpert;
public: 
  adaptor();
  /// This should not be defined by the user...
  ASSEMBLE_FUNCTION;
  
  bool EVAL_RES;
  bool EVAL_MAT;
  /// the reciprocal of the time step. (May be null for steady problems)
  double rec_Dt;
  /// The trapezoidal rule parameter.
  double alpha;
  /// The number of Gauss points
  int npg;
  /// The dimension of the space
  int ndim;
  /// The dimension of the element
  int ndimel;
  /// The element number (may be used for printing errors, for instance)
  int elem;
  /// Use FEM interpolation (shape, dshapexi members)
  int fem_interpolation;
  ///  Shape function at this GP (size #nel#). 
  FastMat2 shape;
  /// Gradient with respect to master element coordinates (size #ndimel x nel#) 
  FastMat2 dshapexi;
  /** Gradient with respect to global coordinates. It is only
      dimensioned but it should be computed by the user. (size #ndim x
      nel#) */
  FastMat2 dshapex;
  /// Gauss points weights (size #nel#)
  FastMat2 wpg;
  /// Parameters passed to the element from the main
  GlobParam *glob_param;
  /** User defined callback function to be defined by the
      user. Called {\bf before} the element loop. May be used for
      initialization operations. 
  */
  virtual void init() {}
  /** User defined callback function to be defined by the
      user. Called {\bf after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  /** Callback that computes the residual and jacobian 
      of the element contribution. 
      @param xloc (input) Coordinates of the nodes.
      @param state_old (input) The state at time $t^n$
      @param state_new (input) The state at time $t^{n+1}$
      @param res (output) residual vector
      @param mat (input) jacobian of the residual with respect to
      $x^{n+1}$
  */
  virtual void element_connector(const FastMat2 &xloc,
				 const FastMat2 &state_old,
				 const FastMat2 &state_new,
				 FastMat2 &res,FastMat2 &mat)=0;
  /** Callback similar to the #element_connector()#, the
      difference cames when the option #use_jacobian_fdj# to
      compute Jacobians of the residual using finite
      differences is activated. The idea is that each of the
      callbacks compute a term in the residual.  In the
      following we use the convention to call #resd#, #matd#
      the values returned by #element_connector()# and
      #resa#, #mata# those returned by
      #element_connector_analytic()#. The residual as seen
      from the FEM program is the sum of the two,
      i.e. #res=resd+resa#.  If the Jacobians are computed
      analytically (i.e.  #use_jacobian_fdj=0#) then the
      Jacobian will be simply #matd+mata#.  If the Jacobians
      are computed by finite differences (i.e.
      #use_jacobian_fdj=1#) then the Jacobian will be
      #matd_fdj+mata#, where #matd_fdj# is the Jacobian of
      the #resd# part computed by perturbation of the
      #element_connector()# analytic function.  Note that
      #matd# is irrelevant when #use_jacobian_fdj=1#.  The
      idea is that the computation of FD Jacobians is
      expensive, so that if some terms of the residual are
      expensive to compute but their Jacobians can be easily
      computed analytically, then it is cheaper to separate
      this part, and leave the FD computation only for the
      rest.  Another optimization may be to compute
      expensive stuff that do not depend on the state only
      once (see doc for #prtb_index()#).
      @param xloc (input) Coordinates of the nodes.
      @param state_old (input) The state at time $t^n$
      @param state_new (input) The state at time $t^{n+1}$
      @param res (output) residual vector
      @param mat (input) jacobian of the residual with respect to
      $x^{n+1}$
  */
  virtual void 
  element_connector_analytic(const FastMat2 &xloc,
                             const FastMat2 &state_old,
                             const FastMat2 &state_new,
                             FastMat2 &resa,FastMat2 &mata) { }
  /** This is called only once for each element after calling 
      initialize().  */ 
  virtual void element_init() { } 

  class ArgHandle {
    friend class adaptor;
  private:
    int index_m;
    ArgHandle(int indx) : index_m(indx) {}
  public:
    ArgHandle() : index_m(-1) {}
    int index() { return index_m; }
    bool operator==(ArgHandle b) const;
  };
  static ArgHandle NullArgHandle;
  size_t nargs() const;

  ArgHandle get_arg_handle(const string &key,
                           const char* errmess=NULL) const;

  void export_vals(ArgHandle h,double *vals,int s=-1); 

  void export_vals(ArgHandle h,FastMat2 &a); 

  void get_vals(ArgHandle h,double *vals,int s=-1); 

  void get_vals(ArgHandle h,FastMat2 &a); 

  /** Returns the stage we are performinf when computing 
      Jacobians by finite differences. The stage number 
      returned may be #j=-1# (no stage), #j=0# (initialization, 
      reference values), #1<=j<=nen# computing residual 
      for the perturbation of the #j# component of 
      the element state vector. 
      @return stage number */ 
  int prtb_index();

  void after_assemble(const char *jobinfo);

  virtual void before_chunk(const char *jobinfo) { }
  virtual void after_chunk(const char *jobinfo) { }
  
  /** Contains constant fields for the element */
  FastMat2 Hloc;
  int nH;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to define elements only by its contributions at Gauss
    points (GP). The user is passed the GP coordinates, the state and
    gradient of state at the GP, for both the state at this time step
    and the next time step, and should return the residual and
    contribution to the Jacobian matrix. One has access also to the
    shape function, This adaptor may be used also for #ndimel<ndim#,
    for instance quad panels in 3D. In that case the gradients of the
    states and shape functions are also of size #ndim*...# and are
    parallel to the surface. Also the user has acces to the normal to
    the element.  */
class adaptor_pg : public adaptor { 
private:
  /** User defined callback function for the `adaptor' class. 
      Implemented in this class. */
  void init();
  /** User defined callback function for the `adaptor' class. 
      Implemented in this class. */
  void clean();
  /** User defined callback function for the `adaptor' class. 
      Implemented in this class. */
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  FastMat2 Jaco,grad_state_new_pg,grad_state_old_pg,
    state_old_pg,state_new_pg,res_pg,mat_pg,
    xpg,normal_m,g,ig,tmp,shape_m,dshapexi_m;
  int compute_H_fields;
public: 
  /** @name Call back functions. */ 
  //@{
  /** Callback hook to be executed before a chunk of elements */ 
  virtual void elemset_init() {};
  /** Callback hook to be executed before a specific element */ 
  virtual void elem_init() {};
  /** Callback hook to be executed after a specific element */ 
  virtual void elem_end() {};
  /** Callback hook to be executed after a chunk of elements */ 
  virtual void elemset_end() {};
  /** Callback function that defines the residual and matrix at the Gauss point. This 
      shouldn't scale by the Gauss weight, neither by the Gauss point volume. 
      Also this function should {\bf not} accumulate on #mat# and #res#.
      It should {\bf set} those variables. (Eventually reset to 0.)
      @param xpg (input) coordinates of the Gauss point (size #ndim#)
      @param state_old_pg (input) state vector at the GP (size #ndof#)
      @param grad_state_old_pg (input) gradient of state vector at the
      GP (size #ndim*ndof#). 
      @param state_new_pg (input) state vector at the GP (size #ndof#)
      @param grad_state_new_pg (input) gradient of state vector at the
      GP (size #ndim*ndof#). 
      @param res_pg (output) Residual computed by the routine. (size #nel*ndof#)
      @param mat_pg (output) Jacobian of the residual (sign reverted, i.e. #-dR/dU#),
      (size #nel*ndof*nel*ndof#). */ 
  virtual void pg_connector(const FastMat2 &xpg,
			    const FastMat2 &state_old_pg,
			    const FastMat2 &grad_state_old_pg,
			    const FastMat2 &state_new_pg,
			    const FastMat2 &grad_state_new_pg,
			    FastMat2 &res_pg,FastMat2 &mat_pg)=0;
  //@}
  /** Returns the normal to the element (in the case #ndimel<ndim#).  */ 
  FastMat2 &normal();
  FastMat2 &shape();
  FastMat2 &dshapexi();
  FastMat2 &dshapex();
  FastMat2 H, grad_H;
};

#define GET_JOBINFO_FLAG_ND(name) name  = !strcmp(jobinfo,#name)

#endif
