// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: advective.h,v 1.69 2004/09/30 16:52:35 mstorti Exp $
 
//#define CHECK_JAC // Computes also the FD Jacobian for debugging
 
#ifndef ADVECTIVE_H
#define ADVECTIVE_H

/** The jacobians of the flux functions. This is an array
    of ndim matrices of ndof x ndof entries each. 
*/
#define AJAC(jd) (*A_jac[(jd)-1])

/// Standard arguments for the typical flux function. Newmat version. 
#define FLUX_FUN_ARGS const RowVector &U,int ndim,const Matrix &iJaco,		\
	      Matrix &H, Matrix &grad_H,					\
	      Matrix &flux,vector<Matrix *> A_jac,				\
	      Matrix &A_grad_U, Matrix &grad_U,					\
	      Matrix &G_source,Matrix &tau_supg, double &delta_sc,		\
	      double &lam_max,							\
	      TextHashTable *thash,						\
              Matrix &nor,Matrix &lambda,Matrix &Vr,Matrix &Vr_inv,		\
              double *propel,							\
	      void *user_data,int options, int &start_chunk, int & ret_options

#define FLUX_FUN_ARGS_GENER(ML_) const ML_ &U,			\
              int ndim,const ML_ &iJaco,			\
	      ML_ &H, ML_ &grad_H,				\
	      ML_ &flux,vector<ML_ *> A_jac,			\
	      ML_ &A_grad_U, ML_ &grad_U,			\
	      ML_ &G_source,ML_ &tau_supg, double &delta_sc,	\
	      double &lam_max,					\
	      TextHashTable *thash,				\
              ML_ &nor,ML_ &lambda,ML_ &Vr,ML_ &Vr_inv,		\
              double *propel,					\
	      void *user_data,int options,			\
              int &start_chunk, int & ret_options

#define AD_FLUX_FUN_ARGS const  FastMat2  &U,					\
	      int ndim, const  FastMat2  &iJaco,  FastMat2  &H,			\
	      FastMat2  &grad_H,  FastMat2  &flux, FastMat2 &fluxd,		\
	      FastMat2  &A_jac,  FastMat2  &A_grad_U,				\
	      FastMat2  &grad_U,  FastMat2  &G_source,				\
              FastMat2  &D_jac, FastMat2 &C_jac,				\
	      FastMat2  &tau_supg, double &delta_sc, double &lam_max,		\
	      TextHashTable *thash,  FastMat2  &nor, FastMat2  &lambda,		\
	      FastMat2  &Vr, FastMat2  &Vr_inv, double *propel,			\
	      void *user_data,int options,int &start_chunk, int & ret_options

#define FLUX_FUN_ARGS_FM FLUX_FUN_ARGS_GENER(FastMat)

/// Standard arguments for the typical flux function (in the function call).
#define FLUX_FUN_CALL_ARGS_GENER(ML_) U##ML_,ndim,				\
              iJaco##ML_,H##ML_,grad_H##ML_,					\
              flux##ML_,A_jac##ML_,A_grad_U##ML_,				\
	      grad_U##ML_,G_source##ML_,					\
              tau_supg##ML_,delta_sc,lam_max,thash,				\
              nor##ML_,lambda##ML_,Vr##ML_,Vr_inv##ML_,propel,user_data,	\
	      options,start_chunk,ret_options

#define FLUX_FUN_CALL_ARGS_FM FLUX_FUN_CALL_ARGS_GENER(_FM)

#define FLUX_FUN_CALL_ARGS_FM2 FLUX_FUN_CALL_ARGS_GENER(_FM2)

#define FLUX_FUN_CALL_ARGS U,ndim,iJaco,H,grad_H,flux,A_jac,A_grad_U,	\
	      grad_U,G_source,tau_supg,delta_sc,lam_max,thash,		\
              nor,lambda,Vr,Vr_inv,propel,user_data,			\
	      options,start_chunk,ret_options
	      
/// Type of flux function. 
typedef int AdvDifFluxFunction (AD_FLUX_FUN_ARGS);

#define ADVDIFFF_ARGS const NewElemset *elemset,const  FastMat2  &U,		\
	   int ndim, const  FastMat2  &iJaco,  FastMat2  &H,			\
	   FastMat2  &grad_H,  FastMat2  &flux, FastMat2 &fluxd,		\
	   FastMat2  &A_jac,  FastMat2  &A_grad_U,				\
	   FastMat2  &grad_U,  FastMat2  &G_source,				\
	   FastMat2  &D_jac,	FastMat2  &C_jac,				\
	   FastMat2  &tau_supg, double &delta_sc, double &lam_max,		\
	   FastMat2  &nor, FastMat2  &lambda,					\
	   FastMat2  &Vr, FastMat2  &Vr_inv, double *propel,			\
	   void *user_data,int options,int &start_chunk, int & ret_options

/// Generic FastMat2 matrix
#if 0
class FastMat2Shell {
public:
  /// Performs A = (*this) * B;
  virtual apply(FastMat2 & A,FastMat2 & B) const;
}
#endif
typedef void FastMat2Shell(FastMat2 & A,FastMat2 & B);

#define COMPUTE_FLUX_ARGS						\
		const FastMat2 &U,					\
		const FastMat2 &iJaco, FastMat2 &H,			\
		FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,	\
		FastMat2 &A_grad_U,					\
		FastMat2 &grad_U, FastMat2 &G_source,			\
		FastMat2 &tau_supg, double &delta_sc,			\
		double &lam_max,					\
		FastMat2 &nor, FastMat2 &lambda,			\
		FastMat2 &Vr, FastMat2 &Vr_inv,				\
		int options

// This is the flux function for a given physical problem. 
class AdvDifFF {
private:
  // The list of variables to be logarithmically transformed
  vector<int> log_vars_v;
public:
  virtual int operator() (ADVDIFFF_ARGS) {
    PetscPrintf(PETSC_COMM_WORLD,"Undefined flux function\n");
    return 0;
  }
  /** Returns the list of variables that are 
      logarithmic transformed.
      @param nlog_vars (input) the number of logarithmic variables
      @param log_vars (input) the list of logarithmic variables
      @author M. Storti
  */ 
  virtual void get_log_vars(const NewElemset *elemset,int &nlog_vars, 
			    const int *& log_vars);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Virtual class thet defines the relation between vector
    state and enthalpy content. 
*/ 
class EnthalpyFun {
public:
  /// The actual state
  FastMat2 UU;
  /** Allows updating the data for the object. 
      @param e (input) cofficients for updating the object
  */ 
  virtual void update(const double *e) {}

  /** Allows setting the state of the object
      @param U (input) the state of the object
  */ 
  virtual void set_state(const FastMat2 &U) { UU.set(U); }

  /** Computes the enthalpy vector from the state vector
      @param H (output) the enthalpy content vector
      @param U (input) the state vector
  */ 
  virtual void enthalpy(FastMat2 &H,const  FastMat2 &U) {
    // fixme:= This may be wrong, we should first
    // `U.resize()' in a `init()' function and then `U.clear()' in a
    // `clear()' function 
    set_state(U);
    enthalpy(H);
  }

  /** Computes the enthalpy vector from the state vector.
      Uses the state previously set with `set_state'
      @param H (output) the enthalpy content vector
  */ 
  virtual void enthalpy(FastMat2 &H)=0;

  /** Computes the product #(W_Cp_N)_(p,mu,q,nu) = W_p N_q Cp_(mu,nu)#
      @param W_Cp_N (output) size #nel# x #ndof# x #nel# x #ndof#
      @param W (input) weight function, size #nel#
      @param N (input) interpolation function, size #nel#
      @param w (input) scalar weight
  */ 
  virtual void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
			   double w)=0;
  /** Computes the product #(P_Cp)_(mu,nu) = (P_supg)_(mu,lambda) 
      Cp_(lambda,nu)#
      @param P_Cp (output) size #ndof# x #ndof#
      @param P_supg (input) matricial weight function, size #ndof# x #ndof#
  */ 
  virtual void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg)=0;
};

/// Constant Cp for all fields
class GlobalScalarEF : public EnthalpyFun {
  /// Aux var. identity of size ndof
  FastMat2 eye_ndof,htmp1,htmp2;

protected:
  /// The actual Cp
  double Cp;

public:
  /// initializes dimensions and sets Cp
  void init(int ndim,int ndof,int nel,double Cp=1.);

  /// Sets Cp from elemset data
  void update(const double *Cp_) {Cp=*Cp_;};

  /// Multiplies U by Cp with 'scale'
  void enthalpy(FastMat2 &H);

  /// Scales at the same time by w*Cp
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);

  /// Scales #P_supg# by #Cp#
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);
};

/** Constant Cp=1 for all the fields. 
    Identity relation between H and T
*/
class IdentityEF : public GlobalScalarEF {
public:
  IdentityEF() { Cp=1.; }
  /// Does nothing
  void update(const double *Cp_) { }
  /// Copies #U# in #H#
  void enthalpy(FastMat2 &H) { H.set(UU); }
  /// Copies #P_supg# in #P_Cp#
  void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg) { P_Cp.set(P_supg); }
};

extern IdentityEF identity_ef;
class NewAdvDif;

/// This is the flux function for a given physical problem. 
class NewAdvDifFF {
private:
  /// The list of variables to be logarithmically transformed
  vector<int> log_vars_v;
public:
  /// The elemset associated with the flux function
  const NewElemset *elemset;
  /** Pointer to NewAdvDif elemset (FIXME:= `elemset' doesn't serve because
      it is a pointer to `NewElemset')
  */
  const NewAdvDif *new_adv_dif_elemset;
  /// The enthalpy function for this flux function
  EnthalpyFun *enthalpy_fun;
  /// Constructor from the elemset
  NewAdvDifFF(const NewElemset *elemset_=NULL);

  /** Define the list of variables that are 
      treated logarithmically. Reads from the options 
      #nlog_vars# and #log_vars#. 
  */
  virtual void get_log_vars(int &nlog_vars,const int *& log_vars);

  virtual void set_profile(FastMat2 &seed);

  /** This is called before any other in a loop and may help in
      optimization 
      @param ret_options (input/output) this is used by the flux
      function writer for returning some options. Currently the only
      option used is #SCALAR_TAU#. This options tells the elemset
      whether the flux function returns a scalar or matrix
      #tau_supg#. 
  */ 
  virtual void start_chunk(int &ret_options) =0;

  /** This is called before entering the Gauss points loop and may
      help in optimization. 
      @param element (input) an iterator on the elemlist. 
  */ 
  virtual void element_hook(ElementIterator &element) =0;

  /** @name Advective jacobians related */
  //@{
  /** Computes the product #(A_grad_N)_(p,mu,nu) = A_(i,mu,nu) (grad_N)_(i,p)#
      @param A_grad_N (output, size #nel# x #nd# x #nd#) 
      @param grad_N (input, size #nel# x #ndof#)
  */ 
  virtual void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) =0;

  /** Computes the product #(A_jac_n)_(mu,nu) = A_(i,mu,nu) normal_i#
      @param A_jac_n (output, size #ndof# x #ndof#) 
      @param normal (input, size #ndim#)
  */ 
  virtual void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) =0;

  /** This sets the local state of the flow and is called before to
      call all the enthalpy functions functions and the like (those
      that don't need the state gradient #grad_U#.  
      @param U (input) the local state of the fluid
  */
  virtual void set_state(const FastMat2 &U) {}

  /** This sets the local state of the flow and is called before to
      call all the `comp_grad_N_D_grad_N' functions and alike.
      @param U (input) the local state of the fluid
      @param grad_U (input) the local gradient of state of the fluid 
  */
  virtual void set_state(const FastMat2 &U,const FastMat2 &grad_U) {}

  /** Computes fluxes, upwind parameters etc...
      fixme:= more doc here ...
  */ 
  virtual void compute_flux(COMPUTE_FLUX_ARGS) =0;
  //@}

  /** @name Diffusive jacobians related */
  //@{
  /** Computes the product #(grad_N_D_grad_N)_(p,mu,q,nu) 
      = D_(i,j,mu,nu) (grad_N)_(i,p) (grad_N)_(j,q)#
      @param grad_N_D_grad_N (output) size #nel# x #ndof# x #nel# x #ndof# 
      @param grad_N (input) size #nel# x #ndof#
  */ 
  virtual void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				    FastMat2 & dshapex,double w) =0 ;
  //@}

  /** @name Reactive jacobians related */
  //@{
  /** Computes the product #(N_N_C)_(p,mu,q,nu) 
      = w C_(mu,nu) N_p N_q#
      @param N_N_C (output, size #nel# x #ndof# x #nel# x #ndof# 
      @param N (input) FEM interpolation function size #nel#
      @param w (input) a scalar coefficient
  */ 
  virtual void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w)=0;

  /** Computes the product #(N_P_C)_(mu,q,nu) 
      = w (P_supg)_(mu,lambda) C_(lambda,nu) N_q #
      @param N_P_C (output) size  #ndof# x #nel# x #ndof# 
      @param P_supg (input) SUPG perturbation function size #ndof# x #ndof#
      @param N (input) FEM interpolation function size #nel#
      @param w (input) a scalar coefficient
  */ 
  virtual void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			  FastMat2 &N,double w)=0;

  /** Computes the SUPG perturbation function from the gradient of the
      shape function.  If you want to define a formulation for the
      SUPG perturbation function that can't be cast in the standard
      form (i.e. it is not of the form #tau * A * grad_N# then you
      should define #tau# as scalar (may be 0) and then compute your
      own expression for #P_supg#. 
      @param P_supg (output) the SUPG perturbation function for all
      nodes at this Gauss point, has dimensions #nel# x #ndof# x
      #ndof# 
*/
  virtual void comp_P_supg(FastMat2 &P_supg);

  /** Returns the dimension of the element (May be different from the
      dimension space). For instance, a river may be a 1D elemset in a
      2D space. If this is equal to the space dimension (the most
      common case) then return -1. 
      @return the dimension of the advective elemset.
  */
  virtual int dim() const { return -1; }

  virtual void get_C(FastMat2 &C);

  virtual void get_Cp(FastMat2 &Cp);

  virtual void compute_delta_sc_v(FastMat2 &delta_sc_v);

  virtual void get_bcconv_factor(FastMat2 &bcconv_factor);

  /** Returns the Riemann Invariants and jacobians 
      for Adv-Diff absorbent condition
  */
  virtual void Riemann_Inv(const FastMat2 &U, const FastMat2 &normaln,
			   FastMat2 &Rie, FastMat2 &drdU, FastMat2 &C_);


  virtual ~NewAdvDifFF()=0;
  //@}
};

/** Generic elemset for advective diffusive problems. 
    Several physical problems may be solved by defining the
    corresponding flux function object. 
*/
class NewAdvDif : public NewElemset { 
protected:
  /** A pointer to the flux function. Different
      physics are obtained by redefining the flux function
  */
  NewAdvDifFF *adv_diff_ff;
  int volume_flag;
  double Volume,rec_Dt_m;
  /// Options returned by the flux function
  int ff_options;
  /// The actual time
  double time_m;
  /// All these are for checking advective jacobians with numerical
  int compute_fd_adv_jacobian, comp_checked, comp_total;
  double compute_fd_adv_jacobian_eps,
    compute_fd_adv_jacobian_rel_err_threshold,
    A_fd_jac_norm_max, A_fd_jac_norm_min, A_jac_norm_max, 
    A_jac_norm_min, A_jac_err_norm_max, A_jac_err_norm_min,
    A_rel_err_min, A_rel_err_max;

public:
  FastMat2 dshapex,Uo,Ao_grad_N,tau_supg,P_supg,grad_Uo;
  friend class NewAdvDifFF;
  /// Contructor from the pointer to the fux function
  NewAdvDif(NewAdvDifFF *adv_diff_ff_=NULL) :
    adv_diff_ff(adv_diff_ff_), volume_flag(0) {};
  /** Destructor. Destroys the flux function object. fixme:= Warning: this is
      not good!! We cannot destroy the flux function object here if it
      is built in the derived class, because it may happen, for
      instance, that we pass a pointer to a global object. To be fixed
      later... 
   */
  ~NewAdvDif() {delete adv_diff_ff;}
  
  /** Prepare variables for report of error on
      flux advective jacobians */
  void before_assemble(arg_data_list &arg_datav,Nodedata *nodedata,
		       Dofmap *dofmap, const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_data);

  /// Report erros on jacobian fluxes
  void after_assemble(const char *jobinfo);

  /// The assemble function for the elemset. 
  NewAssembleFunction new_assemble;
  /// The ask function for the elemset. 
  ASK_FUNCTION;
  double volume() const ;
  const FastMat2 *grad_N() const ;
  double rec_Dt() const { return rec_Dt_m; }
  const FastMat2 &Uold() const { return Uo; }
  const FastMat2 &grad_Uold() const { return grad_Uo; }

  /// axisymmetric key
  int axi;

  /** Returns the actual time. 
      @return actual time
  */ 
  double time() const { return time_m; };

};

/** This is the companion elemset to advdif that computes the boundary
    term when using the weak-form option. 
*/
class NewBcconv : public NewElemset { 
  /** A pointer to the flux function. Different
      physics are obtained by redefining the flux function
  */
  NewAdvDifFF *adv_diff_ff;
public:
  /// Contructor from the pointer to the fux function
  NewBcconv(NewAdvDifFF *adv_diff_ff_=NULL) : 
    adv_diff_ff(adv_diff_ff_) {};
  /// Destructor. Destroys the flux function object. 
  ~NewBcconv() { delete adv_diff_ff; }
  /// The assemble function for the elemset. 
  NewAssembleFunction new_assemble;
  /// The ask function for the elemset. 
  ASK_FUNCTION;
};

/** The class AdvDif is a NewElemset class plus a
    advdif flux function object.
*/
class AdvDif : public NewElemset { 
public:
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
  AdvDifFF *adv_diff_ff;
};

class BcconvAdv : public NewElemset { 
public:
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
  AdvDifFF *adv_diff_ff;
};

enum flux_fun_opt {
  DEFAULT     = 0x0000,
  COMP_SOURCE = 0x0001,
  COMP_UPWIND = 0x0002,
  COMP_EIGENV = 0x0004,
  SCALAR_TAU  = 0x0008,
  COMP_SOURCE_NOLUMPED = 0x0010,
  COMP_SOURCE_LUMPED = 0x0020,
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#define ADVECTIVE_ELEMSET(name)				\
FluxFunction flux_fun_##name;				\
							\
class volume_##name : public Advective {		\
public:							\
  volume_##name() {flux_fun=&flux_fun_##name;};		\
};							\
							\
class absorb_##name : public Absorb {			\
public:							\
  absorb_##name() {flux_fun=&flux_fun_##name;};		\
};							\
							\
class bcconv_adv_##name : public BcconvAdv {		\
public:							\
  bcconv_adv_##name() {flux_fun=&flux_fun_##name;};	\
}

#define ADVECTIVE_ELEMSET_FM2(name)			\
FluxFunctionFM2 flux_fun_##name;			\
							\
class volume_##name : public AdvectiveFM2 {		\
public:							\
  volume_##name() {flux_fun=&flux_fun_##name;};		\
};							\
class bcconv_adv_##name : public BcconvAdvFM2 {		\
public:							\
  bcconv_adv_##name() {flux_fun=&flux_fun_##name;};	\
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Add here declarations for further advective elemsets. 
/// Euler equations for inviscid (Gas dynamics eqs.)

#define ADVDIF_ELEMSET(name)				\
class name##_ff_t : public AdvDifFF {			\
public:							\
  int operator()(ADVDIFFF_ARGS);			\
};							\
							\
class advdif_##name : public AdvDif {			\
public:							\
  advdif_##name() {adv_diff_ff = new name##_ff_t;};	\
};							\
							\
class bcconv_adv_##name : public BcconvAdv {		\
public:							\
  bcconv_adv_##name() {adv_diff_ff = new name##_ff_t;};	\
};

ADVDIF_ELEMSET(advecfm2);	// linear advective diffusive 
ADVDIF_ELEMSET(burgers);	// 1D scalar Burgers equation
ADVDIF_ELEMSET(swfm2t);	        // shallow water 2d turbulent
//ADVDIF_ELEMSET(swfm1t);	        // shallow water 1d supg

class wall_swfm2t : NewElemset { 
public: 
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
};

class GenLoad;

/// Generic surface flux function (film function) element
class HFilmFun {
public:
  const FastMat2 &H,&H_out,&H_in;
  GenLoad * elemset;
  // for one layer
  virtual void q(FastMat2 &uin,FastMat2 &flux,FastMat2 &jacin);
  // for two layers
  virtual void q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		 FastMat2 &jacin,FastMat2 &jacout);
  virtual void init()=0;
  virtual void element_hook(ElementIterator &element)=0;
  HFilmFun(GenLoad *e);
};

/// Generic surface flux element
class GenLoad : public NewElemset { 
  /** These are to pass the state of the `H' quantities on the
      internal and external layers. `Hin' is an alias for `H'
  */
  FastMat2 H_m,H_out_m;
public: 
  const FastMat2 &H,&H_out,&H_in;
  GenLoad() : H_in(H_m), H(H_m), H_out(H_out_m) {}
  HFilmFun *h_film_fun;
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
  virtual ~GenLoad()=0;
};

/// Global parameters passed to the elemsets
struct GlobParam {
  /// trapezoidal temporal integarion rule parameter
  double alpha;
  /// time step
  double Dt;
  /// if set to one then Dt=infty
  int steady;
  /// State vectors;
  Vec x,xold;
  /// Current time
  TimeData *time;
  /// Constructor
  GlobParam() : alpha(0), Dt(1.), steady(0) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Transforms state vector from logarithmic. The indices of fields
    logarithmically tranformed are listed in #log_vars#. 
    @author M. Storti
    CORREGIR:=
    @param true_lstate (output) Transformed from logarithm
    to positive variable. 
    @param lstate (ouput) input state logarithmically transformed
    (only those fields in #log_vars#).
    @param nlog_vars (input) number of fields logarithmically
    transformed
    @param log_vars (input) list of fields logarithmically
    transformed.
*/ 
void log_transf(FastMat2 &true_lstate,const FastMat2 &lstate,
		const int nlog_vars,const int *log_vars);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Sets an error for the negative jacobian error case. 
    @param detjaco (input) the determinant of the jacobian
    @param elem (input) the element number */ 
void detj_error(double &detJaco,int elem);

/** Cutoff function. It is very near to ${\rm ctff(x)}\approx \rm tol$ for
    $x<0$ and ${\rm ctff}(x)=x$ for $x\gg \rm tol$. 
*/ 
double ctff(double x, double & diff_ctff, double tol);

#endif
