// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: stream.h,v 1.19 2002/10/12 15:49:37 mstorti Exp $
#ifndef STREAM_H
#define STREAM_H

#include "advective.h"

#define GETOPT_PROP(type,name,default) elemset->get_prop(name##_prop,#name) //nd

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class AdvDifFFWEnth;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DummyEnthalpyFun : public EnthalpyFun {
  AdvDifFFWEnth *s;
public:
  DummyEnthalpyFun(AdvDifFFWEnth *s_) : s(s_) {}
  void set_state(const FastMat2 &U);
  void enthalpy(FastMat2 &H);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class AdvDifFFWEnth : public NewAdvDifFF {
  DummyEnthalpyFun ef;
  FastMat2 UU;
public:
  friend class DummyEnthalpyFun;
  AdvDifFFWEnth(const NewElemset *elemset_=NULL) 
    : NewAdvDifFF(elemset_), ef(this) { enthalpy_fun = &ef; }
  virtual void enthalpy(FastMat2 &H)=0;
  virtual void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w)=0;
  virtual void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Defines the shape of the channel
class ChannelShape {
protected:
  const NewElemset *elemset;
public:
  ChannelShape(const NewElemset *e) : elemset(e) {}

  /// All objects should be crated with this
  static ChannelShape *factory(const NewElemset *e);
  
  /// Initialize properties (perhaps) from the elemset table 
  virtual void init() {}

  /// Read local element properties
  virtual void element_hook(ElementIterator element) {}
  
  /** For a given water depth (with respect to the bottom of the
      channel) give the fluid area, cross sectional water-line and wet
      channel perimeter.
      @param u (input) the water depth
      @param A (ouput) the fluid cross sectional area
      @param w (ouput) the fluid cross sectional water line
      @param P (ouput) the fluid cross sectional perimeter
  */
  virtual void geometry(double u,double &A,double &w,double &P)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Rectangular shaped channel
class rect_channel : public ChannelShape {
  Property width_prop;
  double width;
public:
  rect_channel(const NewElemset *e) : ChannelShape(e) {}

  /// Initializes the object
  void init() { 
    //o Width of the channel
    GETOPT_PROP(double,width,<required>);
    elemset->get_prop(width_prop,"width"); 
  }

  /// Read local element properties
  void element_hook(ElementIterator element) {
    width = elemset->prop_val(element,width_prop); 
  }
  
  /** For a given water depth (with respect to the bottom of the
      channel) give the fluid area, cross sectional water-line and wet
      channel perimeter.
      @param u (input) the water depth
      @param area (ouput) the fluid cross sectional area
      @param wl_width (ouput) the fluid cross sectional water line
      @param perimeter (ouput) the fluid cross sectional perimeter
  */
  void geometry(double u,double &area,
		double &wl_width,double &perimeter) {
    assert(u>=0.);
    assert(u>=0.);
    area = u*width;
    wl_width = width;
    perimeter = width+2*u;
  }

};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Circular shaped channel
class circular_channel : public ChannelShape {
  Property radius_prop;
  double radius;
public:
  circular_channel(const NewElemset *e) : ChannelShape(e) {}

  /// Initializes the object
  void init() { 
    //o Radius of the channel
    GETOPT_PROP(double,radius,<required>); }

  /// Read local element properties
  void element_hook(ElementIterator element) {
    radius = elemset->prop_val(element,radius_prop); 
  }
  
  /** For a given water depth (with respect to the bottom of the
      channel) give the fluid area, cross sectional water-line and wet
      channel perimeter.
      @param u (input) the water depth
      @param area (ouput) the fluid cross sectional area
      @param wl_width (ouput) the fluid cross sectional water line
      @param perimeter (ouput) the fluid cross sectional perimeter
  */
  void geometry(double u,double &area,
		double &wl_width,double &perimeter) {
    assert(u>=0.);
    assert(u<radius);
    assert(u>=0.);
    double cos_phi = (radius-u)/radius;
    double phi = acos(cos_phi);
    double sin_phi = sin(phi);
    area = radius*radius*(phi-sin_phi*cos_phi);
    wl_width = 2*radius*sin_phi;
    perimeter = 2*phi*radius;
  }

};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Abstract class representing all friction laws
class FrictionLaw {
protected:
  const NewElemset *elemset;
public:
  /// fixme:= This should be private
  FrictionLaw(const NewElemset *e);

  /// All objects should be crated with this
  static FrictionLaw *factory(const NewElemset *e);
  
  /// Initialize properties (perhaps) from the elemset table 
  virtual void init() {}

  /// Read local element properties
  virtual void element_hook(ElementIterator element) {}

  /** Should return the volumetric flow `Q' for a given area `A' and
      the derivative `dQ/dA'. 
      @param area (input) tranversal area
      @param perimeter (input) wetted perimeter
      @param S (input) bottom slope
      @param Q (output) volumetric flow
      @param C (output) derivative of volumetric flow w.r.t. A,
      i.e. #C = dQ/dQ#
  */ 
  virtual void flow(double area,double perimeter,
		    double S,double &Q,double &C) const =0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// This implements the Chezy friction law
class Chezy : public FrictionLaw {
  /// Chezy friction coefficient property
  Property Ch_prop;
  /// Chezy friction coefficient value
  double Ch;
public:
  Chezy(const NewElemset *e) : FrictionLaw(e) {}

  // Initialize properties
  void init() { 
    //o Chezy roughness coefficient
    GETOPT_PROP(double,Ch,<required>); }

  /// Read local element properties
  void element_hook(ElementIterator element) {
    Ch = elemset->prop_val(element,Ch_prop);
  }

  /** Should return the volumetric flow `Q' for a given area `A' and
      the derivative `dQ/dA'. 
      @param area (input) tranversal area
      @param perimeter (input) wetted perimeter
      @param S (input) bottom slope
      @param Q (output) volumetric flow
      @param C (output) derivative of volumetric flow w.r.t. A,
      i.e. #C = dQ/dQ#
  */ 
  void flow(double area,double perimeter,
	    double S,double &Q,double &C) const {
    double m=1.5;
    double gamma = Ch*sqrt(S/perimeter);
    Q = gamma*pow(area,m);
    C = Q*m/area;
  }    
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// This implements the Manning friction law
class Manning : public FrictionLaw {
  /// Manning friction coefficient property
  Property roughness_prop;
  /// Manning friction coefficient value
  double roughness;
  /// Conversion factor
  double a_bar;
public:
  Manning(const NewElemset *e) : FrictionLaw(e),
  a_bar(1.) {}

  // Initialize properties
  void init();

  /// Read local element properties
  void element_hook(ElementIterator element) {
    roughness = elemset->prop_val(element,roughness_prop);
  }

  /** Should return the volumetric flow `Q' for a given area `A' and
      the derivative `dQ/dA'. 
      @param area (input) tranversal area
      @param perimeter (input) wetted perimeter
      @param S (input) bottom slope
      @param Q (output) volumetric flow
      @param C (output) derivative of volumetric flow w.r.t. A,
      i.e. #C = dQ/dQ#
  */ 
  void flow(double area,double perimeter,
	    double S,double &Q,double &C) const {
    double m = 5./3.;
    double gamma = a_bar/roughness*sqrt(S)/pow(perimeter,2./3.);
    Q = gamma*pow(area,m);
    C = Q*m/area;
  }    
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The flux function for flow in a channel with arbitrary shape and
    using the Kinematic Wave Model.
*/ 
class stream_ff : public AdvDifFFWEnth {
  /// The local depth of the fluid
  double u;
  /// The slope of the channel
  // double S;
  /** Local state of the fluid obtained from geometrical parameters
      from the water depth `u'
  */
  double area,wl_width,perimeter;
  /// The local wave velocity C:= dQ/dA
  double C;
  /// Type of friction law used
  FrictionLaw *friction_law;
  /// Pointer to the channel shape object
  ChannelShape *channel;
  /// Properties related to friction
  //   Property slope_prop;
public:
  stream_ff(const NewAdvDif *e);

  ~stream_ff();
  
  /** This is called before any other in a loop and may help in
      optimization 
      @param ret_options (input/output) this is used by the flux
      function writer for returning some options. Currently the only
      option used is #SCALAR_TAU#. This options tells the elemset
      whether the flux function returns a scalar or matrix
      #tau_supg#. 
  */ 
  void start_chunk(int &ret_options);
  
  /** This is called before entering the Gauss points loop and may
      help in optimization. 
      @param element (input) an iterator on the elemlist. 
  */ 
  void element_hook(ElementIterator &element);

  /** Basically stores `U(1)' in the water depth `u' and computes
      geometric parameters of the channel
      @param U (input) the state of the fluid
      @param grad_U (input) gradient of the state of the fluid
  */ 
  void set_state(const FastMat2 &U,const FastMat2 &grad_U);

  /** Basically stores `U(1)' in the water depth `u' and computes
      geometric parameters of the channel
      @param U (input) the state of the fluid
  */ 
  void set_state(const FastMat2 &U);

  /** @name Advective jacobians related */
  //@{
  /** Computes the product #(A_grad_N)_(p,mu,nu) = A_(i,mu,nu) (grad_N)_(i,p)#
      @param A_grad_N (output, size #nel# x #nd# x #nd#) 
      @param grad_N (input, size #nel# x #ndof#)
  */ 
  void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N);

  /** Computes the product #(A_jac_n)_(mu,nu) = A_(i,mu,nu) normal_i#
      @param A_jac_n (output, size #ndof# x #ndof#) 
      @param normal (input, size #ndim#)
  */ 
  void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal);

  /** Computes fluxes, upwind parameters etc...
      fixme:= include more doc here ...
  */ 
  void compute_flux(COMPUTE_FLUX_ARGS);
  //@}

  /** @name Diffusive jacobians related */
  //@{
  /** Computes the product #(grad_N_D_grad_N)_(p,mu,q,nu) 
      = D_(i,j,mu,nu) (grad_N)_(i,p) (grad_N)_(j,q)#
      @param grad_N_D_grad_N (output) size #nel# x #ndof# x #nel# x #ndof# 
      @param grad_N (input) size #nel# x #ndof#
  */ 
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & dshapex,double w);
  //@}

  /** @name Reactive jacobians related */
  //@{
  /** Computes the product #(N_N_C)_(p,mu,q,nu) 
      = w C_(mu,nu) N_p N_q#
      @param N_N_C (output, size #nel# x #ndof# x #nel# x #ndof# 
      @param N (input) FEM interpolation function size #nel#
      @param w (input) a scalar coefficient
  */ 
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);

  /** Computes the product #(N_P_C)_(mu,q,nu) 
      = w (P_supg)_(mu,lambda) C_(lambda,nu) N_q #
      @param N_P_C (output) size  #ndof# x #nel# x #ndof# 
      @param P_supg (input) SUPG perturbation function size #ndof# x #ndof#
      @param N (input) FEM interpolation function size #nel#
      @param w (input) a scalar coefficient
  */ 
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w);

  void enthalpy(FastMat2 &H);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);

  /** This stream elemset is essentially 1D.
      @return the dimension of this element that is 1
  */
  int dim() const { return 1; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// The `stream' (river or channel) element.
class stream : public NewAdvDif {
public:
  stream() :  NewAdvDif(new stream_ff(this)) {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Losses from a stream to the aquifer. Here the `in'
    side of the element represents the stream and
    the `out' side the aquifer. 
*/ 
class StreamLossFilmFun : public HFilmFun {
  /** The property corresponding to the resistance of the stream
      bottom surface
  */
  Property Rf_prop;
  /// The inverse of resistance of the stream bottom surface `k=1/Rf'
  double k;
  /// If set to 1 then Rf = infty
  int impermeable;
public:
  StreamLossFilmFun(GenLoad *e) : HFilmFun(e) {}
  // Only defines the double layer source term
  void q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
	 FastMat2 &jacin,FastMat2 &jacout);

  void init() { 
    //o _T: double
    // _N: Rf
    // _D: 1.
    // _DOC: Resistivity (including perimeter) of the stream to 
    //          loss to the aquifer. 
    // _END
    elemset->get_prop(Rf_prop,"Rf"); 
    int ierr;
    //o Flag whether the element is impermeable ($R_f\to\infty$) or not. 
    EGETOPTDEF_ND(elemset,int,impermeable,0);
    assert(ierr==0);
    if (impermeable) k=0.;
  }

  void element_hook(ElementIterator &element) {
    if (!impermeable) k = 1./elemset->prop_val(element,Rf_prop);
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class stream_loss : public GenLoad { 
public: 
  StreamLossFilmFun stream_loss_film_fun;
  stream_loss() : stream_loss_film_fun(this) {h_film_fun = &stream_loss_film_fun;};
  ~stream_loss() {};
};

#endif
