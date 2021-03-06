// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$
#ifndef PETSCFEM_NSI_TET_H  
#define PETSCFEM_NSI_TET_H

#ifdef USE_ANN
#include <ANN/ANN.h>			// ANN declarations
#endif
#include <vector>

#include <src/secant.h>
#include <src/pfobject.h>

extern int fractional_step;
extern int reuse_mat;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// This is the typical element for volume computations. The `ask'
// function is the same.
class ns_volume_element : public Elemset { 
public: 
  ASK_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
//  class nsi_tet : public ns_volume_element { 
//  public: 
//    ASSEMBLE_FUNCTION;
//  };

//  //-------<*>-------<*>-------<*>-------<*>-------<*>------- 
//  class nsi_tet_les : public ns_volume_element { 
//  public: 
//    ASSEMBLE_FUNCTION;
//  };

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_fm2 : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_comp : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_ther : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_asm : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_asm : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_asm_avgvol : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class ns_gasflow : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_keps : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_rot : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_ns_fm2 : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_nsi_tet_asm : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

#if 1
//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_nsi_tet_asm_avgvol : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};
#endif

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_nsther_fm2 : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_nsasm_fm2 : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_ns_gasflow : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class wall : public Elemset { 
private:
  int ndim;
public: 
  void initialize();
  void before_assemble(arg_data_list &arg_datav,Nodedata *nodedata,
		       Dofmap *dofmap, const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_data);
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

typedef pair<int,Elemset *> ElemToPtr;

typedef vector<ElemToPtr> ElemToPtrV;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class WallData {
private:
#ifdef USE_ANN
  /// The octree
  ANNkd_tree *kd_tree;                 // search structure
  /// The position of the points
  ANNpointArray data_pts;
#endif
  /// Number of points
  int npoints;
  /// number of spacial dimensions
  int ndim;
  /// pointers to the corresponding elemsets
  ElemToPtr *elemset_pointer;
  /// The length of the elemset_pointer array
  int nelemset;
public:
  /// constructor from a vector of coordinates, pointers and dimensions
  WallData();

  /// Dtor. 
  ~WallData();

  void init(vector<double> *data_pts_,vector<ElemToPtr>
	    *elemset_pointer,int ndim_);

  void clear();

#if USE_ANN
  /// find the nearest neighbor
  void nearest(const ANNpoint &point, Elemset *& elemset, int &elem, ANNidx &nn_idx,
	       ANNpoint &nn,ANNdist &dist);
#endif
  /// find the nearest neighbor (short version) only returns the index.
  void nearest(double *point,int &nn);
  /** Give wall element info. 
      Given the nearest neighbor index give the corresponding
      wall element info: elemset, element number, coords of
      center of wall element. 
      @author M. Storti
      @param nn (input) the index of the wall element in
      the WallData structure as returned by nearest()
      @param elemset (output) the `wall' elemset to which the element pertains
      @param elem (output) the wall element index (base 0)
      @param coords (output) the coords of the center of the wall element. 
  */
  void nearest_elem_info(const int nn, Elemset *& elemset, int &elem,
			 const double *& coords);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

void wall_fun(double yp,double &f,double &fprime);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// This is the generic wall function
class WallFun {
public:
  /// Initialize the object
  virtual void init() {};
  /** Compute the function $f = f(y^+)$ and $fp = f'(y^+)$. 
      #yp#=$y^+$ is the nondimesional coordinate normal to the wall. 
  */
  virtual void w(double yp,double &f,double &fp)=0;
  /// The destructor
  virtual ~WallFun()=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The standard wall function, composed of a linear laminar profile, 
    a buffer region and the logarithmic region. 
*/
class WallFunStd : public WallFun {
private:
  /// A pointer to the elemeset in order to get the physical data. 
  Elemset *elemset;
  /// Constants of the law
  double c1,c2;
public:
  /// The specific wall function
  void w(double yp,double &f,double &fprime);
  // assuming a wall law of the form f = 2.5*log(yplus) + 5.5
  // then if it has to be compatible with f = 1/Chi * log(E*yplus)
  // we find: Chi = 0.4; E = 9.025;
  // WallFunStd(double Chi_=0.4, double E_=9.025) : Chi(Chi_), E_star(E_) {};
  // WallFunStd() : Chi(Chi_), E_star(E_) {};
  // WallFunStd(Elemset *e) : elemset(e) {};
  /// Constructor
  WallFunStd(Elemset *e);
  /// Destructor
  ~WallFunStd() {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This class provides the solution of the equation for the friction
    velocity with the secant method. 
*/
class WallFunSecant : public Secant {
public:
  /// viscosity
  double nu;
  /// coordinate normal to the wall
  double y_wall;
  /// Velocity at the point near the wall
  double u;
  /// Density
  double rho;
  /// Specific wall function
  WallFun *wf;
  /// Provides the residual of the function to be solved
  double residual(double ustar,void *user_data=NULL);
  /// Constructor from the wall function
  WallFunSecant(WallFun *wf_) : wf(wf_) {};
  /** Solve the nonlinear equation in the friction velocity
      for a given velocity at the near wall. The wall law 
      $u/u_* = f(y_w u_* /\nu)$ gives a non-linear equation in $u_*$
      for a given $u$. So that we can eliminate $u=u(u_*)$. 
      @param u (input) velocity at the near wall
      @param ustar (output) friction velocity
      @param tau_w (input) traction (friction) at the wall
      @param yplus (output) non-dimensional normal coordinate at
      the near wall
      @param fwall (output) wall function at the near wall
      @param fprime (output) slope of the wall function at the near
      wall
      @param dustar_du (output) is $du/du_*$. 
  */ 
  void solve(double u_,double &ustar,double &tau_w,double &yplus,
	     double &fwall, double &fprime,
	     double &dustar_du);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
#define NonLinearRes NonLinearResNS
class NonLinearRes : public Elemset {
 public:
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
  /// Number of restrictions
  virtual int nres()=0;
  /// Initialize the elemset (maybe reads hash table)
  virtual void init()=0;
  /** Computes the residual and jacobian of the function to be
      imposed. Usually you derive #NonLinearRes# and instantiate this
      function that defines the restriction to be imposed. 
      @param k (input) element number in elemset
      @param U (input) state vector at both nodes
      @param r (output) a vector of length #nres# containing the
      residuals
      @param lambda (input) the state of the Lagrange multipliers 
      @param jac (output) the jacobian of the residuals with respect
      of the node state. (size #nres x ndof#)
  */ 
  virtual void res(int k,FastMat2 &U,FastMat2 & r,
		   FastMat2 & lambda,FastMat2 & jac)=0;
  /// Make it pure virtual. 
  virtual ~NonLinearRes()=0;
};

#define MAXPROP_WLR 10
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This elemset provides the restriction (via Lagrange multipliers)
    of the nonlinear Dirichlet type b.c. at the wall of the form $k_w
    = u_*^2/C_\mu$ and similarly for $\epsilon_w$. 
*/
class wall_law_res : public NonLinearRes {
  /// Index (field numbers) of $k$ and $\epsilon$, number of dimensions
  int nk,ne,ndim;
  /// Normal coordinate where the wall law is imposed (non-dimensional)
  double y_wall_plus;
  /// Normal coordinate where the wall law is imposed (dimensional)
  double y_wall;
  /// Wall function at the near wall. 
  double fwall;
  /// Slope of wall function at the near wall. 
  double fprime;
  /// Turbulence model constant 
  double C_mu;
  /// Turbulence model constant 
  double von_Karman_cnst;
  /// Viscosity
  double viscosity;
  /// Deactivate turbulence
  double turbulence_coef;
  /// Specific wall function
  WallFun *wf;
  /** Wall function with secant solver if use $y=$cnst
      instead of $y^+=$cnst
  */
  WallFunSecant *wfs;
  /// Index for #u_wall# property
  int u_wall_indx;
  /// Number of properties
  int nprops;
  /// Absolute velocity at the near wall 
  FastMat2 u_wall;
  /// Relative  velocity at the near wall 
  FastMat2 du_wall;
  /// Array of property indices 
  int elprpsindx[MAXPROP_WLR]; 
  /// Array of properties
  double propel[MAXPROP_WLR];
public:
  /// Number of restrictions (=2)
  int nres() {return 2;};
  /// Initialize the elemset (maybe reads hash table)
  void init();
  /// computes the residual and jacobian of the function to be imposed
  void res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
	   FastMat2 & jac);
  /// Contructor
  wall_law_res() {wf = new WallFunStd(this); wfs = new WallFunSecant(wf);};
  /// Destructor
  ~wall_law_res() {delete wf; delete wfs;};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Wall elemset. Imposes traction as a function
    of velocity vial universal wall law functions. 
*/ 
class wallke : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

/** Cutoff function used in turbulence calculations. It is very near
    to ${\rm ctff}(x)\approx \rm tol$ for $x<0$ and ${\rm ctff}(x)=x$
    for $x\gg \rm tol$.  
    @param x (input) the argument where to compute the cutoff fuction
    @param (output) the derivative of the cutoff function at $x$
    @param (input) cuttoff scale parameter. 
    @return the cutoff'ed value
*/ 
double ctff(double x, double & diff_ctff, double tol=1e-5);

/** Reads a #FastMat2# matrix from a #TextHashTable# entry. 
    @param thash (input) the options table
    @param s (input) the key in the table where to extract the matrix coefficients
    @param ndof (input) Let #coef[]# be the list of coefficients entered in
    the line and #ncoef# the length of #coef#. 
    If #ncoef=1# then #cond = v[0] * Id(ndof)#, else if #ncoef=ndof# then, 
    #cond = diag(v)#, else if #ncoef=ndof*ndof# then #cond(i,j) = v[ndof*i+j]# 
    (rowwise). Any other causes an error. 
    @param ndof (input) the number of fields
    @param cond (output) the matrix to be read
    @return a reference to the matrix. */ 
void read_cond_matrix(TextHashTable *thash, const char *s,
		      int ndof,FastMat2 &cond);

/** This is the factory of BasicObjects (Objects that are read from
    the user data file), specific for the NS module */ 
BasicObject_factory_t BasicObject_ns_factory;

/// Fixes the jacobian of the element. 
void detj_error(double &detJaco,int elem);

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Creates hooks depending on the name. 
    @param name (input) the name of the hook. 
    @return a pointer to the hook. */ 
Hook *ns_hook_factory(const char *name);
#endif

double smabs(double x);
double smabs(double x,double &dydx);

#endif
