// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: embgath.h,v 1.20 2003/02/25 03:14:01 mstorti Exp $
#ifndef EMBGATH_H
#define EMBGATH_H

#include <src/surf2vol.h>

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to compute integrals, or any kind of function
    that returns a set of scalar values over surface
    elements but depending on the gradient of state function
    normal to the surface. */
class embedded_gatherer : public Elemset { 
  /** Contains the interpolation/integration
      object (defines the element geometry). */
  Surf2Vol *sv_gp_data;
  /// Number of Gauss points. 
  int npg;
  int nel_surf, nel_vol, layer, layers;
public: 
  /** Constructor. Initializes the pointer 
      to interpolation/integration element. */
  embedded_gatherer() : sv_gp_data(NULL) { }
  /** Destructor. Safely destroys #sv_gp_data# because
      initialized to null pointer. */
  ~embedded_gatherer() { delete sv_gp_data; }
  /** Finds the volume element for each surface element and 
      rotates it in order to give the face to standard position. */
  void initialize();
  /** Number of elements in the global vector 
      that the element will contribute. */
  int gather_length;
  /** Loops over elements and Gauss points. Computes state vectors, 
      gradients, shape functions, and passes them to user as
      callback functions in order to compute contributions. */
  ASSEMBLE_FUNCTION;
  /// Only accepts the #gather# jobinfo. 
  ASK_FUNCTION;
  
  /** @name User defined callback functions
      @doc Deriving the class and defining these
      functions you can define several integrators.
  */
  //@{
  /// Called \emph{before} the element loop. 
  virtual void init() { }

  /** Called for each element, \emph{before} the Gauss point loop. 
      @param k (input) element number in the elemeset. */ 
  virtual void element_hook(int k) {}

  /** Called at each Gauss point loop. Compute the contribution at
      each Gauss point. 
      @param pg_values (output) your contributions at the Gauss point
      @param u (input) the field state at time n+1
      @param uold (input) the field state at time n
      @param grad_u (input) the gradient of field state of variables at time n+1
      @param grad_uold (input) the gradient of field state of variables at time n
      @param xpg (input) coordinates of the Gauss point
      @param n (input) normal to the surface (if #ndimel=ndim-1#)
      between real and master coordinates
      @param wpgdet (input) the weight of the Gauss point
      @param time (input) time value ($t^\alpha$) */ 
  virtual void set_pg_values(vector<double> &pg_values,FastMat2 &u,
			     FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold,
			     FastMat2 &xpg,FastMat2 &n, double wpgdet,double time)=0;

  /** Called \textbf{after} the element loop. May be used for
      clean-up operations. */
  virtual void clean() {};
  //@}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Computes the force (viscous and pressure) for Navier-Stokes eqs.
class visc_force_integrator : public embedded_gatherer { 
private:
  /// Auxiliary matrix, contains the force on the Gauss point. 
  FastMat2 force, moment, x_center, dx, strain_rate,
    sigma, Omega, rigid_grad_u;
  /// Flag to compute moments or not, number of dimensions
  int compute_moment, ndim_m;
  /// Viscosity
  double viscosity;
public:
  /// Get dimension number, define #compute_moment#, viscosity
  void init();
  /** Adds force (pressure+viscous) for a given Gauss point. 
   FIXME:= enter doc... */ 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold, 
		     FastMat2 &xpg,FastMat2 &n, double wpgdet,double time);
};

#endif
