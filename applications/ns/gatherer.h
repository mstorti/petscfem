// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: gatherer.h,v 1.11 2002/05/17 00:37:59 mstorti Exp $
#ifndef GATHERER_H
#define GATHERER_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to compute integrals, or any kind of function
    that returns scalar values.
*/
class gatherer : public Elemset { 
public: 
  int gather_length;
  /// This should not be defined by the user
  ASSEMBLE_FUNCTION;
  /// Return which "jobinfos" will process
  ASK_FUNCTION;
  
  /** @name User defined callback functions
      @doc Deriving the class and defining these
      functions you can define several integrators.
  */
  //@{
  /// Called \textbf{before} the element loop. 
  virtual void init() { }

  /** Called for each element, before the Gauss point loop. 
      @param k (input) element number in the elemeset
  */ 
  virtual void element_hook(int k) {}

  /** Called at each Gauss point loop. Compute the contribution at
      each Gauss point. 
      @param pg_values (output) your contributions at the Gauss point
      @param u (input) the state of variables at time n+1
      @param uold (input) the state of variables at time n
      @param xpg (input) coordinates of the Gauss point
      @param n (input) normal to the surface (if #ndimel=ndim-1#)
      between real and master coordinates
      @param wpgdet (input) the weight of the Gauss point
      @param time (input) time value ($t^\alpha$)
      
  */ 
  virtual void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &n,
		     double wpgdet,double time)=0;

  /** Called \textbf{after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  //@}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the force that the container is
    performing on the wall. 
    # force = \int_\surf - p * normal \dS #, where
    #force# is the computed force, #surf# is the
    surface of the elemset, #p# is pressure #normal#
    is the unit vector \emph{exterior} to the fluid.
    #\dS# is the diferential os surface. 
*/ 
class force_integrator : public gatherer {
private:
  /** Auxiliary vectors to store force, moment, the point about
      which to compute moments, the displacemente vector from the
      moment center to the actual integration point .
  */
  FastMat2 force,moment,x_center,dx;
  /// Flag to compute moments or not, number of dimensions
  int compute_moment,ndim_m;
public:
  /// perform several checks and initialization
  void init();
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
  /// perform some cleaning
  void clean();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the flow rate at wall. 
    fixme:= agregar doc.
*/ 
class flow_rate_integrator : public gatherer {
private:
  FastMat2 Q;
  int ndim_m;
public:
  /// perform several checks and initialization
  void init();
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the total height of the free surface (linear free surface
    b.c. version).
    fixme:= agregar doc.
*/ 
class free_surface_level_integrator : public gatherer {
private:
public:
  /// perform several checks and initialization
  void init() { assert(gather_length==1); }
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time) {
    pg_values[0] = wpgdet * u.get(1);
  }
};

#endif
