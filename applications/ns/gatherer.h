// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: gatherer.h,v 1.2 2002/03/17 03:53:40 mstorti Exp $
#ifndef GATHERER_H
#define GATHERER_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to compute integrals, or any kind of function
    that returns scalar values.
*/
class gatherer : public Elemset { 
private:
  // The coordinates of the Gauss point
  FastMat2 xpg; 
public: 
  /// This should not be defined by the user
  ASSEMBLE_FUNCTION;
  /// Return which "jobinfos" will process
  ASK_FUNCTION;
  
  /// Return the coordinate of the actual Gauss point
  const FastMat2 & x_pg() { return xpg; }

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
      @param Jaco (input) Jacobians of coordinates transformation
      between real and master coordinates
      @param wpgdet (input) the weight of the Gauss point
      @param time (input) time value ($t^\alpha$)
      
  */ 
  virtual void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time)=0;

  /** Called \textbf{after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  //@}
};

class force_integrator : public gatherer {
private:
  FastMat2 XCG,normal;
public:
  void init() { XCG.resize(1,3); normal.resize(1,3); }
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

#endif
