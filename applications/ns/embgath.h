// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: embgath.h,v 1.1 2002/08/06 14:28:46 mstorti Exp $
#ifndef EMBGATH_H
#define EMBGATH_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to compute integrals, or any kind of function
    that returns scalar values.
*/
class embedded_gatherer : public Elemset { 
public: 
  void initialize();
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

  /** Called for each element, \emph{before} the Gauss point loop. 
      @param k (input) element number in the elemeset
  */ 
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
      @param time (input) time value ($t^\alpha$)
  */ 
  virtual void set_pg_values(vector<double> &pg_values,FastMat2 &u,
			     FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold, 
			     FastMat2 &xpg,FastMat2 &n,
			     double wpgdet,double time)=0;

  /** Called \textbf{after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  //@}
};

#endif
