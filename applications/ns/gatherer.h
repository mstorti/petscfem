// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: gatherer.h,v 1.1 2002/03/16 21:36:43 mstorti Exp $
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

  /** @name Call back functions
      @doc Deriving the class and defining these
      functions you can define several integrators.
  */
  //@{



  
  /** User defined callback function to be defined by the
      user. Called \textbf{before} the element loop. 
  */
  virtual void init()=0;

  /** Let the user do some things for the element (common for all
      Gauss Points)
      @param k (input) element number in the elemeset
  */ 
  void element_hook(k);

  /** User defined callback function to be defined by the
      user. Called \textbf{after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  //@}
};


#endif
