// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: adaptor.h,v 1.7 2002/08/07 16:17:08 mstorti Exp $
#ifndef ADAPTOR_H
#define ADAPTOR_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Generic class that allows the user not make explicitly the element
    loop, but instead he has to code the `element_connector' function
    tha computes the residual vector and jacobian of the element. 
*/
class adaptor : public ns_volume_element { 
public: 
  /// This should not be defined by the user...
  ASSEMBLE_FUNCTION;
  /// the reciprocal of the time step. (May be null for steady problems)
  double rec_Dt;
  /// The number of Gauss points
  int npg;
  /// The dimension of the element
  int ndim;
  /// The element number (may be used for printing errors, for instance)
  int elem;
  ///  Shape function
  FastMat2 shape;
  /// Gradient with respect to master element coordinates 
  FastMat2 dshapexi;
  /** Gradient with respect to global coordinates. It is only
      dimensioned but it should be computed by the user.
  */
  FastMat2 dshapex;
  /// Gauss points weights
  FastMat2 wpg;
  /// Parameters passed to the element from the main
  GlobParam *glob_param;
  /** User defined callback function to be defined by the
      user. Called \textbf{before} the element loop. 
  */
  virtual void init()=0;
  /** User defined callback function to be defined by the
      user. Called \textbf{after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  /** User defined callback function to be defined by the
      user. Called \textbf{after} the element loop. May be used for
      clean-up operations. 
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
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class adaptor_pg : public adaptor { 
  /** User defined callback function for the `adaptor' class. 
      Implemented in this class.
  */
  void init();
  /** User defined callback function for the `adaptor' class. 
      Implemented in this class.
  */
  void clean();
public: 
  FastMat2 Jaco,iJaco,grad_state_new_pg,grad_state_old_pg,
    state_old_pg,state_new_pg,res_pg,mat_pg,
    xpg;
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  virtual void elemset_init() {};
  virtual void elem_init() {};
  virtual void elemset_end() {};
  virtual void elem_end() {};
  /// Warning: this function should *accumulate* on `mat' and `res'
  virtual void pg_connector(const FastMat2 &xpg,
			    const FastMat2 &state_old_pg,
			    const FastMat2 &grad_state_old_pg,
			    const FastMat2 &state_new_pg,
			    const FastMat2 &grad_state_new_pg,
			    FastMat2 &res_pg,FastMat2 &mat_pg)=0;
};

#endif
