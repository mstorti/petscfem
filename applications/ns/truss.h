// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id mstorti-v6-branch-1.0.2-8-ge5d52ed Sun Oct 14 10:36:32 2007 -0300$

#ifndef TRUSS_H
#define TRUSS_H

class  truss : public adaptor { 
private:
  /// Auxiliary variables
  FastMat2 G;
public: 
  /** Initializes the elemset. Reads parameters, 
      resizes auxiliary matrices. 
  */ 
  void init();
  /** Computes the residual (gradient of the distortion functional) and 
      the jacobian (Hessian, i.e. second derivatives, of the 
      distortion functional. 
      @param xloc (input) element node coordinates
      @param state_old (input) element node displacements (previous step)
      @param state_new (input) element node displacements (actual step)
      @param res (output) gradient of distortion functional
      @param mat (input) Hessian (matrix of second derivatives) 
      of distortion functional 
   */ 
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
