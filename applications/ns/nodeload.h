// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$

#ifndef NODELOAD_H
#define NODELOAD_H

class nodeload : public adaptor { 
private:
  /// Auxiliary variables
  dvector<int> elprpsindx; 
  int nprops;
  dvector<double> propel;
  int load_indx;
  double load_fac;

public: 
  nodeload();
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
