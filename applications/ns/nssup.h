// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: nssup.h,v 1.5 2001/10/16 18:33:45 mstorti Exp $
#ifndef ROCKNS_H
#define ROCKNS_H

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/texthf.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include <src/pfmat.h>
#include <src/iisdmat.h>

#include <applications/ns/nsi_tet.h>

/// Restriction for linearized free surface boundary condition
class ns_sup_res : public LagrangeMult {
private:
  /// The dimension of the problem
  int ndim;
  /// The index (number of degree of freedom) for pressure
  int p_indx;
  /// Gravity acceleration
  double gravity;
  // Density
  double rho;
public:
  /// Number of restrictions
  int nres() {return 1;};
  /** Return the node/dof pair to be used as lagrange multiplier for
      the #jr#-th restriction. 
  */
  void lag_mul_dof(int jr,int &node,int &dof) {assert(jr==1); node=2; dof=2; };
  /// Initialize the elemset (maybe reads hash table)
  void init();
  /// computes the residual and jacobian of the function to be imposed
  void res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
	   FastMat2 & jac);
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class ns_sup : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

/// General restriction elemset
class linear_restriction : public LagrangeMult {
private:
  /// Coefficients of the restrictions. 
  vector<double> coef;
public:
  linear_restriction() : LagrangeMult() { coef=NULL;};
  /// Number of restrictions 
  int nres();
  /** Return the node/dof pair to be used as lagrange multiplier for
      the #jr#-th restriction. 
  */
  void lag_mul_dof(int jr,int &node,int &dof);
  /// Initialize the elemset (maybe reads hash table)
  void init();
  /// computes the residual and jacobian of the function to be imposed
  void res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
	   FastMat2 & jac);
};

#endif
