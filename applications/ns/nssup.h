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

class ns_sup_res : public LagrangeMult {
private:
  double gravity, rho;
public:
  /// Number of restrictions per node pair
  int nres() {return 1;};
  /** Return the node/dof pair to be used as lagrange multiplier for
      the #jr#-th restriction. 
  */
  void lag_mul_dof(int jr,int &node,int &dof) {node=2; dof=2;};
  /// Initialize the elemset (maybe reads hash table)
  void init();
  /// computes the residual and jacobian of the function to be imposed
  void res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
	   FastMat2 & jac);
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class ns_sup : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

#endif
