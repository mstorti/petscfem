// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: condwall.h,v 1.1 2005/03/28 03:29:31 mstorti Exp $
#ifndef PETSCFEM_CONDWALL_H
#define PETSCFEM_CONDWALL_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>
#include "./nslagmul.h"
#include "./nsi_tet.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class cond_wall : public NSLagrangeMult {
private:
  int ndim,nel,ndof;
  int R;			// Resistance of the membrane
  // Property normal_prop;
public:
  // First two nodes are real nodes at both sides of the membrane. 
  // Other two nodes are lagrange multipliers. 
  cond_wall() {} 
  ~cond_wall() { } 
  int nres() { return ndof*2; }
  void lag_mul_dof(int jr,int &node,int &dof) {
    if (jr<=ndof) {
      node = 1; dof=jr;
    } else {
      node = 2; dof=jr-ndof;
    }
  }
  //element init
  void element_hook(ElementIterator &element);
  void lm_initialize();
  void init();
  void res(int k,FastMat2 &U,FastMat2 &r,
	   FastMat2 &w,FastMat2 &jac);
  void close() {}
};

#endif
