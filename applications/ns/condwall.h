// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: condwall.h,v 1.6 2005/04/01 06:25:18 mstorti Exp $
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

struct cond_wall_data_t {
  dvector<double> 
  Rv,u1,u2;			// Resistance vector,
				// velocity on both layers. 
};

typedef map<string,cond_wall_data_t> 
cond_wall_data_map_t;

extern cond_wall_data_map_t cond_wall_data_map;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class cond_wall : public NSLagrangeMult {
private:
  int ndim,nel,ndof;
  int use_vector_resistance;
  double R;			// Resistance of the membrane
  // Property normal_prop;
  FastMat2 U1,U2,u1,u2;
  cond_wall_data_t *data_p;
public:
  // First two nodes are real nodes at both sides of the membrane. 
  // Other two nodes are lagrange multipliers. 
  cond_wall() : data_p(NULL) { } 
  ~cond_wall() { } 
  int nres();
  void lag_mul_dof(int jr,int &node,int &dof);
  //element init
  void element_hook(ElementIterator &element);
  void lm_initialize();
  void init();
  void res(int k,FastMat2 &U,FastMat2 &r,
	   FastMat2 &w,FastMat2 &jac);
  void close() {}
};

#endif
