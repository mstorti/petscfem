// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: condwallpen.h,v 1.3 2005/04/09 15:49:04 mstorti Exp $
#ifndef PETSCFEM_CONDWALLPEN_H
#define PETSCFEM_CONDWALLPEN_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/penalize.h>
#include "./nspenal.h"
#include "./nsi_tet.h"
#include "./condwall.h"

class CondWallRestriction : public Restriction {
private:
  int nel,ndof,ndim;
  int use_vector_resistance;
  double R;			// Resistance of the membrane
  cond_wall_data_t *data_p;
  FastMat2 U1,U2,u1,u2;
public:
  CondWallRestriction() : data_p(NULL) { }
  int init(int nel_a,int ndof_a,
	   TextHashTable *thash,const char *name);
  void res(int k,FastMat2 &U,FastMat2 & r,
	   FastMat2 & w,FastMat2 & jac);
  ~CondWallRestriction() { }
};

class cond_wall_pen : public NSPenalize {
public:
  // First two nodes are real nodes at both sides of the membrane. 
  // Other two nodes are lagrange multipliers. 
  cond_wall_pen() : NSPenalize(new CondWallRestriction){ }
};
  
#endif
