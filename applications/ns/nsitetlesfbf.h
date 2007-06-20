// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id$
#ifndef PETSCFEM_NSITETLESFBF_H
#define PETSCFEM_NSITETLESFBF_H

#include "src/fm2temp.h"
#include "./nsitetlesf.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_full_bf : nsi_tet_les_full { 
public: 
  double param1;
  double param2;

public: 
  double* field1;
  double* field2;
  string pot1name;
  string pot2name;
  vector<double> pot1vec;
  vector<double> pot2vec;
  FastMat2 pot1col;
  FastMat2 pot2col;
  FastMat2 grad_pot1;
  FastMat2 grad_pot2;
  FastMat2Tmp tmp;

public: 
  virtual void bf_init(Nodedata*);
  virtual void bf_eval_el(int iel);
  virtual void bf_eval_pg(FastMat2& shape,FastMat2& dshape,
			  FastMat2& body_force);
};

#endif
