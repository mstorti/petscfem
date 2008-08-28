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
  string potential_field;
  string charge_field;

  double* potn_ptr;
  double* chrg_ptr;
  FastMat2 potncol;
  FastMat2 chrgcol;
  FastMat2 grad_potn;
  FastMat2Tmp tmp;

public: 
  virtual void bf_init(Nodedata*);
  virtual void bf_eval_el(int iel);
  virtual void bf_eval_pg(FastMat2& shape,FastMat2& dshape,
			  FastMat2& body_force);
};

#endif
