// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: fracstep.h,v 1.7 2002/07/17 02:55:01 mstorti Exp $
#ifndef FRACSTEP_H
#define FRACSTEP_H

#include "nsi_tet.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class fracstep : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

#endif
