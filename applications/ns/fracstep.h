// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: fracstep.h,v 1.6.2.1 2002/07/15 23:16:57 mstorti Exp $
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
