// -*- mode: C++ -*- 
// $Id: id.h,v 1.3 2003/12/08 23:24:45 mstorti Exp $
#ifndef ID_H
#define ID_H

#include "advective.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// This is the Identity operator, useful for debugging 
class id : public Elemset { 
public: 
  ASSEMBLE_FUNCTION;
  ASK_FUNCTION;
};

#endif
