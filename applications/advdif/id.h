// -*- mode: C++ -*- 
// $Id: id.h,v 1.2 2003/10/16 19:13:42 mstorti Exp $
#ifndef ID_H
#define ID_H

#include "advective.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// This is the Identity operator, useful for debugging 
class id : public Elemset { 
public: 
  ASSEMBLE_FUNCTION;
};

#endif
