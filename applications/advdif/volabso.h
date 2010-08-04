// -*- mode: C++ -*- 
// $Id: id.h,v 1.4 2005/02/23 01:40:34 mstorti Exp $
#ifndef PETSCFEM_VOLABSO_H
#define PETSCFEM_VOLABSO_H

#include "advective.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// This is the Identity operator, useful for debugging 
class volabso : public Elemset { 
public: 
  ASSEMBLE_FUNCTION;
  ASK_FUNCTION;
};

#endif
