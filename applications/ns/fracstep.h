// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: fracstep.h,v 1.9 2004/12/21 12:20:40 mstorti Exp $
#ifndef FRACSTEP_H
#define FRACSTEP_H

#include "nsi_tet.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class fracstep : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class fracstep_fm2 : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_fstep_fm2 : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

#endif
