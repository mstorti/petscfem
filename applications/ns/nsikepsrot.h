// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: nsikepsrot.h,v 1.1 2002/05/04 23:32:19 mstorti Exp $
#ifndef NSIKEPSROT_H
#define NSIKEPSROT_H

#include "nsi_tet.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_keps_rot : public ns_volume_element { 
public: 
  int real_nodes(int iele,const int *&node);
  ASSEMBLE_FUNCTION;
};

#endif
