// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: nsid.h,v 1.1 2002/10/07 00:26:08 mstorti Exp $
#ifndef NSID_H
#define NSID_H

#include "nsi_tet.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// This is the Identity operator, useful for debugging 
class ns_id : public ns_volume_element { 
private:
  int part_include_fic, real_nodes_m;
public: 
  void initialize();
  int real_nodes(int iele,const int *&node);
  ASSEMBLE_FUNCTION;
};

#endif
