// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: nsikepsrot.h,v 1.2 2002/05/04 23:54:13 mstorti Exp $
#ifndef NSIKEPSROT_H
#define NSIKEPSROT_H

#include "nsi_tet.h"

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_keps_rot : public ns_volume_element { 
private:
  int non_inertial_frame, part_include_fic;
public: 
  void initialize();
  int real_nodes(int iele,const int *&node);
  ASSEMBLE_FUNCTION;
};

#endif
